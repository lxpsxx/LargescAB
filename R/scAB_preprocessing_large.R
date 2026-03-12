choose_scab_graph_mats <- function(
  A,
  prefer_author = TRUE,
  check_numeric_sym = TRUE,
  sym_tol = 1e-12,
  force_symmetrize = FALSE
) {
  stopifnot(inherits(A, "dgCMatrix"))

  sym_fro <- NA_real_
  is_num_sym <- NA
  if (check_numeric_sym) {
    sym_fro <- Matrix::norm(A - Matrix::t(A), "F")
    is_num_sym <- is.finite(sym_fro) && sym_fro <= sym_tol
  }

  A_use <- A
  if (force_symmetrize) {
    A_use <- methods::as((A + Matrix::t(A)) / 2, "dgCMatrix")
  }

  if (prefer_author) {
    tA_use <- Matrix::t(A_use)
    mode <- "author:t(A)"
  } else if (isTRUE(is_num_sym)) {
    tA_use <- A_use
    mode <- "optimized:use_A_as_tA"
  } else {
    tA_use <- Matrix::t(A_use)
    mode <- "fallback:t(A)"
  }

  list(A = A_use, tA = tA_use, sym_fro = sym_fro, is_num_sym = is_num_sym, mode = mode)
}


center_scale_scab_columns <- function(M) {
  if (!is.matrix(M)) {
    M <- as.matrix(M)
  }
  if (nrow(M) < 2) {
    stop("At least 2 genes are required to compute correlation columns.")
  }

  mu <- colMeans(M)
  M <- sweep(M, 2, mu, "-")
  s <- sqrt(colSums(M^2) / (nrow(M) - 1))
  s[!is.finite(s) | s == 0] <- 1
  sweep(M, 2, s, "/")
}


sanitize_scab_block_size <- function(block_size, n_col, fallback = 10000L) {
  if (!is.finite(n_col) || n_col < 0) {
    stop("n_col must be a non-negative finite integer.")
  }
  if (n_col == 0) {
    return(1L)
  }

  if (is.null(block_size) || !is.finite(block_size) || block_size <= 0) {
    block_size <- fallback
  }

  as.integer(max(1L, min(as.integer(block_size), as.integer(n_col))))
}


quantile_normalize_scab_block <- function(M, target) {
  M <- as.matrix(M)
  out <- preprocessCore::normalize.quantiles.use.target(M, target = target)
  out <- matrix(
    as.numeric(out),
    nrow = nrow(M),
    ncol = ncol(M),
    dimnames = dimnames(M)
  )
  out
}


compute_scab_quantile_target <- function(
  bulk_dataset,
  sc_exprs,
  common_genes,
  bulk_block_size = NULL,
  cell_block_size = 10000L
) {
  n_gene <- length(common_genes)
  if (n_gene == 0) {
    stop("common_genes must contain at least one gene.")
  }

  n_bulk <- ncol(bulk_dataset)
  n_cell <- ncol(sc_exprs)
  bulk_block_size <- sanitize_scab_block_size(
    block_size = bulk_block_size,
    n_col = n_bulk,
    fallback = min(1024L, max(1L, n_bulk))
  )
  cell_block_size <- sanitize_scab_block_size(
    block_size = cell_block_size,
    n_col = n_cell,
    fallback = 10000L
  )

  target_sum <- numeric(n_gene)
  total_columns <- 0L

  accumulate_target <- function(block) {
    block <- as.matrix(block)
    if (ncol(block) == 0) {
      return(invisible(NULL))
    }

    block_target <- preprocessCore::normalize.quantiles.determine.target(block)
    target_sum <<- target_sum + as.numeric(block_target) * ncol(block)
    total_columns <<- total_columns + ncol(block)
    invisible(NULL)
  }

  if (n_bulk > 0) {
    bulk_starts <- seq.int(1L, n_bulk, by = bulk_block_size)
    for (st in bulk_starts) {
      ed <- min(st + bulk_block_size - 1L, n_bulk)
      accumulate_target(bulk_dataset[common_genes, seq.int(st, ed), drop = FALSE])
    }
  }

  if (n_cell > 0) {
    cell_starts <- seq.int(1L, n_cell, by = cell_block_size)
    for (st in cell_starts) {
      ed <- min(st + cell_block_size - 1L, n_cell)
      accumulate_target(sc_exprs[common_genes, seq.int(st, ed), drop = FALSE])
    }
  }

  if (total_columns == 0L) {
    stop("No columns were available to determine the quantile target.")
  }

  list(
    target = target_sum / total_columns,
    total_columns = total_columns,
    bulk_block_size = bulk_block_size,
    cell_block_size = cell_block_size
  )
}


materialize_scab_quantile_bulk <- function(
  bulk_dataset,
  common_genes,
  target,
  bulk_block_size = NULL
) {
  n_gene <- length(common_genes)
  n_bulk <- ncol(bulk_dataset)
  bulk_block_size <- sanitize_scab_block_size(
    block_size = bulk_block_size,
    n_col = n_bulk,
    fallback = min(1024L, max(1L, n_bulk))
  )

  bulk_qn <- matrix(
    0,
    nrow = n_gene,
    ncol = n_bulk,
    dimnames = list(common_genes, colnames(bulk_dataset))
  )

  if (n_bulk == 0) {
    return(bulk_qn)
  }

  bulk_starts <- seq.int(1L, n_bulk, by = bulk_block_size)
  for (st in bulk_starts) {
    ed <- min(st + bulk_block_size - 1L, n_bulk)
    idx <- seq.int(st, ed)
    bulk_qn[, idx] <- quantile_normalize_scab_block(
      bulk_dataset[common_genes, idx, drop = FALSE],
      target = target
    )
  }

  bulk_qn
}


split_scab_batches <- function(values, batch_size) {
  if (length(values) == 0) {
    return(list())
  }
  split(values, ceiling(seq_along(values) / batch_size))
}


compute_scab_correlation_block <- function(
  sc_exprs,
  common_genes,
  target,
  Zb,
  n_gene,
  start_col,
  end_col,
  worker_num_threads = NULL
) {
  with_scab_num_threads(worker_num_threads, {
    idx <- seq.int(start_col, end_col)
    Cblk <- quantile_normalize_scab_block(
      sc_exprs[common_genes, idx, drop = FALSE],
      target = target
    )
    Zc <- center_scale_scab_columns(Cblk)

    Xblk <- crossprod(Zb, Zc) / (n_gene - 1)
    Xblk[is.na(Xblk)] <- 0

    list(
      start_col = start_col,
      end_col = end_col,
      Xblk = Xblk,
      fro_sq = sum(Xblk^2)
    )
  })
}


build_scab_correlation_blocks_from_quantile_target <- function(
  Expression_bulk,
  sc_exprs,
  common_genes,
  target,
  block_size = 10000L,
  x_backend = c("hdf5", "memory"),
  x_h5_path = NULL,
  x_h5_dataset = "X",
  parallel = FALSE,
  nworkers = 1L,
  worker_num_threads = NULL
) {
  x_backend <- match.arg(x_backend)

  Expression_bulk <- as.matrix(Expression_bulk)
  n_bulk <- ncol(Expression_bulk)
  n_cell <- ncol(sc_exprs)
  n_gene <- nrow(Expression_bulk)
  block_size <- sanitize_scab_block_size(
    block_size = block_size,
    n_col = n_cell,
    fallback = 10000L
  )

  Zb <- center_scale_scab_columns(Expression_bulk)
  idx_starts <- if (n_cell > 0) seq.int(1L, n_cell, by = block_size) else integer(0)
  fro_sq_total <- 0
  parallel_cfg <- resolve_scab_parallel_config(
    num_threads = worker_num_threads,
    use_parallel = parallel,
    nworkers = nworkers,
    parallel_label = "x_parallel"
  )

  if (x_backend == "memory") {
    X_out <- matrix(0, nrow = n_bulk, ncol = n_cell)
  } else {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      stop("Please install rhdf5 to use x_backend='hdf5'.")
    }
    if (is.null(x_h5_path) || !nzchar(x_h5_path)) {
      x_h5_path <- tempfile(pattern = "scAB_X_", fileext = ".h5")
    }
    if (file.exists(x_h5_path)) {
      unlink(x_h5_path)
    }

    rhdf5::h5createFile(x_h5_path)
    rhdf5::h5createDataset(
      file = x_h5_path,
      dataset = x_h5_dataset,
      dims = c(n_bulk, n_cell),
      storage.mode = "double",
      chunk = c(n_bulk, min(block_size, max(1L, n_cell)))
    )
    X_out <- NULL
  }

  process_block_result <- function(res) {
    idx <- seq.int(res$start_col, res$end_col)
    if (x_backend == "memory") {
      X_out[, idx] <<- res$Xblk
    } else {
      rhdf5::h5write(
        res$Xblk,
        file = x_h5_path,
        name = x_h5_dataset,
        index = list(seq_len(n_bulk), idx)
      )
    }
    fro_sq_total <<- fro_sq_total + res$fro_sq
    invisible(NULL)
  }

  if (parallel_cfg$parallel_enabled) {
    for (batch in split_scab_batches(idx_starts, parallel_cfg$nworkers_effective)) {
      batch_res <- parallel::mclapply(
        batch,
        function(st) {
          ed <- min(st + block_size - 1L, n_cell)
          compute_scab_correlation_block(
            sc_exprs = sc_exprs,
            common_genes = common_genes,
            target = target,
            Zb = Zb,
            n_gene = n_gene,
            start_col = st,
            end_col = ed,
            worker_num_threads = parallel_cfg$num_threads_effective
          )
        },
        mc.cores = parallel_cfg$nworkers_effective
      )
      for (res in batch_res) {
        process_block_result(res)
      }
    }
  } else {
    for (st in idx_starts) {
      ed <- min(st + block_size - 1L, n_cell)
      process_block_result(
        compute_scab_correlation_block(
          sc_exprs = sc_exprs,
          common_genes = common_genes,
          target = target,
          Zb = Zb,
          n_gene = n_gene,
          start_col = st,
          end_col = ed,
          worker_num_threads = parallel_cfg$num_threads_effective
        )
      )
    }
  }

  fro_norm <- sqrt(fro_sq_total)
  if (!is.finite(fro_norm) || fro_norm == 0) {
    stop("Encountered a non-finite Frobenius norm while building X.")
  }

  if (x_backend == "memory") {
    X_out <- X_out / fro_norm
  } else {
    for (st in idx_starts) {
      ed <- min(st + block_size - 1L, n_cell)
      idx <- seq.int(st, ed)
      Xblk <- rhdf5::h5read(
        x_h5_path,
        x_h5_dataset,
        index = list(seq_len(n_bulk), idx)
      )
      Xblk <- Xblk / fro_norm
      rhdf5::h5write(
        Xblk,
        file = x_h5_path,
        name = x_h5_dataset,
        index = list(seq_len(n_bulk), idx)
      )
    }
  }

  list(
    X = X_out,
    X_h5 = if (x_backend == "hdf5") x_h5_path else NULL,
    X_h5_dataset = if (x_backend == "hdf5") x_h5_dataset else NULL,
    X_dim = c(n_bulk, n_cell),
    fro_norm = fro_norm,
    qn_mode = "blockwise_target"
  )
}


build_scab_correlation_blocks <- function(
  Expression_bulk,
  Expression_cell,
  block_size = 10000L,
  x_backend = c("hdf5", "memory"),
  x_h5_path = NULL,
  x_h5_dataset = "X"
) {
  x_backend <- match.arg(x_backend)

  Expression_bulk <- as.matrix(Expression_bulk)
  Expression_cell <- as.matrix(Expression_cell)

  n_bulk <- ncol(Expression_bulk)
  n_cell <- ncol(Expression_cell)
  n_gene <- nrow(Expression_bulk)

  Zb <- center_scale_scab_columns(Expression_bulk)
  idx_starts <- seq.int(1L, n_cell, by = as.integer(block_size))
  fro_sq_total <- 0

  if (x_backend == "memory") {
    X_out <- matrix(0, nrow = n_bulk, ncol = n_cell)
  } else {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      stop("Please install rhdf5 to use x_backend='hdf5'.")
    }
    if (is.null(x_h5_path) || !nzchar(x_h5_path)) {
      x_h5_path <- tempfile(pattern = "scAB_X_", fileext = ".h5")
    }
    if (file.exists(x_h5_path)) {
      unlink(x_h5_path)
    }

    rhdf5::h5createFile(x_h5_path)
    rhdf5::h5createDataset(
      file = x_h5_path,
      dataset = x_h5_dataset,
      dims = c(n_bulk, n_cell),
      storage.mode = "double",
      chunk = c(n_bulk, min(as.integer(block_size), n_cell))
    )
    X_out <- NULL
  }

  for (st in idx_starts) {
    ed <- min(st + as.integer(block_size) - 1L, n_cell)
    Cblk <- Expression_cell[, seq.int(st, ed), drop = FALSE]
    Zc <- center_scale_scab_columns(Cblk)

    Xblk <- crossprod(Zb, Zc) / (n_gene - 1)
    Xblk[is.na(Xblk)] <- 0

    if (x_backend == "memory") {
      X_out[, seq.int(st, ed)] <- Xblk
    } else {
      rhdf5::h5write(
        Xblk,
        file = x_h5_path,
        name = x_h5_dataset,
        index = list(seq_len(n_bulk), seq.int(st, ed))
      )
    }

    fro_sq_total <- fro_sq_total + sum(Xblk^2)
  }

  fro_norm <- sqrt(fro_sq_total)
  if (!is.finite(fro_norm) || fro_norm == 0) {
    stop("Encountered a non-finite Frobenius norm while building X.")
  }

  if (x_backend == "memory") {
    X_out <- X_out / fro_norm
  } else {
    for (st in idx_starts) {
      ed <- min(st + as.integer(block_size) - 1L, n_cell)
      Xblk <- rhdf5::h5read(
        x_h5_path,
        x_h5_dataset,
        index = list(seq_len(n_bulk), seq.int(st, ed))
      )
      Xblk <- Xblk / fro_norm
      rhdf5::h5write(
        Xblk,
        file = x_h5_path,
        name = x_h5_dataset,
        index = list(seq_len(n_bulk), seq.int(st, ed))
      )
    }
  }

  list(
    X = X_out,
    X_h5 = if (x_backend == "hdf5") x_h5_path else NULL,
    X_h5_dataset = if (x_backend == "hdf5") x_h5_dataset else NULL,
    X_dim = c(n_bulk, n_cell),
    fro_norm = fro_norm
  )
}


materialize_scAB_X <- function(Object) {
  row_idx <- if (!is.null(Object$X_row_index)) as.integer(Object$X_row_index) else NULL
  x_scale <- if (!is.null(Object$X_scale)) as.numeric(Object$X_scale[[1]]) else NULL

  if (!is.null(Object$X)) {
    out <- as.matrix(
      Object$X[
        if (is.null(row_idx)) seq_len(nrow(Object$X)) else row_idx,
        ,
        drop = FALSE
      ]
    )
    if (!is.null(x_scale) && is.finite(x_scale) && x_scale != 1) {
      out <- out * x_scale
    }
    return(out)
  }

  if (is.null(Object$X_h5) || is.null(Object$X_h5_dataset)) {
    stop("Object does not contain in-memory X or an HDF5 X reference.")
  }
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Please install rhdf5 to materialize X from HDF5.")
  }

  x_dim <- if (!is.null(Object$X_dim)) as.integer(Object$X_dim) else NULL
  if (is.null(x_dim) || length(x_dim) != 2) {
    stop("Object$X_dim is required to materialize X from HDF5.")
  }

  out <- rhdf5::h5read(
    Object$X_h5,
    Object$X_h5_dataset,
    index = list(
      if (is.null(row_idx)) seq_len(x_dim[1]) else row_idx,
      seq_len(x_dim[2])
    )
  )
  if (!is.null(x_scale) && is.finite(x_scale) && x_scale != 1) {
    out <- out * x_scale
  }
  out
}


read_scAB_X_block <- function(Object, start_col, end_col) {
  row_idx <- if (!is.null(Object$X_row_index)) as.integer(Object$X_row_index) else NULL
  x_scale <- if (!is.null(Object$X_scale)) as.numeric(Object$X_scale[[1]]) else NULL

  if (!is.null(Object$X)) {
    out <- as.matrix(
      Object$X[
        if (is.null(row_idx)) seq_len(nrow(Object$X)) else row_idx,
        seq.int(start_col, end_col),
        drop = FALSE
      ]
    )
    if (!is.null(x_scale) && is.finite(x_scale) && x_scale != 1) {
      out <- out * x_scale
    }
    return(out)
  }

  if (is.null(Object$X_h5) || is.null(Object$X_h5_dataset)) {
    stop("Object does not contain in-memory X or an HDF5 X reference.")
  }
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Please install rhdf5 to read X blocks from HDF5.")
  }

  x_dim <- if (!is.null(Object$X_dim)) Object$X_dim else NULL
  if (is.null(x_dim) || length(x_dim) != 2) {
    stop("Object$X_dim is required to read X blocks from HDF5.")
  }

  out <- rhdf5::h5read(
    Object$X_h5,
    Object$X_h5_dataset,
    index = list(
      if (is.null(row_idx)) seq_len(x_dim[1]) else row_idx,
      seq.int(start_col, end_col)
    )
  )
  if (!is.null(x_scale) && is.finite(x_scale) && x_scale != 1) {
    out <- out * x_scale
  }
  out
}


get_scAB_X_dim <- function(Object) {
  row_idx <- if (!is.null(Object$X_row_index)) as.integer(Object$X_row_index) else NULL

  if (!is.null(Object$X)) {
    x_dim <- dim(Object$X)
    if (is.null(row_idx)) {
      return(x_dim)
    }
    return(c(length(row_idx), x_dim[2]))
  }

  if (!is.null(Object$X_dim)) {
    x_dim <- as.integer(Object$X_dim)
    if (is.null(row_idx)) {
      return(x_dim)
    }
    return(c(length(row_idx), x_dim[2]))
  }

  stop("Object does not contain X dimensions.")
}


materialize_scAB_graph_dense <- function(Object) {
  if (is.null(Object$A)) {
    stop("Object does not contain graph matrix A.")
  }

  A_dense <- as.matrix(Object$A)

  if (is.matrix(Object$D)) {
    D_diag <- diag(Object$D)
    D_dense <- Object$D
  } else {
    D_diag <- as.numeric(Object$D)
    D_dense <- diag(D_diag, nrow = length(D_diag), ncol = length(D_diag))
  }

  L_dense <- if (!is.null(Object$L)) {
    as.matrix(Object$L)
  } else {
    D_dense - A_dense
  }

  list(A = A_dense, D = D_dense, D_diag = D_diag, L = L_dense)
}


create_scAB_large <- function(
  Obejct,
  bulk_dataset,
  phenotype,
  method = c("survival", "binary"),
  assay = "RNA",
  graph_name = "RNA_snn",
  binarize_graph = TRUE,
  block_size = 10000L,
  bulk_block_size = NULL,
  num_threads = NULL,
  x_parallel = FALSE,
  x_nworkers = 1L,
  x_backend = c("hdf5", "memory"),
  x_h5_path = NULL,
  x_h5_dataset = "X",
  verify_equivalence = FALSE
) {
  method <- match.arg(method)
  x_backend <- match.arg(x_backend)
  x_cfg <- resolve_scab_parallel_config(
    num_threads = num_threads,
    use_parallel = x_parallel,
    nworkers = x_nworkers,
    parallel_label = "x_parallel"
  )

  with_scab_num_threads(x_cfg$num_threads_effective, {
    A <- get_scab_graph_sparse(Obejct = Obejct, graph_name = graph_name)
    Matrix::diag(A) <- 0
    A <- Matrix::drop0(A)
    if (binarize_graph && length(A@x) > 0) {
      A@x[A@x != 0] <- 1
    }

    deg <- Matrix::rowSums(A)
    inv_sqrt_deg <- ifelse(deg > 0, 1 / sqrt(deg), 0)

    i_idx <- A@i + 1L
    j_idx <- rep.int(seq_len(ncol(A)), diff(A@p))
    A@x <- A@x * (inv_sqrt_deg[i_idx] * inv_sqrt_deg[j_idx])
    Ahat <- A
    dhat <- ifelse(deg > 0, 1, 0)
    gm <- choose_scab_graph_mats(
      Ahat,
      prefer_author = TRUE,
      check_numeric_sym = TRUE,
      sym_tol = 1e-12,
      force_symmetrize = FALSE
    )

    bulk_dataset <- as.matrix(bulk_dataset)
    sc_exprs <- get_scab_assay_layer(Obejct = Obejct, assay = assay, layer = "data")
    common <- intersect(rownames(bulk_dataset), rownames(sc_exprs))
    if (length(common) == 0) {
      stop("No overlapping genes were found between bulk_dataset and single-cell RNA data.")
    }

    qn_res <- compute_scab_quantile_target(
      bulk_dataset = bulk_dataset,
      sc_exprs = sc_exprs,
      common_genes = common,
      bulk_block_size = bulk_block_size,
      cell_block_size = block_size
    )
    Expression_bulk <- materialize_scab_quantile_bulk(
      bulk_dataset = bulk_dataset,
      common_genes = common,
      target = qn_res$target,
      bulk_block_size = qn_res$bulk_block_size
    )

    if (method == "survival") {
      phenotype <- as.data.frame(phenotype)
      if (is.null(rownames(phenotype))) {
        stop("For method='survival', phenotype must have row names matching bulk samples.")
      }
      ss <- guanrank(phenotype[, c("time", "status")])
      S <- diag(1 - ss[rownames(phenotype), 3])
    } else {
      S <- diag(1 - phenotype)
    }

    cor_res <- build_scab_correlation_blocks_from_quantile_target(
      Expression_bulk = Expression_bulk,
      sc_exprs = sc_exprs,
      common_genes = common,
      target = qn_res$target,
      block_size = qn_res$cell_block_size,
      x_backend = x_backend,
      x_h5_path = x_h5_path,
      x_h5_dataset = x_h5_dataset,
      parallel = x_cfg$parallel_enabled,
      nworkers = x_cfg$nworkers_effective,
      worker_num_threads = x_cfg$num_threads_effective
    )

    if (verify_equivalence) {
      dataset0 <- cbind(
        bulk_dataset[common, , drop = FALSE],
        as.matrix(sc_exprs[common, , drop = FALSE])
      )
      dataset1 <- preprocessCore::normalize.quantiles(as.matrix(dataset0))
      rownames(dataset1) <- rownames(dataset0)
      colnames(dataset1) <- colnames(dataset0)
      n_bulk <- ncol(bulk_dataset)
      Expression_cell <- dataset1[, (n_bulk + 1L):ncol(dataset1), drop = FALSE]

      X_ref <- cor(Expression_bulk, Expression_cell)
      X_ref[is.na(X_ref)] <- 0
      X_ref <- X_ref / norm(X_ref, "F")
      X_cmp <- if (x_backend == "memory") cor_res$X else rhdf5::h5read(cor_res$X_h5, cor_res$X_h5_dataset)
      max_abs <- max(abs(X_ref - X_cmp))
      if (!isTRUE(max_abs <= 1e-8)) {
        stop(sprintf("Blockwise correlation check failed: max_abs_diff=%0.12f", max_abs))
      }
    }

    obj <- list(
      X = cor_res$X,
      X_h5 = cor_res$X_h5,
      X_h5_dataset = cor_res$X_h5_dataset,
      X_dim = cor_res$X_dim,
      S = S,
      A = gm$A,
      tA = gm$tA,
      D = dhat,
      L = NULL,
      phenotype = phenotype,
      method = method,
      graph_info = list(
        graph_name = graph_name,
        sym_fro = gm$sym_fro,
        is_num_sym = gm$is_num_sym,
        transpose_mode = gm$mode,
        zero_degree = sum(deg == 0),
        min_degree = unname(min(deg)),
        max_degree = unname(max(deg)),
        mean_degree = unname(mean(deg))
      ),
      common_genes = common,
      x_backend = x_backend,
      block_size = as.integer(qn_res$cell_block_size),
      x_fro_norm = cor_res$fro_norm,
      qn_info = list(
        mode = cor_res$qn_mode,
        target_length = length(qn_res$target),
        total_columns = qn_res$total_columns,
        bulk_block_size = qn_res$bulk_block_size,
        cell_block_size = qn_res$cell_block_size
      ),
      parallel_info = list(
        num_threads_requested = x_cfg$num_threads_requested,
        num_threads_effective = x_cfg$num_threads_effective,
        x_parallel_requested = isTRUE(x_parallel),
        x_parallel_effective = x_cfg$parallel_enabled,
        x_nworkers_requested = x_cfg$nworkers_requested,
        x_nworkers_effective = x_cfg$nworkers_effective
      )
    )
    class(obj) <- "scAB_data"
    obj
  })
}
