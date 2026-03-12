resolve_scab_seed <- function(method, seed = NULL) {
  if (!is.null(seed)) {
    return(as.integer(seed[[1]]))
  }

  if (!is.null(method) && length(method) == 1 && nzchar(method)) {
    return(if (method == "survival") 7L else 5L)
  }

  NULL
}


initialize_scab_factors <- function(nr, nc, K, method = "", seed = NULL, W0 = NULL, H0 = NULL) {
  if (!is.null(W0) && !is.null(H0)) {
    W <- matrix(as.numeric(W0), nrow = nr, ncol = K)
    H <- matrix(as.numeric(H0), nrow = K, ncol = nc)
    if (!identical(dim(W), c(nr, K)) || !identical(dim(H), c(K, nc))) {
      stop("Provided W0/H0 do not match the requested dimensions.")
    }
    return(list(W = W, H = H, seed_effective = NA_integer_, init_source = "provided"))
  }

  seed_to_use <- resolve_scab_seed(method = method, seed = seed)
  if (!is.null(seed_to_use) && !is.na(seed_to_use)) {
    set.seed(seed_to_use)
  }

  list(
    W = matrix(runif(nr * K), nrow = nr, ncol = K),
    H = matrix(runif(K * nc), nrow = K, ncol = nc),
    seed_effective = seed_to_use,
    init_source = "runif"
  )
}


get_scab_col_block_starts <- function(nc, block_size) {
  block_size <- sanitize_scab_block_size(
    block_size = block_size,
    n_col = nc,
    fallback = 10000L
  )
  seq.int(1L, nc, by = block_size)
}


read_scab_fit_X_block <- function(X, Object, start_col, end_col) {
  if (!is.null(X)) {
    return(as.matrix(X[, seq.int(start_col, end_col), drop = FALSE]))
  }

  read_scAB_X_block(Object, start_col, end_col)
}


multiply_scab_H_tA_block <- function(H, tA, start_col, end_col) {
  idx <- seq.int(start_col, end_col)
  as.matrix(H %*% tA[, idx, drop = FALSE])
}


accumulate_scab_XHt <- function(X, H, Object, nr, K, idx_starts, block_size, nc) {
  out <- matrix(0, nrow = nr, ncol = K)
  for (st in idx_starts) {
    ed <- min(st + block_size - 1L, nc)
    idx <- seq.int(st, ed)
    Xblk <- read_scab_fit_X_block(X, Object, st, ed)
    out <- out + Xblk %*% t(H[, idx, drop = FALSE])
  }
  out
}


accumulate_scab_reconstruction_loss <- function(X, W, H, Object, idx_starts, block_size, nc) {
  total <- 0
  for (st in idx_starts) {
    ed <- min(st + block_size - 1L, nc)
    idx <- seq.int(st, ed)
    Xblk <- read_scab_fit_X_block(X, Object, st, ed)
    Hblk <- H[, idx, drop = FALSE]
    total <- total + sum((Xblk - W %*% Hblk)^2)
  }
  total
}


accumulate_scab_graph_quadratic <- function(H, tA, idx_starts, block_size, nc) {
  total <- 0
  for (st in idx_starts) {
    ed <- min(st + block_size - 1L, nc)
    idx <- seq.int(st, ed)
    HAt_blk <- multiply_scab_H_tA_block(H, tA, st, ed)
    total <- total + sum(H[, idx, drop = FALSE] * HAt_blk)
  }
  total
}


scAB_large <- function(
  Object,
  K,
  alpha = 0.005,
  alpha_2 = 0.005,
  maxiter = 2000,
  convergence_threshold = 1e-5,
  eps = 0,
  seed = NULL,
  guard_nonfinite = FALSE,
  verbose = FALSE,
  check_every = 1L,
  materialize_hdf5 = TRUE,
  x_block_size = NULL,
  num_threads = NULL
) {
  fit_cfg <- resolve_scab_parallel_config(
    num_threads = num_threads,
    use_parallel = FALSE,
    nworkers = 1L,
    parallel_label = "scAB_large"
  )

  with_scab_num_threads(fit_cfg$num_threads_effective, {
    if (is.null(Object$X) && !materialize_hdf5) {
      if (is.null(Object$X_h5) || is.null(Object$X_h5_dataset)) {
        stop("Object$X is NULL and no HDF5 X reference is available for blockwise fitting.")
      }
    }

    X <- if (is.null(Object$X) && !materialize_hdf5) NULL else materialize_scAB_X(Object)
    x_dim <- get_scAB_X_dim(Object)
    nr <- x_dim[1]
    nc <- x_dim[2]
    stream_X <- is.null(X)

    x_block_size <- sanitize_scab_block_size(
      block_size = x_block_size,
      n_col = nc,
      fallback = if (!is.null(Object$block_size)) as.integer(Object$block_size) else 10000L
    )
    idx_starts <- get_scab_col_block_starts(nc = nc, block_size = x_block_size)

    A <- Object$A
    if (is.null(A)) {
      stop("Object does not contain graph matrix A.")
    }
    tA <- Object$tA
    if (is.null(tA)) {
      tA <- Matrix::t(A)
    }

    D_diag <- if (is.matrix(Object$D)) diag(Object$D) else as.numeric(Object$D)
    if (length(D_diag) != nc) {
      stop("Graph degree vector length does not match the number of cells.")
    }

    S <- Object$S
    SS <- S %*% S

    init <- initialize_scab_factors(
      nr = nr,
      nc = nc,
      K = K,
      method = Object$method,
      seed = seed,
      W0 = Object$W0,
      H0 = Object$H0
    )
    W <- init$W
    H <- init$H

    right_multiply_diag <- function(Hblk, d) sweep(Hblk, 2, d, "*")

    loss_func <- function(X, W, H, use_L = !is.null(Object$L)) {
      loss1 <- accumulate_scab_reconstruction_loss(
        X = X,
        W = W,
        H = H,
        Object = Object,
        idx_starts = idx_starts,
        block_size = x_block_size,
        nc = nc
      )
      loss2 <- alpha * (norm(S %*% W, "F")^2)
      if (use_L) {
        loss3 <- alpha_2 * sum(diag(H %*% Object$L %*% t(H)))
      } else {
        HDH <- sum(H * right_multiply_diag(H, D_diag))
        HAH <- accumulate_scab_graph_quadratic(
          H = H,
          tA = tA,
          idx_starts = idx_starts,
          block_size = x_block_size,
          nc = nc
        )
        loss3 <- alpha_2 * (HDH - HAH)
      }
      as.numeric(loss1 + loss2 + loss3)
    }

    old_eucl <- loss_func(X, W, H)

    for (iter in 1:maxiter) {
      WtW <- t(W) %*% W
      H_old <- H
      for (st in idx_starts) {
        ed <- min(st + x_block_size - 1L, nc)
        idx <- seq.int(st, ed)
        Xblk <- read_scab_fit_X_block(X, Object, st, ed)
        H_old_blk <- H_old[, idx, drop = FALSE]
        HAt_blk <- multiply_scab_H_tA_block(H_old, tA, st, ed)
        WtXblk <- t(W) %*% Xblk
        HDt_blk <- right_multiply_diag(H_old_blk, D_diag[idx])

        H_den_blk <- (WtW %*% H_old_blk + alpha_2 * HDt_blk)
        if (eps != 0) {
          H_den_blk <- H_den_blk + eps
        }
        H[, idx] <- as.matrix(
          H_old_blk * (WtXblk + alpha_2 * HAt_blk) / H_den_blk
        )
      }

      HHt <- H %*% t(H)
      XHt <- accumulate_scab_XHt(
        X = X,
        H = H,
        Object = Object,
        nr = nr,
        K = K,
        idx_starts = idx_starts,
        block_size = x_block_size,
        nc = nc
      )
      Pena <- SS %*% W

      W_den <- (W %*% HHt + alpha * Pena)
      if (eps != 0) {
        W_den <- W_den + eps
      }
      W <- W * XHt / W_den

      if (guard_nonfinite) {
        H[!is.finite(H)] <- 0
        H[H < 0] <- 0
        W[!is.finite(W)] <- 0
        W[W < 0] <- 0
      }

      if (iter %% check_every == 0) {
        eucl_dist <- loss_func(X, W, H)
        d_eucl <- abs(eucl_dist - old_eucl)
        if (verbose && iter %% 10 == 0) {
          message("iter=", iter, " loss=", eucl_dist, " d=", d_eucl)
        }

        if (!is.finite(d_eucl) || is.na(d_eucl)) {
          break
        }
        if (d_eucl < convergence_threshold) {
          old_eucl <- eucl_dist
          break
        }
        old_eucl <- eucl_dist
      }

      iter <- iter + 1L
    }

    list(
      W = as.matrix(W),
      H = as.matrix(H),
      iter = iter,
      loss = loss_func(X, W, H),
      method = Object$method,
      seed_effective = init$seed_effective,
      init_source = init$init_source,
      num_threads_requested = fit_cfg$num_threads_requested,
      num_threads_effective = fit_cfg$num_threads_effective
    )
  })
}
