###  Non-negative Matrix Factorization
#' Classical non-negative matrix factorization algorithm
#'
#' @param X  a non-negative matrix
#' @param K  the rank of matrix factorization
#' @param maxiter the maximum number of iterations
#'
#' @return a list with the submatrix and loss value
#' @export
#'
#' @examples
NMF <- function(X, K,maxiter=2000){
  eps=2.2204e-256
  nr = dim(X)[1];
  nc = dim(X)[2];
  W = matrix(runif(nr*K),nrow=nr, ncol=K);
  H = matrix(runif(K*nc),nrow=K, ncol=nc);
  loss_func=function(X,W,H){
    loss<-norm(X-W%*%H,"F")^2
    return(as.numeric(loss))
  }
  
  for (iter in 1:maxiter){
    H = H*(t(W)%*%X)/ ((t(W)%*%W)%*%H )
    W=W*(X%*%t(H))/(W %*% H %*% t(H)  )
    {
      if(iter!=1){
        eucl_dist = loss_func(X,W,H)
        d_eucl=abs(eucl_dist-old_eucl)
        if(d_eucl<10^(-5)) {break;}
        old_eucl=eucl_dist
      }
      else{ old_eucl =loss_func(X,W,H)
      }
    }
    iter=iter+1
  }
  return(list(W=W,H=H,iter=iter,loss=loss_func(X,W,H) ))
}

###  Non-negative Matrix Factorization with phenotype and cell-cell similarity regularization.
#' Non-negative Matrix Factorization with phenotype and cell-cell similarity regularization, for identifing phenotype-associated cell states at different resolutions.
#'
#' @param Object  a scAB_data object
#' @param K  the rank of matrix factorization
#' @param maxiter the maximum number of iterations
#' @param alpha Coefficient of phenotype regularization
#' @param alpha_2 Coefficient of cell-cell similarity regularization
#'
#' @return a list with the submatrix and loss value
#' @export
#'
#' @examples
scAB <- function(Object, K,alpha=0.005,alpha_2=0.005,maxiter=2000,seed=NULL){
  seed_to_use <- NULL
  if (is.null(seed)) {
    if (!is.null(Object$method) && length(Object$method) == 1 && nzchar(Object$method)) {
      seed_to_use <- if (Object$method == "survival") 7L else 5L
    }
  } else {
    seed_to_use <- as.integer(seed[[1]])
  }
  if (!is.null(seed_to_use) && !is.na(seed_to_use)) set.seed(seed_to_use)
  X <- Object$X
  A <- Object$A
  L <- Object$L
  D <- Object$D
  S <- Object$S
  eps=2.2204e-256
  nr = dim(X)[1];
  nc = dim(X)[2];
  W = matrix(runif(nr*K),nrow=nr, ncol=K);
  H = matrix(runif(K*nc),nrow=K, ncol=nc);
  loss_func=function(X,W,H,S,L,alpha,alpha_2){
    loss<-norm(X-W%*%H,"F")^2 +
      alpha*( norm(S%*%W,"F")^2) +
      alpha_2*sum(diag(H %*% L %*% t(H)))
    return(as.numeric(loss))
  }
  tD <- t(D)
  tA <- t(A)
  for (iter in 1:maxiter){
    W_old <- W
    H <- H*(t(W)%*%X+alpha_2*H%*%tA)/ ((t(W)%*%W)%*%H + alpha_2*H%*%tD)
    Pena <- (S%*%S)%*%W
    W <- W*(X%*%t(H))/(W %*% H %*% t(H)  + alpha*Pena )
    {
      if(iter!=1){
        eucl_dist <- loss_func(X,W,H,S,L,alpha,alpha_2)
        d_eucl <- abs(eucl_dist-old_eucl)
        if(d_eucl<10^(-5)) {break;}
        old_eucl <- eucl_dist
      }
      else{ old_eucl <- loss_func(X,W,H,S,L,alpha,alpha_2)
      }
    }
    iter <- iter+1
  }
  return(list(W=W,H=H,iter=iter,loss=loss_func(X,W,H,S,L,alpha,alpha_2),method=Object$method ))
}


###  Create subsets of cross-validation
#'
#' @param k  k-fold cross validation
#' @param datasize  the size of samples
#' @param seed random seed
#'
#' @return a list with subsets of cross-validation
#' @export
#'
#' @examples
CVgroup <- function(k,datasize,seed = 0){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]
  temp <- sample(n,datasize)
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])
  return(cvlist)
}


###  Selection of parameter K
#'
#' @param Object a scAB_data object
#' @param k_max  the maximum value of the rank in the matrix factorization
#' @param repeat_times  the number of repetitions
#' @param seed random seed
#' @param verbose Print output
#' @param use_large Whether to force the large-object path. `NULL` means
#'   auto-detect based on the object structure.
#' @param materialize_hdf5 Whether to materialize HDF5-backed `X` before fitting
#'   when `use_large=TRUE`.
#' @param x_block_size Optional block size for HDF5-backed `X` during
#'   `use_large=TRUE`.
#' @param num_threads Optional BLAS/OMP thread count for `use_large=TRUE`.
#' @param k_parallel Whether to parallelize repeated fits within each candidate
#'   `K` when `use_large=TRUE`.
#' @param k_nworkers Worker count for `k_parallel=TRUE`.
#' @param return_details Whether to return a detail list instead of only the
#'   selected `K`.
#'
#' @return an appropriate value of K
#' @export
#'
#' @examples
select_K_large <- function(
  Object,
  K_max = 20,
  repeat_times = 10,
  maxiter = 2000,
  seed = 0,
  verbose = FALSE,
  materialize_hdf5 = TRUE,
  x_block_size = NULL,
  num_threads = NULL,
  k_parallel = FALSE,
  k_nworkers = 1L,
  return_details = FALSE
) {
  if (K_max < 2) {
    stop("`K_max` must be at least 2.")
  }

  k_cfg <- resolve_scab_parallel_config(
    num_threads = num_threads,
    use_parallel = k_parallel,
    nworkers = k_nworkers,
    parallel_label = "k_parallel"
  )
  K_all <- 2:K_max
  dist_K <- matrix(NA_real_, nrow = length(K_all), ncol = repeat_times)
  eii <- rep(NA_real_, length(K_all))
  completed_K <- integer(0)

  set.seed(seed)
  seed_grid <- matrix(
    sample.int(.Machine$integer.max, length(K_all) * repeat_times, replace = FALSE),
    nrow = length(K_all),
    ncol = repeat_times
  )

  Ki_last <- K_all[[1]]
  for (ki_idx in seq_along(K_all)) {
    Ki <- K_all[[ki_idx]]
    Ki_last <- Ki

    run_one_repeat <- function(Kj) {
      res_ij <- scAB_large(
        Object = Object,
        K = Ki,
        alpha = 0,
        alpha_2 = 0,
        maxiter = maxiter,
        seed = seed_grid[ki_idx, Kj],
        materialize_hdf5 = materialize_hdf5,
        x_block_size = x_block_size,
        num_threads = k_cfg$num_threads_effective
      )
      as.numeric(res_ij$loss)
    }

    if (k_cfg$parallel_enabled) {
      losses_ki <- unlist(
        parallel::mclapply(
          seq_len(repeat_times),
          run_one_repeat,
          mc.cores = k_cfg$nworkers_effective
        ),
        use.names = FALSE
      )
    } else {
      losses_ki <- vapply(seq_len(repeat_times), run_one_repeat, numeric(1))
    }

    dist_K[ki_idx, ] <- losses_ki
    completed_K <- c(completed_K, Ki)

    if (verbose) {
      print(paste0("loss of ", Ki, ": ", mean(dist_K[ki_idx, ], na.rm = TRUE)))
    }

    if (Ki == 2) {
      next
    }

    mean_dist <- rowMeans(dist_K[seq_len(ki_idx), , drop = FALSE], na.rm = TRUE)
    loss_drop <- mean_dist[ki_idx - 1L] - mean_dist[ki_idx]
    denom <- mean_dist[1L] - mean_dist[ki_idx]
    eii[ki_idx] <- if (isTRUE(denom > 0)) loss_drop / denom else NA_real_

    if (!isTRUE(loss_drop > 0)) {
      break
    }
    if (is.finite(eii[ki_idx]) && eii[ki_idx] < 0.05) {
      break
    }
  }

  selected_K <- max(2L, as.integer(Ki_last - 1L))
  detail <- list(
    selected_K = selected_K,
    K_all = K_all,
    completed_K = completed_K,
    dist_K = dist_K,
    mean_loss = rowMeans(dist_K, na.rm = TRUE),
    eii = eii,
    parallel_info = list(
      k_parallel_requested = isTRUE(k_parallel),
      k_parallel_effective = k_cfg$parallel_enabled,
      k_nworkers_requested = k_cfg$nworkers_requested,
      k_nworkers_effective = k_cfg$nworkers_effective,
      num_threads_requested = k_cfg$num_threads_requested,
      num_threads_effective = k_cfg$num_threads_effective
    )
  )

  if (isTRUE(return_details)) {
    return(detail)
  }
  selected_K
}


select_K <- function(
  Object,
  K_max = 20,
  repeat_times = 10,
  maxiter = 2000,
  seed = 0,
  verbose = FALSE,
  use_large = NULL,
  materialize_hdf5 = TRUE,
  x_block_size = NULL,
  num_threads = NULL,
  k_parallel = FALSE,
  k_nworkers = 1L,
  return_details = FALSE
) {
  if (K_max < 2) {
    stop("`K_max` must be at least 2.")
  }

  auto_large <- is.null(use_large) && (
    is.null(Object$X) ||
      !is.null(Object$X_h5) ||
      !is.null(Object$x_backend) ||
      !is.null(Object$block_size)
  )

  if (isTRUE(use_large) || isTRUE(auto_large)) {
    return(select_K_large(
      Object = Object,
      K_max = K_max,
      repeat_times = repeat_times,
      maxiter = maxiter,
      seed = seed,
      verbose = verbose,
      materialize_hdf5 = materialize_hdf5,
      x_block_size = x_block_size,
      num_threads = num_threads,
      k_parallel = k_parallel,
      k_nworkers = k_nworkers,
      return_details = return_details
    ))
  }

  if (is.null(Object$X)) {
    stop("Object$X is NULL; use `select_K_large()` or `select_K(..., use_large=TRUE)` for large objects.")
  }

  X <- Object$X
  set.seed(seed)
  K_all <- 2:K_max
  dist_K <- matrix(NA_real_, nrow = length(K_all), ncol = repeat_times)
  eii <- rep(NA_real_, length(K_all))
  Ki_last <- K_all[[1]]

  for (ki_idx in seq_along(K_all)) {
    Ki <- K_all[[ki_idx]]
    Ki_last <- Ki
    for (Kj in seq_len(repeat_times)) {
      res_ij <- NMF(X = X, K = Ki, maxiter = maxiter)
      dist_K[ki_idx, Kj] <- norm(X - res_ij$W %*% res_ij$H, "F")^2
      if (Kj == repeat_times && verbose) {
        print(paste0("loss of ", Ki, ": ", mean(dist_K[ki_idx, ], na.rm = TRUE)))
      }
    }
    if (Ki == 2) {
      next
    }

    mean_dist <- rowMeans(dist_K[seq_len(ki_idx), , drop = FALSE], na.rm = TRUE)
    loss_drop <- mean_dist[ki_idx - 1L] - mean_dist[ki_idx]
    denom <- mean_dist[1L] - mean_dist[ki_idx]
    eii[ki_idx] <- if (isTRUE(denom > 0)) loss_drop / denom else NA_real_

    if (!isTRUE(loss_drop > 0)) {
      break
    }
    if (is.finite(eii[ki_idx]) && eii[ki_idx] < 0.05) {
      break
    }
  }

  selected_K <- max(2L, as.integer(Ki_last - 1L))
  if (isTRUE(return_details)) {
    return(list(
      selected_K = selected_K,
      K_all = K_all,
      completed_K = K_all[seq_len(ki_idx)],
      dist_K = dist_K,
      mean_loss = rowMeans(dist_K, na.rm = TRUE),
      eii = eii,
      parallel_info = NULL
    ))
  }
  selected_K
}


###  Selection of parameter alpha and alpha_2
#'
#' @param Object a scAB_data object
#' @param K the rank of matrix factorization
#' @param cross_k k-fold cross validation
#' @param seed random seed for the fold split
#' @param alpha1_list candidate values for `alpha`
#' @param alpha2_list candidate values for `alpha_2`
#' @param maxiter maximum iteration count for each fold fit
#' @param model_seed optional seed passed to the model fit inside each fold
#' @param verbose whether to print CV summaries
#' @param use_large whether to force the large-object path. `NULL` means
#'   auto-detect based on the object structure
#' @param materialize_hdf5 whether to materialize HDF5-backed `X`
#' @param x_block_size optional HDF5 column-block size
#' @param num_threads optional BLAS/OMP thread count
#' @param cv_parallel whether to parallelize CV tasks for large objects
#' @param cv_nworkers worker count for `cv_parallel=TRUE`
#' @param return_details whether to return fold-level detail outputs
#'
#' @return a list with alpha, alpha_2 and cross-validation error
#' @importFrom MASS ginv
#' @import survival
#' @export
#'
#' @examples
coerce_scab_cv_phenotype <- function(Object) {
  if (identical(Object$method, "survival")) {
    pheno <- as.data.frame(Object$phenotype)
    if (!all(c("time", "status") %in% colnames(pheno))) {
      stop("Survival phenotype must contain `time` and `status` columns.")
    }
    if (is.null(rownames(pheno))) {
      rownames(pheno) <- seq_len(nrow(pheno))
    }
    return(pheno[, c("time", "status"), drop = FALSE])
  }

  status <- ifelse(Object$phenotype, 1, 0)
  time <- ifelse(Object$phenotype, 1, 100)
  row_ids <- if (!is.null(Object$X)) {
    rownames(Object$X)
  } else if (!is.null(names(Object$phenotype))) {
    names(Object$phenotype)
  } else {
    as.character(seq_along(status))
  }

  data.frame(time = time, status = status, row.names = row_ids)
}


make_scab_cv_S <- function(pheno_sub) {
  ss <- guanrank(pheno_sub[, c("time", "status"), drop = FALSE])
  diag(1 - ss[rownames(pheno_sub), 3], nrow = nrow(pheno_sub), ncol = nrow(pheno_sub))
}


build_scab_cv_result_matrix <- function(cv_summary, alpha1_list, alpha2_list) {
  out <- matrix(
    NA_real_,
    nrow = length(alpha1_list),
    ncol = length(alpha2_list),
    dimnames = list(as.character(alpha1_list), as.character(alpha2_list))
  )

  for (i in seq_len(nrow(cv_summary))) {
    out[as.character(cv_summary$alpha1[[i]]), as.character(cv_summary$alpha2[[i]])] <- cv_summary$cindex[[i]]
  }

  out
}


compose_scab_row_index <- function(Object, row_index) {
  row_index <- as.integer(row_index)
  if (is.null(Object$X_row_index)) {
    return(row_index)
  }
  as.integer(Object$X_row_index)[row_index]
}


make_scab_row_view_object <- function(Object, row_index, phenotype = NULL, S = NULL, x_scale = NULL) {
  out <- Object
  out$X_row_index <- compose_scab_row_index(Object, row_index)

  parent_scale <- if (!is.null(Object$X_scale)) as.numeric(Object$X_scale[[1]]) else 1
  scale_value <- if (is.null(x_scale)) parent_scale else parent_scale * as.numeric(x_scale[[1]])
  if (!is.finite(scale_value) || scale_value <= 0) {
    stop("`x_scale` must be a single positive finite number.")
  }
  if (abs(scale_value - 1) <= .Machine$double.eps^0.5) {
    out$X_scale <- NULL
  } else {
    out$X_scale <- scale_value
  }

  if (!is.null(phenotype)) {
    out$phenotype <- phenotype
  }
  if (!is.null(S)) {
    out$S <- S
  }

  class(out) <- class(Object)
  out
}


compute_scab_object_fro_norm <- function(Object, x_block_size = NULL) {
  x_dim <- get_scAB_X_dim(Object)
  nc <- x_dim[2]
  block_size <- sanitize_scab_block_size(
    block_size = x_block_size,
    n_col = nc,
    fallback = if (!is.null(Object$block_size)) Object$block_size else 10000L
  )
  idx_starts <- get_scab_col_block_starts(nc, block_size)
  fro_sq_total <- 0

  for (st in idx_starts) {
    ed <- min(st + block_size - 1L, nc)
    Xblk <- read_scab_fit_X_block(NULL, Object, st, ed)
    fro_sq_total <- fro_sq_total + sum(Xblk * Xblk)
  }

  sqrt(fro_sq_total)
}


predict_scab_test_W_large <- function(Object, H, x_block_size = NULL) {
  x_dim <- get_scAB_X_dim(Object)
  nr <- x_dim[1]
  nc <- x_dim[2]
  block_size <- sanitize_scab_block_size(
    block_size = x_block_size,
    n_col = nc,
    fallback = if (!is.null(Object$block_size)) Object$block_size else 10000L
  )
  idx_starts <- get_scab_col_block_starts(nc, block_size)
  ginvH <- MASS::ginv(H)
  out <- matrix(0, nrow = nr, ncol = nrow(H))

  for (st in idx_starts) {
    ed <- min(st + block_size - 1L, nc)
    idx <- seq.int(st, ed)
    Xblk <- read_scab_fit_X_block(NULL, Object, st, ed)
    out <- out + Xblk %*% ginvH[idx, , drop = FALSE]
  }

  out
}


run_scab_cv_grid_large <- function(
  Object,
  K,
  alpha1_list,
  alpha2_list,
  cross_k = 5L,
  cv_seed = 0L,
  model_seed = NULL,
  maxiter = 2000L,
  verbose = FALSE,
  materialize_hdf5 = FALSE,
  x_block_size = NULL,
  num_threads = NULL,
  cv_parallel = FALSE,
  cv_nworkers = 1L
) {
  validate_scab_parallel_combination(
    cv_parallel = cv_parallel,
    x_parallel = FALSE,
    context = "select_alpha_large"
  )
  cv_cfg <- resolve_scab_parallel_config(
    num_threads = num_threads,
    use_parallel = cv_parallel,
    nworkers = cv_nworkers,
    parallel_label = "cv_parallel"
  )

  train_phenotype <- coerce_scab_cv_phenotype(Object)
  datasize <- get_scAB_X_dim(Object)[1]
  if (nrow(train_phenotype) != datasize) {
    stop("`nrow(phenotype)` must match the row dimension of `X`.")
  }

  cvlist <- CVgroup(k = cross_k, datasize = datasize, seed = cv_seed)
  row_ids <- seq_len(datasize)

  fold_data_list <- lapply(seq_len(cross_k), function(fold_id) {
    test_id <- cvlist[[fold_id]]
    train_id <- setdiff(row_ids, test_id)

    ph_train <- train_phenotype[train_id, , drop = FALSE]
    ph_test <- train_phenotype[test_id, , drop = FALSE]
    S_train <- make_scab_cv_S(ph_train)

    obj_train_raw <- make_scab_row_view_object(
      Object = Object,
      row_index = train_id,
      phenotype = ph_train,
      S = S_train
    )
    fn <- compute_scab_object_fro_norm(obj_train_raw, x_block_size = x_block_size)
    if (!is.finite(fn) || fn <= 0) {
      fn <- 1
    }

    list(
      train = make_scab_row_view_object(
        Object = Object,
        row_index = train_id,
        phenotype = ph_train,
        S = S_train,
        x_scale = 1 / fn
      ),
      test = make_scab_row_view_object(
        Object = Object,
        row_index = test_id,
        phenotype = ph_test,
        S = NULL,
        x_scale = NULL
      ),
      ph_train = ph_train,
      ph_test = ph_test
    )
  })

  tasks <- expand.grid(
    alpha1 = alpha1_list,
    alpha2 = alpha2_list,
    fold = seq_len(cross_k),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  if (isTRUE(verbose)) {
    message(
      sprintf(
        "select_alpha_large_start tasks=%d cross_k=%d alpha1_n=%d alpha2_n=%d cv_parallel=%s cv_workers=%d num_threads_effective=%d",
        nrow(tasks),
        cross_k,
        length(alpha1_list),
        length(alpha2_list),
        cv_cfg$parallel_enabled,
        cv_cfg$nworkers_effective,
        cv_cfg$num_threads_effective
      )
    )
  }

  run_one_task <- function(i) {
    task <- tasks[i, , drop = FALSE]
    fd <- fold_data_list[[task$fold]]
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "select_alpha_large_task_start alpha1=%s alpha2=%s fold=%d",
          format(task$alpha1, trim = TRUE, scientific = FALSE),
          format(task$alpha2, trim = TRUE, scientific = FALSE),
          task$fold
        )
      )
    }
    fit <- scAB_large(
      Object = fd$train,
      K = K,
      alpha = task$alpha1,
      alpha_2 = task$alpha2,
      maxiter = maxiter,
      seed = model_seed,
      materialize_hdf5 = materialize_hdf5,
      x_block_size = x_block_size,
      num_threads = cv_cfg$num_threads_effective
    )

    new_W <- predict_scab_test_W_large(
      Object = fd$test,
      H = fit$H,
      x_block_size = x_block_size
    )

    df_train <- data.frame(time = fd$ph_train$time, status = fd$ph_train$status, fit$W)
    df_test <- data.frame(time = fd$ph_test$time, status = fd$ph_test$status, new_W)

    cox_fit <- survival::coxph(survival::Surv(time, status) ~ ., data = df_train)
    score_test <- as.numeric(stats::predict(cox_fit, newdata = df_test, type = "lp"))
    cind <- as.numeric(survival::concordance(survival::Surv(df_test$time, df_test$status) ~ score_test)$concordance)

    out <- data.frame(
      alpha1 = task$alpha1,
      alpha2 = task$alpha2,
      fold = task$fold,
      cindex = cind,
      iter = fit$iter,
      loss = fit$loss,
      stringsAsFactors = FALSE
    )
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "select_alpha_large_task_end alpha1=%s alpha2=%s fold=%d cindex=%.6f iter=%d loss=%.10f",
          format(task$alpha1, trim = TRUE, scientific = FALSE),
          format(task$alpha2, trim = TRUE, scientific = FALSE),
          task$fold,
          cind,
          fit$iter,
          fit$loss
        )
      )
    }
    out
  }

  runner <- if (cv_cfg$parallel_enabled) {
    function(ids) parallel::mclapply(ids, run_one_task, mc.cores = cv_cfg$nworkers_effective)
  } else {
    function(ids) lapply(ids, run_one_task)
  }

  timing <- system.time({
    res_list <- runner(seq_len(nrow(tasks)))
  })

  cv_res <- do.call(rbind, res_list)
  cv_res <- cv_res[order(cv_res$alpha1, cv_res$alpha2, cv_res$fold), , drop = FALSE]
  rownames(cv_res) <- NULL

  cv_summary <- aggregate(cindex ~ alpha1 + alpha2, data = cv_res, FUN = mean)
  cv_summary <- cv_summary[order(-cv_summary$cindex, cv_summary$alpha1, cv_summary$alpha2), , drop = FALSE]
  rownames(cv_summary) <- NULL

  if (verbose) {
    print(cv_summary)
  }

  list(
    cv_res = cv_res,
    cv_summary = cv_summary,
    best = cv_summary[1, , drop = FALSE],
    result_cv = build_scab_cv_result_matrix(cv_summary, alpha1_list, alpha2_list),
    timing = timing,
    parallel_info = cv_cfg
  )
}


select_alpha_large <- function(
  Object,
  K,
  cross_k = 5,
  seed = 0,
  alpha1_list = c(0.01, 0.005, 0.001),
  alpha2_list = c(0.01, 0.005, 0.001),
  maxiter = 2000L,
  model_seed = NULL,
  verbose = FALSE,
  materialize_hdf5 = FALSE,
  x_block_size = NULL,
  num_threads = NULL,
  cv_parallel = FALSE,
  cv_nworkers = 1L,
  return_details = FALSE
) {
  cv_out <- run_scab_cv_grid_large(
    Object = Object,
    K = K,
    alpha1_list = alpha1_list,
    alpha2_list = alpha2_list,
    cross_k = cross_k,
    cv_seed = seed,
    model_seed = model_seed,
    maxiter = maxiter,
    verbose = verbose,
    materialize_hdf5 = materialize_hdf5,
    x_block_size = x_block_size,
    num_threads = num_threads,
    cv_parallel = cv_parallel,
    cv_nworkers = cv_nworkers
  )

  out <- list(
    alpha_1 = cv_out$best$alpha1[[1]],
    alpha_2 = cv_out$best$alpha2[[1]],
    result_cv = cv_out$result_cv
  )
  if (isTRUE(return_details)) {
    out <- c(out, cv_out[c("cv_res", "cv_summary", "best", "timing", "parallel_info")])
  }
  out
}


select_alpha <- function(
  Object,
  K,
  cross_k = 5,
  seed = 0,
  alpha1_list = c(0.01, 0.005, 0.001),
  alpha2_list = c(0.01, 0.005, 0.001),
  maxiter = 2000L,
  model_seed = NULL,
  verbose = FALSE,
  use_large = NULL,
  materialize_hdf5 = FALSE,
  x_block_size = NULL,
  num_threads = NULL,
  cv_parallel = FALSE,
  cv_nworkers = 1L,
  return_details = FALSE
) {
  auto_large <- is.null(use_large) && (
    is.null(Object$X) ||
      !is.null(Object$X_h5) ||
      !is.null(Object$x_backend) ||
      !is.null(Object$block_size)
  )

  if (isTRUE(use_large) || isTRUE(auto_large)) {
    return(select_alpha_large(
      Object = Object,
      K = K,
      cross_k = cross_k,
      seed = seed,
      alpha1_list = alpha1_list,
      alpha2_list = alpha2_list,
      maxiter = maxiter,
      model_seed = model_seed,
      verbose = verbose,
      materialize_hdf5 = materialize_hdf5,
      x_block_size = x_block_size,
      num_threads = num_threads,
      cv_parallel = cv_parallel,
      cv_nworkers = cv_nworkers,
      return_details = return_details
    ))
  }

  if (is.null(Object$X)) {
    stop("Object$X is NULL; use `select_alpha_large()` or `select_alpha(..., use_large=TRUE)` for large objects.")
  }

  train_phenotype <- coerce_scab_cv_phenotype(Object)
  train_data <- Object$X
  A_cv <- Object$A
  L_cv <- Object$L
  D_cv <- Object$D
  datasize <- nrow(train_data)
  cvlist <- CVgroup(k = cross_k, datasize = datasize, seed = seed)

  result_cv <- matrix(
    NA_real_,
    nrow = length(alpha1_list),
    ncol = length(alpha2_list),
    dimnames = list(as.character(alpha1_list), as.character(alpha2_list))
  )
  cv_rows <- list()
  times_para <- 0
  pb <- if (!isTRUE(verbose)) txtProgressBar(style = 3) else NULL

  for (para_1 in seq_along(alpha1_list)) {
    for (para_2 in seq_along(alpha2_list)) {
      cv_c <- numeric(cross_k)
      for (cvi in seq_len(cross_k)) {
        train_id <- setdiff(seq_len(datasize), cvlist[[cvi]])
        test_id <- cvlist[[cvi]]
        train <- train_data[train_id, , drop = FALSE]
        test <- train_data[test_id, , drop = FALSE]
        train_c <- train_phenotype[train_id, , drop = FALSE]
        test_c <- train_phenotype[test_id, , drop = FALSE]
        fn <- sqrt(sum(train * train))
        if (!is.finite(fn) || fn <= 0) {
          fn <- 1
        }
        train <- train / fn
        S <- make_scab_cv_S(train_c)
        Object_cv <- list(X = train, S = S, phenotype = train_c, A = A_cv, L = L_cv, D = D_cv, method = Object$method)
        class(Object_cv) <- "scAB_data"
        s_res <- scAB(
          Object = Object_cv,
          K = K,
          alpha = alpha1_list[para_1],
          alpha_2 = alpha2_list[para_2],
          maxiter = maxiter,
          seed = model_seed
        )
        ginvH <- MASS::ginv(s_res$H)
        new_W <- test %*% ginvH
        clin_km <- data.frame(time = train_c$time, status = train_c$status, s_res$W)
        res.cox <- survival::coxph(survival::Surv(time, status) ~ ., data = clin_km)
        pre_test <- as.numeric(stats::predict(res.cox, newdata = data.frame(new_W), type = "lp"))
        cv_c[cvi] <- as.numeric(survival::concordance(survival::Surv(test_c$time, test_c$status) ~ pre_test)$concordance)
        cv_rows[[length(cv_rows) + 1L]] <- data.frame(
          alpha1 = alpha1_list[para_1],
          alpha2 = alpha2_list[para_2],
          fold = cvi,
          cindex = cv_c[cvi],
          iter = s_res$iter,
          loss = s_res$loss,
          stringsAsFactors = FALSE
        )
      }
      result_cv[para_1, para_2] <- mean(cv_c)
      times_para <- times_para + 1L
      if (!is.null(pb)) {
        setTxtProgressBar(pb, times_para / length(alpha1_list) / length(alpha2_list))
      }
    }
  }
  if (!is.null(pb)) {
    close(pb)
  }

  para_index <- as.numeric(which(result_cv == max(result_cv, na.rm = TRUE), arr.ind = TRUE)[1, ])
  alpha_1 <- alpha1_list[para_index[1]]
  alpha_2 <- alpha2_list[para_index[2]]
  out <- list(alpha_1 = alpha_1, alpha_2 = alpha_2, result_cv = result_cv)

  if (isTRUE(return_details)) {
    cv_res <- do.call(rbind, cv_rows)
    cv_summary <- aggregate(cindex ~ alpha1 + alpha2, data = cv_res, FUN = mean)
    cv_summary <- cv_summary[order(-cv_summary$cindex, cv_summary$alpha1, cv_summary$alpha2), , drop = FALSE]
    rownames(cv_summary) <- NULL
    out <- c(out, list(
      cv_res = cv_res,
      cv_summary = cv_summary,
      best = cv_summary[1, , drop = FALSE],
      timing = NULL,
      parallel_info = NULL
    ))
  }

  out
}


###  identification cells above the threshold
#'
#' @param W individual matrix
#' @param H cell matrix
#' @param tred threshold
#'
#' @return a list with cells above the threshold
#' @importFrom  diptest dip.test
#' @importFrom multimode locmodes
#'
#' @export
#'
#' @examples
findModule <- function(H,tred=2,do.dip=FALSE){
  K = dim(H)[1]
  I=length(H)
  module=list()
  meanH=rowMeans(H)
  sdH=apply(H,1,sd)
  for(i in 1:K){
    x <- H[i,]
    if (isTRUE(do.dip) && diptest::dip.test(x)$p.value < 0.05)
    { modes <- multimode::locmodes(x, mod0 = 2)
      module=c(module,list(which( x > modes$locations[2])))}
    else {module=c(module,list( which( H[i,]-meanH[i] > tred*sdH[i] )) )}
  }
  return(module)
}



###  Subsets identification
#'
#' @param Object a Seurat object
#' @param scAB_Object a scAB_data object
#' @param tred threshold
#'
#' @return a Seurat object
#' @import Seurat
#' @export
#'
#' @examples
findSubset <- function(Object, scAB_Object, tred = 2){
  do.dip <- ifelse(scAB_Object$method=="binary",1,0)
  module <- findModule(scAB_Object$H, tred = tred, do.dip = do.dip)
  scAB_index <- unique(unlist(module))
  
  scAB_select <- rep("Other cells",ncol(Object))
  scAB_select[scAB_index] <- "scAB+ cells"
  Object <- Seurat::AddMetaData(Object,metadata = scAB_select, col.name = "scAB_select")
  
  for(i in 1:length(module)){
    M <- rep("Other cells",ncol(Object))
    M[as.numeric(module[[i]])]="scAB+ cells"
    Object <- Seurat::AddMetaData(Object,metadata = M, col.name = paste0("scAB_Subset",i))
    Object <- Seurat::AddMetaData(Object,metadata = scAB_Object$H[i,], col.name = paste0("Subset",i,"_loading"))
  }
  return(Object)
}
