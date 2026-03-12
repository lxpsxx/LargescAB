sanitize_scab_positive_int <- function(x, arg_name) {
  if (length(x) != 1 || is.na(x) || !is.finite(x) || x < 1) {
    stop(sprintf("`%s` must be a single positive integer.", arg_name))
  }
  as.integer(x)
}


restore_scab_envvar <- function(name, value) {
  if (length(value) != 1 || is.na(value) || !nzchar(value)) {
    Sys.unsetenv(name)
  } else {
    do.call(Sys.setenv, stats::setNames(list(value), name))
  }
  invisible(NULL)
}


resolve_scab_parallel_config <- function(
  num_threads = NULL,
  use_parallel = FALSE,
  nworkers = 1L,
  parallel_label = "parallel"
) {
  use_parallel <- isTRUE(use_parallel)
  nworkers_requested <- sanitize_scab_positive_int(nworkers, "nworkers")
  num_threads_requested <- if (is.null(num_threads)) {
    NULL
  } else {
    sanitize_scab_positive_int(num_threads, "num_threads")
  }

  parallel_enabled <- use_parallel && nworkers_requested > 1L
  if (use_parallel && nworkers_requested <= 1L) {
    warning(
      sprintf("`%s=TRUE` but `nworkers<=1`; falling back to serial execution.", parallel_label),
      call. = FALSE
    )
  }

  if (parallel_enabled && .Platform$OS.type != "unix") {
    stop(sprintf("`%s` currently requires a Unix-like OS for mclapply.", parallel_label))
  }

  force_single_thread <- parallel_enabled
  num_threads_effective <- if (force_single_thread) {
    1L
  } else if (is.null(num_threads_requested)) {
    NULL
  } else {
    num_threads_requested
  }

  if (force_single_thread && !is.null(num_threads_requested) && num_threads_requested != 1L) {
    warning(
      sprintf(
        "`%s=TRUE` forces `num_threads=1` inside each worker; requested `%d`, using `%d`.",
        parallel_label,
        num_threads_requested,
        num_threads_effective
      ),
      call. = FALSE
    )
  }

  list(
    parallel_enabled = parallel_enabled,
    nworkers_requested = nworkers_requested,
    nworkers_effective = if (parallel_enabled) nworkers_requested else 1L,
    num_threads_requested = num_threads_requested,
    num_threads_effective = num_threads_effective,
    force_single_thread = force_single_thread,
    parallel_label = parallel_label
  )
}


validate_scab_parallel_combination <- function(
  cv_parallel = FALSE,
  x_parallel = FALSE,
  context = "this workflow"
) {
  if (isTRUE(cv_parallel) && isTRUE(x_parallel)) {
    stop(
      sprintf(
        "Nested parallelism is not supported in %s: `cv_parallel=TRUE` and `x_parallel=TRUE` cannot be enabled together.",
        context
      )
    )
  }
  invisible(TRUE)
}


with_scab_num_threads <- function(num_threads, expr) {
  expr_sub <- substitute(expr)

  if (is.null(num_threads)) {
    return(eval.parent(expr_sub))
  }

  num_threads <- sanitize_scab_positive_int(num_threads, "num_threads")
  env_names <- c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS"
  )
  old_env <- Sys.getenv(env_names, unset = NA_character_)
  old_blas <- NA_integer_
  old_omp <- NA_integer_
  has_blasctl <- requireNamespace("RhpcBLASctl", quietly = TRUE)

  if (has_blasctl) {
    old_blas <- suppressWarnings(as.integer(RhpcBLASctl::blas_get_num_procs()))
    old_omp <- suppressWarnings(as.integer(RhpcBLASctl::omp_get_max_threads()))
  }

  on.exit({
    for (nm in names(old_env)) {
      restore_scab_envvar(nm, old_env[[nm]])
    }
    if (has_blasctl) {
      if (is.finite(old_blas)) {
        RhpcBLASctl::blas_set_num_threads(old_blas)
      }
      if (is.finite(old_omp)) {
        RhpcBLASctl::omp_set_num_threads(old_omp)
      }
    }
  }, add = TRUE)

  thread_value <- as.character(num_threads)
  do.call(
    Sys.setenv,
    stats::setNames(as.list(rep(thread_value, length(env_names))), env_names)
  )

  if (has_blasctl) {
    RhpcBLASctl::blas_set_num_threads(num_threads)
    RhpcBLASctl::omp_set_num_threads(num_threads)
  }

  eval.parent(expr_sub)
}
