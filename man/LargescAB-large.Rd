\name{LargescAB-large}
\alias{create_scAB_large}
\alias{scAB_large}
\alias{select_K_large}
\alias{select_alpha_large}
\title{Large-data entry points in LargescAB}
\description{
The large-data workflow in \code{LargescAB} keeps the original scAB
factorization semantics while replacing dense-memory bottlenecks with sparse,
blockwise, and HDF5-backed implementations.
}
\usage{
create_scAB_large(
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
)

scAB_large(
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
)

select_K_large(
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
)

select_alpha_large(
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
)
}
\arguments{
\item{Obejct}{A Seurat object. The \code{RNA} assay data layer and the graph
\code{RNA_snn} are used by default.}

\item{bulk_dataset}{A numeric bulk matrix in gene-by-sample format.}

\item{phenotype}{For \code{method = "survival"}, a data frame with row names
matching \code{colnames(bulk_dataset)} and columns \code{time} and
\code{status}. For \code{method = "binary"}, a logical or numeric label
vector.}

\item{method}{One of \code{"survival"} or \code{"binary"}.}

\item{assay}{Assay name used to read the single-cell expression layer.}

\item{graph_name}{Graph name used as the cell-cell similarity graph.}

\item{binarize_graph}{Whether to convert graph edge weights to 0/1 before
normalization.}

\item{block_size}{Cell-column block size used during preprocessing and HDF5
construction.}

\item{bulk_block_size}{Optional block size used for bulk-side quantile target
construction.}

\item{num_threads}{Single-process BLAS/OMP thread count. When a task-level
parallel mode is enabled, worker-level \code{num_threads} is forced to 1.}

\item{x_parallel}{Whether to parallelize \code{X} block construction in
\code{create_scAB_large()}.}

\item{x_nworkers}{Worker count for \code{x_parallel = TRUE}.}

\item{x_backend}{Whether \code{X} should be stored in memory or in HDF5.}

\item{x_h5_path}{Optional HDF5 file path for \code{X}.}

\item{x_h5_dataset}{Dataset name used inside the HDF5 file.}

\item{verify_equivalence}{Debug-only check that compares the blockwise
preprocessing result against the dense reference implementation.}

\item{Object}{A \code{scAB_data} object returned by
\code{create_scAB_large()} or a compatible large object.}

\item{K}{Matrix factorization rank / number of inferred programs.}

\item{alpha}{Phenotype regularization coefficient.}

\item{alpha_2}{Graph regularization coefficient.}

\item{maxiter}{Maximum iteration count used by fitting or selection sub-fits.}

\item{convergence_threshold}{Stopping threshold on the loss change.}

\item{eps}{Optional denominator stability constant. The default remains 0 to
stay close to the original implementation.}

\item{seed}{Random seed. In \code{scAB_large()} it controls factor
initialization; in \code{select_K_large()} and \code{select_alpha_large()} it
controls repeated fits or CV fold generation.}

\item{guard_nonfinite}{Whether to zero-out non-finite or negative values
after updates. Intended for rescue/debug use.}

\item{verbose}{Whether to print progress or diagnostic messages.}

\item{check_every}{How often to evaluate convergence.}

\item{materialize_hdf5}{Whether to fully materialize HDF5-backed \code{X}
before fitting. Set \code{FALSE} for large-data streaming.}

\item{x_block_size}{Optional fit-time HDF5 block size.}

\item{K_max}{Largest candidate \code{K} evaluated by \code{select_K_large()}.}

\item{repeat_times}{Number of repeated fits for each candidate \code{K}.}

\item{k_parallel}{Whether to parallelize repeated fits inside each candidate
\code{K}.}

\item{k_nworkers}{Worker count for \code{k_parallel = TRUE}.}

\item{return_details}{Whether to return detailed diagnostics such as per-fold
or per-repeat tables.}

\item{cross_k}{Number of folds used by \code{select_alpha_large()}.}

\item{alpha1_list}{Candidate values for \code{alpha}.}

\item{alpha2_list}{Candidate values for \code{alpha_2}.}

\item{model_seed}{Optional seed passed to the model fit inside each CV fold.}

\item{cv_parallel}{Whether to parallelize the CV task grid in
\code{select_alpha_large()}.}

\item{cv_nworkers}{Worker count for \code{cv_parallel = TRUE}.}
}
\value{
\code{create_scAB_large()} returns a large-data \code{scAB_data} object.
\code{scAB_large()} returns a list containing at least \code{W}, \code{H},
\code{iter}, and \code{loss}. \code{select_K_large()} returns either the
selected \code{K} or a detail list. \code{select_alpha_large()} returns the
selected \code{alpha_1}, \code{alpha_2}, and optionally CV detail tables.
}
\details{
Recommended production workflow:
\enumerate{
  \item Build a large object with \code{create_scAB_large()}.
  \item Select \code{K} with \code{select_K_large()} when needed.
  \item Select \code{alpha/alpha_2} with \code{select_alpha_large()} when needed.
  \item Fit the final model with \code{scAB_large()}.
}

Task-level parallel modes and \code{num_threads} are intentionally separated:
\enumerate{
  \item \code{x_parallel} only affects preprocessing.
  \item \code{k_parallel} only affects repeated \code{K} fits.
  \item \code{cv_parallel} only affects the alpha-selection CV grid.
  \item Any worker process forces \code{num_threads = 1} internally.
}
}
\examples{
\dontrun{
obj <- create_scAB_large(
  Obejct = sc_use,
  bulk_dataset = bulk_use,
  phenotype = Y_cox,
  method = "survival",
  block_size = 5000L,
  x_backend = "hdf5"
)

k_info <- select_K_large(Object = obj, K_max = 8, repeat_times = 3, return_details = TRUE)

alpha_info <- select_alpha_large(
  Object = obj,
  K = k_info$selected_K,
  return_details = TRUE
)

fit <- scAB_large(
  Object = obj,
  K = k_info$selected_K,
  alpha = alpha_info$alpha_1,
  alpha_2 = alpha_info$alpha_2,
  materialize_hdf5 = FALSE
)
}
}
