# LargescAB

LargescAB is a large-scale engineering extension of `scAB`.  
It is designed to preserve the core scAB modeling workflow while improving
compatibility and scalability for modern single-cell datasets.

## What LargescAB Adds

- Seurat v5 / Assay5 compatibility
- Sparse graph handling
- HDF5-backed `X` support
- Blockwise quantile normalization
- Blockwise / streaming fitting
- Large-data `select_K_large()`
- Large-data `select_alpha_large()`
- Unified thread/process control

## Recommended API

Use this four-step workflow for large datasets:

1. `create_scAB_large()`
2. `select_K_large()` (optional if `K` is fixed)
3. `select_alpha_large()` (optional if `alpha`/`alpha_2` are fixed)
4. `scAB_large()`

The original dense interfaces are still available for compatibility:

- `create_scAB()`
- `scAB()`
- `select_K()`
- `select_alpha()`

## Typical Usage

```r
obj <- create_scAB_large(
  Obejct = sc_use,
  bulk_dataset = bulk_use,
  phenotype = Y_cox,
  method = "survival",
  block_size = 5000L,
  num_threads = 8L,
  x_parallel = FALSE,
  x_backend = "hdf5",
  x_h5_path = "scab_X.h5"
)

k_info <- select_K_large(
  Object = obj,
  K_max = 8,
  repeat_times = 3,
  maxiter = 100,
  k_parallel = TRUE,
  k_nworkers = 3L,
  return_details = TRUE
)

alpha_info <- select_alpha_large(
  Object = obj,
  K = k_info$selected_K,
  cross_k = 5,
  alpha1_list = c(0.01, 0.005, 0.001),
  alpha2_list = c(0.01, 0.005, 0.001),
  cv_parallel = TRUE,
  cv_nworkers = 24L,
  return_details = TRUE
)

fit <- scAB_large(
  Object = obj,
  K = k_info$selected_K,
  alpha = alpha_info$alpha_1,
  alpha_2 = alpha_info$alpha_2,
  maxiter = 2000,
  seed = 7L,
  materialize_hdf5 = FALSE,
  num_threads = 48L
)
```

## Practical Parameter Tiers

Frequently tuned:

- `K`
- `alpha`
- `alpha_2`
- `block_size`
- `num_threads`
- `x_parallel` / `x_nworkers`
- `k_parallel` / `k_nworkers`
- `cv_parallel` / `cv_nworkers`
- `materialize_hdf5`

Usually left as default:

- `assay = "RNA"`
- `graph_name = "RNA_snn"`
- `binarize_graph = TRUE`
- `convergence_threshold = 1e-5`
- `check_every = 1`

Debug or rescue only:

- `verify_equivalence`
- `eps`
- `guard_nonfinite`

Detailed parameter notes: [`PARAMETERS.md`](./PARAMETERS.md)

## Parallel Rules

- `num_threads` controls only single-process BLAS/OMP.
- `x_parallel` controls only preprocessing `X` block construction.
- `k_parallel` controls only repeated fits in `select_K_large()`.
- `cv_parallel` controls only CV task parallelism in `select_alpha_large()`.
- When task-level multiprocessing is enabled, worker-level `num_threads` is forced to `1`.
- Nested parallelism is intentionally disallowed in the first stable version.

## Relationship to Original scAB

LargescAB aims to preserve the core scAB methodology on valid inputs while
fixing compatibility and stability issues where needed (for example Seurat v5
access, missing graph safeguards, and zero-degree graph edge cases).

## Suggested GitHub Description

`LargescAB is a scalable implementation of scAB for large single-cell transcriptomic datasets. It is compatible with Seurat v5, reduces dense-memory bottlenecks, and is designed for datasets with large numbers of cells, including datasets with more than 100,000 cells. Its computational workflow preserves the core scAB methodology while introducing engineering changes for scalability and reproducibility.`
