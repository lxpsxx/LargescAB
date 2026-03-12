# LargescAB

LargescAB is a scalable fork of `scAB` for integrating single-cell and bulk
RNA-seq data with phenotype information.

It preserves the core scAB modeling workflow while adding engineering
improvements required for large single-cell datasets, including Seurat v5
compatibility and memory-aware large-matrix handling.

## Key Capabilities

- Seurat v5 / Assay5 compatibility
- Sparse graph support for cell-cell regularization
- HDF5-backed `X` matrix support
- Blockwise preprocessing and fitting paths
- Large-data `K` selection (`select_K_large()`)
- Large-data alpha selection (`select_alpha_large()`)
- Explicit thread/process parallel controls

## Installation

```r
# install.packages("devtools")
devtools::install_github("lxpsxx/LargescAB")
```

Core dependencies are listed in `DESCRIPTION`.

## Main API

Large-data production workflow:

1. `create_scAB_large()`
2. `select_K_large()` (optional)
3. `select_alpha_large()` (optional)
4. `scAB_large()`

Original dense interfaces are retained for compatibility:
`create_scAB()`, `scAB()`, `select_K()`, and `select_alpha()`.

## Input Requirements

- `sc_use`: Seurat object with an expression `data` layer and graph
  (default graph name: `RNA_snn`)
- `bulk_use`: numeric matrix in `gene x sample` format
- `Y_cox`: phenotype table aligned to `colnames(bulk_use)`
  - survival mode: columns `time` and `status`

## Quick Start

```r
obj <- create_scAB_large(
  Obejct = sc_use,
  bulk_dataset = bulk_use,
  phenotype = Y_cox,
  method = "survival",
  block_size = 5000L,
  x_backend = "hdf5",
  x_h5_path = "scab_X.h5",
  num_threads = 8L
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
  materialize_hdf5 = FALSE,
  num_threads = 48L
)
```

For a concise parameter guide, see [`PARAMETERS.md`](./PARAMETERS.md).

## Reproducibility and Parallel Notes

- Use `seed` for controlled initialization/reproducibility.
- `num_threads` controls single-process BLAS/OMP threads.
- `x_parallel`, `k_parallel`, and `cv_parallel` are task-level parallel modes.
- Worker-level `num_threads` is forced to `1` when task-level parallelism is enabled.

## License and Attribution

LargescAB is distributed under `GPL-3` and is based on the original `scAB`
package. Please keep attribution and license terms when redistributing.
