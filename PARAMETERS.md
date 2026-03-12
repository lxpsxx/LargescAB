# LargescAB Parameter Quick Reference

This file summarizes practical parameter usage for the large-data workflow.

## `create_scAB_large()`

Commonly tuned:

- `method`
- `block_size`
- `num_threads`
- `x_parallel`
- `x_nworkers`
- `x_backend`
- `x_h5_path`

Usually left as default:

- `assay = "RNA"`
- `graph_name = "RNA_snn"`
- `binarize_graph = TRUE`
- `bulk_block_size = NULL`
- `x_h5_dataset = "X"`

Debug only:

- `verify_equivalence`

## `scAB_large()`

Commonly tuned:

- `K`
- `alpha`
- `alpha_2`
- `maxiter`
- `seed`
- `materialize_hdf5`
- `num_threads`

Usually left as default:

- `convergence_threshold = 1e-5`
- `x_block_size = NULL`
- `verbose = FALSE`
- `check_every = 1`

Debug or rescue only:

- `eps`
- `guard_nonfinite`

## `select_K_large()`

Commonly tuned:

- `K_max`
- `repeat_times`
- `maxiter`
- `k_parallel`
- `k_nworkers`

Usually left as default:

- `seed = 0`
- `materialize_hdf5`
- `return_details = FALSE`

## `select_alpha_large()`

Commonly tuned:

- `K`
- `cross_k`
- `alpha1_list`
- `alpha2_list`
- `maxiter`
- `cv_parallel`
- `cv_nworkers`

Usually left as default:

- `seed = 0`
- `model_seed = NULL`
- `return_details = FALSE`

## Parallel Rules

- `num_threads` is for single-process BLAS/OMP only.
- `x_parallel`, `k_parallel`, and `cv_parallel` are process-level parallel modes.
- If process-level parallelism is enabled, worker `num_threads` is forced to `1`.
- Nested parallelism is intentionally disabled in the current stable workflow.

## Practical Starter Settings

For checkpoint full runs:

- `block_size = 5000L` or `10000L`
- `x_backend = "hdf5"`
- `materialize_hdf5 = FALSE`
- `num_threads = 24L` or `48L`

For `K` selection:

- `K_max = 8`
- `repeat_times = 3`
- `maxiter = 100`
- `k_parallel = TRUE`

For alpha selection:

- `cross_k = 5`
- `alpha1_list = c(0.01, 0.005, 0.001)`
- `alpha2_list = c(0.01, 0.005, 0.001)`
- `cv_parallel = TRUE`
