# LargescAB

`LargescAB` 是 `scAB` 的大规模工程化扩展版本。它优先解决 Seurat v5 兼容性、
dense 内存瓶颈和大规模单细胞输入的可运行性，同时尽量保留 `scAB` 的核心计算逻辑。

## 当前定位

这个包已经不是单纯的兼容补丁，而是当前项目的正式算法实现，重点包括：

- Seurat v5 / Assay5 兼容
- 稀疏图矩阵
- HDF5-backed `X`
- blockwise quantile normalization
- blockwise / streaming 拟合
- large-data `select_K_large()`
- large-data `select_alpha_large()`
- 线程和多进程规则统一

## 推荐 API

主线入口固定为四个函数：

- `create_scAB_large()`
- `scAB_large()`
- `select_K_large()`
- `select_alpha_large()`

兼容保留的原始入口仍然存在：

- `create_scAB()`
- `scAB()`
- `select_K()`
- `select_alpha()`

但大规模数据不建议再走 dense 主线。

## 最常用调用方式

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

## 参数分层

常改参数：

- `K`
- `alpha`
- `alpha_2`
- `block_size`
- `num_threads`
- `x_parallel / x_nworkers`
- `k_parallel / k_nworkers`
- `cv_parallel / cv_nworkers`
- `materialize_hdf5`

一般不用动：

- `assay = "RNA"`
- `graph_name = "RNA_snn"`
- `binarize_graph = TRUE`
- `convergence_threshold = 1e-5`
- `check_every = 1`

调试时才改：

- `verify_equivalence`
- `eps`
- `guard_nonfinite`

详细参数解释见：

- [`PARAMETERS.md`](./PARAMETERS.md)
- [`man/LargescAB-large.Rd`](./man/LargescAB-large.Rd)

## 并行规则

- `num_threads`
  只控制单进程 BLAS/OMP
- `x_parallel`
  只控制 `X` 分块构建
- `k_parallel`
  只控制 `select_K_large()` 的 repeat 并行
- `cv_parallel`
  只控制 `select_alpha_large()` 的 CV 任务并行

约束：

- 一旦进入多进程 worker，`num_threads` 会被强制改成 `1`
- 第一版不允许嵌套并行
- `cv_parallel` 和 `x_parallel` 不能在同一条工作流里同时启用

## 与原始 scAB 的关系

设计原则是：

- 对正常输入，尽量保持原始算法语义和结果
- 对原包本身有 bug 或边界不稳的场景，修掉错误行为

已经确认修复的典型问题包括：

- Seurat v5 `data` 层访问不兼容
- 缺失 `RNA_snn` 时继续误跑
- 零度节点导致 `Inf/NA`
- 原始 `scAB()` seed 写死且不可控

## GitHub 描述建议

推荐使用下面这段更稳的英文描述：

`LargescAB is a scalable implementation of scAB for large single-cell transcriptomic datasets. It is compatible with Seurat v5, reduces dense-memory bottlenecks, and is designed for datasets with large numbers of cells, including datasets with more than 100,000 cells. Its computational workflow preserves the core scAB methodology while introducing engineering changes for scalability and reproducibility.`

这版比 `strictly follows the original scAB methodology` 更准确，因为当前实现确实做了：

- Seurat v5 兼容层
- HDF5 / blockwise 重写
- 稀疏图表示
- 并行和线程控制

这些都不改变主算法目标，但属于明确的工程层改造。

## 当前工程配套

包本体只负责算法和对象操作。正式项目运行入口见：

- [`../../02脚本/40_run_largerscab_standard.R`](../../02脚本/40_run_largerscab_standard.R)
