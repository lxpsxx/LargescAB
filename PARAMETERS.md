# LargescAB 参数速查

这份文档只整理当前 large 主线真正会用到的参数。

## `create_scAB_large()`

### 实际常改

- `method`
  一般是 `"survival"`；如果是二分类任务才改 `"binary"`
- `block_size`
  单细胞分块大小，直接影响预处理峰值内存和速度
- `num_threads`
  单进程 BLAS/OMP 线程数
- `x_parallel`
  是否对 `X` 分块构建启用多进程
- `x_nworkers`
  `x_parallel=TRUE` 时的 worker 数
- `x_backend`
  `"hdf5"` 或 `"memory"`；大数据推荐 `"hdf5"`
- `x_h5_path`
  HDF5 `X` 的输出位置

### 一般不用动

- `assay = "RNA"`
- `graph_name = "RNA_snn"`
- `binarize_graph = TRUE`
- `bulk_block_size = NULL`
- `x_h5_dataset = "X"`

### 调试时才用

- `verify_equivalence`
  会额外和 dense 参考实现做对照，平时关闭

## `scAB_large()`

### 实际常改

- `K`
  必填，表示分解 rank / 细胞程序数
- `alpha`
  phenotype regularization 系数
- `alpha_2`
  graph regularization 系数
- `maxiter`
  最大迭代次数
- `seed`
  想严格复现时建议显式设置
- `materialize_hdf5`
  `FALSE` 表示按块流式读取 HDF5-backed `X`
- `num_threads`
  单进程 BLAS/OMP 线程数

### 一般不用动

- `convergence_threshold = 1e-5`
- `x_block_size = NULL`
- `verbose = FALSE`
- `check_every = 1`

### 调试时才用

- `eps`
  分母稳定项；默认保持 `0` 是为了贴近原始实现
- `guard_nonfinite`
  数值炸掉时的保底清理开关

## `select_K_large()`

### 实际常改

- `K_max`
  搜索上限
- `repeat_times`
  每个候选 `K` 的重复次数
- `maxiter`
  每次重复拟合的最大迭代
- `k_parallel`
  是否并行每个 `K` 的重复拟合
- `k_nworkers`
  `k_parallel=TRUE` 时 worker 数

### 一般不用动

- `seed = 0`
- `materialize_hdf5 = TRUE/FALSE`
  根据对象大小和内存定
- `return_details = FALSE`

## `select_alpha_large()`

### 实际常改

- `K`
  当前固定要评估的分解 rank
- `cross_k`
  CV 折数
- `alpha1_list`
  `alpha` 候选集合
- `alpha2_list`
  `alpha_2` 候选集合
- `maxiter`
  每个 fold 拟合的最大迭代
- `cv_parallel`
  是否并行 CV 任务
- `cv_nworkers`
  `cv_parallel=TRUE` 时 worker 数

### 一般不用动

- `seed = 0`
- `model_seed = NULL`
- `return_details = FALSE`

## 并行规则

- `num_threads` 只控制单进程 BLAS/OMP
- `x_parallel`、`k_parallel`、`cv_parallel` 都是多进程
- 只要开了多进程，worker 内 `num_threads` 自动强制成 `1`
- 第一版不支持嵌套并行

## 推荐起步配置

### checkpoint 全流程

- `block_size = 5000L` 或 `10000L`
- `x_backend = "hdf5"`
- `materialize_hdf5 = FALSE`
- `num_threads = 24L` 或 `48L`

### `K` 选择

- `K_max = 8`
- `repeat_times = 3`
- `maxiter = 100`
- `k_parallel = TRUE`

### `alpha` 选择

- `cross_k = 5`
- `alpha1_list = c(0.01, 0.005, 0.001)`
- `alpha2_list = c(0.01, 0.005, 0.001)`
- `cv_parallel = TRUE`
