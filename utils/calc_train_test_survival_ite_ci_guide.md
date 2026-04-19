# `calc_train_test_survival_ite_ci()` 函数说明文档

## 📌 文档目的

本文档用于详细说明 `test.val.R` 中最终版本函数 `calc_train_test_survival_ite_ci()` 的设计目标、计算步骤、输入输出、bootstrap 实现逻辑以及推荐使用方式。

该函数的核心目标是：

- 基于 `personalized` 包在训练集上拟合生存结局个体化治疗效应模型
- 同时计算训练集和测试集的逐患者 ITE 点估计
- 利用 bootstrap 重采样对训练集反复重建模型
- 在每一轮 bootstrap 中同时预测原始训练集和固定测试集
- 最终分别给出训练集和测试集逐患者 ITE 的可信区间

这里的 ITE 口径与当前项目原有脚本保持一致，即：

- 将 `personalized::predict(..., type = "benefit.score")` 输出作为逐患者 ITE 的代理量

因此，本函数更准确地说，是用于计算：

- 逐患者 `benefit score` 点估计
- 逐患者 `benefit score` 的 bootstrap 百分位法可信区间

***

## 🧭 函数所在位置

- 主脚本：`test.val.R`
- 结果整理函数：`build_ite_result_table()`
- 主函数：`calc_train_test_survival_ite_ci()`

***

## 🔄 整体计算流程

下面用流程图说明该函数的总体执行逻辑。

```mermaid
flowchart TD
    accTitle: 训练测试集 ITE 计算流程
    accDescr: 该流程展示函数如何读取训练与测试数据，统一构建设计矩阵，在训练集拟合主模型，并通过bootstrap同时得到训练集和测试集的逐患者ITE点估计与可信区间。

    start["开始"]
    check["检查训练集/测试集字段与缺失值"]
    matrix["统一训练集与测试集设计矩阵编码"]
    fit0["完整训练集拟合主模型"]
    pred0["主模型预测训练集与测试集 ite_hat"]
    boot["训练集进行 B 次 bootstrap 重采样"]
    refit["每轮 bootstrap 重拟合 personalized 模型"]
    pred_both["每轮同时预测原始训练集与固定测试集"]
    quantile["对每位患者计算 2.5% 与 97.5% 分位数"]
    output["输出训练集与测试集结果表"]
    end["结束"]

    start --> check
    check --> matrix
    matrix --> fit0
    fit0 --> pred0
    pred0 --> boot
    boot --> refit
    refit --> pred_both
    pred_both --> quantile
    quantile --> output
    output --> end
```

***

## 🧪 为什么要“一次 bootstrap 同时计算训练集和测试集”

你当前的需求不是只看测试集，而是希望：

- 训练集也有逐患者 ITE 点估计和区间
- 测试集也有逐患者 ITE 点估计和区间
- 两者尽量来自同一套 bootstrap 过程，保证分析流程统一

本函数的设计正是围绕这一点展开：

- 主模型只在完整训练集上拟合一次，用于给训练集和测试集分别计算点估计
- bootstrap 每一轮都只对训练集做有放回重采样
- 但每一轮训练出的模型会同时预测：
- 原始训练集
- 固定测试集

这样做的好处是：

- 训练集和测试集区间来自同一套重采样机制
- 训练/测试两套结果口径一致
- 不需要写两套重复逻辑
- 后续画图、排序、筛选显著患者会非常方便

***

## 🧾 输入数据要求

函数要求同时提供：

- `train_data`：训练集数据框
- `test_data`：测试集数据框

训练集和测试集都必须包含以下几类字段：

- 患者 ID 字段，例如 `ID`
- 生存时间字段，例如 `OS`
- 结局状态字段，例如 `OS_CNSR`
- 治疗分组字段，例如 `Trt_ITT`
- 建模协变量字段，即 `feature_cols` 指定的全部变量

当前脚本默认使用的协变量为：

```r
feature_cols_default <- c(
  "Tumor_border_status",
  "CD4_CD8",
  "Peritumoral_enhancement",
  "CRP",
  "CPgrade",
  "Radiologic_morphology",
  "NLR",
  "CAperihepatic_99",
  "Neutro",
  "Cr",
  "PDL1",
  "DBIL",
  "BUN",
  "ALB",
  "TumorDiameterMax",
  "Vascular_invasion",
  "Tumor_deposit"
)
```

另外还要求：

- `OS_CNSR` 必须编码为 `0/1`
- `Trt_ITT` 必须编码为 `0/1`
- 参与建模的关键字段不能有缺失值

***

## 🛠 关键辅助函数说明

### 1. `assert_required_columns()`

作用：

- 检查训练集或测试集是否缺少必要字段

如果缺少列，函数会立即报错，避免后面建模时报出不易定位的错误。

### 2. `assert_no_missing()`

作用：

- 检查关键字段中是否存在缺失值

如果某些建模字段存在缺失值，会直接停止，并指出是哪几列有多少缺失。

### 3. `build_aligned_design_matrices()`

作用：

- 统一训练集和测试集的协变量编码方式

这是非常关键的一步。因为如果训练集和测试集分别做哑变量展开，分类变量的水平可能不一致，导致：

- 训练集矩阵列数与测试集矩阵列数不同
- 某些分类变量在测试集中缺少水平
- 预测时 `newx` 与训练时设计矩阵结构不匹配

本函数的处理方式是：

1. 先取出训练集和测试集的原始协变量
2. 纵向拼接训练集和测试集
3. 对字符变量与逻辑变量统一转为因子
4. 用 `model.matrix(~ . - 1)` 一次性生成设计矩阵
5. 再按行切回训练集矩阵和测试集矩阵

这样可以保证：

- 两套数据的哑变量列完全对齐
- 分类变量的编码方式完全一致

***

## ⚙ 主函数计算步骤详解

## 1. 参数初始化

函数开始时会先对参数做基础规范化：

- `B` 转成整数
- `nfolds` 转成整数
- `prop_nfolds_crossfit` 转成整数

并检查：

- `B >= 2`
- `nfolds >= 2`

***

## 2. 训练集与测试集字段检查

函数会分别对训练集和测试集执行：

- 字段完整性检查
- 缺失值检查
- 治疗变量与状态变量编码检查

这样可以在真正拟合模型前尽早发现问题。

***

## 3. 构建设计矩阵

函数通过 `build_aligned_design_matrices()` 完成：

- 训练集设计矩阵 `x_train`
- 测试集设计矩阵 `x_test`

同时生成：

- 训练集样本量 `n_train`
- 测试集样本量 `n_test`

***

## 4. 构造训练集生存对象

使用：

```r
Surv(time = time_train, event = status_train)
```

构造训练集生存结局对象 `y_train`。

***

## 5. 创建倾向评分函数

函数内部使用：

```r
personalized::create.propensity.function()
```

生成倾向评分函数，并将其传入 `fit.subgroup()`。

当前脚本使用的是：

- `crossfit = TRUE`
- `nfolds.crossfit = 5`
- `cv.glmnet.args = list(type.measure = "auc", nfolds = 5)`

这与 `personalized` 文档的推荐方向是一致的：高维或观察性数据中，倾向评分估计通常建议配合交叉拟合与内部交叉验证。

***

## 6. 在完整训练集上拟合主模型

函数通过内部封装的：

```r
fit_subgroup_once()
```

执行：

```r
personalized::fit.subgroup(
  x = x_train,
  y = y_train,
  trt = trt_train,
  propensity.func = prop_func,
  loss = "cox_loss_lasso",
  nfolds = nfolds
)
```

这里模型只在完整训练集上拟合一次，作为主模型。

***

## 7. 主模型点估计

主模型拟合完成后，分别对：

- 原始训练集 `x_train`
- 固定测试集 `x_test`

进行预测：

```r
predict(object = base_model, newx = x_train, type = "benefit.score")
predict(object = base_model, newx = x_test, type = "benefit.score")
```

得到：

- `ite_hat_train`
- `ite_hat_test`

这两个量就是训练集和测试集的逐患者点估计。

***

## 8. bootstrap 重采样重建模型

函数执行 `B` 次循环。

每一轮流程如下：

```mermaid
flowchart LR
    accTitle: 单轮 Bootstrap 流程
    accDescr: 单轮bootstrap中，先对训练集有放回重采样，再用重采样训练集拟合模型，最后同时预测原始训练集和固定测试集。

    sample["训练集有放回重采样"]
    fit["在bootstrap训练集上重拟合模型"]
    pred_train["预测原始训练集"]
    pred_test["预测固定测试集"]
    store["保存本轮预测结果"]

    sample --> fit
    fit --> pred_train
    fit --> pred_test
    pred_train --> store
    pred_test --> store
```

具体而言：

### 第一步：训练集有放回抽样

从训练集行索引中抽取 `n_train` 个样本，允许重复：

```r
boot_idx <- sample.int(n = n_train, size = n_train, replace = TRUE)
```

### 第二步：构造本轮 bootstrap 生存对象

使用重采样后的训练样本重新构造：

```r
Surv(time = time_train[boot_idx], event = status_train[boot_idx])
```

### 第三步：在 bootstrap 样本上拟合模型

重新调用 `fit_subgroup_once()`。

### 第四步：同时预测原始训练集和固定测试集

如果本轮模型拟合成功，则分别预测：

- 原始训练集 `x_train`
- 固定测试集 `x_test`

### 第五步：保存到两张预测矩阵

分别保存到：

- `boot_pred_train_mat`
- `boot_pred_test_mat`

矩阵结构为：

- 行：患者
- 列：bootstrap 轮次

***

## 9. 为什么预测“原始训练集”而不是“bootstrap 样本自身”

这是本函数的重要设计点。

如果每轮 bootstrap 只预测重采样样本自身，会有几个问题：

- 同一患者在某轮可能出现多次
- 不同轮次的患者集合不完全一致
- 很难对“原始训练集每一位患者”构造稳定的逐患者区间

因此，本函数采用更稳定的做法：

- 每轮 bootstrap 虽然在训练集上重采样建模
- 但模型一旦建好，统一去预测原始训练集所有患者

这样训练集和测试集都具有统一口径：

- 每位患者在每轮尽可能对应一个预测值
- 因而能直接对每位患者做分位数区间计算

***

## 10. 可信区间计算

对于训练集和测试集，函数都调用：

- `build_ite_result_table()`

这一步会逐患者计算：

- `ite_ci_lower`：2.5% 分位数
- `ite_ci_upper`：97.5% 分位数
- `boot_valid_B`：非缺失预测的次数
- `ite_ci_width`：区间宽度

因此当前区间类型是：

- bootstrap 百分位法可信区间

不是：

- 正态近似区间
- BCa 区间
- 偏差校正区间

***

## 11. 结果表的附加解释字段

每张结果表还会额外生成：

### `recommended_trt`

规则为：

```r
ite_hat > 0 -> 1
ite_hat <= 0 -> 0
```

即：

- 预测获益大于 0 时推荐治疗 `1`
- 否则推荐治疗 `0`

### `ci_excludes_zero`

用于判断区间是否不跨 `0`：

- 若 `ite_ci_lower > 0` 或 `ite_ci_upper < 0`，则记为 `TRUE`
- 否则记为 `FALSE`

这通常可以作为“该患者个体化获益方向较稳定”的一个简单标记。

### `ite_rank_desc`

按 `ite_hat` 从大到小排序的排名。

***

## 📤 输出对象说明

函数最终返回一个列表，包含：

- `base_model`：完整训练集拟合得到的主模型
- `train_results_tbl`：训练集结果表
- `test_results_tbl`：测试集结果表
- `boot_pred_train_mat`：训练集 bootstrap 预测矩阵
- `boot_pred_test_mat`：测试集 bootstrap 预测矩阵
- `design_colnames`：最终设计矩阵列名

***

## 📁 文件输出说明

函数支持三类输出：

### 1. 训练集 CSV

由参数：

- `output_train_csv`

控制。

### 2. 测试集 CSV

由参数：

- `output_test_csv`

控制。

### 3. Excel 汇总文件

由参数：

- `output_xlsx`

控制。

Excel 文件中包含两个 sheet：

- `train_ite_results`
- `test_ite_results`

***

## 🧷 当前脚本中的默认调用方式

当前 `test.val.R` 脚本底部已写好一个可直接运行的示例。

主要配置如下：

- 训练集文件：`data/final_data_train.xlsx`
- 测试集文件：`data/final_data_test.xlsx`
- 建模变量：上面列出的 `17` 个变量
- `loss = "cox_loss_lasso"`
- `B = 100`
- `nfolds = 5`
- `prop_crossfit = TRUE`

调用示例为：

```r
train_test_ite_res <- calc_train_test_survival_ite_ci(
  train_data = train_data,
  test_data = test_data,
  feature_cols = feature_cols_default,
  time_col = "OS",
  status_col = "OS_CNSR",
  trt_col = "Trt_ITT",
  id_col = "ID",
  B = 100L,
  nfolds = 5L,
  prop_crossfit = TRUE,
  prop_nfolds_crossfit = 5L,
  prop_cv_metric = "auc",
  loss = "cox_loss_lasso",
  seed = 20260418L,
  output_train_csv = output_train_csv_file,
  output_test_csv = output_test_csv_file,
  output_xlsx = output_xlsx_file,
  verbose = TRUE
)
```

***

## ✅ 当前项目中已生成的结果文件

本次 `100` 次 bootstrap 分析实际生成了以下文件：

- `results/survival_ite_ci_train_results_boot100.csv`
- `results/survival_ite_ci_test_results_boot100.csv`
- `results/survival_ite_ci_train_test_results_boot100.xlsx`

***

## 📌 结果字段解释

训练集和测试集结果表字段完全一致：

| 字段名                | 含义                               |
| ------------------ | -------------------------------- |
| `dataset`          | 数据来源，`train` 或 `test`            |
| `patient_id`       | 患者 ID                            |
| `observed_time`    | 生存时间                             |
| `observed_status`  | 生存结局状态                           |
| `observed_trt`     | 实际治疗分组                           |
| `ite_hat`          | 逐患者 ITE 点估计代理量，即 `benefit.score` |
| `ite_ci_lower`     | 95% 区间下限                         |
| `ite_ci_upper`     | 95% 区间上限                         |
| `boot_valid_B`     | 成功获得预测的 bootstrap 次数             |
| `ite_ci_width`     | 区间宽度                             |
| `recommended_trt`  | 根据 `ite_hat` 得到的推荐治疗             |
| `ci_excludes_zero` | 区间是否不跨 0                         |
| `ite_rank_desc`    | 按 `ite_hat` 从高到低的排名              |

***

## ⚠ 方法学注意事项

### 1. 这里的 ITE 不是严格意义上的“真实个体治疗效应”

当前输出的是：

- `benefit.score`

它是 `personalized` 框架下用于个体化治疗推荐的获益评分，可视为逐患者 ITE 的代理量，但并不等同于真实不可观测的反事实治疗效应。

### 2. 当前区间是 bootstrap 百分位区间

即：

- 直接对每位患者 bootstrap 预测分布取分位数

优点是直观、实现清晰；缺点是：

- 对偏态分布或极端值时可能不如 BCa 等区间稳健

### 3. 当前训练集区间属于“重建模型后回推原始训练集”

因此它反映的是：

- 训练样本在模型不稳定性下的预测波动

而不是：

- 完全独立外部验证下的区间

### 4. 测试集区间更接近外部样本预测不确定性

因为测试集本身不参与每轮 bootstrap 拟合，只作为固定外部样本接受预测。

***

## 🚀 适合的使用场景

该函数适合以下场景：

- 已有明确训练集与测试集划分
- 需要同时得到训练集和测试集逐患者 ITE 点估计
- 需要给每位患者的 ITE 提供一个 bootstrap 区间
- 后续还要进行：
- 排名图绘制
- 高获益人群筛选
- 区间是否跨 0 的分层分析
- 训练/测试结果并排对比

***

## 🧠 推荐后续扩展

如果你后续还要继续完善，可以考虑增加：

- 支持自定义分位数区间，例如 `90%`、`99%`
- 支持 BCa 或 bias-corrected 区间
- 支持并行 bootstrap
- 支持自动绘制训练/测试两套排序图
- 支持将模型参数与运行摘要自动写入单独 sheet

***

## 📚 与 `personalized` 官方验证函数的关系

`personalized::validate.subgroup()` 官方更偏向：

- 亚组层面的治疗效应验证
- 重复训练/测试拆分验证
- bootstrap bias correction

而当前函数更偏向：

- 逐患者层面的 `benefit.score`
- 逐患者点估计
- 逐患者 bootstrap 区间

因此二者是互补关系，不是替代关系。

***

## ✅ 总结

`calc_train_test_survival_ite_ci()` 的本质可以概括为一句话：

- 用完整训练集建立 `personalized` 生存获益模型，再通过训练集 bootstrap 重建模型的不确定性，同时量化训练集和固定测试集中每一位患者的 `benefit score` 点估计与可信区间。

它的优势在于：

- 训练集与测试集口径统一
- 单次 bootstrap 流程同时覆盖两套数据
- 结果结构清晰，便于后续画图、排序和汇报
- 已经在当前项目数据上完成 `100` 次重采样实跑验证，代码可正常工作

