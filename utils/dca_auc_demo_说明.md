# DCA 曲线面积计算工具

本项目提供了两个R脚本，用于计算决策曲线分析（Decision Curve Analysis, DCA）的曲线面积，分别适用于二分类结局和生存结局。

## 功能说明

- **二分类结局**：计算二分类结局的DCA曲线面积，包括正面积、负面积和有符号总面积
- **生存结局**：计算生存结局的DCA曲线面积，同样包括正面积、负面积和有符号总面积

## 依赖包

运行本项目需要安装以下R包：

- `dcurves`：提供DCA主函数
- `tidyverse`：用于数据处理
- `survival`（仅生存分析脚本需要）：提供Surv与Cox拟合工具

## 安装依赖

```R
install.packages(c("dcurves", "tidyverse", "survival"))
```

## 文件结构

- `dca_binary_auc_demo.R`：二分类结局DCA面积计算脚本
- `dca_survival_auc_demo.R`：生存结局DCA面积计算脚本

## 函数说明

### 1. 二分类结局DCA面积计算

#### 函数签名

```R
calc_dca_auc_binary(
  data,           # 数据框
  outcome,        # 结局列名（二分类变量）
  model_cols,     # 模型预测列名向量
  thresholds = seq(0.01, 0.99, by = 0.01)  # 阈值序列
)
```

#### 参数说明

- `data`：包含结局变量和模型预测值的数据框
- `outcome`：二分类结局变量的列名
- `model_cols`：包含一个或多个模型预测值的列名向量
- `thresholds`：决策阈值序列，默认从0.01到0.99，步长0.01

#### 返回值

返回一个数据框，包含以下列：
- `label`：模型标签
- `dca_auc_pos`：正面积（曲线在0线上方的面积）
- `dca_auc_neg`：负面积（曲线在0线下方的面积）
- `dca_auc_sum_signed`：有符号总面积（正面积+负面积）

### 2. 生存结局DCA面积计算

#### 函数签名

```R
calc_dca_auc_survival(
  data,           # 数据框
  time_col,       # 时间列名
  event_col,      # 事件列名
  model_cols,     # 模型预测列名向量
  eval_time,      # 评估时间点
  thresholds = seq(0.01, 0.99, by = 0.01)  # 阈值序列
)
```

#### 参数说明

- `data`：包含时间、事件和模型预测值的数据框
- `time_col`：时间变量的列名
- `event_col`：事件变量的列名
- `model_cols`：包含一个或多个模型预测值的列名向量
- `eval_time`：评估时间点
- `thresholds`：决策阈值序列，默认从0.01到0.99，步长0.01

#### 返回值

返回一个数据框，包含以下列：
- `label`：模型标签
- `dca_auc_pos`：正面积（曲线在0线上方的面积）
- `dca_auc_neg`：负面积（曲线在0线下方的面积）
- `dca_auc_sum_signed`：有符号总面积（正面积+负面积）

## 使用示例

### 1. 二分类结局示例

```R
# 加载脚本
source("dca_binary_auc_demo.R")

# 载入示例数据
data("df_binary", package = "dcurves")

# 拟合逻辑回归模型得到风险预测
fit_glm <- glm(
  formula = cancer ~ age + famhistory + marker,
  data = df_binary,
  family = binomial()
)
risk_glm <- predict(fit_glm, type = "response")
df_use <- df_binary %>% mutate(risk_glm = risk_glm)

# 计算DCA曲线面积
result <- calc_dca_auc_binary(
  data = df_use,
  outcome = "cancer",
  model_cols = c("risk_glm", "cancerpredmarker"),
  thresholds = seq(0.01, 0.99, by = 0.01)
)

# 查看结果
print(result)
```

### 2. 生存结局示例

```R
# 加载脚本
source("dca_survival_auc_demo.R")

# 载入示例数据
data("df_surv", package = "dcurves")

# 设定评估时间点
t0 <- 1

# 拟合Cox模型得到风险预测
fit_cox <- coxph(
  formula = Surv(ttcancer, cancer) ~ age + marker + famhistory,
  data = df_surv,
  x = TRUE
)

# 计算风险预测值
lp <- predict(fit_cox, type = "lp")
bh <- basehaz(fit_cox, centered = TRUE)
H0_t0 <- approx(
  x = bh$time,
  y = bh$hazard,
  xout = t0,
  method = "linear",
  rule = 2
)$y
S0_t0 <- exp(-H0_t0)
risk_cox <- 1 - (S0_t0 ^ exp(lp))
df_use <- df_surv %>% mutate(risk_cox = risk_cox)

# 计算DCA曲线面积
result <- calc_dca_auc_survival(
  data = df_use,
  time_col = "ttcancer",
  event_col = "cancer",
  model_cols = c("risk_cox", "cancerpredmarker"),
  eval_time = t0,
  thresholds = seq(0.01, 0.99, by = 0.01)
)

# 查看结果
print(result)
```

## 直接运行示例

两个脚本都包含了示例代码，当直接执行时会运行一个最小可复现实例：

```bash
# 运行二分类结局示例
Rscript dca_binary_auc_demo.R

# 运行生存结局示例
Rscript dca_survival_auc_demo.R
```

## 输出解释

输出结果按 `dca_auc_sum_signed` 降序排列，值越大表示模型的净获益越好。

- `dca_auc_pos`：表示在所有阈值范围内，模型净获益为正的区域面积
- `dca_auc_neg`：表示在所有阈值范围内，模型净获益为负的区域面积
- `dca_auc_sum_signed`：有符号总面积，综合反映模型在整个阈值范围内的表现

## 注意事项

1. 输入数据必须包含指定的结局列和模型预测列
2. 模型预测值应该是0-1之间的概率值
3. 对于生存分析，需要指定评估时间点，该时间点应该在数据的时间范围内
4. 脚本会自动过滤不存在的模型列，避免因列不存在导致的错误

## 扩展应用

本工具可以用于：
- 比较不同预测模型的临床决策价值
- 评估模型在不同阈值范围内的表现
- 作为模型选择的辅助指标

## 版本信息

- 日期：2026-04-16
- 依赖包版本：
  - dcurves >= 0.3.0
  - tidyverse >= 1.3.0
  - survival >= 3.2.0