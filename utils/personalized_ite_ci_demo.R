# ============================================================
# personalized 包：逐患者 ITE（个体化治疗效应）及其 95% 置信区间演示
# 说明：
# 1) 本脚本严格按 personalized 官网教程的核心分析流程组织：
#    倾向评分函数 -> （可选）增强函数 -> fit.subgroup 拟合 -> 个体层预测
# 2) ITE 置信区间采用非参数 bootstrap 百分位法构建，便于逐患者给区间
# 3) 脚本包含基础可执行测试（assertion），用于快速验证结果结构正确
# ============================================================

# -------------------------------
# 0. 固定随机种子，保证可复现
# -------------------------------
set.seed(20260416)

# ----------------------------------------------------------
# 1. 自动检查并加载依赖包（尽量采用 tidyverse 现代语法流）
# ----------------------------------------------------------
# required_pkgs <- c("personalized", "dplyr", "purrr", "tibble", "readr")

# # 逐个包检查：若未安装则自动安装，然后再加载
# for (pkg in required_pkgs) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg, repos = "https://cloud.r-project.org")
#   }
# }

# 加载本脚本会用到的包
library(personalized)
library(dplyr)
library(purrr)
library(tibble)
library(readr)

# ----------------------------------------------------------
# 2. 按官网示例思路模拟一个“观察性研究”数据集
#    - 处理变量 trt 取值为 {-1, 1}
#    - 结果变量 y 为连续型（高值更好）
# ----------------------------------------------------------
n_obs <- 300L
n_vars <- 12L

# 构造协变量矩阵 x（每行对应一个患者）
x <- matrix(
  data = rnorm(n_obs * n_vars, sd = 2),
  nrow = n_obs,
  ncol = n_vars
)

# 构造非随机治疗分配概率（模拟观察性数据）
xbetat <- 0.3 + 0.4 * x[, 2] - 0.35 * x[, 8]
trt_prob <- exp(xbetat) / (1 + exp(xbetat))
trt01 <- rbinom(n_obs, size = 1, prob = trt_prob)

# 转成 personalized 常用二元处理编码 {-1, 1}
trt <- 2 * trt01 - 1

# 构造真实异质性治疗效应 delta(x)
delta <- 1.8 * (0.4 + x[, 1] - 0.8 * x[, 3] + x[, 2] * x[, 5])

# 构造主效应 g(x)
g_x <- x[, 1] + 0.5 * x[, 4] - 1.3 * x[, 6]^2 + 0.7 * x[, 10]

# 构造连续结局 y（更大更好）
y <- drop(g_x + delta * trt + rnorm(n_obs, sd = 1.8))

# 为每个患者定义唯一 ID，后续 bootstrap 用于身份追踪
patient_id <- seq_len(n_obs)

# ----------------------------------------------------------
# 3. 构建 propensity 与 augmentation 函数（官网推荐流程）
#    注：教程里常用 crossfit=TRUE；演示时折数适中以兼顾速度与稳定性
# ----------------------------------------------------------
prop_func <- create.propensity.function(
  crossfit = TRUE,
  nfolds.crossfit = 5,
  cv.glmnet.args = list(type.measure = "auc", nfolds = 5)
)

aug_func <- create.augmentation.function(
  family = "gaussian",
  crossfit = TRUE,
  nfolds.crossfit = 5,
  cv.glmnet.args = list(type.measure = "mae", nfolds = 5)
)

# ----------------------------------------------------------
# 4. 拟合 subgroup 模型
#    - loss 使用官网常见的连续结局损失 sq_loss_lasso
#    - method 使用 weighting
# ----------------------------------------------------------
subgrp_model <- fit.subgroup(
  x = x,
  y = y,
  trt = trt,
  propensity.func = prop_func,
  augment.func = aug_func,
  loss = "sq_loss_lasso",
  method = "weighting",
  nfolds = 5
)

# ----------------------------------------------------------
# 5. 计算逐患者点估计 ITE
#    这里使用 predict(..., type = "benefit.score")
#    在 personalized 框架中，benefit score 是用于个体化决策的关键量，
#    常作为个体治疗获益（ITE/CATE 的可操作估计）来排序与分层。
# ----------------------------------------------------------
ite_hat <- as.numeric(
  predict(
    object = subgrp_model,
    newx = x,
    type = "benefit.score"
  )
)

# ----------------------------------------------------------
# 6. bootstrap（基于患者 ID 追踪）计算每位患者 ITE 的 95% 置信区间
#    - 每次有放回重采样，并保留该行对应的原始患者 ID（orig_id）
#    - 在重采样数据上拟合模型
#    - 在“本次重采样样本自身”上预测 ITE，并记录 (orig_id, ite_boot)
#    - 最终按 orig_id 合并所有轮次结果，得到每位患者的 ITE 经验分布
#      再计算 2.5% 与 97.5% 分位数作为区间
# ----------------------------------------------------------
B <- 10L

# 用列表收集每轮 bootstrap 的“ID-ITE”明细，最后再拼接
boot_pred_list <- vector("list", length = B)

for (b in seq_len(B)) {
  # 有放回抽样得到 bootstrap 索引
  idx_b <- sample.int(n = n_obs, size = n_obs, replace = TRUE)

  # 取出 bootstrap 样本
  x_b <- x[idx_b, , drop = FALSE]
  y_b <- y[idx_b]
  trt_b <- trt[idx_b]

  # 使用 tryCatch 防止个别重采样失败导致整体中断
  boot_tbl_b <- tryCatch(
    {
      model_b <- fit.subgroup(
        x = x_b,
        y = y_b,
        trt = trt_b,
        propensity.func = prop_func,
        augment.func = aug_func,
        loss = "sq_loss_lasso",
        method = "weighting",
        nfolds = 5
      )

      # 关键：在“本次重采样样本”上预测，并保留每行对应的原始患者 ID
      ite_boot_b <- as.numeric(
        predict(
          object = model_b,
          newx = x_b,
          type = "benefit.score"
        )
      )

      tibble(
        boot_id = b,
        orig_id = patient_id[idx_b],
        ite_boot = ite_boot_b
      )
    },
    error = function(e) {
      # 失败时返回空 tibble，后续 bind_rows 与汇总会自动兼容
      tibble(
        boot_id = integer(),
        orig_id = integer(),
        ite_boot = numeric()
      )
    }
  )

  # 保存本轮 bootstrap 明细结果
  boot_pred_list[[b]] <- boot_tbl_b
}

# ----------------------------------------------------------
# 7. 汇总逐患者 ITE 置信区间结果
# ----------------------------------------------------------
# 将所有 bootstrap 明细合并成一张长表
boot_pred_long <- bind_rows(boot_pred_list)

# 按患者 ID 聚合，得到每位患者的区间与有效样本量
boot_ci_by_id <- boot_pred_long %>%
  group_by(orig_id) %>%
  summarise(
    ite_ci_lower = quantile(ite_boot, probs = 0.025, na.rm = TRUE),
    ite_ci_upper = quantile(ite_boot, probs = 0.975, na.rm = TRUE),
    # 统计该患者在所有 bootstrap 样本中被抽到的总次数
    boot_valid_B = n(),
    .groups = "drop"
  )

ite_ci_tbl <- tibble(
  patient_id = patient_id,
  ite_hat = ite_hat
) %>%
  left_join(
    boot_ci_by_id,
    by = c("patient_id" = "orig_id")
) %>%
  mutate(
    # 区间宽度可用于衡量每位患者估计不确定性
    ite_ci_width = ite_ci_upper - ite_ci_lower
  )

# ----------------------------------------------------------
# 8. 演示测试代码（可执行断言）
# ----------------------------------------------------------
# 测试1：结果行数应与患者数一致
stopifnot(nrow(ite_ci_tbl) == n_obs)

# 测试2：每位患者的下界应不大于上界（允许极端情况下同值）
stopifnot(all(ite_ci_tbl$ite_ci_lower <= ite_ci_tbl$ite_ci_upper, na.rm = TRUE))

# 测试3：每位患者至少应在 bootstrap 汇总中出现 1 次（基于 ID 合并）
stopifnot(all(ite_ci_tbl$boot_valid_B >= 1L))

# ----------------------------------------------------------
# 9. 输出结果到控制台和 CSV 文件
# ----------------------------------------------------------
cat("\n===== ITE 与 95% 置信区间（前 12 行）=====\n")
print(ite_ci_tbl %>% slice_head(n = 12))

# 将完整结果写出，便于下游分析或可视化
output_file <- "ite_ci_results.csv"
write_csv(ite_ci_tbl, file = output_file)

cat("\n脚本运行完成，结果已保存到：", normalizePath(output_file), "\n", sep = "")
