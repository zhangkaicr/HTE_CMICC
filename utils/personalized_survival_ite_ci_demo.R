# ============================================================
# personalized 包：生存结局下逐患者 ITE（benefit score 代理）及 95% 置信区间演示
# 说明：
# 1) 参考 personalized 官网生存结局流程：Surv(...) + fit.subgroup(..., loss = "cox_loss_lasso")
# 2) 区间估计采用“基于患者 ID 追踪”的 bootstrap 百分位法
# 3) 由于个体真实 ITE 通常不可直接观测，这里以 benefit score 作为个体化治疗效应代理量
# ============================================================

# -------------------------------
# 0. 固定随机种子，保证可复现
# -------------------------------
set.seed(20260416)

# ----------------------------------------------------------
# 1. 加载依赖包
# ----------------------------------------------------------
# 如需自动安装可取消下面注释：
# required_pkgs <- c("personalized", "survival", "dplyr", "purrr", "tibble", "readr")
# for (pkg in required_pkgs) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg, repos = "https://cloud.r-project.org")
#   }
# }

# 加载本脚本依赖
library(personalized)
library(survival)
library(dplyr)
library(purrr)
library(tibble)
library(readr)

# ----------------------------------------------------------
# 2. 模拟生存数据（观察性分配 + 生存时间 + 删失）
# ----------------------------------------------------------
n_obs <- 300L
n_vars <- 12L

# 模拟协变量矩阵（每行一个患者）
x <- matrix(
  data = rnorm(n_obs * n_vars, sd = 2),
  nrow = n_obs,
  ncol = n_vars
)

# 模拟观察性治疗分配概率（0/1 编码，生存示例常用该编码）
xbetat <- 0.2 + 0.35 * x[, 2] - 0.25 * x[, 7]
trt_prob <- plogis(xbetat)
trt01 <- rbinom(n_obs, size = 1, prob = trt_prob)

# 构造“个体化治疗异质性”信号：delta_x > 0 时更偏向治疗组获益
delta_x <- 0.5 + 0.6 * x[, 1] - 0.5 * x[, 3] + 0.3 * x[, 1] * x[, 4]

# 构造基线风险线性预测子（值越大，事件风险越高）
lp_base <- 0.2 * x[, 2] - 0.3 * x[, 5] + 0.15 * x[, 9]

# 构造治疗作用到风险上的线性项：
# 这里设置为：治疗组(trt01=1)在 delta_x 较大时风险下降（生存改善）
trt_effect_lp <- -0.45 * delta_x * trt01

# 事件时间服从指数分布：rate 越高事件越早
event_rate <- exp(lp_base + trt_effect_lp)
event_time <- rexp(n_obs, rate = event_rate)

# 独立删失时间（控制删失比例）
censor_time <- rexp(n_obs, rate = 0.12)

# 观测时间与结局指示
time <- pmin(event_time, censor_time)
status <- as.integer(event_time <= censor_time)

# 生存结局对象
y_surv <- Surv(time = time, event = status)

# 患者唯一 ID（用于 bootstrap 身份追踪）
patient_id <- seq_len(n_obs)

# ----------------------------------------------------------
# 3. 构建 propensity 函数并拟合生存子群模型
# ----------------------------------------------------------
prop_func <- create.propensity.function(
  crossfit = TRUE,
  nfolds.crossfit = 5,
  cv.glmnet.args = list(type.measure = "auc", nfolds = 5)
)

# 生存结局采用 cox 损失；方法用 weighting（默认）
subgrp_model <- fit.subgroup(
  x = x,
  y = y_surv,
  trt = trt01,
  propensity.func = prop_func,
  loss = "cox_loss_lasso",
  nfolds = 5
)

# ----------------------------------------------------------
# 4. 计算逐患者 ITE 点估计（benefit score 代理）
# ----------------------------------------------------------
ite_hat <- as.numeric(
  predict(
    object = subgrp_model,
    newx = x,
    type = "benefit.score"
  )
)

# ----------------------------------------------------------
# 5. bootstrap（ID 追踪）计算每位患者 ITE 95% CI
# ----------------------------------------------------------
# 为了演示速度，先设 B=10；正式分析建议增大到 200 或更高
B <- 10L

# 预分配列表收集每轮 bootstrap 的明细结果
boot_pred_list <- vector("list", length = B)

for (b in seq_len(B)) {
  # 有放回重采样行索引
  idx_b <- sample.int(n = n_obs, size = n_obs, replace = TRUE)

  # 取出本轮 bootstrap 样本
  x_b <- x[idx_b, , drop = FALSE]
  y_b <- Surv(time = time[idx_b], event = status[idx_b])
  trt_b <- trt01[idx_b]

  # 拟合与预测失败时用空表兜底，避免中断整个流程
  boot_tbl_b <- tryCatch(
    {
      model_b <- fit.subgroup(
        x = x_b,
        y = y_b,
        trt = trt_b,
        propensity.func = prop_func,
        loss = "cox_loss_lasso",
        nfolds = 5
      )

      # 在“本轮重采样样本自身”上预测，并绑定原始患者 ID
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
      tibble(
        boot_id = integer(),
        orig_id = integer(),
        ite_boot = numeric()
      )
    }
  )

  # 保存本轮结果
  boot_pred_list[[b]] <- boot_tbl_b
}

# 合并所有轮次明细
boot_pred_long <- bind_rows(boot_pred_list)

# 按患者 ID 汇总区间
boot_ci_by_id <- boot_pred_long %>%
  group_by(orig_id) %>%
  summarise(
    ite_ci_lower = quantile(ite_boot, probs = 0.025, na.rm = TRUE),
    ite_ci_upper = quantile(ite_boot, probs = 0.975, na.rm = TRUE),
    boot_valid_B = n(),
    .groups = "drop"
  )

# 生成最终逐患者结果
ite_ci_tbl <- tibble(
  patient_id = patient_id,
  ite_hat = ite_hat
) %>%
  left_join(
    boot_ci_by_id,
    by = c("patient_id" = "orig_id")
  ) %>%
  mutate(
    ite_ci_width = ite_ci_upper - ite_ci_lower
  )

# ----------------------------------------------------------
# 6. 演示测试代码（可执行断言）
# ----------------------------------------------------------
# 测试1：行数必须等于患者总数
stopifnot(nrow(ite_ci_tbl) == n_obs)

# 测试2：区间下界不大于上界
stopifnot(all(ite_ci_tbl$ite_ci_lower <= ite_ci_tbl$ite_ci_upper, na.rm = TRUE))

# 测试3：每位患者至少出现过 1 次 bootstrap 命中
stopifnot(all(ite_ci_tbl$boot_valid_B >= 1L))

# ----------------------------------------------------------
# 7. 打印摘要并写出结果文件
# ----------------------------------------------------------
cat("\n===== 生存结局：ITE（benefit score 代理）与 95% CI（前 12 行）=====\n")
print(ite_ci_tbl %>% slice_head(n = 12))

cat("\n===== 生存数据概况 =====\n")
cat("样本量 n = ", n_obs, "\n", sep = "")
cat("事件率 = ", round(mean(status), 4), "\n", sep = "")
cat("中位随访时间 = ", round(median(time), 4), "\n", sep = "")

# 输出 CSV
output_file <- "survival_ite_ci_results.csv"
write_csv(ite_ci_tbl, file = output_file)

cat("\n脚本运行完成，结果已保存到：", normalizePath(output_file), "\n", sep = "")
