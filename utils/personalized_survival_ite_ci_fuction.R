# ============================================================
# personalized 包：生存结局下逐患者 ITE（benefit score 代理）及 95% 置信区间
# 说明：
# 1) 将核心流程封装为函数，便于复用到真实数据
# 2) 函数保留关键可调参数，含 bootstrap 次数 B
# 3) 演示部分仍可直接运行，验证函数可执行
# ============================================================

# ----------------------------------------------------------
# 0. 加载依赖包
# ----------------------------------------------------------
# 如需自动安装可取消下面注释：
# required_pkgs <- c("personalized", "survival", "dplyr", "tibble", "readr")
# for (pkg in required_pkgs) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg, repos = "https://cloud.r-project.org")
#   }
# }

# 加载本脚本依赖包
library(personalized)
library(survival)
library(dplyr)
library(tibble)
library(readr)

# ----------------------------------------------------------
# 1. 函数：生存结局个体 ITE 及 95% CI 计算
# ----------------------------------------------------------
# 参数说明：
# x: 数值矩阵，行=患者，列=协变量
# time: 生存时间向量
# status: 结局指示（1=事件发生，0=删失）
# trt: 治疗指示（建议 0/1）
# patient_id: 患者唯一 ID；默认自动生成 1:n
# B: bootstrap 重采样次数
# nfolds: fit.subgroup 内部交叉验证折数
# prop_crossfit: 倾向评分函数是否启用 crossfit
# prop_nfolds_crossfit: 倾向评分 crossfit 折数
# prop_cv_metric: 倾向评分内部 cv 指标（如 "auc"）
# seed: 随机种子（保证重现）
# output_file: 若非 NULL，则把结果写出到 csv
# verbose: 是否打印进度信息
calc_survival_ite_ci <- function(
    x,
    time,
    status,
    trt,
    patient_id = NULL,
    B = 200L,
    nfolds = 5L,
    prop_crossfit = TRUE,
    prop_nfolds_crossfit = 5L,
    prop_cv_metric = "auc",
    seed = 20260416L,
    output_file = NULL,
    verbose = TRUE
) {
  # 固定随机种子，确保重复运行结果可复现
  set.seed(seed)

  # 将 x 强制转为矩阵，避免 data.frame 类型导致建模接口不一致
  x <- as.matrix(x)

  # 记录样本量
  n_obs <- nrow(x)

  # 若未提供 patient_id，则自动按行号生成唯一 ID
  if (is.null(patient_id)) {
    patient_id <- seq_len(n_obs)
  }

  # 基础输入检查：维度必须一致
  stopifnot(length(time) == n_obs)
  stopifnot(length(status) == n_obs)
  stopifnot(length(trt) == n_obs)
  stopifnot(length(patient_id) == n_obs)

  # 基础输入检查：B 和 nfolds 至少为 2，避免不稳定或报错
  stopifnot(B >= 2L)
  stopifnot(nfolds >= 2L)

  # status 必须是 0/1，统一转为整数
  status <- as.integer(status)
  stopifnot(all(status %in% c(0L, 1L)))

  # trt 建议为 0/1；若不是则报错提示用户先转换
  trt <- as.integer(trt)
  stopifnot(all(trt %in% c(0L, 1L)))

  # 构建生存结局对象 Surv(time, status)
  y_surv <- survival::Surv(time = time, event = status)

  # 构建 propensity 函数（观察性研究下是关键步骤）
  prop_func <- personalized::create.propensity.function(
    crossfit = prop_crossfit,
    nfolds.crossfit = prop_nfolds_crossfit,
    cv.glmnet.args = list(type.measure = prop_cv_metric, nfolds = nfolds)
  )

  # 在全样本拟合主模型（cox_loss_lasso）
  subgrp_model <- personalized::fit.subgroup(
    x = x,
    y = y_surv,
    trt = trt,
    propensity.func = prop_func,
    loss = "cox_loss_lasso",
    nfolds = nfolds
  )

  # 计算全样本逐患者点估计（benefit score 代理 ITE）
  ite_hat <- as.numeric(
    stats::predict(
      object = subgrp_model,
      newx = x,
      type = "benefit.score"
    )
  )

  # 预分配列表，存储每轮 bootstrap 的 ID-预测值明细
  boot_pred_list <- vector("list", length = B)

  # 循环 bootstrap 次数
  for (b in seq_len(B)) {
    # 每轮有放回重采样索引
    idx_b <- sample.int(n = n_obs, size = n_obs, replace = TRUE)

    # 取本轮重采样数据
    x_b <- x[idx_b, , drop = FALSE]
    y_b <- survival::Surv(time = time[idx_b], event = status[idx_b])
    trt_b <- trt[idx_b]

    # 拟合和预测过程加 tryCatch，防止个别轮次失败影响整体
    boot_tbl_b <- tryCatch(
      {
        # 在 bootstrap 样本上拟合模型
        model_b <- personalized::fit.subgroup(
          x = x_b,
          y = y_b,
          trt = trt_b,
          propensity.func = prop_func,
          loss = "cox_loss_lasso",
          nfolds = nfolds
        )

        # 在“本轮样本自身”上预测，并保留原始患者 ID
        ite_boot_b <- as.numeric(
          stats::predict(
            object = model_b,
            newx = x_b,
            type = "benefit.score"
          )
        )

        # 输出本轮明细：每行对应一次患者命中和其预测值
        tibble::tibble(
          boot_id = b,
          orig_id = patient_id[idx_b],
          ite_boot = ite_boot_b
        )
      },
      error = function(e) {
        # 失败轮次返回空表，后续汇总自动兼容
        tibble::tibble(
          boot_id = integer(),
          orig_id = integer(),
          ite_boot = numeric()
        )
      }
    )

    # 存储本轮结果
    boot_pred_list[[b]] <- boot_tbl_b
  }

  # 合并全部轮次明细
  boot_pred_long <- dplyr::bind_rows(boot_pred_list)

  # 按患者 ID 汇总区间与命中次数
  boot_ci_by_id <- boot_pred_long %>%
    dplyr::group_by(orig_id) %>%
    dplyr::summarise(
      ite_ci_lower = stats::quantile(ite_boot, probs = 0.025, na.rm = TRUE),
      ite_ci_upper = stats::quantile(ite_boot, probs = 0.975, na.rm = TRUE),
      boot_valid_B = dplyr::n(),
      .groups = "drop"
    )

  # 合并点估计和区间估计
  ite_ci_tbl <- tibble::tibble(
    patient_id = patient_id,
    ite_hat = ite_hat
  ) %>%
    dplyr::left_join(
      boot_ci_by_id,
      by = c("patient_id" = "orig_id")
    ) %>%
    dplyr::mutate(
      ite_ci_width = ite_ci_upper - ite_ci_lower
    )

  # 可执行断言：保障输出质量
  stopifnot(nrow(ite_ci_tbl) == n_obs)
  stopifnot(all(ite_ci_tbl$ite_ci_lower <= ite_ci_tbl$ite_ci_upper, na.rm = TRUE))
  stopifnot(all(ite_ci_tbl$boot_valid_B >= 1L))

  # 若传入输出路径，则写出 csv
  if (!is.null(output_file)) {
    readr::write_csv(ite_ci_tbl, file = output_file)
  }

  # 按需打印摘要，便于快速检查
  if (isTRUE(verbose)) {
    cat("\n===== calc_survival_ite_ci 运行完成 =====\n")
    cat("样本量 n = ", n_obs, "\n", sep = "")
    cat("事件率 = ", round(mean(status), 4), "\n", sep = "")
    cat("bootstrap 次数 B = ", B, "\n", sep = "")
    print(ite_ci_tbl %>% dplyr::slice_head(n = 12))
    if (!is.null(output_file)) {
      cat("结果文件：", normalizePath(output_file), "\n", sep = "")
    }
  }

  # 返回列表，含模型对象、明细表和汇总表，便于下游分析
  list(
    model = subgrp_model,
    ite_ci_tbl = ite_ci_tbl,
    boot_pred_long = boot_pred_long
  )
}

# ----------------------------------------------------------
# 2. 演示：模拟生存数据并调用函数（可直接执行）
# ----------------------------------------------------------
# 固定随机种子，保证演示数据可复现
set.seed(20260416)

# 演示数据规模
n_obs <- 300L
n_vars <- 12L

# 生成协变量矩阵
x_demo <- matrix(
  data = rnorm(n_obs * n_vars, sd = 2),
  nrow = n_obs,
  ncol = n_vars
)

# 生成观察性治疗分配
xbetat <- 0.2 + 0.35 * x_demo[, 2] - 0.25 * x_demo[, 7]
trt_prob <- plogis(xbetat)
trt01_demo <- rbinom(n_obs, size = 1, prob = trt_prob)

# 生成个体异质性信号
delta_x <- 0.5 + 0.6 * x_demo[, 1] - 0.5 * x_demo[, 3] + 0.3 * x_demo[, 1] * x_demo[, 4]

# 生成基线风险与治疗效应
lp_base <- 0.2 * x_demo[, 2] - 0.3 * x_demo[, 5] + 0.15 * x_demo[, 9]
trt_effect_lp <- -0.45 * delta_x * trt01_demo

# 生成事件时间与删失时间
event_rate <- exp(lp_base + trt_effect_lp)
event_time <- rexp(n_obs, rate = event_rate)
censor_time <- rexp(n_obs, rate = 0.12)

# 构造观测到的 time/status
time_demo <- pmin(event_time, censor_time)
status_demo <- as.integer(event_time <= censor_time)

# 构造患者 ID
patient_id_demo <- seq_len(n_obs)

# 调用函数进行计算（B 为可调参数）
res_demo <- calc_survival_ite_ci(
  x = x_demo,
  time = time_demo,
  status = status_demo,
  trt = trt01_demo,
  patient_id = patient_id_demo,
  B = 10L,
  nfolds = 5L,
  prop_crossfit = TRUE,
  prop_nfolds_crossfit = 5L,
  prop_cv_metric = "auc",
  seed = 20260416L,
  output_file = "survival_ite_ci_results.csv",
  verbose = TRUE
)
