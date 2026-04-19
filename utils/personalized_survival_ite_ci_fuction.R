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
library(personalized)   # 用于亚组识别与个性化治疗效应估计
library(here)           # 便于相对路径管理
library(tidyverse)      # 数据清洗与可视化
library(openxlsx)       # 读写 Excel 文件
library(survival)       # 生存分析基础函数
library(fastDummies)    # 快速生成哑变量
library(MatchIt)
library("MatchIt")
library(tidyverse)
library(here)
library(openxlsx)
library(cobalt)
library(compareGroups)
library(survival)
library(broom)
library(personalized)   # 用于亚组识别与个性化治疗效应估计
library(here)           # 便于相对路径管理
library(tidyverse)      # 数据清洗与可视化
library(openxlsx)       # 读写 Excel 文件
library(survival)       # 生存分析基础函数
library(fastDummies)    # 快速生成哑变量
library(tidyverse)  # 用于数据处理和可视化的综合工具包
library(PSweight)  # 实现倾向性评分匹配的专用包
library(openxlsx)   # 用于读写Excel格式文件
library(here)       # 用于项目路径的智能管理
library(survival)   # 提供生存分析的核心功能
library(survminer)  # 用于生存分析结果的可视化
library(ggsci)      # 提供科学期刊标准配色方案
library(tidysmd)    # 用于计算标准化均值差异
library(propensity) # 提供倾向性评分分析工具
library(halfmoon)   # 用于倾向性评分的诊断分析     # 用于评估协变量平衡性
library(mlr3verse)  # 机器学习框架包
library(predRupdate)  # 用于预测模型更新
library(broom)
library(jskm)
library(survey)
library(gtsummary)
library(tidyverse)
library(survival)
library(survey)
library(rempsyc)
library(remotes)
library(HTEPredictionMetrics)
library(tidyverse)                                       # 数据清洗与可视化核心包（dplyr/ggplot2/readr 等）
library(grf)                                             # 广义随机森林，用于因果生存分析
library(survival)                                        # 生存分析基础函数（Surv, survfit 等）
library(fastDummies)                                     # 快速将分类变量转为哑变量
library(here)                                            # 以项目根目录为基准的稳健路径解析
library(visdat)                                          # 可视化缺失值分布
library(ggpubr)
library(survminer)
library(personalized)   # 用于亚组识别与个性化治疗效应估计
library(here)           # 便于相对路径管理
library(tidyverse)      # 数据清洗与可视化
library(openxlsx)       # 读写 Excel 文件
library(survival)       # 生存分析基础函数
library(fastDummies)  
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

final_data <- read.xlsx(here("data", "final_data.xlsx"))
head(final_data)

# 2.定义需要的变量---------------------------------
x_matched <- final_data %>%
  select(Tumor_border_status,CD4_CD8,Peritumoral_enhancement, 
  CRP,CPgrade,Radiologic_morphology,NLR,CAperihepatic_99,Neutro,Cr,
  PDL1,DBIL,BUN,ALB,TumorDiameterMax,
  Vascular_invasion,Tumor_deposit,  
  )

# 对所有因子/字符列做哑变量，删除原列，并去掉首个水平以防共线
x_dummy_matched <- fastDummies::dummy_cols(
  x_matched,
  remove_selected_columns = TRUE,
  remove_first_dummy = TRUE
)

# 将哑变量数据框转换为矩阵，便于后续算法使用
x_demo <- as.matrix(x_dummy_matched)

# 提取生存时间
time_demo <- final_data$OS %>% as.numeric()

# 提取生存状态（1 = 死亡，0 = 存活/删失），并转为整数型
status_demo <- final_data$OS_CNSR %>% as.integer()

trt01_demo <- final_data$Trt_ITT %>% as.integer()

patient_id_demo <- final_data$ID



# 调用函数进行计算（B 为可调参数）
res_demo <- calc_survival_ite_ci(
  x = x_demo,
  time = time_demo,
  status = status_demo,
  trt = trt01_demo,
  patient_id = patient_id_demo,
  B = 5000L,
  nfolds = 10L,
  prop_crossfit = TRUE,
  prop_nfolds_crossfit = 5L,
  prop_cv_metric = "auc",
  seed = 20260416L,
  output_file = "survival_ite_ci_results.xlsx",
  verbose = TRUE
)


############绘图函数###########################

library(here)
here()
  df2<- read_csv(here("results", "survival_ite_ci_test_results.csv"))
  df2$ite_hat <- as.numeric(df2$ite_hat)
  df2$ite_ci_lower <- as.numeric(df2$ite_ci_lower)
  df2$ite_ci_upper <- as.numeric(df2$ite_ci_upper)
  head(df2)
  
# 准备数据：创建Rank列
# 按ite_hat从高到低排序，创建Rank列
df2_ranked <- df2 %>%
  dplyr::arrange(desc(ite_hat)) %>%                        # 按ite_hat降序排序（从高到低）  
  dplyr::mutate(Rank = dplyr::row_number())                      # 创建Rank列（1到n，1为最高分）

# 查看数据摘要
cat("\n=== 数据摘要 ===\n")
cat(sprintf("总患者数: %d\n", nrow(df2_ranked)))
cat(sprintf("ite_hat范围: %.6f 到 %.6f\n", 
            min(df2_ranked$ite_hat, na.rm = TRUE),
            max(df2_ranked$ite_hat, na.rm = TRUE)))
cat(sprintf("ite_hat均值: %.6f\n", 
            mean(df2_ranked$ite_hat, na.rm = TRUE)))
cat(sprintf("ite_hat中位数: %.6f\n", 
            median(df2_ranked$ite_hat, na.rm = TRUE)))

# 绘制排名图：使用颜色渐变
# 创建散点图，横轴为Rank，纵轴为ite_hat，颜色按ite_hat渐变
ite_hat_rank_plot <- ggplot2::ggplot(df2_ranked, 
                                      ggplot2::aes(x = Rank, y = ite_hat)) +
  # 绘制置信区间带
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ite_ci_lower, ymax = ite_ci_upper),
    fill = "lightblue",
    alpha = 0.3
  ) +
  # 绘制散点图
  ggplot2::geom_point(
    ggplot2::aes(color = ite_hat),
    size = 3, 
    alpha = 0.7
  ) +                  # 点的大小和透明度
  # 设置颜色渐变方案：从#666666be（灰色，低分）到#2d8737af（绿色，高分）均匀变化
  ggplot2::scale_color_gradient(
    low = "#643292",                                            # 低分颜色：灰色（带透明度）
    high = "#04f87e",                                           # 高分颜色：绿色（带透明度）
    name = "Benefit\nScore",                                      # 图例标题
    guide = ggplot2::guide_colorbar(                              # 使用颜色条图例
      title.position = "top",                                     # 标题位置在顶部
      barwidth = 2,                                               # 颜色条宽度（垂直图例时控制宽度）
      barheight = 20,                                             # 颜色条高度（垂直图例时控制高度，增大以拉长）
      title.hjust = 0.5                                           # 标题水平居中
    )
  ) +
  # 设置坐标轴标签和标题
  ggplot2::labs(
    title = "Predicted individualized treatment effect for each patient",  # 主标题
    x = "Patients ranked by predicted individualized treatment effect",    # x轴标签
    y = "Predicted individualized treatment effect"                      # y轴标签
  ) +
  # 设置主题样式，与之前图表一致
  ggplot2::theme(
    axis.title.x = element_text(size = 24, face = "bold", color = "#000000"),  # x轴标题：24号字体，加粗，黑色
    axis.title.y = element_text(size = 24, face = "bold", color = "#000000"),  # y轴标题：24号字体，加粗，黑色
    axis.text.x = element_text(size = 22, color = "#000000"),                  # x轴刻度标签：22号字体，黑色
    axis.text.y = element_text(size = 22, color = "#000000"),                  # y轴刻度标签：22号字体，黑色
    legend.title = element_text(size = 22, face = "bold", color = "#000000"),  # 图例标题：22号字体，加粗，黑色
    legend.text = element_text(size = 20, color = "#000000"),                  # 图例文字：20号字体，黑色
    legend.position = "right",                                                 # 图例位置：右侧
    plot.title = element_text(size = 24, face = "bold", color = "#000000",     # 主标题：24号字体，加粗，黑色
                              hjust = 0.5),                                     # 标题居中
    plot.subtitle = element_text(size = 22, color = "#000000"),                # 副标题：22号字体，黑色
    # 添加边框线
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.8), # 黑色边框，线宽0.8
    # 设置背景颜色为白色
    panel.background = element_rect(fill = "#FFFFFF"),                         # 白色背景
    # 设置网格线颜色为白色（隐藏网格线）
    panel.grid.major = element_line(color = "#FFFFFF"),                         # 主网格线为白色（不可见）
    panel.grid.minor = element_line(color = "#FFFFFF")                           # 次网格线为白色（不可见）
  )

# 显示图形
print(ite_hat_rank_plot)

# 直接保存到当前目录
cat("Saving PDF to ite_hat_ranked_plot.pdf...\n")
pdf("ite_hat_ranked_plot_test.pdf", width = 12, height = 8)
print(ite_hat_rank_plot)
dev.off()
cat("PDF saved successfully to current directory\n")


