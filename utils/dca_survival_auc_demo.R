# ----------------------------- #
# 生存结局 DCA 面积函数脚本
# ----------------------------- #

# 加载 survival 包（提供 Surv 与 Cox 拟合工具）
library(survival)

# 加载 dcurves 包（提供 dca 主函数）
library(dcurves)

# 加载 tidyverse 包（用于数据处理）
library(tidyverse)

# 定义内部工具函数：计算单条曲线的正面积、负面积与有符号总面积
.calc_auc_signed <- function(threshold, net_benefit) {
  # 将阈值转成数值型，避免类型问题
  x <- as.numeric(threshold)
  # 将净获益转成数值型，确保可积分
  y <- as.numeric(net_benefit)
  # 仅保留成对非缺失值，避免 NA 触发逻辑判断报错
  keep <- !(is.na(x) | is.na(y))
  x <- x[keep]
  y <- y[keep]
  # 按阈值从小到大排序，保证积分方向正确
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  # 若有效点少于 2 个，无法形成线段，直接返回 0 面积
  if (length(x) < 2) {
    return(
      tibble(
        auc_pos = 0,
        auc_neg = 0,
        auc_signed = 0
      )
    )
  }
  # 初始化正面积累加器
  area_pos <- 0
  # 初始化负面积累加器
  area_neg <- 0
  # 逐线段积分
  for (i in seq_len(length(x) - 1)) {
    # 左端点横坐标
    x1 <- x[i]
    # 右端点横坐标
    x2 <- x[i + 1]
    # 左端点纵坐标
    y1 <- y[i]
    # 右端点纵坐标
    y2 <- y[i + 1]
    # 当前线段宽度
    dx <- x2 - x1
    # 若宽度异常则跳过
    if (dx <= 0) next
    # 同号在线上方，计入正面积
    if (y1 >= 0 && y2 >= 0) {
      area_pos <- area_pos + dx * (y1 + y2) / 2
    # 同号在线下方，计入负面积
    } else if (y1 <= 0 && y2 <= 0) {
      area_neg <- area_neg + dx * (y1 + y2) / 2
    # 跨 0 线，切分两段积分
    } else {
      # 线性插值零点
      x0 <- x1 - y1 * (x2 - x1) / (y2 - y1)
      # 左段宽度
      dx1 <- x0 - x1
      # 右段宽度
      dx2 <- x2 - x0
      # 按左端符号分配正负面积
      if (y1 > 0) {
        area_pos <- area_pos + dx1 * (y1 + 0) / 2
        area_neg <- area_neg + dx2 * (0 + y2) / 2
      } else {
        area_neg <- area_neg + dx1 * (y1 + 0) / 2
        area_pos <- area_pos + dx2 * (0 + y2) / 2
      }
    }
  }
  # 返回三种面积指标
  tibble(
    auc_pos = area_pos,
    auc_neg = area_neg,
    auc_signed = area_pos + area_neg
  )
}

# 定义对外函数：计算生存结局的 DCA 曲线面积
calc_dca_auc_survival <- function(
  data,
  time_col,
  event_col,
  model_cols,
  eval_time,
  thresholds = seq(0.01, 0.99, by = 0.01)
) {
  # 校验输入是否为 data.frame
  if (!is.data.frame(data)) stop("data 必须是 data.frame")
  # 校验时间列是否存在
  if (!time_col %in% names(data)) stop("time_col 不存在于 data 中")
  # 校验事件列是否存在
  if (!event_col %in% names(data)) stop("event_col 不存在于 data 中")
  # 校验模型列是否提供
  if (length(model_cols) == 0) stop("model_cols 不能为空")
  # 仅保留实际存在的模型列
  model_cols <- model_cols[model_cols %in% names(data)]
  # 若过滤后为空则报错
  if (length(model_cols) == 0) stop("model_cols 在 data 中均不存在")
  # 构建公式左侧：Surv(time_col, event_col)
  lhs <- paste0("Surv(", time_col, ", ", event_col, ")")
  # 构建公式右侧：多个模型列
  rhs <- paste(model_cols, collapse = " + ")
  # 拼成完整生存 DCA 公式
  fml <- as.formula(paste(lhs, "~", rhs))
  # 运行生存 DCA
  dca_res <- dca(
    formula = fml,
    data = data,
    time = eval_time,
    thresholds = thresholds
  )
  # 计算每条模型曲线面积，排除默认策略线
  auc_df <- dca_res %>%
    as_tibble() %>%
    filter(!str_detect(tolower(label), "all|none")) %>%
    group_by(label) %>%
    arrange(threshold, .by_group = TRUE) %>%
    reframe(.calc_auc_signed(threshold, net_benefit)) %>%
    ungroup() %>%
    mutate(
      dca_auc_pos = auc_pos,
      dca_auc_neg = auc_neg,
      dca_auc_sum_signed = auc_signed
    ) %>%
    select(label, dca_auc_pos, dca_auc_neg, dca_auc_sum_signed) %>%
    arrange(desc(dca_auc_sum_signed))
  # 返回结果数据框
  auc_df
}

# 当脚本被 Rscript 直接执行时，运行一个最小可复现实例
if (sys.nframe() == 0) {
  # 载入 dcurves 自带生存示例数据
  data("df_surv", package = "dcurves")
  # 设定评估时间点（与 ttcancer 单位一致）
  t0 <- 1
  # 拟合 Cox 模型
  fit_cox <- coxph(
    formula = Surv(ttcancer, cancer) ~ age + marker + famhistory,
    data = df_surv,
    x = TRUE
  )
  # 计算线性预测值
  lp <- predict(fit_cox, type = "lp")
  # 提取基线累计风险
  bh <- basehaz(fit_cox, centered = TRUE)
  # 插值获得 t0 的基线累计风险
  H0_t0 <- approx(
    x = bh$time,
    y = bh$hazard,
    xout = t0,
    method = "linear",
    rule = 2
  )$y
  # 计算 t0 的基线生存概率
  S0_t0 <- exp(-H0_t0)
  # 计算个体绝对风险
  risk_cox <- 1 - (S0_t0 ^ exp(lp))
  # 对风险做边界裁剪
  risk_cox <- pmin(pmax(risk_cox, 0), 1)
  # 组装分析数据
  df_use <- df_surv %>%
    mutate(risk_cox = risk_cox)
  # 自动识别可比较模型列
  model_cols <- c("risk_cox", "cancerpredmarker", "test")
  model_cols <- model_cols[model_cols %in% names(df_use)]
  # 调用函数并打印结果
  print(
    calc_dca_auc_survival(
      data = df_use,
      time_col = "ttcancer",
      event_col = "cancer",
      model_cols = model_cols,
      eval_time = t0,
      thresholds = seq(0.01, 0.99, by = 0.01)
    )
  )
}
