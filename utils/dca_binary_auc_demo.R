# ----------------------------- #
# 二分类结局 DCA 面积函数脚本
# ----------------------------- #

# 加载 dcurves 包（提供 dca 主函数）
library(dcurves)

# 加载 tidyverse 包（用于数据处理）
library(tidyverse)

# 定义内部工具函数：计算单条曲线的正面积、负面积与有符号总面积
.calc_auc_signed <- function(threshold, net_benefit) {
  # 将阈值转成数值型，避免因子类型导致运算错误
  x <- as.numeric(threshold)
  # 将净获益转成数值型，确保后续可做插值与积分
  y <- as.numeric(net_benefit)
  # 初始化正面积累加器
  area_pos <- 0
  # 初始化负面积累加器（该值通常为负数）
  area_neg <- 0
  # 逐段遍历相邻阈值点，按线段积分
  for (i in seq_len(length(x) - 1)) {
    # 提取左端点横坐标
    x1 <- x[i]
    # 提取右端点横坐标
    x2 <- x[i + 1]
    # 提取左端点纵坐标
    y1 <- y[i]
    # 提取右端点纵坐标
    y2 <- y[i + 1]
    # 计算当前线段宽度
    dx <- x2 - x1
    # 如果宽度非正则跳过，防御异常输入
    if (dx <= 0) next
    # 若整段在 0 线上方，全部计入正面积
    if (y1 >= 0 && y2 >= 0) {
      area_pos <- area_pos + dx * (y1 + y2) / 2
    # 若整段在 0 线下方，全部计入负面积
    } else if (y1 <= 0 && y2 <= 0) {
      area_neg <- area_neg + dx * (y1 + y2) / 2
    # 若线段跨越 0，则先求零点再拆成两段
    } else {
      # 线性插值计算与 y=0 的交点横坐标
      x0 <- x1 - y1 * (x2 - x1) / (y2 - y1)
      # 左子段宽度
      dx1 <- x0 - x1
      # 右子段宽度
      dx2 <- x2 - x0
      # 若左端为正，则左段计正、右段计负
      if (y1 > 0) {
        area_pos <- area_pos + dx1 * (y1 + 0) / 2
        area_neg <- area_neg + dx2 * (0 + y2) / 2
      # 若左端为负，则左段计负、右段计正
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

# 定义对外函数：计算二分类结局的 DCA 曲线面积
calc_dca_auc_binary <- function(
  data,
  outcome,
  model_cols,
  thresholds = seq(0.01, 0.99, by = 0.01)
) {
  # 校验输入数据是否为 data.frame
  if (!is.data.frame(data)) stop("data 必须是 data.frame")
  # 校验结局列是否存在
  if (!outcome %in% names(data)) stop("outcome 列不存在于 data 中")
  # 校验模型列是否提供
  if (length(model_cols) == 0) stop("model_cols 不能为空")
  # 仅保留真实存在的数据列，避免因缺列导致公式报错
  model_cols <- model_cols[model_cols %in% names(data)]
  # 若过滤后模型列为空，给出明确错误
  if (length(model_cols) == 0) stop("model_cols 在 data 中均不存在")
  # 拼接公式右侧字符串
  rhs <- paste(model_cols, collapse = " + ")
  # 构建 dca 公式（例如：cancer ~ risk_glm + cancerpredmarker）
  fml <- as.formula(paste(outcome, "~", rhs))
  # 运行二分类 DCA
  dca_res <- dca(
    formula = fml,
    data = data,
    thresholds = thresholds
  )
  # 计算每条模型曲线的面积，并排除默认策略线
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
  # 载入 dcurves 自带二分类示例数据
  data("df_binary", package = "dcurves")
  # 拟合逻辑回归模型得到示例风险列
  fit_glm <- glm(
    formula = cancer ~ age + famhistory + marker,
    data = df_binary,
    family = binomial()
  )
  # 计算个体预测概率
  risk_glm <- predict(fit_glm, type = "response")
  # 对概率进行边界裁剪
  risk_glm <- pmin(pmax(risk_glm, 0), 1)
  # 组装分析数据
  df_use <- df_binary %>%
    mutate(risk_glm = risk_glm)
  # 自动识别可比较模型列（兼容是否存在 test 列）
  model_cols <- c("risk_glm", "cancerpredmarker", "test")
  model_cols <- model_cols[model_cols %in% names(df_use)]
  # 调用函数并打印结果
  print(
    calc_dca_auc_binary(
      data = df_use,
      outcome = "cancer",
      model_cols = model_cols,
      thresholds = seq(0.01, 0.99, by = 0.01)
    )
  )
}
