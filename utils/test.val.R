############################################################
# personalized 包：
# 基于训练集拟合生存结局 ITE（以 benefit score 作为逐患者代理值），
# 并在固定测试集上计算点估计与 bootstrap 百分位置信区间
#
# 说明：
# 1. 训练集用于拟合 personalized::fit.subgroup() 模型
# 2. 测试集始终保持固定，只接受每一轮 bootstrap 模型的预测
# 3. 区间估计来自“训练集有放回重采样 -> 重拟合 -> 固定测试集预测”
# 4. 输出的 ite_hat 本质上是 personalized 的 benefit.score，
#    与你原脚本保持同一口径，作为逐患者 ITE 的代理量使用
############################################################

# ----------------------------------------------------------
# 0. 加载依赖包
# ----------------------------------------------------------
# 加载 personalized 包，用于拟合个体化治疗规则/benefit score
library(personalized)

# 加载 survival 包，用于构造生存结局对象 Surv()
library(survival)

# 加载 tidyverse 包，便于使用现代数据处理语法
library(tidyverse)

# 加载 openxlsx 包，用于读取和写出 Excel 文件
library(openxlsx)

# ----------------------------------------------------------
# 1. 工具函数：检查列是否存在
# ----------------------------------------------------------
# 定义一个小工具函数，用于统一检查数据中是否包含指定字段
assert_required_columns <- function(data, required_cols, data_name) {
  # 找出当前数据缺失的列名
  missing_cols <- setdiff(required_cols, names(data))

  # 如果存在缺失列，则直接报错并给出明确提示
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        data_name,
        " 缺少以下字段：",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

# ----------------------------------------------------------
# 2. 工具函数：检查关键字段是否缺失
# ----------------------------------------------------------
# 定义一个工具函数，用于统计并阻止关键字段缺失值进入建模流程
assert_no_missing <- function(data, check_cols, data_name) {
  # 逐列统计缺失值个数
  na_count <- data %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(check_cols), ~ sum(is.na(.x)))) %>%
    as.list()

  # 挑出缺失值个数大于 0 的列
  bad_cols <- names(na_count)[unlist(na_count) > 0]

  # 若存在缺失值，则报错并提示具体列及缺失数
  if (length(bad_cols) > 0) {
    bad_msg <- purrr::map_chr(
      bad_cols,
      ~ paste0(.x, "=", na_count[[.x]])
    ) %>%
      paste(collapse = ", ")

    stop(
      paste0(data_name, " 中关键字段存在缺失值：", bad_msg),
      call. = FALSE
    )
  }
}

# ----------------------------------------------------------
# 3. 工具函数：统一训练集与测试集设计矩阵
# ----------------------------------------------------------
# 该函数把训练集和测试集的协变量先合并后一起做 model.matrix()，
# 这样可以保证分类变量哑变量编码完全一致，避免两套数据列不对齐
build_aligned_design_matrices <- function(train_df, test_df, feature_cols) {
  # 提取训练集协变量原始数据
  train_x_raw <- train_df %>%
    dplyr::select(dplyr::all_of(feature_cols))

  # 提取测试集协变量原始数据
  test_x_raw <- test_df %>%
    dplyr::select(dplyr::all_of(feature_cols))

  # 记录训练集样本量，后面切分矩阵时要用
  n_train <- nrow(train_x_raw)

  # 把训练集和测试集协变量纵向拼接，确保编码规则一致
  combined_x_raw <- dplyr::bind_rows(train_x_raw, test_x_raw) %>%
    # 把字符型变量显式转为因子，便于 model.matrix 自动展开哑变量
    dplyr::mutate(dplyr::across(where(is.character), as.factor)) %>%
    # 把逻辑型变量也转为因子，避免后续隐式类型转换带来不一致
    dplyr::mutate(dplyr::across(where(is.logical), as.factor))

  # 用 model.matrix 生成无截距设计矩阵
  combined_x_mat <- stats::model.matrix(
    object = ~ . - 1,
    data = as.data.frame(combined_x_raw)
  )

  # 按原始训练集行数切回训练集矩阵
  train_x_mat <- combined_x_mat[seq_len(n_train), , drop = FALSE]

  # 按原始训练集行数之后的部分切回测试集矩阵
  test_x_mat <- combined_x_mat[(n_train + 1):nrow(combined_x_mat), , drop = FALSE]

  # 返回对齐后的训练集矩阵、测试集矩阵以及最终列名
  list(
    train_x_mat = train_x_mat,
    test_x_mat = test_x_mat,
    design_colnames = colnames(combined_x_mat)
  )
}

# ----------------------------------------------------------
# 4. 工具函数：把点估计与 bootstrap 预测整理成标准结果表
# ----------------------------------------------------------
# 该函数用于将任意一个数据集（训练集或测试集）的点估计与
# bootstrap 预测矩阵统一整理为结果表，避免同样的代码写两遍
build_ite_result_table <- function(
    data_use,
    id_col,
    time_col,
    status_col,
    trt_col,
    ite_hat,
    boot_pred_mat,
    dataset_label
) {
  # 计算每位患者在 bootstrap 中成功获得预测的次数
  boot_valid_B <- rowSums(!is.na(boot_pred_mat))

  # 计算每位患者的区间下限
  ite_ci_lower <- apply(
    X = boot_pred_mat,
    MARGIN = 1,
    FUN = function(x) {
      if (all(is.na(x))) {
        return(NA_real_)
      }
      as.numeric(stats::quantile(x, probs = 0.025, na.rm = TRUE, names = FALSE))
    }
  )

  # 计算每位患者的区间上限
  ite_ci_upper <- apply(
    X = boot_pred_mat,
    MARGIN = 1,
    FUN = function(x) {
      if (all(is.na(x))) {
        return(NA_real_)
      }
      as.numeric(stats::quantile(x, probs = 0.975, na.rm = TRUE, names = FALSE))
    }
  )

  # 组装结果表并追加常用解释字段
  tibble::tibble(
    dataset = dataset_label,
    patient_id = data_use[[id_col]],
    observed_time = as.numeric(data_use[[time_col]]),
    observed_status = as.integer(data_use[[status_col]]),
    observed_trt = as.integer(data_use[[trt_col]]),
    ite_hat = ite_hat,
    ite_ci_lower = ite_ci_lower,
    ite_ci_upper = ite_ci_upper,
    boot_valid_B = boot_valid_B
  ) %>%
    dplyr::mutate(
      ite_ci_width = ite_ci_upper - ite_ci_lower,
      recommended_trt = dplyr::if_else(ite_hat > 0, 1L, 0L),
      ci_excludes_zero = dplyr::if_else(
        !is.na(ite_ci_lower) & !is.na(ite_ci_upper) &
          (ite_ci_lower > 0 | ite_ci_upper < 0),
        TRUE,
        FALSE,
        missing = FALSE
      ),
      ite_rank_desc = dplyr::min_rank(dplyr::desc(ite_hat))
    )
}

# ----------------------------------------------------------
# 5. 主函数：一次 bootstrap 同时输出训练集与测试集结果
# ----------------------------------------------------------
# 该函数是本脚本的核心：
# 1. 用完整训练集拟合主模型
# 2. 在原始训练集和固定测试集上分别得到点估计
# 3. 每次 bootstrap 都在训练集上有放回重采样重建模型
# 4. 每轮模型同时预测“原始训练集”和“固定测试集”
# 5. 最终一次性返回训练集和测试集两套 ITE 点估计与可信区间
calc_train_test_survival_ite_ci <- function(
    train_data,
    test_data,
    feature_cols,
    time_col = "OS",
    status_col = "OS_CNSR",
    trt_col = "Trt_ITT",
    id_col = "ID",
    B = 500L,
    nfolds = 5L,
    prop_crossfit = TRUE,
    prop_nfolds_crossfit = 5L,
    prop_cv_metric = "auc",
    loss = "cox_loss_lasso",
    seed = 20260418L,
    output_train_csv = NULL,
    output_test_csv = NULL,
    output_xlsx = NULL,
    verbose = TRUE
) {
  # 固定随机种子，确保结果可复现
  set.seed(seed)

  # 把 bootstrap 次数强制转成整数，避免用户传入数值型小数
  B <- as.integer(B)

  # 把交叉验证折数强制转成整数
  nfolds <- as.integer(nfolds)

  # 把倾向评分 crossfit 折数强制转成整数
  prop_nfolds_crossfit <- as.integer(prop_nfolds_crossfit)

  # 对 bootstrap 次数做基本合法性检查
  if (B < 2L) {
    stop("B 至少需要 >= 2。", call. = FALSE)
  }

  # 对建模折数做基本合法性检查
  if (nfolds < 2L) {
    stop("nfolds 至少需要 >= 2。", call. = FALSE)
  }

  # 整理训练集必须包含的字段
  train_required_cols <- c(id_col, time_col, status_col, trt_col, feature_cols)

  # 整理测试集必须包含的字段
  test_required_cols <- c(id_col, time_col, status_col, trt_col, feature_cols)

  # 检查训练集字段是否齐全
  assert_required_columns(
    data = train_data,
    required_cols = train_required_cols,
    data_name = "train_data"
  )

  # 检查测试集字段是否齐全
  assert_required_columns(
    data = test_data,
    required_cols = test_required_cols,
    data_name = "test_data"
  )

  # 检查训练集关键字段是否有缺失
  assert_no_missing(
    data = train_data,
    check_cols = train_required_cols,
    data_name = "train_data"
  )

  # 检查测试集关键字段是否有缺失
  assert_no_missing(
    data = test_data,
    check_cols = test_required_cols,
    data_name = "test_data"
  )

  # 抽取训练集本次需要使用的字段，减少后续对象体积
  train_use <- train_data %>%
    dplyr::select(dplyr::all_of(train_required_cols))

  # 抽取测试集本次需要使用的字段，减少后续对象体积
  test_use <- test_data %>%
    dplyr::select(dplyr::all_of(test_required_cols))

  # 把训练集生存时间转成数值型
  time_train <- as.numeric(train_use[[time_col]])

  # 把训练集事件指标转成整数型
  status_train <- as.integer(train_use[[status_col]])

  # 把训练集治疗指示转成整数型
  trt_train <- as.integer(train_use[[trt_col]])

  # 把测试集生存时间转成数值型，主要用于结果回填和核对
  time_test <- as.numeric(test_use[[time_col]])

  # 把测试集事件指标转成整数型，主要用于结果回填和核对
  status_test <- as.integer(test_use[[status_col]])

  # 把测试集治疗指示转成整数型，主要用于结果回填和核对
  trt_test <- as.integer(test_use[[trt_col]])

  # 检查训练集事件指标是否为 0/1
  if (!all(status_train %in% c(0L, 1L))) {
    stop("训练集结局状态必须编码为 0/1。", call. = FALSE)
  }

  # 检查测试集事件指标是否为 0/1
  if (!all(status_test %in% c(0L, 1L))) {
    stop("测试集结局状态必须编码为 0/1。", call. = FALSE)
  }

  # 检查训练集治疗变量是否为 0/1
  if (!all(trt_train %in% c(0L, 1L))) {
    stop("训练集治疗变量必须编码为 0/1。", call. = FALSE)
  }

  # 检查测试集治疗变量是否为 0/1
  if (!all(trt_test %in% c(0L, 1L))) {
    stop("测试集治疗变量必须编码为 0/1。", call. = FALSE)
  }

  # 统一生成训练集和测试集的设计矩阵
  design_obj <- build_aligned_design_matrices(
    train_df = train_use,
    test_df = test_use,
    feature_cols = feature_cols
  )

  # 取出训练集设计矩阵
  x_train <- design_obj$train_x_mat

  # 取出测试集设计矩阵
  x_test <- design_obj$test_x_mat

  # 记录训练集样本量
  n_train <- nrow(x_train)

  # 记录测试集样本量
  n_test <- nrow(x_test)

  # 构造训练集生存对象
  y_train <- survival::Surv(time = time_train, event = status_train)

  # 创建倾向评分函数
  prop_func <- personalized::create.propensity.function(
    crossfit = prop_crossfit,
    nfolds.crossfit = prop_nfolds_crossfit,
    cv.glmnet.args = list(
      type.measure = prop_cv_metric,
      nfolds = nfolds
    )
  )

  # 封装一次模型拟合，便于主拟合和 bootstrap 重复调用
  fit_subgroup_once <- function(x, y, trt) {
    personalized::fit.subgroup(
      x = x,
      y = y,
      trt = trt,
      propensity.func = prop_func,
      loss = loss,
      nfolds = nfolds
    )
  }

  # 如需打印进度，则先输出主模型拟合提示
  if (isTRUE(verbose)) {
    cat("\n===== 开始拟合主模型（训练集） =====\n")
  }

  # 在完整训练集上拟合主模型
  base_model <- fit_subgroup_once(
    x = x_train,
    y = y_train,
    trt = trt_train
  )

  # 在原始训练集上计算点估计
  ite_hat_train <- as.numeric(
    stats::predict(
      object = base_model,
      newx = x_train,
      type = "benefit.score"
    )
  )

  # 在固定测试集上计算点估计
  ite_hat_test <- as.numeric(
    stats::predict(
      object = base_model,
      newx = x_test,
      type = "benefit.score"
    )
  )

  # 预分配训练集 bootstrap 预测矩阵：
  # 行对应原始训练集患者，列对应 bootstrap 轮次
  boot_pred_train_mat <- matrix(
    data = NA_real_,
    nrow = n_train,
    ncol = B
  )

  # 预分配测试集 bootstrap 预测矩阵：
  # 行对应固定测试集患者，列对应 bootstrap 轮次
  boot_pred_test_mat <- matrix(
    data = NA_real_,
    nrow = n_test,
    ncol = B
  )

  # 若需打印进度，则输出 bootstrap 开始提示
  if (isTRUE(verbose)) {
    cat("\n===== 开始 bootstrap（训练集重采样，同时预测训练集与测试集） =====\n")
  }

  # 循环进行 B 次 bootstrap
  for (b in seq_len(B)) {
    # 在训练集索引上执行有放回重采样
    boot_idx <- sample.int(
      n = n_train,
      size = n_train,
      replace = TRUE
    )

    # 基于本轮重采样索引构造 bootstrap 生存对象
    y_boot <- survival::Surv(
      time = time_train[boot_idx],
      event = status_train[boot_idx]
    )

    # 尝试拟合本轮 bootstrap 模型；
    # 若某一轮因数值问题失败，则保留整列 NA 并继续
    boot_model <- tryCatch(
      {
        fit_subgroup_once(
          x = x_train[boot_idx, , drop = FALSE],
          y = y_boot,
          trt = trt_train[boot_idx]
        )
      },
      error = function(e) {
        NULL
      }
    )

    # 若本轮模型拟合成功，则同时对原始训练集和固定测试集进行预测
    if (!is.null(boot_model)) {
      boot_pred_train_mat[, b] <- tryCatch(
        {
          as.numeric(
            stats::predict(
              object = boot_model,
              newx = x_train,
              type = "benefit.score"
            )
          )
        },
        error = function(e) {
          rep(NA_real_, n_train)
        }
      )

      boot_pred_test_mat[, b] <- tryCatch(
        {
          as.numeric(
            stats::predict(
              object = boot_model,
              newx = x_test,
              type = "benefit.score"
            )
          )
        },
        error = function(e) {
          rep(NA_real_, n_test)
        }
      )
    }

    # 按设定频率打印进度，便于长任务监控
    if (isTRUE(verbose) && (b %% 10L == 0L || b == 1L || b == B)) {
      cat("bootstrap 进度：", b, "/", B, "\n", sep = "")
    }
  }

  # 把训练集点估计与 bootstrap 预测整理为标准结果表
  train_results_tbl <- build_ite_result_table(
    data_use = train_use,
    id_col = id_col,
    time_col = time_col,
    status_col = status_col,
    trt_col = trt_col,
    ite_hat = ite_hat_train,
    boot_pred_mat = boot_pred_train_mat,
    dataset_label = "train"
  )

  # 把测试集点估计与 bootstrap 预测整理为标准结果表
  test_results_tbl <- build_ite_result_table(
    data_use = test_use,
    id_col = id_col,
    time_col = time_col,
    status_col = status_col,
    trt_col = trt_col,
    ite_hat = ite_hat_test,
    boot_pred_mat = boot_pred_test_mat,
    dataset_label = "test"
  )

  # 如果给定训练集 CSV 输出路径，则先确保目录存在再写出
  if (!is.null(output_train_csv)) {
    dir.create(dirname(output_train_csv), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(train_results_tbl, file = output_train_csv)
  }

  # 如果给定测试集 CSV 输出路径，则先确保目录存在再写出
  if (!is.null(output_test_csv)) {
    dir.create(dirname(output_test_csv), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(test_results_tbl, file = output_test_csv)
  }

  # 如果给定 Excel 输出路径，则按两个 sheet 一次性写出
  if (!is.null(output_xlsx)) {
    dir.create(dirname(output_xlsx), recursive = TRUE, showWarnings = FALSE)
    openxlsx::write.xlsx(
      x = list(
        train_ite_results = train_results_tbl,
        test_ite_results = test_results_tbl
      ),
      file = output_xlsx,
      overwrite = TRUE
    )
  }

  # 按需打印运行摘要，便于快速核对结果
  if (isTRUE(verbose)) {
    cat("\n===== 训练集 ITE 结果摘要 =====\n")
    cat("训练集样本量：", n_train, "\n", sep = "")
    cat("测试集样本量：", n_test, "\n", sep = "")
    cat("设计矩阵列数：", ncol(x_train), "\n", sep = "")
    cat("bootstrap 次数：", B, "\n", sep = "")
    cat(
      "训练集 ite_hat 均值：",
      round(mean(train_results_tbl$ite_hat, na.rm = TRUE), 6),
      "\n",
      sep = ""
    )
    cat(
      "训练集 ite_hat 中位数：",
      round(stats::median(train_results_tbl$ite_hat, na.rm = TRUE), 6),
      "\n",
      sep = ""
    )
    cat(
      "训练集 CI 不跨 0 的患者数：",
      sum(train_results_tbl$ci_excludes_zero, na.rm = TRUE),
      "\n",
      sep = ""
    )
    print(train_results_tbl %>% dplyr::slice_head(n = 10))

    cat("\n===== 测试集 ITE 结果摘要 =====\n")
    cat(
      "测试集 ite_hat 均值：",
      round(mean(test_results_tbl$ite_hat, na.rm = TRUE), 6),
      "\n",
      sep = ""
    )
    cat(
      "测试集 ite_hat 中位数：",
      round(stats::median(test_results_tbl$ite_hat, na.rm = TRUE), 6),
      "\n",
      sep = ""
    )
    cat(
      "测试集 CI 不跨 0 的患者数：",
      sum(test_results_tbl$ci_excludes_zero, na.rm = TRUE),
      "\n",
      sep = ""
    )
    print(test_results_tbl %>% dplyr::slice_head(n = 10))
  }

  # 返回完整结果，便于后续在会话中继续分析
  list(
    base_model = base_model,
    train_results_tbl = train_results_tbl,
    test_results_tbl = test_results_tbl,
    boot_pred_train_mat = boot_pred_train_mat,
    boot_pred_test_mat = boot_pred_test_mat,
    design_colnames = design_obj$design_colnames
  )
}

# ----------------------------------------------------------
# 6. 示例：基于当前项目中的训练/测试文件直接运行
# ----------------------------------------------------------
# 指定训练集文件路径
train_file <- file.path("data", "final_data_train.xlsx")

# 指定测试集文件路径
test_file <- file.path("data", "final_data_test.xlsx")

# 指定训练集 CSV 结果输出路径
output_train_csv_file <- file.path("results", "survival_ite_ci_train_results_boot100.csv")

# 指定测试集 CSV 结果输出路径
output_test_csv_file <- file.path("results", "survival_ite_ci_test_results_boot100.csv")

# 指定 Excel 汇总结果输出路径
output_xlsx_file <- file.path("results", "survival_ite_ci_train_test_results_boot100.xlsx")

# 指定与原训练脚本保持一致的 17 个默认建模变量
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

# 读取训练集数据
train_data <- openxlsx::read.xlsx(train_file)

# 读取测试集数据
test_data <- openxlsx::read.xlsx(test_file)

# 调用主函数，直接在当前训练/测试数据上完成训练集与测试集 ITE 点估计及区间计算
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

# 读取训练集结果表到当前对象，便于后续继续查看
train_ite_tbl <- train_test_ite_res$train_results_tbl

# 读取测试集结果表到当前对象，便于后续继续查看
test_ite_tbl <- train_test_ite_res$test_results_tbl

# 再打印一次训练集按 ite_hat 从高到低排列后的前 10 位患者，便于快速检查
cat("\n===== 训练集 ITE 排名前 10 =====\n")
print(
  train_ite_tbl %>%
    dplyr::arrange(dplyr::desc(ite_hat)) %>%
    dplyr::slice_head(n = 10)
)

# 再打印一次测试集按 ite_hat 从高到低排列后的前 10 位患者，便于快速检查
cat("\n===== 测试集 ITE 排名前 10 =====\n")
print(
  test_ite_tbl %>%
    dplyr::arrange(dplyr::desc(ite_hat)) %>%
    dplyr::slice_head(n = 10)
)

# 打印结果文件路径，便于定位输出
cat("\n训练集 CSV 结果已写出到：", normalizePath(output_train_csv_file), "\n", sep = "")
cat("测试集 CSV 结果已写出到：", normalizePath(output_test_csv_file), "\n", sep = "")
cat("训练/测试 Excel 结果已写出到：", normalizePath(output_xlsx_file), "\n", sep = "")
