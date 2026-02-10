library(personalized)   
library(here)          
library(tidyverse)      
library(openxlsx)      
library(survival)      
library(fastDummies)    
library(MatchIt)
library("MatchIt")
library(tidyverse)
library(here)
library(openxlsx)
library(cobalt)
library(compareGroups)
library(survival)
library(broom)
library(personalized)   
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


rm(list = ls())         
df <- read.xlsx(here("data", "test.xlsx"))
head(df)

# ============================================================================
# 分层随机抽样：从df中随机抽取300例数据，保证每个分类变量的类别都包括
# ============================================================================

# 设置随机种子，确保结果可重复
set.seed(123)

# 识别所有分类变量（factor类型或character类型）
# sapply函数遍历df的每一列，检查是否为factor或character类型
categorical_vars <- names(df)[sapply(df, function(x) is.factor(x) | is.character(x))]

# 打印识别到的分类变量名称，便于检查
print(paste("识别到的分类变量:", paste(categorical_vars, collapse = ", ")))

# 检查每个分类变量的类别数量，确保有足够的样本进行分层抽样
# 如果某个类别的样本数太少，可能需要调整抽样策略
for(var in categorical_vars) {
  cat_counts <- table(df[[var]], useNA = "ifany")  # 统计每个类别的频数，包括NA值
  print(paste("变量", var, "的类别分布:"))
  print(cat_counts)
}

# 方法1：使用dplyr进行分层抽样（适用于单个主要分类变量）
# 如果数据中有多个分类变量，可以选择最重要的一个进行分层
# 这里我们使用所有分类变量的组合进行分层

# 创建分层标识：将所有分类变量的值组合成一个字符串，作为分层标识
# 这样可以确保每个分类变量的每个类别组合都有代表
df <- df %>%
  dplyr::mutate(
    # 将所有分类变量的值用"_"连接，创建唯一的分层标识
    stratum_id = paste(!!!syms(categorical_vars), sep = "_")
  )

# 计算每个分层的样本数
stratum_counts <- df %>%
  dplyr::count(stratum_id) %>%
  dplyr::arrange(desc(n))

# 打印每个分层的样本数
print("各分层的样本数:")
print(stratum_counts)

# 计算每个分层应该抽取的样本数（按比例分配）
# 如果某个分层的样本数不足，则全部抽取
total_samples <- 300  # 目标总样本数

# 按比例分配样本数到每个分层
stratum_samples <- df %>%
  dplyr::count(stratum_id) %>%
  dplyr::mutate(
    # 计算每个分层应抽取的样本数（按比例，但至少抽取1个）
    sample_size = pmax(1, round(n * total_samples / nrow(df)))
  ) %>%
  # 如果按比例分配的总数超过300，则按比例缩减
  dplyr::mutate(
    sample_size = if(sum(sample_size) > total_samples) {
      # 按比例缩减，但确保每个分层至少抽取1个
      floor(sample_size * total_samples / sum(sample_size))
    } else {
      sample_size
    }
  )

# 确保每个分层至少抽取1个样本（如果该分层有数据）
stratum_samples <- stratum_samples %>%
  dplyr::mutate(
    sample_size = pmax(1, pmin(sample_size, n))  # 至少1个，最多不超过该分层的总数
  )

# 如果分配的总数不足300，从样本数较多的分层中补充
current_total <- sum(stratum_samples$sample_size)
if(current_total < total_samples) {
  # 计算还需要抽取的样本数
  remaining <- total_samples - current_total
  # 从样本数较多的分层中补充
  stratum_samples <- stratum_samples %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(
      # 按顺序补充，直到达到目标总数
      sample_size = if(row_number() <= remaining) {
        pmin(sample_size + 1, n)  # 增加1个，但不超过该分层的总数
      } else {
        sample_size
      }
    ) %>%
    dplyr::arrange(stratum_id)  # 恢复原始顺序
}

# 执行分层抽样
# 使用purrr::map_dfr来对每个分层分别进行抽样
df_sampled <- purrr::map_dfr(
  .x = unique(df$stratum_id),  # 遍历每个唯一的分层标识
  .f = function(stratum) {
    # 获取当前分层的样本
    stratum_data <- df %>% dplyr::filter(stratum_id == stratum)
    # 获取该分层应该抽取的样本数
    sample_size <- stratum_samples$sample_size[stratum_samples$stratum_id == stratum]
    # 从该分层中随机抽取指定数量的样本
    stratum_data %>%
      dplyr::sample_n(size = min(sample_size, nrow(stratum_data)), replace = FALSE)
  }
)

# 如果抽样后的总数不等于300，进行微调
if(nrow(df_sampled) != total_samples) {
  if(nrow(df_sampled) < total_samples) {
    # 如果样本数不足，从剩余数据中随机补充
    df_remaining <- df %>%
      dplyr::anti_join(df_sampled, by = names(df))  # 获取未抽中的样本
    
    # 随机抽取补充样本
    additional_samples <- total_samples - nrow(df_sampled)
    if(nrow(df_remaining) >= additional_samples) {
      df_additional <- df_remaining %>%
        dplyr::sample_n(additional_samples, replace = FALSE)
      df_sampled <- dplyr::bind_rows(df_sampled, df_additional)
    } else {
      # 如果剩余样本不足，则全部加入
      df_sampled <- dplyr::bind_rows(df_sampled, df_remaining)
    }
  } else {
    # 如果样本数过多，随机删除多余的样本
    df_sampled <- df_sampled %>%
      dplyr::sample_n(total_samples, replace = FALSE)
  }
}

# 移除临时创建的分层标识列
df_sampled <- df_sampled %>%
  dplyr::select(-stratum_id)

# 验证抽样结果：检查每个分类变量的类别是否都包含在抽样数据中
print("\n========== 抽样结果验证 ==========")
all_categories_present <- TRUE  # 标记是否所有类别都被包含

for(var in categorical_vars) {
  # 统计原始数据和抽样数据中各类别的频数（包括NA值）
  cat_counts_sampled <- table(df_sampled[[var]], useNA = "ifany")
  cat_counts_original <- table(df[[var]], useNA = "ifany")
  
  print(paste("\n变量", var, ":"))
  print("原始数据类别分布:")
  print(cat_counts_original)
  print("抽样数据类别分布:")
  print(cat_counts_sampled)
  
  # 检查是否所有类别都被包含
  # 获取原始数据中的所有类别名称（包括NA）
  original_cats <- names(cat_counts_original)
  sampled_cats <- names(cat_counts_sampled)
  
  # 检查缺失的类别
  missing_categories <- setdiff(original_cats, sampled_cats)
  if(length(missing_categories) > 0) {
    warning(paste("警告：变量", var, "的以下类别在抽样数据中缺失:", paste(missing_categories, collapse = ", ")))
    all_categories_present <- FALSE
  } else {
    print(paste("✓ 变量", var, "的所有类别都已包含在抽样数据中"))
  }
}

# 打印总体验证结果
if(all_categories_present) {
  print("\n✓ 所有分类变量的所有类别都已包含在抽样数据中")
} else {
  print("\n⚠ 部分分类变量的某些类别在抽样数据中缺失，请检查警告信息")
}

# 打印抽样结果摘要
print("\n========== 抽样摘要 ==========")
print(paste("原始数据样本数:", nrow(df)))
print(paste("抽样后样本数:", nrow(df_sampled)))
print(paste("目标样本数:", total_samples))
print(paste("抽样比例:", round(nrow(df_sampled) / nrow(df) * 100, 2), "%"))

# 将抽样后的数据赋值回df（如果需要替换原始数据）
# df <- df_sampled

# 或者创建新的变量保存抽样数据（保留原始df）
df <- df_sampled  # 保存300例抽样数据到新变量df_300

exposure.formu <- Trt_ITT ~ sex + age_grade + ecog + cancer_grade + complication + PDL1 + CPgrade + Distant_metastasis + TNM + disease_status                      
matched_data <- matchit(exposure.formu,
                data = df,
                method = "nearest",
                distance = "logit",
                link = "logit",
                estimand = "ATT",
                ratio = 1,
                caliper = 0.05)
summary(matched_data)

final_data <- match.data(matched_data, 
                         group = "all")  
head(final_data)


x_matched <- final_data %>%
  select(Tumor_border_status,CD4_CD8,Peritumoral_enhancement, 
  CRP,CPgrade,Radiologic_morphology,NLR,CA199,Neutro,Cr,
  PDL1,DBIL,BUN,ALB,TumorDiameterMax,
  Vascular_invasion,Satellite_nodules,  
  )

x_dummy_matched <- fastDummies::dummy_cols(
  x_matched,
  remove_selected_columns = TRUE,
  remove_first_dummy = TRUE
)
x_matched <- as.matrix(x_dummy_matched)
y.time.to.event <- final_data$OS %>% as.numeric()
y.status <- final_data$OS_CNSR %>% as.integer()
trt <- final_data$Trt_ITT %>% as.integer()
subclass <- final_data$subclass %>% as.integer()


set.seed(123)
subgrp_model <- fit.subgroup(
    x = x_matched,
    y = Surv(y.time.to.event, y.status),
    trt = trt,
    match.id = subclass,  
    loss = "cox_loss_lasso",
    method = "weighting",
    cutpoint = 0
)
summary(subgrp_model)
plot(subgrp_model)
plot(subgrp_model, type = "interaction")

validation <- validate.subgroup(subgrp_model, 
                                B = 100L, 
                                method = "training_test_replication",
                                train.fraction = 0.80)
validation
plot(validation)
plot(validation, type = "interaction")


trt <-subgrp_model$trts
y <-subgrp_model$y
pi.x <- subgrp_model$pi.x
wts <- 1 / (pi.x * (trt[1]) + (1 - pi.x) * (trt[2]))
resp <-y*wts
surv_obj <- subgrp_model$y
resp <- surv_obj[, 2]  
resp_weighted <- resp * wts
recom<-subgrp_model$recommended.trts
TX<-subgrp_model$trt.received
train_cohort <-data.frame(resp,TX,recom)
anov <- aov(resp ~ TX * recom, data = train_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Discovery p-int : %f',pvalue)
rm(pvalue)



score_pred <- predict(subgrp_model, newx = x_matched)
head(score_pred)
class(score_pred)
score_df <- data.frame(score_pred)
head(score_df)



benefit.scores <- predict(subgrp_model, newx = x_matched, type = "benefit.score")
head(benefit.scores)
data_clean <- final_data %>% dplyr::mutate(ITE = benefit.scores) %>% 
 dplyr::mutate(                                                        
    Benefit_Level = dplyr::case_when(
      ITE < quantile(ITE, 0.33, na.rm=TRUE) ~ "Negative",      
      ITE >= quantile(ITE, 0.33, na.rm=TRUE) &                             
        ITE < quantile(ITE, 0.66, na.rm=TRUE)           ~ "Moderate",
      ITE >= quantile(ITE, 0.66, na.rm=TRUE)            ~ "High",      
      TRUE                                              ~ "Unknown"     
    ))
data_Negative <- data_clean %>% dplyr::filter(Benefit_Level == "Negative")
data_High <- data_clean %>% dplyr::filter(Benefit_Level == "High")
data_Moderate <- data_clean %>% dplyr::filter(Benefit_Level == "Moderate")
group_counts <- data_clean %>% dplyr::count(Benefit_Level)
print(group_counts)
surfit_Negative <- survival::survfit(survival::Surv(OS, OS_CNSR) ~ Trt_ITT, data = data_Negative)
surfit_Moderate <- survival::survfit(survival::Surv(OS, OS_CNSR) ~ Trt_ITT, data = data_Moderate)
surfit_High <- survival::survfit(survival::Surv(OS, OS_CNSR) ~ Trt_ITT, data = data_High)



rec.trt.grp <- predict(subgrp_model, newx = x_matched, type = "trt.group")
head(rec.trt.grp)



save.image("personalized_ICC.RData")
load("personalized_ICC.RData")

