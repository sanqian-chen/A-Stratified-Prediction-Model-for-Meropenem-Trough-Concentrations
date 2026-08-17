rm(list=ls())
setwd("E:/哥/撰稿文章/美罗培南文章/数据/R/训练集和验证集合并")
library(readxl)
data <- read_excel("data.xlsx")

#dataset为1——训练集test_data，dataset为0——外部验证集valid_data


# 变量分组及类别 -----------------------------------------------------------------
# 查看 data 数据框所有变量的类型
str(data)

# 血药浓度分组
data$血药浓度分组 <- ifelse(
  data$浓度 < 2, 1,
  ifelse(data$浓度 <= 8, 2, 3)
)
# 确保血药浓度分组为有序因子
data$血药浓度分组 <- factor(data$血药浓度分组, levels = c(1,2,3), ordered = TRUE)
# 查看结果（确保分组正确）
table(
  数据集 = ifelse(data$dataset == 1, "训练集", "外部验证集"),
  血药浓度分组 = data$血药浓度分组
)


#CKDEPI分组
#严重程度分组：根据是否行CRRT和CKD-EPI指标
library(dplyr)
data <- data %>%
  mutate(
    CKDEPI分组 = case_when(
      是否行CRRT == 1 ~ 5,                # 行CRRT → 极重
      CKDEPI < 10 ~ 4,                    # <10 → 重
      CKDEPI >= 10 & CKDEPI <= 25 ~ 3,    # 10-25 → 中
      CKDEPI > 25 & CKDEPI <= 50 ~ 2,     # 25-50 → 轻
      CKDEPI > 50 ~ 1,                    # >50 → 正常
      TRUE ~ NA_integer_                   # 其他情况（如缺失）设为NA
    ))


# 统计例数
result <- data %>%
  mutate(
    数据集 = ifelse(dataset == 1, "训练集", "外部验证集")
  ) %>%
  group_by(数据集, CKDEPI分组, 血药浓度分组) %>%
  summarise(例数 = n(), .groups = "drop") %>%
  arrange(数据集, CKDEPI分组, 血药浓度分组)

result




# 二分类变量定义为因子factor
binary_vars <- c("性别", "是否有创机械通气", "是否输血", 
                 "是否使用血管活性药物", "是否行ECMO")

# 将二分类变量从数值转换为因子，不添加标签
for(var in binary_vars) {
  if(var %in% names(data)) {
    # 获取变量的唯一值
    unique_vals <- sort(unique(data[[var]]))
    
    # 将变量转换为因子，不添加标签
    data[[var]] <- factor(data[[var]], levels = unique_vals)
  }
}






# 感染病灶列名可能为 "感染病灶"
library(tidyr)
library(dplyr)

# 1. 拆分、去重、清理空格
df_split <- data %>%
  separate_rows(感染病灶, sep = "[、,，]") %>%
  mutate(感染病灶 = trimws(感染病灶))  # 去除首尾空格

# 2. 合并感染病灶类别
df_split <- df_split %>%
  mutate(感染病灶合并 = case_when(
    感染病灶 %in% c("肺部", "肺") ~ "肺部",
    感染病灶 %in% c("血流", "主动脉瓣赘生物") ~ "血流",
    感染病灶 %in% c("腹腔", "腹膜炎", "胆道", "肝脏", "腹部", "腹膜透析中腹腔感染") ~ "腹腔",
    感染病灶 %in% c("泌尿", "泌尿道", "泌尿系", "尿路") ~ "泌尿系统",
    感染病灶 %in% c("中枢", "颅内", "中枢神经系统") ~ "中枢神经系统",
    TRUE ~ "其他"
  ))

# 3. 创建数据集标签
df_split <- df_split %>%
  mutate(数据集 = ifelse(dataset == 1, "训练集", "外部验证集"))

# 4. 分组统计各数据集中每个合并类别的病例数
感染病灶统计 <- df_split %>%
  group_by(数据集, 感染病灶合并) %>%
  summarise(病例数 = n(), .groups = "drop") %>%
  arrange(数据集, 感染病灶合并)

感染病灶统计








# CKDEPI分组变量的Jonckheere-Terpstra趋势检验+Spearman秩相关+百分比条图 ---------------------------------------------------
library(dplyr)
library(clinfun)

# 确保变量为数值型有序等级
data <- data %>%
  mutate(
    CKDEPI分组_num = as.numeric(as.character(CKDEPI分组)),
    血药浓度分组_num = as.numeric(as.character(血药浓度分组))
  )

# 分出训练集和外部验证集
train_data <- data %>% filter(dataset == 1)
valid_data <- data %>% filter(dataset == 0)

#训练集
table(train_data$CKDEPI分组, train_data$血药浓度分组)
#Jonckheere-Terpstra趋势检验
jonckheere.test(
  x = train_data$血药浓度分组_num,
  g = train_data$CKDEPI分组_num,
  alternative = "two.sided",
  nperm = 10000   # 使用置换法得到精确 P 值
)
#Spearman秩相关
cor.test(
  train_data$CKDEPI分组_num,
  train_data$血药浓度分组_num,
  method = "spearman",
  exact = FALSE
)


#验证集
table(valid_data$CKDEPI分组, valid_data$血药浓度分组)
#Jonckheere-Terpstra趋势检验
jonckheere.test(
  x = valid_data$血药浓度分组_num,
  g = valid_data$CKDEPI分组_num,
  alternative = "two.sided",
  nperm = 10000
)
#Spearman秩相关
cor.test(
  valid_data$CKDEPI分组_num,
  valid_data$血药浓度分组_num,
  method = "spearman",
  exact = FALSE
)




#百分比条图
library(dplyr)
library(ggplot2)
library(scales)


plot_data <- data %>%
  mutate(
    数据集 = ifelse(
      dataset == 1,
      "训练集",
      "外部验证集"
    )
  ) %>%
  count(
    数据集,
    CKDEPI分组,
    血药浓度分组
  ) %>%
  group_by(
    数据集,
    CKDEPI分组
  ) %>%
  mutate(
    Percent = n / sum(n)
  )

ggplot(
  plot_data,
  aes(
    x = factor(CKDEPI分组),
    y = Percent,
    fill = factor(血药浓度分组)
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill",
    color = "black"
  ) +
  facet_wrap(
    ~数据集,
    nrow = 1
  ) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  labs(
    x = "CKD-EPI group",
    y = "Percentage",
    fill = "Drug concentration group"
  ) +
  theme_bw(base_size = 14)





####修改成英文格式
#百分比条图
library(dplyr)
library(ggplot2)
library(scales)

plot_data <- data %>%
  mutate(
    数据集 = factor(
      ifelse(
        dataset == 1,
        "Training set",             # dataset == 1 → 训练集
        "External validation set"   # dataset == 0 → 外部验证集
      ),
      levels = c("Training cohort", "Validation cohort")  # 强制Training在左，External在右
    )
  ) %>%
  count(
    数据集,
    CKDEPI分组,
    血药浓度分组
  ) %>%
  group_by(
    数据集,
    CKDEPI分组
  ) %>%
  mutate(
    Percent = n / sum(n)
  )

ggplot(
  plot_data,
  aes(
    x = factor(CKDEPI分组),
    y = Percent,
    fill = factor(血药浓度分组, 
                  levels = c(1, 2, 3),
                  labels = c("Subtherapeutic group", "Therapeutic group", "High-exposure group"))
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill",
    color = "black"
  ) +
  facet_wrap(
    ~数据集,
    nrow = 1
  ) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  labs(
    x = "CKD-EPI group",
    y = "Percentage",
    fill = "Drug concentration group"
  ) +
  theme_bw(base_size = 14)











# 训练集和外部验证集患者的基线和临床特征的单因素分析 -----------------------------------------------
#变量感染病灶合并
# 构建列联表：数据集 × 感染病灶合并
tab <- table(df_split$数据集, df_split$感染病灶合并)
tab
#卡方检验
chi_result <- chisq.test(tab)
chi_result
chi_result$expected


# CKDEPI分组变量比较 
tab_ckd <- table(
  data$数据集,
  data$CKDEPI分组
)
tab_ckd
#卡方检验
chi_ckd <- chisq.test(tab_ckd)
chi_ckd
chi_ckd$expected

#Mann-Whitney U检验
wilcox.test(
  as.numeric(CKDEPI分组) ~ 数据集,
  data = data,
  exact = FALSE
)



#在训练集和外部验证集中，开始对连续型定量变量和二分类变量进行单因素分析
library(dplyr)

# 设置数据集标签
data <- data %>%
  mutate(
    数据集 = factor(
      dataset,
      levels = c(1, 0),
      labels = c("训练集", "外部验证集")
    )
  )

# 定义变量
# 二分类变量
bin_vars <- c(
  "性别",
  "是否有创机械通气",
  "是否输血",
  "是否使用血管活性药物",
  "是否行ECMO",
  "存活死亡",
  "是否行CRRT"
)

# 连续型定量变量
cont_vars <- c(
  "年龄",
  "肌酐",
  "尿素氮",
  "白蛋白",
  "ALT",
  "AST",
  "总胆红素",
  "血红蛋白",
  "血小板",
  "APACHEII",
  "PT",
  "APTT",
  "D二聚体",
  "Fbg",
  "测量血药浓度时肌酐",
  "负平衡量",
  "使用剂量",
  "输注时间",
  "入住ICU时间",
  "无呼吸机时间",
  "无血管活性药物使用时间",
  "无CRRT支持时间"
)

# 连续变量：中位数(Q1,Q3) + Mann-Whitney或t检验
median_iqr <- function(x) {
  x <- x[!is.na(x)]
  med <- median(x)
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  paste0(round(med, 2), " (", round(q1,2), ", ", round(q3,2), ")")
}

continuous_results <- lapply(cont_vars, function(var) {
  
  x_train <- data %>% filter(数据集 == "训练集") %>% pull(var) %>% na.omit()
  x_valid <- data %>% filter(数据集 == "外部验证集") %>% pull(var) %>% na.omit()
  
  # Shapiro-Wilk正态性检验
  shapiro_train <- shapiro.test(x_train)
  shapiro_valid <- shapiro.test(x_valid)
  
  normal_train <- shapiro_train$p.value > 0.05
  normal_valid <- shapiro_valid$p.value > 0.05
  
  if(normal_train & normal_valid) {
    test <- t.test(data[[var]] ~ data$数据集, var.equal = FALSE)
    method <- "两独立样本t检验"
    stat <- unname(test$statistic)
    p <- test$p.value
  } else {
    test <- wilcox.test(data[[var]] ~ data$数据集, exact = FALSE)
    method <- "Mann-Whitney U检验"
    stat <- unname(test$statistic)
    p <- test$p.value
  }
  
  data.frame(
    变量 = var,
    训练集 = median_iqr(x_train),
    外部验证集 = median_iqr(x_valid),
    方法 = method,
    统计量 = stat,
    P值 = p,
    stringsAsFactors = FALSE
  )
})

continuous_results <- do.call(rbind, continuous_results)

# 二分类变量：例数 + Pearson卡方检验
# 二分类变量：例数 + 百分比 + Pearson卡方检验

binary_results <- lapply(bin_vars, function(var) {
  
  # 构建列联表
  tab <- table(data$数据集, data[[var]], useNA = "no")
  
  # Pearson卡方检验
  test <- chisq.test(tab, correct = FALSE)
  
  # 获取变量类别名称
  levels_var <- colnames(tab)
  
  # 训练集和验证集例数
  train_n <- tab["训练集", ]
  valid_n <- tab["外部验证集", ]
  
  # 计算百分比
  train_pct <- round(train_n / sum(train_n) * 100, 1)
  valid_pct <- round(valid_n / sum(valid_n) * 100, 1)
  
  # 拼接输出格式：类别 n (%)
  train_str <- paste0(
    levels_var,
    ": ",
    train_n,
    " (",
    train_pct,
    "%)",
    collapse = "; "
  )
  
  valid_str <- paste0(
    levels_var,
    ": ",
    valid_n,
    " (",
    valid_pct,
    "%)",
    collapse = "; "
  )
  
  
  data.frame(
    变量 = var,
    训练集 = train_str,
    外部验证集 = valid_str,
    方法 = "Pearson卡方检验",
    统计量 = unname(test$statistic),
    P值 = test$p.value,
    stringsAsFactors = FALSE
  )
})


binary_results <- do.call(rbind, binary_results)

# 合并连续变量和二分类变量结果
baseline_results <- bind_rows(
  continuous_results,
  binary_results
)

baseline_results

#导出Excel文件
library(writexl)
write_xlsx(
  baseline_results,
  path = "训练集与外部验证集基线特征比较.xlsx"
)







# 训练集中3个血药浓度分组患者的基线、临床特征及预后变量比较 -------------------

library(dplyr)
library(tidyr)
library(openxlsx)

# 1. 提取训练集 

train_data <- data %>%
  filter(dataset == 1) %>%
  mutate(
    血药浓度分组 = factor(
      血药浓度分组,
      levels = c(1, 2, 3),
      ordered = TRUE
    ),
    血药浓度分组_num = as.numeric(血药浓度分组),
    
    # 存活死亡：0 = 死亡，1 = 存活
    存活死亡 = factor(
      存活死亡,
      levels = c(0, 1),
      labels = c("死亡", "存活")
    )
  )


# 2. 定义变量 

# 二分类变量：加入“存活死亡”
binary_vars <- c(
  "性别",
  "是否有创机械通气",
  "是否输血",
  "是否使用血管活性药物",
  "是否行ECMO",
  "存活死亡",
  "是否行CRRT"
)

# 连续型定量变量：加入4个预后连续变量
cont_vars <- c(
  "年龄",
  "肌酐",
  "尿素氮",
  "白蛋白",
  "ALT",
  "AST",
  "总胆红素",
  "血红蛋白",
  "血小板",
  "APACHEII",
  "PT",
  "APTT",
  "D二聚体",
  "Fbg",
  "负平衡量",
  "使用剂量",
  "输注时间",
  "测量血药浓度时肌酐",
  
  # 新增4个连续型预后变量
  "入住ICU时间",
  "无呼吸机时间",
  "无血管活性药物使用时间",
  "无CRRT支持时间"
)


# 3. 连续变量：中位数(Q1,Q3) + Kruskal-Wallis H检验 

median_iqr <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA)
  }
  
  med <- median(x)
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  
  sprintf("%.2f(%.2f,%.2f)", med, q1, q3)
}


continuous_results <- lapply(cont_vars, function(var) {
  
  if (!var %in% names(train_data)) {
    return(NULL)
  }
  
  # 描述性统计
  desc <- train_data %>%
    group_by(血药浓度分组) %>%
    summarise(
      统计量 = median_iqr(.data[[var]]),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = 血药浓度分组,
      values_from = 统计量,
      names_prefix = "组"
    )
  
  # Kruskal-Wallis H检验
  test <- kruskal.test(
    train_data[[var]] ~ train_data$血药浓度分组
  )
  
  desc %>%
    mutate(
      变量 = var,
      变量类型 = "连续型定量变量",
      检验方法 = "Kruskal-Wallis H检验",
      统计量_H = round(unname(test$statistic), 3),
      P值 = test$p.value
    ) %>%
    select(
      变量,
      变量类型,
      组1,
      组2,
      组3,
      检验方法,
      统计量_H,
      P值
    )
})

continuous_results <- bind_rows(continuous_results)


# 4. 二分类变量：例数 + Mann-Whitney U检验 

binary_results <- lapply(binary_vars, function(var) {
  
  if (!var %in% names(train_data)) {
    return(NULL)
  }
  
  # 根据变量类型设置因子水平和显示方式
  if (var == "性别") {
    
    # 性别：假设1=男，2=女，显示为 男/女
    train_data[[var]] <- factor(
      train_data[[var]],
      levels = c(1, 2),
      labels = c("男", "女")
    )
    
    row_levels <- c("男", "女")
    display_type <- "二分类变量，显示格式：男/女"
    
  } else if (var == "存活死亡") {
    
    # 存活死亡：0=死亡，1=存活，显示为 死亡/存活
    train_data[[var]] <- factor(
      train_data[[var]],
      levels = c("死亡", "存活")
    )
    
    row_levels <- c("死亡", "存活")
    display_type <- "二分类变量，显示格式：死亡/存活"
    
  } else {
    
    # 其他二分类变量：0=否，1=是，显示为 是/否
    train_data[[var]] <- factor(
      train_data[[var]],
      levels = c(0, 1),
      labels = c("否", "是")
    )
    
    row_levels <- c("是", "否")
    display_type <- "二分类变量，显示格式：是/否"
  }
  
  
  # Mann-Whitney U检验
  # 比较不同二分类状态下，血药浓度分组等级分布是否不同
  test <- wilcox.test(
    train_data$血药浓度分组_num ~ train_data[[var]],
    correct = TRUE,
    exact = FALSE
  )
  
  
  # 列联表：行=二分类变量水平，列=血药浓度分组
  tab <- table(
    train_data[[var]],
    train_data$血药浓度分组
  )
  
  group_levels <- c("1", "2", "3")
  
  # 补齐可能缺失的行或列，避免下标报错
  tab_full <- matrix(
    0,
    nrow = length(row_levels),
    ncol = length(group_levels),
    dimnames = list(row_levels, group_levels)
  )
  
  common_rows <- intersect(rownames(tab_full), rownames(tab))
  common_cols <- intersect(colnames(tab_full), colnames(tab))
  
  tab_full[common_rows, common_cols] <- tab[common_rows, common_cols]
  
  
  # 按变量类型生成每组描述
  if (var == "性别") {
    
    group1_str <- paste0(tab_full["男", "1"], "/", tab_full["女", "1"])
    group2_str <- paste0(tab_full["男", "2"], "/", tab_full["女", "2"])
    group3_str <- paste0(tab_full["男", "3"], "/", tab_full["女", "3"])
    
  } else if (var == "存活死亡") {
    
    group1_str <- paste0(tab_full["死亡", "1"], "/", tab_full["存活", "1"])
    group2_str <- paste0(tab_full["死亡", "2"], "/", tab_full["存活", "2"])
    group3_str <- paste0(tab_full["死亡", "3"], "/", tab_full["存活", "3"])
    
  } else {
    
    group1_str <- paste0(tab_full["是", "1"], "/", tab_full["否", "1"])
    group2_str <- paste0(tab_full["是", "2"], "/", tab_full["否", "2"])
    group3_str <- paste0(tab_full["是", "3"], "/", tab_full["否", "3"])
  }
  
  
  data.frame(
    变量 = var,
    变量类型 = display_type,
    组1 = group1_str,
    组2 = group2_str,
    组3 = group3_str,
    检验方法 = "Mann-Whitney U检验",
    统计量_W = round(unname(test$statistic), 3),
    P值 = test$p.value,
    stringsAsFactors = FALSE
  )
})

binary_results <- bind_rows(binary_results)


# 5. 合并连续变量和二分类变量结果 

train_univariate_results <- bind_rows(
  continuous_results,
  binary_results
) %>%
  mutate(
    P值_格式 = ifelse(
      P值 < 0.001,
      "<0.001",
      sprintf("%.3f", P值)
    )
  )


# 6. 查看结果 

print(train_univariate_results)


# 7. 导出 Excel 

write.xlsx(
  train_univariate_results,
  file = "训练集_血药浓度分组单因素分析_含5个预后变量.xlsx",
  rowNames = FALSE
)










#训练集数据中，单独对df_split数据中感染病灶合并变量进行单因素分析
str(df_split)
library(dplyr)
library(tidyr)

# 1. 提取df_split中的训练集
df_train_split <- df_split %>%
  filter(dataset == 1) %>%
  mutate(
    血药浓度分组 = factor(
      血药浓度分组,
      levels = c(1, 2, 3),
      ordered = TRUE
    ),
    血药浓度分组_num = as.numeric(as.character(血药浓度分组)),
    感染病灶合并 = factor(感染病灶合并)
  )

# 2. 统计训练集中三个血药浓度分组下6类感染病灶合并的例数
感染病灶_血药浓度统计 <- df_train_split %>%
  group_by(感染病灶合并, 血药浓度分组) %>%
  summarise(例数 = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = 血药浓度分组,
    values_from = 例数,
    values_fill = 0,
    names_prefix = "组"
  ) %>%
  arrange(感染病灶合并)

感染病灶_血药浓度统计

# 3. Kruskal-Wallis H检验
kw_infection <- kruskal.test(
  血药浓度分组_num ~ 感染病灶合并,
  data = df_train_split
)

kw_infection

# 4. 整理检验结果
感染病灶_单因素结果 <- data.frame(
  变量 = "感染病灶合并",
  变量类型 = "无序多分类变量",
  检验方法 = "Kruskal-Wallis H检验",
  统计量_H = round(unname(kw_infection$statistic), 3),
  自由度 = unname(kw_infection$parameter),
  P值 = kw_infection$p.value,
  P值_格式 = ifelse(
    kw_infection$p.value < 0.001,
    "<0.001",
    sprintf("%.3f", kw_infection$p.value)
  )
)

感染病灶_单因素结果






# 5个预后指标绘制表格 --------------------------------------------------------------
# 根据“存活死亡”计算28天死亡率 

library(dplyr)
library(tidyr)
library(openxlsx)

# 1. 整理变量编码 

data_mortality <- data %>%
  mutate(
    数据集 = factor(
      dataset,
      levels = c(1, 0),
      labels = c("训练集", "外部验证集")
    ),
    
    血药浓度分组 = factor(
      血药浓度分组,
      levels = c(1, 2, 3),
      labels = c("组1", "组2", "组3"),
      ordered = TRUE
    ),
    
    # 存活死亡：0 = 死亡，1 = 存活
    存活死亡 = case_when(
      as.character(存活死亡) %in% c("0", "死亡") ~ "死亡",
      as.character(存活死亡) %in% c("1", "存活") ~ "存活",
      TRUE ~ NA_character_
    )
  )


# 2. 定义死亡率计算函数 

calc_mortality <- function(df, group_var) {
  
  df %>%
    filter(!is.na(存活死亡)) %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      死亡例数 = sum(存活死亡 == "死亡"),
      存活例数 = sum(存活死亡 == "存活"),
      总例数 = 死亡例数 + 存活例数,
      `28天死亡率` = 死亡例数 / 总例数 * 100,
      .groups = "drop"
    ) %>%
    mutate(
      `28天死亡率（死亡/存活例数）` = paste0(
        sprintf("%.2f", `28天死亡率`),
        "%（",
        死亡例数,
        "/",
        存活例数,
        "）"
      )
    )
}


# 3. 训练集和外部验证集总体28天死亡率 

mortality_dataset <- calc_mortality(
  df = data_mortality,
  group_var = "数据集"
)

mortality_dataset


# 4. 训练集中不同血药浓度分组的28天死亡率 

mortality_train_group <- data_mortality %>%
  filter(数据集 == "训练集") %>%
  calc_mortality(
    group_var = "血药浓度分组"
  )

mortality_train_group



# 5. 比较训练集和外部验证集总体28天死亡率差异 
# 行 = 训练集/外部验证集，列 = 死亡/存活
tab_mortality_dataset <- table(
  data_mortality$数据集,
  data_mortality$存活死亡
)

tab_mortality_dataset

chisq.test(
  tab_mortality_dataset,
  correct = FALSE
)


# 6. 比较训练集中不同血药浓度分组的28天死亡率差异 

train_mortality_data <- data_mortality %>%
  filter(数据集 == "训练集")

# 行 = 血药浓度分组，列 = 死亡/存活
tab_mortality_train_group <- table(
  train_mortality_data$血药浓度分组,
  train_mortality_data$存活死亡
)

tab_mortality_train_group

chisq.test(
  tab_mortality_train_group,
  correct = FALSE
)









# 3个预后变量：筛选对应治疗人群后统计n和M(Q1,Q3) --------

library(dplyr)
library(tidyr)
library(openxlsx)

# 1. 定义变量及其筛选条件 

prognosis_info <- data.frame(
  变量 = c(
    "无呼吸机时间",
    "无血管活性药物使用时间",
    "无CRRT支持时间"
  ),
  筛选变量 = c(
    "是否有创机械通气",
    "是否使用血管活性药物",
    "是否行CRRT"
  ),
  纳入对象 = c(
    "是否有创机械通气=1",
    "是否使用血管活性药物=1",
    "是否行CRRT=1"
  ),
  stringsAsFactors = FALSE
)


# 2. 整理数据 

data_prognosis <- data %>%
  mutate(
    数据集 = factor(
      dataset,
      levels = c(1, 0),
      labels = c("训练集", "外部验证集")
    ),
    血药浓度分组 = factor(
      血药浓度分组,
      levels = c(1, 2, 3),
      labels = c("组1", "组2", "组3"),
      ordered = TRUE
    )
  )


# 3. 定义函数：计算n和M(Q1,Q3)

median_iqr_all <- function(x) {
  
  x <- as.numeric(as.character(x))
  
  n <- length(x)
  
  if (n == 0) {
    return("n=0；NA")
  }
  
  med <- median(x)
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  
  paste0(
    "n=", n,
    "；",
    sprintf("%.2f(%.2f,%.2f)", med, q1, q3)
  )
}


# 4. 维度一：训练集与外部验证集比较 

prognosis_dataset_results <- lapply(seq_len(nrow(prognosis_info)), function(i) {
  
  var <- prognosis_info$变量[i]
  filter_var <- prognosis_info$筛选变量[i]
  include_note <- prognosis_info$纳入对象[i]
  
  if (!var %in% names(data_prognosis) | !filter_var %in% names(data_prognosis)) {
    return(NULL)
  }
  
  # 只筛选对应治疗人群，不剔除0值，不剔除缺失值
  temp_data <- data_prognosis %>%
    filter(
      as.numeric(as.character(.data[[filter_var]])) == 1
    ) %>%
    mutate(
      value = as.numeric(as.character(.data[[var]]))
    )
  
  # 描述性统计
  train_stat <- median_iqr_all(
    temp_data$value[temp_data$数据集 == "训练集"]
  )
  
  valid_stat <- median_iqr_all(
    temp_data$value[temp_data$数据集 == "外部验证集"]
  )
  
  # Mann-Whitney U检验
  test <- wilcox.test(
    value ~ 数据集,
    data = temp_data,
    exact = FALSE
  )
  
  data.frame(
    变量 = var,
    纳入对象 = include_note,
    比较内容 = "训练集 vs 外部验证集",
    训练集 = train_stat,
    外部验证集 = valid_stat,
    检验方法 = "Mann-Whitney U检验",
    统计量_W = round(unname(test$statistic), 3),
    P值 = test$p.value,
    stringsAsFactors = FALSE
  )
})

prognosis_dataset_results <- bind_rows(prognosis_dataset_results) %>%
  mutate(
    P值_格式 = ifelse(
      P值 < 0.001,
      "<0.001",
      sprintf("%.3f", P值)
    )
  )

prognosis_dataset_results


# 5. 维度二：训练集中不同血药浓度分组比较 

train_prognosis_group_results <- lapply(seq_len(nrow(prognosis_info)), function(i) {
  
  var <- prognosis_info$变量[i]
  filter_var <- prognosis_info$筛选变量[i]
  include_note <- prognosis_info$纳入对象[i]
  
  if (!var %in% names(data_prognosis) | !filter_var %in% names(data_prognosis)) {
    return(NULL)
  }
  
  # 先筛选训练集，再筛选对应治疗人群，不剔除0值，不剔除缺失值
  temp_train_data <- data_prognosis %>%
    filter(数据集 == "训练集") %>%
    filter(
      as.numeric(as.character(.data[[filter_var]])) == 1
    ) %>%
    mutate(
      value = as.numeric(as.character(.data[[var]]))
    )
  
  # 描述性统计
  group1_stat <- median_iqr_all(
    temp_train_data$value[temp_train_data$血药浓度分组 == "组1"]
  )
  
  group2_stat <- median_iqr_all(
    temp_train_data$value[temp_train_data$血药浓度分组 == "组2"]
  )
  
  group3_stat <- median_iqr_all(
    temp_train_data$value[temp_train_data$血药浓度分组 == "组3"]
  )
  
  # Kruskal-Wallis H检验
  test <- kruskal.test(
    value ~ 血药浓度分组,
    data = temp_train_data
  )
  
  data.frame(
    变量 = var,
    纳入对象 = include_note,
    比较内容 = "训练集中不同血药浓度分组",
    组1 = group1_stat,
    组2 = group2_stat,
    组3 = group3_stat,
    检验方法 = "Kruskal-Wallis H检验",
    统计量_H = round(unname(test$statistic), 3),
    P值 = test$p.value,
    stringsAsFactors = FALSE
  )
})

train_prognosis_group_results <- bind_rows(train_prognosis_group_results) %>%
  mutate(
    P值_格式 = ifelse(
      P值 < 0.001,
      "<0.001",
      sprintf("%.3f", P值)
    )
  )

train_prognosis_group_results


