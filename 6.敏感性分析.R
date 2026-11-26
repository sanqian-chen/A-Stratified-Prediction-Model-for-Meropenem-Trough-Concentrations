# ============================================================================
# 美罗培南谷浓度分组的界值敏感性分析
#
# 主分析：    1 = <2 mg/L；2 = 2–8 mg/L； 3 = >8 mg/L
# 敏感性分析：1 = <2 mg/L；2 = 2–16 mg/L；3 = >16 mg/L
#
# dataset = 1：训练集
# dataset = 0：时间验证集
#
# 输出：分组人数、重分类交叉表、重分类比例、时间验证集性能对照
# ============================================================================


# 1. 读取数据 ----------------------------------------------------------------

# 如首次使用，可先单独运行：
# install.packages(c("readxl", "MASS", "writexl"))

work_dir <- "E:/哥/撰稿文章/美罗培南文章/数据/R/训练集和验证集合并"

dat <- as.data.frame(
  readxl::read_excel(file.path(work_dir, "data.xlsx"))
)

required_vars <- c(
  "dataset", "浓度", "是否行CRRT", "CKDEPI",
  "尿素氮", "血红蛋白", "血小板", "年龄", "APTT", "输注时间"
)

missing_vars <- setdiff(required_vars, names(dat))
if (length(missing_vars) > 0) {
  stop("Excel中缺少以下列：", paste(missing_vars, collapse = "、"))
}

# 将Excel中可能以文本储存的数字转换为数值。
# 如果出现不能识别的非空文本，明确报告列名，避免静默转成缺失值。
for (v in required_vars) {
  x <- trimws(as.character(dat[[v]]))
  x[x == ""] <- NA_character_
  value <- suppressWarnings(as.numeric(x))
  
  if (any(!is.na(x) & !is.finite(value))) {
    stop("变量“", v, "”包含非数值内容，请核对Excel。")
  }
  dat[[v]] <- value
}

if (!all(dat$dataset %in% c(0, 1))) {
  stop("dataset存在缺失或非0/1的取值。")
}
if (!all(dat$是否行CRRT %in% c(0, 1))) {
  stop("是否行CRRT存在缺失或非0/1的取值。")
}
if (anyNA(dat$浓度) || any(dat$浓度 < 0)) {
  stop("浓度存在缺失或负值，请核对首次TDM浓度。")
}


# 2. 建立肾功能分组及两种浓度结局 ------------------------------------------

# 沿用原分析：CRRT优先归第5组；非CRRT者按CKDEPI分组。
# >50：1；>25至50：2；10至25：3；<10：4；CRRT：5。
dat$CKDEPI分组 <- with(dat, ifelse(
  是否行CRRT == 1, 5,
  ifelse(CKDEPI < 10, 4,
         ifelse(CKDEPI <= 25, 3,
                ifelse(CKDEPI <= 50, 2, 1)))
))

make_outcome <- function(concentration, upper_limit) {
  group <- ifelse(concentration < 2, 1L,
                  ifelse(concentration <= upper_limit, 2L, 3L))
  ordered(group, levels = 1:3)
}

dat$分组8  <- make_outcome(dat$浓度, upper_limit = 8)
dat$分组16 <- make_outcome(dat$浓度, upper_limit = 16)

predictors <- c("尿素氮", "血红蛋白", "血小板", "年龄",
                "APTT", "输注时间", "CKDEPI分组")

# 防止模型因缺失值自动删行，导致分析人数与论文不一致
na_counts <- colSums(is.na(dat[, predictors, drop = FALSE]))
if (any(na_counts > 0)) {
  print(na_counts[na_counts > 0])
  stop("模型变量存在缺失。请使用与主分析相同的缺失值处理结果后再运行。")
}

train <- dat[dat$dataset == 1, , drop = FALSE]
valid <- dat[dat$dataset == 0, , drop = FALSE]

# 用已知样本数和原三组人数核对文件及队列划分
if (nrow(train) != 327 || nrow(valid) != 71 ||
    !all(as.integer(table(train$分组8)) == c(91, 94, 142))) {
  stop("样本数与原研究不一致。请核对数据文件、dataset编码和浓度单位。")
}


# 3. 统计分组人数和重分类情况 ------------------------------------------------

count_groups <- function(y, cohort, analysis) {
  n <- as.integer(table(y))
  data.frame(
    数据集 = cohort,
    分析 = analysis,
    暴露组 = c("低暴露", "中间暴露", "高暴露"),
    例数 = n,
    百分比 = round(100 * n / sum(n), 2)
  )
}

group_results <- rbind(
  count_groups(train$分组8,  "训练集", "主分析（2和8）"),
  count_groups(train$分组16, "训练集", "敏感性分析（2和16）"),
  count_groups(valid$分组8,  "时间验证集", "主分析（2和8）"),
  count_groups(valid$分组16, "时间验证集", "敏感性分析（2和16）")
)

reclassification <- function(df, cohort) {
  changed <- df$分组8 != df$分组16
  
  # 只有8 < 浓度 ≤ 16的患者会从原第3组转入新第2组
  stopifnot(all(changed == (df$浓度 > 8 & df$浓度 <= 16)))
  
  data.frame(
    数据集 = cohort,
    总人数 = nrow(df),
    高暴露转入中间暴露例数 = sum(changed),
    重分类百分比 = round(100 * mean(changed), 2)
  )
}

reclass_results <- rbind(
  reclassification(train, "训练集"),
  reclassification(valid, "时间验证集")
)

cross_table <- function(df) {
  tab <- with(df, table(原分组 = 分组8, 新分组 = 分组16))
  data.frame(
    原分组 = c("低暴露", "中间暴露", "高暴露"),
    新低暴露 = as.integer(tab[, "1"]),
    新中间暴露 = as.integer(tab[, "2"]),
    新高暴露 = as.integer(tab[, "3"])
  )
}


# 4. 在训练集分别拟合两个有序Logistic回归模型 -------------------------------

fit_ordinal <- function(outcome_name) {
  if (any(table(train[[outcome_name]]) == 0)) {
    stop(outcome_name, "在训练集中有类别为0例，无法拟合三分类模型。")
  }
  
  fit <- MASS::polr(
    formula = reformulate(predictors, response = outcome_name),
    data = train,
    method = "logistic",
    Hess = TRUE,
    na.action = na.fail
  )
  
  if (fit$convergence != 0 || length(coef(fit)) != length(predictors)) {
    stop(outcome_name, "模型未正常收敛或未完整估计7个预测参数。")
  }
  fit
}

model8  <- fit_ordinal("分组8")
model16 <- fit_ordinal("分组16")


# 5. 获取时间验证集预测结果 --------------------------------------------------

validation_predictions <- function(fit, y) {
  prob <- predict(fit, newdata = valid, type = "probs")
  prob <- prob[, c("1", "2", "3"), drop = FALSE]
  
  # 与原代码一致：各类别预测概率 × 类别等级，得到排序分数
  score <- as.numeric(prob %*% c(1, 2, 3))
  
  # 有序模型线性预测值，用于共同有序校准斜率
  beta <- coef(fit)
  x <- model.matrix(delete.response(terms(fit)), data = valid)
  lp <- as.numeric(x[, names(beta), drop = FALSE] %*% beta)
  
  # 多分类Brier：每位患者三个类别的平方误差之和，再对患者取平均
  observed <- cbind(y == "1", y == "2", y == "3")
  
  data.frame(
    y = y,
    score = score,
    lp = lp,
    p_high = prob[, "3"],
    brier = rowSums((prob - observed)^2)
  )
}

pred8  <- validation_predictions(model8,  valid$分组8)
pred16 <- validation_predictions(model16, valid$分组16)


# 6. 定义有序C-index和共同有序校准斜率 -------------------------------------

ordinal_cindex <- function(y, score) {
  # 仅比较实际结局类别不同的患者对；预测得分相同计0.5
  dy <- outer(as.integer(y), as.integer(y), "-")
  ds <- outer(score, score, "-")
  use <- upper.tri(dy) & dy != 0
  
  if (!any(use)) return(NA_real_)
  
  pair_product <- dy[use] * ds[use]
  mean((pair_product > 0) + 0.5 * (pair_product == 0))
}

ordinal_slope <- function(df) {
  # 某次验证集重抽样可能缺少一个类别，此时不能估计三分类校准斜率
  if (any(table(df$y) == 0) || length(unique(df$lp)) < 2) {
    return(NA_real_)
  }
  
  fit <- tryCatch(
    suppressWarnings(
      MASS::polr(y ~ lp, data = df, method = "logistic", Hess = FALSE)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || fit$convergence != 0) return(NA_real_)
  slope <- unname(coef(fit)["lp"])
  if (is.finite(slope)) slope else NA_real_
}


# 7. 计算时间验证集性能及bootstrap 95% CI ----------------------------------
# 两个训练模型固定。每次按验证患者重抽样，两套界值使用同一批行号。

set.seed(2026)
B <- 1000L
boot_rows <- replicate(B, sample.int(nrow(valid), replace = TRUE))

get_ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  unname(quantile(x, c(0.025, 0.975)))
}

summarise_performance <- function(df, analysis) {
  boot_values <- vapply(seq_len(B), function(b) {
    boot_df <- df[boot_rows[, b], , drop = FALSE]
    c(
      C_index = ordinal_cindex(boot_df$y, boot_df$score),
      slope = ordinal_slope(boot_df)
    )
  }, FUN.VALUE = c(C_index = 0, slope = 0))
  
  ci_c <- get_ci(boot_values["C_index", ])
  ci_s <- get_ci(boot_values["slope", ])
  
  data.frame(
    分析 = analysis,
    验证集人数 = nrow(df),
    高暴露例数 = sum(df$y == "3"),
    C_index = ordinal_cindex(df$y, df$score),
    C_index下限 = ci_c[1],
    C_index上限 = ci_c[2],
    共同有序校准斜率 = ordinal_slope(df),
    斜率下限 = ci_s[1],
    斜率上限 = ci_s[2],
    多分类Brier分数 = mean(df$brier),
    高暴露平均预测百分比 = 100 * mean(df$p_high),
    高暴露实际百分比 = 100 * mean(df$y == "3"),
    C_index有效bootstrap次数 = sum(is.finite(boot_values["C_index", ])),
    斜率有效bootstrap次数 = sum(is.finite(boot_values["slope", ]))
  )
}

performance <- rbind(
  summarise_performance(pred8,  "主分析（2和8）"),
  summarise_performance(pred16, "敏感性分析（2和16）")
)


# 8. 打印并导出结果 ----------------------------------------------------------

print(group_results, row.names = FALSE)
print(reclass_results, row.names = FALSE)
print(performance, digits = 4, row.names = FALSE)

output_file <- file.path(work_dir, "美罗培南谷浓度界值敏感性分析_2和8对比2和16.xlsx")

writexl::write_xlsx(
  list(
    分组人数及百分比 = group_results,
    重分类汇总 = reclass_results,
    训练集重分类交叉表 = cross_table(train),
    时间验证集重分类交叉表 = cross_table(valid),
    时间验证集性能对照 = performance
  ),
  path = output_file
)

message("分析完成，结果已保存至：", output_file)
