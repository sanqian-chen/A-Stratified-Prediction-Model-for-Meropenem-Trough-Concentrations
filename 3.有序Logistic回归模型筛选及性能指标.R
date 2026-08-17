library(dplyr)
library(MASS)
library(brant)
library(writexl)
library(ggplot2)
library(tidyr)

# 1. 准备训练集和外部验证集 -------------------------------------------------

train_data <- data %>%
  dplyr::filter(dataset == 1)

valid_data <- data %>%
  dplyr::filter(dataset == 0)

# 因变量：有序三分类
train_data$血药浓度分组 <- ordered(
  train_data$血药浓度分组,
  levels = c(1, 2, 3)
)

valid_data$血药浓度分组 <- ordered(
  valid_data$血药浓度分组,
  levels = c(1, 2, 3)
)

# CKDEPI分组作为数值型有序变量
train_data$CKDEPI分组 <- as.numeric(as.character(train_data$CKDEPI分组))
valid_data$CKDEPI分组 <- as.numeric(as.character(valid_data$CKDEPI分组))

# 三个等级
class_levels <- levels(train_data$血药浓度分组)

# 空模型预测概率，使用训练集类别构成作为基准概率
null_probs <- prop.table(table(train_data$血药浓度分组))
null_probs <- as.numeric(null_probs[class_levels])
names(null_probs) <- class_levels


# 2. 构建公式函数 -----------------------------------------------------------

make_formula <- function(vars) {
  as.formula(
    paste0(
      "`血药浓度分组` ~ ",
      paste(paste0("`", vars, "`"), collapse = " + ")
    )
  )
}


# 3. 数值保留4位小数函数 ----------------------------------------------------

round_df <- function(df, digits = 4) {
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) round(x, digits) else x
  })
  return(df)
}


# 4. 提取有序Logistic回归结果，包括截距 ------------------------------------

extract_polr_result <- function(model) {
  
  ctable <- coef(summary(model))
  
  p_value <- pnorm(
    abs(ctable[, "t value"]),
    lower.tail = FALSE
  ) * 2
  
  beta <- ctable[, "Value"]
  se <- ctable[, "Std. Error"]
  
  result <- data.frame(
    变量 = rownames(ctable),
    Value = beta,
    Std_Error = se,
    t_value = ctable[, "t value"],
    P值 = p_value,
    OR = exp(beta),
    OR_95CI下限 = exp(beta - 1.96 * se),
    OR_95CI上限 = exp(beta + 1.96 * se),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  result <- round_df(result, 4)
  
  return(result)
}


# 5. 提取自变量P值，不包括截距 ---------------------------------------------

extract_predictor_p <- function(model) {
  
  ctable <- coef(summary(model))
  
  coef_table <- ctable[
    !grepl("\\|", rownames(ctable)),
    ,
    drop = FALSE
  ]
  
  p_value <- pnorm(
    abs(coef_table[, "t value"]),
    lower.tail = FALSE
  ) * 2
  
  p_df <- data.frame(
    变量 = rownames(coef_table),
    P值 = p_value,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  return(p_df)
}


# 6. 计算模型整体似然比检验P值 ---------------------------------------------

get_model_global_p <- function(model, null_model) {
  
  LR <- 2 * (as.numeric(logLik(model)) - as.numeric(logLik(null_model)))
  
  df <- attr(logLik(model), "df") - attr(logLik(null_model), "df")
  
  p <- pchisq(LR, df = df, lower.tail = FALSE)
  
  return(
    data.frame(
      LR_chisq = LR,
      LR_df = df,
      模型整体P值 = p
    )
  )
}


# 7. 平行线假设检验：Brant test --------------------------------------------

extract_brant_result <- function(model) {
  
  br <- tryCatch(
    {
      tmp <- capture.output(
        res <- suppressWarnings(brant(model))
      )
      res
    },
    error = function(e) {
      return(e)
    }
  )
  
  if (inherits(br, "error")) {
    br_df <- data.frame(
      变量 = "Brant检验失败",
      X2 = NA,
      df = NA,
      probability = NA,
      说明 = br$message,
      stringsAsFactors = FALSE
    )
  } else {
    br_df <- as.data.frame(br)
    br_df$变量 <- rownames(br_df)
    br_df <- br_df[, c("变量", setdiff(names(br_df), "变量"))]
  }
  
  br_df <- round_df(br_df, 4)
  
  return(br_df)
}


# 8. 提取预测概率 -----------------------------------------------------------

get_pred_prob <- function(model, newdata) {
  
  pred_prob <- predict(
    model,
    newdata = newdata,
    type = "probs"
  )
  
  if (is.null(dim(pred_prob))) {
    pred_prob <- matrix(pred_prob, nrow = 1)
    colnames(pred_prob) <- levels(newdata$血药浓度分组)
  }
  
  pred_prob <- as.data.frame(pred_prob)
  pred_prob <- pred_prob[, class_levels, drop = FALSE]
  
  return(as.matrix(pred_prob))
}


# 9. 有序结局C-index函数 ----------------------------------------------------

ordinal_cindex <- function(y, score) {
  
  y <- as.numeric(y)
  score <- as.numeric(score)
  
  dy <- outer(y, y, "-")
  ds <- outer(score, score, "-")
  
  use <- upper.tri(dy) & dy != 0
  prod <- dy[use] * ds[use]
  
  n_pairs <- length(prod)
  
  if (n_pairs == 0) {
    return(
      data.frame(
        C_index = NA,
        Dxy = NA,
        可比较对子数 = 0
      )
    )
  }
  
  c_index <- (
    sum(prod > 0) +
      0.5 * sum(prod == 0)
  ) / n_pairs
  
  dxy <- 2 * (c_index - 0.5)
  
  return(
    data.frame(
      C_index = c_index,
      Dxy = dxy,
      可比较对子数 = n_pairs
    )
  )
}


# 10. Multiclass Brier score 和 RPS -----------------------------------------

calc_brier_rps <- function(y_true, pred_prob, null_probs) {
  
  classes <- colnames(pred_prob)
  K <- length(classes)
  
  y_char <- as.character(y_true)
  y_fac <- factor(y_char, levels = classes)
  
  # one-hot矩阵
  y_mat <- model.matrix(~ 0 + y_fac)
  colnames(y_mat) <- classes
  
  # Multiclass Brier score
  brier <- mean(
    rowSums((pred_prob - y_mat)^2)
  )
  
  # 空模型Brier
  null_prob_mat <- matrix(
    rep(null_probs, each = nrow(pred_prob)),
    nrow = nrow(pred_prob),
    byrow = FALSE
  )
  colnames(null_prob_mat) <- classes
  
  null_brier <- mean(
    rowSums((null_prob_mat - y_mat)^2)
  )
  
  scaled_brier <- 1 - brier / null_brier
  
  # RPS：有序多分类概率评分
  pred_cum <- t(
    apply(pred_prob[, 1:(K - 1), drop = FALSE], 1, cumsum)
  )
  
  y_num <- as.numeric(as.character(y_true))
  
  obs_cum <- sapply(1:(K - 1), function(k) {
    as.numeric(y_num <= k)
  })
  
  obs_cum <- as.matrix(obs_cum)
  
  rps_raw <- mean(
    rowSums((pred_cum - obs_cum)^2)
  )
  
  rps_std <- rps_raw / (K - 1)
  
  # 空模型RPS
  null_cum <- cumsum(null_probs)[1:(K - 1)]
  
  null_cum_mat <- matrix(
    rep(null_cum, each = nrow(pred_prob)),
    nrow = nrow(pred_prob),
    byrow = FALSE
  )
  
  null_rps_raw <- mean(
    rowSums((null_cum_mat - obs_cum)^2)
  )
  
  null_rps_std <- null_rps_raw / (K - 1)
  
  scaled_rps <- 1 - rps_raw / null_rps_raw
  
  out <- data.frame(
    Multiclass_Brier = brier,
    Null_Brier = null_brier,
    Scaled_Brier = scaled_brier,
    RPS = rps_raw,
    Standardized_RPS = rps_std,
    Null_RPS = null_rps_raw,
    Null_Standardized_RPS = null_rps_std,
    Scaled_RPS = scaled_rps
  )
  
  return(out)
}


# 11. 计算整体预测性能 ------------------------------------------------------

get_prediction_performance <- function(model, newdata, dataset_name, null_probs) {
  
  if (nrow(newdata) == 0) {
    return(
      data.frame(
        数据集 = dataset_name,
        N = 0,
        Accuracy = NA,
        Misclassification = NA,
        C_index = NA,
        Dxy = NA,
        可比较对子数 = NA,
        Multiclass_Brier = NA,
        Null_Brier = NA,
        Scaled_Brier = NA,
        RPS = NA,
        Standardized_RPS = NA,
        Null_RPS = NA,
        Null_Standardized_RPS = NA,
        Scaled_RPS = NA
      )
    )
  }
  
  pred_prob <- get_pred_prob(model, newdata)
  
  pred_class <- predict(
    model,
    newdata = newdata,
    type = "class"
  )
  
  y_true <- newdata$血药浓度分组
  
  classes <- colnames(pred_prob)
  
  # 预测期望等级，作为有序风险分数
  pred_score <- as.numeric(
    pred_prob %*% as.numeric(classes)
  )
  
  accuracy <- mean(
    as.character(pred_class) == as.character(y_true)
  )
  
  misclass <- 1 - accuracy
  
  cindex_res <- ordinal_cindex(
    y = as.numeric(as.character(y_true)),
    score = pred_score
  )
  
  score_res <- calc_brier_rps(
    y_true = y_true,
    pred_prob = pred_prob,
    null_probs = null_probs
  )
  
  out <- data.frame(
    数据集 = dataset_name,
    N = nrow(newdata),
    Accuracy = accuracy,
    Misclassification = misclass,
    C_index = cindex_res$C_index,
    Dxy = cindex_res$Dxy,
    可比较对子数 = cindex_res$可比较对子数,
    score_res
  )
  
  out <- round_df(out, 4)
  
  return(out)
}


# 12. 平滑分类别校准曲线函数 ------------------------------------------------
# 只绘制 P(Y=1)、P(Y=2)、P(Y=3) 的分类别校准曲线
# 删除累积概率校准曲线

smooth_curve <- function(prob, event, grid, span = 0.8) {
  
  df <- data.frame(
    prob = as.numeric(prob),
    event = as.numeric(event)
  ) %>%
    dplyr::filter(!is.na(prob), !is.na(event))
  
  if (
    nrow(df) < 10 ||
    length(unique(df$prob)) < 4 ||
    length(unique(df$event)) < 2
  ) {
    return(rep(mean(df$event), length(grid)))
  }
  
  fit <- tryCatch(
    loess(
      event ~ prob,
      data = df,
      span = span,
      degree = 1,
      control = loess.control(surface = "direct")
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(rep(mean(df$event), length(grid)))
  }
  
  pred <- as.numeric(
    predict(
      fit,
      newdata = data.frame(prob = grid)
    )
  )
  
  if (any(is.na(pred))) {
    
    not_na <- !is.na(pred)
    
    if (sum(not_na) >= 2) {
      pred <- approx(
        x = grid[not_na],
        y = pred[not_na],
        xout = grid,
        rule = 2
      )$y
    } else {
      pred[is.na(pred)] <- mean(df$event)
    }
  }
  
  pred <- pmin(pmax(pred, 0), 1)
  
  return(pred)
}


calibration_smooth_by_class <- function(model, newdata, model_name, dataset_name,
                                        span = 0.8, n_grid = 100) {
  
  pred_prob <- get_pred_prob(model, newdata)
  y_true <- as.character(newdata$血药浓度分组)
  
  out_list <- list()
  
  for (cls in colnames(pred_prob)) {
    
    prob_cls <- pred_prob[, cls]
    event_cls <- as.numeric(y_true == cls)
    
    grid_min <- max(0.01, min(prob_cls, na.rm = TRUE))
    grid_max <- min(0.99, max(prob_cls, na.rm = TRUE))
    
    if (grid_min >= grid_max) {
      grid_min <- max(0.01, grid_min - 0.01)
      grid_max <- min(0.99, grid_max + 0.01)
    }
    
    grid_cls <- seq(
      from = grid_min,
      to = grid_max,
      length.out = n_grid
    )
    
    observed_smooth <- smooth_curve(
      prob = prob_cls,
      event = event_cls,
      grid = grid_cls,
      span = span
    )
    
    out_list[[cls]] <- data.frame(
      模型 = model_name,
      数据集 = dataset_name,
      校准类型 = "分类别平滑校准",
      类别 = paste0("P(Y=", cls, ")"),
      mean_pred = grid_cls,
      observed = observed_smooth
    )
  }
  
  out <- dplyr::bind_rows(out_list)
  
  return(out)
}


# 13. 非目标浓度DCA函数 -----------------------------------------------------
# 目标浓度 = Y=1
# 非目标浓度 = Y=2 或 Y=3
# DCA事件定义：非目标浓度，即 Y=2 或 Y=3
# 模型预测概率：P(Y=2) + P(Y=3)

dca_nontarget <- function(model, newdata, model_name, dataset_name,
                          thresholds = seq(0.01, 0.80, by = 0.01)) {
  
  pred_prob <- get_pred_prob(model, newdata)
  
  # 非目标浓度 = Y=2 或 Y=3
  prob_nontarget <- pred_prob[, "2"] + pred_prob[, "3"]
  
  y_true <- as.character(newdata$血药浓度分组)
  
  # 事件 = 非目标浓度，即Y=2或Y=3
  event_nontarget <- as.numeric(y_true %in% c("2", "3"))
  
  df <- data.frame(
    y = event_nontarget,
    p = prob_nontarget
  )
  
  N <- nrow(df)
  prevalence <- mean(df$y == 1)
  
  dca_model <- lapply(thresholds, function(pt) {
    
    pred_pos <- df$p >= pt
    
    TP <- sum(pred_pos & df$y == 1)
    FP <- sum(pred_pos & df$y == 0)
    
    net_benefit <- TP / N - FP / N * (pt / (1 - pt))
    
    data.frame(
      threshold = pt,
      net_benefit = net_benefit,
      模型 = model_name,
      数据集 = dataset_name
    )
  })
  
  dca_model <- do.call(rbind, dca_model)
  
  dca_all <- data.frame(
    threshold = thresholds,
    net_benefit = prevalence - (1 - prevalence) * (thresholds / (1 - thresholds)),
    模型 = "All",
    数据集 = dataset_name
  )
  
  dca_none <- data.frame(
    threshold = thresholds,
    net_benefit = 0,
    模型 = "None",
    数据集 = dataset_name
  )
  
  out <- bind_rows(
    dca_model,
    dca_all,
    dca_none
  )
  
  out <- round_df(out, 4)
  
  return(out)
}


# 14. 初始9变量模型 ---------------------------------------------------------

vars_current <- c(
  "尿素氮",
  "年龄",
  "血红蛋白",
  "血小板",
  "APACHEII",
  "APTT",
  "负平衡量",
  "输注时间",
  "CKDEPI分组"
)

alpha <- 0.05

model_list <- list()
result_list <- list()
brant_list <- list()
compare_list <- list()
performance_list <- list()

# 空模型，用于整体似然比检验
fit_null <- polr(
  `血药浓度分组` ~ 1,
  data = train_data,
  Hess = TRUE,
  method = "logistic"
)

step_id <- 1


# 15. 从9变量模型开始，逐个删除P值最大的非显著变量 -------------------------

repeat {
  
  formula_current <- make_formula(vars_current)
  
  fit_current <- polr(
    formula_current,
    data = train_data,
    Hess = TRUE,
    method = "logistic"
  )
  
  model_name <- paste0("模型", step_id, "_", length(vars_current), "变量")
  
  model_list[[model_name]] <- fit_current
  
  result_list[[model_name]] <- extract_polr_result(fit_current)
  
  p_current <- extract_predictor_p(fit_current)
  
  max_p <- max(p_current$P值)
  remove_var <- p_current$变量[which.max(p_current$P值)]
  
  global_p <- get_model_global_p(fit_current, fit_null)
  
  brant_current <- extract_brant_result(fit_current)
  brant_list[[paste0("Brant_", model_name)]] <- brant_current
  
  brant_omnibus_p <- brant_current$probability[
    brant_current$变量 == "Omnibus"
  ]
  
  if (length(brant_omnibus_p) == 0) {
    brant_omnibus_p <- NA
  }
  
  perf_train <- get_prediction_performance(
    model = fit_current,
    newdata = train_data,
    dataset_name = "训练集",
    null_probs = null_probs
  )
  
  perf_valid <- get_prediction_performance(
    model = fit_current,
    newdata = valid_data,
    dataset_name = "外部验证集",
    null_probs = null_probs
  )
  
  perf_current <- rbind(perf_train, perf_valid)
  perf_current$模型 <- model_name
  perf_current <- perf_current[
    ,
    c("模型", setdiff(names(perf_current), "模型"))
  ]
  
  performance_list[[model_name]] <- perf_current
  
  compare_current <- data.frame(
    模型 = model_name,
    变量数 = length(vars_current),
    当前变量 = paste(vars_current, collapse = " + "),
    AIC = AIC(fit_current),
    BIC = BIC(fit_current),
    logLik = as.numeric(logLik(fit_current)),
    LR_chisq = global_p$LR_chisq,
    LR_df = global_p$LR_df,
    模型整体P值 = global_p$模型整体P值,
    Brant_Omnibus_P = brant_omnibus_p,
    最大P值 = max_p,
    最大P值变量 = remove_var,
    stringsAsFactors = FALSE
  )
  
  compare_list[[model_name]] <- compare_current
  
  if (max_p <= alpha) {
    break
  }
  
  if (length(vars_current) <= 1) {
    break
  }
  
  vars_current <- setdiff(vars_current, remove_var)
  
  step_id <- step_id + 1
}


# 16. 合并模型比较和预测性能结果 -------------------------------------------

model_compare <- do.call(rbind, compare_list)
model_compare <- round_df(model_compare, 4)

performance_all <- do.call(rbind, performance_list)
performance_all <- round_df(performance_all, 4)

model_compare
performance_all


# 17. 计算所有模型的平滑分类别校准曲线和非目标浓度DCA -----------------------

cal_class_list <- list()
dca_list <- list()

thresholds <- seq(0.01, 0.80, by = 0.01)

for (model_name in names(model_list)) {
  
  model_i <- model_list[[model_name]]
  
  # 平滑分类别校准：训练集
  cal_class_list[[paste0(model_name, "_训练集_分类别平滑")]] <-
    calibration_smooth_by_class(
      model = model_i,
      newdata = train_data,
      model_name = model_name,
      dataset_name = "训练集",
      span = 0.8,
      n_grid = 100
    )
  
  # 平滑分类别校准：外部验证集
  cal_class_list[[paste0(model_name, "_外部验证集_分类别平滑")]] <-
    calibration_smooth_by_class(
      model = model_i,
      newdata = valid_data,
      model_name = model_name,
      dataset_name = "外部验证集",
      span = 1,
      n_grid = 100
    )
  
  # 非目标浓度DCA：训练集
  dca_list[[paste0(model_name, "_训练集_DCA")]] <-
    dca_nontarget(
      model = model_i,
      newdata = train_data,
      model_name = model_name,
      dataset_name = "训练集",
      thresholds = thresholds
    )
  
  # 非目标浓度DCA：外部验证集
  dca_list[[paste0(model_name, "_外部验证集_DCA")]] <-
    dca_nontarget(
      model = model_i,
      newdata = valid_data,
      model_name = model_name,
      dataset_name = "外部验证集",
      thresholds = thresholds
    )
}

cal_class_all <- bind_rows(cal_class_list)
dca_all_raw <- bind_rows(dca_list)


# All / None 每个数据集中只保留一次 ---------------------------------------

dca_models <- dca_all_raw %>%
  filter(!模型 %in% c("All", "None"))

dca_reference <- dca_all_raw %>%
  filter(模型 %in% c("All", "None")) %>%
  distinct(数据集, threshold, 模型, .keep_all = TRUE)

dca_all <- bind_rows(
  dca_models,
  dca_reference
) %>%
  mutate(
    数据集 = factor(
      数据集,
      levels = c("训练集", "外部验证集")
    )
  )

dca_all <- round_df(dca_all, 4)


# 18. 绘制平滑分类别校准曲线 ------------------------------------------------
# 不绘制累积概率校准曲线

p_cal_class_train <- cal_class_all %>%
  filter(数据集 == "训练集") %>%
  ggplot(
    aes(
      x = mean_pred,
      y = observed,
      color = 模型,
      linetype = 模型,
      group = 模型
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.8
  ) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ 类别, nrow = 1) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    title = "Smoothed Category-specific Calibration Curves in Training Set",
    x = "Predicted probability",
    y = "Observed proportion"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

p_cal_class_valid <- cal_class_all %>%
  filter(数据集 == "外部验证集") %>%
  ggplot(
    aes(
      x = mean_pred,
      y = observed,
      color = 模型,
      linetype = 模型,
      group = 模型
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.8
  ) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ 类别, nrow = 1) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    title = "Smoothed Category-specific Calibration Curves in External Validation Set",
    x = "Predicted probability",
    y = "Observed proportion"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

p_cal_class_train
p_cal_class_valid


# 19. 绘制非目标浓度DCA -----------------------------------------------------
# 目标浓度 = Y=1
# 非目标浓度 = Y=2 或 Y=3
# DCA事件 = 非目标浓度

# 19.1 训练集DCA -----------------------------------------------------------

p_dca_train <- dca_all %>%
  filter(数据集 == "训练集") %>%
  filter(threshold >= 0.01, threshold <= 0.80) %>%
  ggplot(
    aes(
      x = threshold,
      y = net_benefit,
      color = 模型,
      linetype = 模型
    )
  ) +
  geom_line(linewidth = 1.0) +
  scale_x_continuous(
    limits = c(0, 0.80),
    breaks = seq(0, 0.80, 0.10)
  ) +
  scale_color_manual(
    values = c(
      "All" = "#F8766D",
      "None" = "#00BA38",
      "模型1_9变量" = "#00BFC4",
      "模型2_8变量" = "#619CFF",
      "模型3_7变量" = "#C77CFF"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "All" = "solid",
      "None" = "dashed",
      "模型1_9变量" = "dashed",
      "模型2_8变量" = "dashed",
      "模型3_7变量" = "dotted"
    )
  ) +
  labs(
    title = "Decision Curve Analysis for Non-target Concentration in Training Set",
    subtitle = "Target concentration: Y=1; Non-target concentration: Y=2 or Y=3",
    x = "Threshold probability",
    y = "Net benefit"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

p_dca_train


# 19.2 外部验证集DCA -------------------------------------------------------

p_dca_valid <- dca_all %>%
  filter(数据集 == "外部验证集") %>%
  filter(threshold >= 0.01, threshold <= 0.80) %>%
  ggplot(
    aes(
      x = threshold,
      y = net_benefit,
      color = 模型,
      linetype = 模型
    )
  ) +
  geom_line(linewidth = 1.0) +
  scale_x_continuous(
    limits = c(0, 0.80),
    breaks = seq(0, 0.80, 0.10)
  ) +
  scale_color_manual(
    values = c(
      "All" = "#F8766D",
      "None" = "#00BA38",
      "模型1_9变量" = "#00BFC4",
      "模型2_8变量" = "#619CFF",
      "模型3_7变量" = "#C77CFF"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "All" = "solid",
      "None" = "dashed",
      "模型1_9变量" = "dashed",
      "模型2_8变量" = "dashed",
      "模型3_7变量" = "dotted"
    )
  ) +
  labs(
    title = "Decision Curve Analysis for Non-target Concentration in External Validation Set",
    subtitle = "Target concentration: Y=1; Non-target concentration: Y=2 or Y=3",
    x = "Threshold probability",
    y = "Net benefit"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

p_dca_valid


# 19.3 合并绘制训练集和外部验证集DCA --------------------------------------
# 训练集在左，外部验证集在右

p_dca_combined <- dca_all %>%
  filter(threshold >= 0.01, threshold <= 0.80) %>%
  mutate(
    数据集 = factor(
      数据集,
      levels = c("训练集", "外部验证集")
    )
  ) %>%
  ggplot(
    aes(
      x = threshold,
      y = net_benefit,
      color = 模型,
      linetype = 模型
    )
  ) +
  geom_line(linewidth = 1.0) +
  facet_wrap(
    ~ 数据集,
    nrow = 1
  ) +
  scale_x_continuous(
    limits = c(0, 0.80),
    breaks = seq(0, 0.80, 0.10)
  ) +
  scale_color_manual(
    values = c(
      "All" = "#F8766D",
      "None" = "#00BA38",
      "模型1_9变量" = "#00BFC4",
      "模型2_8变量" = "#619CFF",
      "模型3_7变量" = "#C77CFF"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "All" = "solid",
      "None" = "dashed",
      "模型1_9变量" = "dashed",
      "模型2_8变量" = "dashed",
      "模型3_7变量" = "dotted"
    )
  ) +
  labs(
    title = "Decision Curve Analysis for Predicting Non-target Concentration",
    subtitle = "Target concentration: Y=1; Non-target concentration: Y=2 or Y=3",
    x = "Threshold probability",
    y = "Net benefit"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

p_dca_combined


# 20. 查看最终模型 ----------------------------------------------------------

final_model_name <- names(model_list)[length(model_list)]
final_model <- model_list[[final_model_name]]

summary(final_model)
formula(final_model)

final_result <- result_list[[final_model_name]]
final_result

final_brant <- brant_list[[paste0("Brant_", final_model_name)]]
final_brant






# 21. 多重共线性诊断（VIF） ---------------------------------------------
#install.packages("car")
library(car)

# 获取最终模型公式
final_formula <- formula(final_model)

# 提取自变量名称
predictor_names <- attr(
  terms(final_formula),
  "term.labels"
)

predictor_names

# 将有序因子转换为数值型（1,2,3...）
train_data$血药浓度分组_num <- as.numeric(train_data$血药浓度分组)

# 构建VIF辅助模型
vif_formula <- as.formula(
  paste("血药浓度分组_num ~", paste(predictor_names, collapse = " + "))
)
vif_model <- lm(vif_formula, data = train_data)

# 计算VIF
vif_result <- vif(vif_model)
vif_result


vif_table <- data.frame(
  变量 = names(vif_result),
  VIF = round(vif_result, 3)
)
vif_table


# 22. 导出Excel -------------------------------------------------------------

export_list <- c(
  result_list,
  brant_list,
  list(
    模型比较 = model_compare,
    整体预测性能 = performance_all,
    平滑分类别校准曲线数据 = cal_class_all,
    非目标浓度DCA数据 = dca_all
  )
)

write_xlsx(
  export_list,
  path = "有序三分类模型_整体评价_平滑分类别校准曲线_DCA.xlsx"
)


# 23. 保存图片 --------------------------------------------------------------

ggsave(
  "训练集_平滑分类别校准曲线.png",
  p_cal_class_train,
  width = 10,
  height = 4,
  dpi = 300
)

ggsave(
  "外部验证集_平滑分类别校准曲线.png",
  p_cal_class_valid,
  width = 10,
  height = 4,
  dpi = 300
)

ggsave(
  "训练集_非目标浓度DCA.png",
  p_dca_train,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "外部验证集_非目标浓度DCA.png",
  p_dca_valid,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "训练集和外部验证集_非目标浓度DCA.png",
  p_dca_combined,
  width = 11,
  height = 5,
  dpi = 300
)