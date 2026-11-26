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


# 18. 绘制平滑分类别校准曲线 -----------------------------------------------

library(dplyr)
library(ggplot2)

# 如尚未安装字体包，首次运行时取消下一行的注释
# install.packages(c("sysfonts", "showtext"))

# 18.1. 注册 Windows 系统中的 Times New Roman 字体 ----------------------------

font_file <- file.path(
  Sys.getenv("WINDIR", "C:/Windows"),
  "Fonts",
  "times.ttf"
)

if (!file.exists(font_file)) {
  stop("未找到 Times New Roman，请将 font_file 改为本机 times.ttf 的实际路径。")
}

font_family <- "Times New Roman"
sysfonts::font_add(family = font_family, regular = font_file)
showtext::showtext_auto(enable = TRUE)

# 18.2. 设置分面标题、模型颜色、线型及图例文字 -------------------------------

class_labels <- c(
  "P(Y=1)" = "Low",
  "P(Y=2)" = "Intermediate",
  "P(Y=3)" = "High"
)

model_colors <- c(
  "模型1_9变量" = "#F8766D",
  "模型2_8变量" = "#00BA38",
  "模型3_7变量" = "#619CFF"
)

model_linetypes <- c(
  "模型1_9变量" = "solid",
  "模型2_8变量" = "dashed",
  "模型3_7变量" = "dotdash"
)

model_labels <- c(
  "模型1_9变量" = "9-predictor model",
  "模型2_8变量" = "8-predictor model",
  "模型3_7变量" = "7-predictor model"
)

# 18.3. 定义训练集和外部验证集共用的绘图函数 ---------------------------------

plot_class_calibration <- function(data, dataset_name) {
  
  data %>%
    filter(数据集 == dataset_name) %>%
    ggplot(
      aes(
        x = mean_pred,
        y = observed,
        color = 模型,
        linetype = 模型,
        group = 模型
      )
    ) +
    # 理想校准线
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.8
    ) +
    # 三个模型的校准曲线
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    facet_wrap(
      ~ 类别,
      nrow = 1,
      labeller = as_labeller(class_labels)
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    scale_color_manual(
      values = model_colors,
      breaks = names(model_labels),
      labels = unname(model_labels)
    ) +
    scale_linetype_manual(
      values = model_linetypes,
      breaks = names(model_labels),
      labels = unname(model_labels)
    ) +
    labs(
      x = "Predicted probability",
      y = "Observed proportion",
      color = NULL,
      linetype = NULL
    ) +
    # 全图文字统一继承 Times New Roman
    theme_bw(base_size = 16, base_family = font_family) +
    theme(
      # 横纵坐标标题
      axis.title.x = element_text(
        size = 18, color = "black", margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 18, color = "black", margin = margin(r = 8)
      ),
      
      # 横纵坐标刻度数字
      axis.text = element_text(size = 16, color = "black"),
      
      # 顶部分面标题 Low、Intermediate、High
      strip.text.x = element_text(
        size = 16, face = "plain", color = "black",
        margin = margin(t = 5, b = 5)
      ),
      strip.background = element_rect(
        fill = "grey85", color = "grey30", linewidth = 0.8
      ),
      
      # 右侧模型图例
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 16, color = "black"),
      legend.key.width = grid::unit(1.2, "cm"),
      
      # 网格线、面板边框及间距
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black", fill = NA, linewidth = 0.8
      ),
      panel.spacing.x = grid::unit(0.35, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# 18.4. 生成训练集和外部验证集的校准图 ---------------------------------------

p_cal_class_train <- plot_class_calibration(cal_class_all, "训练集")
p_cal_class_valid <- plot_class_calibration(cal_class_all, "外部验证集")

# 显示训练集图形
print(p_cal_class_train)

# 查看外部验证集图形时，运行下一行
# print(p_cal_class_valid)
p_cal_class_valid






# 21. 训练集：9、8、7变量模型的分类别DCA曲线 -------------------------------

# 每个类别均采用 one-vs-rest：
# Y=1：事件为Y=1，非事件为Y=2或Y=3
# Y=2：事件为Y=2，非事件为Y=1或Y=3
# Y=3：事件为Y=3，非事件为Y=1或Y=2


# 21.1 计算单个模型的分类别DCA --------------------------------------------

calculate_dca_by_class <- function(
    model,
    newdata,
    model_name,
    dataset_name = "训练集",
    thresholds = seq(0.01, 0.80, by = 0.01)
) {
  
  pred_prob <- get_pred_prob(
    model = model,
    newdata = newdata
  )
  
  y_true <- as.character(newdata$血药浓度分组)
  
  out_list <- list()
  
  for (cls in class_levels) {
    
    # 当前类别的预测概率
    prob_cls <- as.numeric(pred_prob[, cls])
    
    # one-vs-rest事件定义
    event_cls <- as.numeric(y_true == cls)
    
    n <- length(event_cls)
    prevalence <- mean(event_cls)
    
    dca_cls <- lapply(thresholds, function(pt) {
      
      pred_positive <- prob_cls >= pt
      
      TP <- sum(pred_positive & event_cls == 1)
      FP <- sum(pred_positive & event_cls == 0)
      FN <- sum(!pred_positive & event_cls == 1)
      TN <- sum(!pred_positive & event_cls == 0)
      
      # 当前模型的净获益
      nb_model <- TP / n -
        FP / n * (pt / (1 - pt))
      
      # Treat-all策略的净获益
      nb_all <- prevalence -
        (1 - prevalence) * (pt / (1 - pt))
      
      # Treat-none策略的净获益
      nb_none <- 0
      
      data.frame(
        数据集 = dataset_name,
        类别 = paste0("P(Y=", cls, ")"),
        threshold = pt,
        prevalence = prevalence,
        TP = TP,
        FP = FP,
        FN = FN,
        TN = TN,
        策略 = model_name,
        net_benefit = nb_model,
        stringsAsFactors = FALSE
      )
    }) %>%
      dplyr::bind_rows()
    
    # 添加Treat-all和Treat-none参考策略
    dca_reference <- data.frame(
      数据集 = dataset_name,
      类别 = paste0("P(Y=", cls, ")"),
      threshold = thresholds,
      prevalence = prevalence,
      TP = NA,
      FP = NA,
      FN = NA,
      TN = NA,
      策略 = rep(
        c("Treat all", "Treat none"),
        each = length(thresholds)
      ),
      net_benefit = c(
        prevalence -
          (1 - prevalence) *
          (thresholds / (1 - thresholds)),
        rep(0, length(thresholds))
      ),
      stringsAsFactors = FALSE
    )
    
    out_list[[cls]] <- dplyr::bind_rows(
      dca_cls,
      dca_reference
    )
  }
  
  dplyr::bind_rows(out_list)
}


# 21.2 分别计算9、8、7变量模型 ---------------------------------------------

dca_by_class_train_list <- list()

for (model_name in names(model_list)) {
  
  dca_by_class_train_list[[model_name]] <-
    calculate_dca_by_class(
      model = model_list[[model_name]],
      newdata = train_data,
      model_name = model_name,
      dataset_name = "训练集",
      thresholds = seq(0.01, 0.80, by = 0.01)
    )
}

dca_by_class_train_raw <- dplyr::bind_rows(
  dca_by_class_train_list
)


# Treat-all和Treat-none在每个模型中相同，只保留一次 -------------------------

dca_model_curves <- dca_by_class_train_raw %>%
  dplyr::filter(
    !策略 %in% c("Treat all", "Treat none")
  )

dca_reference_curves <- dca_by_class_train_raw %>%
  dplyr::filter(
    策略 %in% c("Treat all", "Treat none")
  ) %>%
  dplyr::distinct(
    数据集,
    类别,
    threshold,
    策略,
    .keep_all = TRUE
  )

dca_by_class_train <- dplyr::bind_rows(
  dca_reference_curves,
  dca_model_curves
) %>%
  dplyr::mutate(
    类别 = factor(
      类别,
      levels = c(
        "P(Y=1)",
        "P(Y=2)",
        "P(Y=3)"
      )
    ),
    策略 = factor(
      策略,
      levels = c(
        "Treat all",
        "Treat none",
        "模型1_9变量",
        "模型2_8变量",
        "模型3_7变量"
      )
    )
  )


# 21.3 对3条模型曲线进行LOESS平滑 -----------------------------------------

# 仅对9、8、7变量模型曲线进行平滑；
# Treat-all和Treat-none参考线保持原始形式。
# 平滑仅用于图形展示，不改变原始净获益结果。

smooth_dca_curves <- function(dca_data, span = 0.25) {
  
  dca_data %>%
    dplyr::group_by(
      数据集,
      类别,
      策略
    ) %>%
    dplyr::group_modify(function(.x, .y) {
      
      df <- .x %>%
        dplyr::arrange(threshold)
      
      strategy_name <- as.character(.y$策略)
      
      # 参考策略不进行平滑
      if (strategy_name %in% c(
        "Treat all",
        "Treat none"
      )) {
        df$net_benefit_plot <- df$net_benefit
        return(df)
      }
      
      # 数据点不足时不平滑
      if (
        nrow(df) < 10 ||
        length(unique(df$threshold)) < 4 ||
        length(unique(df$net_benefit)) < 4
      ) {
        df$net_benefit_plot <- df$net_benefit
        return(df)
      }
      
      fit <- tryCatch(
        loess(
          net_benefit ~ threshold,
          data = df,
          span = span,
          degree = 1,
          control = loess.control(
            surface = "direct"
          )
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit)) {
        df$net_benefit_plot <- df$net_benefit
        return(df)
      }
      
      pred <- as.numeric(
        predict(
          fit,
          newdata = data.frame(
            threshold = df$threshold
          )
        )
      )
      
      # 平滑结果存在NA时进行插值
      if (anyNA(pred)) {
        
        not_na <- !is.na(pred)
        
        if (sum(not_na) >= 2) {
          pred <- approx(
            x = df$threshold[not_na],
            y = pred[not_na],
            xout = df$threshold,
            rule = 2
          )$y
        } else {
          pred <- df$net_benefit
        }
      }
      
      df$net_benefit_plot <- pred
      
      df
    }) %>%
    dplyr::ungroup()
}


dca_by_class_train_smooth <- smooth_dca_curves(
  dca_data = dca_by_class_train,
  span = 0.25
)


# 21.4 设置字体、分面名称及图例 --------------------------------------------

library(ggplot2)

# 如尚未安装字体包，首次运行时取消下一行的注释
# install.packages(c("sysfonts", "showtext"))

# 注册 Windows 系统中的 Times New Roman 字体
font_family <- "Times New Roman"
font_file <- file.path(
  Sys.getenv("WINDIR", "C:/Windows"),
  "Fonts",
  "times.ttf"
)

if (!file.exists(font_file)) {
  stop("未找到 Times New Roman，请将 font_file 改为本机 times.ttf 的实际路径。")
}

sysfonts::font_add(family = font_family, regular = font_file)
showtext::showtext_auto(enable = TRUE)

# 分面标题
class_labels <- c(
  "P(Y=1)" = "Low",
  "P(Y=2)" = "Intermediate",
  "P(Y=3)" = "High"
)

# 图例文字及排列顺序
strategy_labels <- c(
  "Treat all"   = "Treat all",
  "Treat none"  = "Treat none",
  "模型1_9变量" = "9-predictor model",
  "模型2_8变量" = "8-predictor model",
  "模型3_7变量" = "7-predictor model"
)

# 各策略的颜色
strategy_colors <- c(
  "Treat all"   = "#F8766D",
  "Treat none"  = "#00BA38",
  "模型1_9变量" = "#E76F51",
  "模型2_8变量" = "#00A65A",
  "模型3_7变量" = "#619CFF"
)

# 各策略的线型
strategy_linetypes <- c(
  "Treat all"   = "solid",
  "Treat none"  = "dashed",
  "模型1_9变量" = "solid",
  "模型2_8变量" = "dashed",
  "模型3_7变量" = "dotdash"
)


# 21.5 绘制训练集分类别 DCA ------------------------------------------------

p_dca_by_class_train <- ggplot(
  dca_by_class_train_smooth,
  aes(
    x = threshold,
    y = net_benefit_plot,
    color = 策略,
    linetype = 策略,
    group = 策略
  )
) +
  geom_line(
    linewidth = 1.0,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ 类别,
    nrow = 1,
    labeller = as_labeller(class_labels)
  ) +
  scale_x_continuous(
    limits = c(0.01, 0.80),
    breaks = seq(0.10, 0.80, by = 0.10),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  coord_cartesian(
    ylim = c(-0.10, NA)
  ) +
  scale_color_manual(
    values = strategy_colors,
    breaks = names(strategy_labels),
    labels = unname(strategy_labels)
  ) +
  scale_linetype_manual(
    values = strategy_linetypes,
    breaks = names(strategy_labels),
    labels = unname(strategy_labels)
  ) +
  labs(
    x = "Threshold probability",
    y = "Net benefit",
    color = NULL,
    linetype = NULL
  ) +
  # 所有文字元素统一使用 Times New Roman
  theme_bw(base_size = 16, base_family = font_family) +
  theme(
    # 横纵坐标标题：18 磅
    axis.title.x = element_text(
      size = 18, color = "black",
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      size = 18, color = "black",
      margin = margin(r = 8)
    ),
    
    # 横纵坐标刻度数字：16 磅
    axis.text = element_text(size = 16, color = "black"),
    
    # 顶部分面标题 Low、Intermediate、High：16 磅
    strip.text.x = element_text(
      size = 16, face = "plain", color = "black",
      margin = margin(t = 5, b = 5)
    ),
    strip.background = element_rect(
      fill = "grey85", color = "grey30"
    ),
    
    # 右侧五项图例文字：16 磅
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 16, color = "black"),
    legend.key.width = grid::unit(1.2, "cm"),
    
    # 给放大后的刻度标签留出分面间距
    panel.spacing.x = grid::unit(0.35, "cm"),
    panel.grid.minor = element_blank()
  )

print(p_dca_by_class_train)


ggsave(
  filename = "训练集_9_8_7变量模型_分类别DCA.png",
  plot = p_dca_by_class_train,
  width = 14,
  height = 5.5,
  dpi = 600
)

writexl::write_xlsx(
  list(
    Raw_category_specific_DCA =
      round_df(dca_by_class_train, 4),
    
    Smoothed_category_specific_DCA =
      round_df(dca_by_class_train_smooth, 4)
  ),
  path = "训练集_9_8_7变量模型_分类别DCA结果.xlsx"
)
















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
