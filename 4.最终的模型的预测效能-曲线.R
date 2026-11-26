# 最终7变量有序Logistic回归模型：
#校准曲线（1000次Bootstrap内部验证）和DCA曲线
# 1000次Bootstrap内部验证 + 分类别校准曲线 + 外部验证集分类别校准曲线 + HL检验 ----

library(dplyr)
library(MASS)
library(ggplot2)
library(tidyr)
library(writexl)

# 1. 准备训练集和外部验证集 -------------------------------------------------

train_data <- data %>%
  dplyr::filter(dataset == 1)

valid_data <- data %>%
  dplyr::filter(dataset == 0)

train_data$血药浓度分组 <- ordered(
  train_data$血药浓度分组,
  levels = c(1, 2, 3)
)

valid_data$血药浓度分组 <- ordered(
  valid_data$血药浓度分组,
  levels = c(1, 2, 3)
)

train_data$CKDEPI分组 <- as.numeric(as.character(train_data$CKDEPI分组))
valid_data$CKDEPI分组 <- as.numeric(as.character(valid_data$CKDEPI分组))

class_levels <- levels(train_data$血药浓度分组)


# 2. 最终7变量模型公式 ------------------------------------------------------

final_formula <- 血药浓度分组 ~
  尿素氮 +
  血红蛋白 +
  血小板 +
  年龄 +
  APTT +
  输注时间 +
  CKDEPI分组


# 3. 拟合最终有序Logistic回归模型 ------------------------------------------

fit_final <- MASS::polr(
  final_formula,
  data = train_data,
  Hess = TRUE,
  method = "logistic"
)

summary(fit_final)


# 4. 基础函数 ---------------------------------------------------------------

round_df <- function(df, digits = 4) {
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) round(x, digits) else x
  })
  return(df)
}


format_hl_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ paste0("=", sprintf("%.3f", p))
  )
}


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
        C_index = NA_real_,
        Dxy = NA_real_,
        可比较对子数 = 0
      )
    )
  }
  
  c_index <- (
    sum(prod > 0) +
      0.5 * sum(prod == 0)
  ) / n_pairs
  
  dxy <- 2 * (c_index - 0.5)
  
  data.frame(
    C_index = c_index,
    Dxy = dxy,
    可比较对子数 = n_pairs
  )
}


# 训练集类别构成作为空模型概率
null_probs <- prop.table(table(train_data$血药浓度分组))
null_probs <- as.numeric(null_probs[class_levels])
names(null_probs) <- class_levels


calc_brier_rps <- function(y_true, pred_prob, null_probs) {
  
  classes <- colnames(pred_prob)
  K <- length(classes)
  
  y_char <- as.character(y_true)
  y_fac <- factor(y_char, levels = classes)
  
  y_mat <- model.matrix(~ 0 + y_fac)
  colnames(y_mat) <- classes
  
  # Multiclass Brier score
  brier <- mean(
    rowSums((pred_prob - y_mat)^2)
  )
  
  # Null Brier score
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
  
  # Ranked Probability Score
  pred_cum <- t(
    apply(pred_prob[, 1:(K - 1), drop = FALSE], 1, cumsum)
  )
  
  y_num <- as.numeric(as.character(y_true))
  
  obs_cum <- sapply(1:(K - 1), function(k) {
    as.numeric(y_num <= k)
  })
  
  obs_cum <- as.matrix(obs_cum)
  
  rps <- mean(
    rowSums((pred_cum - obs_cum)^2)
  )
  
  standardized_rps <- rps / (K - 1)
  
  # Null RPS
  null_cum <- cumsum(null_probs)[1:(K - 1)]
  
  null_cum_mat <- matrix(
    rep(null_cum, each = nrow(pred_prob)),
    nrow = nrow(pred_prob),
    byrow = FALSE
  )
  
  null_rps <- mean(
    rowSums((null_cum_mat - obs_cum)^2)
  )
  
  null_standardized_rps <- null_rps / (K - 1)
  
  scaled_rps <- 1 - rps / null_rps
  
  data.frame(
    Multiclass_Brier = brier,
    Null_Brier = null_brier,
    Scaled_Brier = scaled_brier,
    RPS = rps,
    Standardized_RPS = standardized_rps,
    Null_RPS = null_rps,
    Null_Standardized_RPS = null_standardized_rps,
    Scaled_RPS = scaled_rps
  )
}


# 5. 预测效能函数 -----------------------------------------------------------

get_prediction_performance_raw <- function(model, newdata, null_probs) {
  
  pred_prob <- get_pred_prob(model, newdata)
  
  pred_class <- predict(
    model,
    newdata = newdata,
    type = "class"
  )
  
  y_true <- newdata$血药浓度分组
  
  pred_score <- as.numeric(
    pred_prob %*% as.numeric(colnames(pred_prob))
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
  
  data.frame(
    Accuracy = accuracy,
    Misclassification = misclass,
    C_index = cindex_res$C_index,
    Dxy = cindex_res$Dxy,
    Multiclass_Brier = score_res$Multiclass_Brier,
    Scaled_Brier = score_res$Scaled_Brier,
    RPS = score_res$RPS,
    Standardized_RPS = score_res$Standardized_RPS,
    Scaled_RPS = score_res$Scaled_RPS
  )
}


# 6. 校准斜率函数 -----------------------------------------------------------

get_lp_polr <- function(model, newdata) {
  
  mm <- model.matrix(
    delete.response(terms(model)),
    data = newdata
  )
  
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  
  beta <- coef(model)
  mm <- mm[, names(beta), drop = FALSE]
  
  lp <- as.numeric(mm %*% beta)
  
  return(lp)
}


calibration_slope_polr <- function(model, newdata) {
  
  lp <- get_lp_polr(model, newdata)
  
  temp <- data.frame(
    y = newdata$血药浓度分组,
    lp = lp
  )
  
  fit_cal <- tryCatch(
    MASS::polr(
      y ~ lp,
      data = temp,
      Hess = FALSE,
      method = "logistic"
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_cal)) {
    return(NA_real_)
  }
  
  slope <- unname(coef(fit_cal)["lp"])
  
  return(slope)
}


# 7. 训练集表观性能和外部验证性能 ------------------------------------------

apparent_perf <- get_prediction_performance_raw(
  model = fit_final,
  newdata = train_data,
  null_probs = null_probs
)

external_perf <- get_prediction_performance_raw(
  model = fit_final,
  newdata = valid_data,
  null_probs = null_probs
)

apparent_slope <- calibration_slope_polr(
  model = fit_final,
  newdata = train_data
)

external_slope <- calibration_slope_polr(
  model = fit_final,
  newdata = valid_data
)

performance_apparent_external <- dplyr::bind_rows(
  apparent_perf %>%
    dplyr::mutate(
      模型 = "最终模型_7变量",
      数据集 = "训练集表观性能",
      Calibration_slope = apparent_slope
    ),
  external_perf %>%
    dplyr::mutate(
      模型 = "最终模型_7变量",
      数据集 = "外部验证集性能",
      Calibration_slope = external_slope
    )
) %>%
  dplyr::select(
    模型,
    数据集,
    dplyr::everything()
  ) %>%
  round_df(4)

performance_apparent_external


# 8. 平滑校准曲线函数 --------------------------------------------------------

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


# 9. 分类别校准曲线函数 ------------------------------------------------------
# P(Y=1), P(Y=2), P(Y=3)

calibration_smooth_by_class <- function(model, newdata, grid_df, span = 0.8) {
  
  pred_prob <- get_pred_prob(model, newdata)
  y_true <- as.character(newdata$血药浓度分组)
  
  out_list <- list()
  
  for (cls in class_levels) {
    
    grid_cls <- grid_df %>%
      dplyr::filter(类别 == paste0("P(Y=", cls, ")")) %>%
      dplyr::pull(mean_pred)
    
    prob_cls <- pred_prob[, cls]
    event_cls <- as.numeric(y_true == cls)
    
    observed_smooth <- smooth_curve(
      prob = prob_cls,
      event = event_cls,
      grid = grid_cls,
      span = span
    )
    
    out_list[[cls]] <- data.frame(
      类别 = paste0("P(Y=", cls, ")"),
      mean_pred = grid_cls,
      observed = observed_smooth
    )
  }
  
  dplyr::bind_rows(out_list)
}


# 10. 建立训练集分类别校准曲线x轴网格 --------------------------------------

pred_prob_final_train <- get_pred_prob(
  model = fit_final,
  newdata = train_data
)

grid_df <- lapply(class_levels, function(cls) {
  
  prob_cls <- pred_prob_final_train[, cls]
  
  grid_min <- max(0.01, min(prob_cls, na.rm = TRUE))
  grid_max <- min(0.99, max(prob_cls, na.rm = TRUE))
  
  data.frame(
    类别 = paste0("P(Y=", cls, ")"),
    mean_pred = seq(grid_min, grid_max, length.out = 100)
  )
}) %>%
  dplyr::bind_rows()


# 11. 进行1000次Bootstrap内部验证 ------------------------------------------

set.seed(123)

B <- 1000

boot_results_list <- vector("list", B)
curve_optimism_list <- vector("list", B)

for (b in 1:B) {
  
  boot_id <- sample(
    1:nrow(train_data),
    size = nrow(train_data),
    replace = TRUE
  )
  
  boot_data <- train_data[boot_id, ]
  
  boot_data$血药浓度分组 <- ordered(
    boot_data$血药浓度分组,
    levels = c(1, 2, 3)
  )
  
  # 如果Bootstrap样本偶然缺少某个结局类别，则跳过
  if (length(unique(boot_data$血药浓度分组)) < 3) {
    next
  }
  
  fit_boot <- tryCatch(
    suppressWarnings(
      MASS::polr(
        final_formula,
        data = boot_data,
        Hess = FALSE,
        method = "logistic"
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_boot)) {
    next
  }
  
  # Bootstrap样本中的表观性能
  perf_boot_app <- tryCatch(
    get_prediction_performance_raw(
      model = fit_boot,
      newdata = boot_data,
      null_probs = null_probs
    ),
    error = function(e) NULL
  )
  
  # Bootstrap模型回代到原训练集中的性能
  perf_boot_test <- tryCatch(
    get_prediction_performance_raw(
      model = fit_boot,
      newdata = train_data,
      null_probs = null_probs
    ),
    error = function(e) NULL
  )
  
  if (is.null(perf_boot_app) || is.null(perf_boot_test)) {
    next
  }
  
  slope_boot_app <- calibration_slope_polr(
    model = fit_boot,
    newdata = boot_data
  )
  
  slope_boot_test <- calibration_slope_polr(
    model = fit_boot,
    newdata = train_data
  )
  
  boot_results_list[[b]] <- data.frame(
    bootstrap = b,
    
    Accuracy_app = perf_boot_app$Accuracy,
    Accuracy_test = perf_boot_test$Accuracy,
    
    Misclassification_app = perf_boot_app$Misclassification,
    Misclassification_test = perf_boot_test$Misclassification,
    
    C_index_app = perf_boot_app$C_index,
    C_index_test = perf_boot_test$C_index,
    
    Dxy_app = perf_boot_app$Dxy,
    Dxy_test = perf_boot_test$Dxy,
    
    Multiclass_Brier_app = perf_boot_app$Multiclass_Brier,
    Multiclass_Brier_test = perf_boot_test$Multiclass_Brier,
    
    Scaled_Brier_app = perf_boot_app$Scaled_Brier,
    Scaled_Brier_test = perf_boot_test$Scaled_Brier,
    
    RPS_app = perf_boot_app$RPS,
    RPS_test = perf_boot_test$RPS,
    
    Standardized_RPS_app = perf_boot_app$Standardized_RPS,
    Standardized_RPS_test = perf_boot_test$Standardized_RPS,
    
    Scaled_RPS_app = perf_boot_app$Scaled_RPS,
    Scaled_RPS_test = perf_boot_test$Scaled_RPS,
    
    Calibration_slope_app = slope_boot_app,
    Calibration_slope_test = slope_boot_test
  )
  
  
  # Bootstrap分类别校准曲线乐观偏倚
  curve_app <- tryCatch(
    calibration_smooth_by_class(
      model = fit_boot,
      newdata = boot_data,
      grid_df = grid_df,
      span = 0.8
    ) %>%
      dplyr::rename(observed_app = observed),
    error = function(e) NULL
  )
  
  curve_test <- tryCatch(
    calibration_smooth_by_class(
      model = fit_boot,
      newdata = train_data,
      grid_df = grid_df,
      span = 0.8
    ) %>%
      dplyr::rename(observed_test = observed),
    error = function(e) NULL
  )
  
  if (!is.null(curve_app) && !is.null(curve_test)) {
    
    curve_optimism_list[[b]] <- curve_app %>%
      dplyr::left_join(
        curve_test,
        by = c("类别", "mean_pred")
      ) %>%
      dplyr::mutate(
        bootstrap = b,
        optimism = observed_app - observed_test
      )
  }
  
  if (b %% 50 == 0) {
    message("Bootstrap完成次数：", b)
  }
}

boot_results <- dplyr::bind_rows(boot_results_list)

curve_optimism <- dplyr::bind_rows(curve_optimism_list)

boot_results


# 12. 计算Bootstrap乐观校正性能 --------------------------------------------

higher_better <- c(
  "Accuracy",
  "C_index",
  "Dxy",
  "Scaled_Brier",
  "Scaled_RPS"
)

lower_better <- c(
  "Misclassification",
  "Multiclass_Brier",
  "RPS",
  "Standardized_RPS"
)


correct_high <- lapply(higher_better, function(metric) {
  
  app_col <- paste0(metric, "_app")
  test_col <- paste0(metric, "_test")
  
  optimism <- mean(
    boot_results[[app_col]] - boot_results[[test_col]],
    na.rm = TRUE
  )
  
  apparent <- apparent_perf[[metric]]
  
  corrected <- apparent - optimism
  
  data.frame(
    指标 = metric,
    指标方向 = "越大越好",
    训练集表观性能 = apparent,
    平均乐观偏倚 = optimism,
    Bootstrap校正性能 = corrected
  )
})


correct_low <- lapply(lower_better, function(metric) {
  
  app_col <- paste0(metric, "_app")
  test_col <- paste0(metric, "_test")
  
  optimism <- mean(
    boot_results[[test_col]] - boot_results[[app_col]],
    na.rm = TRUE
  )
  
  apparent <- apparent_perf[[metric]]
  
  corrected <- apparent + optimism
  
  data.frame(
    指标 = metric,
    指标方向 = "越小越好",
    训练集表观性能 = apparent,
    平均乐观偏倚 = optimism,
    Bootstrap校正性能 = corrected
  )
})


# 校准斜率
slope_optimism <- mean(
  boot_results$Calibration_slope_app - boot_results$Calibration_slope_test,
  na.rm = TRUE
)

slope_corrected <- apparent_slope - slope_optimism

correct_slope <- data.frame(
  指标 = "Calibration_slope",
  指标方向 = "理想值为1",
  训练集表观性能 = apparent_slope,
  平均乐观偏倚 = slope_optimism,
  Bootstrap校正性能 = slope_corrected
)


internal_validation_results <- dplyr::bind_rows(
  correct_high,
  correct_low,
  correct_slope
) %>%
  dplyr::mutate(
    训练集表观性能 = round(训练集表观性能, 4),
    平均乐观偏倚 = round(平均乐观偏倚, 4),
    Bootstrap校正性能 = round(Bootstrap校正性能, 4)
  )

internal_validation_results


# 13. Bootstrap校正后的训练集分类别校准曲线数据 -----------------------------

apparent_cal_class_curve <- calibration_smooth_by_class(
  model = fit_final,
  newdata = train_data,
  grid_df = grid_df,
  span = 0.8
) %>%
  dplyr::rename(Apparent = observed)


mean_curve_optimism <- curve_optimism %>%
  dplyr::group_by(类别, mean_pred) %>%
  dplyr::summarise(
    平均乐观偏倚 = mean(optimism, na.rm = TRUE),
    .groups = "drop"
  )


bootstrap_corrected_class_curve <- apparent_cal_class_curve %>%
  dplyr::left_join(
    mean_curve_optimism,
    by = c("类别", "mean_pred")
  ) %>%
  dplyr::mutate(
    平均乐观偏倚 = ifelse(is.na(平均乐观偏倚), 0, 平均乐观偏倚),
    `Bias-corrected` = Apparent - 平均乐观偏倚,
    `Bias-corrected` = pmin(pmax(`Bias-corrected`, 0), 1)
  )

bootstrap_corrected_class_curve


# 14. 分类别Hosmer-Lemeshow检验函数 ----------------------------------------
# 说明：
# 有序三分类结局不能直接做一个传统二分类HL检验。
# 这里分别对P(Y=1)、P(Y=2)、P(Y=3)做one-vs-rest形式的HL检验。

hosmer_lemeshow_binary <- function(event, prob, g = 10) {
  
  df <- data.frame(
    event = as.numeric(event),
    prob = as.numeric(prob)
  ) %>%
    dplyr::filter(!is.na(event), !is.na(prob))
  
  # 基本检查
  if (
    nrow(df) < 10 ||
    length(unique(df$event)) < 2 ||
    length(unique(df$prob)) < 3
  ) {
    return(
      data.frame(
        HL_chisq = NA_real_,
        df = NA_integer_,
        P值 = NA_real_,
        分组数 = NA_integer_
      )
    )
  }
  
  g <- min(g, nrow(df))
  
  df <- df %>%
    dplyr::arrange(prob) %>%
    dplyr::mutate(
      group = dplyr::ntile(dplyr::row_number(), g)
    )
  
  hl_tab <- df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      observed = sum(event),
      expected = sum(prob),
      .groups = "drop"
    )
  
  df_hl <- nrow(hl_tab) - 2
  
  if (df_hl <= 0) {
    return(
      data.frame(
        HL_chisq = NA_real_,
        df = df_hl,
        P值 = NA_real_,
        分组数 = nrow(hl_tab)
      )
    )
  }
  
  expected_event <- pmax(hl_tab$expected, 1e-10)
  expected_nonevent <- pmax(hl_tab$n - hl_tab$expected, 1e-10)
  
  hl_stat <- sum(
    (hl_tab$observed - hl_tab$expected)^2 / expected_event +
      ((hl_tab$n - hl_tab$observed) - (hl_tab$n - hl_tab$expected))^2 /
      expected_nonevent
  )
  
  p_value <- pchisq(
    hl_stat,
    df = df_hl,
    lower.tail = FALSE
  )
  
  data.frame(
    HL_chisq = hl_stat,
    df = df_hl,
    P值 = p_value,
    分组数 = nrow(hl_tab)
  )
}


hosmer_lemeshow_by_class <- function(model, newdata, dataset_name, g = 10) {
  
  pred_prob <- get_pred_prob(
    model = model,
    newdata = newdata
  )
  
  y_true <- as.character(newdata$血药浓度分组)
  
  out <- lapply(class_levels, function(cls) {
    
    event_cls <- as.numeric(y_true == cls)
    prob_cls <- pred_prob[, cls]
    
    hl_res <- hosmer_lemeshow_binary(
      event = event_cls,
      prob = prob_cls,
      g = g
    )
    
    data.frame(
      数据集 = dataset_name,
      类别 = paste0("P(Y=", cls, ")"),
      HL_chisq = hl_res$HL_chisq,
      df = hl_res$df,
      P值 = hl_res$P值,
      分组数 = hl_res$分组数
    )
  })
  
  dplyr::bind_rows(out)
}


# 训练集分类别HL检验，训练集用9组
hl_train_by_class <- hosmer_lemeshow_by_class(
  model = fit_final,
  newdata = train_data,
  dataset_name = "训练集",
  g = 9
) %>%
  dplyr::mutate(
    HL_chisq = round(HL_chisq, 4),
    P值_格式 = format_hl_p(P值)
  )

hl_train_by_class


# 外部验证集分类别HL检验，外部验证集样本量较小，建议用5组或6组
hl_valid_by_class <- hosmer_lemeshow_by_class(
  model = fit_final,
  newdata = valid_data,
  dataset_name = "外部验证集",
  g = 6
) %>%
  dplyr::mutate(
    HL_chisq = round(HL_chisq, 4),
    P值_格式 = format_hl_p(P值)
  )

hl_valid_by_class




# 15. 绘制训练集 Bootstrap 校正后的分类别校准曲线 --------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# 如尚未安装字体包，首次运行时取消下一行的注释
# install.packages(c("sysfonts", "showtext"))

# 15.1 注册 Windows 系统中的 Times New Roman 字体 --------------------------

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

# 15.2 设置分面顺序和英文名称 ----------------------------------------------

class_labels <- c(
  "P(Y=1)" = "Low",
  "P(Y=2)" = "Intermediate",
  "P(Y=3)" = "High"
)

curve_order <- c("Ideal", "Apparent", "Bias-corrected")

# 15.3 整理训练集分类别校准曲线数据 ----------------------------------------

# 提取 Apparent 和 Bias-corrected 曲线
plot_boot_class_curve <- bootstrap_corrected_class_curve %>%
  dplyr::select(类别, mean_pred, Apparent, `Bias-corrected`) %>%
  tidyr::pivot_longer(
    cols = c(Apparent, `Bias-corrected`),
    names_to = "曲线类型",
    values_to = "observed"
  ) %>%
  dplyr::mutate(
    类别 = factor(类别, levels = names(class_labels))
  )

# 为每个分面生成 Ideal 对角线
ideal_curve_train <- tidyr::expand_grid(
  类别 = factor(names(class_labels), levels = names(class_labels)),
  mean_pred = c(0, 1)
) %>%
  dplyr::mutate(
    observed = mean_pred,
    曲线类型 = "Ideal"
  )

# 合并三种曲线，并固定图例顺序
plot_boot_class_curve <- dplyr::bind_rows(
  ideal_curve_train,
  plot_boot_class_curve
) %>%
  dplyr::mutate(
    曲线类型 = factor(曲线类型, levels = curve_order)
  )

# 15.4 准备各面板的 Hosmer–Lemeshow P 值标注 -------------------------------
# 曲线和 P 值使用相同的类别顺序，确保标注出现在对应面板

hl_train_label <- hl_train_by_class %>%
  mutate(
    类别 = factor(类别, levels = names(class_labels)),
    label = paste0("Hosmer-Lemeshow P", P值_格式)
  )

# 15.5 绘制训练集分类别校准曲线 --------------------------------------------

p_boot_class_calibration <- ggplot(
  plot_boot_class_curve,
  aes(
    x = mean_pred,
    y = observed,
    color = 曲线类型,
    linetype = 曲线类型
  )
) +
  geom_line(linewidth = 1.1) +
  
  # 面板内的文字需要单独指定字体和字号
  geom_text(
    data = hl_train_label,
    aes(label = label),
    x = 0.04,
    y = 0.95,
    hjust = 0,
    family = font_family,
    size = 14 / ggplot2::.pt,  # 将 14 磅换算为 geom_text 使用的毫米
    color = "black",
    inherit.aes = FALSE
  ) +
  
  facet_wrap(
    ~ 类别,
    nrow = 1,
    labeller = as_labeller(class_labels)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_color_manual(
    breaks = curve_order,
    values = c(
      "Ideal"          = "grey",
      "Apparent"       = "#E69F00",
      "Bias-corrected" = "#1E88E5"
    )
  ) +
  scale_linetype_manual(
    breaks = curve_order,
    values = c(
      "Ideal"          = "dashed",
      "Apparent"       = "solid",
      "Bias-corrected" = "dotdash"
    )
  ) +
  labs(
    x = "Predicted probability",
    y = "Observed proportion",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 16, base_family = font_family) +
  theme(
    # 横纵坐标标题：18 磅
    axis.title.x = element_text(
      size = 18, color = "black", margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      size = 18, color = "black", margin = margin(r = 8)
    ),
    
    # 横纵坐标刻度数字：16 磅
    axis.text = element_text(size = 16, color = "black"),
    
    # 顶部分面标题 Low、Intermediate、High：16 磅
    strip.text.x = element_text(
      size = 16, face = "plain", color = "black",
      margin = margin(t = 5, b = 5)
    ),
    strip.background = element_rect(fill = "grey85", color = "grey30"),
    
    # 右侧图例：16 磅
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 16, color = "black"),
    legend.key.width = grid::unit(1.2, "cm"),
    
    panel.spacing.x = grid::unit(0.35, "cm"),
    panel.grid.minor = element_blank()
  )

print(p_boot_class_calibration)





# 16. 绘制外部验证集分类别校准曲线，并显示HL检验 ---------------------------
# 外部验证集图例也显示：Ideal、Apparent、Bias-corrected三根线
# 注意：
# 外部验证集严格来说不进行Bootstrap内部乐观偏倚校正。
# 这里的Bias-corrected是将训练集Bootstrap得到的平均乐观偏倚插值到外部验证集预测概率网格后得到的校正参考线。
# 论文解释时建议说明：外部验证集主要依据Apparent曲线和HL检验评价模型校准度。
# 外部验证集分类别校准曲线：数据整理与绘图 

library(dplyr)
library(tidyr)
library(ggplot2)

# 16.1 设置 Times New Roman 字体 -------------------------------------------

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

# 16.2 设置分面名称、分面顺序和图例顺序 ------------------------------------

class_labels <- c(
  "P(Y=1)" = "Low",
  "P(Y=2)" = "Intermediate",
  "P(Y=3)" = "High"
)

curve_order <- c("Ideal", "Apparent", "Bias-corrected")

# 16.3 整理外部验证集曲线数据 ----------------------------------------------
# Bias-corrected 为使用训练集 Bootstrap 平均乐观偏倚得到的参考曲线

plot_external_class_curve <- external_corrected_class_curve %>%
  dplyr::select(类别, mean_pred, Apparent, `Bias-corrected`) %>%
  tidyr::pivot_longer(
    cols = c(Apparent, `Bias-corrected`),
    names_to = "曲线类型",
    values_to = "observed"
  )

# 每个分面添加一条从 (0, 0) 到 (1, 1) 的 Ideal 对角线
ideal_curve_valid <- tidyr::expand_grid(
  类别 = names(class_labels),
  mean_pred = c(0, 1)
) %>%
  dplyr::mutate(
    observed = mean_pred,
    曲线类型 = "Ideal"
  )

plot_external_class_curve <- dplyr::bind_rows(
  ideal_curve_valid,
  plot_external_class_curve
) %>%
  dplyr::mutate(
    类别 = factor(类别, levels = names(class_labels)),
    曲线类型 = factor(曲线类型, levels = curve_order)
  )

# 16.4 准备外部验证集各面板的 HL 检验标注 ----------------------------------

hl_valid_label <- hl_valid_by_class %>%
  dplyr::mutate(
    类别 = factor(类别, levels = names(class_labels)),
    label = paste0("Hosmer-Lemeshow P", P值_格式)
  )

# 16.5 绘制外部验证集分类别校准曲线 ----------------------------------------

p_external_class_calibration <- ggplot(
  plot_external_class_curve,
  aes(
    x = mean_pred,
    y = observed,
    color = 曲线类型,
    linetype = 曲线类型
  )
) +
  geom_line(linewidth = 1.1) +
  
  # 各面板的 HL 检验标注：15 磅
  geom_text(
    data = hl_valid_label,
    aes(label = label),
    x = 0.01,
    y = 0.95,
    hjust = 0,
    family = font_family,
    size = 14 / ggplot2::.pt,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  facet_wrap(
    ~ 类别,
    nrow = 1,
    labeller = as_labeller(class_labels)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_color_manual(
    breaks = curve_order,
    values = c(
      "Ideal"          = "grey",
      "Apparent"       = "#E69F00",
      "Bias-corrected" = "#1E88E5"
    )
  ) +
  scale_linetype_manual(
    breaks = curve_order,
    values = c(
      "Ideal"          = "dashed",
      "Apparent"       = "solid",
      "Bias-corrected" = "dotdash"
    )
  ) +
  labs(
    x = "Predicted probability",
    y = "Observed proportion",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 16, base_family = font_family) +
  theme(
    # 横纵坐标标题：18 磅
    axis.title.x = element_text(
      size = 18, color = "black", margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      size = 18, color = "black", margin = margin(r = 8)
    ),
    
    # 横纵坐标刻度数字：16 磅
    axis.text = element_text(size = 16, color = "black"),
    
    # 顶部分面标题：16 磅
    strip.text.x = element_text(
      size = 16, face = "plain", color = "black",
      margin = margin(t = 5, b = 5)
    ),
    strip.background = element_rect(fill = "grey85", color = "grey30"),
    
    # 右侧图例：16 磅
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 16, color = "black"),
    legend.key.width = grid::unit(1.2, "cm"),
    
    panel.spacing.x = grid::unit(0.35, "cm"),
    panel.grid.minor = element_blank()
  )

print(p_external_class_calibration)








# 18. 分类别DCA曲线 ---------------------------------------------------------
# 方式一：两张图，每张3个子图
# 图1：训练集分类别DCA，分别为 P(Y=1)、P(Y=2)、P(Y=3)
# 图2：外部验证集分类别DCA，分别为 P(Y=1)、P(Y=2)、P(Y=3)
#
# 每个子图均采用 one-vs-rest 方式：
# P(Y=1): 事件 = Y=1，非事件 = Y=2或Y=3，预测概率 = P(Y=1)
# P(Y=2): 事件 = Y=2，非事件 = Y=1或Y=3，预测概率 = P(Y=2)
# P(Y=3): 事件 = Y=3，非事件 = Y=1或Y=2，预测概率 = P(Y=3)


# 18.1 分类别DCA计算函数 ---------------------------------------------------

calculate_dca_by_class <- function(model, newdata, dataset_name,
                                   thresholds = seq(0.01, 0.80, by = 0.01)) {
  
  pred_prob <- get_pred_prob(
    model = model,
    newdata = newdata
  )
  
  y_true <- as.character(newdata$血药浓度分组)
  
  out_list <- list()
  
  for (cls in class_levels) {
    
    # 当前类别的预测概率
    prob_cls <- pred_prob[, cls]
    
    # one-vs-rest事件定义：
    # 当前类别为事件，其余类别为非事件
    event_cls <- ifelse(y_true == cls, 1, 0)
    
    n <- length(event_cls)
    prevalence <- mean(event_cls)
    
    dca_cls <- lapply(thresholds, function(pt) {
      
      pred_positive <- ifelse(prob_cls >= pt, 1, 0)
      
      TP <- sum(pred_positive == 1 & event_cls == 1)
      FP <- sum(pred_positive == 1 & event_cls == 0)
      FN <- sum(pred_positive == 0 & event_cls == 1)
      TN <- sum(pred_positive == 0 & event_cls == 0)
      
      # Model策略净获益
      nb_model <- TP / n - FP / n * (pt / (1 - pt))
      
      # All策略净获益
      nb_all <- prevalence - (1 - prevalence) * (pt / (1 - pt))
      
      # None策略净获益
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
        Model = nb_model,
        All = nb_all,
        None = nb_none
      )
    }) %>%
      dplyr::bind_rows()
    
    out_list[[cls]] <- dca_cls
  }
  
  out_wide <- dplyr::bind_rows(out_list)
  
  out_long <- out_wide %>%
    dplyr::select(
      数据集,
      类别,
      threshold,
      prevalence,
      Model,
      All,
      None
    ) %>%
    tidyr::pivot_longer(
      cols = c(Model, All, None),
      names_to = "策略",
      values_to = "net_benefit"
    ) %>%
    dplyr::mutate(
      策略 = factor(
        策略,
        levels = c("All", "None", "Model")
      ),
      类别 = factor(
        类别,
        levels = paste0("P(Y=", class_levels, ")")
      ),
      数据集 = factor(
        数据集,
        levels = c("训练集", "外部验证集")
      )
    )
  
  return(out_long)
}


# 18.2 对Model曲线进行loess平滑 --------------------------------------------
# 说明：
# DCA中Model曲线因为TP、FP随阈值概率跳跃变化，尤其在外部验证集样本量较小时会呈锯齿状。
# 这里仅对Model曲线做loess平滑，All和None保留原始曲线。
# 平滑只用于图形展示，不改变原始DCA计算结果。

smooth_dca_model_curve <- function(dca_data, span = 0.25) {
  
  out <- dca_data %>%
    dplyr::group_by(数据集, 类别, 策略) %>%
    dplyr::group_modify(function(.x, .y) {
      
      df <- .x %>%
        dplyr::arrange(threshold)
      
      if (as.character(.y$策略) != "Model") {
        df$net_benefit_plot <- df$net_benefit
        return(df)
      }
      
      # 如果点数太少或净获益变化太少，则不平滑
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
          control = loess.control(surface = "direct")
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
          newdata = data.frame(threshold = df$threshold)
        )
      )
      
      # 若平滑后出现NA，则用线性插值补齐
      if (any(is.na(pred))) {
        
        not_na <- !is.na(pred)
        
        if (sum(not_na) >= 2) {
          pred <- approx(
            x = df$threshold[not_na],
            y = pred[not_na],
            xout = df$threshold,
            rule = 2
          )$y
        } else {
          pred[is.na(pred)] <- df$net_benefit[is.na(pred)]
        }
      }
      
      df$net_benefit_plot <- pred
      
      return(df)
    }) %>%
    dplyr::ungroup()
  
  return(out)
}


# 18.3 计算训练集和外部验证集分类别DCA数据 ---------------------------------

dca_by_class_train <- calculate_dca_by_class(
  model = fit_final,
  newdata = train_data,
  dataset_name = "训练集",
  thresholds = seq(0.01, 0.80, by = 0.01)
)

dca_by_class_valid <- calculate_dca_by_class(
  model = fit_final,
  newdata = valid_data,
  dataset_name = "外部验证集",
  thresholds = seq(0.01, 0.80, by = 0.01)
)

dca_by_class_results <- dplyr::bind_rows(
  dca_by_class_train,
  dca_by_class_valid
)

# 对Model曲线进行平滑
# 训练集样本量较大，可用span=0.25
# 外部验证集样本量较小，可用span=0.35，使曲线更平滑

dca_by_class_train_smooth <- smooth_dca_model_curve(
  dca_data = dca_by_class_train,
  span = 0.25
)

dca_by_class_valid_smooth <- smooth_dca_model_curve(
  dca_data = dca_by_class_valid,
  span = 0.35
)

dca_by_class_results_smooth <- dplyr::bind_rows(
  dca_by_class_train_smooth,
  dca_by_class_valid_smooth
)

dca_by_class_results_smooth


# 18.4 设置绘图格式并绘制训练集分类别 DCA 曲线 ----------------------------

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

# 分面名称与显示顺序
class_labels <- c(
  "P(Y=1)" = "Low",
  "P(Y=2)" = "Intermediate",
  "P(Y=3)" = "High"
)

# 训练集和外部验证集共用同一个绘图函数
plot_dca_by_class <- function(dca_data) {
  
  dca_data$类别 <- factor(dca_data$类别, levels = names(class_labels))
  
  ggplot(
    dca_data,
    aes(
      x = threshold,
      y = net_benefit_plot,
      color = 策略,
      linetype = 策略
    )
  ) +
    geom_line(linewidth = 1.0) +
    facet_wrap(
      ~ 类别,
      nrow = 1,
      labeller = as_labeller(class_labels)
    ) +
    scale_x_continuous(
      limits = c(0.01, 0.80),
      breaks = seq(0, 0.80, by = 0.10),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    coord_cartesian(
      ylim = c(-0.10, 0.35)
    ) +
    scale_color_manual(
      breaks = c("All", "None", "Model"),
      values = c(
        "All"   = "#F8766D",
        "None"  = "#00BA38",
        "Model" = "#619CFF"
      )
    ) +
    scale_linetype_manual(
      breaks = c("All", "None", "Model"),
      values = c(
        "All"   = "solid",
        "None"  = "dashed",
        "Model" = "dashed"
      )
    ) +
    labs(
      x = "Threshold probability",
      y = "Net benefit",
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 16, base_family = font_family) +
    theme(
      # 横纵坐标标题：18 磅
      axis.title.x = element_text(
        size = 18, color = "black", margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 18, color = "black", margin = margin(r = 8)
      ),
      
      # 横纵坐标刻度数字：16 磅
      axis.text = element_text(size = 16, color = "black"),
      
      # 分面标题 Low、Intermediate、High：16 磅
      strip.text.x = element_text(
        size = 16, face = "plain", color = "black",
        margin = margin(t = 5, b = 5)
      ),
      strip.background = element_rect(fill = "grey85", color = "grey30"),
      
      # 右侧图例 All、None、Model：16 磅
      legend.position = "right",
      legend.text = element_text(size = 16, color = "black"),
      legend.key.width = grid::unit(1.1, "cm"),
      
      panel.spacing.x = grid::unit(0.35, "cm"),
      panel.grid.minor = element_blank()
    )
}

# 训练集
p_dca_by_class_train <- plot_dca_by_class(dca_by_class_train_smooth)
print(p_dca_by_class_train)


# 18.5 绘制外部验证集分类别 DCA 曲线 --------------------------------------

p_dca_by_class_valid <- plot_dca_by_class(dca_by_class_valid_smooth)
print(p_dca_by_class_valid)




# 19. 导出结果 --------------------------------------------------------------

write_xlsx(
  list(
    训练集表观与外部验证性能 = performance_apparent_external,
    Bootstrap内部验证校正结果 = internal_validation_results,
    Bootstrap每次重抽样结果 = boot_results,
    训练集Bootstrap校正分类别校准曲线数据 = bootstrap_corrected_class_curve,
    外部验证分类别校准曲线数据 = external_corrected_class_curve,
    训练集分类别HL检验 = hl_train_by_class,
    外部验证集分类别HL检验 = hl_valid_by_class
  ),
  path = "最终7变量有序Logistic模型_1000次Bootstrap内部验证_分类别校准曲线_HL检验.xlsx"
)


# 20. 保存图片 --------------------------------------------------------------

ggsave(
  "最终7变量模型_训练集Bootstrap校正分类别校准曲线_HL检验.png",
  p_boot_class_calibration,
  width = 11,
  height = 4,
  dpi = 300
)

ggsave(
  "最终7变量模型_外部验证集分类别校准曲线_HL检验.png",
  p_external_class_calibration,
  width = 11,
  height = 4,
  dpi = 300
)


# 21. 查看主要结果 ----------------------------------------------------------

performance_apparent_external
internal_validation_results
hl_train_by_class
hl_valid_by_class

p_boot_class_calibration
p_external_class_calibration









# 12A. C-index及95%置信区间 -------------------------------------------------
# 说明：
# 1）训练集和外部验证集C-index的95%CI：
#    固定最终模型fit_final，对相应数据集进行个体层面Bootstrap重抽样，
#    每次计算一次C-index，取2.5%和97.5%分位数作为95%CI。
#
# 2）Bootstrap内部验证校正C-index的95%CI：
#    每次Bootstrap中：
#    optimism_b = C_index_app_b - C_index_test_b
#    corrected_b = apparent_C_index - optimism_b
#    对corrected_b取2.5%和97.5%分位数作为校正C-index的95%CI。

# 12A.1 提取每个个体的有序预测分数 ----------------------------------------

get_pred_score_polr <- function(model, newdata) {
  
  pred_prob <- get_pred_prob(
    model = model,
    newdata = newdata
  )
  
  pred_score <- as.numeric(
    pred_prob %*% as.numeric(colnames(pred_prob))
  )
  
  return(pred_score)
}


# 12A.2 固定模型，Bootstrap计算C-index的95%CI -------------------------------

bootstrap_cindex_ci_fixed_model <- function(model, newdata, dataset_name,
                                            B = 1000, seed = 123) {
  
  set.seed(seed)
  
  y <- as.numeric(as.character(newdata$血药浓度分组))
  score <- get_pred_score_polr(
    model = model,
    newdata = newdata
  )
  
  # 原始数据C-index
  cindex_original <- ordinal_cindex(
    y = y,
    score = score
  )$C_index
  
  boot_cindex <- rep(NA_real_, B)
  
  for (b in 1:B) {
    
    id <- sample(
      seq_len(nrow(newdata)),
      size = nrow(newdata),
      replace = TRUE
    )
    
    y_b <- y[id]
    score_b <- score[id]
    
    # 如果Bootstrap样本中结局类别不足，可能没有可比较对子，返回NA
    boot_cindex[b] <- ordinal_cindex(
      y = y_b,
      score = score_b
    )$C_index
  }
  
  boot_cindex <- boot_cindex[!is.na(boot_cindex)]
  
  out <- data.frame(
    数据集 = dataset_name,
    C_index = cindex_original,
    lower_95 = as.numeric(quantile(boot_cindex, 0.025, na.rm = TRUE)),
    upper_95 = as.numeric(quantile(boot_cindex, 0.975, na.rm = TRUE)),
    Bootstrap有效次数 = length(boot_cindex)
  )
  
  return(out)
}


# 12A.3 训练集表观C-index及95%CI -------------------------------------------

cindex_ci_train <- bootstrap_cindex_ci_fixed_model(
  model = fit_final,
  newdata = train_data,
  dataset_name = "训练集表观性能",
  B = 1000,
  seed = 123
)


# 12A.4 外部验证集C-index及95%CI -------------------------------------------

cindex_ci_valid <- bootstrap_cindex_ci_fixed_model(
  model = fit_final,
  newdata = valid_data,
  dataset_name = "外部验证集性能",
  B = 1000,
  seed = 456
)


# 12A.5 Bootstrap内部验证校正C-index及95%CI --------------------------------
# 这里使用你前面Bootstrap循环中已经得到的：
# C_index_app  = Bootstrap模型在Bootstrap样本中的C-index
# C_index_test = Bootstrap模型回代原始训练集中的C-index

apparent_cindex <- apparent_perf$C_index

boot_cindex_corrected_df <- boot_results %>%
  dplyr::filter(
    !is.na(C_index_app),
    !is.na(C_index_test)
  ) %>%
  dplyr::mutate(
    C_index_optimism = C_index_app - C_index_test,
    C_index_corrected_each = apparent_cindex - C_index_optimism
  )

cindex_ci_bootstrap_corrected <- boot_cindex_corrected_df %>%
  dplyr::summarise(
    数据集 = "Bootstrap内部验证校正性能",
    C_index = apparent_cindex - mean(C_index_optimism, na.rm = TRUE),
    lower_95 = as.numeric(quantile(C_index_corrected_each, 0.025, na.rm = TRUE)),
    upper_95 = as.numeric(quantile(C_index_corrected_each, 0.975, na.rm = TRUE)),
    Bootstrap有效次数 = dplyr::n()
  )


# 12A.6 合并C-index结果 -----------------------------------------------------

cindex_summary <- dplyr::bind_rows(
  cindex_ci_train,
  cindex_ci_valid,
  cindex_ci_bootstrap_corrected
) %>%
  dplyr::mutate(
    C_index = round(C_index, 4),
    lower_95 = round(lower_95, 4),
    upper_95 = round(upper_95, 4),
    C_index_95CI = paste0(
      sprintf("%.4f", C_index),
      " (",
      sprintf("%.4f", lower_95),
      "-",
      sprintf("%.4f", upper_95),
      ")"
    )
  )

cindex_summary





# 补充时间外部验证队列的三个指标（总体校准斜率和三分类 Brier 分数，校准截距） ---------------------------------------------------------
# 时间验证队列的分类别校准截距和斜率
prob_valid <- get_pred_prob(fit_final, valid_data)
y_valid <- as.character(valid_data$血药浓度分组)

stopifnot(
  nrow(prob_valid) == length(y_valid),
  all(y_valid %in% colnames(prob_valid)),
  all(is.finite(prob_valid))
)

category_labels <- c(
  "1" = "Subtherapeutic",
  "2" = "Therapeutic",
  "3" = "High exposure"
)

eps <- 1e-8

calibration_valid <- do.call(rbind, lapply(colnames(prob_valid), function(cls) {
  cal_data <- data.frame(
    event = as.integer(y_valid == cls),
    logit_p = qlogis(pmin(pmax(prob_valid[, cls], eps), 1 - eps))
  )
  
  intercept_fit <- glm(
    event ~ 1 + offset(logit_p),
    family = binomial(),
    data = cal_data
  )
  
  slope_fit <- glm(
    event ~ logit_p,
    family = binomial(),
    data = cal_data
  )
  
  data.frame(
    Category = unname(category_labels[cls]),
    Events = sum(cal_data$event),
    Calibration_intercept = unname(coef(intercept_fit)["(Intercept)"]),
    Calibration_slope = unname(coef(slope_fit)["logit_p"])
  )
}))

print(round_df(calibration_valid, 3), row.names = FALSE)

# 原脚本已计算的时间验证队列总体指标
print(round(data.frame(
  Common_ordinal_calibration_slope = external_slope,
  Multiclass_Brier = external_perf$Multiclass_Brier
), 3), row.names = FALSE)


# 每个类别的 one-vs-rest Brier 分数，均使用全部 71 例
brier_by_class <- vapply(
  colnames(prob_valid),
  function(cls) {
    observed <- as.integer(y_valid == cls)
    mean((prob_valid[, cls] - observed)^2)
  },
  FUN.VALUE = numeric(1)
)

# 核对分类别 Brier 分数之和是否等于原脚本的三分类 Brier 分数
brier_overall <- external_perf$Multiclass_Brier
stopifnot(isTRUE(all.equal(sum(brier_by_class), brier_overall)))

# 按类别名称对齐校准结果和 Brier 分数
class_codes <- names(category_labels)[
  match(calibration_valid$Category, unname(category_labels))
]
stopifnot(all(class_codes %in% names(brier_by_class)))

temporal_validation_table <- rbind(
  data.frame(
    Scope = "Overall temporal validation cohort",
    Patients_or_events = length(y_valid),
    Calibration_intercept = "Reported by category",
    Calibration_slope = round(external_slope, 3),
    Brier_score = round(brier_overall, 3),
    Brier_definition = "Multiclass"
  ),
  data.frame(
    Scope = calibration_valid$Category,
    Patients_or_events = calibration_valid$Events,
    Calibration_intercept = sprintf(
      "%.3f", calibration_valid$Calibration_intercept
    ),
    Calibration_slope = round(calibration_valid$Calibration_slope, 3),
    Brier_score = round(unname(brier_by_class[class_codes]), 3),
    Brier_definition = "One-vs-rest"
  )
)

print(temporal_validation_table, row.names = FALSE)

writexl::write_xlsx(
  list("Temporal validation" = temporal_validation_table),
  path = "Temporal_Validation_Metrics_Updated.xlsx"
)
