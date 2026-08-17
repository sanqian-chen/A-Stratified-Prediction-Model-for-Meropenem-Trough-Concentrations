library(dplyr)
library(rms)

# 1. 提取训练集 -------------------------------------------------------------

train_data <- data %>%
  dplyr::filter(dataset == 1)

# 2. 设置因变量 -------------------------------------------------------------

train_data$血药浓度分组 <- ordered(
  train_data$血药浓度分组,
  levels = c(1, 2, 3)
)

# 3. CKDEPI分组作为数值型有序变量 ------------------------------------------

train_data$CKDEPI分组 <- as.numeric(as.character(train_data$CKDEPI分组))


# 4. 只保留建模需要的变量 ---------------------------------------------------

nom_data <- train_data %>%
  dplyr::select(
    血药浓度分组,
    尿素氮,
    血红蛋白,
    血小板,
    年龄,
    APTT,
    输注时间,
    CKDEPI分组
  )


# 5. 检查是否仍有常量变量 ---------------------------------------------------
# 如果所有变量unique数都大于1，就不会再出现constant警告

sapply(nom_data, function(x) length(unique(x)))


# 6. 设置rms建模环境 --------------------------------------------------------

dd <- datadist(nom_data)
options(datadist = "dd")


# 7. 拟合有序Logistic回归模型 ----------------------------------------------

fit_nom <- lrm(
  血药浓度分组 ~
    尿素氮 +
    血红蛋白 +
    血小板 +
    年龄 +
    APTT +
    输注时间 +
    CKDEPI分组,
  data = nom_data,
  x = TRUE,
  y = TRUE
)

fit_nom


# 8. 提取截距 ---------------------------------------------------------------

kint <- fit_nom$non.slopes
alpha <- coef(fit_nom)[1:kint]

alpha


# 9. 定义三个预测概率函数 ---------------------------------------------------
# lrm输出为 y>=2、y>=3，因此使用以下函数

p_y1 <- function(lp) {
  1 - plogis(alpha[1] + lp)
}

p_y2 <- function(lp) {
  plogis(alpha[1] + lp) - plogis(alpha[2] + lp)
}

p_y3 <- function(lp) {
  plogis(alpha[2] + lp)
}


# 10. 绘制列线图 ------------------------------------------------------------

nom <- nomogram(
  fit_nom,
  fun = list(
    p_y1,
    p_y2,
    p_y3
  ),
  funlabel = c(
    "Pr(Y=1)",
    "Pr(Y=2)",
    "Pr(Y=3)"
  ),
  fun.at = c(
    0.05, 0.10, 0.20, 0.30, 0.40,
    0.50, 0.60, 0.70, 0.80, 0.90
  ),
  lp = FALSE
)
nom


plot(
  nom,
  xfrac = 0.35,
  cex.var = 1.0,
  cex.axis = 0.8,
  lmgp = 0.25
)







####把血药浓度的3个组别变成英文
nom <- nomogram(
  fit_nom,
  fun = list(
    p_y1,
    p_y2,
    p_y3
  ),
  funlabel = c(
    "Pr(Subtherapeutic group)",
    "Pr(Therapeutic group)",
    "Pr(High-exposure group)"
  ),
  fun.at = c(
    0.05, 0.10, 0.20, 0.30, 0.40,
    0.50, 0.60, 0.70, 0.80, 0.90
  ),
  lp = FALSE
)

plot(
  nom,
  xfrac = 0.35,
  cex.var = 1.0,
  cex.axis = 0.8,
  lmgp = 0.25
)









# #添加风险评分表 ----------------------------------------------------------------
alpha
alpha[1]
alpha[2]
fit_nom$linear.predictors
range(fit_nom$linear.predictors)
total_score <- seq(0,160,10)
total_score <- seq(0,280,10)
prob_table <- data.frame(
  Total_points=seq(0,280,10)
)
str(nom)
names(nom)



nom$total.points
nom$total.points$x
###################################################
# 更接近nomogram实际计算
###################################################

sc <- attr(nom,"info")$sc


# nomogram最低LP
lp_min <- min(fit_nom$linear.predictors)


total_points <- seq(0,220,10)


# Points = (LP-lp_min)*sc
lp <- total_points/sc + lp_min



risk_table <- data.frame(
  
  Total_points = total_points,
  
  `Pr(<2 μg/mL)` =
    1-plogis(alpha[1]+lp),
  
  `Pr(2-8 μg/mL)` =
    plogis(alpha[1]+lp)-
    plogis(alpha[2]+lp),
  
  `Pr(>8 μg/mL)` =
    plogis(alpha[2]+lp)
  
)


risk_table[,2:4] <-
  round(risk_table[,2:4]*100,2)


risk_table
