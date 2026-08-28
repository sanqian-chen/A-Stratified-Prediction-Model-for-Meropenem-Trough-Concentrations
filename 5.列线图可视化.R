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









# =============================================================================
# 美罗培南血药浓度分组电子计算器——最终修正版
# 适用于rms::lrm建立的三分类有序Logistic回归模型fit_nom
# =============================================================================

# 首次使用时安装
# install.packages("openxlsx")

library(openxlsx)


# =============================================================================
# 一、提取模型参数
# =============================================================================

var_names <- c(
  "尿素氮",
  "血红蛋白",
  "血小板",
  "年龄",
  "APTT",
  "输注时间",
  "CKDEPI分组"
)

# 检查模型是否存在
if (!exists("fit_nom")) {
  stop("未检测到fit_nom，请先运行有序Logistic回归模型代码。")
}

# 检查建模数据是否存在
if (!exists("nom_data")) {
  stop("未检测到nom_data，请先运行建模数据整理代码。")
}

# 提取全部模型参数
coef_all <- coef(fit_nom)

# 有序三分类结局应当有两个截距
kint <- fit_nom$non.slopes

if (kint != 2) {
  stop("当前模型截距数量不是2，请检查结局是否为三分类有序变量。")
}

# 两个截距
alpha <- coef_all[seq_len(kint)]

# 7个预测变量的回归系数
beta <- coef_all[var_names]

if (any(is.na(beta))) {
  stop(
    paste0(
      "以下变量未能从fit_nom中正确提取：",
      paste(var_names[is.na(beta)], collapse = "、")
    )
  )
}

# 显示模型参数
model_parameters <- data.frame(
  参数 = c(names(alpha), names(beta)),
  估计值 = as.numeric(c(alpha, beta)),
  row.names = NULL
)

print(model_parameters)


# =============================================================================
# 二、验证手动概率公式与fit_nom是否完全一致
# =============================================================================

complete_index <- which(
  complete.cases(nom_data[, var_names])
)[1]

if (is.na(complete_index)) {
  stop("nom_data中没有7个变量均完整的患者，无法验证概率公式。")
}

test_patient <- nom_data[
  complete_index,
  var_names,
  drop = FALSE
]

test_values <- as.numeric(
  test_patient[1, var_names]
)

# 线性预测值
lp_test <- sum(
  test_values * as.numeric(beta)
)

# rms::lrm模型的两个累计概率
p_ge2_test <- plogis(
  alpha[1] + lp_test
)

p_ge3_test <- plogis(
  alpha[2] + lp_test
)

# 3个类别的概率
p_manual <- c(
  "亚治疗组" = 1 - p_ge2_test,
  "治疗范围组" = p_ge2_test - p_ge3_test,
  "高暴露组" = p_ge3_test
)

# fit_nom直接预测结果
p_model <- as.numeric(
  predict(
    fit_nom,
    newdata = test_patient,
    type = "fitted.ind"
  )
)

names(p_model) <- c(
  "亚治疗组",
  "治疗范围组",
  "高暴露组"
)

validation_result <- data.frame(
  血药浓度组别 = names(p_manual),
  手动公式计算概率 = as.numeric(p_manual),
  fit_nom预测概率 = as.numeric(p_model),
  绝对差值 = abs(
    as.numeric(p_manual) -
      as.numeric(p_model)
  )
)

print(validation_result)

if (max(validation_result$绝对差值) > 1e-8) {
  stop("概率公式与fit_nom预测结果不一致，请停止生成计算器。")
} else {
  message("概率公式验证通过，可以生成Excel计算器。")
}


# =============================================================================
# 三、提取训练集观察范围
# =============================================================================

get_range_text <- function(x, digits = 1) {
  
  x_range <- range(
    x,
    na.rm = TRUE
  )
  
  paste0(
    round(x_range[1], digits),
    "～",
    round(x_range[2], digits)
  )
}

input_ranges <- c(
  get_range_text(nom_data$尿素氮, 2),
  get_range_text(nom_data$血红蛋白, 1),
  get_range_text(nom_data$血小板, 1),
  get_range_text(nom_data$年龄, 0),
  get_range_text(nom_data$APTT, 1),
  get_range_text(nom_data$输注时间, 1),
  "1～5"
)

# 用于Excel输入限制
range_list <- list(
  range(nom_data$尿素氮, na.rm = TRUE),
  range(nom_data$血红蛋白, na.rm = TRUE),
  range(nom_data$血小板, na.rm = TRUE),
  range(nom_data$年龄, na.rm = TRUE),
  range(nom_data$APTT, na.rm = TRUE)
)


# =============================================================================
# 四、创建Excel工作簿
# =============================================================================

wb <- createWorkbook(
  creator = "Meropenem TDM Prediction Study"
)

addWorksheet(
  wb,
  sheetName = "临床计算器",
  gridLines = FALSE
)

addWorksheet(
  wb,
  sheetName = "模型参数",
  gridLines = FALSE
)

# 用于保存下拉菜单选项
addWorksheet(
  wb,
  sheetName = "下拉选项",
  gridLines = FALSE,
  visible = FALSE
)


# =============================================================================
# 五、创建下拉菜单数据源
# =============================================================================

# 输注时间
writeData(
  wb,
  sheet = "下拉选项",
  x = data.frame(
    输注时间 = c(1, 2, 3)
  ),
  startCol = 1,
  startRow = 1,
  colNames = FALSE
)

# CKD-EPI分组
ckdepi_options <- c(
  "1: >50 mL/min",
  "2: >25至≤50 mL/min",
  "3: ≥10至≤25 mL/min",
  "4: <10 mL/min",
  "5: CRRT"
)

writeData(
  wb,
  sheet = "下拉选项",
  x = data.frame(
    CKDEPI分组 = ckdepi_options
  ),
  startCol = 2,
  startRow = 1,
  colNames = FALSE
)


# =============================================================================
# 六、设置Excel样式
# =============================================================================

title_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 16,
  fontColour = "#FFFFFF",
  fgFill = "#1F4E78",
  textDecoration = "bold",
  halign = "center",
  valign = "center"
)

subtitle_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 11,
  fontColour = "#FFFFFF",
  fgFill = "#5B9BD5",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "TopBottomLeftRight",
  borderColour = "#FFFFFF"
)

header_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 11,
  fontColour = "#FFFFFF",
  fgFill = "#4472C4",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "TopBottomLeftRight"
)

input_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 11,
  fgFill = "#FFF2CC",
  halign = "center",
  valign = "center",
  border = "TopBottomLeftRight",
  borderColour = "#A6A6A6",
  locked = FALSE
)

normal_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 11,
  valign = "center",
  border = "TopBottomLeftRight",
  borderColour = "#D9D9D9"
)

parameter_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 10,
  numFmt = "0.000000000",
  border = "TopBottomLeftRight",
  borderColour = "#D9D9D9"
)

output_label_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 11,
  fontColour = "#FFFFFF",
  fgFill = "#70AD47",
  textDecoration = "bold",
  halign = "left",
  valign = "center",
  border = "TopBottomLeftRight"
)

probability_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 12,
  fgFill = "#E2F0D9",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  numFmt = "0.00%",
  border = "TopBottomLeftRight"
)

result_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 13,
  fontColour = "#C00000",
  fgFill = "#FCE4D6",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "TopBottomLeftRight"
)

note_style <- createStyle(
  fontName = "微软雅黑",
  fontSize = 10,
  fontColour = "#666666",
  wrapText = TRUE,
  valign = "top"
)


# =============================================================================
# 七、创建“模型参数”工作表
# =============================================================================

mergeCells(
  wb,
  sheet = "模型参数",
  cols = 1:3,
  rows = 1
)

writeData(
  wb,
  sheet = "模型参数",
  x = "有序Logistic回归模型参数",
  startCol = 1,
  startRow = 1
)

addStyle(
  wb,
  sheet = "模型参数",
  style = title_style,
  rows = 1,
  cols = 1:3,
  gridExpand = TRUE
)

# 7个回归系数
parameter_table <- data.frame(
  预测变量 = var_names,
  回归系数 = as.numeric(beta),
  单位或编码 = c(
    "mmol/L",
    "g/L",
    "×10^9/L",
    "岁",
    "s",
    "h",
    "1～5级"
  )
)

writeData(
  wb,
  sheet = "模型参数",
  x = parameter_table,
  startCol = 1,
  startRow = 2,
  colNames = TRUE,
  headerStyle = header_style
)

# 参数位于B3:B9
addStyle(
  wb,
  sheet = "模型参数",
  style = parameter_style,
  rows = 3:9,
  cols = 1:3,
  gridExpand = TRUE
)

# 两个截距
intercept_table <- data.frame(
  参数 = c("y>=2", "y>=3"),
  估计值 = as.numeric(alpha),
  说明 = c(
    "P(Y≥2)的截距",
    "P(Y≥3)的截距"
  )
)

writeData(
  wb,
  sheet = "模型参数",
  x = intercept_table,
  startCol = 1,
  startRow = 11,
  colNames = TRUE,
  headerStyle = header_style
)

# 截距位于B12:B13
addStyle(
  wb,
  sheet = "模型参数",
  style = parameter_style,
  rows = 12:13,
  cols = 1:3,
  gridExpand = TRUE
)

# CKD-EPI分组说明
ckdepi_table <- data.frame(
  CKDEPI分组 = 1:5,
  定义 = c(
    ">50 mL/min",
    ">25至≤50 mL/min",
    "≥10至≤25 mL/min",
    "<10 mL/min",
    "接受CRRT"
  ),
  模型编码 = 1:5
)

writeData(
  wb,
  sheet = "模型参数",
  x = ckdepi_table,
  startCol = 1,
  startRow = 16,
  colNames = TRUE,
  headerStyle = header_style
)

addStyle(
  wb,
  sheet = "模型参数",
  style = normal_style,
  rows = 17:21,
  cols = 1:3,
  gridExpand = TRUE
)

writeData(
  wb,
  sheet = "模型参数",
  x = paste0(
    "模型建立日期：",
    format(Sys.Date(), "%Y-%m-%d"),
    "。模型参数由fit_nom自动提取，请勿修改。"
  ),
  startCol = 1,
  startRow = 24
)

addStyle(
  wb,
  sheet = "模型参数",
  style = note_style,
  rows = 24,
  cols = 1:3,
  gridExpand = TRUE
)


# =============================================================================
# 八、创建“临床计算器”工作表
# =============================================================================

# 标题
mergeCells(
  wb,
  sheet = "临床计算器",
  cols = 1:4,
  rows = 1
)

writeData(
  wb,
  sheet = "临床计算器",
  x = "美罗培南首次血药浓度分组预测计算器",
  startCol = 1,
  startRow = 1
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = title_style,
  rows = 1,
  cols = 1:4,
  gridExpand = TRUE
)

setRowHeights(
  wb,
  sheet = "临床计算器",
  rows = 1,
  heights = 30
)

# 副标题
mergeCells(
  wb,
  sheet = "临床计算器",
  cols = 1:4,
  rows = 2
)

writeData(
  wb,
  sheet = "临床计算器",
  x = "适用于重症成人患者首次美罗培南TDM血药浓度分组预测",
  startCol = 1,
  startRow = 2
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = subtitle_style,
  rows = 2,
  cols = 1:4,
  gridExpand = TRUE
)

# 输入表头
input_header <- data.frame(
  预测变量 = "预测变量",
  患者数值 = "请输入患者数值",
  单位 = "单位",
  训练集范围 = "训练集观察范围"
)

writeData(
  wb,
  sheet = "临床计算器",
  x = input_header,
  startCol = 1,
  startRow = 4,
  colNames = FALSE
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = header_style,
  rows = 4,
  cols = 1:4,
  gridExpand = TRUE
)

# 输入项目
input_table <- data.frame(
  预测变量 = c(
    "尿素氮",
    "血红蛋白",
    "血小板",
    "年龄",
    "APTT",
    "输注时间",
    "CKD-EPI分组"
  ),
  患者数值 = rep("", 7),
  单位 = c(
    "mmol/L",
    "g/L",
    "×10^9/L",
    "岁",
    "s",
    "h",
    "等级"
  ),
  训练集观察范围 = input_ranges,
  stringsAsFactors = FALSE
)

writeData(
  wb,
  sheet = "临床计算器",
  x = input_table,
  startCol = 1,
  startRow = 5,
  colNames = FALSE
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = normal_style,
  rows = 5:11,
  cols = c(1, 3, 4),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = input_style,
  rows = 5:11,
  cols = 2,
  gridExpand = TRUE
)


# =============================================================================
# 九、添加输入范围限制和下拉菜单
# =============================================================================

# 连续变量：
# B5尿素氮、B6血红蛋白、B7血小板、B8年龄、B9 APTT
for (i in seq_along(range_list)) {
  
  dataValidation(
    wb,
    sheet = "临床计算器",
    cols = 2,
    rows = 4 + i,
    type = "decimal",
    operator = "between",
    value = range_list[[i]],
    allowBlank = TRUE,
    showInputMsg = FALSE,
    showErrorMsg = TRUE
  )
}

# B10：输注时间下拉菜单
dataValidation(
  wb,
  sheet = "临床计算器",
  cols = 2,
  rows = 10,
  type = "list",
  value = "'下拉选项'!$A$1:$A$3",
  allowBlank = TRUE,
  showInputMsg = FALSE,
  showErrorMsg = TRUE
)

# B11：CKD-EPI分组下拉菜单
dataValidation(
  wb,
  sheet = "临床计算器",
  cols = 2,
  rows = 11,
  type = "list",
  value = "'下拉选项'!$B$1:$B$5",
  allowBlank = TRUE,
  showInputMsg = FALSE,
  showErrorMsg = TRUE
)


# =============================================================================
# 十、建立线性预测值和三分类概率公式
# =============================================================================

writeData(
  wb,
  sheet = "临床计算器",
  x = "模型线性预测值（LP）",
  startCol = 1,
  startRow = 14
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = normal_style,
  rows = 14,
  cols = 1,
  gridExpand = TRUE
)

# 模型参数工作表：
# B3:B8为前6个连续变量系数
# B9为CKD-EPI分组系数
lp_formula <- paste0(
  '=IF(OR(COUNT(B5:B10)<6,B11=""),"",',
  'SUMPRODUCT(B5:B10,',
  "'模型参数'!$B$3:$B$8)",
  '+VALUE(LEFT(B11,1))*',
  "'模型参数'!$B$9)"
)

writeFormula(
  wb,
  sheet = "临床计算器",
  x = lp_formula,
  startCol = 2,
  startRow = 14
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = normal_style,
  rows = 14,
  cols = 2,
  gridExpand = TRUE
)

# 概率标题
mergeCells(
  wb,
  sheet = "临床计算器",
  cols = 1:2,
  rows = 16
)

writeData(
  wb,
  sheet = "临床计算器",
  x = "三类血药浓度分组预测概率",
  startCol = 1,
  startRow = 16
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = header_style,
  rows = 16,
  cols = 1:2,
  gridExpand = TRUE
)

group_labels <- c(
  "亚治疗组（<2 μg/mL）",
  "治疗范围组（2～8 μg/mL）",
  "高暴露组（>8 μg/mL）"
)

writeData(
  wb,
  sheet = "临床计算器",
  x = group_labels,
  startCol = 1,
  startRow = 18,
  colNames = FALSE
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = output_label_style,
  rows = 18:20,
  cols = 1,
  gridExpand = TRUE
)

# 模型参数工作表：
# B12为y>=2截距
# B13为y>=3截距

# P(Y=1)
formula_p1 <- paste0(
  '=IF($B$14="","",',
  '1-1/(1+EXP(-(',
  "'模型参数'!$B$12+$B$14)))",
  ')'
)

# P(Y=2)
formula_p2 <- paste0(
  '=IF($B$14="","",',
  '1/(1+EXP(-(',
  "'模型参数'!$B$12+$B$14)))",
  '-1/(1+EXP(-(',
  "'模型参数'!$B$13+$B$14)))",
  ')'
)

# P(Y=3)
formula_p3 <- paste0(
  '=IF($B$14="","",',
  '1/(1+EXP(-(',
  "'模型参数'!$B$13+$B$14)))",
  ')'
)

writeFormula(
  wb,
  sheet = "临床计算器",
  x = formula_p1,
  startCol = 2,
  startRow = 18
)

writeFormula(
  wb,
  sheet = "临床计算器",
  x = formula_p2,
  startCol = 2,
  startRow = 19
)

writeFormula(
  wb,
  sheet = "临床计算器",
  x = formula_p3,
  startCol = 2,
  startRow = 20
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = probability_style,
  rows = 18:20,
  cols = 2,
  gridExpand = TRUE
)

# 概率合计
writeData(
  wb,
  sheet = "临床计算器",
  x = "概率合计",
  startCol = 1,
  startRow = 21
)

writeFormula(
  wb,
  sheet = "临床计算器",
  x = '=IF($B$14="","",SUM(B18:B20))',
  startCol = 2,
  startRow = 21
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = normal_style,
  rows = 21,
  cols = 1,
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = probability_style,
  rows = 21,
  cols = 2,
  gridExpand = TRUE
)


# =============================================================================
# 十一、输出最终预测组别
# =============================================================================

writeData(
  wb,
  sheet = "临床计算器",
  x = "模型预测组别",
  startCol = 1,
  startRow = 23
)

predicted_group_formula <- paste0(
  '=IF($B$14="","",',
  'INDEX($A$18:$A$20,',
  'MATCH(MAX($B$18:$B$20),$B$18:$B$20,0)))'
)

writeFormula(
  wb,
  sheet = "临床计算器",
  x = predicted_group_formula,
  startCol = 2,
  startRow = 23
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = result_style,
  rows = 23,
  cols = 1:2,
  gridExpand = TRUE
)

writeData(
  wb,
  sheet = "临床计算器",
  x = "最高预测概率",
  startCol = 1,
  startRow = 24
)

writeFormula(
  wb,
  sheet = "临床计算器",
  x = '=IF($B$14="","",MAX($B$18:$B$20))',
  startCol = 2,
  startRow = 24
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = normal_style,
  rows = 24,
  cols = 1,
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = probability_style,
  rows = 24,
  cols = 2,
  gridExpand = TRUE
)


# =============================================================================
# 十二、添加使用说明
# =============================================================================

instructions <- paste0(
  "使用说明：\n",
  "1. 仅在黄色单元格中输入或选择患者的7项预测指标。\n",
  "2. CKD-EPI分组：1=>50 mL/min；2=>25至≤50 mL/min；",
  "3=≥10至≤25 mL/min；4=<10 mL/min；5=CRRT。\n",
  "3. 接受CRRT的患者，无论CKD-EPI数值如何，均选择第5组。\n",
  "4. 计算器自动输出3个血药浓度组别的预测概率，",
  "概率最高的组别作为模型预测组别。\n",
  "5. 三个预测概率之和应为100%。\n",
  "6. 本工具不能替代实际治疗药物监测及临床判断。\n",
  "7. 超出训练集观察范围的患者属于模型外推，",
  "预测结果应谨慎解释。"
)

mergeCells(
  wb,
  sheet = "临床计算器",
  cols = 1:4,
  rows = 27:33
)

writeData(
  wb,
  sheet = "临床计算器",
  x = instructions,
  startCol = 1,
  startRow = 27
)

addStyle(
  wb,
  sheet = "临床计算器",
  style = note_style,
  rows = 27:33,
  cols = 1:4,
  gridExpand = TRUE
)

setRowHeights(
  wb,
  sheet = "临床计算器",
  rows = 27:33,
  heights = 22
)


# =============================================================================
# 十三、调整版式
# =============================================================================

setColWidths(
  wb,
  sheet = "临床计算器",
  cols = 1,
  widths = 28
)

setColWidths(
  wb,
  sheet = "临床计算器",
  cols = 2,
  widths = 28
)

setColWidths(
  wb,
  sheet = "临床计算器",
  cols = 3,
  widths = 14
)

setColWidths(
  wb,
  sheet = "临床计算器",
  cols = 4,
  widths = 20
)

setColWidths(
  wb,
  sheet = "模型参数",
  cols = 1,
  widths = 22
)

setColWidths(
  wb,
  sheet = "模型参数",
  cols = 2,
  widths = 24
)

setColWidths(
  wb,
  sheet = "模型参数",
  cols = 3,
  widths = 30
)

freezePane(
  wb,
  sheet = "临床计算器",
  firstActiveRow = 4
)


# =============================================================================
# 十四、保存Excel文件
# =============================================================================

output_file <- file.path(
  getwd(),
  "美罗培南血药浓度分组电子计算器_最终版.xlsx"
)

saveWorkbook(
  wb,
  file = output_file,
  overwrite = TRUE
)

message("Excel电子计算器已成功生成：")
message(
  normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )
)
