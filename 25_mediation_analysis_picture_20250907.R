
# no_source() 

no_source() 
rm(list = ls())
# setwd(r4projects::get_project_wd())
source("R/100-tools.R")

library(tidyverse)
library(tidymass)



pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts,scales)


library(ggplot2)
library(dplyr)




####-------------回归系数的95%CI置信区间-------------############
# 
# 
# df_esti <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/20250714_final/median_estimate_20250711.xlsx")
# 
# 
# df_esti 
# 
# 
# df_esti$CI2.5 <- df_esti$Estimate-1.96*df_esti$Std.Error
# df_esti$CI97.5 <- df_esti$Estimate+1.96*df_esti$Std.Error
# 
# df_esti 
# 
# df_esti$CI95<- paste0(round(df_esti$Estimate,digits = 2)," [",round(df_esti$CI2.5,2),', ',round(df_esti$CI97.5,2),"]")
# 
# df_esti
# 
# 
# 
# 
# library(dplyr)
# library(tidyr)
# 
# # 假设你的数据框叫 df
# # 核心列: value...5 = group, value...6 = metabolite, Estimate, CI95
# 
# df_wide <- df_esti %>%
#   # 选择关键列
#   select(value...6, value...5, Estimate, CI95) %>%
#   # 转换为宽数据
#   pivot_wider(
#     names_from = value...5,      # a, b, c 作为列名
#     values_from = c(Estimate, CI95),  # 同时展开 Estimate 和 CI95
#     names_sep = "_"
#   )
# 
# df_wide
# 
# 
# 
# write.xlsx(df_wide,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/20250714_final/median_estimate_20250711_wide.xlsx")
# 
# 






# 数据
df <- data.frame(
  metabolite = c("12-Oxo-20-carboxy-leukotriene B4",
                 "(3E,5Z,11Z)-Pentadeca-3,5,11-trienedioylcarnitine",
                 "9-Hydroxy-16-oxooctadeca-10,12,14-trienoylcarnitine",
                 "3-Hydroxyhex-4-enoylcarnitine",
                 "13,14-dihydro-15-keto-tetranor Prostaglandin F1beta",
                 "Testosterone",
                 "L-Acetylcarnitine",
                 "O-Propanoylcarnitine",
                 "4-Hydroxy-4-methyl-2-oxoadipate",
                 "Prostaglandin E1",
                 "2-Hydroxyestradiol",
                 "Tetrahydrocortisol"),
  hr  = c(0.08815545, 0.15146368, 0.12579598, 0.09211156, 
          0.07847621, 0.08138687, 0.13714610, 0.10855098,
          0.06858654, 0.05627941, 0.11491812, 0.09433149),
  lhr = c(-0.01527781, 0.04178313, 0.00459490, 0.01626829,
          -0.01489335, -0.00738831, -0.00466834, 0.01007618,
          -0.00527702, 0.00123371, 0.02341064, 0.00476567),
  uhr = c(0.25247879, 0.40604362, 0.36776816, 0.27718537,
          0.23350071, 0.23327026, 0.43436928, 0.29072731,
          0.23377014, 0.16990541, 0.33047679, 0.25596688)
)

# 添加一列标签文本 (HR [95% CI])
df <- df %>%
  mutate(label = sprintf("%.2f [%.2f, %.2f]", hr, lhr, uhr))

# 变量顺序翻转（让第一个在上方）
df$metabolite <- factor(df$metabolite, levels = rev(df$metabolite))

# 绘图
ggplot(df, aes(x = hr, y = metabolite)) +
  geom_point(shape = 18, size = 3, color = "steelblue") +   # 方块点
  geom_errorbarh(aes(xmin = lhr, xmax = uhr), height = 0.2, color = "black") + # 水平误差棒
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + # 基准线 (logHR=0)
  geom_text(aes(label = label, x = uhr + 0.05), hjust = 0, size = 3.5,colour = "black",family = "Arial",face = "bold") + # 在右侧显示 HR (95% CI)
  labs(x = "Proportion of mediation(95% CI)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 9),       # y 轴代谢物名称
    plot.title = element_text(hjust = 0.5,colour = "black",family = "Arial",face = "bold"), # 标题居中
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
    text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
    # legend.position = c(0.25,0.15),
    legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
    axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
    # axis.line = element_line(linewidth = 0.9),
    # axis.ticks = element_line(linewidth = 0.9),
    axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
    axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold")
    # panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
    # panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
    # panel.border = element_rect(size = 0.8, colour = "black")
  ) +
  xlim(min(df$lhr) - 0.1, max(df$uhr) + 0.3)   # 给右侧留空间放标签





# 绘图
ggplot(df, aes(x = hr, y = metabolite)) +
  geom_point(shape = 18, size = 6, color = "steelblue") +   # 方块点扩大一倍
  geom_errorbarh(aes(xmin = lhr, xmax = uhr), height = 0.2, color = "black", size = 1.2) + # 水平误差棒加粗
  geom_vline(xintercept = 0, linetype = "dashed", color = "red",size=1.2) + # 基准线 (logHR=0)
  geom_text(aes(label = label, x = uhr + 0.05), 
            hjust = 0, size = 4, colour = "black", family = "Arial", face = "bold") + # 标签加粗
  labs(x = "Proportion of mediation (95% CI)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 9),       
    plot.title = element_text(hjust = 0.5, colour = "black", family = "Arial", face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    legend.title = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold")
  ) +
  xlim(min(df$lhr) - 0.1, max(df$uhr) + 0.3)   # 给右侧留空间放标签






library(stringr)

df$metabolite <- str_wrap(df$metabolite, width = 20)  # 25 表示每行最多 25 个字符，自动换行



# 固定标签位置（比最大 uhr 稍微往右）
label_x <- max(df$uhr) + 0.01  

# 绘图
p1 <- ggplot(df, aes(x = hr, y = metabolite)) +
      geom_point(shape = 18, size = 6, color = "steelblue") +   # 方块点扩大一倍
     geom_errorbarh(aes(xmin = lhr, xmax = uhr), height = 0.2, color = "black", size = 1.2) + # 水平误差棒加粗
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1.2) + # 基准线 (logHR=0)
     geom_text(aes(label = label), 
            x = label_x, hjust = 0, size = 6, 
            colour = "black", family = "Arial", face = "bold") + # 标签固定一列并加粗
  labs(x = "Proportion of mediation (95% CI)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, colour = "black", family = "Arial", face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    axis.text.y = element_text(
      size = 14, colour = "black", family = "Arial", face = "bold",
      margin = margin(r = -75)  # 15 pt，向右偏移
    ),
    legend.title = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    # axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
    axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold")
  ) +
  xlim(min(df$lhr) - 0.1, label_x + 0.3)   # 给右侧留足空间放标签


p1

# 
# 
# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/20250714_final/mediation_20250907.pdf",width = 12, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p1 #出图

par(opar)
dev.off()




################################################################################
###########                                                 ####################
###########   中介分析重新绘制森林图                        ####################
###########                                                 ####################
################################################################################









library(ggplot2)
library(dplyr)
library(readr)

# 1. 读取数据
# 假设是制表符分隔
df <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/20250714_final/mediation_20250908.xlsx")
# names(df) <- c("metabolite", "acme", "lacme", "uacme", "hr", "lhr", "uhr")

# 2. 合并 acme 列，用于标签显示

format_num <- function(x, digits = 2) {
  rounded <- round(x, digits)
  ifelse(rounded == 0 & x < 0,
         paste0("-0.", strrep("0", digits)),
         sprintf(paste0("%.", digits, "f"), rounded))
}


df <- df %>%
  mutate(
    hr_label = paste0(
      format_num(hr), " [", format_num(lhr), ", ", format_num(uhr), "]"
    ),
    acme_label = paste0(
      format_num(acme), " [", format_num(lacme), ", ", format_num(uacme), "]"
    )
  )


# 
# df <- df %>%
#   mutate(acme_label = paste0(round(acme, 2), " [", round(lacme, 2), ", ", round(uacme, 2), "]"))
# 
# df <- df %>%
#   mutate(
#     acme_label = paste0(
#       ifelse(round(acme, 2) == 0 & acme < 0,
#              "-0.00",
#              sprintf("%.2f", round(acme, 2))),
#       " [",
#       ifelse(round(lacme, 2) == 0 & lacme < 0,
#              "-0.00",
#              sprintf("%.2f", round(lacme, 2))),
#       ", ",
#       ifelse(round(uacme, 2) == 0 & uacme < 0,
#              "-0.00",
#              sprintf("%.2f", round(uacme, 2))),
#       "]"
#     )
#   )
# 
# # 3. 创建 HR+CI 标签列
# df <- df %>%
#   mutate(hr_label = paste0(round(hr, 2), " [", round(lhr, 2), ", ", round(uhr, 2), "]"))

# 4. 合并显示：HR+CI | ACME
df <- df %>%
  mutate(label = paste0(hr_label, " | ACME: ", acme_label," | ACME: ",a," | ACME: ",b," | ACME: ",c))

# 5. 设置 y 轴顺序（按 hr 从小到大排序）
df <- df %>%
  arrange(hr) %>%
  mutate(metabolite = factor(metabolite, levels = metabolite))

# 6. 标签固定在右侧
# label_x <- max(df$uhr) + 0.05

# 7. 绘图
# ggplot(df, aes(x = hr, y = metabolite)) +
#   geom_point(shape = 18, size = 4, color = "steelblue") +
#   geom_errorbarh(aes(xmin = lhr, xmax = uhr), height = 0.2, color = "black", size = 1) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#   geom_text(aes(label = label), x = label_x, hjust = 0, size = 4, family = "Arial", face = "bold") +
#   labs(x = "Effect size / HR (95% CI)", y = NULL) +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.y = element_text(size = 12, face = "bold", margin = margin(r = 15)),
#     axis.text.x = element_text(size = 12),
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   ) +
#   xlim(min(df$lhr) - 0.05, label_x + 0.2)

# 
# ggplot(df, aes(x = hr, y = metabolite)) +
#   geom_point(shape = 18, size = 6, color = "steelblue") +   # 方块点扩大一倍
#   geom_errorbarh(aes(xmin = lhr, xmax = uhr), height = 0.2, color = "black", size = 1.2) + # 水平误差棒加粗
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1.2) + # 基准线 (logHR=0)
#   geom_text(aes(label = label), 
#             x = label_x, hjust = 0, size = 6, 
#             colour = "black", family = "Arial", face = "bold") + # 标签固定一列并加粗
#   labs(x = "Proportion of mediation (95% CI)", y = NULL) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5, colour = "black", family = "Arial", face = "bold"),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.y = element_blank(),
#     legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.text.y = element_text(
#       size = 14, colour = "black", family = "Arial", face = "bold",
#       margin = margin(r = -75)  # 15 pt，向右偏移
#     ),
#     legend.title = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     # axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold")
#   ) +
#   xlim(min(df$lhr) - 0.1, label_x + 0.3)   # 给右侧留足空间放标签
# 
# 




# 
# # 1. 定义两列的x位置（右侧对齐）
# label_x1 <- max(df$uhr) + 0.05   # HR 列
# label_x2 <- label_x1 + 0.55      # ACME 列（在HR右边留间距）
# label_x3 <- label_x2 + 0.55      # ACME 列（在HR右边留间距）
# label_x4 <- label_x3 + 0.55      # ACME 列（在HR右边留间距）
# label_x5 <- label_x4 + 0.55      # ACME 列（在HR右边留间距）
# 
# 
# 
# 
# # 2. 绘图
# ggplot(df, aes(x = hr, y = metabolite)) +
#   geom_point(shape = 18, size = 6, color = "steelblue") +   # 方块点
#   geom_errorbarh(aes(xmin = lhr, xmax = uhr), height = 0.2, color = "black", size = 1.2) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1.2) +
#   
#   # 左列：HR
#   geom_text(aes(label = hr_label), 
#             x = label_x1, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   
#   # 右列：ACME
#   geom_text(aes(label = acme_label), 
#             x = label_x2, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   # 右列：ACME
#   geom_text(aes(label = a), 
#             x = label_x3, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   # 右列：ACME
#   geom_text(aes(label = b), 
#             x = label_x4, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   # 右列：ACME
#   geom_text(aes(label = c), 
#             x = label_x5, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   
#   # 可选：列标题
#   annotate("text", x = label_x1, y = length(df$metabolite) + 1,
#            label = "HR (95% CI)", hjust = 0, size = 6,
#            family = "Arial", face = "bold") +
#   annotate("text", x = label_x2, y = length(df$metabolite) + 1,
#            label = "ACME (95% CI)", hjust = 0, size = 6,
#            family = "Arial", face = "bold") +
#   
#   labs(x = "Proportion of mediation (95% CI)", y = NULL) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5, colour = "black", family = "Arial", face = "bold"),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.y = element_blank(),
#     text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.text.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold",
#                                margin = margin(r = -20)),  # 适当调整
#     axis.text.x = element_text(size = 14, family = "Arial", face = "bold"),
#     axis.title.x = element_text(size = 14, family = "Arial", face = "bold")
#   ) +
#   xlim(min(df$lhr) - 0.1, label_x2 + 0.4) +  # 给右边两列留空间
#   xlim(min(df$label_x1) - 0.1, label_x3 + 0.4)   # 给右边两列留空间
# 
# 



# 
# 
# ggplot(df, aes(x = hr, y = metabolite)) +
#   geom_point(shape = 18, size = 6, color = "steelblue") +   # 方块点扩大一倍
#   geom_errorbarh(aes(xmin = lhr, xmax = uhr), height = 0.2, color = "black", size = 1.2) + # 水平误差棒加粗
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1.2) + # 基准线 (logHR=0)
#   geom_text(aes(label = label), 
#             x = label_x, hjust = 0, size = 6, 
#             colour = "black", family = "Arial", face = "bold") + # 标签固定一列并加粗
#   labs(x = "Proportion of mediation (95% CI)", y = NULL) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5, colour = "black", family = "Arial", face = "bold"),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.y = element_blank(),
#     legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.text.y = element_text(
#       size = 14, colour = "black", family = "Arial", face = "bold",
#       margin = margin(r = -75)  # 15 pt，向右偏移
#     ),
#     legend.title = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     # axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold")
#   ) +
#   xlim(min(df$lhr) - 0.1, label_x + 0.3)   # 给右侧留足空间放标签
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 1. 定义 5 列的 x 位置
# label_x1 <- max(df$uhr) + 0.05   # HR
# label_x2 <- label_x1 + 0.95      # ACME
# label_x3 <- label_x2 + 0.95      # a
# label_x4 <- label_x3 + 0.95      # b
# label_x5 <- label_x4 + 0.95      # c
# 
# # 2. 绘图
# ggplot(df, aes(x = hr, y = metabolite)) +
#   geom_point(shape = 18, size = 6, color = "steelblue") +
#   geom_errorbarh(aes(xmin = lhr, xmax = uhr), height = 0.2,
#                  color = "black", size = 1.2) +
#   geom_vline(xintercept = 0, linetype = "dashed",
#              color = "red", size = 1.2) +
#   
#   # 五列文字
#   geom_text(aes(label = hr_label), x = label_x1, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   geom_text(aes(label = acme_label), x = label_x2, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   geom_text(aes(label = a), x = label_x3, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   geom_text(aes(label = b), x = label_x4, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   geom_text(aes(label = c), x = label_x5, hjust = 0, size = 6,
#             colour = "black", family = "Arial", face = "bold") +
#   
#   # 列标题
#   annotate("text", x = label_x1, y = length(df$metabolite) + 1,
#            label = "Proportion of mediation \nMasking effect (95% CI)", hjust = 0, size = 6,
#            family = "Arial", face = "bold") +
#   annotate("text", x = label_x2, y = length(df$metabolite) + 1,
#            label = "ACME (95% CI)", hjust = 0, size = 6,
#            family = "Arial", face = "bold") +
#   annotate("text", x = label_x3, y = length(df$metabolite) + 1,
#            label = "a", hjust = 0, size = 6,
#            family = "Arial", face = "bold") +
#   annotate("text", x = label_x4, y = length(df$metabolite) + 1,
#            label = "b", hjust = 0, size = 6,
#            family = "Arial", face = "bold") +
#   annotate("text", x = label_x5, y = length(df$metabolite) + 1,
#            label = "c", hjust = 0, size = 6,
#            family = "Arial", face = "bold") +
#   
#   labs(x = "Proportion of mediation (95% CI)", y = NULL) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5, colour = "black",
#                               family = "Arial", face = "bold"),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.y = element_blank(),
#     text = element_text(size = 14, colour = "black",
#                         family = "Arial", face = "bold"),
#     axis.text.y = element_text(size = 14, colour = "black",
#                                family = "Arial", face = "bold",
#                                margin = margin(r = -20)),
#     axis.text.x = element_text(size = 14, family = "Arial", face = "bold"),
#     axis.title.x = element_text(size = 14, family = "Arial", face = "bold")
#   ) +
#   # 👇 注意：只保留一个 xlim，覆盖到最后一列
#   xlim(min(df$lhr) - 0.1, label_x5 + 0.4)
# 
# 
# 
# 
# 
# 




# 定义 5 列的等距位置
# label_positions <- seq(max(df$uhr) + 0.6, max(df$uhr) + 4.8, by = 1.0)  
# 你可以调节起始值 0.8 和步长 1.0 来控制整体位置和间距

offset <- 0.70   # 调小越往左
spacing <- 1.5  # 列之间的间距
label_positions <- seq(max(df$uhr) + offset,
                       max(df$uhr) + offset + spacing * 4,
                       by = spacing)



# 绘图
ggplot(df, aes(x = hr, y = metabolite)) +
  geom_point(shape = 18, size = 6, color = "steelblue") +
  geom_errorbarh(aes(xmin = lhr, xmax = uhr), height = 0.2,
                 color = "black", size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "red", size = 1.2) +
  
  # 五列文字（改用向量 & 居中对齐）
  geom_text(aes(label = hr_label), x = label_positions[1], hjust = 0.5, size = 6,
            colour = "black", family = "Arial", face = "bold") +
  geom_text(aes(label = acme_label), x = label_positions[2], hjust = 0.5, size = 6,
            colour = "black", family = "Arial", face = "bold") +
  geom_text(aes(label = a), x = label_positions[3], hjust = 0.5, size = 6,
            colour = "black", family = "Arial", face = "bold") +
  geom_text(aes(label = b), x = label_positions[4], hjust = 0.5, size = 6,
            colour = "black", family = "Arial", face = "bold") +
  geom_text(aes(label = c), x = label_positions[5], hjust = 0.5, size = 6,
            colour = "black", family = "Arial", face = "bold") +
  
  # 列标题（同样居中） Proportion of mediation \n Masking effect (95% CI)
  annotate("text", x = label_positions[1], y = length(df$metabolite) + 1,
           label = "(95% CI)", hjust = 0.5, size = 6,
           family = "Arial", face = "bold") +
  annotate("text", x = label_positions[2], y = length(df$metabolite) + 1,
           label = "ACME (95% CI)", hjust = 0.5, size = 6,
           family = "Arial", face = "bold") +
  annotate("text", x = label_positions[3], y = length(df$metabolite) + 1,
           label = "a", hjust = 0.5, size = 6,
           family = "Arial", face = "bold") +
  annotate("text", x = label_positions[4], y = length(df$metabolite) + 1,
           label = "b", hjust = 0.5, size = 6,
           family = "Arial", face = "bold") +
  annotate("text", x = label_positions[5], y = length(df$metabolite) + 1,
           label = "c", hjust = 0.5, size = 6,
           family = "Arial", face = "bold") +
  
  # labs(x = "Proportion of mediation (95% CI)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, colour = "black",
                              family = "Arial", face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    text = element_text(size = 14, colour = "black",
                        family = "Arial", face = "bold"),
    axis.text.y = element_text(size = 14, colour = "black",
                               family = "Arial", face = "bold",
                               margin = margin(r = -35)),
    axis.text.x = element_text(size = 14, family = "Arial", face = "bold"),
    axis.title.x = element_text(size = 14, family = "Arial", face = "bold")
  ) +
  # 👇 覆盖到最后一列
  xlim(min(df$lhr) - 0.1, label_positions[5] + 0.4)  -> p2



p2

# 
# 
# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/20250714_final/mediation_20250908_final.pdf",width = 18, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p2 #出图

par(opar)
dev.off()




