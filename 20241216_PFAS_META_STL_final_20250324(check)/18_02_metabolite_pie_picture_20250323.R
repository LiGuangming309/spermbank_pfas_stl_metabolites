

no_source()
rm(list = ls())
# setwd(r4projects::get_project_wd())

source("R/100-tools.R")

library(tidyverse)
library(tidymass)



pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts)
# get_project_wd()
# masstools::setwd_project()

getwd() 




####>>>>>>>>>>>>>>>>>>>>>>>>>>读入胎粪代谢组学<<<<<<<<<<<<<<<<<<<<<<<


####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2

str(df_unique_2)   ###一共有357个确定的代谢物


######>>>>>>删除ICC<0.3的代谢物  

df_icc <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc_20250222_number.xlsx")
####重复测量的32个同样样本的代谢物特征的一致性

df_icc %>% 
  filter(icc>=0.3) -> df_icc_1   ###筛选出icc>0.3的代谢物特征


df_icc_1$metabolites


df_unique_2 %>% 
  filter(metabolite %in% df_icc_1$metabolites) -> df_unique_3


str(df_unique_3)

#####0.1.1大类中各个代谢成分的占比
# 使用 table() 计算每个类别的频次
freq_table <- table(df_unique_3$class_material)

# 计算每个类别的占比
percentage_Classification <- freq_table / sum(freq_table) * 100

# 将频次和占比组合成一个数据框
result_Classification <- data.frame(class_material = names(freq_table), 
                                    Frequency = as.integer(freq_table), 
                                    Percentage = round(percentage_Classification, 2))



library(ggplot2)
# install.packages("omu")
library(omu) 
library(scales)

show_col(ggsci::pal_nejm()(n = 9))

result_Classification %>%
  arrange(desc(Percentage.Freq)) -> result_Classification_order 


result_Classification_order %>% 
    dplyr::mutate(class_material=case_when(class_material=="脂质"~"Lipids", 
                                         class_material=="氨基酸"~ "Amino acids and their derivatives",
                                         class_material=="核苷酸"~"Nucleotides and their derivatives",
                                         class_material=="其他代谢物"~"Other metabolites",
                                         class_material=="糖类"~ "Carbohydrates",
                                         class_material=="有机酸及其衍生物"~"Organic acids and their derivatives"),
                  Percentage.Var1=case_when(Percentage.Var1=="脂质"~"Lipids", 
                                            Percentage.Var1=="氨基酸"~ "Amino acids and their derivatives",
                                            Percentage.Var1=="核苷酸"~"Nucleotides and their derivatives",
                                            Percentage.Var1=="其他代谢物"~"Other metabolites",
                                            Percentage.Var1=="糖类"~ "Carbohydrates",
                                            Percentage.Var1=="有机酸及其衍生物"~"Organic acids and their derivatives")
  ) -> result_Classification_order_2



# 绘制饼图
p_pie <- ggplot(result_Classification_order_2, aes(x = "", y = order(Percentage.Freq), fill = class_material)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # 将条形图转换为饼图
  scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF","#6F99ADFF")) +
  labs(title = "Pie Chart of Categories") +
  geom_text(aes(label = paste(class_material, "\n", Frequency, " (",Percentage.Freq, "%)", sep = "")), 
            position = position_stack(vjust = 0.6)) +  # 添加标签
  theme_void()  # 去除坐标轴和背景




library(ggplot2)

ggplot(result_Classification_order_2, aes(x = "", y = order(Percentage.Freq), fill = class_material)) +
  geom_bar(stat = "identity", width = 1,color = "black") +  # 添加黑色分割线
  coord_polar(theta = "y") +  # 将条形图转换为饼图
  scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF","#6F99ADFF")) +
  labs(title = "Pie Chart of Categories") +
  geom_text(aes(label = paste(class_material, "\n", Frequency, " (", Percentage.Freq, "%)", sep = "")),
            position = position_stack(vjust = 0.6), 
            family = "Arial", size = 4, color = "black") +  # 设置字体为Arial
  theme_void() +  # 去除坐标轴和背景
  theme(text = element_text(family = "Arial"))  # 全局字体设置为Arial
 



library(ggplot2)
library(ggrepel)  # 加载ggrepel包

p_pie <- ggplot(result_Classification_order_2, aes(x = "", y = order(Percentage.Freq), fill = class_material)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # 添加黑色分割线
  coord_polar(theta = "y") +  # 将条形图转换为饼图
  scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF","#6F99ADFF")) +
  # labs(title = "Pie Chart of Categories") +
  geom_text_repel(
    aes(label = paste(class_material, "\n", Frequency, " (", Percentage.Freq, "%)", sep = "")),
    position = position_stack(vjust = 0.5),  # 调整标签位置
    family = "Arial", size = 4, color = "black",  # 设置字体为Arial
    box.padding = 0.5,  # 调整标签间距
    point.padding = 0.5,  # 调整标签与点的间距
    max.overlaps = Inf  # 允许标签重叠时强制调整
  ) +
  theme_void() +  # 去除坐标轴和背景
  theme(text = element_text(family = "Arial"))  # 全局字体设置为Arial





p_pie <- ggplot(result_Classification_order_2, aes(x = "", y = order(Percentage.Freq), fill = class_material)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # 添加黑色分割线
  coord_polar(theta = "y") +  # 将条形图转换为饼图
  scale_fill_manual(values = c("Amino acids and their derivatives"="#ff7f00", 
                               "Carbohydrates"="#984ea3",
                               "Lipids"="#e41a1c", 
                               "Nucleotides and their derivatives"="#4daf4a", 
                               "Organic acids and their derivatives"= "#377eb8",
                               "Other metabolites"="#f781bf")) +
  # labs(title = "Pie Chart of Categories") +
  geom_text_repel(
    aes(label = paste(class_material, "\n", Frequency, " (", Percentage.Freq, "%)", sep = "")),
    position = position_stack(vjust = 0.5),  # 调整标签位置
    family = "Arial", size = 4, color = "black",  # 设置字体为Arial
    box.padding = 0.5,  # 调整标签间距
    point.padding = 0.5,  # 调整标签与点的间距
    max.overlaps = Inf  # 允许标签重叠时强制调整
  ) +
  theme_void() +  # 去除坐标轴和背景
  theme(text = element_text(family = "Arial"))  # 全局字体设置为Arial


p_pie

library(showtext)
library(sysfonts)

font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置



cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/describe_pie/class_metabolite.pdf",width = 12, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p_pie #出图

par(opar)
dev.off()

