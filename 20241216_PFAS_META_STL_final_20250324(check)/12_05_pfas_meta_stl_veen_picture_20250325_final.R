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
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts)
# get_project_wd()
# masstools::setwd_project()

getwd() 


###====================读入代谢组学数据===============================

####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2

str(df_unique_2)   ###一共有378个确定的代谢物

######>>>>>>删除ICC<0.3的代谢物  

df_icc <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc_20250222_number.xlsx")
####重复测量的32个同样样本的代谢物特征的一致性

df_icc %>% 
  filter(icc>=0.3) -> df_icc_1   ###筛选出icc>0.3的代谢物特征


df_icc_1$metabolites


df_unique_2 %>% 
  filter(metabolite %in% df_icc_1$metabolites) -> df_unique_3



####>>>>>>>>读入混合物计算后的P值------------------------------
# write.xlsx(results_df_mix_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/results_df_mix_p_value.xlsx")
# 
# pfas_meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q_PFAS.xlsx")

###取出精浆代谢与端粒长度相关的代谢物
meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q.xlsx")

meta_STL %>% 
  filter(q_value<0.2) -> meta_STL_1

meta_STL_1

str(meta_STL_1)  ###  62

meta_STL_1$meta -> STL

###取出PFAS混合暴露与精液代谢相关的代谢物
pfas_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value_mix_liner.xlsx")

pfas_meta %>% 
  filter(q_value<0.2) -> pfas_meta_1 

str(pfas_meta_1)   ###79
names(meta_STL)
names(pfas_meta)

pfas_meta_1$metabolite -> PFAS






pfas_meta_1 %>% 
  left_join(meta_STL_1,by=c("metabolite"="meta")) -> pfas_meta_stl

colSums(is.na(pfas_meta_stl))

pfas_meta_stl %>% 
  filter(!Estimate=="NA") -> pfas_meta_stl_all 



str(pfas_meta_stl_all)    ###18

names(pfas_meta_stl_all)



pfas_meta_stl_all$metabolite















library(ggvenn) #加载ggvenn包
# # install.packages("ggvenn") #安装ggvenn包
# library(ggvenn) #加载ggvenn包
# # gene_set1 <- c("BIRC3","JUN", "FOS", "MAP3K14",
# #                "MAP3K8","CREB5", "CCL2", "CCL20",
# #                "CXCL1", "CXCL3","CXCL5") #输入第一个基因集
# # 
# # gene_set2 <- c("FOS","FOSB", "JUN", "JUND",
# #                "FOSL1","TNFAIP3", "CXCL1", "CXCL3",
# #                "CXCL5", "CXCL8","CXCL10") #输入第二个基因集  
# base::intersect(gene_set1, gene_set2) #获取两个基因集的交集
# #[1] "JUN"   "FOS"   "CXCL1" "CXCL3""CXCL5"
# union(gene_set1, gene_set2) #获取两个基因集的并集
# #[1] "BIRC3"   "JUN"     "FOS"     "MAP3K14""MAP3K8"  "CREB5"   "CCL2"  
# #[8] "CCL20"   "CXCL1"   "CXCL3"   "CXCL5"   "FOSB"    "JUND"    "FOSL1" 
# #[15] "TNFAIP3""CXCL8"   "CXCL10"
# setdiff(gene_set1, gene_set2) #获取gene_set1减去gene_set2的差集
# #[1] "BIRC3"   "MAP3K14" "MAP3K8"  "CREB5"   "CCL2"    "CCL20"
# setdiff(gene_set2, gene_set1) #获取gene_set2减去gene_set1的差集
# #[1] "FOSB"    "JUND"    "FOSL1"   "TNFAIP3" "CXCL8"   "CXCL10"
a <- list(`PFAS` = PFAS,
          `STL` = STL) #将基因集变成列表变量


p1 <- ggvenn(a, show_elements = TRUE,fill_color = c("blue", "green","red"),
             label_sep = "\n", stroke_size = 1.5,set_name_size = 5,
             text_size = 3)
p1



p2 <- ggvenn(a, show_elements = FALSE,fill_color = c("blue", "green","red"),
             label_sep = "\n", stroke_size = 1.5,set_name_size = 5,
             text_size = 3)
p2



library(ggvenn)

# 加载 VennDiagram 包
library(VennDiagram)

library(ggplot2)
library(ggforce)



ggplot() +
  # 第一个椭圆（旋转 -45°）
  geom_ellipse(aes(x0 = 2.4, y0 = 2, a = 2.5, b = 1.5, angle = -45), 
               fill = "blue", alpha = 0.5,size=0.2) +
  # 第二个椭圆（旋转 45°）
  geom_ellipse(aes(x0 = 4, y0 = 2, a = 2.5, b = 1.5, angle = 45), 
               fill = "red", alpha = 0.5,size=0.2) +
  # 添加集合名称
  annotate("text", x = 1.4, y = 4.5, label = "PFAS", size = 5, color = "black",family = "Arial",face = "bold") +
  annotate("text", x = 4.8, y = 4.5, label = "STL", size = 5, color = "black",family = "Arial",face = "bold") +
  # annotate("text", x = 3.5, y = 2.8, label = "Intersection", size = 5, color = "black") +
  
  annotate("text", x = 1.6, y = 2.8, label = "61 \n (49.6%)", size = 5, color = "black",family = "Arial",face = "bold") +
  annotate("text", x = 4.6, y = 2.8, label = "44 \n (35.8%)", size = 5, color = "black",family = "Arial",face = "bold") +
  
  # 添加重叠部分的标签
  annotate("text", x = 3.2, y = 1.5, label = "18 \n (14.6%)", size = 5, color = "black",family = "Arial",face = "bold")+
  theme_void()+
  theme(
    text = element_text(family = "Arial",face = "bold"),  # 设置所有字体为 Arial
    plot.margin = margin(10, 10, 10, 10)    # 调整边距
  )  -> VEEN

VEEN
  
  # 
  # 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/VEEN_meeting_mediator.pdf",width = 7, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
VEEN #出图

par(opar)
dev.off()












