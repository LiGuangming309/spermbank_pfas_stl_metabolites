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

# ###取出精浆代谢与端粒长度相关的代谢物
# meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q.xlsx")
# 
# meta_STL %>% 
#   filter(q_value<0.2) -> meta_STL_1
# 
# meta_STL_1
# 
# names(meta_STL_1)


###取出PFOA混合暴露与精液代谢相关的代谢物
pfas_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_only_p_value_mix_liner.xlsx")

pfas_meta



# PFOS.beta , PFOA.beta, PFDA.beta, PFUdA.beta, L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PFOA.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
# names(meta_STL_1)
names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))

library("magrittr")

#>>>>>>>>>>>>>>>>>>>>>>>>>>所有化合物的火山图<<<<<<<<<<<<<<<<<<<<<<<<<<<


marker_color <-
  c(
    "Up" = ggsci::pal_futurama()(n = 9)[6],
    "Down" = ggsci::pal_futurama()(n = 9)[3],
    "No" = ggsci::pal_futurama()(n = 9)[9]
  )

show_col(ggsci::pal_futurama()(n = 9))


pfas_meta_STL_p_kegg


names(pfas_meta_STL_p_kegg)

library(showtext)
library(sysfonts)

font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置



####利用全部的代谢化合物绘制火山图
pfas_meta %>% 
  filter(var.names.x == "PFOA.beta") -> pfas_meta_stl_all_2 


pfas_meta_stl_all_2  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_2


pfas_meta_STL_p_kegg_2


names(pfas_meta_STL_p_kegg_2)


# 
# 
# ###（1）取出STL数据 
# pfas_meta_STL_p_kegg_2 %>% 
#   select(metabolite,Estimate,CI2.5,CI97.5,name,p.value,q_value.y) %>% 
#   mutate(group=rep("STL",271))-> pfas_meta_STL_p_kegg_2_stl
# names(pfas_meta_STL_p_kegg_2_stl)



###（2）取出PFAS数据
pfas_meta_STL_p_kegg_2 %>% 
  select(metabolite,estimate.x,estimate_lcl.x,estimate_ucl.x,name,p_value,q_value) -> pfas_meta_STL_p_kegg_2_pfas

names(pfas_meta_STL_p_kegg_2_pfas) 




pfas_meta_STL_p_kegg_2_pfas %>% 
  # group_by(group) %>% 
  mutate(marker=case_when(
    p_value < 0.05 & estimate.x  > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    p_value < 0.05 & estimate.x < 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate.x, y = -log10(p_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(p_value, 10),color = marker),alpha = 0.7) +
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2,size = 1) +
  labs(x = "Differential metabolites", y = "-log10(p value)")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )+
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE)  -> volcano_plot_all




volcano_plot_all

# dev.off()


# 
# 
# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/volcano_pfas_meta_STL_merger_line_PFOS.pdf",width = 7, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_all #出图

par(opar)
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




#################################PFOS#########################################


# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PFOS.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
names(meta_STL_1)
names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))


library("magrittr")

#>>>>>>>>>>>>>>>>>>>>>>>>>>所有化合物的火山图<<<<<<<<<<<<<<<<<<<<<<<<<<<


marker_color <-
  c(
    "Up" = ggsci::pal_futurama()(n = 9)[6],
    "Down" = ggsci::pal_futurama()(n = 9)[3],
    "No" = ggsci::pal_futurama()(n = 9)[9]
  )

show_col(ggsci::pal_futurama()(n = 9))


pfas_meta_STL_p_kegg


names(pfas_meta_STL_p_kegg)

library(showtext)
library(sysfonts)

font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置



####利用全部的代谢化合物绘制火山图
pfas_meta %>% 
  filter(var.names.x == "PFOS.beta") -> pfas_meta_stl_all_2 


pfas_meta_stl_all_2  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_2


pfas_meta_STL_p_kegg_2


names(pfas_meta_STL_p_kegg_2)


# 
# 
# ###（1）取出STL数据 
# pfas_meta_STL_p_kegg_2 %>% 
#   select(metabolite,Estimate,CI2.5,CI97.5,name,p.value,q_value.y) %>% 
#   mutate(group=rep("STL",271))-> pfas_meta_STL_p_kegg_2_stl
# names(pfas_meta_STL_p_kegg_2_stl)



###（2）取出PFAS数据
pfas_meta_STL_p_kegg_2 %>% 
  select(metabolite,estimate.x,estimate_lcl.x,estimate_ucl.x,name,p_value,q_value) -> pfas_meta_STL_p_kegg_2_pfas

names(pfas_meta_STL_p_kegg_2_pfas) 




pfas_meta_STL_p_kegg_2_pfas %>% 
  # group_by(group) %>% 
  mutate(marker=case_when(
    p_value < 0.05 & estimate.x  > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    p_value < 0.05 & estimate.x < 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate.x, y = -log10(p_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(p_value, 10),color = marker),alpha = 0.7) +
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2,size = 1) +
  labs(x = "Differential metabolites", y = "-log10(p value)")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )+
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE)  -> volcano_plot_all




volcano_plot_all

# dev.off()


# 
# 
# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/volcano_pfas_meta_STL_merger_line_PFOS.pdf",width = 7, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_all #出图

par(opar)
dev.off()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>





#################################PFDA.beta#########################################


# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PFDA.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
names(meta_STL_1)
names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))



library("magrittr")

#>>>>>>>>>>>>>>>>>>>>>>>>>>所有化合物的火山图<<<<<<<<<<<<<<<<<<<<<<<<<<<


marker_color <-
  c(
    "Up" = ggsci::pal_futurama()(n = 9)[6],
    "Down" = ggsci::pal_futurama()(n = 9)[3],
    "No" = ggsci::pal_futurama()(n = 9)[9]
  )

show_col(ggsci::pal_futurama()(n = 9))


pfas_meta_STL_p_kegg


names(pfas_meta_STL_p_kegg)

library(showtext)
library(sysfonts)

font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置



####利用全部的代谢化合物绘制火山图
pfas_meta %>% 
  filter(var.names.x == "PFDA.beta") -> pfas_meta_stl_all_2 


pfas_meta_stl_all_2  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_2


pfas_meta_STL_p_kegg_2


names(pfas_meta_STL_p_kegg_2)


# 
# 
# ###（1）取出STL数据 
# pfas_meta_STL_p_kegg_2 %>% 
#   select(metabolite,Estimate,CI2.5,CI97.5,name,p.value,q_value.y) %>% 
#   mutate(group=rep("STL",271))-> pfas_meta_STL_p_kegg_2_stl
# names(pfas_meta_STL_p_kegg_2_stl)



###（2）取出PFAS数据
pfas_meta_STL_p_kegg_2 %>% 
  select(metabolite,estimate.x,estimate_lcl.x,estimate_ucl.x,name,p_value,q_value) -> pfas_meta_STL_p_kegg_2_pfas

names(pfas_meta_STL_p_kegg_2_pfas) 




pfas_meta_STL_p_kegg_2_pfas %>% 
  # group_by(group) %>% 
  mutate(marker=case_when(
    p_value < 0.05 & estimate.x  > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    p_value < 0.05 & estimate.x < 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate.x, y = -log10(p_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(p_value, 10),color = marker),alpha = 0.7) +
  scale_color_manual(values=c("#525252","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2,size = 1) +
  labs(x = "Differential metabolites", y = "-log10(p value)")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )+
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE)  -> volcano_plot_all




volcano_plot_all

# dev.off()


# 
# 
# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/volcano_pfas_meta_STL_merger_line_PFDA.pdf",width = 7, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_all #出图

par(opar)
dev.off()





#################################PFUdA.beta#########################################


# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PFUdA.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
names(meta_STL_1)
names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))



library("magrittr")

#>>>>>>>>>>>>>>>>>>>>>>>>>>所有化合物的火山图<<<<<<<<<<<<<<<<<<<<<<<<<<<


marker_color <-
  c(
    "Up" = ggsci::pal_futurama()(n = 9)[6],
    "Down" = ggsci::pal_futurama()(n = 9)[3],
    "No" = ggsci::pal_futurama()(n = 9)[9]
  )

show_col(ggsci::pal_futurama()(n = 9))


pfas_meta_STL_p_kegg


names(pfas_meta_STL_p_kegg)

library(showtext)
library(sysfonts)

font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置



####利用全部的代谢化合物绘制火山图
pfas_meta %>% 
  filter(var.names.x == "PFUdA.beta") -> pfas_meta_stl_all_2 


pfas_meta_stl_all_2  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_2


pfas_meta_STL_p_kegg_2


names(pfas_meta_STL_p_kegg_2)


# 
# 
# ###（1）取出STL数据 
# pfas_meta_STL_p_kegg_2 %>% 
#   select(metabolite,Estimate,CI2.5,CI97.5,name,p.value,q_value.y) %>% 
#   mutate(group=rep("STL",271))-> pfas_meta_STL_p_kegg_2_stl
# names(pfas_meta_STL_p_kegg_2_stl)



###（2）取出PFAS数据
pfas_meta_STL_p_kegg_2 %>% 
  select(metabolite,estimate.x,estimate_lcl.x,estimate_ucl.x,name,p_value,q_value) -> pfas_meta_STL_p_kegg_2_pfas

names(pfas_meta_STL_p_kegg_2_pfas) 




pfas_meta_STL_p_kegg_2_pfas %>% 
  # group_by(group) %>% 
  mutate(marker=case_when(
    p_value < 0.05 & estimate.x  > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    p_value < 0.05 & estimate.x < 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate.x, y = -log10(p_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(p_value, 10),color = marker),alpha = 0.7) +
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2,size = 1) +
  labs(x = "Differential metabolites", y = "-log10(p value)")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )+
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE)  -> volcano_plot_all




volcano_plot_all

# dev.off()


# 
# 
# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/volcano_pfas_meta_STL_merger_line_PFUdA.pdf",width = 7, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_all #出图

par(opar)
dev.off()









#################################L_PFHxS.beta#########################################


# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "L_PFHxS.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
names(meta_STL_1)
names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))



library("magrittr")

#>>>>>>>>>>>>>>>>>>>>>>>>>>所有化合物的火山图<<<<<<<<<<<<<<<<<<<<<<<<<<<


marker_color <-
  c(
    "Up" = ggsci::pal_futurama()(n = 9)[6],
    "Down" = ggsci::pal_futurama()(n = 9)[3],
    "No" = ggsci::pal_futurama()(n = 9)[9]
  )

show_col(ggsci::pal_futurama()(n = 9))


pfas_meta_STL_p_kegg


names(pfas_meta_STL_p_kegg)

library(showtext)
library(sysfonts)

font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置



####利用全部的代谢化合物绘制火山图
pfas_meta %>% 
  filter(var.names.x == "L_PFHxS.beta") -> pfas_meta_stl_all_2 


pfas_meta_stl_all_2  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_2


pfas_meta_STL_p_kegg_2


names(pfas_meta_STL_p_kegg_2)


# 
# 
# ###（1）取出STL数据 
# pfas_meta_STL_p_kegg_2 %>% 
#   select(metabolite,Estimate,CI2.5,CI97.5,name,p.value,q_value.y) %>% 
#   mutate(group=rep("STL",271))-> pfas_meta_STL_p_kegg_2_stl
# names(pfas_meta_STL_p_kegg_2_stl)



###（2）取出PFAS数据
pfas_meta_STL_p_kegg_2 %>% 
  select(metabolite,estimate.x,estimate_lcl.x,estimate_ucl.x,name,p_value,q_value) -> pfas_meta_STL_p_kegg_2_pfas

names(pfas_meta_STL_p_kegg_2_pfas) 




pfas_meta_STL_p_kegg_2_pfas %>% 
  # group_by(group) %>% 
  mutate(marker=case_when(
    p_value < 0.05 & estimate.x  > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    p_value < 0.05 & estimate.x < 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate.x, y = -log10(p_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(p_value, 10),color = marker),alpha = 0.7) +
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2,size = 1) +
  labs(x = "Differential metabolites", y = "-log10(p value)")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )+
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE)  -> volcano_plot_all




volcano_plot_all

# dev.off()


# 
# 
# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/volcano_pfas_meta_STL_merger_line_L_PFHxS.pdf",width = 7, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_all #出图

par(opar)
dev.off()







#################################PF3ONS_9CL.beta#########################################


# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PF3ONS_9CL.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
names(meta_STL_1)
names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))



library("magrittr")

#>>>>>>>>>>>>>>>>>>>>>>>>>>所有化合物的火山图<<<<<<<<<<<<<<<<<<<<<<<<<<<


marker_color <-
  c(
    "Up" = ggsci::pal_futurama()(n = 9)[6],
    "Down" = ggsci::pal_futurama()(n = 9)[3],
    "No" = ggsci::pal_futurama()(n = 9)[9]
  )

show_col(ggsci::pal_futurama()(n = 9))


pfas_meta_STL_p_kegg


names(pfas_meta_STL_p_kegg)

library(showtext)
library(sysfonts)

font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置



####利用全部的代谢化合物绘制火山图
pfas_meta %>% 
  filter(var.names.x == "PF3ONS_9CL.beta") -> pfas_meta_stl_all_2 


pfas_meta_stl_all_2  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_2


pfas_meta_STL_p_kegg_2


names(pfas_meta_STL_p_kegg_2)


# 
# 
# ###（1）取出STL数据 
# pfas_meta_STL_p_kegg_2 %>% 
#   select(metabolite,Estimate,CI2.5,CI97.5,name,p.value,q_value.y) %>% 
#   mutate(group=rep("STL",271))-> pfas_meta_STL_p_kegg_2_stl
# names(pfas_meta_STL_p_kegg_2_stl)



###（2）取出PFAS数据
pfas_meta_STL_p_kegg_2 %>% 
  select(metabolite,estimate.x,estimate_lcl.x,estimate_ucl.x,name,p_value,q_value) -> pfas_meta_STL_p_kegg_2_pfas

names(pfas_meta_STL_p_kegg_2_pfas) 




pfas_meta_STL_p_kegg_2_pfas %>% 
  # group_by(group) %>% 
  mutate(marker=case_when(
    p_value < 0.05 & estimate.x  > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    p_value < 0.05 & estimate.x < 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate.x, y = -log10(p_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(p_value, 10),color = marker),alpha = 0.7) +
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2,size = 1) +
  labs(x = "Differential metabolites", y = "-log10(p value)")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )+
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE)  -> volcano_plot_all




volcano_plot_all

# dev.off()


# 
# 
# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/volcano_pfas_meta_STL_merger_line_PF3ONS_9CL.pdf",width = 7, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_all #出图

par(opar)
dev.off()





