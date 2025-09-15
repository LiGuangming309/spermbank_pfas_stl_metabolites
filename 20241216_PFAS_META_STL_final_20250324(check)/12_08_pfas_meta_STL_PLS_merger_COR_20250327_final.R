
# no_source()
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







####>>>>>>>>>>>>>>>>>>>一、读入所有变量PFAS和精浆代谢组学关联数据进行合并绘图<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 


read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PFOA.xlsx") %>% 
  dplyr::rename("PFAS"="PFOA") -> temp_data_PFOA

read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PFOS.xlsx") %>% 
  dplyr::rename("PFAS"="PFOS") -> temp_data_PFOS

read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PFDA.xlsx")%>% 
  dplyr::rename("PFAS"="PFDA") -> temp_data_PFDA

read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PFUdA.xlsx")%>% 
  dplyr::rename("PFAS"="PFUdA") -> temp_data_PFUdA


read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_L_PFHxS.xlsx")%>% 
  dplyr::rename("PFAS"="L_PFHxS") -> temp_data_L_PFHxS


read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PF3ONS_9CL.xlsx")%>% 
  dplyr::rename("PFAS"="PF3ONS_9CL") -> temp_data_PF3ONS_9CL

read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_mix_PFAS_data.xlsx")%>% 
  dplyr::rename("PFAS"="PFAS") -> temp_data_PFAS


# 
# write.xlsx(temp_data,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_mix_PFAS_data.xlsx")



rbind(temp_data_PFOA,
      temp_data_PFOS,
      temp_data_PFDA,
      temp_data_PFUdA,
      temp_data_L_PFHxS,
      temp_data_PF3ONS_9CL,
      temp_data_PFAS) ->temp_data_all 

temp_data_all$class

RColorBrewer::brewer.pal(n = 12, name = "Set3")


omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "6:2 Cl-PFESA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6],
    "PFAS"=RColorBrewer::brewer.pal(n = 12, name = "Set3")[7]
  )



# ,
# gut_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[7],
# skin_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[8],
# oral_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[10],
# nasal_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[11]

temp_data_all$class



temp_data_all  %>%
  mutate(class=case_when(
    class=="PFOA"~"PFOA",
    class=="PFOS"~"PFOS",
    class=="PFDA"~"PFDA",
    class=="PFUdA"~"PFUdA",
    class=="L_PFHxS"~"PFHxS",
    class=="PF3ONS_9CL"~"6:2 Cl-PFESA",
    class=="PFAS"~"PFAS")) %>%  # 倒序排序class
  ggplot(aes(PFAS,PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  theme_bw() +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  scale_color_manual(values = omics_color) +
  labs(x = "PFAS", y = "PC1")+
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.5, linetype = "dashed"),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.5, linetype = "dashed") # 次网格线粗细,
        # panel.border = element_rect(size = 1, colour = "black"), # 设置边框线的粗细和颜色
        # plot.margin = unit(c(1, 1, 1, 1), "cm")  # 设置整个图形的外边距
  )+
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 2) -> p_all

p_all


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

###导出保存为pdf
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_pfas_meta_cor_2.pdf",width =12, height =8, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p_all
# par(opar)
dev.off()


temp_data_all %>%
  group_by(class) %>%
  dplyr::summarise(cor = cor(PFAS,PC, method = "spearman"))




temp_data_all %>% 
  filter(class=="PFOA") -> temp_data_all_pfoa
  
cor.test(temp_data_all_pfoa$PFAS, temp_data_all_pfoa$PC, method = "spearman")

summary(lm(temp_data_all_pfoa$PC~temp_data_all_pfoa$PFAS))


####绝对值回归系数与p值
library(plyr)
temp_data_all %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PFAS, x$PC, method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()





########将代谢组学与端粒长度进行关联分析 

read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/temp_data_STL.xlsx")  -> temp_data_meta_stl


temp_data_meta_stl



omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "6:2 Cl-PFESA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6],
    "STL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[7]
  )


###回归系数
lm(STL~PC,data=temp_data_meta_stl) -> line1

summary(line1)

lm(PC~STL,data=temp_data_meta_stl) -> line1

summary(line1)

1/0.46669


temp_data_meta_stl  %>%
  ggplot(aes(STL,PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  theme_bw() +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  scale_color_manual(values = omics_color) +
  labs(x = "PC1", y = "STL" )+
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.5, linetype = "dashed"),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.5, linetype = "dashed") # 次网格线粗细,
        # panel.border = element_rect(size = 1, colour = "black"), # 设置边框线的粗细和颜色
        # plot.margin = unit(c(1, 1, 1, 1), "cm")  # 设置整个图形的外边距
  )+
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 1) -> p_stl


p_stl


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

###导出保存为pdf
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_meta_STL_cor_2.pdf",width =4, height =4, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p_stl
# par(opar)
dev.off()


temp_data_meta_stl  %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(STL,PC, method = "spearman"))



####绝对值回归系数与p值
library(plyr)

temp_data_meta_stl %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$STL,x$PC,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()





########<<<<<<<<<<<<<<<<<<<合并PFAS和STL与精浆代谢组学数据>>>>>>>>>>>>>>>>>>>
library(ggplot2)
library(gridExtra)

# # 使用 grid.arrange() 拼接
# grid.arrange(p_all, p_stl, ncol = 2)  # ncol = 2 表示水平排列
# 
# 





temp_data_meta_stl %>% 
  dplyr::rename("PFAS"="STL") %>% 
  rbind(temp_data_all) -> temp_data_pfas_meta_stl



temp_data_pfas_meta_stl %>%
  filter(class=="STL") %>% 
  ggplot(aes(PC,PFAS)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  theme_bw() +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  scale_color_manual(values = omics_color) +
  labs(x = "PC1", y = "STL" ) +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.5, linetype = "dashed"),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.5, linetype = "dashed") # 次网格线粗细,
        # panel.border = element_rect(size = 1, colour = "black"), # 设置边框线的粗细和颜色
        # plot.margin = unit(c(1, 1, 1, 1), "cm")  # 设置整个图形的外边距
  )+
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 1)


# 
# 
# omics_color <-
#   c(
#     "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
#     "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
#     "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
#     "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
#     "PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
#     "6:2 Cl-PFESA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6],
#    
#   )


omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "6:2 Cl-PFESA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6],
    "PFAS"=RColorBrewer::brewer.pal(n = 12, name = "Set3")[7],
    "STL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[8]
  )



# ,
# gut_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[7],
# skin_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[8],
# oral_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[10],
# nasal_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[11]

temp_data_pfas_meta_stl$class



temp_data_pfas_meta_stl %>%
  mutate(class=case_when(
    class=="PFOA"~"PFOA",
    class=="PFOS"~"PFOS",
    class=="PFDA"~"PFDA",
    class=="PFUdA"~"PFUdA",
    class=="L_PFHxS"~"PFHxS",
    class=="PF3ONS_9CL"~"6:2 Cl-PFESA",
    class=="PFAS"~"PFAS",
    class=="STL"~"STL"
  )) %>% 
  ggplot(aes(PFAS,PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  scale_color_manual(values = omics_color) +
  theme_bw() +
  labs(x = "Spearman corrlation", y = "PC1")+
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.5, linetype = "dashed"),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.5, linetype = "dashed") # 次网格线粗细,
        # panel.border = element_rect(size = 1, colour = "black"), # 设置边框线的粗细和颜色
        # plot.margin = unit(c(1, 1, 1, 1), "cm")  # 设置整个图形的外边距
  ) +
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 2) -> p_pfas_meta_stl

p_pfas_meta_stl


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

###导出保存为pdf
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_pfas_meta_stl_cor.pdf",width =16, height =8, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p_pfas_meta_stl
# par(opar)
dev.off()





temp_data_pfas_meta_stl %>%
  group_by(class) %>%
  dplyr::summarise(cor = cor(PFAS,PC, method = "spearman"))





####绝对值回归系数与p值
library(plyr) 

temp_data_pfas_meta_stl %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PFAS, x$PC, method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()


####>>>相反方向的关联，主要做第一组分对端粒关联分析


temp_data_pfas_meta_stl %>%
  mutate(class=case_when(
    class=="PFOA"~"PFOA",
    class=="PFOS"~"PFOS",
    class=="PFDA"~"PFDA",
    class=="PFUdA"~"PFUdA",
    class=="L_PFHxS"~"PFHxS",
    class=="PF3ONS_9CL"~"6:2 Cl-PFESA",
    class=="PFAS"~"PFAS",
    class=="STL"~"STL"
  )) %>% 
  ggplot(aes(PC,PFAS)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  scale_color_manual(values = omics_color) +
  theme_bw() +
  labs(x = "Spearman corrlation", y = "PC1")+
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.5, linetype = "dashed"),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.5, linetype = "dashed") # 次网格线粗细,
        # panel.border = element_rect(size = 1, colour = "black"), # 设置边框线的粗细和颜色
        # plot.margin = unit(c(1, 1, 1, 1), "cm")  # 设置整个图形的外边距
  ) +
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 2) -> p_pfas_meta_stl

p_pfas_meta_stl


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

###导出保存为pdf
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_pfas_meta_stl_cor_reverse.pdf",width =16, height =8, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p_pfas_meta_stl
# par(opar)
dev.off()





temp_data_pfas_meta_stl %>%
  group_by(class) %>%
  dplyr::summarise(cor = cor(PC,PFAS, method = "spearman"))




####绝对值回归系数与p值
library(plyr) 

temp_data_pfas_meta_stl %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PC, x$PFAS,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()




###########################二、读入所有PFAS和精浆代谢数据绘制混合暴露的火山图###############################


read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PFOA_meta.xlsx") %>% 
  cbind(rep("PFOA",271) %>% as_tibble() %>% dplyr::rename("PFAS"="value")) -> temp_data_cor_PFOA_meta

read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PFOS_meta.xlsx") %>% 
  # dplyr::rename("PFAS"="PFOS")  %>% 
  cbind(rep("PFOS",271) %>% as_tibble() %>% dplyr::rename("PFAS"="value"))-> temp_data_cor_PFOS_meta


read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PFDA_meta.xlsx")%>% 
  # dplyr::rename("PFAS"="PFDA") %>% 
  cbind(rep("PFDA",271) %>% as_tibble() %>% dplyr::rename("PFAS"="value")) -> temp_data_cor_PFDA_meta


read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PFUdA_meta.xlsx")%>% 
  # dplyr::rename("PFAS"="PFUdA")%>%
  cbind(rep("PFUdA",271) %>% as_tibble() %>% dplyr::rename("PFAS"="value")) -> temp_data_cor_PFUdA_meta



read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_L_PFHxS_meta.xlsx")%>% 
  # dplyr::rename("PFAS"="L_PFHxS")  %>% 
  cbind(rep("L_PFHxS",271) %>% as_tibble() %>% dplyr::rename("PFAS"="value")) -> temp_data_cor_L_PFHxS_meta


read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PF3ONS_9CL_meta.xlsx")%>% 
  # dplyr::rename("PFAS"="PF3ONS_9CL") %>% 
  cbind(rep("PF3ONS_9CL",271) %>% as_tibble() %>% dplyr::rename("PFAS"="value"))  -> temp_data_cor_PF3ONS_9CL_meta



read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_mix_pfas_meta.xlsx") %>% 
  # dplyr::rename("PFAS"="PF3ONS_9CL") %>% 
  cbind(rep("PFAS",271) %>% as_tibble() %>% dplyr::rename("PFAS"="value"))  -> temp_data_cor_PFAS_meta




read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/cor_meta_stl.xlsx")%>% 
  # dplyr::rename("PFAS"="PF3ONS_9CL") %>% 
  cbind(rep("STL",271) %>% as_tibble() %>% dplyr::rename("PFAS"="value"))  -> temp_data_cor_meta_stl

###将代谢组学对精液质量影响的所有第一组成分进行回归分析

rbind(temp_data_cor_PFOA_meta,
      temp_data_cor_PFOS_meta,
      temp_data_cor_PFDA_meta,
      temp_data_cor_PFUdA_meta,
      temp_data_cor_L_PFHxS_meta,
      temp_data_cor_PF3ONS_9CL_meta,
      temp_data_cor_PFAS_meta,
      temp_data_cor_meta_stl) ->temp_data_all_cor 

temp_data_all_cor$PFAS

RColorBrewer::brewer.pal(n = 12, name = "Set3")

omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "6:2 Cl-PFESA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6],
    "PFAS"=RColorBrewer::brewer.pal(n = 12, name = "Set3")[7],
    "STL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[8]
  )


# ,
# 
# skin_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[8],
# oral_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[10],
# nasal_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[11]



library(dplyr)
library(ggplot2)

names(temp_data_all_cor)

temp_data_all_cor %>%
  mutate(PFAS=case_when(
    PFAS=="PFOA"~"PFOA",
    PFAS=="PFOS"~"PFOS",
    PFAS=="PFDA"~"PFDA",
    PFAS=="PFUdA"~"PFUdA",
    PFAS=="L_PFHxS"~"PFHxS",
    PFAS=="PF3ONS_9CL"~"6:2 Cl-PFESA",
    PFAS== "PFAS"~"PFAS",
    PFAS=="STL"~"STL"
  )) %>% 
  mutate(
    marker = case_when(
      cor_p < 0.05 & spearman_cor > 0 ~ "Up",
      cor_p < 0.05 & spearman_cor < 0 ~ "Down",
      TRUE ~ "No"),
    marker_pfas_combo = paste(marker, PFAS, sep = " ")  # Combine marker and sperm_names2
  ) %>%
  ggplot(aes(spearman_cor, -log(cor_p, 10))) +
  geom_point(aes(
    size = -log(cor_p, 10),
    color = marker_pfas_combo  # Use the combined variable for color
  ), alpha = 0.7) +
  scale_color_manual(values = c(
    "Up PFOS" = "#8DD3C7",
    "Up PFOA" = "#FFED6F",
    "Up PFDA" = "#BEBADA",
    "Up PFUdA" = "#FB8072",
    "Up PFHxS" = "#80B1D3",
    "Up 6:2 Cl-PFESA" = "#FDB462",
    "Up PFAS" = "#B3DE69",
    "Up STL" = "#FCCDE5",
    "Down PFOS" = "#8DD3C7",
    "Down PFOA" = "#FFED6F",
    "Down PFDA" = "#BEBADA",
    "Down PFUdA" = "#FB8072",
    "Down PFHxS" = "#80B1D3",
    "Down 6:2 Cl-PFESA" = "#FDB462",
    "Down PFAS" = "#B3DE69",
    "Down STL" = "#FCCDE5",
    "No PFOS" = "#3D3B25FF",
    "No PFOA" = "#3D3B25FF",
    "No PFDA" = "#3D3B25FF",
    "No PFUdA" = "#3D3B25FF",
    "No PFHxS" = "#3D3B25FF",
    "No 6:2 Cl-PFESA" = "#3D3B25FF",
    "No PFAS" = "#3D3B25FF",
    "No STL" = "#3D3B25FF"
  )) +
  geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2,size = 1) +
  labs(x = "Spearman correlation", y = "-log10(P value)")+
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
        panel.grid.minor = element_line(size = 0.8, linetype = "dashed") # 次网格线粗细,
        # panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  ) ->volcano_plot_2   

# scale_size_continuous(range = c(3, 10))  # Adjust point size range

# c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF")





volcano_plot_2

library(showtext)
library(sysfonts)


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/volcano_correlation_meta.pdf",width = 16, height = 7, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_2

par(opar)
dev.off()



##############################画每个变量的分类火山图###########################
###############################################################################
###############################################################################


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


str(df_unique_3)


names(temp_data_all_cor)

temp_data_all_cor %>% 
  left_join(df_unique_3,by=c("variable_id"="metabolite")) -> temp_data_all_cor_2








temp_data_all_cor_2 %>%
  mutate(PFAS=case_when(
    PFAS=="PFOA"~"PFOA",
    PFAS=="PFOS"~"PFOS",
    PFAS=="PFDA"~"PFDA",
    PFAS=="PFUdA"~"PFUdA",
    PFAS=="L_PFHxS"~"PFHxS",
    PFAS=="PF3ONS_9CL"~"6:2 Cl-PFESA",
    PFAS== "PFAS"~"PFAS",
    PFAS=="STL"~"STL"
  )) %>% 
  mutate(
    marker = case_when(
      cor_p < 0.05 & spearman_cor > 0 ~ "Up",
      cor_p < 0.05 & spearman_cor < 0 ~ "Down",
      TRUE ~ "No"),
    marker_pfas_combo = paste(marker, PFAS, sep = " ")  # Combine marker and sperm_names2
  ) %>%
  ggplot(aes(spearman_cor, -log(cor_p, 10))) +
  geom_point(aes(
    size = -log(cor_p, 10),
    color = marker_pfas_combo  # Use the combined variable for color
  ), alpha = 1) +
  scale_color_manual(values = c(
    "Up PFOS" = "#BC3C29FF",
    "Up PFOA" = "#BC3C29FF",
    "Up PFDA" = "#BC3C29FF",
    "Up PFUdA" = "#BC3C29FF",
    "Up PFHxS" = "#BC3C29FF",
    "Up 6:2 Cl-PFESA" = "#BC3C29FF",
    "Up PFAS" = "#BC3C29FF",
    "Up STL" = "#BC3C29FF",
    "Down PFOS" = "#0072B5FF",
    "Down PFOA" = "#0072B5FF",
    "Down PFDA" = "#0072B5FF",
    "Down PFUdA" = "#0072B5FF",
    "Down PFHxS" = "#0072B5FF",
    "Down 6:2 Cl-PFESA" = "#0072B5FF",
    "Down PFAS" = "#0072B5FF",
    "Down STL" = "#0072B5FF",
    "No PFOS" = "#3D3B25FF",
    "No PFOA" = "#3D3B25FF",
    "No PFDA" = "#3D3B25FF",
    "No PFUdA" = "#3D3B25FF",
    "No PFHxS" = "#3D3B25FF",
    "No 6:2 Cl-PFESA" = "#3D3B25FF",
    "No PFAS" = "#3D3B25FF",
    "No STL" = "#3D3B25FF"
  ), guide = "none") +
  facet_wrap(facets = vars(PFAS),
             scales = "free",
             nrow = 2)+
  geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2,size = 1) +
  labs(x = "Spearman correlation", y = "-log10(P value)")+
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
        panel.grid.major = element_line(size = 0.5, linetype = "dashed"),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.5, linetype = "dashed") # 次网格线粗细,
        # panel.border = element_rect(size = 1, colour = "black"), # 设置边框线的粗细和颜色
        # plot.margin = unit(c(1, 1, 1, 1), "cm")  # 设置整个图形的外边距
  ) +
  geom_text(aes(label = name), 
            hjust = 0.5, 
            vjust = -1.5, 
            size = 4,  # Set font size to 12
            family = "Arial",  # Set font family to Arial
            check_overlap = TRUE, 
            show.legend = FALSE,
            position = position_jitter(width = 0.05, height = 0.05)) -> volcano_plot_3


volcano_plot_3

# 
# 
# library(dplyr)
# library(ggplot2)
# library(ggrepel)
# 
# df_plot <- temp_data_all_cor_2 %>%
#   mutate(PFAS=case_when(
#     PFAS=="PFOA"~"PFOA",
#     PFAS=="PFOS"~"PFOS",
#     PFAS=="PFDA"~"PFDA",
#     PFAS=="PFUdA"~"PFUdA",
#     PFAS=="L_PFHxS"~"PFHxS",
#     PFAS=="PF3ONS_9CL"~"6:2 Cl-PFESA",
#     PFAS== "PFAS"~"PFAS",
#     PFAS=="STL"~"STL"
#   )) %>% 
#   mutate(
#     marker = case_when(
#       cor_p < 0.05 & spearman_cor > 0 ~ "Up",
#       cor_p < 0.05 & spearman_cor < 0 ~ "Down",
#       TRUE ~ "No"),
#     marker_pfas_combo = paste(marker, PFAS, sep = " ")
#   )
# 
# df_plot
# 
# 
# # 每个 PFAS 选出 p 值最小的前 5 个点
# df_labels <- df_plot %>%
#   group_by(PFAS) %>%
#   slice_min(order_by = cor_p, n = 3, with_ties = FALSE) %>%
#   ungroup()
# 
# ggplot(df_plot, aes(spearman_cor, -log10(cor_p))) +
#   geom_point(aes(
#     size = -log10(cor_p),
#     color = marker_pfas_combo
#   ), alpha = 1) +
#   scale_color_manual(values = c(
#     "Up PFOS" = "#BC3C29FF",
#     "Up PFOA" = "#BC3C29FF",
#     "Up PFDA" = "#BC3C29FF",
#     "Up PFUdA" = "#BC3C29FF",
#     "Up PFHxS" = "#BC3C29FF",
#     "Up 6:2 Cl-PFESA" = "#BC3C29FF",
#     "Up PFAS" = "#BC3C29FF",
#     "Up STL" = "#BC3C29FF",
#     "Down PFOS" = "#0072B5FF",
#     "Down PFOA" = "#0072B5FF",
#     "Down PFDA" = "#0072B5FF",
#     "Down PFUdA" = "#0072B5FF",
#     "Down PFHxS" = "#0072B5FF",
#     "Down 6:2 Cl-PFESA" = "#0072B5FF",
#     "Down PFAS" = "#0072B5FF",
#     "Down STL" = "#0072B5FF",
#     "No PFOS" = "#3D3B25FF",
#     "No PFOA" = "#3D3B25FF",
#     "No PFDA" = "#3D3B25FF",
#     "No PFUdA" = "#3D3B25FF",
#     "No PFHxS" = "#3D3B25FF",
#     "No 6:2 Cl-PFESA" = "#3D3B25FF",
#     "No PFAS" = "#3D3B25FF",
#     "No STL" = "#3D3B25FF"
#   ), guide = "none") +
#   facet_wrap(vars(PFAS), scales = "free", nrow = 2) +
#   geom_vline(xintercept = 0, linetype = 2, size=1, color="#e41a1c") +
#   geom_hline(yintercept = -log10(0.05), linetype = 2, size = 1) +
#   labs(x = "Spearman correlation", y = "-log10(P value)") +
#   theme_bw() +
#   theme(
#     legend.text=element_text(size =14, colour = "black", family = "Arial", face = "bold"),
#     text = element_text(size =14, colour = "black", family = "Arial", face = "bold"),
#     legend.title = element_text(size =14, colour = "black", family = "Arial", face = "bold"),
#     axis.text=element_text(size =14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     panel.grid.major = element_line(size = 0.5, linetype = "dashed"),
#     panel.grid.minor = element_line(size = 0.5, linetype = "dashed")
#   ) +
# geom_text_repel(
#   data = df_labels,
#   aes(x = spearman_cor, y = -log10(cor_p), label = name),  # 显式指定x和y
#   inherit.aes = FALSE,   # 避免继承全局映射
#   size = 3,
#   family = "Arial",
#   box.padding = -0.2,
#   point.padding = 0.5,
#   max.overlaps = Inf,
#   segment.color = "black",
#   segment.size = 0.3,
#   force = 8,
#   show.legend = FALSE
# )















library(ggsci)
library(scales)
RColorBrewer::brewer.pal(n = 12, name = "Set3")

show_col(pal_nejm("default")(8))

volcano_plot_3




library(showtext)
library(sysfonts)





# font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/volcano_correlation_meta_P_value.pdf",width = 16, height = 8, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

volcano_plot_3

par(opar)
dev.off()



