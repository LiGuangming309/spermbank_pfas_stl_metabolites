# no_source() 
# 

no_source()
rm(list = ls())
# setwd(r4projects::get_project_wd())

source("R/100-tools.R")

library(tidyverse)
library(tidymass)

getwd()

pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts)
# get_project_wd()
# masstools::setwd_project()

getwd() 



###====================读入一般资料和THM的信息===============================

####基线样本
###----------读入插补后一般信息和DNA、STL合并的数据-----------
df_background_DNA_LT <- read_excel("result/20240407_PFOA_sperm/df_background_DNA_noNA_final.xlsx")
str(df_background_DNA_LT)

df_background_DNA_LT %>% 
  dplyr::select(number:showinter_class,mtDNAcn,STL,abstinence,season) ->df_background_DNA_LT 

###将禁欲时间划分为几类
str(df_background_DNA_LT)

df_background_DNA_LT %>%   ###要利用检测mtDNA和TL的禁欲时间
  mutate(
    abstinence_class=case_when(
      abstinence>2 & abstinence<7 ~1,
      abstinence<=2~2,
      abstinence>=7~3
    )
  )  -> df_background_DNA_LT


df_background_DNA_LT %>% 
  dplyr::count(abstinence_class)





median(df_background_DNA_LT$age)  ###中位年龄是27

mean(df_background_DNA_LT$age)  ###中位年龄是27

df_background_DNA_LT%>%   ###要利用检测mtDNA和TL的禁欲时间
  mutate(
    age_class_2=case_when(
      age<28 ~0,
      age>=28~1   ###将年龄以平均年龄划分为大于等于28岁，和小于28岁
      
    )
  )  -> df_background_DNA_LT



df_background_DNA_LT %>% 
  dplyr::count(BMI_class)

df_background_DNA_LT %>% 
  dplyr::count(smk0)

433+86+645

86+645

731/1164  ###现在从不吸烟




df_background_DNA_LT %>% 
  dplyr::count(drk)

12+720+298

1030/1164   ###现在不饮酒


df_background_DNA_LT %>% 
  dplyr::count( age_class_2)

# 
df_background_DNA_LT%>%   ###要利用检测mtDNA和TL的禁欲时间
  mutate(
    smk_class_2=case_when(
      smk0==3~0,  ###从不吸烟的为0
      smk0==1|smk0==2~1  ###现在和之前吸烟的为1
    ),
    drk_class_2=case_when(
      drk==4~1,  ###绝不饮酒
      drk==1|drk==2|drk==3~2  ###现在或之前饮酒
      
    ),
    marriage_class_2=case_when(
      marriage==1~1,  ###未婚
      marriage==2|marriage==4~2),  ###已婚
    
    abstinence_class_2=case_when(
      abstinence<7 ~1,    ###小于7天禁欲时间
      abstinence>=7~2    ###≥ 7天的禁欲时间
    ),
    BMI_class_2=case_when(BMI_class==1|BMI_class==2~1, 
                          BMI_class==3~2 )
  )  -> df_background_DNA_LT

df_background_DNA_LT %>% 
  dplyr::count(marriage)

df_background_DNA_LT %>% 
  dplyr::count(marriage_class_2)

df_background_DNA_LT %>% 
  dplyr::count(season)


df_background_DNA_LT$marriage <- factor(df_background_DNA_LT$marriage,levels = c(1,2,4),labels = c("Unmarried","Married","Divorced"))

df_background_DNA_LT$marriage_class_2 <- factor(df_background_DNA_LT$marriage_class_2,levels = c(1,2),labels = c("Unmarried","Married"))

df_background_DNA_LT$age_class <- factor(df_background_DNA_LT$age_class,levels = c(1,2,3),labels = c("<25","25-30","30"))

df_background_DNA_LT$age_class_2 <- factor(df_background_DNA_LT$age_class_2,levels = c(0,1),labels = c("low","High"))

df_background_DNA_LT$BMI_class <- factor(df_background_DNA_LT$BMI_class,levels = c(1,2,3),labels = c("18.5-24","<18.5","≥24"))

df_background_DNA_LT$BMI_class_2 <- factor(df_background_DNA_LT$BMI_class_2,
                                           levels = c(1,2),labels = c("<24","≥24"))

df_background_DNA_LT$education <- factor(df_background_DNA_LT$education,levels = c(1,2,3),labels = c("High school","College","Undergraduate and above"))

df_background_DNA_LT$childn_class <- factor(df_background_DNA_LT$childn_class,levels = c(0,1),labels = c("No","Yes"))

df_background_DNA_LT$abstinence_class <- factor(df_background_DNA_LT$abstinence_class,levels = c(1,2,3),labels = c("2.1-7","≤2","≥7"))

df_background_DNA_LT$abstinence_class_2 <- factor(df_background_DNA_LT$abstinence_class_2,levels = c(1,2),labels = c("<7","≥7"))

df_background_DNA_LT$income <- factor(df_background_DNA_LT$income,levels = c(1,2,3),labels = c("<4000","4000-8000",">8000"))

df_background_DNA_LT$smk0 <- factor(df_background_DNA_LT$smk0,levels = c(1,2,3),labels = c("Current","Former","Never"))

df_background_DNA_LT$smk_class_2 <- factor(df_background_DNA_LT$smk_class_2,levels = c(0,1),labels = c("Never","Former/Current"))

df_background_DNA_LT$drk <- factor(df_background_DNA_LT$drk,levels = c(1,2,3,4),labels = c("Current","Former","Occasional","Never"))

df_background_DNA_LT$drk_class_2 <- factor(df_background_DNA_LT$drk_class_2,levels = c(1,2),labels = c("Never","Former/Current"))

df_background_DNA_LT$season<- factor(df_background_DNA_LT$season,levels = c(0,1,2,3),labels = c("Spring","Summer","Autumn","Winter"))


df_background_DNA_LT$waterv_class <- factor(df_background_DNA_LT$waterv_class,levels = c(0,1,2),labels = c("0L","1-1000",">1000"))

df_background_DNA_LT$showinter_class <- factor(df_background_DNA_LT$showinter_class,levels = c(0,1,2),labels = c("<12","12.1-24",">24"))

str(df_background_DNA_LT)


df_background_DNA_LT %>% 
  dplyr::count(marriage_class_2)


df_background_DNA_LT %>% 
  dplyr::count(BMI_class_2)


df_background_DNA_LT %>% 
  dplyr::count(income)

names(df_background_DNA_LT)

df_background_DNA_LT %>% 
  dplyr::select(number,age,BMI,age_class_2,BMI_class_2,education,marriage_class_2,income,smk_class_2,drk_class_2,season,abstinence_class_2,mtDNAcn,STL) -> df_background_DNA_LT_tidy
###将基线中一般资料、精液质量参数、DNA和血压参数挑选出来

str(df_background_DNA_LT_tidy)


# write.xlsx(df_background_DNA_LT_tidy,"result/20240407_PFOA_sperm/df_background_DNA_LT_tidy.xlsx")

str(df_background_DNA_LT_tidy)



###对精液质量参数取对数
df_background_DNA_LT_tidy %>% 
  dplyr::mutate(lg_DNA=log(mtDNAcn),
                lg_TL=log(STL)) ->df_background_DNA_LT__tidy_lg 

# write.xlsx(df_background_DNA_LT_all_Hyper_tidy_lg,"result/20240106_Hyperten_DNA_sperm/df_background_DNA_LT_all_Hyper_tidy_lg.xlsx")

# write.xlsx(df_background_DNA_LT_all_Hyper_tidy_lg,"result/20240106_Hyperten_DNA_sperm/df_background_DNA_LT_all_Hyper_tidy_lg.xlsx")

str(df_background_DNA_LT__tidy_lg)



# 
# ###将BMI转换为24二分类
# df_background_DNA_LT__tidy_lg %>% 
#   mutate(BMI_class_2=case_when(
#     BMI_class=="18.5-24"|BMI_class=="<18.5"~1,
#     BMI_class=="≥24"~2
#   )) -> df_background_DNA_LT__tidy_lg_2
# 


df_background_DNA_LT__tidy_lg %>% 
  dplyr::count(BMI_class_2)

# df_background_DNA_LT__tidy_lg_2$BMI_class_2 <- factor(df_background_DNA_LT__tidy_lg_2$BMI_class_2)


str(df_background_DNA_LT__tidy_lg)


###>>>>>>>>>>>>>>>>>>>>>>>>>一、画第一层
library(circlize)


df <-
  data.frame(
    factors = df_background_DNA_LT__tidy_lg$number,
    x = 1,
    y = 1,
    df_background_DNA_LT__tidy_lg,
    stringsAsFactors = TRUE
  ) %>%
  dplyr::arrange(age) %>%
  dplyr::mutate(factors = factor(factors, levels = factors))

# 设置全局字体为 Arial
par(family = "Arial")

circos.par(
  "track.height" = 0.2,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df) - 1), 90),
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = df$factors,
                  x = df$x,
                  xlim = c(0.8,3))
# Randomize label positions to avoid overlap
# random_angles <- runif(length(df$factors), min = 0, max = 2 * pi)  # Random angle for each sector

##1.age
range(df$age, na.rm = TRUE)
temp_value <- df$age

circos.track(
  factors = df$factors,
  # x = df$x,
  y = temp_value,
  ylim = c(0.5 * min(temp_value), 1.1 * max(temp_value, na.rm = TRUE)),
  bg.border = "white",
  # bg.col = NA,
  track.height = 0.2,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = c(0.8 * min(temp_value),
             round((
               min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)
             ) / 2, 2),
             round(max(
               temp_value, na.rm = TRUE
             ), 2)),
      sector.index = get.all.sector.index()[1],
      labels.cex = 2,labels.font =2,tick = 2,tick.length = 6,lwd = 2,
      labels.niceFacing = FALSE
    )
    
    circos.lines(
      x = mean(xlim, na.rm = TRUE),
      y =  temp_value[i],
      pch = 16,
      cex = 0.8,
      type = "h",
      col = ggsci::pal_aaas()(n = 10)[4],
      lwd = 0.8
    )
    # Randomized label placement to avoid overlap
    # random_angle <- random_angles[i]  # Get the random angle for this sector
    #plot country labels
    circos.text(
      x = 1,
      y = 50,
      labels = name,
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.3,family="Arial"
      # adj = aa
    )
    
    circos.points(
      x = mean(xlim),
      y =  temp_value[i],
      pch = 16,
      cex = 0.2,
      col = ggsci::pal_aaas()(n = 10)[4]
    )
  }
)

library(scales)

library(ggsci)
# 
# show_col(ggsci::pal_aaas()(n = 10))
# 


##2.BMI
range(df$BMI, na.rm = TRUE)
temp_value <- df$BMI

circos.track(
  factors = df$factors,
  # x = df$x,
  y = temp_value,
  ylim = c(
    0.5 * min(temp_value, na.rm = TRUE),
    1.1 * max(temp_value, na.rm = TRUE)
  ),
  bg.border = "white",
  # bg.col = NA,
  track.height = 0.2,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = c(
        round(0.8 * min(temp_value, na.rm = TRUE),2),
        round((
          min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)
        ) / 2, 2),
        round(max(temp_value, na.rm = TRUE), 2)
      ),
      sector.index = get.all.sector.index()[1],
      labels.cex = 2,labels.font =2,tick = 2,tick.length = 6,lwd = 2,
      labels.niceFacing = FALSE
    )
    
    circos.lines(
      x = mean(xlim, na.rm = TRUE),
      y =  temp_value[i],
      pch = 16,
      cex = 0.8,
      type = "h",
      col = ggsci::pal_tron()(n = 10)[1],
      lwd = 0.8
    )
    
    circos.points(
      x = mean(xlim),
      y =  temp_value[i],
      pch = 16,
      cex = 0.2,
      col = ggsci::pal_tron()(n = 10)[1]
    )
  }
)



# install.packages("wesanderson")
library(wesanderson)
###3.education
education_color <-
  c(
    "High school" = wesanderson::wes_palettes$Rushmore1[5],
    "College"=wesanderson::wes_palettes$Rushmore1[3],
    "Undergraduate and above" = wesanderson::wes_palettes$Rushmore1[1]
  )

# df$age_class_2
## sex
temp_education <- as.character(df$education)   ###去除数据框格式

names(df)

as.vector(df$age_class_2)
temp_education[is.na(temp_education)] <- "grey"
temp_education[temp_education =="High school"] <- education_color["High school"]
temp_education[temp_education =="College"] <- education_color["College"]
temp_education[temp_education =="Undergraduate and above"] <- education_color["Undergraduate and above"]
temp_education
# age_color[temp_age]
# 
# temp_age %>% 
#   mutate(value=case_when(
#     temp_age=="low"~age_color["low"],
#     temp_age=="High"~age_color["High"]
#   ))



circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "white",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_education[i],
      bg.border = "white"
    )
  }
)


###4.income

income_color <-
  c(
    "<4000" =wesanderson::wes_palettes$Darjeeling2[1],
    "4000-8000"= wesanderson::wes_palettes$Darjeeling2[2],
    ">8000" = wesanderson::wes_palettes$Darjeeling2[3]
  )

df$income
## sex
temp_income <- as.character(df$income)   ###去除数据框格式

names(df)

as.vector(df$income)
temp_income[is.na(temp_income)] <- "grey"
temp_income[temp_income =="<4000"] <- income_color["<4000"]
temp_income[temp_income =="4000-8000"] <- income_color["4000-8000"]
temp_income[temp_income ==">8000"] <- income_color[">8000"]
temp_income
# age_color[temp_age]
# 
# temp_age %>% 
#   mutate(value=case_when(
#     temp_age=="low"~age_color["low"],
#     temp_age=="High"~age_color["High"]
#   ))



circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "white",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_income[i],
      bg.border = "white"
    )
  }
)



###5.abstinence

abstinence_color <-
  c(
    "<7" =RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "≥7" = wesanderson::wes_palettes$Zissou1[1]
  )

## sex
temp_abstinence <- as.character(df$abstinence_class_2)   ###去除数据框格式

names(df)

# as.vector(df$income)
temp_abstinence[is.na(temp_abstinence)] <- "grey"
temp_abstinence[temp_abstinence =="<7"] <- abstinence_color["<7"]
temp_abstinence[temp_abstinence =="≥7"] <- abstinence_color["≥7"]

temp_abstinence
# age_color[temp_age]
# 
# temp_age %>% 
#   mutate(value=case_when(
#     temp_age=="low"~age_color["low"],
#     temp_age=="High"~age_color["High"]
#   ))



circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "white",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_abstinence[i],
      bg.border = "white"
    )
  }
)

# # 
# ###6.smoke
# 
# # show_col(RColorBrewer::brewer.pal(n = 12, name = "Set3"))
# smk_color <-
#   c(
#     "Never" =RColorBrewer::brewer.pal(n = 12, name = "Set3")[7],
#     "Former/Current" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[8]
#   )
# 
# df$smk_class_2
# ## sex
# temp_smk <- as.character(df$smk_class_2)   ###去除数据框格式
# 
# names(df)
# 
# # as.vector(df$income)
# temp_smk[is.na(temp_smk)] <- "grey"
# temp_smk[temp_smk =="Never"] <- smk_color["Never"]
# temp_smk[temp_smk =="Former/Current"] <- smk_color["Former/Current"]
# 
# temp_abstinence
# # age_color[temp_age]
# #
# # temp_age %>%
# #   mutate(value=case_when(
# #     temp_age=="low"~age_color["low"],
# #     temp_age=="High"~age_color["High"]
# #   ))
# 
# 
# 
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = df$y,
#   ylim = c(0, 1),
#   bg.border = "white",
#   # bg.col = NA,
#   track.height = 0.1,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
# 
#     #text direction (dd) and adjusmtents (aa)
#     theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
#     dd <-
#       ifelse(theta < 90 ||
#                theta > 270, "clockwise", "reverse.clockwise")
#     aa = c(0.5, 1)
#     # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
# 
#     circos.rect(
#       xleft = xlim[1],
#       ybottom = ylim[1],
#       xright = xlim[2],
#       ytop = ylim[2],
#       col = temp_smk[i],
#       bg.border = "white"
#     )
#   }
# )
# 
# 
# ###7.drk
# 
# # show_col(RColorBrewer::brewer.pal(n = 12, name = "Set3"))
# drk_color <-
#   c(
#     "Never" =RColorBrewer::brewer.pal(n = 12, name = "Set3")[11],
#     "Former/Current" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12]
#   )
# 
# df$drk_class_2
# ## sex
# temp_drk <- as.character(df$drk_class_2)   ###去除数据框格式
# 
# names(df)
# 
# # as.vector(df$income)
# temp_drk[is.na(temp_smk)] <- "grey"
# temp_drk[temp_drk =="Never"] <- smk_color["Never"]
# temp_drk[temp_drk =="Former/Current"] <- smk_color["Former/Current"]
# 
# temp_drk
# # age_color[temp_age]
# #
# # temp_age %>%
# #   mutate(value=case_when(
# #     temp_age=="low"~age_color["low"],
# #     temp_age=="High"~age_color["High"]
# #   ))
# 
# 
# 
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = df$y,
#   ylim = c(0, 1),
#   bg.border = "white",
#   # bg.col = NA,
#   track.height = 0.1,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
# 
#     #text direction (dd) and adjusmtents (aa)
#     theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
#     dd <-
#       ifelse(theta < 90 ||
#                theta > 270, "clockwise", "reverse.clockwise")
#     aa = c(0.5, 1)
#     # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
# 
#     circos.rect(
#       xleft = xlim[1],
#       ybottom = ylim[1],
#       xright = xlim[2],
#       ytop = ylim[2],
#       col = temp_drk[i],
#       bg.border = "white"
#     )
#   }
# )
# 



# 
# ggsave(plot_bmi,
#        filename = "plot_bmi.pdf",
#        width = 3,
#        height = 10)






#########


###age
#####age
age <-df$age

IQR(age)

mean(age)
median(age)
round(quantile(age,0.10),digits = 6)
round(quantile(age,0.25),digits = 6)
round(quantile(age,0.50),digits = 6)
round(quantile(age,0.75),digits = 6)
round(quantile(age,0.90),digits = 6)


library(gghalves)
plot_age <-
  age %>%
  data.frame(class = "class", value = .) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(outlier.shape = NA,size = 1, color = "black") +
  geom_dotplot(
    binaxis = "y",
    color = ggsci::pal_aaas()(n = 10)[4],
    fill = ggsci::pal_aaas()(n = 10)[4],
    shape = 16,
    binwidth = 0.15,
    stackdir = "center"
  ) +
 
  theme_bw() +
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black", size = 1.5)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
  )+
  labs(x = "", y = "") 


plot_age



###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/age_describe.pdf",width = 4, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_age

par(opar)
dev.off()






###BMI
#####age
BMI <-df$BMI
IQR(BMI)

mean(BMI)
median(BMI)
round(quantile(BMI,0.10),digits = 6)
round(quantile(BMI,0.25),digits = 6)
round(quantile(BMI,0.50),digits = 6)
round(quantile(BMI,0.75),digits = 6)
round(quantile(BMI,0.90),digits = 6)

library(gghalves)
plot_BMI <-
  BMI %>%
  data.frame(class = "class", value = .) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(outlier.shape = NA,size = 1, color = "black") +
  geom_dotplot(
    binaxis = "y",
    color = ggsci::pal_tron()(n = 10)[1],
    fill = ggsci::pal_tron()(n = 10)[1],
    shape = 16,
    binwidth = 0.15,
    stackdir = "center"
  ) +
  
  theme_bw() +
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black", size = 1.5)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
  )+
  labs(x = "", y = "") 


plot_BMI



###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/BMI_describe.pdf",width =4, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_BMI

par(opar)
dev.off()



###education
education_color <-
  c(
    "High school" = wesanderson::wes_palettes$Rushmore1[5],
    "College"=wesanderson::wes_palettes$Rushmore1[3],
    "Undergraduate and above" = wesanderson::wes_palettes$Rushmore1[1]
  )
education <-
  df$education

plot_education <-
  education %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("High school", "College","Undergraduate and above"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = education_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_education


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/education_describe.pdf",width =2.35, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_education

par(opar)
dev.off()



###education图例
education_color <-
  c(
    "High school" = wesanderson::wes_palettes$Rushmore1[5],
    "College"=wesanderson::wes_palettes$Rushmore1[3],
    "Undergraduate and above" = wesanderson::wes_palettes$Rushmore1[1]
  )
education <-
  df$education

plot_education_legend <-
  education %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("High school", "College","Undergraduate and above"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = TRUE,
    width = 2
  ) +
  scale_fill_manual(values = education_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1),
    legend.title = element_text(size = 14, family = "Arial", face = "bold", color = "black"),  # 图例标题样式
    legend.text = element_text(size = 12, family = "Arial", color = "black"),  # 图例文字样式
    legend.position = "right"  # 图例位置
  )

plot_education_legend 


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/education_describe_legend.pdf",width =2.35, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_education_legend 

par(opar)
dev.off()



###income

show_col(wesanderson::wes_palettes$Darjeeling2)

income_color <-
  c(
    "<4000" =wesanderson::wes_palettes$Darjeeling2[1],
    "4000-8000"= wesanderson::wes_palettes$Darjeeling2[2],
    ">8000" = wesanderson::wes_palettes$Darjeeling2[3]
  )

income <-
  df$income

plot_income <-
  income %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("<4000","4000-8000",">8000"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = income_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_income


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/income_describe.pdf",width =2.2, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_income

par(opar)
dev.off()


###income图列 

plot_income_legend <-
  income %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("<4000","4000-8000",">8000"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = TRUE,
    width = 2
  ) +
  scale_fill_manual(values = income_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1),
    legend.title = element_text(size = 14, family = "Arial", face = "bold", color = "black"),  # 图例标题样式
    legend.text = element_text(size = 12, family = "Arial", color = "black"),  # 图例文字样式
    legend.position = "right"  # 图例位置
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_income_legend


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/income_describe_legend.pdf",width =2.2, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_income_legend

par(opar)
dev.off()






###income

show_col(wesanderson::wes_palettes$Rushmore1)

show_col(RColorBrewer::brewer.pal(n = 12, name = "Set3"))

show_col(wesanderson::wes_palettes$Zissou1)

abstinence_color <-
  c(
    "<7" =RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "≥7" = wesanderson::wes_palettes$Zissou1[1]
  )

abstinence <-
  df$abstinence_class_2

plot_abstinence <-
  abstinence %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("<7","≥7"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = abstinence_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_abstinence 


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/abstinence _describe.pdf",width =2.2, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_abstinence 

par(opar)
dev.off()






###income图例

show_col(wesanderson::wes_palettes$Rushmore1)

show_col(RColorBrewer::brewer.pal(n = 12, name = "Set3"))

show_col(wesanderson::wes_palettes$Zissou1)

abstinence_color <-
  c(
    "<7" =RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "≥7" = wesanderson::wes_palettes$Zissou1[1]
  )

abstinence <-
  df$abstinence_class_2

plot_abstinence_legend <-
  abstinence %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("<7","≥7"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = TRUE,
    width = 2
  ) +
  scale_fill_manual(values = abstinence_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1),
    legend.title = element_text(size = 14, family = "Arial", face = "bold", color = "black"),  # 图例标题样式
    legend.text = element_text(size = 12, family = "Arial", color = "black"),  # 图例文字样式
    legend.position = "right"  # 图例位置
  )

plot_abstinence 


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/abstinence _describe_legend.pdf",width =2.2, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_abstinence_legend

par(opar)
dev.off()



###########>>>>>>>>>>>>>>>>>>>>fig1_小图<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


###>>>>>>>>>>>>>>>>>>>>>>>>>一、画第一层
library(circlize)


df <-
  data.frame(
    factors = df_background_DNA_LT__tidy_lg$number,
    x = 1,
    y = 1,
    df_background_DNA_LT__tidy_lg,
    stringsAsFactors = TRUE
  ) %>%
  dplyr::arrange(age) %>%
  dplyr::mutate(factors = factor(factors, levels = factors))

# 设置全局字体为 Arial
par(family = "Arial")

circos.par(
  "track.height" = 0.2,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df) - 1), 90),
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = df$factors,
                  x = df$x,
                  xlim = c(0.8,2))
# Randomize label positions to avoid overlap
# random_angles <- runif(length(df$factors), min = 0, max = 2 * pi)  # Random angle for each sector


# # 
###6.smoke
library(scales)
library(ggsci)
# show_col(RColorBrewer::brewer.pal(n = 12, name = "Set3"))

# show_col(pal_ucscgb("default", alpha = 1)(26))

smk_color <-
  c(
    "Never" =pal_ucscgb("default", alpha = 1)(26)[19],
    "Former/Current" = pal_ucscgb("default", alpha = 1)(26)[20]
  )

df$smk_class_2
## sex
temp_smk <- as.character(df$smk_class_2)   ###去除数据框格式

names(df)

# as.vector(df$income)
temp_smk[is.na(temp_smk)] <- "grey"
temp_smk[temp_smk =="Never"] <- smk_color["Never"]
temp_smk[temp_smk =="Former/Current"] <- smk_color["Former/Current"]

temp_smk
# age_color[temp_age]
#
# temp_age %>%
#   mutate(value=case_when(
#     temp_age=="low"~age_color["low"],
#     temp_age=="High"~age_color["High"]
#   ))



circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "white",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")

    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)

    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_smk[i],
      bg.border = "white"
    )
  }
)


###7.drk

# show_col(RColorBrewer::brewer.pal(n = 12, name = "Set3"))
# 
# show_col(RColorBrewer::brewer.pal(n = 12, name = "Set3"))

# show_col(pal_ucscgb("default", alpha = 1)(26))
# drk_color <-
#   c(
#     "Never" =pal_ucscgb("default", alpha = 1)(26)[13],
#     "Former/Current" = pal_ucscgb("default", alpha = 1)(26)[10]
#   )

drk_color <-
  c(
    "Never" ="#2CCDF8",
    "Former/Current" = "#FF699C"
  )

df$drk_class_2
## sex
temp_drk <- as.character(df$drk_class_2)   ###去除数据框格式

names(df)

# as.vector(df$income)
temp_drk[is.na(temp_smk)] <- "grey"
temp_drk[temp_drk =="Never"] <- drk_color["Never"]
temp_drk[temp_drk =="Former/Current"] <- drk_color["Former/Current"]

temp_drk
# age_color[temp_age]
#
# temp_age %>%
#   mutate(value=case_when(
#     temp_age=="low"~age_color["low"],
#     temp_age=="High"~age_color["High"]
#   ))



circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "white",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")

    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)

    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_drk[i],
      bg.border = "white"
    )
  }
)







###7.season

# show_col(RColorBrewer::brewer.pal(n = 12, name = "Set3"))
# 
# show_col(RColorBrewer::brewer.pal(n = 12, name = "Set3"))

# show_col(pal_ucscgb("default", alpha = 1)(26))
# drk_color <-
#   c(
#     "Never" =pal_ucscgb("default", alpha = 1)(26)[13],
#     "Former/Current" = pal_ucscgb("default", alpha = 1)(26)[10]
#   )

season_color <-
  c(
    "Spring" = "#4D388F",
    "Summer" = "#EF99FC",
    "Autumn" = "#9185B5",
    "Winter" =  "white"
    
  )

df$drk_class_2
## sex
temp_season <- as.character(df$season)   ###去除数据框格式

names(df)

# as.vector(df$income)
temp_season[is.na(temp_season)] <- "grey"
temp_season[temp_season =="Spring"] <- season_color["Spring"]
temp_season[temp_season =="Summer"] <- season_color["Summer"]
temp_season[temp_season =="Autumn"] <- season_color["Autumn"]
temp_season[temp_season =="Winter"] <- season_color["Winter"]
temp_season
# age_color[temp_age]
#
# temp_age %>%
#   mutate(value=case_when(
#     temp_age=="low"~age_color["low"],
#     temp_age=="High"~age_color["High"]
#   ))



circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "white",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_season[i],
      bg.border = "white"
    )
  }
)






###绘制图例




###smk


# temp_smk[temp_smk =="Never"] <- smk_color["Never"]
# temp_smk[temp_smk =="Former/Current"] <- smk_color["Former/Current"]


smk_color <-
  c(
    "Never" =pal_ucscgb("default", alpha = 1)(26)[19],
    "Former/Current" = pal_ucscgb("default", alpha = 1)(26)[20]
  )

smk <- df$smk_class_2

plot_smk <-
  smk %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("Never","Former/Current"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = smk_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_smk  


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/smk _describe.pdf",width =1.5, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_smk 

par(opar)
dev.off()






###smk图例


smk_color <-
  c(
    "Never" =pal_ucscgb("default", alpha = 1)(26)[19],
    "Former/Current" = pal_ucscgb("default", alpha = 1)(26)[20]
  )

smk <- df$smk_class_2

plot_smk_legend <-
  smk %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("Never","Former/Current"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = TRUE,
    width = 2
  ) +
  scale_fill_manual(values = smk_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_smk_legend 



###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/smk _describe_legend.pdf",width =1.5, height = 9, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_smk_legend

# par(opar)
dev.off()










###drk


# temp_smk[temp_smk =="Never"] <- smk_color["Never"]
# temp_smk[temp_smk =="Former/Current"] <- smk_color["Former/Current"]


drk_color <-
  c(
    "Never" ="#2CCDF8",
    "Former/Current" = "#FF699C"
  )

drk <- df$drk_class_2

plot_drk <-
  drk%>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("Never","Former/Current"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = drk_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_drk 


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/drk _describe.pdf",width =1.5, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_drk 

par(opar)
dev.off()






###smk图例


drk_color <-
  c(
    "Never" ="#2CCDF8",
    "Former/Current" = "#FF699C"
  )

drk <- df$drk_class_2

plot_drk_legend <-
  drk%>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("Never","Former/Current"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = TRUE,
    width = 2
  ) +
  scale_fill_manual(values = drk_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_drk_legend


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/drk _describe_legend.pdf",width =1.5, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_drk_legend 

par(opar)
dev.off()









###season


# temp_smk[temp_smk =="Never"] <- smk_color["Never"]
# temp_smk[temp_smk =="Former/Current"] <- smk_color["Former/Current"]



season_color <-
  c(
    "Spring" = "#4D388F",
    "Summer" = "#EF99FC",
    "Autumn" = "#9185B5",
    "Winter" =  "white"
    
  )

season <- df$season

plot_season <-
  season%>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("Spring","Summer","Autumn","Winter"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = season_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_season 


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/season _describe.pdf",width =1.5, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_season 

par(opar)
dev.off()









season_color <-
  c(
    "Spring" = "#4D388F",
    "Summer" = "#EF99FC",
    "Autumn" = "#9185B5",
    "Winter" =  "white"
    
  )

season <- df$season

plot_season_legend <-
  season%>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("Spring","Summer","Autumn","Winter"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = TRUE,
    width = 2
  ) +
  scale_fill_manual(values = season_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_text(size = 14, family = "Arial", face = "bold", color = "black"),
    # axis.text = element_text(family = "Arial",size = 14,colour = "black"),
    panel.border = element_rect(color = "black",size=1)
    # 修改纵横坐标轴线的粗细
    # axis.line = element_line(color = "black", size = 1),  # 纵横坐标轴线粗细
    # axis.ticks = element_line(color = "black", size = 1)  # 坐标轴刻度线粗细
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

plot_season_legend


###添加竖直的参考线
cairo_pdf( "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/Fifure1/season _describe_legend.pdf",width =1.5, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
plot_season_legend

par(opar)
dev.off()





























# 
# 
# 
# ##========================读取精浆PFOA和血浆PFOA数据=================================
# 
# 
# ###-----------------读入精浆PFOA的信息-----------------------------
# 
# df_sperm_PFOA <- read_excel("raw_data/Sperm_data_share_lgm/血浆-精浆PFAS数据-final-20240524.xlsx",sheet = "精浆PFAS")
# 
# str(df_sperm_PFOA)
# 
# ###删除前2行数据
# df_sperm_PFOA %>% 
#   slice(-c(1:2)) ->df_sperm_PFOA_tidy 
# 
# ###将编码变量重新命名
# df_sperm_PFOA_tidy  %>% 
#   dplyr::rename(number=`Semen PFAS`) ->df_sperm_PFOA_tidy_2
# str(df_sperm_PFOA_tidy_2)
# 
# df_sperm_PFOA_tidy_2$number <- as.numeric(df_sperm_PFOA_tidy_2$number)
# 
# 
# str(df_sperm_PFOA_tidy_2)
# 
# names(df_sperm_PFOA_tidy_2)
# 
# df_sperm_PFOA_tidy_2 %>% 
#   dplyr::rename(L_PFHxS="L-PFHxS",
#                 PF3ONS_9CL="9CL-PF3ONS",
#                 L_PFBS="L-PFBS",
#                 L_PFHpS="L-PFHpS",
#                 PF3OUdS_11CL="11CL-PF3OUdS") -> df_sperm_PFOA_tidy_3
# 
# 
# 
# 
# #### 将PFOA检测限LOD/根号2替代缺失值
# df_sperm_PFOA_tidy_3 %>% 
#   mutate(
#     PFOS=case_when(
#       PFOS<0.0516440006885867~0.0516440006885867/sqrt(2),
#       TRUE~ PFOS),
#     PFOA=case_when(
#       PFOA<0.00675113081441142~0.00675113081441142/sqrt(2),
#       TRUE~ PFOA),
#     PFDA=case_when(
#       PFDA<0.00397719740156436~0.00397719740156436/sqrt(2),
#       TRUE~ PFDA),
#     PFUdA=case_when(
#       PFUdA<0.00351658656663932~0.00351658656663932/sqrt(2),
#       TRUE~ PFUdA),
#     L_PFHxS=case_when(
#       L_PFHxS<0.00338104361546264~0.00338104361546264/sqrt(2),
#       TRUE~ L_PFHxS),
#     PF3ONS_9CL=case_when(
#       PF3ONS_9CL<0.000503330369276714~0.000503330369276714/sqrt(2),
#       TRUE~ PF3ONS_9CL),
#     PFBA=case_when(
#       PFBA<0.0319216854649925~0.0319216854649925/sqrt(2),
#       TRUE~ PFBA),
#     PFNA=case_when(
#       PFNA<0.045045045045045~0.045045045045045/sqrt(2),
#       TRUE~ PFNA),
#     PFDoA=case_when(
#       PFDoA<0.00255493101686255~0.00255493101686255/sqrt(2),
#       TRUE~ PFDoA),
#     L_PFBS=case_when(
#       L_PFBS<0.00494967827091239~0.00494967827091239/sqrt(2),
#       TRUE~ L_PFBS),
#     L_PFHpS=case_when(
#       L_PFHpS<0.00917992656058752~0.00917992656058752/sqrt(2),
#       TRUE~ L_PFHpS),
#     PF3OUdS_11CL=case_when(
#       PF3OUdS_11CL<0.000540365286933967~0.000540365286933967/sqrt(2),
#       TRUE~ PF3OUdS_11CL)
#   ) -> df_sperm_PFOA_tidy_4
# 
# str(df_sperm_PFOA_tidy_4)
# 
# df_background_DNA_LT__tidy_lg %>% 
#   left_join(df_sperm_PFOA_tidy_4,by="number") ->df_background_DNA_LT_sperm_PFOA 
# 
# str(df_background_DNA_LT_sperm_PFOA)
# 
# summary(df_background_DNA_LT_sperm_PFOA)
# colSums(is.na(df_background_DNA_LT_sperm_PFOA)) 
# 
# 
# 
# 
# ####删除缺失值
# 
# df_background_DNA_LT_sperm_PFOA %>% 
#   dplyr::filter(!PFOA=="NA") ->df_background_DNA_LT_sperm_PFOA_noNA  ###得到的精浆和人群基线数据一样的数据
# 
# str(df_background_DNA_LT_sperm_PFOA_noNA)   ###836人
# 
# length(unique(df_background_DNA_LT_sperm_PFOA_noNA$number))
# 
# colSums(is.na(df_background_DNA_LT_sperm_PFOA_noNA)) 
# 
# 
# 
# ####最终数据
# 
# str(df_background_DNA_LT_sperm_PFOA_noNA)
# 
# ###将数据重新排序整理
# df_background_DNA_LT_sperm_PFOA_noNA %>% 
#   select(PFOS:PF3OUdS_11CL,mtDNAcn,STL,everything()) -> df_background_DNA_LT_sperm_PFOA_tidy
# 
# 
# str(df_background_DNA_LT_sperm_PFOA_tidy)
# 
# 
# ####对PFAS取对数
# lapply(df_background_DNA_LT_sperm_PFOA_tidy[,c(1:14)],log) ->df_background_DNA_LT_sperm_PFOA_tidy_lg 
# 
# ### PFAS和DNA均取对数了
# df_background_DNA_LT_sperm_PFOA_tidy %>% 
#   dplyr::select(-c(1:14,25,26)) %>% 
#   cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 
# 
# str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###====================读入代谢组学数据===============================
# 
# df_meta<- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/240823_最终确定版data.xlsx",sheet = "已核对")
# str(df_meta) 
# 
# 
# df_meta %>% 
#   dplyr::rename(number="ID") -> df_meta_1
# 
# str(df_meta_1)
# 
# names(df_meta_1)
# 
# ####>>>>>>>删除外源或二肽物质
# df_out <- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/外源_20241020.xlsx",sheet = "外源和重复")
# 
# df_out %>% 
#   dplyr::rename(metabolite="META")  -> df_out_2
# 
# str(df_out_2)   ###一共有266个二肽或外源物质
# 
# 
# df_out_2$metabolite
# 
# 
# 
# df_meta_1 %>% 
#   select(-c(df_out_2$metabolite))  ->  df_meta_2  ####最终得到414个代谢物质
# 
# str(df_meta_2)
# 
# 794-437
# 
# 437-23
# 
# 
# ###将剩余的物质与PFAS和STL进行合并
# 
# df_background_DNA_LT_sperm_PFOA_tidy_final_lg  %>% 
#   left_join(df_meta_2,by="number") -> df_background_DNA_LT_sperm_PFOA_meta_tidy
# 
# colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy))
# 
# 
# 
# 
# ####删除缺失值
# 
# df_background_DNA_LT_sperm_PFOA_meta_tidy  %>% 
#   dplyr::filter(!sampleID=="NA")  %>% 
#   as_tibble()->df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA  ###得到的精浆和人群基线数据一样的数据
# 
# 
# 
# colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA))
# 
# df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA
# 
# names(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA)
