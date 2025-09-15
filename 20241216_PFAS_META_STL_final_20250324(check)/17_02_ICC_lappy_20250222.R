# no_source() 
# 

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

df_background_DNA_LT %>% 
  dplyr::count(drk)

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

df_background_DNA_LT$age_class_2 <- factor(df_background_DNA_LT$age_class_2,levels = c(0,1),labels = c("<28","≥28"))

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



df_background_DNA_LT %>% 
  dplyr::select(number,age_class_2,BMI_class_2,education,marriage_class_2,income,smk_class_2,drk_class_2,season,abstinence_class_2,mtDNAcn,STL) -> df_background_DNA_LT_tidy
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


# write.xlsx(df_background_DNA_LT__tidy_lg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_repeatsemen_final_1164.xlsx")






##========================读取精浆PFOA和血浆PFOA数据=================================


###-----------------读入精浆PFOA的信息-----------------------------

df_sperm_PFOA <- read_excel("raw_data/Sperm_data_share_lgm/血浆-精浆PFAS数据-final-20240524.xlsx",sheet = "精浆PFAS")

str(df_sperm_PFOA)

###删除前2行数据
df_sperm_PFOA %>% 
  slice(-c(1:2)) ->df_sperm_PFOA_tidy 

###将编码变量重新命名
df_sperm_PFOA_tidy  %>% 
  dplyr::rename(number=`Semen PFAS`) ->df_sperm_PFOA_tidy_2
str(df_sperm_PFOA_tidy_2)

df_sperm_PFOA_tidy_2$number <- as.numeric(df_sperm_PFOA_tidy_2$number)


str(df_sperm_PFOA_tidy_2)

names(df_sperm_PFOA_tidy_2)

df_sperm_PFOA_tidy_2 %>% 
  dplyr::rename(L_PFHxS="L-PFHxS",
                PF3ONS_9CL="9CL-PF3ONS",
                L_PFBS="L-PFBS",
                L_PFHpS="L-PFHpS",
                PF3OUdS_11CL="11CL-PF3OUdS") -> df_sperm_PFOA_tidy_3




#### 将PFOA检测限LOD/根号2替代缺失值
df_sperm_PFOA_tidy_3 %>% 
  mutate(
    PFOS=case_when(
      PFOS<0.0516440006885867~0.0516440006885867/sqrt(2),
      TRUE~ PFOS),
    PFOA=case_when(
      PFOA<0.00675113081441142~0.00675113081441142/sqrt(2),
      TRUE~ PFOA),
    PFDA=case_when(
      PFDA<0.00397719740156436~0.00397719740156436/sqrt(2),
      TRUE~ PFDA),
    PFUdA=case_when(
      PFUdA<0.00351658656663932~0.00351658656663932/sqrt(2),
      TRUE~ PFUdA),
    L_PFHxS=case_when(
      L_PFHxS<0.00338104361546264~0.00338104361546264/sqrt(2),
      TRUE~ L_PFHxS),
    PF3ONS_9CL=case_when(
      PF3ONS_9CL<0.000503330369276714~0.000503330369276714/sqrt(2),
      TRUE~ PF3ONS_9CL),
    PFBA=case_when(
      PFBA<0.0319216854649925~0.0319216854649925/sqrt(2),
      TRUE~ PFBA),
    PFNA=case_when(
      PFNA<0.045045045045045~0.045045045045045/sqrt(2),
      TRUE~ PFNA),
    PFDoA=case_when(
      PFDoA<0.00255493101686255~0.00255493101686255/sqrt(2),
      TRUE~ PFDoA),
    L_PFBS=case_when(
      L_PFBS<0.00494967827091239~0.00494967827091239/sqrt(2),
      TRUE~ L_PFBS),
    L_PFHpS=case_when(
      L_PFHpS<0.00917992656058752~0.00917992656058752/sqrt(2),
      TRUE~ L_PFHpS),
    PF3OUdS_11CL=case_when(
      PF3OUdS_11CL<0.000540365286933967~0.000540365286933967/sqrt(2),
      TRUE~ PF3OUdS_11CL)
  ) -> df_sperm_PFOA_tidy_4

str(df_sperm_PFOA_tidy_4)


write.xlsx(df_sperm_PFOA_tidy_4,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_sperm_PFOA_tidy_4_1053.xlsx")

df_background_DNA_LT__tidy_lg %>% 
  left_join(df_sperm_PFOA_tidy_4,by="number") ->df_background_DNA_LT_sperm_PFOA 

str(df_background_DNA_LT_sperm_PFOA)

summary(df_background_DNA_LT_sperm_PFOA)
colSums(is.na(df_background_DNA_LT_sperm_PFOA)) 




####删除缺失值

df_background_DNA_LT_sperm_PFOA %>% 
  dplyr::filter(!PFOA=="NA") ->df_background_DNA_LT_sperm_PFOA_noNA  ###得到的精浆和人群基线数据一样的数据

str(df_background_DNA_LT_sperm_PFOA_noNA)   ###836人

length(unique(df_background_DNA_LT_sperm_PFOA_noNA$number))

colSums(is.na(df_background_DNA_LT_sperm_PFOA_noNA)) 



####最终数据

str(df_background_DNA_LT_sperm_PFOA_noNA)

###将数据重新排序整理
df_background_DNA_LT_sperm_PFOA_noNA %>% 
  select(PFOS:PF3OUdS_11CL,mtDNAcn,STL,everything()) -> df_background_DNA_LT_sperm_PFOA_tidy


str(df_background_DNA_LT_sperm_PFOA_tidy)

###先取对数，再求平均值
####对PFAS取对数
lapply(df_background_DNA_LT_sperm_PFOA_tidy[,c(1:14)],log2) ->df_background_DNA_LT_sperm_PFOA_tidy_lg 

### PFAS和DNA均取对数了
df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,25,26)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

length(unique(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number))

length(unique(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number))

which(duplicated(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number))

# write.xlsx(df_background_DNA_LT_sperm_PFOA_tidy_final_lg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_DNA_LT_sperm_pfas_tidy_final_lg.xlsx")  ###此处的pfas和端粒长度均取对数了，并未进行标化或归一化







####>>>>>>>>读入代谢物匹配的代谢通路变量名称#############################
# 
# 
# df_kegg<- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/240823_最终确定版data_lgm.xlsx",sheet = "合并")
# str(df_kegg) 
# 
# 
# 
# df_kegg %>% 
#   dplyr::rename(metabolite="META") -> df_kegg_1
# 
# str(df_kegg_1)
# 
# colSums(is.na(df_kegg_1))

###====================读入代谢组学数据===============================

df_meta<- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/240823_最终确定版data_lgm.xlsx",sheet = "已核对")
str(df_meta) 


df_meta %>% 
  dplyr::rename(number="ID") -> df_meta_1

str(df_meta_1)

names(df_meta_1)

####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)

df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2


str(df_unique_2)   ###一共有357个确定的代谢物


df_unique_2$metabolite


df_meta_1

names(df_meta_1)


###筛选出最终的代谢物
df_meta_1 %>% 
  select(number,df_unique_2$metabolite)  ->  df_meta_2  ####最终得到378个代谢物质

str(df_meta_2)

794-378

437-23

df_meta_2


which(duplicated(df_meta_2$number))   ###重复的有335个样本

which(!duplicated(df_meta_2$number))   ###不重复的有966个样本,去除重复的，并统计唯一个数



library(dplyr)
# 
# # 创建数据框
# df <- data.frame(number = c(1, 2, 2, 3,5,6,6,6,7,7))
# 
# # 去重并统计唯一变量的个数
# df_meta_2 %>%
#   distinct(number) %>%
#   nrow()
# 
# df  %>%
#   distinct(number) %>%
#   nrow()
# 



df_meta_2 %>%
  group_by(number) %>%
  summarise(count = n()) %>%
  filter(count == 1) %>%
  nrow()

####计算精浆代谢组学测量的研究天数及研究区间

###在本研究中计算研究区间的与精液质量的计算不一样
names(df_meta_1)   ###这里的dateaq是基线时期的日期，sampledaate是测量的样品收集的日期

df_meta_1 %>% 
  mutate(datesq=lubridate::ymd(sampledate)) %>% 
  arrange(number, datesq) %>% 
  group_by(number) %>% 
  mutate(firstday = lubridate::ymd(dateaq)) %>% 
  mutate(days = datesq - firstday) %>%  #计算每个ID下每次随访时间与基线时间间隔的天数
  ungroup()  -> df_dayds

df_dayds$days <- as.numeric(df_dayds$days)
as.numeric(df_dayds$days)

df_dayds$days


# 
# df_meta_1%>% 
#   mutate(datesq=lubridate::ymd(sampledate)) %>% 
#   arrange(number, datesq) %>% 
#   group_by(number) %>% 
#   mutate(firstday = min(datesq)) %>% 
#   mutate(days = datesq - firstday) %>%  #计算每个ID下每次随访时间与基线时间间隔的天数
#   ungroup()  -> df_background_repeatsemen_days   #得到含有时间天数的变量
# 

# write.xlsx(df_background_repeatsemen_days,"result/20240628_Hyperten_sperm_2th/20240628_result/df_background_repeatsemen_days.xlsx")


# str(df_background_repeatsemen_days)

# as.numbers(df_background_repeatsemen_days$days)

# df_background_repeatsemen_days$days <- as.numeric(df_background_repeatsemen_days$days) #将时间变量转换为数字变量


####根据随访与基线时间间隔的天数，按精子发生时间，划分间隔
df_dayds %>%  
  mutate(
    spermatogenesis=case_when(
      days>=0 & days<=9 ~0,
      days>=10 & days<=14 ~1,
      days>=15 & days<=69 ~2,
      days>=70 &days<=90 ~3
    )
  ) -> df_background_repeatsemen_final1 #生成最终的样本


df_background_repeatsemen_final1$spermatogenesis

table(df_background_repeatsemen_final1$spermatogenesis)

# write.xlsx(df_background_repeatsemen_final1,"result/20240628_Hyperten_sperm_2th/20240628_result/df_background_repeatsemen_final1.xlsx")

####第二种划分方法，根据随访与基线时间间隔的天数，按调查研究天数划分 


df_background_repeatsemen_final1 %>%  
  mutate(
    interval=case_when(
      days==0~0,
      days>=1 & days<=15 ~1,
      days>=16 & days<=31 ~2,
      days>=32 &days<=63 ~3,
      days>=64~4
    )
  ) -> df_background_repeatsemen_final2 #生成最终的样本

str(df_background_repeatsemen_final2)
df_background_repeatsemen_final2$interval

table(df_background_repeatsemen_final2$interval)

# 
# 0    1    2    3    4 
# 1164  838  800 1115 1814 
# write.xlsx(df_background_repeatsemen_final2,"result/20240318_Hyperten_sperm/df_background_repeatsemen_final2.xlsx")


colSums(is.na(df_background_repeatsemen_final2)) 




# write.xlsx(df_background_repeatsemen_final2,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_repeatsemen_final2.xlsx")












###将剩余的物质与PFAS和STL进行合并

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

df_background_DNA_LT_sperm_PFOA_tidy_final_lg  %>% 
  left_join(df_meta_2,by="number") -> df_background_DNA_LT_sperm_PFOA_meta_tidy

colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy))




####删除缺失值

df_background_DNA_LT_sperm_PFOA_meta_tidy  %>% 
  dplyr::filter(!M33=="NA")  %>% 
  as_tibble()->df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA  ###得到的精浆和人群基线数据一样的数据



colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA))

df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA


names(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA)###得到849个重复测量的样本



####计算重复变量的类内相关系数 

df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA

names(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA)


data <- df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA %>% 
  select(number,M642:M625) %>% 
  as.data.frame()

data

# write.xlsx(data,"D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data-meta.xlsx")







############1.读入stata计算的ICC结果——623个个体836次精浆样本测量

library(haven)

###回归系数
data_A <- read_dta("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_meta/A.dta")
data_A %>% 
  as_tibble()


###组间方差 Between-person σ2 (%)
data_B <- read_dta("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_meta/B.dta")

data_B %>% 
  as_tibble()

data_B %>% 
  dplyr::rename(Between_person="b") %>% 
  select(newv,Between_person) -> data_B_1

data_B %>% 
  left_join(data_C,by="newv")



###组内方差  Within-person σ2 (%)
data_C <- read_dta("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_meta/C.dta")
data_C %>% 
  as_tibble()


data_C %>% 
  dplyr::rename(Within_person="b") %>% 
  select(newv,Within_person) -> data_C_1


data_B_1 %>% 
  left_join(data_C_1,by="newv") -> data_b_c



###ICC   ICC (95% CI)
data_D <- read_dta("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_meta/D.dta")


data_D %>% 
  as_tibble() %>% 
  mutate(ICC_95=paste0(round(icc,digits = 2)," (",round(l,2),', ',round(u,2),")")) %>% 
  select(newv,icc,ICC_95) -> data_D_1


data_b_c %>% 
  left_join(data_D_1,by="newv") ->data_icc 



# write.xlsx(data_icc ,"D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc .xlsx")

# .0000159/(.0000159+245.8936)




################2.匹配32个精浆样本的2次重复测量个体的样本的一致性###########


###
data_number <- read_excel("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_number/data_number.xlsx")

data_number %>% 
  select(names(data)) -> data_number_meta

# write.xlsx(data_number_meta,"D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_number/data_number_meta.xlsx")




###组间方差 Between-person σ2 (%)
data_B <- read_dta("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_number/B.dta")

data_B %>% 
  as_tibble()

data_B %>% 
  dplyr::rename(Between_person="b") %>% 
  select(newv,Between_person) -> data_B_1

data_B %>% 
  left_join(data_C,by="newv")



###组内方差  Within-person σ2 (%)
data_C <- read_dta("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_number/C.dta")
data_C %>% 
  as_tibble()


data_C %>% 
  dplyr::rename(Within_person="b") %>% 
  select(newv,Within_person) -> data_C_1


data_B_1 %>% 
  left_join(data_C_1,by="newv") -> data_b_c



###ICC   ICC (95% CI)
data_D <- read_dta("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_number/D.dta")


data_D %>% 
  slice(-1,-355) %>% 
  as_tibble() -> data_D_0 

data_D_0$icc <- as.numeric(data_D_0$icc)  

data_D_0
data_D_0  %>% 
dplyr::mutate(ICC_95=paste0(round(icc,digits = 2)," (",round(l,2),', ',round(u,2),")")) %>% 
  select(newv,icc,ICC_95) -> data_D_1


data_b_c %>% 
  left_join(data_D_1,by="newv") ->data_icc 



# write.xlsx(data_icc ,"D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/repeated_number/data_icc.xlsx")  

# .0000159/(.0000159+245.8936)



#####3.计算ICC的范围 


###3.1计算number
data_number <- read_excel("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc_20250222_number.xlsx")

data_number

range(data_number$icc)
median(data_number$icc)


data_number %>% 
  dplyr::filter(icc>=0.5)


216/357     ####60.50%  超过60.50%的代谢物质重测的一致性大于或等于0.5

###3.1计算repeated_meta

data_repeated_meta <- read_excel("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc_20250222_repeated_meta.xlsx")

data_repeated_meta

data_repeated_meta$icc


data_repeated_meta %>% 
  dplyr::filter(icc>=0.5)

218/357     ###61.06% 大于0.5,超过61.06%代谢物在样本的重复测量收集中类内相关系数大于或等于0.5，表明重复测量能够有效的降低由于每次个体精液物质不同，引起的精浆代谢物的不同








