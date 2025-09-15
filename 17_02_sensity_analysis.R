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


# write.xlsx(df_sperm_PFOA_tidy_4,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_sperm_PFOA_tidy_4_1053.xlsx")

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


###将精浆中相同样本测量的第一次的精浆代谢物浓度替换为第二次的测量的精浆代谢物浓度   

df_meta_sample_repeated <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/240823_最终确定版data_lgm.xlsx",sheet = "重复测量-可做敏感性分析")
str(df_meta_sample_repeated) 


df_meta_sample_repeated %>% 
  dplyr::rename(number="ID") -> df_meta_sample_repeated_1


df_meta_sample_repeated_1


df_meta_1 %>% 
  # select(number,sampleID) %>% 
  left_join(df_meta_sample_repeated_1[,c(1,12)],by="sampleID") %>% 
   dplyr::filter(is.na(number.y)) -> df_meta_2     ###这些样本是第一次测量中单独有的样本

names(df_meta_2)

df_meta_2 %>% 
  select(number.x,sampleID,M1:M794) -> df_meta_3    ####重新提取精浆代谢物的主要变量



df_meta_sample_repeated_1 %>% 
  # select(number,sampleID) %>% 
  left_join(df_meta_1[,c(1,23)],by="sampleID")  %>% 
  dplyr::filter(is.na(number.y)) -> df_meta_sample_repeated_2     ###这些样本是第二次测量中单独有的样本



names(df_meta_sample_repeated_2)


df_meta_sample_repeated_2 %>% 
  select(number.x,sampleID,M1:M794) -> df_meta_sample_repeated_3 ####重新提取精浆代谢物的主要变量


names(df_meta_sample_repeated_3)

rbind(df_meta_3,df_meta_sample_repeated_3) -> df_meta_4     ###重新合成的代谢组学重复测的数据

df_meta_4 %>% 
  dplyr::rename(number=number.x) -> df_meta_5

df_meta_5



####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)

df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2


str(df_unique_2)   ###一共有357个确定的代谢物


df_unique_2$metabolite

# 
# df_meta_1
# 
# names(df_meta_1)


###筛选出最终的代谢物
df_meta_5 %>% 
  select(number,df_unique_2$metabolite)  ->  df_meta_6  ####最终得到378个代谢物质

str(df_meta_6)




794-378

437-23



which(duplicated(df_meta_6$number))   ###重复的有335个样本

which(!duplicated(df_meta_6$number))   ###不重复的有966个样本,去除重复的，并统计唯一个数







###将剩余的物质与PFAS和STL进行合并

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

df_background_DNA_LT_sperm_PFOA_tidy_final_lg  %>% 
  left_join(df_meta_6,by="number") -> df_background_DNA_LT_sperm_PFOA_meta_tidy

colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy))




####删除缺失值

df_background_DNA_LT_sperm_PFOA_meta_tidy  %>% 
  dplyr::filter(!M33=="NA")  %>% 
  as_tibble()->df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA  ###得到的精浆和人群基线数据一样的数据



colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA))

df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA


names(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA)###得到851个重复测量的样本


# length(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$number)
df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA[duplicated(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$number), ]  ###重复的行数

####>>>>>>>>>>>>>>>>>>>>>>一、导出重复的代谢组学数据，进行混合线性模型分析<<<<<<<<<<< 
# 
# write.xlsx(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_meta_mix_rawdata.xlsx")





####>>>>>>>>>>>>>>>>>>>>>>二、将重复的代谢组学数据，取对数进行归一化，然后再取对数，进行线性模型分析<<<<<<<<<<< 
###1.先取对数
names(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA)

df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA %>%
  select(M642:M625) %>%
  `+`(1) %>%
  log2()%>%
  apply(2, function(x) {     ###对每个代谢物进行z-score标化
    (x - mean(x)) / sd(x) 
  }) -> df1




str(df1)
# log(62.8,2)
# 
# 

# ###标化后的数据
df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA  %>%
  select(-c(M642:M625))%>%
  cbind(df1) -> df2_meta_stl







###2.再求均值

str(df2_meta_stl)
names(df2_meta_stl)


df2_meta_stl$STL
library(dplyr)
library(rlang)



# 假设 df_background_DNA_LT 是你的数据框
res2 <- lapply(names(df2_meta_stl[, c(11:381)]), function(x) {    ###这里对PFAS和代谢数据均求均值了，但是PFAS数据并未进行标化
  # 获取列名
  result <- df2_meta_stl %>%
    dplyr::group_by(number) %>%
    # summarise(mean_value = mean(!!sym(x), na.rm = TRUE)) # 使用 !!sym 来动态引用列
    summarise(mean_value = mean(!!sym(x), na.rm = TRUE), .groups = "drop") %>%
    mutate(variable = x)
  return(result)
})




res2 %>% 
  bind_rows() -> df_log_mean



# 3. 使用 pivot_wider() 将长数据转换为宽数据
wide_data <- df_log_mean %>%
  pivot_wider(names_from = variable, values_from = mean_value) %>% 
  dplyr::rename("number"="number")

names(wide_data)

# wide_data$number <- as.character(wide_data$number)

colSums(is.na(wide_data))

str(wide_data)   ###实际测量精液样本的人数只有623人

## 4.将基线数据与代谢数据进行合并  

# df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number <- as.character(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number)

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)   ###实际只有836人同时有端粒和PFAS测量数据




df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  select(-c(PFOS:STL)) %>%
  left_join(wide_data,by="number") ->df_background_DNA_LT_meta_mean





str(df_background_DNA_LT_meta_mean)




head(df_background_DNA_LT_meta_mean)

names(df_background_DNA_LT_meta_mean)
class(df_background_DNA_LT_meta_mean$M1)

colSums(is.na(df_background_DNA_LT_meta_mean))

df_background_DNA_LT_meta_mean %>% 
  # mutate(M1 = as.character(M1)) %>%  # 转换为字符类型
  filter(!is.na(M1))  -> df_background_DNA_LT_meta_mean_final  # 过滤掉 NA

str(df_background_DNA_LT_meta_mean_final)  ###实际只有623人有端粒和PFAS和精液代谢组学数据
names(df_background_DNA_LT_meta_mean_final)




# write.xlsx(df_background_DNA_LT_meta_mean_final,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/result_pfas_rawdata_repeated_623.xlsx")  ###此处的pfas和端粒长度均取对数了，并未进行标化或归一化





####>>>>>>>>>>>>>>>二、读入重复测量求均值模型的pfas_meta数据<<<<<<<<<<<<<<<#####


load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/coefs_PFAS_meta_repeated_sensity.Rdata")

g_coefs_PFAS_meta_repeated_sensity

lapply(seq_along(g_coefs_PFAS_meta_repeated_sensity), function(i) {
  g_coefs_PFAS_meta_repeated_sensity[[i]]  %>% 
    as.data.frame() 
  # res1 %>% 
  #   dplyr::filter(var.names=="pfos.beta"|var.names=="pfhxs.beta"|var.names=="pfhps.beta"|var.names=="pfoa.beta"|var.names=="pfna.beta"|var.names=="pfda.beta") -> res2
  # BCI_CI95 <- paste0(round(res1$Mean,2)," (",round(res1$X2.5.,2),', ',round(res1$X97.5.,2),")")
  # Uni_glm_model <- data.frame(res1,
  #                             "BCI_CI95"=BCI_CI95
  #                           )
  # return(Uni_glm_model)
}) %>% 
  bind_rows() %>%
  write.csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/g_coefs_PFAS_meta_repeated_sensity.csv")






#将导出的数据重新读入数据  

lapply(seq_along(g_coefs_PFAS_meta_repeated_sensity), function(i) {
  g_coefs_PFAS_meta_repeated_sensity[[i]]  %>% 
    as.data.frame() 
  # res1 %>% 
  #   dplyr::filter(var.names=="pfos.beta"|var.names=="pfhxs.beta"|var.names=="pfhps.beta"|var.names=="pfoa.beta"|var.names=="pfna.beta"|var.names=="pfda.beta") -> res2
  # BCI_CI95 <- paste0(round(res1$Mean,2)," (",round(res1$X2.5.,2),', ',round(res1$X97.5.,2),")")
  # Uni_glm_model <- data.frame(res1,
  #                             "BCI_CI95"=BCI_CI95
  #                           )
  # return(Uni_glm_model)
}) %>% 
  bind_rows() %>% 
  as_tibble() -> coefs_tran

# coefs_tran <- read_csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/coefs_pfas_meta.csv")





####将生成的代谢物传递给新的数据
coefs_tran ->coefs_tran_tidy_0

str(coefs_tran_tidy_0)

length(coefs_tran_tidy_0$metabolite)
table(coefs_tran_tidy_0$metabolite) 
coefs_tran_tidy_0$metabolite

coefs_tran_tidy_0[coefs_tran_tidy_0$metabolite %in% coefs_tran_tidy_0$metabolite[duplicated(coefs_tran_tidy_0$metabolite)], ]


coefs_tran_tidy_0$metabolite


###取出代谢物质单独的变量名称
coefs_tran_tidy_0 %>%
  group_by(metabolite) %>%
  dplyr::summarise(count=n())  -> coefs_tran_tidy_unique

coefs_tran_tidy_unique$metabolite

names(coefs_tran_tidy_unique)

str(coefs_tran_tidy_unique)

###将代谢物重新传递给coefs_tran_tidy

coefs_tran_tidy_0 ->coefs_tran_tidy 

coefs_tran_tidy

results_df <- coefs_tran_tidy  %>% 
  # select(-1) %>% 
  dplyr::mutate(p_value_estimate = na_if(coefs_tran_tidy$p.val, NA) %>% as.numeric)   ###去除p值中含有缺失的值

results_df



## Subset eta ------------------------------------
results_df_2 <-results_df %>%
  # dplyr::filter(str_detect(var.names,"eta")) %>%
  dplyr::rename(estimate = Mean,sd=SD,estimate_lcl=X2.5., estimate_ucl=X97.5.)
# mixtures_eta

# estimate, sd, lcl, ucl, p_value

results_df_2$var.names

## Pivot results wider ---------------------------------


# ## Pivot results wider ---------------------------------
# results_df_w <- results_df_2 %>% 
#   filter(var.names != "eta") %>%
#   tidylog::pivot_wider(id_cols = c(cohort, mixture, metabolite, exposure), 
#                        names_from = term, 
#                        values_from = c(estimate, sd, lcl, ucl, p_value)) %>% 
#   dplyr::select(cohort, mixture, metabolite, exposure, 
#                 contains("beta"), contains("pip")) %>% 
#   dplyr::select(-p_value_pip) %>% 
#   dplyr::rename(p_value = p_value_beta)
# 
# 


###--------转换BCI CI95-----------------------


results_df_3 <- results_df_2 %>% 
  dplyr::filter(var.names %in% c("eta.high","eta.low","psi"))  ###提取出eta.high、eta.low、psi三个重要变量

results_df_3

results_df_3 %>%  
  filter(var.names=="psi") -> results_df_4

paste0(round(results_df_4$estimate,2)," (",round(results_df_4$estimate_lcl,2),', ',round(results_df_4$estimate_ucl,2),")") %>% 
  as_tibble() %>% 
  dplyr::rename(BCI_CI95=value)  %>% 
  cbind(results_df_4)-> results_df_BCI_95CI  ###根据需要将想要的变量提取处理 

str(results_df_BCI_95CI)

str(results_df_2) 
results_df_2

results_df_2 %>% 
  dplyr::filter(!var.names %in% c("eta.high","eta.low","psi")) %>% 
  dplyr::select(metabolite,everything())->results_df_5   ###过滤掉不是eta.high、eta.low、psi三个重要变量


###将转换后的BCI_95与不含eta.high、eta.low、psi三个变量进行合并
results_df_5

results_df_5 %>% 
  left_join(results_df_BCI_95CI,by="metabolite") %>% 
  dplyr::filter(var.names.x %in% c("PFOS.beta","PFOA.beta","PFDA.beta","PFUdA.beta","L_PFHxS.beta","PF3ONS_9CL.beta"))-> results_df_6   ###合并后的数据  ###只保留beta变量,每个化合物的beta和混合物beta

results_df_6 
names(results_df_6)


###-------------只保留后验概率PIP数值-------------------------
results_df_5 %>% 
  dplyr::filter(var.names %in% c("PFOS.gamma","PFOA.gamma","PFDA.gamma","PFUdA.gamma","L_PFHxS.gamma","PF3ONS_9CL.gamma"))  %>% 
  tidylog::pivot_wider( id_cols = c(metabolite), 
                        names_from = var.names,
                        values_from = c(estimate)) -> results_df_7  ###保存每个化合物的后验概率

results_df_7   




# Calculate new p values


# Calculate new p values  计算单一物质的wald检验
names(results_df_6)
results_df_6 %>% 
  dplyr::mutate(wald = (estimate.x/sd.x), 
                # p = 2*(1-pnorm(abs(wald),0,1)),
                p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
                p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
                neg_log_p = -log10(p_value))  %>% 
  group_by(var.names.x) %>% 
  dplyr::mutate(q_value=p.adjust(p_value,method = "fdr"),
                significance = if_else(p_value < 0.05,
                                       "p < 0.05", 
                                       "Not Sig."), 
                significancefdr = if_else(q_value < 0.05,
                                          "q < 0.05",
                                          "Not Sig.")) %>% 
  ungroup() -> results_df_only_p_value

results_df_only_p_value

names(results_df_only_p_value)

write.xlsx(results_df_only_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/results_df_only_p_value.xlsx") 




# result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_pfas_meta





########---------------计算混合物质的wald检验---------------------
# Calculate new p values   
# 
results_df_BCI_95CI
names(results_df_BCI_95CI)

results_df_BCI_95CI %>% 
  dplyr::mutate(wald = (estimate/sd), 
                # p = 2*(1-pnorm(abs(wald),0,1)),
                p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
                p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
                neg_log_p = -log10(p_value)) %>% 
  # group_by(metabolite) %>% 
  dplyr::mutate(q_value = p.adjust(p_value,
                                   method = "fdr"), 
                significance=if_else(p_value < 0.05,
                                     "p < 0.05", 
                                     "Not Sig."), 
                significancefdr=if_else(q_value < 0.05,
                                        "q < 0.05",
                                        "Not Sig.")) %>% 
  ungroup() -> results_df_mix_p_value

results_df_mix_p_value

names(results_df_mix_p_value)

write.xlsx(results_df_mix_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/results_df_all_p_value.xlsx")


# 
# q_value = p.adjust(c(0.896,0.928,0.868),
#                    method = "fdr")



###--------------将BCI_CI95与后验概率PIP进行合并-------------


results_df_mix_p_value %>% 
  left_join(results_df_7 ,by="metabolite") %>% 
  dplyr::select(metabolite,BCI_CI95,p_value,q_value,PFOS.gamma:PF3ONS_9CL.gamma,everything())-> results_df_final   ###得到复现文章的Supplemental Excel file文件的结果

write.xlsx(results_df_final,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/results_df_final.xlsx")





names(results_df_only_p_value)


results_df_only_p_value %>% 
  mutate(associated_bci = if_else(estimate_lcl.x > 0 | estimate_ucl.x < 0,
                                  "Associated", 
                                  "Not Associated")) %>% 
  group_by(metabolite, var.names.x) %>% 
  dplyr::summarise(n_features = length(metabolite), 
                   percent_significant_p05 = jag2::npct(significance, 
                                                        "p < 0.05", 
                                                        n.digits = 2), 
                   percent_significant_q05 = jag2::npct(significancefdr, 
                                                        "q < 0.05", 
                                                        n.digits = 2), 
                   pct_associated = jag2::npct(associated_bci, 
                                               "Associated", 
                                               n.digits = 2))  -> mwas_summary


mwas_summary

# pivot summary wider   
mwas_summary_w <- pivot_wider(mwas_summary, 
                              id_cols = metabolite,
                              names_from = var.names.x, 
                              values_from = c(percent_significant_p05, 
                                              percent_significant_q05, 
                                              pct_associated)) 
mwas_summary_w 
# select(exposure, mixture, contains("SOL"), contains("CHS"), contains("Pooled"))  

# write.xlsx(mwas_summary_w,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/results_p_summary_final.xlsx")

###============================画火山图=======================================
####1.###可视化每个化合物有统计学意义的代谢物的火山图
# "Summary of Mixtures MWAS Results_hyper_g_v3.csv")

# -----------SOLAR Volcano Plots ----------------------------------------------
# mwas_results <- read_csv(
#   file = fs::path(dir_results_mixtures, 
#                   "all_pfas_mixtures_results_hyper_g.csv")) %>% 
#   as_tibble()

# Get reduced dataset   ####筛选限定值范围的差异话代谢物
###1.方式画火山图 

# results_df_only_p_value


####>>>>>>>>>>>>>读入代谢物匹配的代谢通路变量名称>>>>>>>>>>>>###################


colSums(is.na(results_df_only_p_value))
range(results_df_only_p_value$q_value)

mwas_reduced <- results_df_only_p_value %>% 
  dplyr::filter(q_value  < 1)
names(mwas_reduced)


####>>>>>>>将代谢物的名称进行合并
####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2

str(df_unique_2)   ###一共有378个确定的代谢物


df_unique_2$metabolite

mwas_reduced$metabolite


mwas_reduced  %>% 
  left_join(df_unique_2,by="metabolite")  -> mwas_reduced_2

colSums(is.na(mwas_reduced_2))




# dev.off()

###2.方式画火山图

library(showtext)
library(sysfonts)

font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

names(results_df_only_p_value)
results_df_only_p_value$var.names




###由于有些化合物可能没有统计学意义，会导致有些化合物标示
mwas_reduced_2 %>% 
  mutate(marker=case_when(
    q_value< 0.05 & estimate.x > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.05 & estimate.x< 0 ~ "Down" ,   ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  ))  -> res1

mwas_reduced_2 %>% 
  mutate(marker=case_when(
    q_value< 0.05 & estimate.x > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.05 & estimate.x< 0 ~ "Down" ,   ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>% 
  mutate(var.names.x=case_when(
    var.names.x=="PFOS.beta"~"PFOS",
    var.names.x=="PFOA.beta"~"PFOA",
    var.names.x=="PFDA.beta"~"PFDA",
    var.names.x=="L_PFHxS.beta"~"PFHxS",
    var.names.x=="PFUdA.beta"~"PFUdA",
    var.names.x=="PF3ONS_9CL.beta"~"6:2 Cl-PFESA"
  )) %>% 
  ggplot(aes(x =estimate.x, y = -log2(q_value),color= marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
  # xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线 
  geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
  # facet_wrap(~var.names, scales = "free") +
  # labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  # ggtitle("单组火山图") + #标题
  labs(x = "Estimate", y = "-log2(q_value)")+
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
  facet_wrap(~var.names.x, scales = "free") +
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) -> volcano_plot_only

volcano_plot_only




# font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/volcano_only_PFAS.pdf",width = 16, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_only #出图
par(opar)
dev.off()



####2.可视化混合化合物有统计学意义的代谢物的火山图

results_df_mix_p_value

names(results_df_mix_p_value)

results_df_mix_p_value$p_value





mwas_reduced_mix <-results_df_mix_p_value %>% 
  dplyr::filter(q_value  < 1)
names(mwas_reduced_mix)
# Volcano Plot     ####画火山图
# # 
# ggplot(mwas_reduced,aes(x = estimate, y = -log10(p_value),color=significance)) +
#   geom_point(size = 1, alpha = 0.5) + 
#   geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
#   # geom_hline(yintercept = -log10(q_value), lty=4,col="black",lwd=0.8)+
#   facet_wrap(~var.names, scales = "free") +
#   ylab("-log10 p") +
#   xlab("Effect Estimate")

# dev.off() 
# 


####>>>>>>>匹配变量名称<<<<<<<<<<<<<<<<<<###########################
####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2

str(df_unique_2)   ###一共有378个确定的代谢物


df_unique_2$metabolite

mwas_reduced$metabolite



mwas_reduced_mix  %>% 
  left_join(df_unique_2,by="metabolite")  -> mwas_reduced_mix_2

# 
library(showtext)
library(sysfonts)
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置


names(mwas_reduced_mix_2)


volcano_plot_mix <-
  mwas_reduced_mix_2 %>% 
  mutate(marker=case_when(
    q_value < 0.05 & estimate > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.05 & estimate< 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate, y = -log(q_value,10),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, s llllllllize=2) +  #点的透明度、大小
  geom_point(aes(size = -log(q_value,10),color = marker),alpha = 0.7) +
  # ylim(c(0,75))+
  # xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线 
  # geom_hline(yintercept = 0, linetype = 2,size=1) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
  # facet_wrap(~var.names, scales = "free") +
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  # ggtitle("单组火山图") + #标题
  labs(x = "PFAS and Metabolism correlation", y = "-log10(q_value)")+
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
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) 


volcano_plot_mix





font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/volcano_all_PFAS.pdf",width = 12, height = 10, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_mix #出图

par(opar)
dev.off()


####二、对端粒进行标化之后再取均值



###2.再求均值

str(df2_meta_stl)
names(df2_meta_stl)


df2_meta_stl$STL
library(dplyr)
library(rlang)





###1.1对端粒进行归一化---在原始的变量上求均值 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  # select(mtDNAcn,STL) %>%
  select(mtDNAcn,STL) %>%
  apply(2, function(x) {     ###对端粒和线粒体拷贝数进行z-score标化
    (x - mean(x)) / sd(x)
  }) %>%
  as_tibble() -> df2_meta_stl_zscore



###1.2.将数据进行整合

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  select(-c(mtDNAcn,STL)) %>% 
  cbind(df2_meta_stl_zscore) %>% 
  select(number:PF3OUdS_11CL,mtDNAcn,STL) -> df_background_pfas_stl_zscore   ###PFAS和端粒长度均进行标化了




df_background_pfas_stl_zscore %>% 
  left_join(df2_meta_stl[,-c(2:24)],by="number") -> df_background_pfas_stl_zscore_final

str(df_background_pfas_stl_zscore_final)

colSums(is.na(df_background_pfas_stl_zscore_final))

df_background_pfas_stl_zscore_final %>% 
  dplyr::filter(!M33=="NA")  %>% 
  as_tibble()->df_background_pfas_stl_zscore_final_tidy  ###得到的精浆和人群基线数据一样的数据

str(df_background_pfas_stl_zscore_final_tidy)


names(df_background_pfas_stl_zscore_final_tidy)







# 假设 df_background_DNA_LT 是你的数据框
res10<- lapply(names(df2_meta_stl[, c(11:381)]), function(x) {    ###这里对PFAS和代谢数据均求均值了，但是PFAS数据并未进行标化
  # 获取列名
  result <- df2_meta_stl %>%
    dplyr::group_by(number) %>%
    # summarise(mean_value = mean(!!sym(x), na.rm = TRUE)) # 使用 !!sym 来动态引用列
    summarise(mean_value = mean(!!sym(x), na.rm = TRUE), .groups = "drop") %>%
    mutate(variable = x)
  return(result)
})




res10 %>% 
  bind_rows() -> df_log_mean



# 3. 使用 pivot_wider() 将长数据转换为宽数据
wide_data <- df_log_mean %>%
  pivot_wider(names_from = variable, values_from = mean_value) %>% 
  dplyr::rename("number"="number")

names(wide_data)

# wide_data$number <- as.character(wide_data$number)

colSums(is.na(wide_data))

str(wide_data)   ###实际测量精液样本的人数只有623人

## 4.将基线数据与代谢数据进行合并  

# df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number <- as.character(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number)

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)   ###实际只有836人同时有端粒和PFAS测量数据




df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  select(-c(PFOS:STL)) %>%
  left_join(wide_data,by="number") ->df_background_DNA_LT_meta_mean





str(df_background_DNA_LT_meta_mean)


head(df_background_DNA_LT_meta_mean)

names(df_background_DNA_LT_meta_mean)
class(df_background_DNA_LT_meta_mean$M1)

colSums(is.na(df_background_DNA_LT_meta_mean))

df_background_DNA_LT_meta_mean %>% 
  # mutate(M1 = as.character(M1)) %>%  # 转换为字符类型
  filter(!is.na(M1))  -> df_background_DNA_LT_meta_mean_final_zscore  # 过滤掉 NA

str(df_background_DNA_LT_meta_mean_final_zscore)  ###实际只有623人有端粒和PFAS和精液代谢组学数据
names(df_background_DNA_LT_meta_mean_final_zscore)









###>>>>>>>>>>>>>>>>>>>>>>三、读入经过BHRMA筛选后的pfas混合物---混合线性模型的<<<<<<<<<<<<<<<<<<

df_pfas_repeated_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/results_df_all_p_value.xlsx")


###1.提取经过BH校正之后q_value<0.2的代谢物 
###
names(df_pfas_repeated_meta)

df_pfas_repeated_meta %>% 
  filter(q_value<0.2) -> df_pfas_repeated_meta_q   ###总共获得95个代谢物质   


df_pfas_repeated_meta_q

str(df_pfas_repeated_meta_q)
names(df_pfas_repeated_meta_q)



# 
# ####>>>>>>>>读入代谢物匹配的代谢通路变量名称#############################


###>>>>>>>>>>>>>>>>四、过滤出与PFAS有相关性的代谢物质，校正后q<0.20的物质，再进行回归分析<<<<<<<<<<<<<<<<<#########

names(df_background_DNA_LT_meta_mean_final_zscore)

df_background_DNA_LT_meta_mean_final_zscore



df_background_DNA_LT_meta_mean_final_zscore  %>% 
  select(number:STL,df_pfas_repeated_meta_q$metabolite) ->df_background_DNA_LT_q   ###过滤出来只与PFAS有关连的代谢物


str(df_background_DNA_LT_q)

names(df_background_DNA_LT_q)



lapply(df_background_DNA_LT_q[,c(25:89)], function(x){
  lm(STL~x+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,data=df_background_DNA_LT_q) -> result
  cox <- summary(result)
  cbind(cox$coefficients[,c(1,2,3,4)],confint(result)) -> result2
  β<-round(cox$coefficients[,1],digits = 2)
  CI95<-paste0(β,"(",round(result2[,5],2),',',round(result2[,6],2),")")
  percent_β <- (exp(cox$coefficients[,1])-1)*100
  percent_CI2.5 <- (exp(result2[,5])-1)*100 
  percent_CI97.5 <- (exp(result2[,6])-1)
  percent_CI95 <- paste0(round(percent_β,digits = 2)," (",round((exp(result2[,5])-1)*100,2),', ',round((exp(result2[,6])-1)*100,2),")")
  res <- result %>% tidy %>% data.frame()   ####自己添加的，将其他输出的变量也添加到包里
  Uni_glm_model <- data.frame('Estimate'=result2[,1],
                              "Std. Error"=result2[,2],
                              't value' = result2[,3],
                              'p value' = result2[,4],
                              'CI2.5' = result2[,5],
                              'CI97.5' = result2[,6],
                              'CI95' = CI95,
                              "percent_β"=percent_β,
                              "percent_CI2.5"=percent_CI2.5,
                              "percent_CI97.5"=percent_CI97.5,
                              "percent_CI95"=percent_CI95,
                              res)
  #返回循环函数继续上述操作                     
  return(Uni_glm_model)}) -> res9  ###将跑的循环中y生成的list结果保存为res1






res9 %>% 
  bind_rows()%>%
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/result_meta_STL_zscore_liner_PFAS_q.xlsx")
# 
# colnames(df_background_DNA_LT_q[,c(47:297)]) %>%
#   as_tibble()  %>% 
#   dplyr::rename("meta"="value") -> name_meta





colnames(df_background_DNA_LT_q[,c(25:89)]) %>%
  as_tibble()  %>% 
  dplyr::rename("meta"="value") -> name_meta


###导出所有的代谢物
# res1 %>% bind_rows()%>%as_tibble()%>%
#   write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/result_PFAS.xlsx")



###代谢物与端粒之间的混合线性模型结果
res9 %>% 
  bind_rows()%>%
  as_tibble() %>% 
  filter(term=="x") %>% 
  cbind(name_meta) -> df_meta_stl



p.adjust(df_meta_stl$p.value,method = "BH") %>% 
  as_tibble() %>% 
  dplyr::rename("q_value"="value")-> q_value

q_value
###代谢物与端粒之间的混合线性模型结果
df_meta_stl %>% 
  cbind(q_value) %>%
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/result_meta_STL_zscore_liner_q_PFAS.xlsx")


df_meta_stl %>% 
  cbind(q_value) -> result_meta_STL_zscore_liner_q_PFAS   #得到即是PFAS 又是端粒的代谢物质






#####======================进行通路富集分析=================  



###====================读入代谢组学数据===============================

####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2

str(df_unique_2)   ###一共有378个确定的代谢物



result_meta_STL_zscore_liner_q_PFAS %>% 
  left_join(df_unique_2,by=c("meta"="metabolite")) -> pfas_meta_STL_p_kegg

pfas_meta_STL_p_kegg

write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/result_meta_STL_zscore_liner_q_PFAS_kegg.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))






####>>>>>1.利用重构的kegg_pathway信号通路进行富集分析<<<<<<<<################


###>>>>>>1.1利用massdatabase下载、读取和转换kegg数据库
library(massdatabase)
library(metpath)

###下载
# download_kegg_pathway(path = "kegg_human_pathway",
#                       sleep = 1,
#                       organism = "hsa")


###读取
data <- 
  read_kegg_pathway(path = "kegg_human_pathway")







####这个函数不起作用的主要原因是pathway_id和pathway_name这里是list,不是向量

convert_kegg2metpath <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    pathway <-
      new(
        Class = "pathway_database",
        database_info = list(source = "KEGG",
                             version = as.character(Sys.Date())),
        pathway_id = unlist(
          lapply(data, function(x) {
            if (is.null(x$pathway_id) || length(x$pathway_id) == 0 || x$pathway_id == "") {
              return(NA)  # 如果 pathway_id 缺失，返回 NA
            } else {
              return(x$pathway_id)  # 否则返回实际值
            }
          }), use.names = FALSE),
        pathway_name =  unlist(
          lapply(data, function(x) {
            if (is.null(x$pathway_name) || length(x$pathway_name) == 0 || x$pathway_name == "") {
              return(NA)  # 如果 pathway_name 缺失，返回 NA
            } else {
              return(x$pathway_name)  # 否则返回实际值
            }
          }), use.names = FALSE),
        describtion = lapply(data, function(x)
          x$describtion),
        pathway_class = lapply(data, function(x)
          x$pathway_class),
        gene_list = lapply(data, function(x)
          x$gene_list),
        compound_list = lapply(data, function(x)
          x$compound_list),
        protein_list = list(),
        reference_list = list(),
        related_disease = lapply(data, function(x)
          x$related_disease),
        related_module = lapply(data, function(x)
          x$related_module)
      )
    
    if (length(pathway@gene_list) == 0) {
      pathway@gene_list <-
        vector(mode = "list",
               length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }
    
    if (length(pathway@compound_list) == 0) {
      pathway@compound_list <-
        vector(mode = "list",
               length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }
    
    if (length(pathway@protein_list) == 0) {
      pathway@protein_list <-
        vector(mode = "list",
               length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }
    kegg_pathway <- pathway
    rm(list = c("pathway"))
    save(kegg_pathway, file = file.path(path, "kegg_pathway"))
    invisible(kegg_pathway)
  }




kegg_human_pathway <-
  convert_kegg2metpath(data = data,
                       path = "kegg_human_pathway",
                       threads = 14)
kegg_human_pathway



pathway_class <-
  metpath::pathway_class(kegg_human_pathway)
# pathway_class

# get_pathway_class(result_case_control)


remain_idx <-
  pathway_class %>%
  unlist() %>% 
  stringr::str_detect("Disease") %>%
  `!`() %>%
  which()

pathway_database_lgm <-
  kegg_human_pathway[remain_idx]


###>>>>>>1.2导入待分析数据合并之后的kegg——id

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

kegg_id_1<- pfas_meta_STL_p_kegg%>% 
  select(id_kegg)%>%
  dplyr::filter(!id_kegg=="NA") %>% 
  as.data.frame()
colnames(kegg_id_1) <- NULL


kegg_id_2 <- kegg_id_1[!is.na(kegg_id_1)]
kegg_id_2 

result_case_control <- 
  enrich_kegg(query_id = kegg_id_2, 
              query_type = "compound", 
              id_type = "KEGG",
              pathway_database = pathway_database_lgm, 
              p_cutoff = 0.05, 
              p_adjust_method = "BH", 
              threads = 32)


result_case_control


library(clusterProfiler)


enrich_bar_plot(object = result_case_control,
                x_axis = "p_value",
                cutoff = 0.05) +
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )-> p1

p1




result_case_control@result %>% 
    as_tibble() %>% 
  write.csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/sensity/sensity_repeated/pfas_meta_STL_p_kegg_pathway.csv")








































