
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


log2(0.05)

log10(0.07)

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

# write.xlsx(df_background_DNA_LT_sperm_PFOA_tidy_final_lg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/df_background_DNA_LT_sperm_PFOA_tidy_final_lg.xlsx")  ###此处的pfas和端粒长度均取对数了，并未进行标化或归一化



####>>>>>>>>读入代谢物匹配的代谢通路变量名称#############################



df_meta<- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/240823_最终确定版data_lgm.xlsx",sheet = "已核对")
str(df_meta) 


df_meta %>% 
  dplyr::rename(number="ID") -> df_meta_1

str(df_meta_1)

names(df_meta_1)


df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2


######>>>>>>删除ICC<0.3的代谢物  

df_icc <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc_20250222_number.xlsx")
####重复测量的32个同样样本的代谢物特征的一致性

df_icc %>% 
  filter(icc>=0.3) -> df_icc_1   ###筛选出icc>0.3的代谢物特征


df_icc_1$metabolites


df_unique_2 %>% 
  filter(metabolite %in% df_icc_1$metabolites) -> df_unique_3




df_meta_1

names(df_meta_1)


###筛选出最终的代谢物
df_meta_1 %>% 
  select(number,df_unique_3$metabolite)  ->  df_meta_2  ####最终得到378个代谢物质

str(df_meta_2)

794-378

437-23

df_meta_2


which(duplicated(df_meta_2$number))
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

names(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA) ###得到849个重复测量的样本


# length(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$number)
# df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA[duplicated(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$number), ]  ###重复的行数


#########代谢物和端粒、线粒体拷贝数的整理
#### 1.先取对数
df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA %>%
  select(M661:M625) %>%
  `+`(1) %>%
  log2() %>%
  apply(2, function(x) {     ###对每个代谢物进行z-score标化
    (x - mean(x)) / sd(x)
  }) -> df1

# log(62.8,2)


#### 标化后的数据
df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA  %>%
  select(-c(M661:M625))%>%
  cbind(df1) -> df2_meta_stl     ###此处的代谢物是已经进行log2对数转化和标化后的物质，可以用作混合线性模型分析

str(df2_meta_stl)

df2_meta_stl   ###



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
  select(number:PF3OUdS_11CL,mtDNAcn,STL) -> df_background_pfas_stl_zscore   ###混合线性回归分析

df_background_pfas_stl_zscore$mtDNAcn

df_background_pfas_stl_zscore$STL

names(df_background_pfas_stl_zscore)   ###已经取对数和归一化之后的端粒

str(df_background_pfas_stl_zscore)

###将zscore变换之后的基线和端粒数据与zscore变换之后的代谢数据进行合并
str(df2_meta_stl)
names(df2_meta_stl)

df_background_pfas_stl_zscore %>% 
  left_join(df2_meta_stl[,-c(2:24)],by="number") -> df_background_pfas_stl_zscore_final

str(df_background_pfas_stl_zscore_final)

colSums(is.na(df_background_pfas_stl_zscore_final))

df_background_pfas_stl_zscore_final %>% 
  dplyr::filter(!M33=="NA")  %>% 
  as_tibble()->df_background_pfas_stl_zscore_final_tidy  ###得到的精浆和人群基线数据一样的数据

str(df_background_pfas_stl_zscore_final_tidy)

# write.xlsx(df_background_pfas_stl_zscore_final_tidy,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy.xlsx")



#######2.再求均值

str(df2_meta_stl)
names(df2_meta_stl)


df2_meta_stl$STL

library(dplyr)
library(rlang)


# 假设 df_background_DNA_LT 是你的数据框
res2 <- lapply(names(df2_meta_stl[, c(25:295)]), function(x) {
  # 获取列名
  result <- df2_meta_stl%>%
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

## df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number <- as.character(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number)

str(df_background_pfas_stl_zscore)   ###实际只有836人同时有端粒和PFAS测量数据


df_background_pfas_stl_zscore %>% 
  left_join(wide_data,by="number") ->df_background_DNA_LT_meta_mean

str(df_background_DNA_LT_meta_mean)

head(df_background_DNA_LT_meta_mean)

names(df_background_DNA_LT_meta_mean)
class(df_background_DNA_LT_meta_mean$M1)

df_background_DNA_LT_meta_mean %>% 
  # mutate(M1 = as.character(M1)) %>%  # 转换为字符类型
  filter(!is.na(M33))  -> df_background_DNA_LT_meta_mean_final  # 最终进行线性回归分析的只有623个变量



str(df_background_DNA_LT_meta_mean_final)  ###实际只有623人有端粒和PFAS和精液代谢组学数据
names(df_background_DNA_LT_meta_mean_final) 

df_background_DNA_LT_meta_mean_final
str(df_background_DNA_LT_meta_mean_final)


# 
# write.xlsx(df_background_DNA_LT_meta_mean_final,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy_623.xlsx")



# write.xlsx(df_background_pfas_stl_zscore_final_tidy,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy_623.xlsx")


# lm(STL~M33+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,data=df_background_DNA_LT_meta_mean_final) -> result 

# summary(result)




########>>>>>>>>>>>>>>>>>>一、对协变量进行调整<<<<<<<<<<<<<<<<<<<<<############

###利用取对数和归一化之后的端粒和代谢物，进行混合线性回归分析
df_background_DNA_LT_meta_mean_final    ###混合线性回归分析


names(df_background_DNA_LT_meta_mean_final)


# 
# ###4.进行混合线性回归分析
# 
# df2_meta_stl_tidy   ####此处的数据集，没有对端粒和线粒体拷贝数进行z-score变换
# names(df2_meta_stl_tidy)
# 
# 
# library(lmer)
# 
# str(df2_meta_stl_tidy)
# 
# # 检查共线性
# library(car)
# vif(lm(STL ~ age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + drk_class_2 + season + abstinence_class_2, data = df2_meta_stl_tidy))
# 
# ####教育与其他协变量可能存在多重共线性
# 
# 
# 
# 
# lapply(df2_meta_stl_tidy[,c(25:402)], function(x){
#   lmer(STL~x+age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + drk_class_2 + season + abstinence_class_2+(1|number),data=df2_meta_stl_tidy) -> result
#   summary(result) -> res1
#   res1$coefficients -> res2
#   β <- res2[,1]
#   SE <- res2[,2]
#   CI2.5 <- β-1.96*SE
#   CI97.5  <- β+1.96*SE
#   CI95<-paste0(round(β,2)," (",round(CI2.5,2),', ',round(CI97.5,2),")")
#   percent_β <- (exp(β)-1)*100
#   percent_CI2.5 <- (exp(β-1.96*SE)-1)*100
#   percent_CI97.5 <- (exp(β+1.96*SE)-1)*100
#   percent_CI95 <- paste0(round(percent_β,2)," (",round(percent_CI2.5,2),', ',round(percent_CI97.5,2),")")
#   res2 %>% 
#     rownames() %>% 
#     as_tibble() -> res ####自己添加的，将其他输出的变量也添加到包里
#   Uni_glm_model <- data.frame('Estimate'=res2[,1],
#                               "Std. Error"=res2[,2],
#                               'df' = res2[,3],
#                               't value' = res2[,4],
#                               'p value' = res2[,5],
#                               'CI2.5' = CI2.5,
#                               'CI97.5' = CI97.5,
#                               'CI95' = CI95,
#                               "percent_β"=percent_β,
#                               "percent_CI2.5"=percent_CI2.5,
#                               "percent_CI97.5"=percent_CI97.5,
#                               "percent_CI95"=percent_CI95,
#                               res)
#   return(Uni_glm_model)}) -> res7  ###将跑的循环中y生成的list结果保存为res1
# 


# 
# res7 %>% 
#   bind_rows()%>%
#   write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_meta_STL_log_mean_all_meta_final_20241201_zscore.xlsx")
# # 
# 




##构建object文件 
###1.1首先构建expression_data 
###
names(df_background_DNA_LT_meta_mean_final)   ###精液质量和代谢组学数据已经进行取对数和标化
str(df_background_DNA_LT_meta_mean_final)

df_background_DNA_LT_meta_mean_final[,c(1,25:295)] %>% 
  t() -> expression_data_1


##第一行作为列名
colnames(expression_data_1) <- expression_data_1[1,]

##删除第一行

expression_data_1[-1,]  %>% 
  as.data.frame()-> expression_data

names(df_background_DNA_LT_meta_mean_final)

###1.2首先构建sample_info
sample_info <- data.frame(df_background_DNA_LT_meta_mean_final[,c(1:24)] %>% 
                            dplyr::rename(sample_id="number"),class = "Subject",
                          subject_id = "Subject")

###1.3其次构建variable_info
variable_info <- rownames(expression_data) %>% as.data.frame() %>% 
  dplyr::rename(variable_id=".") 


object <- 
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

object






##find marker which are change according to aging

###linear regression
library(tidyverse)
library(ggpubr)
library(rstatix)


expression_data <-
  extract_expression_data(object)  %>% 
  # `+`(1) %>%
  # log(2) %>%
  # apply(1, function(x) {
  #   (x - mean(x)) / sd(x)
  # }) %>%
  # t() %>%
  as.data.frame()

library(plyr)

sample_info <-
  object@sample_info

expression_data


colSums(is.na(expression_data))


lm(M33 ~ age_class_2 + BMI_class_2 + education + 
     marriage_class_2 + income + smk_class_2 + drk_class_2 +
     season + abstinence_class_2 , data = df_background_DNA_LT_meta_mean_final) %>%
  residuals()


class(expression_data)

library(tidyverse)
library(tidymass)
str(sample_info)

expression_data %>% as.data.frame() %>% rownames()

library(plyr)
# as.numeric(expression_data[M33,])

# expression_data[1,]
# length(rownames(expression_data))
###探究协变量对代谢物质的影响

lm_adjust <- function(expression_data,
                      sample_info,
                      threads = 5) {
  library(future)
  library(furrr)
  # plan(strategy = multisession(workers = threads))
  lapply(seq_along(rownames(expression_data)), function(i){
    # cat(name, " ")
    x = as.numeric(expression_data[i,])
    temp_data =
      data.frame(x = x, sample_info)
    # 
    # temp_data$age_class_2[temp_data$age_class_2 =='<28'] = 0
    # temp_data$age_class_2[temp_data$age_class_2 =='≥28'] = 1
    # temp_data$age_class_2 = as.factor(temp_data$age_class_2)
    # 
    # temp_data$BMI_class_2[temp_data$BMI_class_2 == '<24'] = 1
    # temp_data$BMI_class_2[temp_data$BMI_class_2 == '≥24'] = 2
    # temp_data$BMI_class_2 = as.factor(temp_data$BMI_class_2)
    # 
    # temp_data$education[temp_data$education == 'High school'] = 1
    # temp_data$education[temp_data$education == 'College'] = 2
    # temp_data$education[temp_data$education == 'Undergraduate and above'] = 3
    # temp_data$education = as.factor(temp_data$education)
    # 
    # temp_data$marriage_class_2[temp_data$marriage_class_2 == 'Unmarried'] = 1
    # temp_data$marriage_class_2[temp_data$marriage_class_2 == 'Married'] = 2
    # temp_data$marriage_class_2 =as.factor(temp_data$marriage_class_2)
    # 
    # temp_data$income[temp_data$income == '<4000'] = 1
    # temp_data$income[temp_data$income == '4000-8000'] = 2
    # temp_data$income[temp_data$income == '>8000'] = 3
    # temp_data$income = as.factor(temp_data$income)
    # 
    # temp_data$smk_class_2[temp_data$smk_class_2 == 'Never'] = 0
    # temp_data$smk_class_2[temp_data$smk_class_2 == 'Former/Current'] = 1
    # temp_data$smk_class_2 = as.factor(temp_data$smk_class_2)
    # 
    # temp_data$drk_class_2[temp_data$drk_class_2 == 'Never'] = 1
    # temp_data$drk_class_2[temp_data$drk_class_2 == 'Former/Current'] = 2
    # temp_data$drk_class_2 = as.factor(temp_data$drk_class_2)
    # 
    # temp_data$season_class[temp_data$season_class == 'Spring'] = 1
    # temp_data$season_class[temp_data$season_class == 'Summer'] = 2
    # temp_data$season_class[temp_data$season_class == 'Autumn'] = 3
    # temp_data$season_class[temp_data$season_class == 'Winter'] = 4
    # temp_data$season_class = as.factor(temp_data$season_class)
    # 
    # temp_data$abstinence_class_2[temp_data$abstinence_class_2 == '<7'] = 1
    # temp_data$abstinence_class_2[temp_data$abstinence_class_2 == '≥7'] = 2
    # temp_data$abstinence_class_2 = as.factor(temp_data$abstinence_class_2)
    # # 
    
    adjusted_x <-
      lm(x ~ age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2, data = temp_data) %>%
      residuals()
    adjusted_x
  }) %>%
    bind_rows() %>%
    as.data.frame() ->  new_expression_data
  
  colnames(new_expression_data) <-
    colnames(expression_data)
  
  rownames(new_expression_data) <-
    rownames(expression_data)
  new_expression_data
}


#######adjust BMI, sex, and IRIS, ethnicity
expression_data_2 <-
  lm_adjust(expression_data = expression_data,
            sample_info = object@sample_info,
            threads = 16)

temp_object <- object

temp_object@expression_data <- expression_data_2

temp_object@expression_data

save(temp_object, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/meta_cor_residuals.RData")
# save(p_fc, file = "p_fc")

temp_object@sample_info

names(temp_object@sample_info)

load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/meta_cor_residuals.RData")

temp_object
# ##correlation
# ##cor
# 
# 

# cor.test(as.numeric(expression_data[1,]),temp_object@sample_info[,67], method = "pearson",conf.level = 0.95) -> cor_result
# 
# cor_result$conf.int  %>% as_tibble() %>% t() -> ci
# 
# data.frame(variable_id = rownames(temp_object)[1],
#            cor_p = cor_result$p.value,
#            spearman_cor = cor_result$estimate,
#             CI2.5=ci[,1],
#            CI92.5=ci[,2])


####1.进行pearson回归分析
# cor_data <-
#   lapply(temp_object@sample_info[,c(67:71)], function(y){
#     lapply(seq_along(rownames(temp_object)), function(i){
#       # cat(name, " ")
#     value <- as.numeric(expression_data[i,])
#     cor_result <-
#       cor.test(value, y, method = "pearson")
#     cor_result$conf.int  %>% as_tibble() %>% t() -> ci
#     data.frame(variable_id = rownames(temp_object)[i],
#                cor_p = cor_result$p.value,
#                pearson_cor = cor_result$estimate,
#                CI2.5=ci[,1],
#                CI97.5=ci[,2])
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#   })

as.numeric(expression_data[1,])

length(as.numeric(expression_data[1,]))

length(temp_object@sample_info[,c(24)])

cor.test(as.numeric(expression_data[1,]), temp_object@sample_info[,c(24)], method = "spearman")


####>>>>>>>>>>>>>>2.对精子端粒进行spearman回归分析<<<<<<<<<<<<<<<<<<<<<<<
cor_data <-
    lapply(seq_along(rownames(temp_object)), function(i){
      value <- as.numeric(expression_data[i,]) 
      cor_result <- cor.test(value,temp_object@sample_info[,c(24)], method = "spearman")

      data.frame(variable_id = rownames(temp_object)[i],
                 cor_p = cor_result$p.value,
                 spearman_cor = cor_result$estimate)
    }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



names(temp_object@sample_info)



cor_data %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_meta_stl_cor.xlsx")



###>>>>>>>>>>>>>>>>>>>3.体积的回归系数及校正p值<<<<<<<<<<<<<<<<<<<<<<===
##### lg_volume  lg_density lg_sperm_count lg_total_motility lg_progressive_motility
##### 
# ###permutation to get the p value   根据组合生成p值
# temp_object@sample_info[,c(67:71)
# dir.create("permutaton_cor_data")
for (idx in 1:271) {
  cat(idx, " ")
  permutation_cor_data <-
    seq_len(nrow(temp_object)) %>%
    purrr::map(function(i) {
      if((i %% 271) == 0)
        cat(i, " ")
      value <-
        as.numeric(unlist(temp_object[i, , drop = TRUE]))
      cor_result <-
        cor.test(value, sample(temp_object@sample_info$STL,
                               replace = FALSE), method = "spearman")
      data.frame(
        variable_id = rownames(temp_object)[i],
        cor_p = cor_result$p.value,
        spearman_cor = cor_result$estimate
      )
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  save(permutation_cor_data,
       file = file.path("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/",
                        paste0("permutation_cor_data_", idx)), compress = "xz")
  
}

load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/permutation_cor_data_1")

permutaton_cor_data_all <-
  permutation_cor_data %>%
  dplyr::select(-cor_p)

for (idx in 2:271) {
  cat(idx, " ")
  load(file.path(
    "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/",
    paste0("permutation_cor_data_", idx)
  ))
  permutaton_cor_data_all <-
    cbind(permutaton_cor_data_all, permutation_cor_data[, 3, drop = FALSE])
}

rownames(permutaton_cor_data_all) <- NULL

permutaton_cor_data_all <-
  permutaton_cor_data_all %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "variable_id")



###提取出精液体系的相关系数


permutated_p_value <-
  1:nrow(cor_data) %>%
  purrr::map(function(idx) {
    cat(idx, " ")
    original_cor <-
      cor_data$spearman_cor[idx]
    permutation_cor <-
      sample(as.numeric(permutaton_cor_data_all[idx,]), 10000, replace = TRUE)
    if (original_cor > 0) {
      sum(permutation_cor > original_cor) / 10000
    } else{
      sum(permutation_cor < original_cor) / 10000
    }
    
  }) %>%
  unlist()

cor_data$permutated_p_value <-
  permutated_p_value

cor_data

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p = cor_data$cor_p,
    permutated_p_value = cor_data$permutated_p_value,
    spearman_cor = cor_data$spearman_cor
  )

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p_adjust = p.adjust(cor_p, method = "fdr"),
    permutated_p_adjust = p.adjust(permutated_p_value, method = "fdr")
  )

temp_object@variable_info

plot(cor_data$cor_p,
     cor_data$permutated_p_value)

sum(temp_object@variable_info$cor_p < 0.05)
sum(temp_object@variable_info$cor_p_adjust < 0.05)
sum(temp_object@variable_info$permutated_p_value < 0.05)
sum(temp_object@variable_info$permutated_p_adjust < 0.05)


temp_object@sample_info



######linear mixed model

# 
linear_model_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    lm_result <-
      glm(formula = STL ~ value, data = temp_data)
    
    lm_result <-
      lm_result %>%
      broom::tidy()
    data.frame(lm_p = lm_result$p.value[2],
               coefficient = lm_result$estimate[2])
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()


save(linear_model_data, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/linear_model_data")



temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p = linear_model_data$lm_p,
                coefficient = linear_model_data$coefficient)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p_adjust = p.adjust(lm_p, method = "fdr"))

temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  ggplot(aes(x = spearman_cor, coefficient)) +
  geom_point() +
  # theme_base +
  labs(x = "Spearman correlation",
       y = "lm beta")

temp_object@variable_info

write.xlsx(temp_object@variable_info,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/cor_meta_stl.xlsx")

# load("result/20241222_PFAS_META_sperm/data_preparation/cor_meta_volume.xlsx")
###volcano plot  火山图
variable_info <-
  extract_variable_info(temp_object)

top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(10) 
# dplyr::pull(Compound.name)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) 
# dplyr::pull(Compound.name)

volcano_plot <-
  variable_info %>%
  mutate(marker = case_when(
    cor_p_adjust < 0.2 & spearman_cor > 0 ~ "Up",
    cor_p_adjust < 0.2 &
      spearman_cor < 0 ~ "Down",
    TRUE ~ "No"
  )) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)") 
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_up_marker_name,
#                  Compound.name, NA)
# ), size = 3) +
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_down_marker_name,
#                  Compound.name, NA)
# ), size = 3)

volcano_plot



####过滤出p<0.05或校正后的p<0.2，过滤出与精液质量相关的代谢物质
variable_info

volume_markers <-
  variable_info %>%
  dplyr::filter(cor_p < 0.05)

temp_object@variable_info$variable_id
###进行主成分分析

pca_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% volume_markers$variable_id) %>%
  massstat::run_pca()

temp_object@sample_info$lg_volume

###画出主成分的散点图
plot <-
  massstat::pca_score_plot(
    object = temp_object,
    pca_object = pca_object,
    color_by = "STL",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
# scale_x_continuous(limits = c(-0.2, 0.2)) +
# scale_y_continuous(limits = c(-0.2, 0.2))

plot

subject_info <-
  sample_info %>%
  dplyr::select(sample_id, STL) %>%
  dplyr::filter(!is.na(STL)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

###提取出第一主成分


lapply(1:72, function(i){
data.frame(
  subject_id = rownames(pca_object$x),
  PC = apply(pca_object$x[,i, drop = FALSE],1, mean),
  class = i
) %>%
  as_tibble()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()




data.frame(
  subject_id = rownames(pca_object$x),
  PC = apply(pca_object$x[, 1, drop = FALSE], 1, mean),
  class = "STL"
) %>%
  as_tibble() -> df_pca


df_pca$subject_id <- as.numeric(df_pca$subject_id)

df_pca %>% 
  dplyr::left_join(subject_info,by = c("subject_id"="sample_id"))  %>% 
  dplyr::filter(!is.na(STL)) -> temp_data 

write.xlsx(temp_data,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/temp_data_STL.xlsx")

###第一组组成分与精液质量的回归系数及直线图


omics_color <-
  c(
    "STL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "Concentration" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "Sperm count" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "Total motility" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "Progressive motility" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5]
  )


temp_data %>%
  ggplot(aes(PC,STL)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  # theme_base +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  # facet_wrap(facets = vars(class),
  #            scales = "free",
  #            nrow = 2) +
  scale_color_manual(values = omics_color) +
  labs(x = "PC1", y = "STL")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        legend.position = c("none"),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 1, colour = "black")  # 设置边框线的粗细和颜色
  )+
  annotate("text", 
           x = Inf, y = Inf, label = "Abs(Spearman cor) = 0.11; P Value < 0.01", 
           hjust =2, vjust = 4, size = 5, color = "black", fontface = "bold") -> p2

p2



cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/pfas_meta_STL_merger_line.pdf",width = 12, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p2   # 使用大写罗马表示
par(opar)
dev.off()


temp_data %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(PC, STL, method = "spearman"))



cor.test(temp_data$PC, temp_data$STL,  method = "spearman")

####绝对值回归系数与p值
library(plyr)
 temp_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PC, x$STL,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()






###############################################################################
####################对PFAS与精浆代谢组进行spearman相关分析#####################
##############################################################################
####>>>>>>>>>>>>>>2.对PFOS进行spearman回归分析<<<<<<<<<<<<<<<<<<<<<<<
temp_object@sample_info

names(temp_object@sample_info)

cor_data <-
  lapply(seq_along(rownames(temp_object)), function(i){
    value <- as.numeric(expression_data[i,]) 
    cor_result <- cor.test(value,temp_object@sample_info[,c(11)], method = "spearman")
    
    data.frame(variable_id = rownames(temp_object)[i],
               cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



names(temp_object@sample_info)



cor_data %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/result_PFOS_meta_cor.xlsx")

###>>>>>>>>>>>>>>>>>>>回归系数及校正p值<<<<<<<<<<<<<<<<<<<<<<===
##### lg_volume  lg_density lg_sperm_count lg_total_motility lg_progressive_motility
##### 
# ###permutation to get the p value   根据组合生成p值
# temp_object@sample_info[,c(67:71)
# dir.create("permutaton_cor_data")
for (idx in 1:271) {
  cat(idx, " ")
  permutation_cor_data <-
    seq_len(nrow(temp_object)) %>%
    purrr::map(function(i) {
      if((i %% 271) == 0)
        cat(i, " ")
      value <-
        as.numeric(unlist(temp_object[i, , drop = TRUE]))
      cor_result <-
        cor.test(value, sample(temp_object@sample_info$PFOS,
                               replace = FALSE), method = "spearman")
      data.frame(
        variable_id = rownames(temp_object)[i],
        cor_p = cor_result$p.value,
        spearman_cor = cor_result$estimate
      )
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  save(permutation_cor_data,
       file = file.path("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
                        paste0("permutation_cor_data_", idx)), compress = "xz")
  
}

load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/permutation_cor_data_1")

permutaton_cor_data_all <-
  permutation_cor_data %>%
  dplyr::select(-cor_p)

for (idx in 2:271) {
  cat(idx, " ")
  load(file.path(
    "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
    paste0("permutation_cor_data_", idx)
  ))
  permutaton_cor_data_all <-
    cbind(permutaton_cor_data_all, permutation_cor_data[, 3, drop = FALSE])
}

rownames(permutaton_cor_data_all) <- NULL

permutaton_cor_data_all <-
  permutaton_cor_data_all %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "variable_id")



###提取出精液体系的相关系数


permutated_p_value <-
  1:nrow(cor_data) %>%
  purrr::map(function(idx) {
    cat(idx, " ")
    original_cor <-
      cor_data$spearman_cor[idx]
    permutation_cor <-
      sample(as.numeric(permutaton_cor_data_all[idx,]), 10000, replace = TRUE)
    if (original_cor > 0) {
      sum(permutation_cor > original_cor) / 10000
    } else{
      sum(permutation_cor < original_cor) / 10000
    }
    
  }) %>%
  unlist()

cor_data$permutated_p_value <-
  permutated_p_value

cor_data

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p = cor_data$cor_p,
    permutated_p_value = cor_data$permutated_p_value,
    spearman_cor = cor_data$spearman_cor
  )

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p_adjust = p.adjust(cor_p, method = "fdr"),
    permutated_p_adjust = p.adjust(permutated_p_value, method = "fdr")
  )

temp_object@variable_info

plot(cor_data$cor_p,
     cor_data$permutated_p_value)


temp_object@sample_info



######linear mixed model

# 
linear_model_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    lm_result <-
      glm(formula = PFOS ~ value, data = temp_data)
    
    lm_result <-
      lm_result %>%
      broom::tidy()
    data.frame(lm_p = lm_result$p.value[2],
               coefficient = lm_result$estimate[2])
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()


save(linear_model_data, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/linear_model_data")



temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p = linear_model_data$lm_p,
                coefficient = linear_model_data$coefficient)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p_adjust = p.adjust(lm_p, method = "fdr"))

temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  ggplot(aes(x = spearman_cor, coefficient)) +
  geom_point() +
  # theme_base +
  labs(x = "Spearman correlation",
       y = "lm beta")

temp_object@variable_info



write.xlsx(temp_object@variable_info,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PFOS_meta.xlsx")

# load("result/20241222_PFAS_META_sperm/data_preparation/cor_meta_volume.xlsx")
###volcano plot  火山图
variable_info <-
  extract_variable_info(temp_object)

top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(10) 
# dplyr::pull(Compound.name)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) 
# dplyr::pull(Compound.name)

volcano_plot <-
  variable_info %>%
  mutate(marker = case_when(
    cor_p_adjust < 0.2 & spearman_cor > 0 ~ "Up",
    cor_p_adjust < 0.2 &
      spearman_cor < 0 ~ "Down",
    TRUE ~ "No"
  )) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)") 
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_up_marker_name,
#                  Compound.name, NA)
# ), size = 3) +
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_down_marker_name,
#                  Compound.name, NA)
# ), size = 3)

volcano_plot

sum(temp_object@variable_info$cor_p < 0.05)
sum(temp_object@variable_info$cor_p_adjust < 0.05)
sum(temp_object@variable_info$permutated_p_value < 0.05)
sum(temp_object@variable_info$permutated_p_adjust < 0.05)
sum(temp_object@variable_info$lm_p < 0.05)
sum(temp_object@variable_info$lm_p_adjust < 0.05)

####过滤出p<0.05或校正后的p<0.2，过滤出与精液质量相关的代谢物质
variable_info

PFAS_markers <-
  variable_info %>%
  dplyr::filter(cor_p < 0.05)

temp_object@variable_info$variable_id
###进行主成分分析

pca_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% PFAS_markers$variable_id) %>%
  massstat::run_pca()

temp_object@sample_info

###画出主成分的散点图
plot <-
  massstat::pca_score_plot(
    object = temp_object,
    pca_object = pca_object,
    color_by = "PFOS",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
# scale_x_continuous(limits = c(-0.2, 0.2)) +
# scale_y_continuous(limits = c(-0.2, 0.2))

plot

subject_info <-
  sample_info %>%
  dplyr::select(sample_id, PFOS) %>%
  dplyr::filter(!is.na(PFOS)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)



###提取出所有的主成分
lapply(1:135, function(i){
  data.frame(
    subject_id = rownames(pca_object$x),
    PC = apply(pca_object$x[,i, drop = FALSE],1, mean),
    class = i
  ) %>%
    as_tibble()
}) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



###提取出第一主成分
data.frame(
  subject_id = rownames(pca_object$x),
  PC = apply(pca_object$x[, 1, drop = FALSE], 1, mean),
  class = "PFOS"
) %>%
  as_tibble() -> df_pca


df_pca$subject_id <- as.numeric(df_pca$subject_id)

df_pca %>% 
  dplyr::left_join(subject_info,by = c("subject_id"="sample_id"))  %>% 
  dplyr::filter(!is.na(PFOS)) -> temp_data 

write.xlsx(temp_data,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PFOS.xlsx")

###第一组组成分与精液质量的回归系数及直线图

names(temp_object@sample_info)

omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "L_PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "PF3ONS_9CL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6]
  )


temp_data %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(PFOS, PC, method = "spearman"))



cor.test(temp_data$PFOS,temp_data$PC,   method = "spearman")

####绝对值回归系数与p值
library(plyr)
temp_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PFOS, x$PC,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()



temp_data %>%
  ggplot(aes(PFOS,PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  # theme_base +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  # facet_wrap(facets = vars(class),
  #            scales = "free",
  #            nrow = 2) +
  scale_color_manual(values = omics_color) +
  labs(x = "PFOS", y = "PC1")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        legend.position = c("none"),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 1, colour = "black")  # 设置边框线的粗细和颜色
  )+
  annotate("text", 
           x = Inf, y = Inf, label = "Abs(Spearman cor) = 0.23; P Value < 0.001", 
           hjust =1.8, vjust = 4, size = 5, color = "black", fontface = "bold") -> p2

p2



cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/pfas_meta_PFOS_merger_line.pdf",width = 12, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p2   # 使用大写罗马表示
par(opar)
dev.off()



temp_data

temp_data %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(PC, PFOS, method = "spearman"))



cor.test(temp_data$PC, temp_data$PFOS,  method = "spearman")

####绝对值回归系数与p值
library(plyr)
temp_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PC, x$PFOS,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()









####>>>>>>>>>>>>>>2.对PFOA进行spearman回归分析<<<<<<<<<<<<<<<<<<<<<<<
temp_object@sample_info

names(temp_object@sample_info)

cor_data <-
  lapply(seq_along(rownames(temp_object)), function(i){
    value <- as.numeric(expression_data[i,]) 
    cor_result <- cor.test(value,temp_object@sample_info[,c(12)], method = "spearman")
    
    data.frame(variable_id = rownames(temp_object)[i],
               cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



names(temp_object@sample_info)



cor_data %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/result_PFOA_meta_cor.xlsx")

###>>>>>>>>>>>>>>>>>>>回归系数及校正p值<<<<<<<<<<<<<<<<<<<<<<===
##### lg_volume  lg_density lg_sperm_count lg_total_motility lg_progressive_motility
##### 
# ###permutation to get the p value   根据组合生成p值
# temp_object@sample_info[,c(67:71)
# dir.create("permutaton_cor_data")
for (idx in 1:271) {
  cat(idx, " ")
  permutation_cor_data <-
    seq_len(nrow(temp_object)) %>%
    purrr::map(function(i) {
      if((i %% 271) == 0)
        cat(i, " ")
      value <-
        as.numeric(unlist(temp_object[i, , drop = TRUE]))
      cor_result <-
        cor.test(value, sample(temp_object@sample_info$PFOA,
                               replace = FALSE), method = "spearman")
      data.frame(
        variable_id = rownames(temp_object)[i],
        cor_p = cor_result$p.value,
        spearman_cor = cor_result$estimate
      )
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  save(permutation_cor_data,
       file = file.path("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
                        paste0("permutation_cor_data_", idx)), compress = "xz")
  
}

load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/permutation_cor_data_1")

permutaton_cor_data_all <-
  permutation_cor_data %>%
  dplyr::select(-cor_p)

for (idx in 2:271) {
  cat(idx, " ")
  load(file.path(
    "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
    paste0("permutation_cor_data_", idx)
  ))
  permutaton_cor_data_all <-
    cbind(permutaton_cor_data_all, permutation_cor_data[, 3, drop = FALSE])
}

rownames(permutaton_cor_data_all) <- NULL

permutaton_cor_data_all <-
  permutaton_cor_data_all %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "variable_id")



###提取出精液体系的相关系数


permutated_p_value <-
  1:nrow(cor_data) %>%
  purrr::map(function(idx) {
    cat(idx, " ")
    original_cor <-
      cor_data$spearman_cor[idx]
    permutation_cor <-
      sample(as.numeric(permutaton_cor_data_all[idx,]), 10000, replace = TRUE)
    if (original_cor > 0) {
      sum(permutation_cor > original_cor) / 10000
    } else{
      sum(permutation_cor < original_cor) / 10000
    }
    
  }) %>%
  unlist()

cor_data$permutated_p_value <-
  permutated_p_value

cor_data

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p = cor_data$cor_p,
    permutated_p_value = cor_data$permutated_p_value,
    spearman_cor = cor_data$spearman_cor
  )

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p_adjust = p.adjust(cor_p, method = "fdr"),
    permutated_p_adjust = p.adjust(permutated_p_value, method = "fdr")
  )

temp_object@variable_info

plot(cor_data$cor_p,
     cor_data$permutated_p_value)


temp_object@sample_info



######linear mixed model

# 
linear_model_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    lm_result <-
      glm(formula = PFOA ~ value, data = temp_data)
    
    lm_result <-
      lm_result %>%
      broom::tidy()
    data.frame(lm_p = lm_result$p.value[2],
               coefficient = lm_result$estimate[2])
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()


save(linear_model_data, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/linear_model_data")



temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p = linear_model_data$lm_p,
                coefficient = linear_model_data$coefficient)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p_adjust = p.adjust(lm_p, method = "fdr"))

temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  ggplot(aes(x = spearman_cor, coefficient)) +
  geom_point() +
  # theme_base +
  labs(x = "Spearman correlation",
       y = "lm beta")

temp_object@variable_info



write.xlsx(temp_object@variable_info,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PFOA_meta.xlsx")

# load("result/20241222_PFAS_META_sperm/data_preparation/cor_meta_volume.xlsx")
###volcano plot  火山图
variable_info <-
  extract_variable_info(temp_object)

top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(10) 
# dplyr::pull(Compound.name)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) 
# dplyr::pull(Compound.name)

volcano_plot <-
  variable_info %>%
  mutate(marker = case_when(
    cor_p_adjust < 0.2 & spearman_cor > 0 ~ "Up",
    cor_p_adjust < 0.2 &
      spearman_cor < 0 ~ "Down",
    TRUE ~ "No"
  )) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)") 
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_up_marker_name,
#                  Compound.name, NA)
# ), size = 3) +
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_down_marker_name,
#                  Compound.name, NA)
# ), size = 3)

volcano_plot

sum(temp_object@variable_info$cor_p < 0.05)
sum(temp_object@variable_info$cor_p_adjust < 0.05)
sum(temp_object@variable_info$permutated_p_value < 0.05)
sum(temp_object@variable_info$permutated_p_adjust < 0.05)
sum(temp_object@variable_info$lm_p < 0.05)
sum(temp_object@variable_info$lm_p_adjust < 0.05)

####过滤出p<0.05或校正后的p<0.2，过滤出与精液质量相关的代谢物质
variable_info

PFAS_markers <-
  variable_info %>%
  dplyr::filter(cor_p < 0.05)

temp_object@variable_info$variable_id
###进行主成分分析

pca_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% PFAS_markers$variable_id) %>%
  massstat::run_pca()

temp_object@sample_info

###画出主成分的散点图
plot <-
  massstat::pca_score_plot(
    object = temp_object,
    pca_object = pca_object,
    color_by = "PFOS",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
# scale_x_continuous(limits = c(-0.2, 0.2)) +
# scale_y_continuous(limits = c(-0.2, 0.2))

plot

subject_info <-
  sample_info %>%
  dplyr::select(sample_id, PFOA) %>%
  dplyr::filter(!is.na(PFOA)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)



###提取出所有的主成分
lapply(1:110, function(i){
  data.frame(
    subject_id = rownames(pca_object$x),
    PC = apply(pca_object$x[,i, drop = FALSE],1, mean),
    class = i
  ) %>%
    as_tibble()
}) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



###提取出第一主成分
data.frame(
  subject_id = rownames(pca_object$x),
  PC = apply(pca_object$x[, 1, drop = FALSE], 1, mean),
  class = "PFOA"
) %>%
  as_tibble() -> df_pca


df_pca$subject_id <- as.numeric(df_pca$subject_id)

df_pca %>% 
  dplyr::left_join(subject_info,by = c("subject_id"="sample_id"))  %>% 
  dplyr::filter(!is.na(PFOA)) -> temp_data 

write.xlsx(temp_data,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PFOA.xlsx")

###第一组组成分与精液质量的回归系数及直线图

names(temp_object@sample_info)

omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "L_PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "PF3ONS_9CL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6]
  )


temp_data %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(PFOA, PC, method = "spearman"))



cor.test(temp_data$PFOA,temp_data$PC,   method = "spearman")



####绝对值回归系数与p值
library(plyr)
temp_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PFOA, x$PC,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()



temp_data %>%
  ggplot(aes(PFOA,PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  # theme_base +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  # facet_wrap(facets = vars(class),
  #            scales = "free",
  #            nrow = 2) +
  scale_color_manual(values = omics_color) +
  labs(x = "PFOA", y = "PC1")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        legend.position = c("none"),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 1, colour = "black")  # 设置边框线的粗细和颜色
  )+
  annotate("text", 
           x = Inf, y = Inf, label = "Abs(Spearman cor) = 0.20; P Value < 0.001", 
           hjust =1.8, vjust = 4, size = 5, color = "black", fontface = "bold") -> p4

p4



cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/pfas_meta_PFOA_merger_line.pdf",width = 12, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p4   # 使用大写罗马表示
par(opar)
dev.off()





####>>>>>>>>>>>>>>4.对PFDA进行spearman回归分析<<<<<<<<<<<<<<<<<<<<<<<
temp_object@sample_info

names(temp_object@sample_info)

cor_data <-
  lapply(seq_along(rownames(temp_object)), function(i){
    value <- as.numeric(expression_data[i,]) 
    cor_result <- cor.test(value,temp_object@sample_info[,c(13)], method = "spearman")
    
    data.frame(variable_id = rownames(temp_object)[i],
               cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



names(temp_object@sample_info)



cor_data %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/result_PFDAPFOA_meta_cor.xlsx")

###>>>>>>>>>>>>>>>>>>>回归系数及校正p值<<<<<<<<<<<<<<<<<<<<<<===
##### lg_volume  lg_density lg_sperm_count lg_total_motility lg_progressive_motility
##### 
# ###permutation to get the p value   根据组合生成p值
# temp_object@sample_info[,c(67:71)
# dir.create("permutaton_cor_data")
for (idx in 1:271) {
  cat(idx, " ")
  permutation_cor_data <-
    seq_len(nrow(temp_object)) %>%
    purrr::map(function(i) {
      if((i %% 271) == 0)
        cat(i, " ")
      value <-
        as.numeric(unlist(temp_object[i, , drop = TRUE]))
      cor_result <-
        cor.test(value, sample(temp_object@sample_info$PFDA,
                               replace = FALSE), method = "spearman")
      data.frame(
        variable_id = rownames(temp_object)[i],
        cor_p = cor_result$p.value,
        spearman_cor = cor_result$estimate
      )
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  save(permutation_cor_data,
       file = file.path("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
                        paste0("permutation_cor_data_", idx)), compress = "xz")
  
}

load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/permutation_cor_data_1")

permutaton_cor_data_all <-
  permutation_cor_data %>%
  dplyr::select(-cor_p)

for (idx in 2:271) {
  cat(idx, " ")
  load(file.path(
    "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
    paste0("permutation_cor_data_", idx)
  ))
  permutaton_cor_data_all <-
    cbind(permutaton_cor_data_all, permutation_cor_data[, 3, drop = FALSE])
}

rownames(permutaton_cor_data_all) <- NULL

permutaton_cor_data_all <-
  permutaton_cor_data_all %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "variable_id")



###提取出精液体系的相关系数


permutated_p_value <-
  1:nrow(cor_data) %>%
  purrr::map(function(idx) {
    cat(idx, " ")
    original_cor <-
      cor_data$spearman_cor[idx]
    permutation_cor <-
      sample(as.numeric(permutaton_cor_data_all[idx,]), 10000, replace = TRUE)
    if (original_cor > 0) {
      sum(permutation_cor > original_cor) / 10000
    } else{
      sum(permutation_cor < original_cor) / 10000
    }
    
  }) %>%
  unlist()

cor_data$permutated_p_value <-
  permutated_p_value

cor_data

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p = cor_data$cor_p,
    permutated_p_value = cor_data$permutated_p_value,
    spearman_cor = cor_data$spearman_cor
  )

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p_adjust = p.adjust(cor_p, method = "fdr"),
    permutated_p_adjust = p.adjust(permutated_p_value, method = "fdr")
  )

temp_object@variable_info

plot(cor_data$cor_p,
     cor_data$permutated_p_value)


temp_object@sample_info



######linear mixed model

# 
linear_model_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    lm_result <-
      glm(formula = PFDA ~ value, data = temp_data)
    
    lm_result <-
      lm_result %>%
      broom::tidy()
    data.frame(lm_p = lm_result$p.value[2],
               coefficient = lm_result$estimate[2])
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()


save(linear_model_data, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/linear_model_data")



temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p = linear_model_data$lm_p,
                coefficient = linear_model_data$coefficient)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p_adjust = p.adjust(lm_p, method = "fdr"))

temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  ggplot(aes(x = spearman_cor, coefficient)) +
  geom_point() +
  # theme_base +
  labs(x = "Spearman correlation",
       y = "lm beta")

temp_object@variable_info



write.xlsx(temp_object@variable_info,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PFDA_meta.xlsx")

# load("result/20241222_PFAS_META_sperm/data_preparation/cor_meta_volume.xlsx")
###volcano plot  火山图
variable_info <-
  extract_variable_info(temp_object)

top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(10) 
# dplyr::pull(Compound.name)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) 
# dplyr::pull(Compound.name)

volcano_plot <-
  variable_info %>%
  mutate(marker = case_when(
    cor_p_adjust < 0.2 & spearman_cor > 0 ~ "Up",
    cor_p_adjust < 0.2 &
      spearman_cor < 0 ~ "Down",
    TRUE ~ "No"
  )) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)") 
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_up_marker_name,
#                  Compound.name, NA)
# ), size = 3) +
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_down_marker_name,
#                  Compound.name, NA)
# ), size = 3)

volcano_plot

sum(temp_object@variable_info$cor_p < 0.05)
sum(temp_object@variable_info$cor_p_adjust < 0.05)
sum(temp_object@variable_info$permutated_p_value < 0.05)
sum(temp_object@variable_info$permutated_p_adjust < 0.05)
sum(temp_object@variable_info$lm_p < 0.05)
sum(temp_object@variable_info$lm_p_adjust < 0.05)

####过滤出p<0.05或校正后的p<0.2，过滤出与精液质量相关的代谢物质
variable_info

PFAS_markers <-
  variable_info %>%
  dplyr::filter(cor_p < 0.05)

temp_object@variable_info$variable_id
###进行主成分分析

pca_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% PFAS_markers$variable_id) %>%
  massstat::run_pca()

temp_object@sample_info

###画出主成分的散点图
plot <-
  massstat::pca_score_plot(
    object = temp_object,
    pca_object = pca_object,
    color_by = "PFDA",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
# scale_x_continuous(limits = c(-0.2, 0.2)) +
# scale_y_continuous(limits = c(-0.2, 0.2))

plot

subject_info <-
  sample_info %>%
  dplyr::select(sample_id, PFDA) %>%
  dplyr::filter(!is.na(PFDA)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)



###提取出所有的主成分
lapply(1:106, function(i){
  data.frame(
    subject_id = rownames(pca_object$x),
    PC = apply(pca_object$x[,i, drop = FALSE],1, mean),
    class = i
  ) %>%
    as_tibble()
}) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



###提取出第一主成分
data.frame(
  subject_id = rownames(pca_object$x),
  PC = apply(pca_object$x[, 1, drop = FALSE], 1, mean),
  class = "PFDA"
) %>%
  as_tibble() -> df_pca


df_pca$subject_id <- as.numeric(df_pca$subject_id)

df_pca %>% 
  dplyr::left_join(subject_info,by = c("subject_id"="sample_id"))  %>% 
  dplyr::filter(!is.na(PFDA)) -> temp_data 

write.xlsx(temp_data,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PFDA.xlsx")

###第一组组成分与精液质量的回归系数及直线图

names(temp_object@sample_info)

omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "L_PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "PF3ONS_9CL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6]
  )


temp_data %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(PFDA, PC, method = "spearman"))



cor.test(temp_data$PFDA,temp_data$PC,   method = "spearman")



####绝对值回归系数与p值
library(plyr)
temp_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PFDA, x$PC,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()



temp_data %>%
  ggplot(aes(PFDA,PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  # theme_base +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  # facet_wrap(facets = vars(class),
  #            scales = "free",
  #            nrow = 2) +
  scale_color_manual(values = omics_color) +
  labs(x = "PFDA", y = "PC1")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        legend.position = c("none"),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 1, colour = "black")  # 设置边框线的粗细和颜色
  )+
  annotate("text", 
           x = Inf, y = Inf, label = "Abs(Spearman cor) = 0.20; P Value < 0.001", 
           hjust =1.8, vjust = 4, size = 5, color = "black", fontface = "bold") -> p5

p5



cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/pfas_meta_PFDA_merger_line.pdf",width = 12, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p5   # 使用大写罗马表示
par(opar)
dev.off()






####>>>>>>>>>>>>>>4.对PFUdA进行spearman回归分析<<<<<<<<<<<<<<<<<<<<<<<
temp_object@sample_info

names(temp_object@sample_info)

cor_data <-
  lapply(seq_along(rownames(temp_object)), function(i){
    value <- as.numeric(expression_data[i,]) 
    cor_result <- cor.test(value,temp_object@sample_info[,c(14)], method = "spearman")
    
    data.frame(variable_id = rownames(temp_object)[i],
               cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



names(temp_object@sample_info)



cor_data %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/result_PFUdA_meta_cor.xlsx")

###>>>>>>>>>>>>>>>>>>>回归系数及校正p值<<<<<<<<<<<<<<<<<<<<<<===
##### lg_volume  lg_density lg_sperm_count lg_total_motility lg_progressive_motility
##### 
# ###permutation to get the p value   根据组合生成p值
# temp_object@sample_info[,c(67:71)
# dir.create("permutaton_cor_data")
for (idx in 1:271) {
  cat(idx, " ")
  permutation_cor_data <-
    seq_len(nrow(temp_object)) %>%
    purrr::map(function(i) {
      if((i %% 271) == 0)
        cat(i, " ")
      value <-
        as.numeric(unlist(temp_object[i, , drop = TRUE]))
      cor_result <-
        cor.test(value, sample(temp_object@sample_info$PFUdA,
                               replace = FALSE), method = "spearman")
      data.frame(
        variable_id = rownames(temp_object)[i],
        cor_p = cor_result$p.value,
        spearman_cor = cor_result$estimate
      )
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  save(permutation_cor_data,
       file = file.path("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
                        paste0("permutation_cor_data_", idx)), compress = "xz")
  
}

load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/permutation_cor_data_1")

permutaton_cor_data_all <-
  permutation_cor_data %>%
  dplyr::select(-cor_p)

for (idx in 2:271) {
  cat(idx, " ")
  load(file.path(
    "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
    paste0("permutation_cor_data_", idx)
  ))
  permutaton_cor_data_all <-
    cbind(permutaton_cor_data_all, permutation_cor_data[, 3, drop = FALSE])
}

rownames(permutaton_cor_data_all) <- NULL

permutaton_cor_data_all <-
  permutaton_cor_data_all %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "variable_id")



###提取出精液体系的相关系数


permutated_p_value <-
  1:nrow(cor_data) %>%
  purrr::map(function(idx) {
    cat(idx, " ")
    original_cor <-
      cor_data$spearman_cor[idx]
    permutation_cor <-
      sample(as.numeric(permutaton_cor_data_all[idx,]), 10000, replace = TRUE)
    if (original_cor > 0) {
      sum(permutation_cor > original_cor) / 10000
    } else{
      sum(permutation_cor < original_cor) / 10000
    }
    
  }) %>%
  unlist()

cor_data$permutated_p_value <-
  permutated_p_value

cor_data

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p = cor_data$cor_p,
    permutated_p_value = cor_data$permutated_p_value,
    spearman_cor = cor_data$spearman_cor
  )

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p_adjust = p.adjust(cor_p, method = "fdr"),
    permutated_p_adjust = p.adjust(permutated_p_value, method = "fdr")
  )

temp_object@variable_info

plot(cor_data$cor_p,
     cor_data$permutated_p_value)


temp_object@sample_info



######linear mixed model

# 
linear_model_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    lm_result <-
      glm(formula = PFUdA ~ value, data = temp_data)
    
    lm_result <-
      lm_result %>%
      broom::tidy()
    data.frame(lm_p = lm_result$p.value[2],
               coefficient = lm_result$estimate[2])
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()


save(linear_model_data, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/linear_model_data")



temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p = linear_model_data$lm_p,
                coefficient = linear_model_data$coefficient)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p_adjust = p.adjust(lm_p, method = "fdr"))

temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  ggplot(aes(x = spearman_cor, coefficient)) +
  geom_point() +
  # theme_base +
  labs(x = "Spearman correlation",
       y = "lm beta")

temp_object@variable_info



write.xlsx(temp_object@variable_info,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PFUdA_meta.xlsx")

# load("result/20241222_PFAS_META_sperm/data_preparation/cor_meta_volume.xlsx")
###volcano plot  火山图
variable_info <-
  extract_variable_info(temp_object)

top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(10) 
# dplyr::pull(Compound.name)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) 
# dplyr::pull(Compound.name)

volcano_plot <-
  variable_info %>%
  mutate(marker = case_when(
    cor_p_adjust < 0.2 & spearman_cor > 0 ~ "Up",
    cor_p_adjust < 0.2 &
      spearman_cor < 0 ~ "Down",
    TRUE ~ "No"
  )) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)") 
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_up_marker_name,
#                  Compound.name, NA)
# ), size = 3) +
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_down_marker_name,
#                  Compound.name, NA)
# ), size = 3)

volcano_plot

sum(temp_object@variable_info$cor_p < 0.05)
sum(temp_object@variable_info$cor_p_adjust < 0.05)
sum(temp_object@variable_info$permutated_p_value < 0.05)
sum(temp_object@variable_info$permutated_p_adjust < 0.05)
sum(temp_object@variable_info$lm_p < 0.05)
sum(temp_object@variable_info$lm_p_adjust < 0.05)

####过滤出p<0.05或校正后的p<0.2，过滤出与精液质量相关的代谢物质
variable_info

PFAS_markers <-
  variable_info %>%
  dplyr::filter(cor_p < 0.05)

temp_object@variable_info$variable_id
###进行主成分分析

pca_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% PFAS_markers$variable_id) %>%
  massstat::run_pca()

temp_object@sample_info

###画出主成分的散点图
plot <-
  massstat::pca_score_plot(
    object = temp_object,
    pca_object = pca_object,
    color_by = "PFUdA",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
# scale_x_continuous(limits = c(-0.2, 0.2)) +
# scale_y_continuous(limits = c(-0.2, 0.2))

plot

subject_info <-
  sample_info %>%
  dplyr::select(sample_id, PFUdA) %>%
  dplyr::filter(!is.na(PFUdA)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)



###提取出所有的主成分
lapply(1:149, function(i){
  data.frame(
    subject_id = rownames(pca_object$x),
    PC = apply(pca_object$x[,i, drop = FALSE],1, mean),
    class = i
  ) %>%
    as_tibble()
}) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



###提取出第一主成分
data.frame(
  subject_id = rownames(pca_object$x),
  PC = apply(pca_object$x[, 1, drop = FALSE], 1, mean),
  class = "PFUdA"
) %>%
  as_tibble() -> df_pca


df_pca$subject_id <- as.numeric(df_pca$subject_id)

df_pca %>% 
  dplyr::left_join(subject_info,by = c("subject_id"="sample_id"))  %>% 
  dplyr::filter(!is.na(PFUdA)) -> temp_data 

write.xlsx(temp_data,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PFUdA.xlsx")

###第一组组成分与精液质量的回归系数及直线图

names(temp_object@sample_info)

omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "L_PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "PF3ONS_9CL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6]
  )


temp_data %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(PFUdA, PC, method = "spearman"))



cor.test(temp_data$PFUdA,temp_data$PC,   method = "spearman")



####绝对值回归系数与p值
library(plyr)
temp_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PFUdA, x$PC,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()



temp_data %>%
  ggplot(aes(PFUdA,PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  # theme_base +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  # facet_wrap(facets = vars(class),
  #            scales = "free",
  #            nrow = 2) +
  scale_color_manual(values = omics_color) +
  labs(x = "PFDA", y = "PC1")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        legend.position = c("none"),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 1, colour = "black")  # 设置边框线的粗细和颜色
  )+
  annotate("text", 
           x = Inf, y = Inf, label = "Abs(Spearman cor) = 0.23; P Value < 0.001", 
           hjust =1.8, vjust = 4, size = 5, color = "black", fontface = "bold") -> p6

p6



cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/pfas_meta_PFUdA_merger_line.pdf",width = 12, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p6   # 使用大写罗马表示
par(opar)
dev.off()







####>>>>>>>>>>>>>>4.对L_PFHxS进行spearman回归分析<<<<<<<<<<<<<<<<<<<<<<<
temp_object@sample_info

names(temp_object@sample_info)

cor_data <-
  lapply(seq_along(rownames(temp_object)), function(i){
    value <- as.numeric(expression_data[i,]) 
    cor_result <- cor.test(value,temp_object@sample_info[,c(15)], method = "spearman")
    
    data.frame(variable_id = rownames(temp_object)[i],
               cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



names(temp_object@sample_info)



cor_data %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/result_L_PFHxS_meta_cor.xlsx")

###>>>>>>>>>>>>>>>>>>>回归系数及校正p值<<<<<<<<<<<<<<<<<<<<<<===
##### lg_volume  lg_density lg_sperm_count lg_total_motility lg_progressive_motility
##### 
# ###permutation to get the p value   根据组合生成p值
# temp_object@sample_info[,c(67:71)
# dir.create("permutaton_cor_data")
for (idx in 1:271) {
  cat(idx, " ")
  permutation_cor_data <-
    seq_len(nrow(temp_object)) %>%
    purrr::map(function(i) {
      if((i %% 271) == 0)
        cat(i, " ")
      value <-
        as.numeric(unlist(temp_object[i, , drop = TRUE]))
      cor_result <-
        cor.test(value, sample(temp_object@sample_info$L_PFHxS,
                               replace = FALSE), method = "spearman")
      data.frame(
        variable_id = rownames(temp_object)[i],
        cor_p = cor_result$p.value,
        spearman_cor = cor_result$estimate
      )
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  save(permutation_cor_data,
       file = file.path("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",paste0("permutation_cor_data_", idx)), compress = "xz")
  
}

load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/permutation_cor_data_1")

permutaton_cor_data_all <-
  permutation_cor_data %>%
  dplyr::select(-cor_p)

for (idx in 2:271) {
  cat(idx, " ")
  load(file.path(
    "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
    paste0("permutation_cor_data_", idx)
  ))
  permutaton_cor_data_all <-
    cbind(permutaton_cor_data_all, permutation_cor_data[, 3, drop = FALSE])
}

rownames(permutaton_cor_data_all) <- NULL

permutaton_cor_data_all <-
  permutaton_cor_data_all %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "variable_id")



###提取出精液体系的相关系数


permutated_p_value <-
  1:nrow(cor_data) %>%
  purrr::map(function(idx) {
    cat(idx, " ")
    original_cor <-
      cor_data$spearman_cor[idx]
    permutation_cor <-
      sample(as.numeric(permutaton_cor_data_all[idx,]), 10000, replace = TRUE)
    if (original_cor > 0) {
      sum(permutation_cor > original_cor) / 10000
    } else{
      sum(permutation_cor < original_cor) / 10000
    }
    
  }) %>%
  unlist()

cor_data$permutated_p_value <-
  permutated_p_value

cor_data

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p = cor_data$cor_p,
    permutated_p_value = cor_data$permutated_p_value,
    spearman_cor = cor_data$spearman_cor
  )

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p_adjust = p.adjust(cor_p, method = "fdr"),
    permutated_p_adjust = p.adjust(permutated_p_value, method = "fdr")
  )

temp_object@variable_info

plot(cor_data$cor_p,
     cor_data$permutated_p_value)


temp_object@sample_info



######linear mixed model

# 
linear_model_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    lm_result <-
      glm(formula = L_PFHxS ~ value, data = temp_data)
    
    lm_result <-
      lm_result %>%
      broom::tidy()
    data.frame(lm_p = lm_result$p.value[2],
               coefficient = lm_result$estimate[2])
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()


save(linear_model_data, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/linear_model_data")



temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p = linear_model_data$lm_p,
                coefficient = linear_model_data$coefficient)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p_adjust = p.adjust(lm_p, method = "fdr"))

temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  ggplot(aes(x = spearman_cor, coefficient)) +
  geom_point() +
  # theme_base +
  labs(x = "Spearman correlation",
       y = "lm beta")

temp_object@variable_info



write.xlsx(temp_object@variable_info,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_L_PFHxS_meta.xlsx")

# load("result/20241222_PFAS_META_sperm/data_preparation/cor_meta_volume.xlsx")
###volcano plot  火山图
variable_info <-
  extract_variable_info(temp_object)

top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(10) 
# dplyr::pull(Compound.name)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) 
# dplyr::pull(Compound.name)

volcano_plot <-
  variable_info %>%
  mutate(marker = case_when(
    cor_p_adjust < 0.2 & spearman_cor > 0 ~ "Up",
    cor_p_adjust < 0.2 &
      spearman_cor < 0 ~ "Down",
    TRUE ~ "No"
  )) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)") 
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_up_marker_name,
#                  Compound.name, NA)
# ), size = 3) +
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_down_marker_name,
#                  Compound.name, NA)
# ), size = 3)

volcano_plot

sum(temp_object@variable_info$cor_p < 0.05)
sum(temp_object@variable_info$cor_p_adjust < 0.05)
sum(temp_object@variable_info$permutated_p_value < 0.05)
sum(temp_object@variable_info$permutated_p_adjust < 0.05)
sum(temp_object@variable_info$lm_p < 0.05)
sum(temp_object@variable_info$lm_p_adjust < 0.05)

####过滤出p<0.05或校正后的p<0.2，过滤出与精液质量相关的代谢物质
variable_info

PFAS_markers <-
  variable_info %>%
  dplyr::filter(cor_p < 0.05)

temp_object@variable_info$variable_id
###进行主成分分析

pca_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% PFAS_markers$variable_id) %>%
  massstat::run_pca()

temp_object@sample_info

###画出主成分的散点图
plot <-
  massstat::pca_score_plot(
    object = temp_object,
    pca_object = pca_object,
    color_by = "L_PFHxS",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
# scale_x_continuous(limits = c(-0.2, 0.2)) +
# scale_y_continuous(limits = c(-0.2, 0.2))

plot

subject_info <-
  sample_info %>%
  dplyr::select(sample_id, L_PFHxS) %>%
  dplyr::filter(!is.na(L_PFHxS)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)



###提取出所有的主成分
lapply(1:88, function(i){
  data.frame(
    subject_id = rownames(pca_object$x),
    PC = apply(pca_object$x[,i, drop = FALSE],1, mean),
    class = i
  ) %>%
    as_tibble()
}) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



###提取出第一主成分
data.frame(
  subject_id = rownames(pca_object$x),
  PC = apply(pca_object$x[, 1, drop = FALSE], 1, mean),
  class = "L_PFHxS"
) %>%
  as_tibble() -> df_pca


df_pca$subject_id <- as.numeric(df_pca$subject_id)

df_pca %>% 
  dplyr::left_join(subject_info,by = c("subject_id"="sample_id"))  %>% 
  dplyr::filter(!is.na(L_PFHxS)) -> temp_data 

write.xlsx(temp_data,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_L_PFHxS.xlsx")

###第一组组成分与精液质量的回归系数及直线图

names(temp_object@sample_info)

omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "L_PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "PF3ONS_9CL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6]
  )


temp_data %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(L_PFHxS, PC, method = "spearman"))



cor.test(temp_data$L_PFHxS,temp_data$PC,   method = "spearman")



####绝对值回归系数与p值
library(plyr)
temp_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$L_PFHxS, x$PC,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()



temp_data %>%
  ggplot(aes(L_PFHxS,PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  # theme_base +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  # facet_wrap(facets = vars(class),
  #            scales = "free",
  #            nrow = 2) +
  scale_color_manual(values = omics_color) +
  labs(x = "L_PFHxS", y = "PC1")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        legend.position = c("none"),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 1, colour = "black")  # 设置边框线的粗细和颜色
  )+
  annotate("text", 
           x = Inf, y = Inf, label = "Abs(Spearman cor) = 0.17; P Value < 0.001", 
           hjust =1.8, vjust = 4, size = 5, color = "black", fontface = "bold") -> p7

p7



cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/pfas_meta_L_PFHxS_merger_line.pdf",width = 12, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p7   # 使用大写罗马表示
par(opar)
dev.off()








####>>>>>>>>>>>>>>4.对PF3ONS_9CL进行spearman回归分析<<<<<<<<<<<<<<<<<<<<<<<
temp_object@sample_info

names(temp_object@sample_info)

cor_data <-
  lapply(seq_along(rownames(temp_object)), function(i){
    value <- as.numeric(expression_data[i,]) 
    cor_result <- cor.test(value,temp_object@sample_info[,c(16)], method = "spearman")
    
    data.frame(variable_id = rownames(temp_object)[i],
               cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



names(temp_object@sample_info)



cor_data %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/result_PF3ONS_9CL_meta_cor.xlsx")

###>>>>>>>>>>>>>>>>>>>回归系数及校正p值<<<<<<<<<<<<<<<<<<<<<<===
##### lg_volume  lg_density lg_sperm_count lg_total_motility lg_progressive_motility
##### 
# ###permutation to get the p value   根据组合生成p值
# temp_object@sample_info[,c(67:71)
# dir.create("permutaton_cor_data")
for (idx in 1:271) {
  cat(idx, " ")
  permutation_cor_data <-
    seq_len(nrow(temp_object)) %>%
    purrr::map(function(i) {
      if((i %% 271) == 0)
        cat(i, " ")
      value <-
        as.numeric(unlist(temp_object[i, , drop = TRUE]))
      cor_result <-
        cor.test(value, sample(temp_object@sample_info$PF3ONS_9CL,
                               replace = FALSE), method = "spearman")
      data.frame(
        variable_id = rownames(temp_object)[i],
        cor_p = cor_result$p.value,
        spearman_cor = cor_result$estimate
      )
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  save(permutation_cor_data,
       file = file.path("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",paste0("permutation_cor_data_", idx)), compress = "xz")
  
}

load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/permutation_cor_data_1")

permutaton_cor_data_all <-
  permutation_cor_data %>%
  dplyr::select(-cor_p)

for (idx in 2:271) {
  cat(idx, " ")
  load(file.path(
    "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/",
    paste0("permutation_cor_data_", idx)
  ))
  permutaton_cor_data_all <-
    cbind(permutaton_cor_data_all, permutation_cor_data[, 3, drop = FALSE])
}

rownames(permutaton_cor_data_all) <- NULL

permutaton_cor_data_all <-
  permutaton_cor_data_all %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "variable_id")



###提取出精液体系的相关系数


permutated_p_value <-
  1:nrow(cor_data) %>%
  purrr::map(function(idx) {
    cat(idx, " ")
    original_cor <-
      cor_data$spearman_cor[idx]
    permutation_cor <-
      sample(as.numeric(permutaton_cor_data_all[idx,]), 10000, replace = TRUE)
    if (original_cor > 0) {
      sum(permutation_cor > original_cor) / 10000
    } else{
      sum(permutation_cor < original_cor) / 10000
    }
    
  }) %>%
  unlist()

cor_data$permutated_p_value <-
  permutated_p_value

cor_data

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p = cor_data$cor_p,
    permutated_p_value = cor_data$permutated_p_value,
    spearman_cor = cor_data$spearman_cor
  )

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p_adjust = p.adjust(cor_p, method = "fdr"),
    permutated_p_adjust = p.adjust(permutated_p_value, method = "fdr")
  )

temp_object@variable_info

plot(cor_data$cor_p,
     cor_data$permutated_p_value)


temp_object@sample_info



######linear mixed model

# 
linear_model_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    lm_result <-
      glm(formula = PF3ONS_9CL ~ value, data = temp_data)
    
    lm_result <-
      lm_result %>%
      broom::tidy()
    data.frame(lm_p = lm_result$p.value[2],
               coefficient = lm_result$estimate[2])
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()


save(linear_model_data, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/linear_model_data")



temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p = linear_model_data$lm_p,
                coefficient = linear_model_data$coefficient)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p_adjust = p.adjust(lm_p, method = "fdr"))

temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  ggplot(aes(x = spearman_cor, coefficient)) +
  geom_point() +
  # theme_base +
  labs(x = "Spearman correlation",
       y = "lm beta")

temp_object@variable_info



write.xlsx(temp_object@variable_info,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/cor_PF3ONS_9CL_meta.xlsx")

# load("result/20241222_PFAS_META_sperm/data_preparation/cor_meta_volume.xlsx")
###volcano plot  火山图
variable_info <-
  extract_variable_info(temp_object)

top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(10) 
# dplyr::pull(Compound.name)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) 
# dplyr::pull(Compound.name)

volcano_plot <-
  variable_info %>%
  mutate(marker = case_when(
    cor_p_adjust < 0.2 & spearman_cor > 0 ~ "Up",
    cor_p_adjust < 0.2 &
      spearman_cor < 0 ~ "Down",
    TRUE ~ "No"
  )) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)") 
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_up_marker_name,
#                  Compound.name, NA)
# ), size = 3) +
# ggrepel::geom_text_repel(aes(
#   label = ifelse(Compound.name %in% top10_down_marker_name,
#                  Compound.name, NA)
# ), size = 3)

volcano_plot

sum(temp_object@variable_info$cor_p < 0.05)
sum(temp_object@variable_info$cor_p_adjust < 0.05)
sum(temp_object@variable_info$permutated_p_value < 0.05)
sum(temp_object@variable_info$permutated_p_adjust < 0.05)
sum(temp_object@variable_info$lm_p < 0.05)
sum(temp_object@variable_info$lm_p_adjust < 0.05)

####过滤出p<0.05或校正后的p<0.2，过滤出与精液质量相关的代谢物质
variable_info

PFAS_markers <-
  variable_info %>%
  dplyr::filter(cor_p < 0.05)

temp_object@variable_info$variable_id
###进行主成分分析

pca_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% PFAS_markers$variable_id) %>%
  massstat::run_pca()

temp_object@sample_info

###画出主成分的散点图
plot <-
  massstat::pca_score_plot(
    object = temp_object,
    pca_object = pca_object,
    color_by = "PF3ONS_9CL",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
# scale_x_continuous(limits = c(-0.2, 0.2)) +
# scale_y_continuous(limits = c(-0.2, 0.2))

plot

subject_info <-
  sample_info %>%
  dplyr::select(sample_id,PF3ONS_9CL) %>%
  dplyr::filter(!is.na(PF3ONS_9CL)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)



###提取出所有的主成分
lapply(1:163, function(i){
  data.frame(
    subject_id = rownames(pca_object$x),
    PC = apply(pca_object$x[,i, drop = FALSE],1, mean),
    class = i
  ) %>%
    as_tibble()
}) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



###提取出第一主成分
data.frame(
  subject_id = rownames(pca_object$x),
  PC = apply(pca_object$x[, 1, drop = FALSE], 1, mean),
  class = "PF3ONS_9CL"
) %>%
  as_tibble() -> df_pca


df_pca$subject_id <- as.numeric(df_pca$subject_id)

df_pca %>% 
  dplyr::left_join(subject_info,by = c("subject_id"="sample_id"))  %>% 
  dplyr::filter(!is.na(PF3ONS_9CL)) -> temp_data 

write.xlsx(temp_data,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/temp_data_PF3ONS_9CL.xlsx")

###第一组组成分与精液质量的回归系数及直线图

names(temp_object@sample_info)

omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "L_PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "PF3ONS_9CL" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6]
  )


temp_data %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(PF3ONS_9CL, PC, method = "spearman"))



cor.test(temp_data$PF3ONS_9CL,temp_data$PC,   method = "spearman")



####绝对值回归系数与p值
library(plyr)
temp_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$PF3ONS_9CL, x$PC,  method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()



temp_data %>%
  ggplot(aes(PF3ONS_9CL,PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  # theme_base +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = TRUE) +
  # facet_wrap(facets = vars(class),
  #            scales = "free",
  #            nrow = 2) +
  scale_color_manual(values = omics_color) +
  labs(x = "PF3ONS_9CL", y = "PC1")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        legend.position = c("none"),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 1, colour = "black")  # 设置边框线的粗细和颜色
  )+
  annotate("text", 
           x = Inf, y = Inf, label = "Abs(Spearman cor) = 0.25; P Value < 0.001", 
           hjust =1.8, vjust = 4, size = 5, color = "black", fontface = "bold") -> p8

p8



cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/pfas_meta_PF3ONS_9CL_merger_line.pdf",width = 12, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p8   # 使用大写罗马表示
par(opar)
dev.off()


