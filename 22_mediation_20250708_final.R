
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


str(df_background_DNA_LT_sperm_PFOA_tidy)
df_background_DNA_LT_sperm_PFOA_tidy



df_background_DNA_LT_sperm_PFOA_tidy %>% 
  mutate(PFAS_mix=PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL) -> df_background_DNA_LT_sperm_PFOA_tidy_1


str(df_background_DNA_LT_sperm_PFOA_tidy_1)


df_background_DNA_LT_sperm_PFOA_tidy_1  %>% 
  left_join(df_meta_2,by="number") -> df_background_DNA_LT_sperm_PFOA_meta_tidy

colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy))




####删除缺失值

df_background_DNA_LT_sperm_PFOA_meta_tidy  %>% 
  dplyr::filter(!M33=="NA")  %>% 
  as_tibble()->df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA  ###得到的精浆和人群基线数据一样的数据



colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA))

df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA

names(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA) ###得到849个重复测量的样本

str(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA)

# length(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$number)
# df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA[duplicated(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$number), ]  ###重复的行数


#########代谢物和端粒、线粒体拷贝数的整理
#### 1.先取对数
df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA %>%
  select(PFAS_mix,M661:M625) %>%
  `+`(1) %>%
  log2() %>%
  apply(2, function(x) {     ###对每个代谢物进行z-score标化
    (x - mean(x)) / sd(x)
  }) -> df1

# log(62.8,2)


#### 标化后的数据
df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA  %>%
  select(-c(PFAS_mix,M661:M625))%>%
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

df_background_pfas_stl_zscore[,-c(2:22)] %>% 
  left_join(df2_meta_stl,by="number") -> df_background_pfas_stl_zscore_final

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
res2 <- lapply(names(df2_meta_stl[, c(27:298)]), function(x) {
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

names(df_background_pfas_stl_zscore)

df_background_pfas_stl_zscore[,-c(11:22)] %>% 
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



####----------取出pfas与STL

###取出精浆代谢与端粒长度相关的代谢物
meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q.xlsx")

meta_STL %>%
  filter(q_value<0.2) -> meta_STL_1

meta_STL_1

###取出PFAS混合暴露与精液代谢相关的代谢物
pfas_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value_mix_liner.xlsx")

pfas_meta %>%
  filter(q_value<0.2) ->pfas_meta_1

names(meta_STL)
names(pfas_meta)

pfas_meta_1 %>%
  left_join(meta_STL_1,by=c("metabolite"="meta")) -> pfas_meta_stl

colSums(is.na(pfas_meta_stl))

pfas_meta_stl %>%
  filter(!Estimate=="NA") -> pfas_meta_stl_all



pfas_meta_stl_all$metabolite







###匹配出暴露和代谢物质中包含df6中的代谢物

str(df_background_DNA_LT_meta_mean_final)

names(df_background_DNA_LT_meta_mean_final)

df_background_DNA_LT_meta_mean_final %>%
  select(number:PFAS_mix,matches(paste(pfas_meta_stl_all$metabolite, collapse = "|"))) -> df_background_DNA_LT_meta_mean_final_pfas



# 
names(df_background_DNA_LT_meta_mean_final_pfas)
# 
str(df_background_DNA_LT_meta_mean_final_pfas)
# 
names(df_background_DNA_LT_meta_mean_final_pfas)
# 
# 
















####--------------------第一种 中介分析------------------------------------------

df_background_DNA_LT_meta_mean_final_pfas -> df_clean

names(df_clean)

str(df_clean)



fita <- lm(M558 ~ PFAS_mix + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
             drk_class_2 + season + abstinence_class_2, data=df_clean)

summary(fita)
# 第二阶段回归模型
fitb <- lm(STL ~M558 + PFAS_mix+ age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
             drk_class_2 + season + abstinence_class_2, data=df_clean)

summary(fitb)



fita <- lm(M378 ~ PFAS_mix + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
             drk_class_2 + season + abstinence_class_2, data=df_clean)

summary(fita)
# 第二阶段回归模型
fitb <- lm(STL ~M378 + PFAS_mix+ age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
             drk_class_2 + season + abstinence_class_2, data=df_clean)

summary(fitb)


names(df_clean)

library(mediation)
####*****进行中介分析，将x和m对应的变量，同时纳入模型进行中介分析<<<<<<<<######

lapply(seq(14,31), function(b) {
     x <- b
  # 
  # # 确保自变量和中介变量列没有 NA 值
  # df_clean <- na.omit(df8)
  
  lapply(df_clean[[names(df_clean)[x]]] %>% as_tibble, function(m1) {   ####中介变量
    
    # lapply(df_clean[[names(df_clean)[m]]] %>% as_tibble, function(m1) {    ###因变量
      
      # 第一阶段回归模型
      fita <- lm(m1 ~ PFAS_mix + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
                   drk_class_2 + season + abstinence_class_2, data=df_clean)
      
      # 第二阶段回归模型
      fitb <- lm(STL ~ m1 + PFAS_mix+ age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
                   drk_class_2 + season + abstinence_class_2, data=df_clean)

      
      library(mediation)
      set.seed(2022)
      
      result <- mediate(fita, fitb, treat='PFAS_mix', mediator="m1", boot=TRUE, sims=1000)
      
      
      return_result <- data.frame(
        term=names(df_clean)[x],
        acme = result$d1,
        acme_ci_lower = result$d1.ci[1],
        acme_ci_upper = result$d1.ci[2],
        acme_p =formatC(result$d1.p, format = "f", digits = 10),
        ade = result$z1,
        ade_ci_lower = result$z1.ci[1],
        ade_ci_upper = result$z1.ci[2],
        ade_p = formatC(result$z1.p, format = "f", digits = 10),
        total_effect = result$tau.coef,
        total_effect_ci_lower = result$tau.ci[1],
        total_effect_ci_upper = result$tau.ci[2],
        total_effect_p =formatC(result$tau.p, format = "f", digits = 10), 
        prop_mediate = result$n1,
        prop_mediate_ci_lower = result$n1.ci[1],
        prop_mediate_ci_upper = result$n1.ci[2],
        prop_mediate_p = formatC(result$n1.p, format = "f", digits = 10)
      ) 
      
      return(return_result)
  #   })
  })
}) -> res4


res4 %>% 
  bind_rows()-> res5

res5$value

# write.csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/result_median_20250711.csv")
# 
# # 假设你有多个类似的嵌套结构，使用 lapply 逐个提取
# df_result <- lapply(res4, function(x) as.data.frame(x$value$value))
# df_result <- bind_rows(df_result)  # 合并成一个数据框
# 
# 
# names(df_clean)[seq(43, 60)] %>% as_tibble() -> names_m
# 
# 
# df_result %>% 
#   cbind(names_m) -> result_median
# 
# result_median


# 
# acme = result$d1,
# acme_ci_lower = result$d1.ci[1],
# acme_ci_upper = result$d1.ci[2],
# acme_p =formatC(result$d1.p, format = "f", digits = 10),
# ade = result$z1,
# ade_ci_lower = result$z1.ci[1],
# ade_ci_upper = result$z1.ci[2],
# ade_p = formatC(result$z1.p, format = "f", digits = 10),
# total_effect = result$tau.coef,
# total_effect_ci_lower = result$tau.ci[1],
# total_effect_ci_upper = result$tau.ci[2],
# total_effect_p =formatC(result$tau.p, format = "f", digits = 10), 
# prop_mediate = result$n1,
# prop_mediate_ci_lower = result$n1.ci[1],
# prop_mediate_ci_upper = result$n1.ci[2],
# prop_mediate_p = formatC(result$n1.p, format = "f", digits = 10)


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


df_unique_2$metabolite

# mwas_reduced$metabolite

###对中介变量中的每个数字附上代谢物质名称

res5$value %>%
  left_join(df_unique_3,by=c("term"="metabolite")) -> result_median_names


write.xlsx(result_median_names,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/result_median_20250711_final.xlsx")






###提取出中介中的a\b\c通路上点系数和p值


library(tibble)
library(dplyr)

# 提取并整理 fitM 模型结果的函数


# 

df_clean


extract_fitM <- function(PFAS_mix,m1, data) {
  # fitM <- lm(m1 ~ x2 + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
  #              drk_class_2 + season + abstinence_class_2, data = data)
  
  fitM <- lm(m1 ~ PFAS_mix + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
               drk_class_2 + season + abstinence_class_2, data=df_clean)
  
  
  summary_fitM <- summary(fitM)[[4]][2, ]
  
  data_M <- data.frame(
    Estimate = summary_fitM["Estimate"],
    Std.Error = summary_fitM["Std. Error"],
    t.value = summary_fitM["t value"],
    P.value = summary_fitM["Pr(>|t|)"]
  )
  
  a <- c("a") %>% as_tibble()
  
  # a %>%
  #   mutate(combined = paste0(value, "_", x_name)) -> a_1
  
  data_M %>% cbind(a) %>% as_tibble() -> medi_m
  
  return(medi_m)
}


# 提取并整理 fitY 模型结果的函数
extract_fitY <- function(PFAS_mix, m1, data) {
  # fitY <- lm(STL ~ m1 +x2+age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
  #              drk_class_2 + season + abstinence_class_2, data = data)
  
  # 第二阶段回归模型
  fitY <- lm(STL ~ m1 + PFAS_mix+ age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
               drk_class_2 + season + abstinence_class_2, data=df_clean)
  
  
  summary_fitY <- summary(fitY)[[4]][c(2, 3), ] %>% as_tibble()
  
  data_Y <- data.frame(
    Estimate = summary_fitY$Estimate,
    Std.Error = summary_fitY$`Std. Error`,
    t.value = summary_fitY$`t value`,
    P.value = summary_fitY$`Pr(>|t|)`
  ) %>% as_tibble()
  
  bc <- c("b", "c") %>% as_tibble()
  
  # bc %>%
  #   mutate(combined = paste0(value, "_", x_name)) -> bc_2
  
  data_Y %>%
    cbind(bc) %>% as_tibble() -> medi_Y
  
  medi_Y %>%
    dplyr::rename(Std.Error = `Std.Error`, P.value = `P.value`) -> medi_Y_2
  
  return(medi_Y_2)
}




# 主要 lapply 函数

lapply(seq(14, 31), function(b) {  
  # x <- b
  m <- b
  # 
  # # 确保自变量和中介变量列没有 NA 值
  # df_clean <- na.omit(df8)
  
  # lapply(df_clean[[names(df_clean)[x]]] %>% as_tibble, function(x_name) {   ####中介变量
  
  # lapply(df_clean[[names(df_clean)[27]]], function(m1) {    ###因变量
  
  # x2 <- df_clean[[x]]
  m1 <- df_clean[[m]]
  medi_m_3 <- extract_fitM(PFAS_mix,m1,df_clean)
  medi_Y_3 <- extract_fitY(PFAS_mix,m1,df_clean)
  
  names <- rep(names(df_clean)[m],times=3) %>% as_tibble()
  combined_result <- rbind(medi_m_3,medi_Y_3) %>% cbind(names)
  
  return(combined_result)
  
  
  # combined_result <- rbind(medi_m_3, medi_Y_3)
  # return(combined_result)
  # })
  # })
  
}) -> res6

# names(df_clean)[27]



res6 %>% 
  bind_rows() %>% 
  as_tibble()-> res7


res7 %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/median_estimate_20250711.xlsx")






####----------------------------第二种混合中介分析----------------------------------


df_background_pfas_stl_zscore_final_tidy

names(df_background_pfas_stl_zscore_final_tidy)


df_background_pfas_stl_zscore[,-c(11:22)] %>% 
  left_join(df_background_pfas_stl_zscore_final_tidy[,-c(2:28)],by="number") ->df_background_DNA_LT_meta_mean

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



####----------取出pfas与STL

###取出精浆代谢与端粒长度相关的代谢物
meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q.xlsx")

meta_STL %>%
  filter(q_value<0.2) -> meta_STL_1

meta_STL_1

###取出PFAS混合暴露与精液代谢相关的代谢物
pfas_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value_mix_liner.xlsx")

pfas_meta %>%
  filter(q_value<0.2) ->pfas_meta_1

names(meta_STL)
names(pfas_meta)

pfas_meta_1 %>%
  left_join(meta_STL_1,by=c("metabolite"="meta")) -> pfas_meta_stl

colSums(is.na(pfas_meta_stl))

pfas_meta_stl %>%
  filter(!Estimate=="NA") -> pfas_meta_stl_all



pfas_meta_stl_all$metabolite







###匹配出暴露和代谢物质中包含df6中的代谢物

str(df_background_DNA_LT_meta_mean_final)

names(df_background_DNA_LT_meta_mean_final)

df_background_DNA_LT_meta_mean_final %>%
  select(number:PFAS_mix,matches(paste(pfas_meta_stl_all$metabolite, collapse = "|"))) -> df_background_DNA_LT_meta_mean_final_pfas



# 
names(df_background_DNA_LT_meta_mean_final_pfas)
# 
str(df_background_DNA_LT_meta_mean_final_pfas)
# 
names(df_background_DNA_LT_meta_mean_final_pfas)
# 
# 


df_background_DNA_LT_meta_mean_final_pfas -> df_clean

names(df_clean)

# library(lmerTest)
unloadNamespace("lmerTest")  ###不加载该包
# detach("package:lmerTest", unload = TRUE)
# detach(lmerTest)
med.fit_1 <- lmer(M558 ~ PFAS_mix + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
                    drk_class_2 + season + abstinence_class_2+(1|number), data=df_clean)
out.fit_1 <- lmer(STL ~ M558 + PFAS_mix+(1|number), data=df_clean)


med.out_1 <- mediate(med.fit_1, out.fit_1, treat = "PFAS_mix", mediator = "M558", sims = 1000)

summary(med.out_1)    ###混合中介做出来的结果与之前的结果相反，不能使用混合中介进行计算









































