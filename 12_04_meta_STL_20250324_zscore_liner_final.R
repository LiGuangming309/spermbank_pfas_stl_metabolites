
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

write.xlsx(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA,"result/20240407_PFOA_sperm/df_background_DNA_LT_tidy_sperm_metabolites.xlsx")


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

write.xlsx(df_background_pfas_stl_zscore_final_tidy,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy.xlsx")



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

write.xlsx(df_background_DNA_LT_meta_mean_final,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy_623.xlsx")


# df_background_pfas_stl_zscore_final_tidy

# df2_meta_stl_tidy

########>>>>>>>>>>>>一、混合线性回归分析中<<<<<<<<<<<<<<<<<<<<<<<##############

# ###利用取对数和归一化之后的端粒和代谢物，进行混合线性回归分析
# df_background_pfas_stl_zscore_final_tidy    ###混合线性回归分析
# 
# df_background_pfas_stl_zscore_final_tidy$mtDNAcn
# 
# df_background_pfas_stl_zscore_final_tidy$STL
# 
# names(df_background_pfas_stl_zscore_final_tidy)
# 
# 
# 
# 
# ###3.进行贝叶斯和蒙特卡洛回归分析
# library(brms)
# ?brm()
# 
# # 定义模型公式
# formula <- bf(STL~ M33+age_class_2 + BMI_class_2 + education + 
#                 marriage_class_2 + income + smk_class_2 + drk_class_2 +
#                 season + abstinence_class_2+(1|number))
# 
# 
# # install.packages("brms")
# 
# 
# 
# 
# # 运行贝叶斯混合效应模型
# fit01 <- brms::brm(
#   formula = formula,
#   data = df_background_pfas_stl_zscore_final_tidy,
#   family = gaussian,
#   chains = 4,
#   iter = 500,
#   warmup = 100,
#   control = list(adapt_delta = 0.95),
#   backend = "rstan"
# )
# 
# summary(fit01) -> fit1
# 
# fit1$fixed -> result
# 
# result
# 
# names()
# 
# # 使用 lapply 对第 25 到 26 列进行循环
# 
# lapply(names(df2_meta_stl_tidy[,c(25:402)]), function(col_name) {
#   
#   # 定义动态模型公式，插入当前列名
#   formula_str <- paste0("STL ~ ", col_name, " + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + drk_class_2 + season + abstinence_class_2+(1|number)")
#   
#   # 加载 brms 包
#   library(brms)
#   
#   # 运行贝叶斯混合效应模型
#   fit01 <- brms::brm(
#     formula = as.formula(formula_str),
#     data = df2_meta_stl_tidy,
#     family = gaussian,
#     chains = 4,
#     iter = 500,
#     warmup = 100,
#     control = list(adapt_delta = 0.95),
#     backend = "rstan"
#   )
#   
#   # 获取模型摘要
#   summary(fit01) -> fit1
#   fit1$fixed -> result
#   
#   result %>% 
#     mutate(p_value=case_when(
#       `l-95% CI` < 0 & `u-95% CI` > 0 ~ ">0.05",  # 置信区间跨越零，p值大于0.05
#       TRUE ~ "<0.05"  # 否则 p值小于0.05
#     )) -> result
#   # 返回模型结果
#   return(result)
#   
# }) -> res1  ###将循环中每个y生成的list结果保存为res1
# 
# 
# 
# res1 -> res5
# 
# res5[[1]]
# rownames(res5[[1]])
# 
# seq_along(res5)
# 
# # lapply(seq_along(res5), function)
# 
# res5[[1]] %>% cbind(rownames(res5[[1]]))
# 
# ###然后再将生成的list重新保存为独立的excel文件
# ###给每个list添加行变量名称
# lapply(seq_along(res5), function(i) {
#   res5[[i]] %>% cbind(rownames(res5[[i]]))
# }) -> res6
# 
# 
# # 
# # res6 %>% 
# #   bind_rows() %>% 
# #   write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/result_meta_STL_log_mean_all_meta_final_20241201_bayes.xlsx")
# # 
# # 
# # 
# 
# ###4.进行混合线性回归分析
# 
# df_background_pfas_stl_zscore_final_tidy   ####此处的数据集，没有对端粒和线粒体拷贝数进行z-score变换
# names(df_background_pfas_stl_zscore_final_tidy)
# 
# library(lme4)
# 
# str(df_background_pfas_stl_zscore_final_tidy)
# 
# # 检查共线性
# library(car)
# vif(lm(STL ~ age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + drk_class_2 + season + abstinence_class_2, data = df_background_pfas_stl_zscore_final_tidy))
# 
# ####教育与其他协变量可能存在多重共线性
# 
# 
# 
# 
# lapply(df_background_pfas_stl_zscore_final_tidy[,c(25:402)], function(x){
#     lmer(STL~x+age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + drk_class_2 + season + abstinence_class_2+(1|number),data=df_background_pfas_stl_zscore_final_tidy) -> result
#     summary(result) -> res1
#     res1$coefficients -> res2
#     β <- res2[,1]
#     SE <- res2[,2]
#     CI2.5 <- β-1.96*SE
#     CI97.5  <- β+1.96*SE
#     CI95<-paste0(round(β,2)," (",round(CI2.5,2),', ',round(CI97.5,2),")")
#     percent_β <- (exp(β)-1)*100
#     percent_CI2.5 <- (exp(β-1.96*SE)-1)*100
#     percent_CI97.5 <- (exp(β+1.96*SE)-1)*100
#     percent_CI95 <- paste0(round(percent_β,2)," (",round(percent_CI2.5,2),', ',round(percent_CI97.5,2),")")
#     res2 %>% 
#       rownames() %>% 
#       as_tibble() -> res ####自己添加的，将其他输出的变量也添加到包里
#     Uni_glm_model <- data.frame('Estimate'=res2[,1],
#                                 "Std. Error"=res2[,2],
#                                 'df' = res2[,3],
#                                 't value' = res2[,4],
#                                 'p value' = res2[,5],
#                                 'CI2.5' = CI2.5,
#                                 'CI97.5' = CI97.5,
#                                 'CI95' = CI95,
#                                 "percent_β"=percent_β,
#                                 "percent_CI2.5"=percent_CI2.5,
#                                 "percent_CI97.5"=percent_CI97.5,
#                                 "percent_CI95"=percent_CI95,
#                                 res)
#     return(Uni_glm_model)}) -> res7  ###将跑的循环中y生成的list结果保存为res1
# 
# 
#   
# 
# # res7 %>% 
# #   bind_rows()%>%
# #   write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_meta_STL_log_mean_all_meta_final_20241201_zscore.xlsx")
# # # 
#   




########>>>>>>>>>>>>二、线性回归分析中<<<<<<<<<<<<<<<<<<<<<<<##############

df_background_DNA_LT_meta_mean_final 
names(df_background_DNA_LT_meta_mean_final)



lapply(df_background_DNA_LT_meta_mean_final[,c(25:295)], function(x){
  lm(STL~x+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,data=df_background_DNA_LT_meta_mean_final) -> result
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
  return(Uni_glm_model)}) -> res8  ###将跑的循环中y生成的list结果保存为res1




res8 %>% 
  bind_rows()%>%
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner.xlsx")
# 
# colnames(df_background_DNA_LT_q[,c(47:297)]) %>%
#   as_tibble()  %>% 
#   dplyr::rename("meta"="value") -> name_meta





colnames(df_background_DNA_LT_meta_mean_final[,c(25:295)]) %>%
  as_tibble()  %>% 
  dplyr::rename("meta"="value") -> name_meta


###导出所有的代谢物
# res1 %>% bind_rows()%>%as_tibble()%>%
#   write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/result_PFAS.xlsx")



###代谢物与端粒之间的混合线性模型结果
res8 %>% 
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
  cbind(q_value)%>%
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q.xlsx")

# filter(q_value<0.05) -> df_meta_stl_q




##

#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
####以上分析均是针对代谢物质和端粒之间的分析



###>>>>>>>>>>>>>>>>>>>>>>三、读入经过BHRMA筛选后的pfas混合物---混合线性模型的<<<<<<<<<<<<<<<<<<

df_pfas_mix_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value_mix_liner.xlsx")


###1.提取经过BH校正之后q_value<0.2的代谢物 
###
names(df_pfas_mix_meta)

df_pfas_mix_meta %>% 
  filter(q_value<0.2) -> df_pfas_mix_meta_q   ###总共获得79个代谢物质   


df_pfas_mix_meta_q 

str(df_pfas_mix_meta_q)
names(df_pfas_mix_meta_q)

# 
# ####>>>>>>>>读入代谢物匹配的代谢通路变量名称#############################
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
# 
# df_pfas_mix_meta_q %>% 
#   left_join(df_kegg_1,by="metabolite") -> df_pfas_mix_meta_q_kegg
# 
# 
# names(df_pfas_mix_meta_q_kegg)
# 
# 
# colSums(is.na(df_pfas_mix_meta_q_kegg))
# 
# str(df_pfas_mix_meta_q_kegg)   ####PFAS相关的只有251个代谢物   
# 





###>>>>>>>>>>>>>>>>四、过滤出与PFAS有相关性的代谢物质，校正后q<0.20的物质，再进行回归分析<<<<<<<<<<<<<<<<<#########

names(df_background_DNA_LT_meta_mean_final)


df_background_DNA_LT_meta_mean_final  %>% 
  select(number:STL,df_pfas_mix_meta_q$metabolite) -> df_background_DNA_LT_q   ###过滤出来只与PFAS有关连的代谢物


str(df_background_DNA_LT_q)

names(df_background_DNA_LT_q)



lapply(df_background_DNA_LT_q[,c(25:103)], function(x){
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
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_PFAS.xlsx")
# 
# colnames(df_background_DNA_LT_q[,c(47:297)]) %>%
#   as_tibble()  %>% 
#   dplyr::rename("meta"="value") -> name_meta





colnames(df_background_DNA_LT_q[,c(25:103)]) %>%
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
     write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q_PFAS.xlsx")
















