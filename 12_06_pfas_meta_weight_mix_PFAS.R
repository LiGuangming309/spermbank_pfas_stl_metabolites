rm(list = ls())

pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts)
# get_project_wd()
# masstools::setwd_project()
rm(list = ls())
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


write.xlsx(df_background_DNA_LT_tidy,"result/20240407_PFOA_sperm/df_background_DNA_LT_tidy.xlsx")

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











###>>>>>>>>>>>>>>>>>>>>>>>>>>>1.读入单个PFAS暴露系数<<<<<<<<<<<<<<<<<<<<<<<<####

df1 <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_only_p_value_mix_liner.xlsx") 

df1


###2.根据回归系数计算权重
# # df_weighted <- df1 %>%
#   # 计算每个代谢物质中PFAS的绝对回归系数
#   group_by(metabolite) %>%
#   mutate(abs_estimate = abs(estimate.x)) %>%
#   # 计算每个代谢物质中的总权重
#   mutate(total_weight = sum(abs_estimate)) %>%
#   # 计算每个PFAS的相对权重
#   mutate(weight = abs_estimate / total_weight) %>%
#   select(metabolite, var.names.x, estimate.x, weight)



###将回归系数直接作为权重
df_weighted <- df1 %>%
  # # 计算每个代谢物质中PFAS的绝对回归系数
  # group_by(metabolite) %>%
  # mutate(abs_estimate = abs(estimate.x)) %>%
  # # 计算每个代谢物质中的总权重
  # mutate(total_weight = sum(abs_estimate)) %>%
  # # 计算每个PFAS的相对权重
  # mutate(weight = abs_estimate / total_weight) %>%
  select(metabolite, var.names.x, estimate.x)


df_weighted

library(stringr)

# 使用 str_replace() 替换掉 ".beta" 后缀
df_weighted$var.names.x <- str_replace(df_weighted$var.names.x, ".beta$", "")

df_weighted


###根据每个代谢物对应的PFAS暴露权重，计算每个代谢物对应的PFAS混合暴露量


###3.读入人群和代谢数据
df2 <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_meta_mix_rawdata.xlsx")  ###此处的pfas和端粒长度均取对数了，并未进行标化或归一化



names(df2)
str(df2)

df2$ID <- paste0("id", seq(1, 849))

df2$ID

####4.根据代谢物将PFAS变量转化为长数据
df2 %>%
  pivot_longer(
    cols = M642:M625,  # 选择所有化合物列
    names_to = "compound",    # 新列用于存储化合物名称
    values_to = "concentration"  # 新列用于存储浓度
  ) %>%
  pivot_longer(
    cols =  PFOS:PF3ONS_9CL,           # 选择 M1, M2, M3 列
    names_to = "variable",   # 新列用于存储变量（M1, M2, M3）
    values_to = "value"      # 新列用于存储对应的值
  )-> df_long





library(dplyr)

##### 5.将化合物和系数权重进行合并，根据 compound 和 variable 进行合并
df_merged <- df_long %>%
  left_join(df_weighted, by = c("compound" = "metabolite", "variable" = "var.names.x"))

# 查看合并后的结果
head(df_merged)

df_merged$estimate.x


str(df_long)

names(df_merged)

####6.生成每个化合物对应的PFAS的混合暴露浓度  

df_merged$value

df_merged %>% 
  mutate(PFAS_MIX=value*estimate.x) -> df_merged_final

names(df_merged_final)


df_merged_final$PFAS_MIX

type(df_merged_final$PFAS_MIX)

df_merged_final$compound

df_merged_final$ID
####7.对每个代谢物对应的化合物的PFAS进行求和，算出混合暴露的估计的PFAS暴露剂量

df_merged_final %>%
  dplyr::group_by(compound,ID) %>%       # 根据 sampleID 和 variable 分组
  summarise(total_PFAS_MIX = sum(PFAS_MIX, na.rm = TRUE)) %>%  # 对 PFAS_MIX 求和，na.rm = TRUE 忽略 NA 值
  ungroup()  ->  df_merged_final_pfas      


df_merged_final_pfas

# 修改 compound 列并添加 .PFAS 后缀
df_merged_final_pfas_2 <- df_merged_final_pfas %>%
  mutate(compound = paste0(compound, ".PFAS"))






###8。将长数据转为宽数据，将计算出代谢物对应的混合暴露的PFAS剂量 

df_merged_final_pfas_2 %>%
  pivot_wider(names_from = compound, values_from = total_PFAS_MIX) ->df_merged_final_pfas_mix 


names(df_merged_final_pfas_mix)


###9.将人群和代谢数据与混合暴露的PFAS进行合并

df2 %>% 
  left_join(df_merged_final_pfas_mix,by="ID") ->  df_merged_final_pfas_mix_meta


str(df_merged_final_pfas_mix_meta)
names(df_merged_final_pfas_mix_meta)



#########9.这里只需要对代谢物和端粒、线粒体拷贝数进行取对数和标准化
#### 1.先取对数
df_merged_final_pfas_mix_meta %>%
  select(M642:M625) %>%
  `+`(1) %>%
  log2() %>%
  apply(2, function(x) {     ###对每个代谢物进行z-score标化
    (x - mean(x)) / sd(x)
  }) -> df3

# log(62.8,2)


#### 标化后的数据
df_merged_final_pfas_mix_meta  %>%
  select(-c(M642:M625))%>%
  cbind(df3) -> df3_pfas_meta_stl     ###此处的代谢物是已经进行log2对数转化和标化后的物质，可以用作混合线性模型分析

str(df3_pfas_meta_stl)

df3_pfas_meta_stl   ###
names(df3_pfas_meta_stl)


###1.1对端粒进行归一化---在原始的变量上求均值 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  # select(mtDNAcn,STL) %>%
  select(mtDNAcn,STL) %>%
  apply(2, function(x) {     ###对端粒和线粒体拷贝数进行z-score标化
    (x - mean(x)) / sd(x)
  }) %>%
  as_tibble() -> df2_meta_stl_zscore

str(df2_meta_stl_zscore)



df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number[duplicated(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$number)]# 返回所有重复的行
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

df_background_pfas_stl_zscore %>% 
  left_join(df3_pfas_meta_stl[,-c(2:25)],by="number") -> df_background_pfas_stl_zscore_final

str(df_background_pfas_stl_zscore_final)

colSums(is.na(df_background_pfas_stl_zscore_final))

df_background_pfas_stl_zscore_final %>% 
  dplyr::filter(!M33=="NA")  %>% 
  as_tibble()->df_background_pfas_stl_zscore_final_tidy  ###得到的精浆和人群基线数据一样的数据





#######2.再求均值

str(df3_pfas_meta_stl)
names(df3_pfas_meta_stl)


df3_pfas_meta_stl$STL

library(dplyr)
library(rlang)


# 假设 df_background_DNA_LT 是你的数据框
res2 <- lapply(names(df3_pfas_meta_stl[, c(26:739)]), function(x) {
  # 获取列名
  result <- df3_pfas_meta_stl%>%
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

write.xlsx(df_background_DNA_LT_meta_mean_final,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/Mix_PFAS_RAW_DATA.xlsx")








