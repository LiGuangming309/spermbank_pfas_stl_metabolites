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
    cols = M661:M625,  # 选择所有化合物列
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
  select(M661:M625) %>%
  `+`(1) %>%
  log2() %>%
  apply(2, function(x) {     ###对每个代谢物进行z-score标化
    (x - mean(x)) / sd(x)
  }) -> df3

# log(62.8,2)


#### 标化后的数据
df_merged_final_pfas_mix_meta  %>%
  select(-c(M661:M625))%>%
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
res2 <- lapply(names(df3_pfas_meta_stl[, c(26:567)]), function(x) {
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

# write.xlsx(df_background_pfas_stl_zscore_final_tidy,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy_623.xlsx")

# lm(STL~M33+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,data=df_background_DNA_LT_meta_mean_final) -> result
# summary(result)



###10.根据计算的混合暴露浓度进行中介分析


# 
# ###读入PFAS与代谢有统计学意义的代谢物，这些代谢物又与端粒有统计学意义的代谢物
# 
# df5<- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q_PFAS.xlsx")
# 
# 
# ###筛选FDR校正之后小于0.2
# df5%>% 
#   filter(q_value<0.2) -> df6
# 
# names(df6)
# 
# df6$meta
# 
# names(df_background_DNA_LT_meta_mean_final)



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
df_background_DNA_LT_meta_mean_final %>% 
  select(number:STL,matches(paste(pfas_meta_stl_all$metabolite, collapse = "|"))) ->df_background_DNA_LT_meta_mean_final_pfas
  



names(df_background_DNA_LT_meta_mean_final_pfas)

str(df_background_DNA_LT_meta_mean_final_pfas)

names(df_background_DNA_LT_meta_mean_final_pfas)




###中介分析



# 需要的包
library(dplyr)
library(mediation)





df_background_DNA_LT_meta_mean_final_pfas -> df7

names(df7)
str(df7)

###X变量
# 选择 df7 中的列并获取列名
selected_columns <- df7[,c(25:42)]

# 获取列名并对其进行排序
sorted_column_names <- sort(names(selected_columns))

# 根据排序后的列名重新排序 df7 中的列
sorted_df <- selected_columns[, sorted_column_names]
sorted_df

str(sorted_df)
names(sorted_df)

###m变量
# 选择 df7 中的列并获取列名
selected_columns_2 <- df7[,c(43:60)]

# 获取列名并对其进行排序
sorted_column_names_2 <- sort(names(selected_columns_2))

# 根据排序后的列名重新排序 df7 中的列
sorted_df_2 <- selected_columns_2[, sorted_column_names_2]


names(sorted_df_2)



names(df7)

df7 %>% 
  select(number:STL)%>%
  cbind(sorted_df,sorted_df_2) -> df8

names(df8)
# 
# df8%>%
#   select(25+123)

names(df8) 
names(df8)[seq(25,60)]


###重新赋值给整理好的数据框
df8 -> df_clean

str(df_clean)
colSums(is.na(df_clean))

####对混合线性模型求均值

# 
# ###读入pfas、stl和代谢的原始物质，其中对端粒和代谢物质已经进行了标化
# df2 <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy_623.xlsx")
# 
# 
# 

df_clean


# 第一阶段回归模型
fita <- lm(M103 ~ M103.PFAS + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
             drk_class_2 + season + abstinence_class_2, data=df_clean)

# 第二阶段回归模型
fitb <- lm(STL ~ M103.PFAS + M103+ age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
             drk_class_2 + season + abstinence_class_2, data=df_clean)
# 定义新的print.summary.mediate函数


library(mediation)


set.seed(2022)

result <- mediate(fita, fitb, treat='M103.PFAS', mediator="M103", boot=TRUE, sims=100000)


return_result <- data.frame(
  acme = result$d1,
  acme_ci_lower = result$d1.ci[1],
  acme_ci_upper = result$d1.ci[2],
  acme_p =result$d1.p,
  ade = result$z1,
  ade_ci_lower = result$z1.ci[1],
  ade_ci_upper = result$z1.ci[2],
  ade_p = result$z1.p,
  total_effect = result$tau.coef,
  total_effect_ci_lower = result$tau.ci[1],
  total_effect_ci_upper = result$tau.ci[2],
  total_effect_p =result$tau.p, 
  prop_mediate = result$n1,
  prop_mediate_ci_lower = result$n1.ci[1],
  prop_mediate_ci_upper = result$n1.ci[2],
  prop_mediate_p = result$n1.p
) 

return_result <- data.frame(
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



return_result



names(df_clean) 

####*****进行中介分析，将x和m对应的变量，同时纳入模型进行中介分析<<<<<<<<######

lapply(seq(25,42), function(b) {  
  x <- b
  m <- b + 18
  # 
  # # 确保自变量和中介变量列没有 NA 值
  # df_clean <- na.omit(df8)
  
  lapply(df_clean[[names(df_clean)[x]]] %>% as_tibble, function(x1) {   ####中介变量
    
    lapply(df_clean[[names(df_clean)[m]]] %>% as_tibble, function(m1) {    ###因变量
      
      # 第一阶段回归模型
      fita <- lm(m1 ~ x1 + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
                   drk_class_2 + season + abstinence_class_2, data=df_clean)
      
      # 第二阶段回归模型
      fitb <- lm(STL ~ m1 + x1+ age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
                   drk_class_2 + season + abstinence_class_2, data=df_clean)
      
      trace(mediation:::print.summary.mediate, 
            at = 11,
            tracer = quote({
              printCoefmat <- function(x, digits) {
                p <- x[, 4] #p-values seem to be stored rounded
                x[, 1:3] <- sprintf("%.10f", x[, 1:3])
                x[, 4] <- sprintf("%.10f", p)
                print(x, quote = FALSE, right = TRUE)
              } 
            }),
            print = FALSE)
      
      library(mediation)
      set.seed(2022)
      
      result <- mediate(fita, fitb, treat='x1', mediator="m1", boot=TRUE, sims=1000)
      
      
      return_result <- data.frame(
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
    })
  })
}) -> res4


res4[[2]]



# 
# 
# lapply(seq(26, 53), function(b) {  
#   m <- b
#   x <- b + 28
#   # 
#   # # 确保自变量和中介变量列没有 NA 值
#   # df_clean <- na.omit(df8)
#   
#   lapply(df_clean[[names(df_clean)[x]]] %>% as_tibble, function(x1) {   ####中介变量
#     
#     lapply(df_clean[[names(df_clean)[m]]] %>% as_tibble, function(m1) {    ###因变量
#       
#       # 第一阶段回归模型
#       fita <- lm(m1 ~ x1 + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
#                    drk_class_2 + season + abstinence_class_2, data=df_clean)
#       
#       # 第二阶段回归模型
#       fitb <- lm(STL ~ m1 + x1 + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
#                    drk_class_2 + season + abstinence_class_2, data=df_clean)
#       
#       library(mediation)
#       set.seed(2022)
#       result <- mediate(fita, fitb, treat='x1', mediator="m1", boot=TRUE, sims=1000)
#       
#       return_result <- data.frame(
#         acme = result$d1,
#         acme_ci_lower = result$d1.ci[1],
#         acme_ci_upper = result$d1.ci[2],
#         acme_p = result$d1.p,
#         ade = result$z1,
#         ade_ci_lower = result$z1.ci[1],
#         ade_ci_upper = result$z1.ci[2],
#         ade_p = result$z1.p,
#         total_effect = result$tau.coef,
#         total_effect_ci_lower = result$tau.ci[1],
#         total_effect_ci_upper = result$tau.ci[2],
#         total_effect_p = result$tau.p,
#         prop_mediate = result$n1,
#         prop_mediate_ci_lower = result$n1.ci[1],
#         prop_mediate_ci_upper = result$n1.ci[2],
#         prop_mediate_p = result$n1.p
#       ) 
#       return(return_result)
#     })
#   })
# }) -> res4
# 



# 假设你有多个类似的嵌套结构，使用 lapply 逐个提取
df_result <- lapply(res4, function(x) as.data.frame(x$value$value))
df_result <- bind_rows(df_result)  # 合并成一个数据框


names(df_clean)[seq(43, 60)] %>% as_tibble() -> names_m


df_result %>% 
  cbind(names_m) -> result_median

result_median






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

result_median %>% 
  left_join(df_unique_3,by=c("value"="metabolite")) -> result_median_names


write.xlsx(result_median_names,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/result_median_20250326_final.xlsx")


####计算中介分析的系数   

# 
# lapply(seq(26, 53), function(b) {  
#   m <- b
#   x <- b + 28
#   # 
#   # # 确保自变量和中介变量列没有 NA 值
#   # df_clean <- na.omit(df8)
#   
#   lapply(df_clean[[names(df_clean)[x]]] %>% as_tibble, function(x1) {   ####中介变量
#     
#     lapply(df_clean[[names(df_clean)[m]]] %>% as_tibble, function(m1) {    ###因变量
#       
#       # 第一阶段回归模型
#       fita <- lm(m1 ~ x1 + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
#                    drk_class_2 + season + abstinence_class_2, data=df_clean)
#       
#       # 第二阶段回归模型
#       fitb <- lm(STL ~ m1 + x1 + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
#                    drk_class_2 + season + abstinence_class_2, data=df_clean)
#       
#       library(mediation)
#       set.seed(2022)
#       result <- mediate(fita, fitb, treat='x1', mediator="m1", boot=TRUE, sims=1000)
#       
#       return_result <- data.frame(
#         acme = result$d1,
#         acme_ci_lower = result$d1.ci[1],
#         acme_ci_upper = result$d1.ci[2],
#         acme_p = result$d1.p,
#         ade = result$z1,
#         ade_ci_lower = result$z1.ci[1],
#         ade_ci_upper = result$z1.ci[2],
#         ade_p = result$z1.p,
#         total_effect = result$tau.coef,
#         total_effect_ci_lower = result$tau.ci[1],
#         total_effect_ci_upper = result$tau.ci[2],
#         total_effect_p = result$tau.p,
#         prop_mediate = result$n1,
#         prop_mediate_ci_lower = result$n1.ci[1],
#         prop_mediate_ci_upper = result$n1.ci[2],
#         prop_mediate_p = result$n1.p
#       ) 
#       return(return_result)
#     })
#   })
# }) -> res4
# 




###提取出中介中的a\b\c通路上点系数和p值


library(tibble)
library(dplyr)

# 提取并整理 fitM 模型结果的函数


# 

df_clean


extract_fitM <- function(x2,m1, data) {
  fitM <- lm(m1 ~ x2 + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
               drk_class_2 + season + abstinence_class_2, data = data)
  
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
extract_fitY <- function(x2, m1, data) {
  fitY <- lm(STL ~ m1 +x2+age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + 
               drk_class_2 + season + abstinence_class_2, data = data)
  
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

lapply(seq(25, 42), function(b) {  
  x <- b
  m <- b + 18
  # 
  # # 确保自变量和中介变量列没有 NA 值
  # df_clean <- na.omit(df8)

  # lapply(df_clean[[names(df_clean)[x]]] %>% as_tibble, function(x_name) {   ####中介变量

    # lapply(df_clean[[names(df_clean)[27]]], function(m1) {    ###因变量

      x2 <- df_clean[[x]]
      m1 <- df_clean[[m]]
      medi_m_3 <- extract_fitM(x2,m1,df_clean)
      medi_Y_3 <- extract_fitY(x2,m1,df_clean)
      
      names <- rep(names(df_clean)[m],times=3) %>% as_tibble()
      combined_result <- rbind(medi_m_3,medi_Y_3) %>% cbind(names)
      
      return(combined_result)


    # combined_result <- rbind(medi_m_3, medi_Y_3)
    # return(combined_result)
  # })
# })
  
}) -> res3

# names(df_clean)[27]



res3 %>% 
  bind_rows() %>% 
  as_tibble()-> res5


res5 %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/median_estimate_20250326.xlsx")



# 
# 
# #####>>>>>>>>>>>>>>>>>>>>>>>>>>>混合中介分析<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 
# library(CMAverse)
# df_clean
# 
# Y <- df_clean$STL
# 
# A <-  df_clean$M189.PFAS
# 
# names(df_clean)
# 
# M189 <- df_clean$M189             
# M252 <- df_clean$M252              
# M259 <- df_clean$M259             
# M261 <- df_clean$M261             
# M310 <- df_clean$M310             
# M352 <- df_clean$M352           
# M355 <- df_clean$M355             
# M374 <- df_clean$M374            
# M376 <- df_clean$M376              
# M378 <- df_clean$M378             
# M400 <- df_clean$M400              
# M405  <-df_clean$M405            
# M423 <- df_clean$M423             
# M428 <- df_clean$M428              
# M450 <- df_clean$M450             
# M520 <- df_clean$M520              
# M567 <- df_clean$M567              
# M581 <- df_clean$M581              
# M607 <- df_clean$M607              
# M609 <- df_clean$M609             
# M642 <- df_clean$M642             
# M661 <- df_clean$M661              
# M668 <- df_clean$M668              
# M731 <- df_clean$M731            
# M87 <- df_clean$M87              
# 
# 
# 
# 
# as.numeric(df_clean$age_class_2)
# 
# # C1 <- as.numeric(df_clean$age_class_2)
# # C2 <- as.numeric(df_clean$BMI_class_2) 
# # C3 <- as.numeric(df_clean$education)
# # C4 <- as.numeric(df_clean$marriage_class_2)
# # C5 <- as.numeric(df_clean$income)
# # C6 <- as.numeric(df_clean$smk_class_2)
# # C7 <- as.numeric(df_clean$drk_class_2)
# # C8 <- as.numeric(df_clean$season)
# # C9 <- as.numeric(df_clean$abstinence_class_2)
# 
# 
# 
# 
# C1 <- df_clean$age_class_2
# C2 <- df_clean$BMI_class_2
# C3 <- df_clean$education
# C4 <- df_clean$marriage_class_2
# C5 <- df_clean$income
# C6 <- df_clean$smk_class_2
# C7 <- df_clean$drk_class_2
# C8 <- df_clean$season
# C9 <- df_clean$abstinence_class_2
# 
# 
# # "M189","M252","M259","M261","M310","M352","M355","M374","M376","M378", "M400","M405","M423","M428","M450","M520","M567","M581","M607","M609","M642","M661","M668","M731","M87"             
# 
# names(df_clean)
# 
# data <- data.frame(A, "M189","M252","M259","M261","M310","M352","M355","M374","M376","M378", "M400","M405","M423","M428","M450","M520","M567","M581","M607","M609","M642","M661","M668","M731","M87",Y,C1,C2,C3,C4,C5, C6,C7,C8,C9)
# 
# 
# cmdag(outcome = "Y", exposure = "A", mediator = c("M189","M252","M259","M261","M310","M352","M355","M374","M376","M378", "M400","M405","M423","M428","M450","M520","M567","M581","M607","M609","M642","M661","M668","M731","M87"),
#       basec = c("C1", "C2","C3", "C4","C5", "C6","C7", "C8","C9"), postc = NULL, node = TRUE, text_col = "white")
# 
# # mediator = c("M103","M196","M234","M248","M251","M252","M378","M400","M405","M423","M428","M445","M450","M520","M558","M581","M607","M609","M642","M661","M668","M731","M786","M787","M87","M91")
# # 
# # res_rb <- cmest(data =data, model = "rb", outcome = "Y", exposure = "A",mediator = c("M103","M196","M234","M248","M251","M252","M378","M400","M405","M423","M428","M445","M450","M520","M558","M581","M607","M609","M642","M661","M668","M731","M786","M787","M87","M91"),  basec = c("C1", "C2","C3", "C4","C5", "C6","C7", "C8","C9"), EMint = TRUE, mreg = list(rep("multinomial", length(mediator))), yreg = "liner",
# #                 astar = 0, a = 1,mval = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
# #                 estimation = "imputation", inference = "bootstrap", nboot = 10,boot.ci.type = "bca")
# # 
# # 
# # ?cmest
# # 
# 
# 
# 
# ?cmest
# 
# 
# 
# 
# ####>>>>>>>>>>>>>>构建循环函数，做每个暴露与混合中介的关系<<<<<<<<<<<<<<<<<<<<#####
# 
# 
# df_clean
# 
# lapply(seq(25, 49), function(i) {  
#   
# Y <- df_clean$STL
# 
# A <-  df_clean[[i]]   ###每次取出单独的暴露变量
# 
# names(df_clean)
# 
# M189 <- df_clean$M189             
# M252 <- df_clean$M252              
# M259 <- df_clean$M259             
# M261 <- df_clean$M261             
# M310 <- df_clean$M310             
# M352 <- df_clean$M352           
# M355 <- df_clean$M355             
# M374 <- df_clean$M374            
# M376 <- df_clean$M376              
# M378 <- df_clean$M378             
# M400 <- df_clean$M400              
# M405  <-df_clean$M405            
# M423 <- df_clean$M423             
# M428 <- df_clean$M428              
# M450 <- df_clean$M450             
# M520 <- df_clean$M520              
# M567 <- df_clean$M567              
# M581 <- df_clean$M581              
# M607 <- df_clean$M607              
# M609 <- df_clean$M609             
# M642 <- df_clean$M642             
# M661 <- df_clean$M661              
# M668 <- df_clean$M668              
# M731 <- df_clean$M731            
# M87 <- df_clean$M87              
# 
# 
# 
# C1 <- df_clean$age_class_2
# C2 <- df_clean$BMI_class_2
# C3 <- df_clean$education
# C4 <- df_clean$marriage_class_2
# C5 <- df_clean$income
# C6 <- df_clean$smk_class_2
# C7 <- df_clean$drk_class_2
# C8 <- df_clean$season
# C9 <- df_clean$abstinence_class_2
# 
# data <- data.frame(A,M189,M252,M259,M261,M310,M352,M355,M374,M376,M378,M400,M405,M423,M428,M450,M520,M567,M581,M607,M609,M642,M661,M668,M731,M87,Y,C1,C2,C3,C4,C5, C6,C7,C8,C9)
# 
# 
# res_rb <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
#                 mediator = c("M189","M252","M259","M261","M310","M352","M355","M374","M376","M378", "M400","M405","M423","M428","M450","M520","M567","M581","M607","M609","M642","M661","M668","M731","M87"),
#                 basec = c("C1","C2","C3","C4","C5","C6","C7","C8","C9"),
#                 EMint = TRUE, 
#                 mreg = list("linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear","linear"),
#                 yreg = "linear", 
#                 astar = 0, 
#                 a = 1,
#                 mval = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
#                 estimation = "imputation", 
#                 inference = "bootstrap", 
#                 nboot = 200,
#                 boot.ci.type = "bca")
# 
# summary(res_rb)$summarydf-> result   ###取出中介系数
# 
# rownames(result) %>% as_tibble() %>% 
#   cbind(result %>% as_tibble())  -> result_2 ###将中介结果的行名生成单独的一列变量
# 
# rep(names(df_clean)[i],15) %>% as_tibble() %>% 
#   dplyr::rename("term"="value") %>% cbind(result_2) -> result_3  ###将暴露变量的名称与中介结果结合在一起，整理成一个结果
# 
# 
# return(result_3)
# 
# 
# }) -> res13
# 
# 
# 
# res13 %>% 
#   bind_rows() %>% 
#   as_tibble()-> res14
# 
# 
# res14 %>% 
#   write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/only_exposure_mix_median_20241219_200.xlsx")
# 
# 
# 
# 
# # 
# # 
# # 
# # 1. CDE: 控制直接效应 (Controlled Direct Effect)
# # 定义: 控制直接效应 (CDE) 是指在控制了中介变量（如 M）的影响后，暴露（X）对结果（Y）的直接效应。也就是说，暴露对结果的影响在固定中介变量的情况下计算。
# # 解释: 量化在其他变量固定不变的情况下，暴露如何直接影响结果。
# # 2. PNDE: 纯自然直接效应 (Pure Natural Direct Effect)
# # 定义: 纯自然直接效应 (PNDE) 是指暴露对结果的直接效应，但不包括通过中介的间接效应。它代表了不经过中介的暴露与结果之间的直接关系。
# # 解释: 表示暴露对结果的直接影响，排除了通过中介的间接效应。
# # 3. TNDE: 总自然直接效应 (Total Natural Direct Effect)
# # 定义: 总自然直接效应 (TNDE) 是指通过中介变量发生的暴露对结果的直接效应，包括所有通过中介的直接和间接效应。
# # 解释: 量化暴露如何通过中介直接或间接影响结果。
# # 4. PNIE: 纯自然间接效应 (Pure Natural Indirect Effect)
# # 定义: 纯自然间接效应 (PNIE) 是指暴露通过中介变量对结果的间接效应，不包括任何交互作用。它是暴露与结果之间的间接路径，但排除了暴露与中介之间的交互作用。
# # 解释: 表示暴露如何通过中介变量间接影响结果，不考虑交互作用。
# # 5. TNIE: 总自然间接效应 (Total Natural Indirect Effect)
# # 定义: 总自然间接效应 (TNIE) 是指暴露通过中介变量对结果的所有间接效应，包括暴露与中介之间的交互作用。
# # 解释: 包括所有暴露通过中介的间接效应，同时考虑了暴露和中介之间的交互作用。
# # 6. TE: 总效应 (Total Effect)
# # 定义: 总效应 (TE) 是暴露对结果的整体影响，包含直接效应和所有间接效应（即直接效应和所有中介效应的总和）。
# # 解释: 量化暴露对结果的总影响，包括通过中介的间接效应和暴露的直接效应。
# # 7. Intref: 参考交互作用 (Reference Interaction)
# # 定义: 参考交互作用 (Intref) 是指暴露和中介之间的交互效应，通常是在一个参考水平下（如中介的基准值）计算的交互作用。
# # 解释: 研究暴露和中介之间交互作用对结果的影响，通常是在其他协变量固定的情况下进行分析。
# # 8. Intmed: 中介交互作用 (Mediated Interaction)
# # 定义: 中介交互作用 (Intmed) 是指暴露和中介变量之间的交互作用如何通过中介传递，从而影响结果。
# # 解释: 研究暴露和中介之间的交互作用如何影响结果，考虑了中介的作用。
# # 9. CDE(prop): CDE的比例 (Proportion CDE)
# # 定义: CDE的比例 (CDE[prop]) 是指暴露对结果的总效应中，由控制直接效应（CDE）所占的比例。
# # 解释: 量化暴露对结果的总效应中，多少部分是由控制直接效应所解释的。
# # 10. Intref(prop): Intref的比例 (Proportion Intref)
# # 定义: Intref的比例 (Intref[prop]) 是指暴露与中介的参考交互作用对总效应的贡献比例。
# # 解释: 量化暴露与中介的参考交互作用对结果的总影响的贡献。
# 11. Intmed(prop): Intmed的比例 (Proportion Intmed)
# 定义: Intmed的比例 (Intmed[prop]) 是指通过中介的交互作用对总效应的贡献比例。
# 解释: 量化暴露与中介之间的交互作用通过中介影响结果的比例。
# 12. PNIE(prop): PNIE的比例 (Proportion PNIE)
# 定义: PNIE的比例 (PNIE[prop]) 是指暴露通过中介对结果的纯自然间接效应（PNIE）所占总效应的比例。
# 解释: 量化暴露通过中介对结果的间接效应，在总效应中所占的比例。
# 13. PM: 总体中介比例 (Overall Proportion Mediated)
# 定义: 总体中介比例 (PM) 是指通过中介的效应占总效应的比例，包括所有直接和间接效应。
# 解释: 量化总效应中有多少部分是通过中介传递的。
# 14. Int: 总体交互作用比例 (Overall Proportion Attributable to Interaction)
# 定义: 总体交互作用比例 (Int) 是指暴露和中介之间的交互作用在总效应中的比例。
# 解释: 量化暴露与中介之间的交互作用在总效应中所占的比例。
# 15. PE: 总体消除比例 (Overall Proportion Eliminated)
# 定义: 总体消除比例 (PE) 是指如果暴露效应被消除或修改，总效应中会被消除的部分所占的比例。
# 解释: 量化如果暴露不存在或被减少到零，总效应将减少多少。
# 总结：
# 直接效应：CDE、PNDE、TNDE 分别表示暴露对结果的不同类型的直接效应。
# 间接效应：PNIE 和 TNIE 分别表示暴露通过中介对结果的不同类型的间接效应。
# 总效应：TE 是暴露对结果的总体效应，包含直接和间接效应。
# 交互作用效应：Intref 和 Intmed 表示暴露与中介之间不同形式的交互作用效应。
# 比例效应：如 CDE(prop)、Intref(prop) 和 PM 分别表示各个效应在总效应中所占的比例。
# PE 表示暴露效应消除后的总体变化。
# 


# 
# 
# res_rb <- cmest(data = df_clean, model = "rb", outcome = "Y", exposure = "A",
#                 mediator = c("M103", "M196"),
#                 basec = c("C1", "C2"),
#                 EMint = TRUE, mreg = list("multinomial", "multinomial"),
#                 yreg = "linear", astar = 0, a = 1,
#                 mval = list(0,0),
#                 estimation = "imputation", inference = "bootstrap",
#                 nboot = 10, boot.ci.type = "bca")
# 
# 
# 
# res_rb <- cmest(data = df_clean, model = "rb", outcome = "Y", exposure = "A",
#                 mediator = c("M103", "M196"),
#                 basec = c("C1", "C2"),
#                 EMint = TRUE, mreg = list("linear","linear"),  # Use linear regression for mediators
#                 yreg = "linear", astar = 0, a = 1,
#                 mval = list(0, 0),
#                 estimation = "imputation", inference = "bootstrap",
#                 nboot = 10, boot.ci.type = "bca")
# 
# 


# 打印结果以查看
# print(result)

# 
# ###然后再将生成的list重新保存为独立的excel文件
# lapply(seq_along(res3), function(i) {
#   write.xlsx(x = res3[[i]] %>% bind_rows(), file = paste0("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/median_estimate_",names(res3) %>% as_tibble() %>% slice(i), ".xlsx"), row.names = FALSE)
# })

# 
# 
# lapply(seq(26, 53), function(b) {
#   m <- 27
#   x <- 27 + 28
#   
#   # 检查索引是否有效
#   if (x <= length(df_clean)) {
#     print(paste("Processing column:", names(df_clean)[x]))
#   } else {
#     stop("Invalid index for x")
#   }
#   
#   # 使用 as_tibble 之前，确保列的索引正确
#   lapply(df_clean[[names(df_clean)[x]]], function(x_name) {   #### 中介变量
#     if (length(x_name) == 0) {
#       stop("Invalid column for x_name")
#     }
#     
#     lapply(df_clean[[names(df_clean)[27]]], function(m1) {   ### 因变量
#       if (length(m1) == 0) {
#         stop("Invalid column for m1")
#       }
#       
#       x <- df_clean[[x_name]]
#       medi_m_3 <- extract_fitM(x %>% as_tibble, df_clean, x_name %>% as_tibble)
#       medi_Y_3 <- extract_fitY(x %>% as_tibble, df_clean[[names(df_clean)[27]]] %>% as_tibble, df_clean, x_name %>% as_tibble)
#       
#       combined_result <- rbind(medi_m_3, medi_Y_3)
#       return(combined_result)
#     })
#   })
# })
# 
# 
# 
# 

# 
# 
# lapply(df_clean[[names(df_clean)[55]]] %>% as_tibble, function(x_name) {   ####中介变量
#   
#   lapply(df_clean[[names(df_clean)[27]]] %>% as_tibble, function(m1) {    ###因变量
#     
#     x <- df_clean[[55]]
#     medi_m_3 <- extract_fitM( df_clean[[55]],df_clean,df_clean[[names(df_clean)[55]]] %>% as_tibble)
#     medi_Y_3 <- extract_fitY(df_clean[[55]]%>% as_tibble, df_clean[[names(df_clean)[27]]] ,df_clean, df_clean[[names(df_clean)[55]]] %>% as_tibble)
#     
#     combined_result <- rbind(medi_m_3, medi_Y_3)
#     return(combined_result)
# 


