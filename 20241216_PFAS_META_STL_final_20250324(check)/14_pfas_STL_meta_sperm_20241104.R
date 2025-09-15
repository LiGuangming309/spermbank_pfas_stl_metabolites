
no_source()
rm(list = ls())
# setwd(r4projects::get_project_wd())

# devtools::install_github("jaspershen/r4projects")


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

df_background_DNA_LT

df_background_DNA_LT %>% 
  dplyr::select(number,age,age_class_2,BMI_class_2,education,marriage_class_2,income,smk_class_2,drk_class_2,season,abstinence_class_2,mtDNAcn,STL) -> df_background_DNA_LT_tidy
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


####对PFAS取对数
lapply(df_background_DNA_LT_sperm_PFOA_tidy[,c(1:14)],log) ->df_background_DNA_LT_sperm_PFOA_tidy_lg 

names(df_background_DNA_LT_sperm_PFOA_tidy_lg) 

### PFAS和DNA均取对数了
df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,26,27)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)









###====================读入代谢组学数据===============================

df_meta<- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/240823_最终确定版data.xlsx",sheet = "已核对")
str(df_meta) 


df_meta %>% 
  dplyr::rename(number="ID") -> df_meta_1

str(df_meta_1)

names(df_meta_1)

####>>>>>>>删除外源或二肽物质
df_out <- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/外源_20241020.xlsx",sheet = "外源和重复")

df_out %>% 
  dplyr::rename(metabolite="META")  -> df_out_2

str(df_out_2)   ###一共有266个二肽或外源物质


df_out_2$metabolite



df_meta_1 %>% 
  select(-c(df_out_2$metabolite))  ->  df_meta_2  ####最终得到414个代谢物质

str(df_meta_2)

794-437

437-23


###将剩余的物质与PFAS和STL进行合并

df_background_DNA_LT_sperm_PFOA_tidy_final_lg  %>% 
  left_join(df_meta_2,by="number") -> df_background_DNA_LT_sperm_PFOA_meta_tidy

colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy))




####删除缺失值

df_background_DNA_LT_sperm_PFOA_meta_tidy  %>% 
  dplyr::filter(!sampleID=="NA")  %>% 
  as_tibble()->df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA  ###得到的精浆和人群基线数据一样的数据



colSums(is.na(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA))

df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA

names(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA)

str(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA)

# write.xlsx(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/result_pfas_rawdata.xlsx")  ###此处的pfas和端粒长度均取对数了，并未进行标化或归一化



####>>>>>>>>>>>>>>>提取PFAS和STL都是P值<0.05,出匹配后的代谢物质<<<<<<<<<<<<<<<<<<<###########################




####>>>>>>>>读入代谢物匹配的代谢通路变量名称#############################


df_kegg<- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/240823_最终确定版data_lgm.xlsx",sheet = "合并")
str(df_kegg) 



df_kegg %>% 
  dplyr::rename(metabolite="META") -> df_kegg_1

str(df_kegg_1)

colSums(is.na(df_kegg_1))

####>>>>>>>>读入混合物计算后的P值------------------------------
# write.xlsx(results_df_mix_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/results_df_mix_p_value.xlsx")
pfas_mix_p<- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/pfas_STL_merge/results_df_mix_p_value.xlsx")
pfas_mix_p

pfas_mix_p%>% 
  left_join(df_kegg_1,by="metabolite") -> pfas_mix_p_kegg


str(pfas_mix_p_kegg)

names(pfas_mix_p_kegg)

# ####>>>>>>>删除外源或二肽物质
# df_out <- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/外源_20241020.xlsx",sheet = "Sheet2")
# 
# df_out %>% 
#   dplyr::rename(metabolite="META")  -> df_out_2
# 
# df_out_2$group
# 
# pfas_mix_p_kegg %>% 
#   left_join(df_out_2,by="metabolite") %>% 
#   dplyr::filter(is.na(group)) ->pfas_mix_p_kegg_final 
# 
# # pfas_mix_p_kegg_final$group
# names(pfas_mix_p_kegg_final)
# str(pfas_mix_p_kegg_final)


###提取出混合暴露物质与代谢物有关联的物质
pfas_mix_p_kegg %>% 
  dplyr::filter(!q_value.x=="NA")  -> pfas_mix_p_kegg_meta 

pfas_mix_p_kegg_meta$q_value.x
# 
# write.csv(pfas_mix_p_kegg_meta,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/pfas_STL_merge/pfas_mix_meta_fdr_20241019.csv")

names(pfas_mix_p_kegg_meta)


pfas_mix_p_kegg_meta$metabolite



df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA %>% 
  select(number:sampleID,pfas_mix_p_kegg_meta$metabolite) -> df_pfas_stl_p

names(df_pfas_stl_p)

str(df_pfas_stl_p)
# 
# df_pfas_stl_p%>% 
#   select(sampleID,M10:M794) -> df_meta_final


# df_meta_final



###>>>>>>>>>>>>>>>>>>>读入精液质量参数<<<<<<<<<<<<<<<<<<<<<<<<################# 



###=================读入重复测量的精液质量参数=======================##########


df2 <- read_excel("raw_data/semen_quality_parameters.xlsx",sheet = "匹配表")



df2 %>%  
  dplyr::select(1:29) -> df2_1  ##每一次精液质量参数与最后一次精液参数匹配后的数据
str(df2_1)
length(unique(df2_1$number...1 ))   #1487份人群精液样本，重复测量6608次


# df2_1[complete.cases(df2_1),] -> df2_1_1  #删除缺失值
# write.xlsx(df2_1,"result/20231211_result/20240213result/semen_quality_parameters_mtDNAcn_match.xlsx") #导出筛选的变量  


###将重复测量中精液质量参数异常结果剔除
df2_1 %>%  
  dplyr::mutate(group_semen=case_when(
    motilityd...13==100|motilityd...13<0~1,   ###1代表精液质量数据有异常的数据，0代表精液质量数据无异常值
    TRUE~0
  )) -> df2_1_1

# write.xlsx(df2_1_1,"result/20231211_result/20240213result/semen_quality_parameters(将有问题的精液质量筛检生成分组).xlsx") ##导出筛选的变量  

str(df2_1_1)

df2_1_1 %>% 
  dplyr::filter(group_semen==0) -> df2_1_2  ###将

str(df2_1_2) 

df2_1_2$density...8[which(df2_1_2$density...8 ==0)] <- 0.01

# write.xlsx(df2_1_2,"result/20231211_result/20240213result/semen_quality_parameters(将有问题的生成分组).xlsx") ##导出筛选的变量  

str(df2_1_2)
length(unique(df2_1_2$number...1))
# colSums()
# 
# 
# 
# df2 %>%  
#   dplyr::select(42:59)  -> df2_2  ##最后一次检测mtDNAcn的精液质量参数
# 
# df2_2[complete.cases(df2_2),] -> df2_2_1  #删除缺失值
# str(df2_2_1)            #1164份含有mtDNA样本
# # write.xlsx(df2_2_1,"result/20231211_result/20240213result/semen_quality_parameters_mtDNAcn.xlsx") #导出筛选的变量  
# 


###=================数据集合并以及变量筛选================================

###==将一般资料的数据与精液质量参数进行匹配
df2_1_2 %>% 
  dplyr::rename(number=number...1) -> df2_1_3 
str(df2_1_3)

colSums(is.na(df2_1_3))

str(df_pfas_stl_p)

length(unique(df_pfas_stl_p$number)) 

df_pfas_stl_p %>%
  left_join(df2_1_3,by="number")  -> df_background_repeatsemen

df_background_repeatsemen 
names( df_background_repeatsemen)
# colSums(is.na(df_background_repeatsemen)) 
str(df_background_repeatsemen)
# 
length(unique(df_background_repeatsemen$mtDNAcn...17)) #匹配得到1164个男性，重复捐精5739个样本
# 

length(unique(df_background_repeatsemen$mtDNAcn...17)) 


df_background_repeatsemen %>% 
  dplyr::filter(!mtDNAcn...17=="NA") ->df_background_repeatsemen_noNA 

str(df_background_repeatsemen_noNA)   ###1164份

length(unique(df_background_repeatsemen_noNA$number))

colSums(is.na(df_background_repeatsemen_noNA)) 





# 
# write.xlsx(df_background_repeatsemen,"result/df_background_repeatsemen.xlsx")

# 
# df_background_repeatsemen %>% 
#   group_by(number) %>%
#   plot_time_series(date, confirmed, .facet_ncol = 2, .interactive = FALSE) + 
#   theme_tq(base_family = cnfont)
# 
# 

# df <- read_excel("raw_data/data1.xlsx",sheet = "匹配表")

#将样本中重复的dateaq重新命名

df_background_repeatsemen_noNA %>% 
  mutate(datesq=lubridate::ymd(datesq)) %>% 
  arrange(number, datesq) %>% 
  group_by(number) %>% 
  mutate(firstday = min(datesq)) %>% 
  mutate(days = datesq - firstday) %>%  #计算每个ID下每次随访时间与基线时间间隔的天数
  ungroup()  -> df_background_repeatsemen_days   #得到含有时间天数的变量


# write.xlsx(df_background_repeatsemen_days,"result/20231211_result/20240213result/df_background_repeatsemen_days.xlsx")


str(df_background_repeatsemen_days)

# as.numbers(df_background_repeatsemen_days$days)

df_background_repeatsemen_days$days <- as.numeric(df_background_repeatsemen_days$days) #将时间变量转换为数字变量


####根据随访与基线时间间隔的天数，按精子发生时间，划分间隔
df_background_repeatsemen_days %>%  
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

# write.xlsx(df_background_repeatsemen_final1,"result/20231211_result/20240213result/df_background_repeatsemen_final1.xlsx")

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
# write.xlsx(df_background_repeatsemen_final2,"result/20231211_result/20240213result/df_background_repeatsemen_final2.xlsx")


colSums(is.na(df_background_repeatsemen_final2)) 

#### ==============提取出重复测量的精液质量参数数据的基线期数据================

###提取基线数据  
df_background_repeatsemen_final2 %>% 
  dplyr::filter(interval==0) -> df_background_repeatsemen_base_lifestyle

str(df_background_repeatsemen_base_lifestyle)

colSums(is.na(df_background_repeatsemen_base_lifestyle))
table(df_background_repeatsemen_base_lifestyle$interval)

# write.xlsx(df_background_repeatsemen_base_lifestyle,"result/20231211_result/20240213result/df_background_repeatsemen_base_lifestyle.xlsx")


# str(df_background_repeatsemen_base_lifestyle)


df_background_repeatsemen_final2


names(df_background_repeatsemen_final2)


###生成基线精液质量参数
df_background_repeatsemen_final2 %>% 
  dplyr::mutate(sperm_count=volume...7*density...8,
                Total_motility=motilitya...10+motilityb...11+motilityc...12,
                Progressive_motility=motilitya...10+motilityb...11) -> df_background_pfas_stl_sperm

###对精液质量参数取对数
df_background_pfas_stl_sperm %>% 
  dplyr::mutate(
                lg_volume=log(volume...7),
                lg_density=log(density...8),
                lg_sperm_count=log(sperm_count),
                lg_total_motility=log(Total_motility),
                lg_progressive_motility=log(Progressive_motility)) ->df_background_pfas_stl_sperm_lg 

str(df_background_pfas_stl_sperm_lg)

names(df_background_pfas_stl_sperm_lg)


# adjusted_y =
#   lme4::lmer(
#     formula = y ~ age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2+ Ethnicity + (1 | number),
#     data = temp_data
#   )



lapply(df_background_pfas_stl_sperm_lg[,c(213:217)], function(y){
  lapply(df_background_pfas_stl_sperm_lg[,c(48:176)], function(x){
    lmer(y~x+age_class_2+BMI_class_2+education+marriage_class_2+income.x+smk_class_2+drk_class_2+season.x+abstinence_class_2+(1|number),data=df_background_pfas_stl_sperm_lg) -> result
    summary(result) -> res1
    res1$coefficients -> res2
    β <- res2[,1]
    SE <- res2[,2]
    CI2.5 <- β-1.96*SE
    CI97.5  <- β+1.96*SE
    CI95<-paste0(round(β,2)," (",round(CI2.5,2),', ',round(CI97.5,2),")")
    percent_β <- (exp(β)-1)*100
    percent_CI2.5 <- (exp(β-1.96*SE)-1)*100
    percent_CI97.5 <- (exp(β+1.96*SE)-1)*100
    percent_CI95 <- paste0(round(percent_β,2)," (",round(percent_CI2.5,2),', ',round(percent_CI97.5,2),")")
    res2 %>% 
      rownames() %>% 
      as_tibble() -> res ####自己添加的，将其他输出的变量也添加到包里
    Uni_glm_model <- data.frame('Estimate'=res2[,1],
                                "Std. Error"=res2[,2],
                                'df' = res2[,3],
                                't value' = res2[,4],
                                'p value' = res2[,5],
                                'CI2.5' = CI2.5,
                                'CI97.5' = CI97.5,
                                'CI95' = CI95,
                                "percent_β"=percent_β,
                                "percent_CI2.5"=percent_CI2.5,
                                "percent_CI97.5"=percent_CI97.5,
                                "percent_CI95"=percent_CI95,
                                res)
    return(Uni_glm_model)})
}) -> res1  ###将跑的循环中y生成的list结果保存为res1



lapply(res1, bind_rows)

lapply(res1, bind_rows) -> df_list



###然后再将生成的list重新保存为独立的excel文件
lapply(seq_along(df_list), function(i) {
  write.xlsx(x =df_list[[i]], file = paste0("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/metabolomic_sperm/result_meta_sperm_0628_MI_",names(res1) %>% as_tibble() %>% slice(i), ".xlsx"), row.names = FALSE)
})






###############>>>>>>>>1.提取精液体积<<<<<<<<<<<<<<<<<<<<<<#######################
####1.提取精液体积
df_list$lg_volume -> df_volume

df_volume %>% 
  filter(value=="x") -> df_volume



colnames(df_background_pfas_stl_sperm_lg[,c(48:176)]) %>%   ###提取y变量的名称，并重新转换为数据框
  as_tibble() %>% 
  dplyr::rename("variable"="value")  -> df_volume_1

df_volume_1 %>% 
  cbind(df_density)  %>% 
  select(-value) %>% 
  cbind(rep("Volume",129) %>% as_tibble()) %>% 
  dplyr::rename("term"="value")   -> df_volume_2

str(df_volume_2)

###BH校正的p值
p.adjust(df_volume_2$p.value,method = "BH") %>% 
  as_tibble()  %>% 
  dplyr::rename("p_adjust"="value") ->df_volume_3 

df_volume_2 %>% 
  cbind(df_volume_3) -> df_volume_final



####2.提取精液体积
df_list$lg_density-> df_concentration

df_density %>% 
  filter(value=="x") -> df_concentration



colnames(df_background_pfas_stl_sperm_lg[,c(48:176)]) %>%   ###提取y变量的名称，并重新转换为数据框
  as_tibble() %>% 
  dplyr::rename("variable"="value")  -> df_concentration_1

df_concentration_1 %>% 
  cbind(df_density)  %>% 
  select(-value) %>% 
  cbind(rep("Concentration",129) %>% as_tibble()) %>% 
  dplyr::rename("term"="value")   -> df_concentration_2 

str(df_concentration_2)

###BH校正的p值
p.adjust(df_concentration_2$p.value,method = "BH") %>% 
  as_tibble()  %>% 
  dplyr::rename("p_adjust"="value") ->df_concentration_3 

df_concentration_2 %>% 
  cbind(df_concentration_3) -> df_concentration_final





####3.提取精子计数
df_list$lg_sperm_count -> df_sperm_count

df_sperm_count %>% 
  filter(value=="x") -> df_sperm_count



colnames(df_background_pfas_stl_sperm_lg[,c(48:176)]) %>%   ###提取y变量的名称，并重新转换为数据框
  as_tibble() %>% 
  dplyr::rename("variable"="value")  -> df_sperm_count_1

df_sperm_count_1 %>% 
  cbind(df_sperm_count)  %>% 
  select(-value) %>% 
  cbind(rep("Sperm count",129) %>% as_tibble()) %>% 
  dplyr::rename("term"="value")   -> df_sperm_count_2 

str(df_sperm_count_2)

###BH校正的p值
p.adjust(df_sperm_count_2$p.value,method = "BH") %>% 
  as_tibble()  %>% 
  dplyr::rename("p_adjust"="value") ->df_sperm_count_3 

df_sperm_count_2 %>% 
  cbind(df_sperm_count_3) -> df_sperm_count_final





####4.提取精液总活力
df_list$lg_total_motility -> df_total_motility

df_total_motility %>% 
  filter(value=="x") -> df_total_motility




colnames(df_background_pfas_stl_sperm_lg[,c(48:176)]) %>%   ###提取y变量的名称，并重新转换为数据框
  as_tibble() %>% 
  dplyr::rename("variable"="value")  -> df_total_motility_1

df_total_motility_1 %>% 
  cbind(df_total_motility)  %>% 
  select(-value) %>% 
  cbind(rep("Total motility",129) %>% as_tibble()) %>% 
  dplyr::rename("term"="value")   -> df_total_motility_2 


str(df_total_motility_2)

###BH校正的p值
p.adjust(df_total_motility_2$p.value,method = "BH") %>% 
  as_tibble()  %>% 
  dplyr::rename("p_adjust"="value") ->df_total_motility_3 

df_total_motility_2 %>% 
  cbind(df_total_motility_3) -> df_total_motility_final




####5.提取精液前向运动活力

df_list$lg_progressive_motility -> df_progressive_motility

df_progressive_motility %>% 
  filter(value=="x") -> df_progressive_motility



colnames(df_background_pfas_stl_sperm_lg[,c(48:176)]) %>%   ###提取y变量的名称，并重新转换为数据框
  as_tibble() %>% 
  dplyr::rename("variable"="value")  -> df_progressive_motility_1
  
df_progressive_motility_1 %>% 
  cbind(df_progressive_motility)  %>% 
  select(-value) %>% 
  cbind(rep("Progressive motility",129) %>% as_tibble()) %>% 
  dplyr::rename("term"="value")   -> df_progressive_motility_2


str(df_progressive_motility_2)

###BH校正的p值
p.adjust(df_progressive_motility_2$p.value,method = "BH") %>% 
  as_tibble()  %>% 
  dplyr::rename("p_adjust"="value") ->df_progressive_motility_3 

df_progressive_motility_2 %>% 
  cbind(df_progressive_motility_3) -> df_progressive_motility_final





####>6.将最终的数据全部合并起来

rbind(df_volume_final,df_concentration_final,df_sperm_count_final,df_total_motility_final,df_progressive_motility_final)  -> df_meta_sperm


####>>>>>7.读入代谢物分类图 


df_meta_class<- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/代谢产物的分类.xlsx",sheet = "Sheet2")
str(df_meta_class) 

names(df_meta_class)

names(pfas_mix_p_kegg)
pfas_mix_p_kegg

colSums(is.na(df_meta_class))

###添加代谢物的分类
colSums(is.na(pfas_mix_p_kegg))
pfas_mix_p_kegg %>% 
  left_join(df_meta_class,by="name") -> pfas_mix_p_kegg_class
 
colSums(is.na(pfas_mix_p_kegg_class))

names(pfas_mix_p_kegg_class)


write.xlsx(pfas_mix_p_kegg_class,"raw_data/Matebolomics_sperm/20240823_haixia/代谢产物的分类_2tidy.xlsx")
# 使用 ggalluvial 绘制桑基图
library(ggplot2)
library(ggalluvial)
library(hrbrthemes)

# 设置字体
library(showtext) 
showtext_auto(enable = TRUE) 
# font_add("songti", 
#          regular = "song.otf",
#          bold = "song.otf",
#          italic = "song.otf",
#          bolditalic = "song.otf")

# 设置 ggplot2 绘图主题
# theme_set(theme_ipsum(base_family = "songti"))


library(ggsci)
library(tidyverse)
library(cowplot)
library("scales") 

show_col(pal_nejm("default", alpha = 1)(8))  ####查看NEJM杂志的配色，并进行着色 

pfas_mix_p_kegg_class %>% 
  dplyr::rename(variable="metabolite") -> pfas_mix_p_kegg_meta

df_meta_sperm %>% 
  left_join(pfas_mix_p_kegg_meta,by="variable") -> df_meta_sperm_2 

names(df_meta_sperm_2)

colSums(is.na(df_meta_sperm_2))

pg <-df_meta_sperm_2 %>% 
  dplyr::count(name, term, class,Estimate)

df_meta_sperm_2$term

ggplot(pg, aes(
  axis1 =name, axis2 = class,
  axis3 = term, y = Estimate
), size = 0.001) +
  geom_stratum(width = 0.5, alpha = 0.2, size = 0.1) + 
  geom_alluvium(aes(fill = term), width = 0.5)+
  scale_fill_manual(values = c(
    "Volume" = "#7876B1FF",
    "Concentration" = "#20854EFF",
    "Sperm count" = "#E18727FF",
    "Progressive motility" = "#0072B5FF",
    "Total motility" = "#BC3C39FF"
    )) + 
  geom_text(
    stat = "stratum",
    infer.label = TRUE,
    family = "Arial", size = 3.5,
    color = "#2A2A2A"
  ) + 
  scale_x_continuous(
    breaks = 1:3,
    labels = c("Metabolite", "Class", "Sperm")
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) 




ggplot(pg, aes(
  axis1 =name, axis2 = class,
  axis3 = term, y = Estimate
), size = 0.001) +
  geom_stratum(width = 0.5, alpha = 0.2, size = 0.1) + 
  geom_alluvium(aes(fill = term), width = 0.5)+
  scale_fill_manual(values = c(
    "Volume" = "#7876B1FF",
    "Concentration" = "#20854EFF",
    "Sperm count" = "#E18727FF",
    "Progressive motility" = "#0072B5FF",
    "Total motility" = "#BC3C39FF"
  )) +
  ggfittext::geom_fit_text(
    stat = "stratum", 
    infer.label = TRUE, 
    family = "Arial", min.size = 0.1,
    color = "#2A2A2A"
  ) +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("Metabolite", "Class", "Sperm")
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) 
  
p_meta
scale_y_continuous(breaks = breaks)


library(ggplot2)
library(ggalluvial)
library(ggfittext)

ggplot(pg, aes(
  axis1 = name, axis2 = class,
  axis3 = as.factor(term), y = Estimate
)) +
  geom_stratum(width = 0.6, alpha = 0.6, fill = "#EAEAEA") + 
  geom_alluvium(aes(fill = as.factor(term)), width = 0.6, alpha = 0.9) +
  scale_fill_manual(values = c(
    "Volume" = "#7876B1FF",
    "Concentration" = "#20854EFF",
    "Sperm count" = "#E18727FF",
    "Progressive motility" = "#0072B5FF",
    "Total motility" = "#BC3C39FF"
  )) +
  ggfittext::geom_fit_text(
    stat = "stratum", 
    infer.label = TRUE, 
    family = "Arial", min.size = 0.1,
    color = "#2A2A2A"
  ) +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("Metabolite", "Class", "Sperm")
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(), # 移除网格线
    axis.title = element_blank(), # 去掉轴标题
    axis.text.y = element_blank(), # 去掉y轴刻度
    axis.ticks.y = element_blank(), # 去掉y轴刻度线
    # legend.title = element_text(size = 14, family = "Arial", face = "bold"),
    legend.text = element_text(size = 14, family = "Arial", face = "bold"),
    axis.text.x = element_text(size = 14, family = "Arial", face = "bold", color = "black"), # 设置x轴文字大小和颜色
    legend.position = "bottom", # 移动图例到底部
    legend.title = element_blank(), # 去掉图例标题
    plot.background = element_rect(fill = "#F5F5F5", color = NA), # 设置背景颜色
    panel.background = element_rect(fill = "#F5F5F5", color = NA) # 面板背景色
  ) -> p_meta

p_meta


cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/metabolomic_sperm/pfas_meta_class.pdf",width = 16, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

p_meta # 使用大写罗马表示
par(opar)
dev.off()






####






# # no_source()
# rm(list = ls())
# setwd(r4projects::get_project_wd())
# source("R/100-tools.R")
# 
# library(tidyverse)
# library(tidymass)
# 
# 
# 
# pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
#                tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
#                rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
#                ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts)
# # get_project_wd()
# # masstools::setwd_project()
# 
# getwd() 
# 
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
# ####>>>>>>>>读入混合物计算后的P值------------------------------
# # write.xlsx(results_df_mix_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/results_df_mix_p_value.xlsx")
# pfas_mix_p<- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/pfas_STL_merge/results_df_mix_p_value.xlsx")
# pfas_mix_p
# 
# pfas_mix_p%>% 
#   left_join(df_kegg_1,by="metabolite") -> pfas_mix_p_kegg
# 
# 
# str(pfas_mix_p_kegg)
# 
# names(pfas_mix_p_kegg)
# 
# # ####>>>>>>>删除外源或二肽物质
# # df_out <- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/外源_20241020.xlsx",sheet = "Sheet2")
# # 
# # df_out %>% 
# #   dplyr::rename(metabolite="META")  -> df_out_2
# # 
# # df_out_2$group
# # 
# # pfas_mix_p_kegg %>% 
# #   left_join(df_out_2,by="metabolite") %>% 
# #   dplyr::filter(is.na(group)) ->pfas_mix_p_kegg_final 
# # 
# # # pfas_mix_p_kegg_final$group
# # names(pfas_mix_p_kegg_final)
# # str(pfas_mix_p_kegg_final)
# 
# 
# ###提取出混合暴露物质与代谢物有关联的物质
# pfas_mix_p_kegg %>% 
#   dplyr::filter(!q_value.x=="NA")  -> pfas_mix_p_kegg_meta 
# 
# pfas_mix_p_kegg_meta$q_value.x
# 
# write.csv(pfas_mix_p_kegg_meta,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/pfas_STL_merge/pfas_mix_meta_fdr_20241019.csv")
# 
# names(pfas_mix_p_kegg_meta)
# 
# 
# pfas_mix_p_kegg_meta$metabolite
# 
# pfas_mix_p_kegg_meta
# 
# 
# 
# 
# 










