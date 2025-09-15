# no_source() 
pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel)
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




###------------------读入血浆PFOA的信息-----------------------

df_blood_PFOA <- read_excel("raw_data/Sperm_data_share_lgm/血浆-精浆PFAS数据-final-20240524.xlsx",sheet = "血浆PFAS")

str(df_blood_PFOA)

###删除前2行数据
df_blood_PFOA %>% 
  slice(-c(1:2)) ->df_blood_PFOA_tidy 

###将编码变量重新命名
df_blood_PFOA_tidy  %>% 
  dplyr::rename(number=`Plasma PFAS`) ->df_blood_PFOA_tidy_2
str(df_blood_PFOA_tidy_2)

df_blood_PFOA_tidy_2$number <- as.numeric(df_blood_PFOA_tidy_2$number)


# str(df_blood_PFOA_tidy_2)

names(df_blood_PFOA_tidy_2)

df_blood_PFOA_tidy_2 %>% 
  dplyr::rename(L_PFHxS="L-PFHxS",
                diPAP_6_2="6:2diPAP",
                PF3OUdS_11CL="11Cl-PF3OUdS",
                PF3ONS_9CL="9Cl-PF3ONS",
                L_PFHpS="L-PFHpS",
                L_PFBS="L-PFBS",
                L_PFDS="L-PFDS",
                FOSA_I="FOSA.I"
  ) -> df_blood_PFOA_tidy_3


#### 将PFOA检测限LOD/根号2替代缺失值
df_blood_PFOA_tidy_3 %>% 
  mutate(
    PFOA=case_when(
      PFOA<0.00675113081441142~0.00675113081441142/sqrt(2),
      TRUE~ PFOA),
    PFOS=case_when(
      PFOS<0.0516440006885867~0.0516440006885867/sqrt(2),
      TRUE~ PFOS),
    PFNA=case_when(
      PFNA<0.045045045045045~0.045045045045045/sqrt(2),
      TRUE~ PFNA),
    PFDA=case_when(
      PFDA<0.00397719740156436~0.00397719740156436/sqrt(2),
      TRUE~ PFDA),
    PFUdA=case_when(
      PFUdA<0.00351658656663932~0.00351658656663932/sqrt(2),
      TRUE~ PFUdA),
    L_PFHxS=case_when(
      L_PFHxS<0.00338104361546264~0.00338104361546264/sqrt(2),
      TRUE~ L_PFHxS),
    diPAP_6_2=case_when(
      diPAP_6_2<0.00221141088014153~0.00221141088014153/sqrt(2),
      TRUE~ diPAP_6_2),
    PF3OUdS_11CL=case_when(
      PF3OUdS_11CL<0.000540365286933967~0.000540365286933967/sqrt(2),
      TRUE~ PF3OUdS_11CL),
    PF3ONS_9CL=case_when(
      PF3ONS_9CL<0.000503330369276714~0.000503330369276714/sqrt(2),
      TRUE~ PF3ONS_9CL),
    PFDoA=case_when(
      PFDoA<0.00255493101686255~0.00255493101686255/sqrt(2),
      TRUE~ PFDoA),
    L_PFHpS=case_when(
      L_PFHpS<0.00917992656058752~0.00917992656058752/sqrt(2),
      TRUE~ L_PFHpS),
    PFTrDA=case_when(
      PFTrDA<0.00586969281940912~0.00586969281940912/sqrt(2),
      TRUE~ PFTrDA),
    L_PFBS=case_when(
      L_PFBS<0.00494967827091239~0.00494967827091239/sqrt(2),
      TRUE~ L_PFBS),
    PFBA=case_when(
      PFBA<0.0319216854649925~0.0319216854649925/sqrt(2),
      TRUE~ PFBA),
    L_PFDS=case_when(
      L_PFDS<0.00312044934470564~0.00312044934470564/sqrt(2),
      TRUE~ L_PFDS),
    FOSA_I=case_when(
      FOSA_I<0.00151095441954168~0.00151095441954168/sqrt(2),
      TRUE~ FOSA_I)
  )-> df_blood_PFOA_tidy_4

str(df_blood_PFOA_tidy_4)


####将PFOS数据与基础数据合并
df_background_DNA_LT__tidy_lg  %>% 
  left_join(df_blood_PFOA_tidy_4,by="number") -> df_background_DNA_LT_blood_PFOA 

str(df_background_DNA_LT_blood_PFOA)

summary(df_background_DNA_LT_blood_PFOA)
colSums(is.na(df_background_DNA_LT_blood_PFOA)) 


####删除缺失值

df_background_DNA_LT_blood_PFOA %>% 
  dplyr::filter(!PFOA=="NA") -> df_background_DNA_LT_blood_PFOA_noNA 

str(df_background_DNA_LT_blood_PFOA_noNA)   ##

length(unique(df_background_DNA_LT_blood_PFOA_noNA$number))  ###最终得到1021人

colSums(is.na(df_background_DNA_LT_blood_PFOA_noNA)) 


####最终数据

str(df_background_DNA_LT_blood_PFOA_noNA)






###-----将数据重新排序整理---------------
df_background_DNA_LT_blood_PFOA_noNA %>% 
  select(PFOA:FOSA_I,mtDNAcn,STL,everything()) -> df_background_DNA_LT_blood_PFOA_tidy



str(df_background_DNA_LT_blood_PFOA_tidy)




####对PFAS取对数
lapply(df_background_DNA_LT_blood_PFOA_tidy[,c(1:18)],log) ->df_background_DNA_LT_blood_PFOA_tidy_lg 

str(df_background_DNA_LT_blood_PFOA_tidy_lg)

### PFAS和DNA均取对数了
df_background_DNA_LT_blood_PFOA_tidy %>% 
  dplyr::select(-c(1:18,29,30)) %>% 
  cbind(df_background_DNA_LT_blood_PFOA_tidy_lg) ->df_background_DNA_LT_blood_PFOA_tidy_final_lg 


str(df_background_DNA_LT_blood_PFOA_tidy_final_lg)






######################<<<<<<<<<<<<<<<<2.精浆PFAS>>>>##############################



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


####对PFAS取对数
lapply(df_background_DNA_LT_sperm_PFOA_tidy[,c(1:14)],log) ->df_background_DNA_LT_sperm_PFOA_tidy_lg 

### PFAS和DNA均取对数了
df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,25,26)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)







#########>>>>>>>>>>>>>>>>>>>>>>血浆和精浆PFAS之间的相关性<<<<<<<<<<<<<<<######## 

names(df_background_DNA_LT_blood_PFOA_tidy_final_lg)

df_background_DNA_LT_blood_PFOA_tidy_final_lg %>% 
  select(number,PFOA:FOSA_I) %>% 
  left_join(df_background_DNA_LT_sperm_PFOA_tidy_final_lg,by="number")  -> df_merger_blood_sperm_pfas
  
df_background_DNA_LT_sperm_PFOA_tidy_final_lg

df_merger_blood_sperm_pfas %>% 
  dplyr::filter(!PFOS.y=="NA") -> df_merger_blood_sperm_pfas_noNA

str(df_merger_blood_sperm_pfas_noNA)  ###血浆和精浆中PFAS配对的有796个


names(df_merger_blood_sperm_pfas_noNA)

####相关性分析    

df_merger_blood_sperm_pfas_noNA %>% 
  select(PFOA.x,PFOS.x,PFDA.x,L_PFHxS.x,PFUdA.x,PF3ONS_9CL.x,PFOA.y,PFOS.y,PFDA.y,L_PFHxS.y,PFUdA.y,PF3ONS_9CL.y
  ) -> df_all_match_pfas

str(df_all_match_pfas)
names(df_all_match_pfas)

cor(df_all_match_pfas)


# 将宽数据转换为长数据
long_df <- df_all_match_pfas %>%
  pivot_longer(
    cols = everything(),  # 选择所有列
    names_to = c("Variable", "Type"),  # 将列名拆分为 "Variable" 和 "Type"
    names_sep = "\\.",  # 以 "." 分隔列名
    values_to = "Value"  # 值列命名为 "Value"
  ) %>%
  pivot_wider(
    names_from = Type,  # 将 "Type" 列（x 和 y）作为新列
    values_from = Value  # 值从 "Value" 列中提取
  )

unlist(long_df)

long_df$Variable

# List of variables
variables <- long_df$Variable

# Repeat each variable 796 times
repeated_variables <- rep(variables, each = 796) %>% as_tibble() %>% 
                      rename(class="value")

long_df$x %>% unlist() %>% 
  as_tibble() %>% 
  rename(x="value")-> x


long_df$y %>% unlist() %>% 
  as_tibble() %>% 
  rename(y="value")-> y

cbind(repeated_variables,x,y) -> df_match_pfas_final





RColorBrewer::brewer.pal(n = 12, name = "Set3")

omics_color <-
  c(
    "PFOS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    "PFOA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    "PFDA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    "PFUdA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    "PFHxS" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    "6:2 Cl-PFESA" = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6]
  )


temp_data %>%
  # group_by(class) %>%
  dplyr::summarise(cor = cor(PF3ONS_9CL, PC, method = "spearman"))



cor.test(temp_data$PF3ONS_9CL,temp_data$PC,   method = "spearman")



####绝对值回归系数与p值
library(plyr)

library(dplyr)
library(tidyr)
library(purrr)
library(broom)

df_match_pfas_final %>%
  group_by(class) %>%  # Group by 'class'
  nest()  -> a  # Nest the data into a list column

a$data


####相关系数
a %>%
  mutate(correlation = map_dbl(data, ~ cor(.$x, .$y, method = "spearman")))


names(df_match_pfas_final)

df_match_pfas_final$class


####PFOA
df_match_pfas_final %>% 
  filter(class=="PFOA") -> df_match_pfas_final_pfoa



df_match_pfas_final_pfoa %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(df_match_pfas_final_pfoa$x, df_match_pfas_final_pfoa$y, method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()





####PFOS
df_match_pfas_final %>% 
  filter(class=="PFOS") -> df_match_pfas_final_PFOS



df_match_pfas_final_PFOS %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(df_match_pfas_final_PFOS$x, df_match_pfas_final_PFOS$y, method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()


# [1] "PFOA"       "PFOS"       "PFDA"       "L_PFHxS"    "PFUdA"      "PF3ONS_9CL"




####PFDA
df_match_pfas_final %>% 
  filter(class=="PFDA") -> df_match_pfas_final_PFDA



df_match_pfas_final_PFDA %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(df_match_pfas_final_PFDA$x, df_match_pfas_final_PFDA$y, method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()




####L_PFHxS
df_match_pfas_final %>% 
  filter(class=="L_PFHxS") -> df_match_pfas_final_L_PFHxS



df_match_pfas_final_L_PFHxS %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(df_match_pfas_final_L_PFHxS$x, df_match_pfas_final_L_PFHxS$y, method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()





####PFUdA
df_match_pfas_final %>% 
  filter(class=="PFUdA") -> df_match_pfas_final_PFUdA



df_match_pfas_final_PFUdA %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(df_match_pfas_final_PFUdA$x, df_match_pfas_final_PFUdA$y, method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()




####PF3ONS_9CL
df_match_pfas_final %>% 
  filter(class=="PF3ONS_9CL") -> df_match_pfas_final_PF3ONS_9CL



df_match_pfas_final_PF3ONS_9CL %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(df_match_pfas_final_PF3ONS_9CL$x, df_match_pfas_final_PF3ONS_9CL$y, method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()




df_match_pfas_final %>%
  mutate(PFAS=case_when(
    class=="PFOA"~"PFOA",
    class=="PFOS"~"PFOS",
    class=="PFDA"~"PFDA",
    class=="PFUdA"~"PFUdA",
    class=="L_PFHxS"~"PFHxS",
    class=="PF3ONS_9CL"~"6:2 Cl-PFESA"
  )) %>% 
  ggplot(aes(x,y)) +
  geom_point(aes(color = PFAS), show.legend = FALSE) +
  # theme_base +
  geom_smooth(aes(color = PFAS), method = "lm",
              show.legend = TRUE) +
  # facet_wrap(facets = vars(PFAS),
  #            scales = "free",
  #            nrow = 2) +
  scale_color_manual(values = omics_color) +
  theme_bw() +
  labs(x = "Plasma PFAS concentrations (ng/mL)", y = "Seminal PFAS concentraions (ng/mL)")+
  theme(legend.position = "none",
       legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
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
  facet_wrap(facets = vars(PFAS),
             scales = "free",
             nrow = 2) -> PFAS_COR


PFAS_COR

###导出保存为pdf
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/temp_data_blood_plasma_pfas_cor.pdf",width =12, height =8, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
PFAS_COR
# par(opar)
dev.off()




# 
# cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/pfas_meta_PF3ONS_9CL_merger_line.pdf",width = 12, height = 9, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# p8   # 使用大写罗马表示
# par(opar)
# dev.off()
# 
# 

