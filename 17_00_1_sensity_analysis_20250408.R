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
  dplyr::count(BMI_class_2)


# df_background_DNA_LT$marriage <- factor(df_background_DNA_LT$marriage,levels = c(1,2,4),labels = c("Unmarried","Married","Divorced"))
# 
# df_background_DNA_LT$marriage_class_2 <- factor(df_background_DNA_LT$marriage_class_2,levels = c(1,2),labels = c("Unmarried","Married"))
# 
# df_background_DNA_LT$age_class <- factor(df_background_DNA_LT$age_class,levels = c(1,2,3),labels = c("<25","25-30","30"))
# 
# df_background_DNA_LT$age_class_2 <- factor(df_background_DNA_LT$age_class_2,levels = c(0,1),labels = c("<28","≥28"))
# 
# df_background_DNA_LT$BMI_class <- factor(df_background_DNA_LT$BMI_class,levels = c(1,2,3),labels = c("18.5-24","<18.5","≥24"))
# 
# df_background_DNA_LT$BMI_class_2 <- factor(df_background_DNA_LT$BMI_class_2,
#                                            levels = c(1,2),labels = c("<24","≥24"))
# 
# df_background_DNA_LT$education <- factor(df_background_DNA_LT$education,levels = c(1,2,3),labels = c("High school","College","Undergraduate and above"))
# 
# df_background_DNA_LT$childn_class <- factor(df_background_DNA_LT$childn_class,levels = c(0,1),labels = c("No","Yes"))
# 
# df_background_DNA_LT$abstinence_class <- factor(df_background_DNA_LT$abstinence_class,levels = c(1,2,3),labels = c("2.1-7","≤2","≥7"))
# 
# df_background_DNA_LT$abstinence_class_2 <- factor(df_background_DNA_LT$abstinence_class_2,levels = c(1,2),labels = c("<7","≥7"))
# 
# df_background_DNA_LT$income <- factor(df_background_DNA_LT$income,levels = c(1,2,3),labels = c("<4000","4000-8000",">8000"))
# 
# df_background_DNA_LT$smk0 <- factor(df_background_DNA_LT$smk0,levels = c(1,2,3),labels = c("Current","Former","Never"))
# 
# df_background_DNA_LT$smk_class_2 <- factor(df_background_DNA_LT$smk_class_2,levels = c(0,1),labels = c("Never","Former/Current"))
# 
# df_background_DNA_LT$drk <- factor(df_background_DNA_LT$drk,levels = c(1,2,3,4),labels = c("Current","Former","Occasional","Never"))
# 
# df_background_DNA_LT$drk_class_2 <- factor(df_background_DNA_LT$drk_class_2,levels = c(1,2),labels = c("Never","Former/Current"))
# 
# df_background_DNA_LT$season<- factor(df_background_DNA_LT$season,levels = c(0,1,2,3),labels = c("Spring","Summer","Autumn","Winter"))
# 
# 
# df_background_DNA_LT$waterv_class <- factor(df_background_DNA_LT$waterv_class,levels = c(0,1,2),labels = c("0L","1-1000",">1000"))
# 
# df_background_DNA_LT$showinter_class <- factor(df_background_DNA_LT$showinter_class,levels = c(0,1,2),labels = c("<12","12.1-24",">24"))

df_background_DNA_LT$marriage_class_2 <- factor(df_background_DNA_LT$marriage_class_2)
df_background_DNA_LT$age_class_2 <- factor(df_background_DNA_LT$age_class_2)

df_background_DNA_LT$BMI_class_2 <- factor(df_background_DNA_LT$BMI_class_2)

df_background_DNA_LT$education <- factor(df_background_DNA_LT$education)

df_background_DNA_LT$abstinence_class_2 <- factor(df_background_DNA_LT$abstinence_class_2)

df_background_DNA_LT$smk_class_2 <- factor(df_background_DNA_LT$smk_class_2)

df_background_DNA_LT$drk_class_2 <- factor(df_background_DNA_LT$drk_class_2) 

df_background_DNA_LT$season<- factor(df_background_DNA_LT$season)

df_background_DNA_LT$income <- factor(df_background_DNA_LT$income)



str(df_background_DNA_LT)


df_background_DNA_LT %>% 
  dplyr::count(marriage_class_2)


df_background_DNA_LT %>% 
  dplyr::count(BMI_class_2)


df_background_DNA_LT %>% 
  dplyr::count(age_class_2)

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






###========================读取精浆PFOA和血浆PFOA数据===========================


###-----------------------------读入精浆PFOA的信息-----------------------------


df_sperm_PFOA <- read_excel("raw_data/guangming/血浆-精浆PFAS数据-final-20240524.xlsx",sheet = "精浆PFAS")

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

############################混合暴露模型########################################


####=============================端粒==================================


####--------------------------qgcomp模型-----------------------------------
####
library(qgcomp)
str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)


Xnm<-c("PFOS","PFOA","PFDA","PFUdA","L_PFHxS","PF3ONS_9CL")


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_sperm_PFOA_tidy_final_lg, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3

# devtools::install_github("alexpkeil1/qgcompint",force = TRUE)

library(qgcompint)

# qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_sperm_PFOA_tidy_final_lg, family=gaussian())
# 
# qgcomp.noboot(f=STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_sperm_PFOA_tidy_final_lg, family=gaussian()) 
# 
# 
###-----------------------------年龄分层---------------------------------------
####
####
### PFAS和DNA均取对数了 
df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,25,26)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  dplyr::count(age_class_2)

###将年龄转换为连续变量做交互效应


df_background_DNA_LT_sperm_PFOA_tidy_final_lg$age_class_2<- as.numeric(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$age_class_2)

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  dplyr::count(age_class_2)

qgcomp.emm.noboot(f=STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="age_class_2", 
                  expnms = c("PFOS","PFOA","PFDA","PFUdA","L_PFHxS","PF3ONS_9CL"), 
                  data=df_background_DNA_LT_sperm_PFOA_tidy_final_lg, 
                  q=4,
                  family=gaussian())




###将年龄<28岁的取出来  

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(age_class_2==1) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_age_1

qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_sperm_PFOA_tidy_final_lg_age_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3




###将年龄≥28岁的取出来 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(age_class_2==2) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_age_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_sperm_PFOA_tidy_final_lg_age_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3



###-----------------------------BMI分层---------------------------------------
####
####
### PFAS和DNA均取对数了 
df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,25,26)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  dplyr::count(BMI_class_2)

###将年龄转换为连续变量做交互效应


df_background_DNA_LT_sperm_PFOA_tidy_final_lg$BMI_class_2<- as.numeric(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$BMI_class_2)

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  dplyr::count(BMI_class_2)

qgcomp.emm.noboot(f=STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="BMI_class_2", 
                  expnms = c("PFOS","PFOA","PFDA","PFUdA","L_PFHxS","PF3ONS_9CL"), 
                  data=df_background_DNA_LT_sperm_PFOA_tidy_final_lg, 
                  q=4,
                  family=gaussian())





###将BMI<24的取出来  

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(BMI_class_2==1) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_BMI_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_sperm_PFOA_tidy_final_lg_BMI_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3




###将BMI≥24的取出来 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(BMI_class_2==2) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_BMI_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_sperm_PFOA_tidy_final_lg_BMI_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3






###-----------------------------吸烟分层---------------------------------------
####
####
### PFAS和DNA均取对数了 

df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,25,26)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  dplyr::count(smk_class_2)



###将年龄转换为连续变量做交互效应

df_background_DNA_LT_sperm_PFOA_tidy_final_lg$smk_class_2<- as.numeric(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$smk_class_2)

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  dplyr::count(smk_class_2)




qgcomp.emm.noboot(f=STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="smk_class_2", 
                  expnms = c("PFOS","PFOA","PFDA","PFUdA","L_PFHxS","PF3ONS_9CL"), 
                  data=df_background_DNA_LT_sperm_PFOA_tidy_final_lg, 
                  q=4,
                  family=gaussian())






###将never的取出来  

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)


df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(smk_class_2==1) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_smk_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df_background_DNA_LT_sperm_PFOA_tidy_final_lg_smk_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3





###将former/current的取出来 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(smk_class_2==2) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_smk_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_sperm_PFOA_tidy_final_lg_smk_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3





###-----------------------------饮酒分层---------------------------------------
####
####
### PFAS和DNA均取对数了 

df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,25,26)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  dplyr::count(drk_class_2)



###将年龄转换为连续变量做交互效应

df_background_DNA_LT_sperm_PFOA_tidy_final_lg$drk_class_2<- as.numeric(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$drk_class_2)

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)



qgcomp.emm.noboot(f=STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="drk_class_2", 
                  expnms = c("PFOS","PFOA","PFDA","PFUdA","L_PFHxS","PF3ONS_9CL"), 
                  data=df_background_DNA_LT_sperm_PFOA_tidy_final_lg, 
                  q=4,
                  family=gaussian())






###将never的取出来  

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)


df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(drk_class_2==1) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_drk_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+season+abstinence_class_2,expnms=Xnm, data = df_background_DNA_LT_sperm_PFOA_tidy_final_lg_drk_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3




###将former/current的取出来 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(drk_class_2==2) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_smk_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_sperm_PFOA_tidy_final_lg_smk_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3








###-----------------------------禁欲时间分层---------------------------------------
####
####
### PFAS和DNA均取对数了 

df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,25,26)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  dplyr::count(abstinence_class_2)



###将年龄转换为连续变量做交互效应

df_background_DNA_LT_sperm_PFOA_tidy_final_lg$abstinence_class_2<- as.numeric(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$abstinence_class_2)

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)
df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  dplyr::count(abstinence_class_2)


qgcomp.emm.noboot(f=STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="abstinence_class_2", 
                  expnms = c("PFOS","PFOA","PFDA","PFUdA","L_PFHxS","PF3ONS_9CL"), 
                  data=df_background_DNA_LT_sperm_PFOA_tidy_final_lg, 
                  q=4,
                  family=gaussian())






###将<7的取出来  

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)


df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(abstinence_class_2==1) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_abstinence_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season,expnms=Xnm, data = df_background_DNA_LT_sperm_PFOA_tidy_final_lg_abstinence_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3




###将≥7的取出来 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  filter(abstinence_class_2==2) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_abstinence_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season,expnms=Xnm, data = df_background_DNA_LT_sperm_PFOA_tidy_final_lg_abstinence_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3








###---------------------读入睡眠质量的数据----------------------------------- 


df_sleep <- read_excel("raw_data/睡眠_CE有的人取均数_842人_chen.xls",sheet = "Sheet1")

str(df_sleep)

###将原始样本与睡眠数据进行匹配
df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,25,26)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)


df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  left_join(df_sleep,by="number") -> df_background_DNA_LT_sperm_PFOA_tidy_sleep

str(df_background_DNA_LT_sperm_PFOA_tidy_sleep)

###删除睡眠缺失的样本
df_background_DNA_LT_sperm_PFOA_tidy_sleep %>% 
  dplyr::filter(!slp_score_mean=="NA") -> df_background_DNA_LT_sperm_PFOA_tidy_sleep_noNA


###最终只剩下完整含有睡眠的人群
length(unique(df_background_DNA_LT_sperm_PFOA_tidy_sleep_noNA$number))   ###仅剩833个人群

str(df_background_DNA_LT_sperm_PFOA_tidy_sleep_noNA)



###提取出分层变量
df_background_DNA_LT_sperm_PFOA_tidy_sleep_noNA %>% 
  dplyr::select(number:STL,slp_score_mean) -> df_background_DNA_LT_sperm_PFOA_tidy_sleep_aj


df_background_DNA_LT_sperm_PFOA_tidy_sleep_aj%>% 
  mutate(slp_class=case_when(
    slp_score_mean>5~1,     ####  1代表睡眠质量差
    slp_score_mean<=5~0    ####  0代表睡眠好
  )) -> df_sperm_DNA_Sleep_aj

str(df_sperm_DNA_Sleep_aj)

df_sperm_DNA_Sleep_aj%>% 
  dplyr::count(slp_class)



df_sperm_DNA_Sleep_aj$slp_class <- factor(df_sperm_DNA_Sleep_aj$slp_class)



df_sperm_DNA_Sleep_aj %>% 
  count(slp_class)


str(df_sperm_DNA_Sleep_aj)


###将睡眠转换为连续变量做交互效应

df_sperm_DNA_Sleep_aj$slp_class<- as.numeric(df_sperm_DNA_Sleep_aj$slp_class)

str(df_sperm_DNA_Sleep_aj)


df_sperm_DNA_Sleep_aj %>% 
  count(slp_class)    ###2代表质量差，1代表质量好


qgcomp.emm.noboot(f=STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2+slp_class,
                  emmvar="slp_class", 
                  expnms = c("PFOS","PFOA","PFDA","PFUdA","L_PFHxS","PF3ONS_9CL"), 
                  data=df_sperm_DNA_Sleep_aj, 
                  q=4,
                  family=gaussian())




###将睡眠质量好的取出来  

str(df_sperm_DNA_Sleep_aj)


df_sperm_DNA_Sleep_aj %>% 
  dplyr::filter(slp_class==1) -> df_sperm_DNA_Sleep_aj_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df_sperm_DNA_Sleep_aj_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3



###将睡眠质量差的取出来  

str(df_sperm_DNA_Sleep_aj)


df_sperm_DNA_Sleep_aj %>% 
  dplyr::filter(slp_class==2) -> df_sperm_DNA_Sleep_aj_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df_sperm_DNA_Sleep_aj_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3



###=========================读取体力活动数据======================= 
###
###

# str(df_background_DNA_LT_all_Hyper_tidy_lg)


df_physical_activity <- read_excel("raw_data/体力活动数据.xls",sheet = "Sheet")

str(df_physical_activity)

# df_physical_activity[duplicated(df_physical_activity$number), ]


df_background_DNA_LT_sperm_PFOA_tidy %>% 
  dplyr::select(-c(1:14,25,26)) %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)



df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  left_join(df_physical_activity,by="number") ->df_background_DNA_LT_physical 

str(df_background_DNA_LT_physical)

summary(df_background_DNA_LT_physical)



###删除身体活动缺失的样本
df_background_DNA_LT_physical  %>% 
  dplyr::filter(!mod_vig=="NA") -> df6_1


###最终只剩下完整含有身体活动的人群
length(unique(df6_1$number))   ###仅剩833个人群



df6_1 %>% 
  mutate(
    mod_vig_class=case_when(
      mod_vig<210~1,  ###
      mod_vig>=210~2
    )
  ) -> df6_2



str(df6_2)
# 
# df6_2 %>% 
#   dplyr::select(-1) -> df7  ###得到最终用于分层数据

df6_2 -> df7

df7$mod_vig_class <- factor(df7$mod_vig_class)


df7 %>% 
  count(mod_vig_class)    ###1代表<210，2代表>210



###将身体活动转换为连续变量做交互效应

df7$mod_vig_class<- as.numeric(df7$mod_vig_class)



df7 %>% 
  count(mod_vig_class)     ###1代表<210，2代表>210

qgcomp.emm.noboot(f=STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2+mod_vig_class,
                  emmvar="mod_vig_class", 
                  expnms = c("PFOS","PFOA","PFDA","PFUdA","L_PFHxS","PF3ONS_9CL"), 
                  data=df7, 
                  q=4,
                  family=gaussian())




###将身体活动<210的取出来  

str(df7)


df7 %>% 
  dplyr::filter(mod_vig_class==1) -> df7_1




qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df7_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3



###将身体活动>210的取出来  

str(df7)


df7 %>% 
  dplyr::filter(mod_vig_class==2) -> df7_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOS+PFOA+PFDA+PFUdA+L_PFHxS+PF3ONS_9CL+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df7_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3


























###------------------读入血浆PFOA的信息-----------------------

df_blood_PFOA <- read_excel("raw_data/guangming/血浆-精浆PFAS数据-final-20240524.xlsx",sheet = "血浆PFAS")

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


###----------------单因素线性回归分析--------------------------------




####------------将LOD大于80%的做多元混合线性模型----------------------

####---------------四分为间距--------------------------

str(df_background_DNA_LT_blood_PFOA_tidy_final_lg)

# PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS

####利用循环函数，批量分段重组

lapply(df_background_DNA_LT_blood_PFOA_tidy_final_lg[,c("PFOA","PFOS","PFNA","PFDA","PFUdA","L_PFHxS","diPAP_6_2","PF3OUdS_11CL","PF3ONS_9CL","PFDoA","L_PFHpS")], function(x){ 
  cut(x,breaks =quantile(x),include.lowest=T,
      labels = c(1,2,3,4)) %>% 
    as_tibble() 
}) %>% 
  bind_cols()  %>%    ####将所有列数据合并在一个数据框中
  set_names(paste(data.frame(colnames(df_background_DNA_LT_blood_PFOA_tidy_final_lg))[c(11,12,13,14,15,16,17,18,19,20,21),],"_","class",sep = "")) -> df_background_DNA_LT_blood_PFOA_tidy_final_lg_class  ####对生成的数据框根据列变量进行重新命名
str(df_background_DNA_LT_blood_PFOA_tidy_final_lg_class)


####将原始数据与新生成的分类数据进行合并 

df_background_DNA_LT_blood_PFOA_tidy_final_lg %>% 
  cbind(df_background_DNA_LT_blood_PFOA_tidy_final_lg_class) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)


####------------将LOD大于80%的做多元混合线性模型----------------------

###精浆检出率大于80%的有PFOS+PFOA+PFUdA+L_PFHxS+PF3ONS_9CL
##11,12,14,






####=============================端粒==================================


####--------------------------qgcomp模型-----------------------------------
####
library(qgcomp)
str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

# names(df_background_DNA_LT_PFOA_tidy)[1:13]

Xnm<-c("PFOA","PFOS","PFNA","PFDA","PFUdA","L_PFHxS","diPAP_6_2","PF3OUdS_11CL","PF3ONS_9CL","PFDoA","L_PFHpS")



A6<- qgcomp.glm.noboot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2, expnms=Xnm, data =df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class, q=4,family=gaussian())
A6
head(A6$qx)
plot(A6)
# dev.off()

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

summary(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2, expnms=Xnm, data =df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class,family=gaussian(), q=4, B=100,seed=125)

###混合物中系数
summary(qcboot.fit3$fit)$coefficients

qcboot.fit3



# 
###-----------------------------年龄分层---------------------------------------
####
####
### PFAS和DNA均取对数了 
df_background_DNA_LT_blood_PFOA_tidy_final_lg %>% 
  cbind(df_background_DNA_LT_blood_PFOA_tidy_final_lg_class) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class%>% 
  dplyr::count(age_class_2)

###将年龄转换为连续变量做交互效应


df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$age_class_2<- as.numeric(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$age_class_2)

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  dplyr::count(age_class_2)

qgcomp.emm.noboot(f=STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="age_class_2", 
                  expnms = c("PFOA","PFOS","PFNA","PFDA","PFUdA","L_PFHxS","diPAP_6_2","PF3OUdS_11CL","PF3ONS_9CL","PFDoA","L_PFHpS"), 
                  data=df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class, 
                  q=4,
                  family=gaussian())




###将年龄<28岁的取出来  

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  filter(age_class_2==1) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_age_1

qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_age_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3




###将年龄≥28岁的取出来 

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  filter(age_class_2==2) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_age_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_age_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3



###-----------------------------BMI分层---------------------------------------
####
####
### PFAS和DNA均取对数了 
df_background_DNA_LT_blood_PFOA_tidy_final_lg %>% 
  cbind(df_background_DNA_LT_blood_PFOA_tidy_final_lg_class) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  dplyr::count(BMI_class_2)

###将年龄转换为连续变量做交互效应


df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$BMI_class_2<- as.numeric(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$BMI_class_2)

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  dplyr::count(BMI_class_2)


qgcomp.emm.noboot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="BMI_class_2", 
                  expnms = c("PFOA","PFOS","PFNA","PFDA","PFUdA","L_PFHxS","diPAP_6_2","PF3OUdS_11CL","PF3ONS_9CL","PFDoA","L_PFHpS"), 
                  data=df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class, 
                  q=4,
                  family=gaussian())





###将BMI<24的取出来  

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  filter(BMI_class_2==1) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_BMI_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_BMI_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3




###将BMI≥24的取出来 

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  filter(BMI_class_2==2) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_BMI_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_BMI_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3






###-----------------------------吸烟分层---------------------------------------
####
####
### PFAS和DNA均取对数了 

df_background_DNA_LT_blood_PFOA_tidy_final_lg %>% 
  cbind(df_background_DNA_LT_blood_PFOA_tidy_final_lg_class) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)


df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  dplyr::count(smk_class_2)



###将年龄转换为连续变量做交互效应

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$smk_class_2<- as.numeric(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$smk_class_2)

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  dplyr::count(smk_class_2)




qgcomp.emm.noboot(f=STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="smk_class_2", 
                  expnms = c("PFOA","PFOS","PFNA","PFDA","PFUdA","L_PFHxS","diPAP_6_2","PF3OUdS_11CL","PF3ONS_9CL","PFDoA","L_PFHpS"), 
                  data=df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class, 
                  q=4,
                  family=gaussian())






###将never的取出来  

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)


df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class%>% 
  filter(smk_class_2==1) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_smk_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_smk_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3





###将former/current的取出来 

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  filter(smk_class_2==2) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_smk_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+drk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_smk_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3





###-----------------------------饮酒分层---------------------------------------
####
####
### PFAS和DNA均取对数了 


df_background_DNA_LT_blood_PFOA_tidy_final_lg %>% 
  cbind(df_background_DNA_LT_blood_PFOA_tidy_final_lg_class) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  dplyr::count(drk_class_2)



###将年龄转换为连续变量做交互效应

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$drk_class_2<- as.numeric(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$drk_class_2)

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)



qgcomp.emm.noboot(f=STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="drk_class_2", 
                  expnms = c("PFOA","PFOS","PFNA","PFDA","PFUdA","L_PFHxS","diPAP_6_2","PF3OUdS_11CL","PF3ONS_9CL","PFDoA","L_PFHpS"), 
                  data=df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class, 
                  q=4,
                  family=gaussian())




###将never的取出来  

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)


df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  filter(drk_class_2==1) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_drk_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+season+abstinence_class_2,expnms=Xnm, data = df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_drk_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3




###将former/current的取出来 

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class%>% 
  filter(drk_class_2==2) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_drk_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+season+abstinence_class_2,expnms=Xnm, data =df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_drk_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3








###-----------------------------禁欲时间分层---------------------------------------
####
####
### PFAS和DNA均取对数了 


df_background_DNA_LT_blood_PFOA_tidy_final_lg %>% 
  cbind(df_background_DNA_LT_blood_PFOA_tidy_final_lg_class) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)



df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  dplyr::count(abstinence_class_2)



###将年龄转换为连续变量做交互效应

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$abstinence_class_2<- as.numeric(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class$abstinence_class_2)

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  dplyr::count(abstinence_class_2)


qgcomp.emm.noboot(f=STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,
                  emmvar="abstinence_class_2", 
                  expnms = c("PFOA","PFOS","PFNA","PFDA","PFUdA","L_PFHxS","diPAP_6_2","PF3OUdS_11CL","PF3ONS_9CL","PFDoA","L_PFHpS"), 
                  data=df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class, 
                  q=4,
                  family=gaussian())






###将<7的取出来  

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)


df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  filter(abstinence_class_2==1) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_abstinence_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season,expnms=Xnm, data = df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_abstinence_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3




###将≥7的取出来 

df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  filter(abstinence_class_2==2) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_abstinence_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season,expnms=Xnm, data =df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class_abstinence_2, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3








###---------------------读入睡眠质量的数据----------------------------------- 


df_sleep <- read_excel("raw_data/睡眠_CE有的人取均数_842人_chen.xls",sheet = "Sheet1")

str(df_sleep)

###将原始样本与睡眠数据进行匹配
df_background_DNA_LT_blood_PFOA_tidy_final_lg %>% 
  cbind(df_background_DNA_LT_blood_PFOA_tidy_final_lg_class) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)



df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  left_join(df_sleep,by="number") -> df_background_DNA_LT_blood_PFOA_tidy_sleep

str(df_background_DNA_LT_blood_PFOA_tidy_sleep)

###删除睡眠缺失的样本
df_background_DNA_LT_blood_PFOA_tidy_sleep %>% 
  dplyr::filter(!slp_score_mean=="NA") -> df_background_DNA_LT_blood_PFOA_tidy_sleep_noNA


###最终只剩下完整含有睡眠的人群
length(unique(df_background_DNA_LT_blood_PFOA_tidy_sleep_noNA$number))   ###仅剩833个人群

str(df_background_DNA_LT_blood_PFOA_tidy_sleep_noNA)



###提取出分层变量
df_background_DNA_LT_blood_PFOA_tidy_sleep_noNA %>% 
  dplyr::select(number:L_PFHpS_class,slp_score_mean) -> df_background_DNA_LT_blood_PFOA_tidy_sleep_aj


df_background_DNA_LT_blood_PFOA_tidy_sleep_aj%>% 
  mutate(slp_class=case_when(
    slp_score_mean>5~1,     ####  1代表睡眠质量差
    slp_score_mean<=5~0    ####  0代表睡眠好
  )) -> df_blood_DNA_Sleep_aj

str(df_blood_DNA_Sleep_aj)

df_blood_DNA_Sleep_aj%>% 
  dplyr::count(slp_class)



df_blood_DNA_Sleep_aj$slp_class <- factor(df_blood_DNA_Sleep_aj$slp_class)





df_blood_DNA_Sleep_aj%>% 
  count(slp_class)


str(df_blood_DNA_Sleep_aj)


###将睡眠转换为连续变量做交互效应

df_blood_DNA_Sleep_aj$slp_class<- as.numeric(df_blood_DNA_Sleep_aj$slp_class)

str(df_blood_DNA_Sleep_aj)


df_blood_DNA_Sleep_aj %>% 
  count(slp_class)    ###2代表质量差，1代表质量好



qgcomp.emm.noboot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2+slp_class,
                  emmvar="slp_class", 
                  expnms = c("PFOA","PFOS","PFNA","PFDA","PFUdA","L_PFHxS","diPAP_6_2","PF3OUdS_11CL","PF3ONS_9CL","PFDoA","L_PFHpS"), 
                  data=df_blood_DNA_Sleep_aj, 
                  q=4,
                  family=gaussian())




###将睡眠质量好的取出来  

str(df_blood_DNA_Sleep_aj)


df_blood_DNA_Sleep_aj %>% 
  dplyr::filter(slp_class==1) -> df_blood_DNA_Sleep_aj_1


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df_blood_DNA_Sleep_aj_1, family=gaussian(), q=4, B=100,seed=125)

qcboot.fit3



###将睡眠质量差的取出来  

str(df_blood_DNA_Sleep_aj)


df_blood_DNA_Sleep_aj %>% 
  dplyr::filter(slp_class==2) -> df_blood_DNA_Sleep_aj_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df_blood_DNA_Sleep_aj_2, family=gaussian(), q=4, B=10000,seed=125)

qcboot.fit3





###=========================读取体力活动数据======================= 
###
###

# str(df_background_DNA_LT_all_Hyper_tidy_lg)


df_physical_activity <- read_excel("raw_data/体力活动数据.xls",sheet = "Sheet")

str(df_physical_activity)

# df_physical_activity[duplicated(df_physical_activity$number), ]

df_background_DNA_LT_blood_PFOA_tidy_final_lg %>% 
  cbind(df_background_DNA_LT_blood_PFOA_tidy_final_lg_class) -> df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class

str(df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class)




df_background_DNA_LT_blood_PFOA_tidy_final_LOD_Class %>% 
  left_join(df_physical_activity,by="number") ->df_background_DNA_LT_physical 

str(df_background_DNA_LT_physical)

summary(df_background_DNA_LT_physical)



###删除身体活动缺失的样本
df_background_DNA_LT_physical  %>% 
  dplyr::filter(!mod_vig=="NA") -> df6_1


###最终只剩下完整含有身体活动的人群
length(unique(df6_1$number))   ###仅剩833个人群



df6_1 %>% 
  mutate(
    mod_vig_class=case_when(
      mod_vig<210~1,  ###
      mod_vig>=210~2
    )
  ) -> df6_2



str(df6_2)

# 
# df6_2 %>% 
#   dplyr::select(-1) -> df7  ###得到最终用于分层数据

df6_2 -> df7

df7$mod_vig_class <- factor(df7$mod_vig_class)


df7 %>% 
  count(mod_vig_class)    ###1代表<210，2代表>210



###将身体活动转换为连续变量做交互效应

df7$mod_vig_class<- as.numeric(df7$mod_vig_class)

str(df7)


df7 %>% 
  count(mod_vig_class)     ###1代表<210，2代表>210

qgcomp.emm.noboot(f=STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2+mod_vig_class,
                  emmvar="mod_vig_class", 
                  expnms = c("PFOA","PFOS","PFNA","PFDA","PFUdA","L_PFHxS","diPAP_6_2","PF3OUdS_11CL","PF3ONS_9CL","PFDoA","L_PFHpS"), 
                  data=df7, 
                  q=4,
                  family=gaussian())





###将身体活动<210的取出来  

str(df7)


df7 %>% 
  dplyr::filter(mod_vig_class==1) -> df7_1




qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df7_1, family=gaussian(), q=4, B=10000,seed=125)

qcboot.fit3



###将身体活动>210的取出来  

str(df7)


df7 %>% 
  dplyr::filter(mod_vig_class==2) -> df7_2


qcboot.fit3 <- qgcomp.glm.boot(STL ~ PFOA+PFOS+PFNA+PFDA+PFUdA+L_PFHxS+diPAP_6_2+PF3OUdS_11CL+PF3ONS_9CL+PFDoA+L_PFHpS+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,expnms=Xnm, data = df7_2, family=gaussian(), q=4, B=10000,seed=125)

qcboot.fit3







