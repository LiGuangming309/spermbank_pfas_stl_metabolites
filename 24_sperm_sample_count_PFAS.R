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




















###============2.读入重复测量的精液质量参数=============================

df2 <- read_excel("raw_data/semen_quality_parameters.xlsx",sheet = "匹配表")

df2 %>%  
  dplyr::select(1:29) -> df2_1  ##每一次精液质量参数与最后一次精液参数匹配后的数据
str(df2_1)
length(unique(df2_1$number...1 ))   #1487份人群精液样本，重复测量6608次


# df2_1[complete.cases(df2_1),] -> df2_1_1  #删除缺失值
# write.xlsx(df2_1,"result/20240628_Hyperten_sperm_2th/20240628_result/semen_quality_parameters_mtDNAcn_match.xlsx") #导出筛选的变量  

df2_1[!is.na(df2_1$mtDNAcn...17),] -> df3
length(unique(df3$number...1))   ###共有1164人 

1487-1164

# write.xlsx(df3,"result/20240628_Hyperten_sperm_2th/20240628_result/semen_quality_parameters_mtDNAcn_20240509.xlsx") #导出筛选的变量

# 
# df2 %>%  
#   dplyr::select(42:59) -> df4  ##每一次精液质量参数与最后一次精液参数匹配后的数据
# 
# str(df4)
# 
# df4_1 <- na.omit(df4)
# 
# str(df4_1)


###或者本研究选择纳入到2018年4月30号
###去掉2018年5月1号之后到2018年7月15号数据,即保留2018年4月30号之前的数据

df2_1 %>% 
  dplyr::rename(number=number...1) -> df2_2
str(df2_2)


###在这一步生成日期以及带有mtDNA的数据是为了进一步比较和排除日期数据是否合理

df2_2 %>% 
  mutate(datemq=lubridate::ymd(datesq)) %>% 
  arrange(number, datemq) %>% 
  group_by(number) %>% 
  mutate(month = month(datemq)) ->df2_3_month   



df2_3_month %>% 
  mutate(season_class=case_when(
    month==3|month==4|month==5~1,  ###春季
    month==6|month==7|month==8~2,  ###夏季
    month==9|month==10|month==11~3,  ###秋季
    month==12|month==1|month==2~4  ###秋季
  )) -> df2_3_month_class

# 
# df2_3_month_class %>% 
#   dplyr::count(season_class)

table(df2_3_month_class$season_class)

# write.xlsx(df2_3_month_class,"result/20240628_Hyperten_sperm_2th/20240628_result/df_background_repeatsemen_month.xlsx")


str(df2_3_month_class)

###最大和最小日期 

min(df2_3_month_class$datemq)
max(df2_3_month_class$datemq)


# 设定日期区间
start_date <- as.Date("2018-05-01")
end_date <- as.Date("2018-07-15")

# ###---------不排除2018年4月30号之后的随访数据----------
# 
# # 筛选出日期在指定区间内的行，即2018年4月30号之前的数据
# 
# df2_3_month %>% 
#   dplyr::filter(datemq<ymd(20180501)) -> df2_4
# 
# str(df2_4)
# 
# 
# length(unique(df2_4$number))   ###共有1487人 
# 
# 
# ##剔除2018年5月1日之后的数据，共有多少人群样本
# df2_3_month %>% 
#   dplyr::filter(datemq>=ymd(20180501)) -> df2_5
# 
# str(df2_5)
# 
# length(unique(df2_5$number))    ####总共有174人
# 


###将重复测量中精液质量参数异常结果剔除
df2_3_month_class %>%  
  dplyr::mutate(group_semen=case_when(
    motilityd...13==100|motilityd...13<0~1,   ###1代表精液质量数据有异常的数据，0代表精液质量数据无异常值
    TRUE~0
  )) -> df2_6

# write.xlsx(df2_6,"result/20240628_Hyperten_sperm_2th/20240628_result/semen_quality_parameters(将有问题的精液质量筛检生成分组).xlsx") ##导出筛选的变量  

str(df2_6)

length(unique(df2_6$number))    ####总共有1487人


df2_6 %>% 
  dplyr::filter(group_semen==0) -> df2_7  ###将死精或无精的样本删除

str(df2_7)

length(unique(df2_7$number))    ####总共有1485人

6608-6596  ###12份精液质量有问题

df2_6 %>% dplyr::filter(group_semen==1) -> df2_8   ###

length(unique(df2_8$number))   #11人份精液样本有问题，总共12个重复精液样本，但是这些人可能还有其他重复测量精液样本

str(df2_8)
# view(df2_1_3)
str(df2_8)

df2_7$density...8[which(df2_7$density...8 ==0)] <- 0.01

# write.xlsx(df2_7,"result/20240628_Hyperten_sperm_2th/20240628_result/semen_quality_parameters(将有问题的生成分组).xlsx") ##导出筛选的变量  

str(df2_7)

length(unique(df2_7$number))   #1487份人群精液样本，重复测量6608次




###=================数据集合并以及变量筛选================================

###==将一般资料的数据与精液质量参数进行匹配


df2_7%>% 
  dplyr::select(number:QJ...15,datemq,month,season_class) -> df2_9  ###挑选出重复测量精液质量的数据

str(df2_9)

length(unique(df2_9$number))

summary(df2_9)



# str(df1_7)
# length(unique(df1_7$number))


df2_9  -> df_background_repeatsemen

# colSums(is.na(df_background_repeatsemen)) 
str(df_background_repeatsemen)
# 


length(unique(df_background_repeatsemen$number))  #匹配得到1383个男性，重复捐精6443个样本

df_background_repeatsemen[is.na(df_background_repeatsemen$motilitya...10),] -> df_background_repeatsemen_NA 

str(df_background_repeatsemen_NA)

# write.xlsx(df_background_repeatsemen_NA,"result/20240628_Hyperten_sperm_2th/20240628_result/df_no_sperm_20240509.xlsx")




df_background_repeatsemen %>% 
  dplyr::filter(!motilitya...10=="NA") ->df_background_repeatsemen_noNA 

str(df_background_repeatsemen_noNA)   ###1164份

length(unique(df_background_repeatsemen_noNA$number))

colSums(is.na(df_background_repeatsemen_noNA)) 




# 
# df_background_repeatsemen_noNA %>% 
#   mutate(datemq=lubridate::ymd(datesq)) %>% 
#   arrange(number, datemq) %>% 
#   group_by(number) %>% 
#   mutate(month = month(datemq)) ->df_background_repeatsemen_month 
# 
# write.xlsx(df_background_repeatsemen_month,"result/20240318_Hyperten_sperm/df_background_repeatsemen_month.xlsx")
# 

df_background_repeatsemen_noNA %>% 
  mutate(datesq=lubridate::ymd(datesq)) %>% 
  arrange(number, datesq) %>% 
  group_by(number) %>% 
  mutate(firstday = min(datesq)) %>% 
  mutate(days = datesq - firstday) %>%  #计算每个ID下每次随访时间与基线时间间隔的天数
  ungroup()  -> df_background_repeatsemen_days   #得到含有时间天数的变量


# write.xlsx(df_background_repeatsemen_days,"result/20240628_Hyperten_sperm_2th/20240628_result/df_background_repeatsemen_days.xlsx")


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

# str(df_background_repeatsemen_final2)
# write.xlsx(df_background_repeatsemen_final2,"result/20240628_Hyperten_sperm_2th/20240628_result/df_background_repeatsemen_final2.xlsx")

#### ==============提取出重复测量的精液质量参数数据的基线期数据================

###提取基线数据  
df_background_repeatsemen_final2 %>% 
  dplyr::filter(interval==0) -> df_background_repeatsemen_base_lifestyle

str(df_background_repeatsemen_base_lifestyle)

colSums(is.na(df_background_repeatsemen_base_lifestyle))
table(df_background_repeatsemen_base_lifestyle$interval)




df_background_DNA_LT__tidy_lg %>% 
  left_join(df_background_repeatsemen_final2,by="number") -> df_backgroud_STL_sperm



colSums(is.na(df_backgroud_STL_sperm))






###生成基线精液质量参数
df_backgroud_STL_sperm  %>% 
  dplyr::mutate(sperm_count=volume...7*density...8,
                Total_motility=motilitya...10+motilityb...11+motilityc...12,
                Progressive_motility=motilitya...10+motilityb...11) -> df_background_DNA_LT_all_Hyper_tidy_2

###对精液质量参数取对数
df_background_DNA_LT_all_Hyper_tidy_2 %>% 
  dplyr::mutate(
    lg_volume=log(volume...7),
    lg_density=log(density...8),
    lg_sperm_count=log(sperm_count),
    lg_total_motility=log(Total_motility),
    lg_progressive_motility=log(Progressive_motility)) ->df_background_DNA_LT_all_Hyper_tidy_lg 


df_background_DNA_LT_all_Hyper_tidy_lg  -> df_background_DNA_LT_all_Hyper_tidy_2_final




names(df_background_DNA_LT_all_Hyper_tidy_2_final)









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



df_background_DNA_LT_all_Hyper_tidy_2_final %>% 
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


df_background_DNA_LT_sperm_PFOA_noNA$interval

table(df_background_DNA_LT_sperm_PFOA_noNA$interval)


df_background_DNA_LT_sperm_PFOA_noNA %>% 
  select(interval,number)



library(dplyr)

df_background_DNA_LT_sperm_PFOA_noNA %>%
  group_by(number) %>%
  summarise(interval_count = n())   # 统计每个个体的行数（即interval次数）




library(dplyr)


# 计算number中每个个体，在每个interval中出现的次数

df_background_DNA_LT_sperm_PFOA_noNA  %>%
  group_by(number, interval) %>%
  summarise(count = n(), .groups = "drop") -> df_background_DNA_LT_sperm_PFOA_noNA_2



df_background_DNA_LT_sperm_PFOA_noNA_2  


####每个样本对应的区间个数，由于每个区间只取第一次采集的精液样本进行等体积的精液样本混合，所以这里只取第一次的就行

df_background_DNA_LT_sperm_PFOA_noNA_2  %>%
  group_by(number) %>%
  summarise(interval_count = n()) -> df_background_DNA_LT_sperm_PFOA_noNA_3   # 统计每个个体的行数（即interval次数）

####每个个体平均3.03±1.37个精液样本

mean(df_background_DNA_LT_sperm_PFOA_noNA_3$interval_count)
sd(df_background_DNA_LT_sperm_PFOA_noNA_3$interval_count)

















# 
# ###将数据重新排序整理
# df_background_DNA_LT_sperm_PFOA_noNA %>% 
#   select(PFOS:PF3OUdS_11CL,mtDNAcn,STL,everything()) -> df_background_DNA_LT_sperm_PFOA_tidy
# 
# 
# str(df_background_DNA_LT_sperm_PFOA_tidy)
# 
# # 
# # ####对PFAS取对数
# lapply(df_background_DNA_LT_sperm_PFOA_tidy[,c(1:14)],log) ->df_background_DNA_LT_sperm_PFOA_tidy_lg 
# 
# ### PFAS和DNA均取对数了
# df_background_DNA_LT_sperm_PFOA_tidy %>% 
#   dplyr::select(-c(1:14,25,26)) %>% 
#   cbind(df_background_DNA_LT_sperm_PFOA_tidy_lg) ->df_background_DNA_LT_sperm_PFOA_tidy_final_lg 
# 
# str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)
# 
# 
# 
# 
# 
# 
# 








