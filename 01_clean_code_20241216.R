
pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel,scales,compareGroups)

# get_project_wd()

rm(list = ls())
getwd()

# library(tidymass)
###读取精子库一般基本信息数据
df1 <- read_excel("raw_data/background.xlsx",sheet = "Sheet")
str(df1)   #总共收集1487份样本

df1  %>% 
  dplyr::select(number,dateaq,age,education,marriage,childn,income,smk0,drk,height,weight,sbp,dbp,waterv,showinter)  -> df1_4  #筛选出基本信息的变量

str(df1_4)
summary(df1_4)

df1_4$education <- factor(df1_4$education,levels = c("专科","中专","大专","高中","本科","硕士","研究生","博士"),labels = c(1,1,1,1,2,3,3,3))

# df1$education <- factor(df1$education,levels = c("专科","中专","大专","高中","本科","硕士","研究生","博士"),labels = c(1,1,1,1,2,3,3,3))


df1_4$education
table(df1_4$education)
str(df1_4)
colSums(is.na(df1_4))   #查看缺失值

str(df1_4)

2/1487

# df1_4 %>% 
#   dplyr::filter(!sbp=="NA") -> df1_4_non_NA
# 
# colSums(is.na(df1_4_non_NA))   #查看缺失值
# str(df1_4_non_NA)  ###删除2个高血压缺失



# # df1_4$education <- as.numeric(df1_4$education)
# 
# ?mice
# 
# # write.xlsx(df1_4,"result/screen_background.xlsx") #导出筛选的变量  
# 
# ###=========对一般资料进行数据插补===================
# 
# library(mice)
# 
# md.pattern(df1_4_non_NA)   ##查看数据缺失情况 
# 
# #### 利用mice函数对数值变量插补法
# colSums(is.na(df1_4))   #查看缺失值
# 

df1_4 %>%
  dplyr::select(age,height,weight,waterv) -> df1_4_continuous
?mice()
str(df1_4_continuous)
mice(df1_4_continuous,m=5,maxit = 50,seed = 1000,method="rf") -> df1_4_continuous_imp
summary(df1_4_continuous_imp)

complete(df1_4_continuous_imp,action = 5) -> df1_4_continuous_imp_final

colSums(is.na(df1_4_continuous_imp_final))   #查看缺失值

####利用mice函数对因变量进行插补

df1_4 %>%
  dplyr::select(education,childn,income) -> df1_4_class

mutate_if(df1_4_class, is.numeric, as.factor) -> df1_4_class_factor  ####将数值变量转换为因变量


mice(df1_4_class_factor,m=5,maxit = 50,seed = 1000,method="polyreg") -> df1_4_class_factor_imp
summary(df1_4_class_factor_imp)

complete(df1_4_class_factor_imp,action = 5) -> df1_4_class_factor_imp_final

colSums(is.na(df1_4_class_factor_imp_final))   #查看缺失值


str(df1_4_class_factor_imp_final)


####将插补的后的变量与之前的数据进行合并

df1_4 %>%
  dplyr::select(-c(age,height,weight,waterv,childn,education,income)) %>%
  bind_cols(df1_4_continuous_imp_final,df1_4_class_factor_imp_final)  -> df1_5   #最终插补后的数据

str(df1_5)





df1_5$smk0 <- factor(df1_5$smk0)
df1_5 %>% 
  mutate(
    income=case_when(
      income==1 |income==2~1,
      income==3 |income==4~2,
      income==5 |income==6~3
    )
  )  -> df1_5





df1_5$childn <- factor(df1_5$childn)
df1_5$drk <- factor(df1_5$drk)

df1_5$income<- factor(df1_5$income)
# df1_5$occupation_exposure<- factor(df1_5$occupation_exposure)
# df1_5$tea<- factor(df1_5$tea)
str(df1_5)

df1_5 %>% 
  mutate(BMI=weight/(height/100)^2)  -> df1_5
df1_5$BMI



df1_5 %>% 
  mutate(
    BMI_class=case_when(
      BMI>=18.5&BMI<24 ~1,
      BMI<18.5 ~2,
      BMI>=24 ~3
    )
  ) -> df1_5

df1_5$BMI_class


df1_5 %>% 
  mutate(
    age_class=case_when(
      age>=25&age<=30 ~2,
      age<25 ~1,
      age>30 ~3
    )
  ) -> df1_5




df1_5 %>% 
  mutate(
    childn_class=case_when(
      childn==0 ~0,
      .default = 1
    )
  ) -> df1_5

?case_when
str(df1_5)


df1_5 %>% 
  mutate(
    waterv_class=case_when(
      waterv==0 ~0,
      waterv>=1& waterv<=1000~1,
      waterv>1000~2
    ),
    showinter_class=case_when(
      showinter<=12 ~0,
      showinter>=12.1 & showinter<=24~1,
      showinter>24~2
    )
  ) ->df1_5 


summary(df1_5)

str(df1_5)
df1_5$showinter_class


mean(df1_5$age)
sd(df1_5$age)

# write.xlsx(df1_5,"result/20240106_Hyperten_DNA_sperm/interference_background.xlsx") #导出插补后的数据


###=======读入基线有疾病史标记的102人=============================

###读取精子库102例有问题的数据
df_102_disease <- read_excel("raw_data/102人有疾病人群.xlsx",sheet = "Sheet1")
str(df_102_disease)   #总共收集1487份样本

df1_5 %>% 
  left_join(df_102_disease,by="number") -> df1_6   ###匹配后数据


str(df1_6)


# df1_6$情况 <- as.numeric(df1_6$情况)
# 
# df1_6 %>% 
#    dplyr::filter(disease==NA)
# 
# summary(df1_6)

###筛选出疾病中有缺失值的数据
df1_6[is.na(df1_6$disease),] ->df1_7 


str(df1_7)




# write.xlsx(df1_5,"result/20240106_Hyperten_DNA_sperm/interference_background.xlsx") #导出插补后的数据





###======读入重复测量的精液质量参数

df2 <- read_excel("raw_data/semen_quality_parameters.xlsx",sheet = "匹配表")


df2 %>%  
  dplyr::select(42:59)  -> df2_2  ##最后一次检测mtDNAcn的精液质量参数

str(df2_2)


df2_2[complete.cases(df2_2),] -> df2_3  #删除缺失值
str(df2_3)            #1164份含有mtDNA样本

df2_3

df2_3 %>% 
  dplyr::rename(number=number...42) -> df2_4 

###将基线与DNA数据合并
df1_7%>%
  left_join(df2_4,by="number")  -> df_background_DNA


# colSums(is.na(df_background_DNA)) 
str(df_background_DNA)
# 
length(unique(df_background_DNA$mtDNAcn...47)) #匹配得到1165个男性
# 


###得到1164个基线和精子DNA和端粒测量的数据
df_background_DNA %>% 
  dplyr::filter(!mtDNAcn...47=="NA") ->df_background_DNA_noNA 

str(df_background_DNA_noNA)   ###1164份

length(unique(df_background_DNA_noNA$mtDNAcn...47))

write.xlsx(df_background_DNA_noNA,"result/20240407_PFOA_sperm/df_background_DNA_noNA_final.xlsx") ##导出筛选的变量











##################精液质量参数##############################################



df2 %>%  
  dplyr::select(1:29) -> df2_1  ##每一次精液质量参数与最后一次精液参数匹配后的数据



###或者本研究选择纳入到2018年4月30号
###去掉2018年5月1号之后到2018年7月15号数据,即保留2018年4月30号之前的数据

df2_1 %>% 
  dplyr::rename(number=number...1) -> df2_1_2
str(df2_1_2)


###在这一步生成日期以及带有mtDNA的数据是为了进一步比较和排除日期数据是否合理

df2_1_2 %>% 
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
  )) -> df2_1_4

# write.xlsx(df2_6,"result/20240628_Hyperten_sperm_2th/20240628_result/semen_quality_parameters(将有问题的精液质量筛检生成分组).xlsx") ##导出筛选的变量  

str(df2_1_4)

length(unique(df2_1_4$number))    ####总共有1487人


df2_1_4 %>% 
  dplyr::filter(group_semen==0) -> df2_1_5  ###将死精或无精的样本删除

str(df2_1_5)

length(unique(df2_1_5$number))    ####总共有1485人

6608-6596  ###12份精液质量有问题

df2_1_4  %>% 
  dplyr::filter(group_semen==1) -> df2_1_6   ###

length(unique(df2_1_5$number))   #11人份精液样本有问题，总共12个重复精液样本，但是这些人可能还有其他重复测量精液样本

str(df2_1_6)
# view(df2_1_3)
str(df2_1_6)

df2_1_5$density...8[which(df2_1_5$density...8 ==0)] <- 0.01

# write.xlsx(df2_7,"result/20240628_Hyperten_sperm_2th/20240628_result/semen_quality_parameters(将有问题的生成分组).xlsx") ##导出筛选的变量  

str(df2_1_5)

length(unique(df2_1_5$number))   #1487份人群精液样本，重复测量6608次




###=================数据集合并以及变量筛选================================

###==将一般资料的数据与精液质量参数进行匹配


df2_1_5  %>% 
  dplyr::select(number:QJ...15,datemq,month,season_class) -> df2_1_7  ###挑选出重复测量精液质量的数据

str(df2_1_7)

length(unique(df2_1_6$number))

summary(df2_1_7)



df_1164 <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_repeatsemen_final_1164.xlsx")

str(df_1164)

names(df2_1_7)


###1164匹配的精液样本数量
df2_1_7 %>% 
  left_join(df_1164,by="number")  -> df_merger


names(df_merger)


###1164匹配的精液样本数据，删除缺失值之后的数量
df_merger %>% 
  dplyr::filter(!mtDNAcn=="NA") -> df_merger_noNA


  
  
  
  
  
df_merger_noNA %>% 
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

  


####精浆中PFAS的pool样本 

df_1053 <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_sperm_PFOA_tidy_4_1053.xlsx")

colSums(is.na(df_1053))

df_background_repeatsemen_final2 %>% 
  left_join(df_1053,by="number") -> df_1053_merger
  


df_1053_merger %>% 
  dplyr::filter(!PFOA=="NA")  -> df_1053_merger_noNA

table(df_1053_merger_noNA$interval)
