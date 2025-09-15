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





###导出DNA和TL长度数据
lapply(c(names(df_background_DNA_LT__tidy_lg[,11:12])),function(x){
  df_background_DNA_LT__tidy_lg %>% 
    # dplyr::group_by(!!sym(x)) %>%
    summarise(Geo.Mean.Days.Stay=round(exp(mean(log(get(x)))),digits = 6),
              Geo.SD.Days.Stay=round(exp(sd(log(get(x)))),digits = 6),
              ari.mean=round(mean(get(x)),digits = 6),
              ari.sd=round(sd(get(x)),digits = 6),
              round(quantile(get(x),0.50),digits = 6),
              round(quantile(get(x),0.25),digits = 6),
              round(quantile(get(x),0.50),digits = 6),
              round(quantile(get(x),0.75),digits = 6),
              round(quantile(get(x),0.95),digits = 6),
              Maxnumber=max(get(x)),
              Medi=round(median(get(x)),digits = 6),
              IQR=round(IQR(get(x)),digits = 6)) %>% 
    dplyr::mutate(geom_Mean_SD=paste(Geo.Mean.Days.Stay,"(",Geo.SD.Days.Stay,")",sep=""),
                  ari_mean_sd=paste(ari.mean,"(",ari.sd,")",sep=""),
                  medi_IQR=paste(Medi,"(",IQR,")",sep=""))  -> result1
  as.data.frame(rownames(result1)) %>%   ###添加变量名称
    cbind(as.data.frame(result1))   ###将list转化为数据框
})  %>% 
  bind_rows() %>% 
  cbind(c(names(df_background_DNA_LT__tidy_lg[,11:12]))) %>%  ###粘贴行变量名
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/table2_DNA_TL.xlsx")






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


str(df_sperm_PFOA_tidy_3)

df_background_DNA_LT__tidy_lg %>% 
  left_join(df_sperm_PFOA_tidy_3,by="number") ->df_background_DNA_LT_sperm_PFOA 

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


###PFOA,PFOS,PFNA,PFDA,PFUdA,L_PFHxS,PF3OUdS_11CL,PF3ONS_9CL,PFDoA,L_PFHpS,L_PFBS,PFBA


df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(PFOA>=0.00675113081441142)  %>% 
  count()/836*100    ###高于于LOD的99.28




df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(PFOS>=0.0516440006885867)  %>% 
      count()/836*100    ###高于于LOD的90.19



df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(PFNA>=0.045045045045045)  %>% 
  count()/836*100    ###高于于LOD的19.74



df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(PFDA>=0.00397719740156436)  %>% 
  count()/836*100    ###高于于LOD的76.08



df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(PFUdA>=0.00351658656663932)  %>% 
  count()/836*100    ###高于于LOD的84.57


df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter( L_PFHxS>=0.00338104361546264)  %>% 
  count()/836*100    ###高于于LOD的98.21



df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter( PF3OUdS_11CL>=0.000540365286933967)  %>% 
  count()/836*100    ###高于于LOD的57.89



df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(PF3ONS_9CL>=0.000503330369276714)  %>% 
  count()/836*100    ###高于于LOD的99.40


df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(PFDoA>=0.00255493101686255)  %>% 
  count()/836*100    ###高于于LOD的18.42


df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(L_PFHpS>=0.00917992656058752)  %>% 
  count()/836*100    ###高于于LOD的9.57



df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(L_PFBS>=0.00494967827091239)  %>% 
  count()/836*100    ###高于于LOD的27.15


df_background_DNA_LT_sperm_PFOA_noNA %>% 
  dplyr::filter(PFBA >=0.0319216854649925)  %>% 
  count()/836*100    ###高于于LOD的18.42


#### 将PFOA检测限LOD/根号2替代缺失值
df_background_DNA_LT_sperm_PFOA_noNA %>% 
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
  ) -> df_background_DNA_LT_sperm_PFOA_noNA_lod

str(df_background_DNA_LT_sperm_PFOA_noNA_lod)

# lapply(df_background_DNA_LT_sperm_PFOA_noNA_lod[,c(15:26)],.*1000)
# 
# lapply(df_background_DNA_LT_sperm_PFOA_noNA_lod[,c(15:26)], function(x){
#    res1  <- x*1000
#    return(res1)}) -> a
# 
# str(a)



###将数据重新排序整理
df_background_DNA_LT_sperm_PFOA_noNA_lod %>% 
  select(PFOA,PFOS,PFNA,PFDA,PFUdA,L_PFHxS,PF3OUdS_11CL,PF3ONS_9CL,PFDoA,L_PFHpS,L_PFBS,PFBA,mtDNAcn,STL,everything()) -> df_background_DNA_LT_sperm_PFOA_tidy


str(df_background_DNA_LT_sperm_PFOA_tidy)

write.xlsx(df_background_DNA_LT_sperm_PFOA_tidy,"result/20240407_PFOA_sperm/df_background_DNA_LT_tidy_sperm_PFAS.xlsx")


# 
# ####--------------绘制箱线图------------------------
# ####1.将PFAS宽数据转换为长数据
# 
# df_background_DNA_LT_sperm_PFOA_tidy %>% 
#   dplyr::select(PFOA:PFBA) -> df_background_DNA_LT_sperm_PFOA_tidy_0
# 
# 
# 
# df_background_DNA_LT_sperm_PFOA_tidy_0 %>% 
#   pivot_longer(1:12,names_to = "PFAS", values_to = "value")  -> df_background_DNA_LT_sperm_PFOA_tidy_0_1
# 
# 
# 
# library(dplyr)  
# 
# df_filtered <- df_background_DNA_LT_sperm_PFOA_tidy_0_1 %>%
#   group_by(PFAS) %>%
#   mutate(
#     Q1 = quantile(value, 0.25, na.rm = TRUE),
#     Q3 = quantile(value, 0.75, na.rm = TRUE),
#     IQR = Q3 - Q1,
#     lower_bound = Q1 - 1.5 * IQR,
#     upper_bound = Q3 + 1.5 * IQR
#   ) %>%
#   filter(value >= lower_bound & value <= upper_bound) %>%  # 剔除异常值
#   ungroup()
# 
# 
# 
# df_filtered %>% 
#   # group_by(PFAS) %>%
#   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
#   geom_boxplot(position = position_dodge(0.1),outlier.shape = NA) +
#   coord_flip() +
#    theme_bw() +
#   labs(x = "Plasma PFAS concentrations (ng/mL)", y = "Term")+
#   theme(legend.position = "none",
#         legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         # legend.position = c(0.25,0.15),
#         legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         # axis.line = element_line(linewidth = 0.9),
#         # axis.ticks = element_line(linewidth = 0.9),
#         axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
#         axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
#         panel.grid.major = element_line(size = 0.5, linetype = "dashed"),  # 主网格线粗细 , colour = "red"
#         panel.grid.minor = element_line(size = 0.5, linetype = "dashed") # 次网格线粗细,
#         # panel.border = element_rect(size = 1, colour = "black"), # 设置边框线的粗细和颜色
#         # plot.margin = unit(c(1, 1, 1, 1), "cm")  # 设置整个图形的外边距
#   ) 
# 
#   # # labs(title = "Trait Response to Warming", y = "Response to warming (%)", x = "") +
#   # theme(legend.position = "top",
#   #       axis.text.y = element_text(size = 8))
# 
# 




####求均值正负标准差以及几何均数 




lapply(c(names(df_background_DNA_LT_sperm_PFOA_tidy[,1:14])),function(x){
  df_background_DNA_LT_sperm_PFOA_tidy %>% 
    # dplyr::group_by(!!sym(x)) %>%
    summarise(Geo.Mean.Days.Stay=round(exp(mean(log(get(x)))),digits = 6),
              Geo.SD.Days.Stay=round(exp(sd(log(get(x)))),digits = 6),
              ari.mean=round(mean(get(x)),digits = 6),
              ari.sd=round(sd(get(x)),digits = 6),
              round(quantile(get(x),0.50),digits = 6),
              round(quantile(get(x),0.25),digits = 6),
              round(quantile(get(x),0.50),digits = 6),
              round(quantile(get(x),0.75),digits = 6),
              round(quantile(get(x),0.95),digits = 6),
              Maxnumber=max(get(x)),
              Medi=round(median(get(x)),digits = 6),
              IQR=round(IQR(get(x)),digits = 6)) %>% 
    dplyr::mutate(geom_Mean_SD=paste(Geo.Mean.Days.Stay,"(",Geo.SD.Days.Stay,")",sep=""),
                  ari_mean_sd=paste(ari.mean,"(",ari.sd,")",sep=""),
                  medi_IQR=paste(Medi,"(",IQR,")",sep=""))  -> result1
  as.data.frame(rownames(result1)) %>%   ###添加变量名称
    cbind(as.data.frame(result1))   ###将list转化为数据框
})  %>% 
  bind_rows() %>% 
  cbind(c(names(df_background_DNA_LT_sperm_PFOA_tidy[,1:14]))) %>%  ###粘贴行变量名
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/table2_sperm_PFAS.xlsx")



###------------------读入精浆PFOA的信息-----------------------

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


####将PFOS数据与基础数据合并
df_background_DNA_LT__tidy_lg  %>% 
  left_join(df_blood_PFOA_tidy_3,by="number") -> df_background_DNA_LT_blood_PFOA 

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



df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PFOA>=0.00675113081441142)  %>% 
  count()/1021*100    ###高于于LOD的99.80



df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PFOS>=0.0516440006885867)  %>% 
  count()/1021*100    ###高于于LOD的99.51



df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PFNA>=0.045045045045045)  %>% 
  count()/1021*100    ###高于于LOD的19.74



df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PFDA>=0.00397719740156436)  %>% 
  count()/1021*100    ###高于于LOD的76.08



df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PFUdA>=0.00351658656663932)  %>% 
  count()/1021*100    ###高于于LOD的84.57


df_background_DNA_LT_blood_PFOA_noNA%>% 
  dplyr::filter( L_PFHxS>=0.00338104361546264)  %>% 
  count()/1021*100    ###高于于LOD的98.21

df_background_DNA_LT_blood_PFOA_noNA%>% 
  dplyr::filter(diPAP_6_2>=0.00221141088014153)  %>% 
  count()/1021*100    ###高于于LOD的98.21



df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PF3OUdS_11CL>=0.000540365286933967)  %>% 
  count()/1021*100    ###高于于LOD的57.89



df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PF3ONS_9CL>=0.000503330369276714)  %>% 
  count()/1021*100    ###高于于LOD的99.40


df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PFDoA>=0.00255493101686255)  %>% 
  count()/1021*100    ###高于于LOD的18.42


df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(L_PFHpS>=0.00917992656058752)  %>% 
  count()/1021*100    ###高于于LOD的9.57


df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PFTrDA>=0.00586969281940912)  %>% 
  count()/1021*100    ###高于于LOD的9.57



df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(L_PFBS>=0.00494967827091239)  %>% 
  count()/1021*100    ###高于于LOD的27.15


df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(PFBA>=0.0319216854649925)  %>% 
  count()/1021*100    ###高于于LOD的18.42


df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(L_PFDS>=0.00312044934470564)  %>% 
  count()/1021*100    ###高于于LOD的18.42


df_background_DNA_LT_blood_PFOA_noNA %>% 
  dplyr::filter(FOSA_I>=0.00151095441954168)  %>% 
  count()/1021*100    ###高于于LOD的18.42












#### 将PFOA检测限LOD/根号2替代缺失值
df_background_DNA_LT_blood_PFOA_noNA %>% 
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
  )-> df_background_DNA_LT_blood_PFOA_noNA_lod

str(df_background_DNA_LT_blood_PFOA_noNA_lod)




###-----将数据重新排序整理---------------
df_background_DNA_LT_blood_PFOA_noNA_lod %>% 
  select(PFOA:FOSA_I,mtDNAcn,STL,everything()) -> df_background_DNA_LT_blood_PFOA_tidy



str(df_background_DNA_LT_blood_PFOA_tidy)

write.xlsx(df_background_DNA_LT_blood_PFOA_tidy,"result/20240407_PFOA_sperm/df_background_DNA_LT_tidy_blood_PFAS.xlsx")


# 
# 
# 
# ####--------------绘制箱线图------------------------
# ####1.将PFAS宽数据转换为长数据
# 
# df_background_DNA_LT_blood_PFOA_tidy %>% 
#   dplyr::select(PFOA:FOSA_I) -> df_background_DNA_LT_blood_PFOA_tidy_0
# 
# 
# 
# df_background_DNA_LT_blood_PFOA_tidy_0 %>% 
#   pivot_longer(1:16,names_to = "PFAS", values_to = "value")  -> df_background_DNA_LT_blood_PFOA_tidy_0_1
# 
# 
# df_background_DNA_LT_blood_PFOA_tidy_0_1
# 
# # 
# # 
# # library(circlize)
# # circos.par("track.height" = 0.1)
# # circos.initialize(df_background_DNA_LT_blood_PFOA_tidy_0_1$PFAS, x = df_background_DNA_LT_blood_PFOA_tidy_0_1$value)
# # 
# # 
# # 
# # circos.track(df$sectors, y = df$y,
# #              panel.fun = function(x, y) {
# #                circos.text(CELL_META$xcenter, 
# #                            CELL_META$cell.ylim[2] + mm_y(5), 
# #                            CELL_META$sector.index)
# #                circos.axis(labels.cex = 0.6)
# #              })
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(dplyr)  
# 
# df_filtered <- df_background_DNA_LT_blood_PFOA_tidy_0_1 %>%
#   group_by(PFAS) %>%
#   mutate(
#     Q1 = quantile(value, 0.25, na.rm = TRUE),
#     Q3 = quantile(value, 0.75, na.rm = TRUE),
#     IQR = Q3 - Q1,
#     lower_bound = Q1 - 1.5 * IQR,
#     upper_bound = Q3 + 1.5 * IQR
#   ) %>%
#   filter(value >= lower_bound & value <= upper_bound) %>%  # 剔除异常值
#   ungroup()
# 
# 
# 
# df_background_DNA_LT_blood_PFOA_tidy_0_1 %>% 
#   # group_by(PFAS) %>%
#   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
#   geom_boxplot(position = position_dodge(0.1),outlier.shape = NA,coef = 1.5) +
#   coord_flip() +
#   theme_bw() +
#   labs(x = "Plasma PFAS concentrations (ng/mL)", y = "Term")+
#   theme(legend.position = "none",
#         legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         # legend.position = c(0.25,0.15),
#         legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         # axis.line = element_line(linewidth = 0.9),
#         # axis.ticks = element_line(linewidth = 0.9),
#         axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
#         axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
#         panel.grid.major = element_line(size = 0.5, linetype = "dashed"),  # 主网格线粗细 , colour = "red"
#         panel.grid.minor = element_line(size = 0.5, linetype = "dashed") # 次网格线粗细,
#         # panel.border = element_rect(size = 1, colour = "black"), # 设置边框线的粗细和颜色
#         # plot.margin = unit(c(1, 1, 1, 1), "cm")  # 设置整个图形的外边距
#   ) 
# 
# # # labs(title = "Trait Response to Warming", y = "Response to warming (%)", x = "") +
# # theme(legend.position = "top",
# #       axis.text.y = element_text(size = 8))
# 
# 
# 
# 
# # install.packages("ggbreak")
# library(ggbreak)
# 
# library(ggplot2)
# library(ggbreak)
# 
# df_filtered %>%
#   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
#   geom_boxplot(position = position_dodge(0.1), outlier.shape = NA, coef = 1.5) +
#   scale_y_break(c(1, 5), scales = 0.5, ticklabels = c(0, 0.5, 1, "//", 5, 10, 15, 20)) +
#   coord_flip() +
#   theme_bw() +
#   labs(x = "Plasma PFAS concentrations (ng/mL)", y = "Term") +
#   theme(
#     legend.position = "none",
#     legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     panel.grid.major = element_line(size = 0.5, linetype = "dashed"),
#     panel.grid.minor = element_line(size = 0.5, linetype = "dashed")
#   )
# 
# 
# 
# 
# df_filtered %>%
#   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
#   geom_boxplot(position = position_dodge(0.1), outlier.shape = NA, coef = 1.5) +
#   scale_y_break(
#     c(0.8,0.81),  # 断裂区间略大一点确保数据包含
#     scales = 0.5
#   ) +
#   coord_flip() +
#   theme_bw() +
#   labs(x = "Plasma PFAS concentrations (ng/mL)", y = "Term") +
#   theme(
#     legend.position = "none",
#     legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     panel.grid.major = element_line(size = 0.5, linetype = "dashed"),
#     panel.grid.minor = element_line(size = 0.5, linetype = "dashed"),
#     panel.border = element_blank(),                # 
#     panel.spacing = unit(0, "lines")               #
#   )
# 
# 
# library(dplyr)
# 
# # 计算每组中位数
# medians <- df_filtered %>%
#   group_by(PFAS) %>%
#   summarise(med = median(value))
# # 
# # 
# # df_filtered %>%
# #   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
# #   geom_boxplot(position = position_dodge(0.1), outlier.shape = NA, coef = 1.5) +
# #   geom_text(
# #             aes(x = PFAS, y = med, label = round(med, 2)),
# #             inherit.aes = FALSE,
# #             hjust = -0.2,
# #             size = 4,
# #             fontface = "bold") +
# #   scale_y_break(
# #     c(1,1),
# #     scales = 0.5) +
# #   coord_flip() +
# #   theme_bw() +
# #   labs(x = "Plasma PFAS concentrations (ng/mL)", y = "Term") +
# #   theme(
# #     legend.position = "none",
# #     legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     panel.grid.major = element_line(size = 0.5, linetype = "dashed"),
# #     panel.grid.minor = element_line(size = 0.5, linetype = "dashed"),
# #     panel.border = element_blank(),
# #     panel.spacing = unit(0, "lines")
# #   )
# # 
# # 
# # 
# # 
# # df_filtered %>%
# #   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
# #   geom_boxplot(position = position_dodge(0.1), outlier.shape = NA, coef = 1.5) +
# #   geom_text(data = medians,
# #             aes(x = PFAS, y = med, label = round(med, 2)),
# #             inherit.aes = FALSE,
# #             hjust = -0.6,
# #             size = 4,
# #             fontface = "bold") +
# #   scale_y_break(
# #     c(1, 1),
# #     scales = 0.5
# #   ) +
# #   coord_flip() +
# #   theme_bw() +
# #   labs(x = "Plasma PFAS concentrations (ng/mL)", y = "Term") +
# #   theme(
# #     legend.position = "none",
# #     legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     panel.grid.major = element_line(size = 0.5, linetype = "dashed"),
# #     panel.grid.minor = element_line(size = 0.5, linetype = "dashed"),
# #     panel.border = element_blank(),
# #     panel.spacing = unit(0, "lines")
# #   )
# # 
# # 
# # 
# # 
# # 
# # df_filtered %>%
# #   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
# #   geom_boxplot(position = position_dodge(0.1), outlier.shape = NA, coef = 1.5) +
# #   geom_text(data = medians,
# #             aes(x = PFAS, y = med, label = round(med, 2)),
# #             inherit.aes = FALSE,
# #             hjust = -0.6,
# #             size = 4,
# #             fontface = "bold") +
# #   scale_y_break(
# #     c(1, 1),  # ⚠️ 建议改为 c(1, 3) 避免出错
# #     scales = 0.8
# #   ) +
# #   coord_flip() +
# #   theme_bw() +
# #   labs(x = "Plasma PFAS concentrations (ng/mL)", y = "Term") +
# #   theme(
# #     legend.position = "none",
# #     legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
# #     panel.grid.major = element_blank(),      # ✅ 取消主网格线
# #     panel.grid.minor = element_blank(),      # ✅ 取消次网格线
# #     panel.border = element_blank(),
# #     panel.spacing = unit(0, "lines")
# #   ) -> p2
# # 
# 
# 
# 
# 
# library(dplyr)
# 
# # 重新排序 PFAS 因子顺序（根据中位数）
# ordered_levels <- medians %>%
#   arrange(med) %>%
#   pull(PFAS)
# 
# # 将 df_filtered 中 PFAS 设置为有序因子
# df_filtered <- df_filtered %>%
#   mutate(PFAS = factor(PFAS, levels = ordered_levels))
# 
# # 也更新 extremes 和 medians 中的因子
# medians <- medians %>%
#   mutate(PFAS = factor(PFAS, levels = ordered_levels))
# 
# extremes <- extremes %>%
#   mutate(PFAS = factor(PFAS, levels = ordered_levels))
# 
# 
# 
# df_filtered %>%
#   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
#   stat_boxplot(geom="errorbar",width=0.2,size=0.8)+
#   # 箱线图
#   geom_boxplot(position = position_dodge(0.1), outlier.shape = NA, coef = 1.5) +
#   
#   # 中位数文本标注
#   geom_text(data = medians,
#             aes(x = PFAS, y = med, label = round(med, 2)),
#             inherit.aes = FALSE,
#             hjust = -0.2,
#             size = 4,
#             fontface = "bold") +
#   
#   # 最小值线
#   geom_segment(data = extremes,
#                aes(x = PFAS, xend = PFAS, y = min_val, yend = min_val),
#                inherit.aes = FALSE,
#                color = "black",
#                linewidth = 0.8) +
#   # 对数刻度
#   # scale_y_log10() +  
#   
#   # # 最大值线
#   # geom_segment(data = extremes,
#   #              aes(x = PFAS, xend = PFAS, y = max_val, yend = max_val),
#   #              inherit.aes = FALSE,
#   #              color = "black",
#   #              linewidth = 0.8) +
#   
#   # 断轴（注意断轴范围别用 c(1,1)）
#   # scale_y_break(c(1, 1), scales = 0.8) +  
#   
#   # 翻转坐标轴
#   coord_flip() +
#   theme_bw() +
#   
#   # 主题美化
#   labs(x = "Plasma PFAS concentrations (ng/mL)", y = "Term") +
#   theme(
#     legend.position = "none",
#     legend.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.text = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.y = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     axis.title.x = element_text(size = 14, colour = "black", family = "Arial", face = "bold"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     panel.spacing = unit(0, "lines")
#   )
# 
# 
# 
# 
# 
# p2
# 
# ###导出保存为pdf
# cairo_pdf("article/Third_artical_PFAS_DNA/最终结果_20250109/20250203_result/main_result/fig2_20250704/result_PFAS_LOD_fina.pdf",width =16 , height =6, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# p2
# par(opar)
# dev.off()
# 
# 
# 
# 
# 
# 
# library(dplyr)
# 
# # 计算每个PFAS的中位数
# PFAS_medians <- df_filtered %>%
#   group_by(PFAS) %>%
#   summarise(med = median(value)) %>%
#   mutate(position = ifelse(med < 0.8, "left", "right"))
# 
# 
# 
# df_plot <- df_filtered %>%
#   left_join(PFAS_medians, by = "PFAS")
# 
# 
# 
# 
# library(ggplot2)
# library(patchwork)
# 
# # 左侧图：中位数 < 0.8
# p1 <- df_plot %>%
#   filter(position == "left") %>%
#   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
#   geom_boxplot(outlier.shape = NA, coef = 1.5) +
#   scale_y_break(c(0.8, 0.81), scales = 0.5) +  # 你原本的断轴逻辑
#   coord_flip() +
#   theme_bw() +
#   labs(x = NULL, y = "Plasma PFAS (ng/mL)") +
#   theme(legend.position = "none")
# 
# # 右侧图：中位数 >= 0.8
# p2 <- df_plot %>%
#   filter(position == "right") %>%
#   ggplot(aes(x = PFAS, y = value, fill = PFAS)) +
#   geom_boxplot(outlier.shape = NA, coef = 1.5) +
#   scale_y_break(c(0.8, 0.81), scales = 0.5) +  # 你原本的断轴逻辑
#   coord_flip() +
#   theme_bw() +
#   labs(x = NULL, y = NULL) +
#   theme(legend.position = "none",
#         axis.text.y = element_blank(),  # 不显示重复的 y 标签
#         axis.ticks.y = element_blank(),
#         panel.border = element_blank(),
#         panel.spacing = unit(0, "lines"))
# 
# # 拼接图形
# p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 





###




lapply(c(names(df_background_DNA_LT_blood_PFOA_tidy[,1:18])),function(x){
  df_background_DNA_LT_blood_PFOA_tidy %>% 
    # dplyr::group_by(!!sym(x)) %>%
    summarise(Geo.Mean.Days.Stay=round(exp(mean(log(get(x)))),digits = 6),
              Geo.SD.Days.Stay=round(exp(sd(log(get(x)))),digits = 6),
              ari.mean=round(mean(get(x)),digits = 6),
              ari.sd=round(sd(get(x)),digits = 6),
              round(quantile(get(x),0.50),digits = 6),
              round(quantile(get(x),0.25),digits = 6),
              round(quantile(get(x),0.50),digits = 6),
              round(quantile(get(x),0.75),digits = 6),
              round(quantile(get(x),0.95),digits = 6),
              Maxnumber=max(get(x)),
              Medi=round(median(get(x)),digits = 6),
              IQR=round(IQR(get(x)),digits = 6)) %>% 
    dplyr::mutate(geom_Mean_SD=paste(Geo.Mean.Days.Stay,"(",Geo.SD.Days.Stay,")",sep=""),
                  ari_mean_sd=paste(ari.mean,"(",ari.sd,")",sep=""),
                  medi_IQR=paste(Medi,"(",IQR,")",sep=""))  -> result1
  as.data.frame(rownames(result1)) %>%   ###添加变量名称
    cbind(as.data.frame(result1))   ###将list转化为数据框
})  %>% 
  bind_rows() %>% 
  cbind(c(names(df_background_DNA_LT_blood_PFOA_tidy[,1:18]))) %>%  ###粘贴行变量名
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/table2_blood_PFAS.xlsx")






