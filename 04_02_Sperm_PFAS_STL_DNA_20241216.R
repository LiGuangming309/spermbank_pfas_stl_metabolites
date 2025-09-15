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


###-------1.---------单因素线性回归分析--------------------------------
# 
# ###-------------------------粗模型Crude-------------------
# 
# str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)
# 
# lapply(df_background_DNA_LT_sperm_PFOA_tidy_final_lg[,c(23:24)], function(y){
#   lapply(df_background_DNA_LT_sperm_PFOA_tidy_final_lg[,c(11:22)], function(x){
#     lm(y~x,data=df_background_DNA_LT_sperm_PFOA_tidy_final_lg) -> result
#     cox <- summary(result)
#     cbind(cox$coefficients[,c(1,2,3,4)],confint(result)) -> result2
#     β<-round(cox$coefficients[,1],digits = 2)
#     CI95<-paste0(β,"(",round(result2[,5],2),',',round(result2[,6],2),")")
#     percent_β <- (exp(cox$coefficients[,1])-1)*100
#     percent_CI2.5 <- (exp(result2[,5])-1)*100 
#     percent_CI97.5 <- (exp(result2[,6])-1)
#     percent_CI95 <- paste0(round(percent_β,digits = 2)," (",round((exp(result2[,5])-1)*100,2),', ',round((exp(result2[,6])-1)*100,2),")")
#     res <- result %>% tidy %>% data.frame()   ####自己添加的，将其他输出的变量也添加到包里
#     Uni_glm_model <- data.frame('Estimate'=result2[,1],
#                                 "Std. Error"=result2[,2],
#                                 't value' = result2[,3],
#                                 'p value' = result2[,4],
#                                 'CI2.5' = result2[,5],
#                                 'CI97.5' = result2[,6],
#                                 'CI95' = CI95,
#                                 "percent_β"=percent_β,
#                                 "percent_CI2.5"=percent_CI2.5,
#                                 "percent_CI97.5"=percent_CI97.5,
#                                 "percent_CI95"=percent_CI95,
#                                 res)
#     #返回循环函数继续上述操作                     
#     return(Uni_glm_model)})
# }) -> res1  ###将跑的循环中y生成的list结果保存为res1
# 
# 
# 
# ###然后再将生成的list重新保存为独立的excel文件
# lapply(seq_along(res1), function(i) {
#   write.xlsx(x = res1[[i]], file = paste0("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/Crude_model/result_Sperm_DNA_",names(res1) %>% as_tibble() %>% slice(i), ".xlsx"), row.names = FALSE)
# })
# 


####----------------校正模型adjust——model-------------------

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)

lapply(df_background_DNA_LT_sperm_PFOA_tidy_final_lg[,c(23:24)], function(y){
  lapply(df_background_DNA_LT_sperm_PFOA_tidy_final_lg[,c(11:22)], function(x){
    lm(y~x+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,data=df_background_DNA_LT_sperm_PFOA_tidy_final_lg) -> result
    cox <- summary(result)
    cbind(cox$coefficients[,c(1,2,3,4)],confint(result)) -> result2
    β<-round(cox$coefficients[,1],digits = 2)
    CI95<-paste0(β,"(",round(result2[,5],2),',',round(result2[,6],2),")")
    percent_β <- (exp(cox$coefficients[,1])-1)*100
    percent_CI2.5 <- (exp(result2[,5])-1)*100 
    percent_CI97.5 <- (exp(result2[,6])-1)*100    ###这里少乘了100
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
    return(Uni_glm_model)})
}) -> res1  ###将跑的循环中y生成的list结果保存为res1

# names(df_background_DNA_LT_sperm_PFOA_tidy_final_lg[,c(11:22)]) %>% 
# 
# res1[[1]]{res1[[1]]}%>% 
#   
#   slice(2)
# 
# 
# names(res1[[2]]) %>% 
#   as_tibble()

lapply(seq_along(res1), function(i){
    res1[[i]] %>% 
      bind_rows() %>%    ###将list数据生成为提取出为data数据
      filter(term=="x") %>%     ###取出 PFOS等变量 
      cbind(names(res1[[1]]) %>%   ###结合上pfos等变量名称
              as_tibble()) -> Uni_glm_model
    return(Uni_glm_model)}) -> res2




###然后再将生成的list重新保存为独立的excel文件
lapply(seq_along(res1), function(i) {
  write.xlsx(x = res1[[i]] %>% 
               bind_rows() %>%    ###将list数据生成为提取出为data数据
               filter(term=="x") %>%     ###取出 PFOS等变量 
               cbind(names(res1[[1]]) %>%   ###结合上pfos等变量名称
                       as_tibble()), file = paste0("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/adjusted_model/20240525_Sperm_PFAS_DNA/result_Sperm_DNA_aj_2024_",names(res1) %>% as_tibble() %>% slice(i), ".xlsx"), row.names = FALSE)
})


###生成组变量

res2[[1]] %>% 
  cbind(rep(1,times=12) %>% tibble() %>% dplyr::rename(group=".")) -> res_DNA  ###重复生成一列数据



res2[[2]] %>% 
  cbind(rep(2,times=12) %>% tibble() %>% dplyr::rename(group=".")) -> res_STL  ###重复生成一列数据

###将两个数据合并 
 
res_DNA %>% 
  rbind(res_STL) -> res_DNA_STL_meg


### 根据estimate正负号生成一列新变量，并赋值为positive 和negative

# 根据 estimate 列生成新变量

res_DNA_STL_meg %>% 
mutate(est_group=ifelse(res_DNA_STL_meg$estimate > 0, "positive", "negative")) -> res_DNA_STL_meg_1



### 过滤掉LOD太低的数值
res_DNA_STL_meg_1  %>% 
  filter(value %in% c("PFOS","PFOA","PFDA","PFUdA","L_PFHxS","PF3ONS_9CL")) -> res_DNA_STL_meg_2

res_DNA_STL_meg_2$group


res_DNA_STL_meg_2 %>% 
  mutate(group_sp=case_when(
   group== 1~"mtDNAcn",
    group==2~"STL"
  )) -> res_DNA_STL_meg_3


res_DNA_STL_meg_3

# Now plot them
# res_DNA_STL_meg_2 %>% 
 # dplyr::group_by(group) %>% 


# 查绘图颜色

library(ggsci)
library(tidyverse)
library(cowplot)
library("scales") 

show_col(pal_nejm("default", alpha = 1)(8))  ####查看NEJM杂志的配色，并进行着色

names(res_DNA_STL_meg_3)

###绘制以变量为分组，DNA和STL为分组的图形

p1 <- res_DNA_STL_meg_3 %>% 
  mutate(value_term=case_when(
    value=="PFOS" ~"PFOS",
    value=="PFOA" ~"PFOA",
    value=="PFDA" ~"PFDA",
    value=="PFUdA" ~"PFUdA",
    value=="L_PFHxS" ~"PFHxS",
    value=="PF3ONS_9CL" ~"6:2 Cl-PFESA"   ###修改为文章的变量名称
  ))  %>% 
    # filter(value_term=="6:2 Cl-PFESA") %>% 
  ggplot(aes(x=reorder(value_term,percent_β), y=percent_β,group=group_sp,shape=group_sp)) +
  geom_point(aes(color=group_sp),size=5,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin=percent_CI2.5,ymax=percent_CI97.5,color=group_sp), 
                width = 0.2,size  = 1.5,position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed",color = "red", size = 1) +
  # scale_shape_manual(values=c(16,17))+
  
  # scale_y_continuous(limits = c(-11,10)) +
  # coord_flip()+
  # facet_grid(group ~ ., scales = "free_x")+
  theme(legend.text = element_text(size =16,colour = "black",family = "Arial",face = "bold"),
        legend.title = element_text(size =16,colour = "black",family = "Arial",face = "bold"),
        panel.background = element_rect(fill="grey95"),
        # legend.position = "none",    ###去除图列
        legend.direction = "vertical",
        # strip.background = element_rect(),
        # strip.background = element_blank(),
        strip.placement = "outside",  # 将分面标题放置在外部
        strip.text.y = element_text(angle = -90, hjust = 0.5,size =16,colour = "black",family = "Arial",face = "bold"),  ###设置每个分面标题的字体及格式
        text = element_text(size = 16),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        panel.grid = element_blank(),
        axis.ticks.x  = element_blank(),  ###去除X轴坐标
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text=element_text(size =16,colour = "black",family = "Arial",face = "bold"),
        axis.title.y = element_text(size =16,colour = "black",family = "Arial",face = "bold")
        # axis.text.x = element_text(colour =c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF"),fill=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")  ###给X轴字体填充不同的颜色
  )+
  labs(x = "", y = "Percent change [βs (95%)]",size =16,colour = "black",family = "Arial",face = "bold")+
  # scale_color_manual(name = "Correlation", values = c("#BC3C29FF", "#0072B5FF")) +   ###修改颜色分组变量的图列标题 
  scale_color_manual(name="Group", values = c("STL" = "#BC3C29FF", "mtDNAcn" = "#0072B5FF"))+
  scale_size_manual(values = c(2.5, 2.5), name = NULL) +
  scale_shape_manual(name = "Group", values = c(16, 17))            ###修改形状分组的图列标题

p1



###导出保存为pdf
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/adjusted_model/20240525_Sperm_PFAS_DNA/result_Sperm_DNA_aj_20241019_fina.pdf",width =16 , height =6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p1
par(opar)
dev.off()









# dev.off()
###绘图

  res_DNA_STL_meg_3
  
res_DNA_STL_meg_3 %>% 
  mutate(value_term=case_when(
    value=="PFOS" ~"PFOS",
    value=="PFOA" ~"PFOA",
    value=="PFDA" ~"PFDA",
    value=="PFUdA" ~"PFUdA",
    value=="L_PFHxS" ~"PFHxS",
    value=="PF3ONS_9CL" ~"6:2 Cl-PFESA"   ###修改为文章的变量名称
  ))  %>% 
ggplot(aes(x=reorder(value_term,percent_β), y=percent_β,group=group_sp,color=est_group)) +
  geom_point(aes(shape=group_sp),size=5,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=percent_CI2.5,ymax=percent_CI97.5), 
                width = 0.2,size  = 1.5,position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed",color = "red", size = 1) +
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values = c("positive" = "#BC3C29FF", "negative" = "#0072B5FF"))+
  scale_size_manual(values = c(2.5, 2.5), name = NULL) +
  # scale_y_continuous(limits = c(-11,7)) +
  # coord_flip()+
  facet_grid(value_term ~ ., scales = "free_x")+
   theme(legend.text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
         legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        panel.background = element_rect(fill="grey95"),
        # legend.position = "none",    ###去除图列
        legend.direction = "vertical",
        strip.background = element_rect(),
        # strip.background = element_blank(),
        strip.placement = "outside",  # 将分面标题放置在外部
        strip.text.y = element_text(angle = -90, hjust = 0.5,size =14,colour = "black",family = "Arial",face = "bold"),  ###设置每个分面标题的字体及格式
        text = element_text(size = 14),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        panel.grid = element_blank(),
        axis.ticks.x  = element_blank(),  ###去除X轴坐标
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.title.y = element_text(size =14,colour = "black",family = "Arial",face = "bold"),

        axis.text.x = element_blank())+
  labs(x = "", y = "Percent change [βs (95%)]",size =14,colour = "black",family = "Arial",face = "bold") -> p
  
p
###修改分面中每个标题的颜色
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
} 

# tiff
grid.draw(g) 


###导出保存为pdf
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/adjusted_model/20240525_Sperm_PFAS_DNA/result_Sperm_DNA_aj_20241019.pdf",width = 12, height = 16, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g) 

par(opar)
dev.off()



# windowsFonts(A = windowsFont("Calibri")) 
###绘制剂量反应图
tiff("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/adjusted_model/20240525_Sperm_PFAS_DNA/result_Sperm_DNA_aj_20241019.tiff",height = 1280,width =960)
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g) 


par(opar)
dev.off()





####>>>>>>>>>>>>>>>>>提取——校正后的四分位绘图<<<<<<<<<<<<<<<<<<<<<<<<<###################



####---------------四分为间距--------------------------

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg)
# PFOS+PFOA+PFUdA+L_PFHxS+PF3ONS_9CL

####利用循环函数，批量分段重组

lapply(df_background_DNA_LT_sperm_PFOA_tidy_final_lg[,c("PFOS","PFOA","PFUdA","L_PFHxS","PF3ONS_9CL")], function(x){ 
  cut(x,breaks =quantile(x),include.lowest=T,
      labels = c(1,2,3,4)) %>% 
    as_tibble() 
}) %>% 
  bind_cols()  %>%    ####将所有列数据合并在一个数据框中
  set_names(paste(data.frame(colnames(df_background_DNA_LT_sperm_PFOA_tidy_final_lg))[c(11,12,14,15,16),],"_","class",sep = "")) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_class  ####对生成的数据框根据列变量进行重新命名
str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg_class)


cut(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$PF3OUdS_11CL,breaks =quantile(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$PF3OUdS_11CL,probs=c(0.0,0.5,0.75,1.00)),include.lowest=T,
    labels = c(1,2,3))%>% 
  bind_cols() %>%    ####将所有列数据合并在一个数据框中
  set_names(paste(data.frame(colnames(df_background_DNA_LT_sperm_PFOA_tidy_final_lg))[c(22),],"_","class",sep = "")) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_class_2



cut(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$PFDA,breaks =quantile(df_background_DNA_LT_sperm_PFOA_tidy_final_lg$PFDA,probs=c(0.0,0.33,0.66,1.00)),include.lowest=T,
    labels = c(1,2,3))%>% 
  bind_cols() %>%    ####将所有列数据合并在一个数据框中
  set_names(paste(data.frame(colnames(df_background_DNA_LT_sperm_PFOA_tidy_final_lg))[c(13),],"_","class",sep = "")) -> df_background_DNA_LT_sperm_PFOA_tidy_final_lg_class_3

str(df_background_DNA_LT_sperm_PFOA_tidy_final_lg_class_3)


####将原始数据与新生成的分类数据进行合并 

df_background_DNA_LT_sperm_PFOA_tidy_final_lg %>% 
  cbind(df_background_DNA_LT_sperm_PFOA_tidy_final_lg_class,df_background_DNA_LT_sperm_PFOA_tidy_final_lg_class_2, df_background_DNA_LT_sperm_PFOA_tidy_final_lg_class_3) %>% 
  dplyr::select(number:STL,PFOS_class,PFOA_class,PFDA_class,PFUdA_class,L_PFHxS_class,PF3ONS_9CL_class,PF3OUdS_11CL_class)-> df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class

str(df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class)




###------用原始的PFAS计算每个四分位区间大于等于LOD的样本数-----------
df_background_DNA_LT_sperm_PFOA_tidy%>% 
  dplyr::filter(PFOA>=0.00675113081441142)  %>% 
  count()/836*100    ###高于于LOD的99.28



###将低于检测线的转换为0或1分类变量
df_background_DNA_LT_sperm_PFOA_tidy %>% 
  mutate(
    PFOS_lod=case_when(
      PFOS < 0.0516440006885867 ~ 0,
      PFOS >= 0.0516440006885867~ 1
    ),
    PFOA_lod=case_when(
      PFOA < 0.00675113081441142~0,
      PFOA>= 0.00675113081441142~1
    ),
    PFDA_lod=case_when(
      PFDA< 0.00397719740156436~0,
      PFDA>=0.00397719740156436~1
    ),
    PFUdA_lod=case_when(
      PFUdA< 0.00351658656663932~0,
      PFUdA>= 0.00351658656663932~1
    ),
    L_PFHxS_lod=case_when(
      L_PFHxS< 0.00338104361546264~0,
      L_PFHxS>= 0.00338104361546264~1
    ),
    PF3ONS_9CL_lod=case_when(
      PF3ONS_9CL< 0.000503330369276714~0,
      PF3ONS_9CL >= 0.000503330369276714~1
    ),
    PF3OUdS_11CL_lod=case_when(
      PF3OUdS_11CL< 0.000540365286933967~0,
      PF3OUdS_11CL>= 0.000540365286933967~1
    )
  ) -> df_background_DNA_LT_sperm_PFOA_tidy_less_LOD

str(df_background_DNA_LT_sperm_PFOA_tidy_less_LOD)

summary(df_background_DNA_LT_sperm_PFOA_tidy_less_LOD)

####将生成的LOD转换为因变量并连接到原来的数据框

df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class %>% 
  cbind(
    lapply(df_background_DNA_LT_sperm_PFOA_tidy_less_LOD[,c(27:33)], factor)  %>% 
      bind_cols()  %>%    ####将所有列数据合并在一个数据框中
      set_names(paste(data.frame(colnames(df_background_DNA_LT_sperm_PFOA_tidy_less_LOD))[c(27:33),],"_","factor",sep = "")) 
  ) -> df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2

summary(df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2)

library(dplyr)
library(tidyr)

str(df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2)
####将生成每个四分位间距中大于LOD的样本数
# ###PFOS
# df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
#   dplyr::group_by(PFOS_class,PFOS_lod_factor) %>% 
#   tally() 
#   
# 836-82
# 
# ###PFOA
# df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
#   dplyr::group_by(PFOA_class,PFOA_lod_factor) %>% 
#   tally() 
# 
# 836-6
# 
# 
# ###PFDA
# df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
#   dplyr::group_by(PFDA_class,PFDA_lod_factor) %>% 
#   tally() 
# 
# 836-200
# 
# 
# 
# ###PFUdA
# df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
#   dplyr::group_by(PFUdA_class,PFUdA_lod_factor) %>% 
#   tally() 
# 
# 836-129
# 
# ###L_PFHxS
# df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
#   dplyr::group_by(L_PFHxS_class,L_PFHxS_lod_factor) %>% 
#   tally() 
# 836-15
# 
# 
# 
# ###PF3ONS_9CL
# df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
#   dplyr::group_by(PF3ONS_9CL_class,PF3ONS_9CL_lod_factor) %>% 
#   tally() 
# 
# 836-5
# 
# ###PF3OUdS_11CL
# df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
#   dplyr::group_by(PF3OUdS_11CL_class,PF3OUdS_11CL_lod_factor) %>% 
#   tally() 
# 
# 836-352
# 




####将生成每个四分位间距中大于LOD的样本数
###PFOS
df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
  dplyr::group_by(PFOS_class) %>% 
  tally() 

836-82

###PFOA
df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
  dplyr::group_by(PFOA_class) %>% 
  tally() 

836-6


###PFDA
df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
  dplyr::group_by(PFDA_class) %>% 
  tally() 

836-200



###PFUdA
df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
  dplyr::group_by(PFUdA_class) %>% 
  tally() 

836-129

###L_PFHxS
df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
  dplyr::group_by(L_PFHxS_class) %>% 
  tally() 
836-15



###PF3ONS_9CL
df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
  dplyr::group_by(PF3ONS_9CL_class) %>% 
  tally() 

836-5

###PF3OUdS_11CL
df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class_final_2 %>% 
  dplyr::group_by(PF3OUdS_11CL_class) %>% 
  tally() 

836-352







####------------将LOD大于80%的做多元混合线性模型----------------------

###精浆检出率大于80%的有PFOS+PFOA+PFUdA+L_PFHxS+PF3ONS_9CL
##11,12,14,



####----------------校正模型adjust——model-------------------

str(df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class)

lapply(df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class[,c(23:24)], function(y){
  lapply(df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class[,c(25:31)], function(x){
    lm(y~x+age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2,data=df_background_DNA_LT_sperm_PFOA_tidy_final_LOD_Class) -> result
    cox <- summary(result)
    cbind(cox$coefficients[,c(1,2,3,4)],confint(result)) -> result2
    β<-round(cox$coefficients[,1],digits = 2)
    CI95<-paste0(β,"(",round(result2[,5],2),',',round(result2[,6],2),")")
    percent_β <- (exp(cox$coefficients[,1])-1)*100
    percent_CI2.5 <- (exp(result2[,5])-1)*100 
    percent_CI97.5 <- (exp(result2[,6])-1)*100
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
    return(Uni_glm_model)})
}) -> res1  ###将跑的循环中y生成的list结果保存为res1




res1[[1]]



lapply(seq_along(res1), function(i) {
  lapply(seq_along(res1[[i]]), function(j) {
    res1[[i]][[j]] %>%
      as_tibble() %>%
      mutate(dataframe_name = names(res1[[i]])[j])  # 添加数据框名称列
  }) %>% bind_rows()  %>% 
    mutate(list_name = names(res1)[i])  # 添加名称列# 合并内层数据框
}) %>% bind_rows()  -> res_all  # 合并所有外层数据框

write.xlsx(res_all,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/adjusted_model/20240525_Sperm_PFAS_DNA/result_Sperm_DNA_aj_2024_IQR.xlsx")

res_all$dataframe_name
res_all$list_name




###合并外层数据

lapply(seq_along(res1), function(i) {
  res1[[i]] %>% 
    bind_rows() %>% 
    as_tibble() %>% 
    mutate(list_name = names(res1)[i])  # 添加名称列
}) %>% 
  bind_rows()


# 
# ###然后再将生成的list重新保存为独立的excel文件
# lapply(seq_along(res1), function(i) {
#   write.xlsx(x = res1[[i]], file = paste0("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/adjusted_model/20240525_Sperm_PFAS_DNA/四分位间距/result_Sperm_DNA_4_class_aj_",names(res1) %>% as_tibble() %>% slice(i), ".xlsx"), row.names = FALSE)
# })
# 



res_all$dataframe_name

names(res_all)

res_all$list_name

str(res_all)

res_all %>% 
   dplyr::filter(term %in% c("x2", "x3", "x4")) -> res_all_1   ###取出含有x2\x3\x4

names(res_all_1)

str(res_all_1)
###根据dataframe_name给每组添加一行X1


res_all_1 %>%
  group_by(dataframe_name) %>%
  do({
    df <- .  # 当前组的数据框
    zero_row <- df[1, ] %>%
      mutate(across(everything(), ~ 0)) %>%
      mutate(term = "x1")  # 修改 term 为 "x1"
    
    # 返回原数据框与新行的合并
    rbind(df, zero_row)
  }) %>%
  ungroup()  # 取消分组


###>>>>>>>#根据两个分组变量，对每个组的每个类编分别生成一个x1,合并所有组的结果,并且每组的填充相应类别的组值

res_all_1 %>%
  group_by(dataframe_name, list_name) %>%
  do({
    df <- .  # 当前组的数据框
    
    # 创建一行所有变量均为0的x1行，填充 dataframe_name 和 list_name
    zero_row <- df[1, ] %>%
      mutate(across(everything(), ~ 0)) %>%
      mutate(term = "x1",  # 修改 term 为 "x1"
             dataframe_name = unique(df$dataframe_name),  # 填充 dataframe_name
             list_name = unique(df$list_name))  # 填充 list_name
    
    # 返回原数据框与新行的合并
    rbind(df, zero_row)
  }) %>%
  ungroup()  -> res_all_2  # 取消分组


res_all_2$dataframe_name
res_all_2$list_name
# 
# res_all_1 %>%
#   group_by(dataframe_name, list_name) %>%
#   do({
#     df <- .  # 当前组的数据框
#     zero_row <- df[1, ] %>%
#       mutate(across(everything(), ~ 0)) %>%
#       mutate(term = "x1")  # 修改 term 为 "x1"
#     
#     # 返回原数据框与新行的合并
#     rbind(df, zero_row)
#   }) %>%
#   ungroup()  -> res_all_2  # 取消分组


write.xlsx(res_all_2,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/adjusted_model/20240525_Sperm_PFAS_DNA/result_Sperm_DNA_aj_2024_IQR_final.xlsx")


# 查绘图颜色

library(ggsci)
library(tidyverse)
library(cowplot)
library("scales") 

show_col(pal_nejm("default", alpha = 1)(8))  ####查看NEJM杂志的配色，并进行着色



###绘制以变量为分组，DNA和STL为分组的图形
####>>>>>>对变量重新赋值及命名 
names(res_all_2) 
str(res_all_2) 



res_all_2$dataframe_name

###>>>>>删除dataframe中PF3OUdS_11CL_class ，主要结果为显示这个分组变量
res_all_2 %>% 
  filter(!dataframe_name=="PF3OUdS_11CL_class") -> res_all_3


res_all_3  %>% 
  dplyr::mutate(value_term=case_when(
    dataframe_name=="PFOS_class"~"PFOS",
    dataframe_name=="PFOA_class"~"PFOA",
    dataframe_name=="PFDA_class"~"PFDA",
    dataframe_name=="PFUdA_class"~"PFUdA",
    dataframe_name=="L_PFHxS_class"~"PFHxS",
    dataframe_name=="PF3ONS_9CL_class"~"6:2 Cl-PFESA"   ###修改为文章的变量名称
  ),
  term_Q=case_when(
    term=="x1"~ "Q1",
    term=="x2"~ "Q2",
    term=="x3"~ "Q3",
    term=="x4"~ "Q4"
  )) %>% 
  dplyr::rename("group"="list_name")-> res_all_4

res_all_4$dataframe_name
res_all_4$value_term

res_all_4$term_Q
res_all_4 %>% 
  as.data.frame() -> res_all_5

names(res_all_5)

res_all_5$group
res_all_5$value_term
res_all_5$term_Q
  # filter(value_term=="6:2 Cl-PFESA") %>% 
ggplot(res_all_5,aes(x=value_term, y=percent_β,group=term_Q,color=term_Q)) +
  geom_point(aes(shape=term_Q),size=5,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin=percent_CI2.5,ymax=percent_CI97.5), 
                width = 0.2,size  = 1.5,position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed",color = "red", size = 1) +
  facet_grid(group ~ ., scales = "free_y")+
  scale_shape_manual(values=c(16,17,18,19))+
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"))+
  scale_size_manual(values = c(2.5, 2.5,2.5,2.5), name = NULL) +
  # scale_y_continuous(limits = c(-11,10)) +
  # coord_flip()+
 
  theme(legend.text = element_text(size =16,colour = "black",family = "Arial",face = "bold"),
        legend.title = element_text(size =16,colour = "black",family = "Arial",face = "bold"),
        panel.background = element_rect(fill="grey95"),
        # legend.position = "none",    ###去除图列
        legend.direction = "vertical",
        # strip.background = element_rect(),
        # strip.background = element_blank(),
        strip.placement = "outside",  # 将分面标题放置在外部
        strip.text.y = element_text(angle = -90, hjust = 0.5,size =16,colour = "black",family = "Arial",face = "bold"),  ###设置每个分面标题的字体及格式
        text = element_text(size = 16),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        panel.grid = element_blank(),
        axis.ticks.x  = element_blank(),  ###去除X轴坐标
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text=element_text(size =16,colour = "black",family = "Arial",face = "bold"),
        axis.title.y = element_text(size =16,colour = "black",family = "Arial",face = "bold")
        # axis.text.x = element_text(colour =c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF"),fill=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")  ###给X轴字体填充不同的颜色
  )+
  labs(x = "", y = "Percent change [βs (95%)]",size =14,colour = "black",family = "Arial",face = "bold")+
  scale_color_manual(name = "Quantile",values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF"))+   ###修改颜色分组变量的图列标题
  scale_shape_manual(name = "Quantile", values =c(16,17,18,19))   -> p3          ###修改形状分组的图列标题


p3












###导出保存为pdf
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/adjusted_model/20240525_Sperm_PFAS_DNA/result_Sperm_DNA_aj_20241019_IQR_fina.pdf",width =16 , height =9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p3
par(opar)
dev.off()










p
###修改分面中每个标题的颜色
g <- ggplot_gtable(ggplot_build(p3))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#BC3C29FF","#0072B5FF")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
} 

# tiff
grid.draw(g) 






# windowsFonts(A = windowsFont("Calibri")) 
###绘制剂量反应图
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/adjusted_model/20240525_Sperm_PFAS_DNA/result_Sperm_DNA_aj_20241019_IQR_fina_color.pdf",width =16 , height =9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
###修改分面中每个标题的颜色
g <- ggplot_gtable(ggplot_build(p3))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#BC3C29FF","#0072B5FF")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
} 

# tiff
grid.draw(g) 


par(opar)
dev.off()








