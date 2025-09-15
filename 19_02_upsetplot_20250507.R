# no_source() 

no_source() 
rm(list = ls())
# setwd(r4projects::get_project_wd())
source("R/100-tools.R")

library(tidyverse)
library(tidymass)

library(scales)



pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts,scales)

# get_project_wd()
# masstools::setwd_project()

getwd() 


###====================读入代谢组学数据===============================

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



####>>>>>>>>读入混合物计算后的P值------------------------------
# write.xlsx(results_df_mix_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/results_df_mix_p_value.xlsx")
# 
# pfas_meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q_PFAS.xlsx")

# ###取出精浆代谢与端粒长度相关的代谢物
# meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q.xlsx")
# 
# meta_STL %>% 
#   filter(q_value<0.2) -> meta_STL_1
# 
# meta_STL_1
# 
# names(meta_STL_1)


###取出PFOA混合暴露与精液代谢相关的代谢物
pfas_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_only_p_value_mix_liner.xlsx")

pfas_meta

271*6

###########################---------------PFOA.beta------------################
# PFOS.beta , PFOA.beta, PFDA.beta, PFUdA.beta, L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PFOA.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
# names(meta_STL_1)
names(pfas_meta_1)


pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_PFOA     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg_PFOA)

names(pfas_meta_STL_p_kegg_PFOA)

colSums(is.na(pfas_meta_STL_p_kegg_PFOA))




#################################PFOS#########################################


# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PFOS.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 

names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_PFOS     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg_PFOS)

names(pfas_meta_STL_p_kegg_PFOS)

colSums(is.na(pfas_meta_STL_p_kegg_PFOS))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>





#################################PFDA.beta#########################################


# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PFDA.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_PFDA    ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg_PFDA)

names(pfas_meta_STL_p_kegg_PFDA)

colSums(is.na(pfas_meta_STL_p_kegg_PFDA))



library("magrittr")


#################################PFUdA.beta#########################################

# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PFUdA.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_PFUdA     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg_PFUdA)

names(pfas_meta_STL_p_kegg_PFUdA)

colSums(is.na(pfas_meta_STL_p_kegg_PFUdA))



#################################L_PFHxS.beta#########################################


# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "L_PFHxS.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
# names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_L_PFHxS     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg_L_PFHxS)

names(pfas_meta_STL_p_kegg_L_PFHxS)

colSums(is.na(pfas_meta_STL_p_kegg_L_PFHxS))



#################################PF3ONS_9CL.beta#########################################


# PFOS.beta , PFOA.beta, PFDA.beta , PFUdA.beta , L_PFHxS.beta, PF3ONS_9CL.beta 

pfas_meta %>% 
  filter(var.names.x == "PF3ONS_9CL.beta") %>%      # PFUdA.beta  7个代谢物；  L_PFHxS.beta### 只有一个代谢物  PF3ONS_9CL.beta  4个代谢物
  filter(p_value <0.05) ->  pfas_meta_1 
# 
# names(meta_STL_1)
names(pfas_meta_1)



pfas_meta_1  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_PF3ONS_9CL    ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg_PF3ONS_9CL)

names(pfas_meta_STL_p_kegg_PF3ONS_9CL)

colSums(is.na(pfas_meta_STL_p_kegg_PF3ONS_9CL))



pfas_meta_STL_p_kegg_PFOA %>% 
  rbind(pfas_meta_STL_p_kegg_PFOA,pfas_meta_STL_p_kegg_PFOS,pfas_meta_STL_p_kegg_PFDA,pfas_meta_STL_p_kegg_PFUdA,pfas_meta_STL_p_kegg_L_PFHxS,
        pfas_meta_STL_p_kegg_PF3ONS_9CL) -> pfas_meta_STL_p_kegg_ALL

pfas_meta_STL_p_kegg_ALL








###取出精浆代谢与端粒长度相关的代谢物
meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q.xlsx")

meta_STL %>% 
  filter(q_value<0.2) -> meta_STL_1

meta_STL_1 %>% 
  mutate(STL=rep("STL",62)) -> meta_STL_2

names(meta_STL_2)

meta_STL_2 %>% 
  dplyr::select(meta,STL) -> meta_STL_STL




###取出PFAS混合暴露与精液代谢相关的代谢物
pfas_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value_mix_liner.xlsx")

pfas_meta %>% 
  filter(q_value<0.2) -> pfas_meta_1 

names(meta_STL)
names(pfas_meta)

pfas_meta_1%>% 
  mutate(PFAS=rep("PFAS",79)) ->  pfas_meta_2  


names(pfas_meta_2)

pfas_meta_2 %>% 
  dplyr::select(metabolite,PFAS) -> pfas_meta_PFAS





####将长数据转换为宽数据


pfas_meta_STL_p_kegg_ALL %>% 
   dplyr::select(metabolite,var.names.x) %>% 
  dplyr::rename(PFAS=var.names.x) %>% 
  rbind(pfas_meta_PFAS) %>%     ###结合上混合暴露的PFAS
  rbind(meta_STL_STL %>% 
          dplyr::rename(metabolite=meta,PFAS=STL)) %>% 
  mutate(value=rep(1,208)) %>% 
  as_tibble()-> pfas_meta_STL_p_kegg_all_long
  

pfas_meta_STL_p_kegg_all_long



  
#   
# pivot_wider(pfas_meta_STL_p_kegg_all_long, names_from = var.names.x, values_from = value)  -> pfas_meta_STL_p_kegg_all_wide
#   
#   
#   
# pfas_meta_STL_p_kegg_all_wide %>% 
#   as_tibble()
#   





###将两列数据转换为一列作为行变量，另一列作为列变量的数据框
library(reshape2)

dcast(pfas_meta_STL_p_kegg_all_long, metabolite ~PFAS, value.var = "value")  -> pfas_meta_STL_p_kegg_all_wide 

pfas_meta_STL_p_kegg_all_wide 

pfas_meta_STL_p_kegg_all_wide %>% 
  dplyr::rename(PFHxS=L_PFHxS.beta,"6:2 Cl-PFESA"=PF3ONS_9CL.beta,PFOA=PFOA.beta,PFOS=PFOS.beta,PFUdA=PFUdA.beta) -> pfas_meta_STL_p_kegg_all_wide_final
  
type(pfas_meta_STL_p_kegg_all_wide_final)
names(pfas_meta_STL_p_kegg_all_wide_final)







#####读入脂质代谢组学数据  


####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2

str(df_unique_2)   ###一共有357个确定的代谢物


######>>>>>>删除ICC<0.3的代谢物  

df_icc <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc_20250222_number.xlsx")
####重复测量的32个同样样本的代谢物特征的一致性

df_icc %>% 
  filter(icc>=0.3) -> df_icc_1   ###筛选出icc>0.3的代谢物特征


df_icc_1$metabolites


df_unique_2 %>% 
  filter(metabolite %in% df_icc_1$metabolites) -> df_unique_3


str(df_unique_3)


pfas_meta_STL_p_kegg_all_wide_final


pfas_meta_STL_p_kegg_all_wide_final %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_all_wide_final_2



pfas_meta_STL_p_kegg_all_wide_final_2 %>% 
  dplyr::mutate(class_material=case_when(class_material=="脂质"~"Lipids", 
                                         class_material=="氨基酸"~ "Amino acids and their derivatives",
                                         class_material=="核苷酸"~"Nucleotides and their derivatives",
                                         class_material=="其他代谢物"~"Other metabolites",
                                         class_material=="糖类"~ "Carbohydrates",
                                         class_material=="有机酸及其衍生物"~"Organic acids and their derivatives")
  ) -> pfas_meta_STL_p_kegg_all_wide_final_3


pfas_meta_STL_p_kegg_all_wide_final_3


names(pfas_meta_STL_p_kegg_all_wide_final_3)




################################################################################

library(ggplot2)
library(ComplexUpset)
# install.packages("ggplot2movies")
library(ggplot2movies)




movies = as.data.frame(pfas_meta_STL_p_kegg_all_wide_final_3)
head(movies, 3)

movies$mpaa

genres = colnames(movies)[2:8]
genres


movies[genres] = movies[genres] == 1

movies

t(head(movies[genres], 3))



set_size(8, 3)
upset(movies, genres, name='genre', width_ratio=0.1)



  
  


set_size(8, 3)
upset(
  movies, genres,
  width_ratio=0.1,
  min_size=10,
  mode='inclusive_union',
  base_annotations=list('Size'=(intersection_size(counts=FALSE, mode='inclusive_union'))),
  intersections='all',
  max_degree=3
)



set_size(8, 6)

upset(
  movies, genres, name='genre', width_ratio=0.1, min_size=10,
  queries=list(
    upset_query(
      intersect=c('PFAS', 'STL'),
      color='red',
      fill='red',
      only_components=c( 'Intersection size')
    ),
    upset_query(
      set='PFAS',
      fill='blue'
    )
  )
)

set_size(5, 3)

upset(
  movies, genres,
  # min_size=2,
  width_ratio=0.3,
  queries=list(upset_query(intersect=c('PFAS', 'STL'), color='red')),
  mode='inclusive_intersection',
  base_annotations=list('Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))),
  set_sizes=(
    upset_set_size(
      position='left'
    )
  ),
  # moves legends over the set sizes
  guides='over'
)


genres

upset(
  movies,c("STL","PFAS","PFOA","PFOS","PFHxS","6:2 Cl-PFESA","PFUdA"),
  # min_size=2,
  width_ratio=0.3,
  queries=list(upset_query(intersect=c('PFAS', 'STL'), color='red',fill="red")),
  mode='inclusive_intersection',
  base_annotations=list('Size'=(intersection_size(counts=TRUE, mode='inclusive_intersection'))),
  set_sizes=(
    upset_set_size(
      position='left'
    )
  ),
  # moves legends over the set sizes
  guides='over'
)


# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
movies

names(movies)

set_size(8, 5)

upset(
  movies, genres,
  name='genre',
  # group_by='sets',
  width_ratio=0.3,
  height_ratio = 0.7,
  mode='inclusive_intersection',
  # queries=list(upset_query(intersect=c('PFAS', 'STL'), color='red')),
  queries=list(
    upset_query(
      intersect=c('PFAS', 'STL'),
      color='red',
      fill='red',
      only_components=c('intersections_matrix', 'Intersection size')
    )
  ),
  base_annotations=list('Size'=(intersection_size(counts=TRUE, mode='inclusive_intersection', mapping=aes(fill=class_material))) + 
                                  scale_fill_manual(
                                    values = c("Amino acids and their derivatives"="#ff7f00", 
                                               "Carbohydrates"="#984ea3",
                                               "Lipids"="#e41a1c", 
                                               "Nucleotides and their derivatives"="#4daf4a", 
                                               "Organic acids and their derivatives"= "#377eb8",
                                               "Other metabolites"="#f781bf"))
                                ),
  set_sizes=(
    upset_set_size( )
    + geom_text(aes(label=..count..), hjust=1, stat='count')
    # you can also add annotations on top of bars:
    + annotate(geom='text', label='@', x='PFAS', y=80, color='white', size=3)
    + expand_limits(y=95)
    + theme(axis.text.x=element_text(angle=0))
  ),
  # themes=upset_default_themes(text=element_text(color='red')),
  themes=upset_modify_themes(
    list(
      'intersections_matrix'=theme(text=element_text(size=14,colour = "black",family = "Arial",face = "bold")),
      'overall_sizes'=theme(axis.text.x=element_text(size=14,colour = "black",family = "Arial",face = "bold"))
    ))
) -> p_upset


p_upset




cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/upset_plot/upset_plot.pdf",width = 16, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p_upset # 使用大写罗马表示
par(opar)
dev.off()





