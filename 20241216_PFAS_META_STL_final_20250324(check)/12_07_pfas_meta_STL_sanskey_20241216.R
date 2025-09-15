

pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts)
# get_project_wd()
# masstools::setwd_project()
rm(list = ls())
getwd()


###>>>>>>>>>>>>>>>>>>>>>>1.读入经过BHRMA筛选后的pfas混合物<<<<<<<<<<<<<<<<<<

df_pfas_mix_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value_mix_liner.xlsx")


###1.提取经过BH校正之后q_value<0.2的代谢物 
###
names(df_pfas_mix_meta)

df_pfas_mix_meta %>% 
  filter(q_value<0.2) -> df_pfas_mix_meta_q   ###总共获得95个代谢物质   


df_pfas_mix_meta_q 



####>>>>>>>>读入代谢物匹配的代谢通路变量名称#############################


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


df_unique_3$metabolite

# mwas_reduced$metabolite


df_pfas_mix_meta_q %>% 
  left_join(df_unique_3,by="metabolite") -> df_pfas_mix_meta_q_kegg


names(df_pfas_mix_meta_q_kegg)


colSums(is.na(df_pfas_mix_meta_q_kegg))

str(df_pfas_mix_meta_q_kegg)   ####PFAS相关的只有79个代谢物   




###>>>>>>>>>>>>>>>>>>>>>>2.读入代谢与端粒有统计学意义的化合物<<<<<<<<<<<<<<<<<<


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

###取出精浆代谢与端粒长度相关的代谢物
meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q.xlsx")

meta_STL %>% 
  filter(q_value<0.2) -> meta_STL_1

meta_STL_1

###取出PFAS混合暴露与精液代谢相关的代谢物
pfas_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value_mix_liner.xlsx")

pfas_meta %>% 
  filter(q_value<0.2) -> pfas_meta_1 

# names(meta_STL)
names(pfas_meta_1)

pfas_meta_1 %>% 
  left_join(meta_STL_1,by=c("metabolite"="meta")) -> pfas_meta_stl

colSums(is.na(pfas_meta_stl))

pfas_meta_stl %>% 
  filter(!Estimate=="NA") -> pfas_meta_stl_all 



str(pfas_meta_stl_all)

names(pfas_meta_stl_all)

# 
# 
pfas_meta_stl_all  %>%
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg     ####取交集之后的最终筛出的代谢物质
# 
# meta_STL_p_kegg_meta$q_value
# 
# # write.csv(STL_meta_p_kegg_meta,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/STL_meta_fdr_20241122.csv")
# 
# names(meta_STL_p_kegg_meta)
# 
# 
# meta_STL_p_kegg_meta$metabolite
# 
# meta_STL_p_kegg_meta$name

#######根据相关性系数生成正负两组与端粒有关的的化合物

pfas_meta_STL_p_kegg  %>% 
  mutate(group_STL=case_when(
    Estimate<0~"Negative",
    Estimate>0~"Positive"
  )) -> meta_STL_p_kegg_meta_2

str(meta_STL_p_kegg_meta_2)

names(meta_STL_p_kegg_meta_2)


###>>>>>>>3.根据均有统计学意义的指标，反推有统计意义的化合物对应的每一个pfas混合物的值<<< 
# 
####对所有化合物对应的PFAS进行桑基图
pfas_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_only_p_value_mix_liner.xlsx")
names(pfas_meta)

pfas_meta$var.names.x
pfas_meta$meta
meta_STL_p_kegg_meta_2$metabolite

length(pfas_meta$metabolite)



meta_STL_p_kegg_meta_2 %>%
  left_join(pfas_meta,by=c("metabolite"="metabolite")) -> pfas_meta_STL_p_kegg_meta_3


# write.xlsx(pfas_meta_STL_p_kegg_meta_3,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/sankey_all.xlsx")

# 
####只筛出有意义的进行桑基图
pfas_meta_STL_p_kegg_meta_3

names(pfas_meta_STL_p_kegg_meta_3)

pfas_meta_STL_p_kegg_meta_3$p_value_estimate.x

pfas_meta_STL_p_kegg_meta_3 %>%
  filter(p_value_estimate.x<0.2) ->pfas_meta_STL_p_kegg_meta_3_1

pfas_meta_STL_p_kegg_meta_3_1

pfas_meta_STL_p_kegg_meta_3_1$var.names.x



# write.xlsx(pfas_meta_STL_p_kegg_meta_3_1,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/sankey.xlsx")



###将生成混合PFAS数据，结合到桑基图中

meta_STL_p_kegg_meta_2

names(meta_STL_p_kegg_meta_2)

str(meta_STL_p_kegg_meta_2)  ###根据行变量个数，查询行数

meta_STL_p_kegg_meta_2 %>%
  mutate(var.names.x=rep("PFAS_mixexposure",18)) -> pfas_meta_mix_2

names(pfas_meta_mix_2)
str(pfas_meta_mix_2)


pfas_meta_mix_2 %>% 
  select(var.names.x,name,group_STL) %>% 
  dplyr::rename("PFAS"="var.names.x","Metabolite"="name","STL"="group_STL") -> df1


df1


names(pfas_meta_STL_p_kegg_meta_3_1)
str(pfas_meta_STL_p_kegg_meta_3_1)



pfas_meta_STL_p_kegg_meta_3_1 %>% 
  select(var.names.x,name,group_STL) %>% 
  dplyr::rename("PFAS"="var.names.x","Metabolite"="name","STL"="group_STL")-> df2



###
df1 %>% 
  rbind(df2) -> df3

df3$PFAS
df3 %>% 
dplyr::mutate(PFAS=case_when(
  PFAS=="PFOS.beta"~"PFOS",
  PFAS=="PFOA.beta"~"PFOA",
  PFAS=="PFDA.beta"~"PFDA",
  PFAS=="L_PFHxS.beta"~"PFHxS",
  PFAS=="PFUdA.beta"~"PFUdA",
  PFAS=="PF3ONS_9CL.beta"~ "6:2 Cl-PFESA",
  PFAS=="PFAS_mixexposure"~"PFAS mixture"
)) -> df4

df4



# # 
# # pfas_meta_mix <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/results_df_mix_p_value_mix.xlsx")
# # names(pfas_meta_mix)
# # str(pfas_meta_mix)
# # pfas_meta_mix$var.names
# 
# pfas_meta_mix %>%
#   mutate(var.names.x=rep("PFAS_mixexposure",414)) -> pfas_meta_mix_2
# 
# 
# meta_STL_p_kegg_meta_2 %>%
#   left_join(pfas_meta_mix_2,by=c("metabolite"="metabolite")) -> pfas_meta_STL_p_kegg_meta_4
# 
# names(pfas_meta_STL_p_kegg_meta_4)
# 
# pfas_meta_STL_p_kegg_meta_4$q_value.y
# 
# # 
# 

###>>>>>>>4.绘制桑基图<<< 

library(scales)
library(ggsci)
library(sankeywheel)

install.packages("ggsankey")

# install.packages("devtools")
# devtools::install_github("davidsjoberg/ggsankey")

###
library(ggsankey)
library(dplyr)
library(ggplot2)

df4$PFAS
pg2 <- df4 %>%
  make_long("PFAS","Metabolite","STL")

phenotype_class_color = c(
  "Positive" = "#BC3C29FF",
  "Negative" = "#0072B5FF")


df4$PFAS
custom_palette_2=c( "PFAS mixture" = ggsci::pal_jama()(n=7)[4],
  "PFOS" = ggsci::pal_jama()(n=7)[1],
  "PFOA" = ggsci::pal_jama()(n=7)[2],
  "PFDA" = ggsci::pal_jama()(n=7)[3],
  "PFUdA" = ggsci::pal_jama()(n=7)[5],
  "PFHxS" = ggsci::pal_jama()(n=7)[6],
  "6:2 Cl-PFESA" = ggsci::pal_jama()(n=7)[7]
                   )


show_col(ggsci::pal_jama()(n=7))
show_col(ggsci::pal_nejm()(n=7))

show_col(ggsci::pal_d3()(n=7))


# df4$Metabolite[1:33] %in% c(viridis(33))[2]

library(viridis)
df4$Metabolite

omics_color = c(
  "Bisnorcholic acid" = c(viridis(25))[1],
  "tetracosatetraenoate" = c(viridis(25))[2],
  "Tetracosapentaenoic acid (24:5n-3)" = c(viridis(25))[3],
  
  "12-Oxo-20-carboxy-leukotriene B4" = c(viridis(25))[4],
  "(3E,5Z,11Z)-Pentadeca-3,5,11-trienedioylcarnitine" = c(viridis(25))[5],
  "9-Hydroxy-16-oxooctadeca-10,12,14-trienoylcarnitine" = c(viridis(25))[6],
  
  "3-Hydroxyhex-4-enoylcarnitine" = c(viridis(25))[7],
  "13,14-dihydro-15-keto-tetranor Prostaglandin F1beta" = c(viridis(25))[8],
  "Testosterone" = c(viridis(25))[9],
  
  "Linoleic acid" = c(viridis(25))[10],
  "L-Acetylcarnitine" = c(viridis(25))[11],
  "O-Propanoylcarnitine" = c(viridis(25))[12],
  
  "4-Hydroxy-4-methyl-2-oxoadipate" = c(viridis(25))[13],
  "Prostaglandin E1" = c(viridis(25))[14],
  "2-Hydroxyestradiol" = c(viridis(25))[15],
  
  "Tetrahydrocortisol" = c(viridis(25))[16],
  "Eicosapentaenoic acid" = c(viridis(25))[17],
  "3-Hydroxyhexadecanoylcarnitine" = c(viridis(25))[18]
  # 
  # "2-Methoxy-17beta-estradiol" = c(viridis(25))[19],
  # "Tetrahydrocortisol" = c(viridis(25))[20],
  # "Eicosapentaenoic acid" = c(viridis(25))[21],
  # 
  # "Docosahexaenoic acid" = c(viridis(25))[22],
  # "3-Hydroxyhexadecanoylcarnitine" = c(viridis(25))[23],
  # "19-Noraldosterone" = c(viridis(25))[24],
  # 
  # "10-Hydroperoxydocosahexaenoic acid" = c(viridis(25))[25]
 
)


# "L-Acetylcarnitine" = c(viridis(28))[1],omics_color= c(wes_palette("Zissou1", 28, type = "continuous"))
# 如果未安装，先安装包
if (!requireNamespace("wesanderson", quietly = TRUE)) {
  install.packages("wesanderson")
}

# 加载包
library(wesanderson)
# 获取Zissou1色带，连续的颜色
colors <- c(wes_palette("Zissou1",18, type = "continuous"))

# 给代谢物变量赋予颜色

viridis(33) ## viridis 主题中提取10个颜色
inferno(28) ## inferno 主题中提取10个颜色


plot(1:28, col=viridis(10),cex=4,pch=20)

barplot(rep(28,28),col=viridis(28),space=0,border = viridis(28))




# "Metabolite" = ggsci::pal_jama()(n=7)[1],
# "Lipidome" = ggsci::pal_aaas()(n=7)[2],

# df5$`Exogenous metals`
# 
# custom_palette_2=c(
#   "Br Intensity"="#00A1D5FF",
#   "Rh Intensity"="#00A1D5FF",
#   "Hg Intensity"="#00A1D5FF",
#   "Rb Intensity"="#00A1D5FF",
#   "Sr Intensity"="#00A1D5FF",
#   "Y Intensity"="#00A1D5FF",
#   "Zr Intensity"="#00A1D5FF",
#   "Sn Intensity"="#00A1D5FF",
#   "Sb Intensity"="#00A1D5FF",
#   "W Intensity"="#00A1D5FF",
#   "Pb Intensity"="#00A1D5FF",
#   "Bi Intensity"="#00A1D5FF",
#   "Rb Percentage"="#FF7F0EFF",
#   "Y Percentage" ="#FF7F0EFF",
#   "Zr Percentage"="#FF7F0EFF",
#   "Mo Percentage"="#FF7F0EFF",
#   "Rh Percentage" ="#FF7F0EFF",
#   "Ag Percentage"="#FF7F0EFF",
#   "Cd Percentage"="#FF7F0EFF",
#   "Sn Percentage"="#FF7F0EFF",
#   "Sb Percentage"="#FF7F0EFF",
#   "W Percentage"="#FF7F0EFF",
#   "Re Percentage"="#FF7F0EFF",
#   "Hg Percentage"="#FF7F0EFF",
#   "Tl Percentage"="#FF7F0EFF",
#   "Pb Percentage"="#FF7F0EFF",
#   "Bi Percentage"="#FF7F0EFF"
# )

show_col(ggsci::pal_jama()(n=7))
# show_col(ggsci::pal_nejm()(n=7))
show_col(ggsci::pal_d3()(n=7))
# omics_color = c(
#   "Metabolite" = ggsci::pal_jama()(n=7)[1],
#   "Lipidome" = ggsci::pal_aaas()(n=7)[2],
#   "Stool microbiome" = ggsci::pal_jama()(n=7)[2],
#   "Skin microbiome" = ggsci::pal_jama()(n=7)[3],
#   "Oral microbiome" = ggsci::pal_jama()(n=7)[4],
#   "Nasal microbiome" = ggsci::pal_jama()(n=7)[5],
#   "Cytokine" = ggsci::pal_jama()(n=7)[6],
#   "Exposome" = ggsci::pal_jama()(n=7)[7],
#   "Proteome" = ggsci::pal_d3()(n=7)[2]
# )

require(wesanderson)
names(wes_palettes)

col = wes_palette("Zissou1", 28, type = "continuous")

show_col(col)

names(col)
# 
# custom_palette_2=c(
#   "Br Intensity"="#3A9AB2",
#   "Rh Intensity"="#4EA3B7",
#   "Hg Intensity"="#62ACBD",
#   "Rb Intensity"="#74B3BF",
#   "Sr Intensity"="#81B6BB",
#   "Y Intensity"="#8EB9B6",
#   "Zr Intensity"="#97BCB0",
#   "Sn Intensity"="#9EBFA8",
#   "Sb Intensity"="#A6C2A0",
#   "W Intensity"="#B0C493",
#   "Pb Intensity"="#B9C786",
#   "Bi Intensity"="#C4C875",
#   "Rb Percentage"="#D0C961",
#   "Y Percentage" ="#DCCB4E",
#   "Zr Percentage"="#DEC336",
#   "Mo Percentage"="#E1BB1E",
#   "Rh Percentage" ="#E3B20E",
#   "Ag Percentage"="#E5A60A",
#   "Cd Percentage"="#E69A05",
#   "Sn Percentage"="#E88E05",
#   "Sb Percentage"="#EA8305",
#   "W Percentage"="#EC7704",
#   "Re Percentage"="#ED6904",
#   "Hg Percentage"="#EE5C03",
#   "Tl Percentage"="#EF4902",
#   "Pb Percentage"="#F03201",
#   "Bi Percentage"="#F11B00"
# )






ggplot(pg2, aes( x = x,
                 next_x = next_x,
                 node = node,
                 next_node = next_node,
                 fill = factor(node),
                 label = node)) +
  geom_sankey(flow.alpha = 1,
              node.color = "gray30") +
  # geom_sankey_label(size = 3.5, color = "black", fill = "white",family = "Arial", face = "bold",position = "identity",nudge_x = 1,nudge_y = 1) +
  scale_fill_manual(values = c(custom_palette_2,omics_color,phenotype_class_color)) +
  # geom_sankey_text(size = 5, vjust = 0.5, hjust = -1,
  #                  color = "black", fontface = "bold") +
  # ggsci::scale_fill_aaas() +
  # theme_sankey(base_size = 18) +
  # labs(x = NULL) +
  # theme(legend.position = "none",
  #       plot.title = element_text(hjust = .5)) +
  # ggtitle("Car features")

  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(
    panel.grid = element_blank(), # 移除网格线
    axis.title = element_blank(), # 去掉轴标题
    axis.text.y = element_blank(), # 去掉y轴刻度
    axis.ticks.y = element_blank(), # 去掉y轴刻度线
    # legend.title = element_text(size = 14, family = "Arial", face = "bold"),
    legend.text = element_text(size = 14, family = "Arial", face = "bold"),
    axis.text.x = element_text(size = 14, family = "Arial", face = "bold", color = "black"), # 设置x轴文字大小和颜色
    legend.position = "bottom" # 移动图例到底部
    # legend.title = element_blank() # 去掉图例标题
    # plot.background = element_rect(fill = "#F5F5F5", color = NA), # 设置背景颜色
    # panel.background = element_rect(fill = "#F5F5F5", color = NA) # 面板背景色
  ) -> p2

p2

# dev.off()

# 
# ###添加竖直的参考线
# cairo_pdf( "result/liu_nian_202400606/mediation/sankey_picture/pool_cell_sperm_tidy_20241111_merge_lgm_p2_legend.pdf",width = 12, height = 8, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# # par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# p2
# 
# par(opar)
# dev.off()




















ggplot(pg2, aes( x = x,
                 next_x = next_x,
                 node = node,
                 next_node = next_node,
                 fill = factor(node),
                 label = node)) +
  geom_sankey(flow.alpha = 0.8,
              node.color = "gray30") +
  geom_sankey_label(size = 3.5, color = "black", fill = "white",family = "Arial", face = "bold",position = "identity",nudge_x = 1,nudge_y = 1) +
  # geom_sankey_text(size = 5, vjust = 0.5, hjust = -1,
  #                  color = "black", fontface = "bold") +
  scale_fill_manual(values = c(custom_palette_2,omics_color, phenotype_class_color)) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(
    panel.grid = element_blank(), # 移除网格线
    axis.title = element_blank(), # 去掉轴标题
    axis.text.y = element_blank(), # 去掉y轴刻度
    axis.ticks.y = element_blank(), # 去掉y轴刻度线
    # legend.title = element_text(size = 14, family = "Arial", face = "bold"),
    legend.text = element_text(size = 14, family = "Arial", face = "bold"),
    axis.text.x = element_text(size = 14, family = "Arial", face = "bold", color = "black"), # 设置x轴文字大小和颜色
    legend.position = " ", # 移动图例到底部
    legend.title = element_blank() # 去掉图例标题
    # plot.background = element_rect(fill = "#F5F5F5", color = NA), # 设置背景颜色
    # panel.background = element_rect(fill = "#F5F5F5", color = NA) # 面板背景色
  ) -> p3


p3



# 
# 
# ###添加竖直的参考线
# cairo_pdf( "result/liu_nian_202400606/mediation/sankey_picture/pool_cell_sperm_tidy_20241111_merge_lgm.pdf",width = 12, height = 8, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# # par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# p3
# 
# par(opar)
# dev.off()
# 
# 
# 
# 










df4 
 # write.xlsx(df4,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/sankey_links.xlsx")

 
 
 
 library(sankeywheel)
library(tidyr)

# 假设数据框是 data
links<- df4 %>%
  pivot_longer(cols = c(PFAS,Metabolite,STL), 
               names_to = "source", 
               values_to = "target") 

###计算权重值值
links %>%
  group_by(source, target) %>%
  summarise(weight = n()) -> links

# write.xlsx(links,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/sankey_links.xlsx")
# 查看转换后的数据
head(links)


nodes<- data.frame(name = c(as.character(links$source),
                            
                            as.character(links$target)) %>% unique)



links$IDsource <- match(links$source, nodes$name) -1

links$IDtarget <- match(links$target, nodes$name) -1

head(nodes)

head(links)

sankeyNetwork(Links = links,
                   
                  Nodes= nodes,
                   
                   Source= "IDsource",
                   
                   Target= "IDtarget",
                   
                   Value= "weight",
                   
                   NodeID= "name",
                   
                   sinksRight= FALSE,
                   
                   nodeWidth= 20,
                   
                   fontSize= 13,
                   
                   nodePadding= 15)








# 创造示例数据-----------------
## 第一层数据
data12 <-  df4 

layer12<- data12 %>%
  select(PFAS,Metabolite) %>% 
  dplyr::rename("source"="PFAS","target"="Metabolite") %>% 
  group_by(source, target) %>%
  summarise(weight = n())




## 第二层数据
# data23 <- data.frame(source=data12$target,
#                      target=df4$Metabolite[1:28])

layer23<- data12 %>%
  select(Metabolite,STL) %>% 
  dplyr::rename("source"="Metabolite","target"="STL") %>% 
  group_by(source, target) %>%
  summarise(weight = n())

# ## 第三层数据
# data34 <- data.frame(source=data23$target,
#                      target=sample(paste0('layer4_',1:4),100,prob = c(0.4,0.2,0.2,0.2),replace = T))
# layer34<- data34 %>%
#   group_by(source, target) %>%
#   summarise(weight = n())

## 合并三层数据
pdata <- rbind(layer12,layer23)
# pdata <- rbind(pdata,layer34)

# 保存示例数据
# write.csv(pdata,file = 'pdata.csv')



# 修改为符合networkD3包要求的数据格式-------------------
# 创建节点名称数据框
nodes <- data.frame(name = c(as.character(pdata$source), 
                             as.character(pdata$target)) %>% unique())
# 把source、target转换为数字
pdata$IDsource = match(pdata$source, nodes$name)-1
pdata$IDtarget = match(pdata$target, nodes$name)-1
head(pdata)


# 加载包
library(readxl)
library(tidyverse)
library(networkD3) #用于绘图
library(webshot) #用于图片格式转换
# 正式绘图---------------------------------------
p1 <- sankeyNetwork(Links = pdata, # 输入数据1
                    Nodes = nodes, # 输入数据2
                    Source = "IDsource", # 来源变量
                    Target = "IDtarget", # 接受变量
                    Value = "weight", # 关系权重
                    NodeID = "name", #节点名称
                    LinkGroup = 'source', # 颜色分组
                    sinksRight = FALSE, # 设置最后一层标签位置在左/右
                    nodeWidth = 15, #节点格子宽度
                    fontSize = 18, #文本标签字体的大小
                    fontFamily = "Arial",
                    width = 2480,
                    height = 1240,
                    nodePadding = 8) #节点格子间空隙宽度
p1

?sankeyNetwork


# webshot::install_phantomjs()
# 保存   
saveNetwork(p1,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey_20250326.html")
webshot("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey.html", "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey_20250326.png")
webshot("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey.html","result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey_20250326.pdf",zoom = 0.5)


?webshot


?sankeyNetwork

###############################################################################
###################筛选出PFAS单独有意义的进行桑葚图#################################
###############################################################################

####只筛出有意义的进行桑基图
pfas_meta_STL_p_kegg_meta_3

names(pfas_meta_STL_p_kegg_meta_3)

pfas_meta_STL_p_kegg_meta_3$p_value_estimate.x

pfas_meta_STL_p_kegg_meta_3 %>% 
  filter(p_value_estimate.x<0.2) ->pfas_meta_STL_p_kegg_meta_3_1 

pfas_meta_STL_p_kegg_meta_3_1

pfas_meta_STL_p_kegg_meta_3_1$var.names.x



# write.xlsx(pfas_meta_STL_p_kegg_meta_3_1,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/sankey.xlsx")



###将生成混合PFAS数据，结合到桑基图中

meta_STL_p_kegg_meta_2

names(meta_STL_p_kegg_meta_2)

str(meta_STL_p_kegg_meta_2)  ###根据行变量个数，查询行数

meta_STL_p_kegg_meta_2 %>%
  mutate(var.names.x=rep("PFAS_mixexposure",18)) -> pfas_meta_mix_2

names(pfas_meta_mix_2)
str(pfas_meta_mix_2)


pfas_meta_mix_2 %>% 
  select(var.names.x,name,group_STL) %>% 
  dplyr::rename("PFAS"="var.names.x","Metabolite"="name","STL"="group_STL") -> df1




pfas_meta_STL_p_kegg_meta_3_1  %>% 
  select(var.names.x,name,group_STL) %>% 
  dplyr::rename("PFAS"="var.names.x","Metabolite"="name","STL"="group_STL")-> df2



###
df1 %>% 
  rbind(df2) -> df3

df3$PFAS
df3 %>% 
  dplyr::mutate(PFAS=case_when(
    PFAS=="PFOS.beta"~"PFOS",
    PFAS=="PFOA.beta"~"PFOA",
    PFAS=="PFDA.beta"~"PFDA",
    PFAS=="L_PFHxS.beta"~"PFHxS",
    PFAS=="PFUdA.beta"~"PFUdA",
    PFAS=="PF3ONS_9CL.beta"~ "6:2 Cl-PFESA",
    PFAS=="PFAS_mixexposure"~"PFAS mixture"
  )) -> df4

df4



# # 
# # pfas_meta_mix <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/results_df_mix_p_value_mix.xlsx")
# # names(pfas_meta_mix)
# # str(pfas_meta_mix)
# # pfas_meta_mix$var.names
# 
# pfas_meta_mix %>%
#   mutate(var.names.x=rep("PFAS_mixexposure",414)) -> pfas_meta_mix_2
# 
# 
# meta_STL_p_kegg_meta_2 %>%
#   left_join(pfas_meta_mix_2,by=c("metabolite"="metabolite")) -> pfas_meta_STL_p_kegg_meta_4
# 
# names(pfas_meta_STL_p_kegg_meta_4)
# 
# pfas_meta_STL_p_kegg_meta_4$q_value.y
# 
# # 
# 

###>>>>>>>4.绘制桑基图<<< 

library(scales)
library(ggsci)
library(sankeywheel)

###
library(ggsankey)
library(dplyr)
library(ggplot2)

df4$PFAS
pg2 <- df4 %>%
  make_long("PFAS","Metabolite","STL")

phenotype_class_color = c(
  "Positive" = "#BC3C29FF",
  "Negative" = "#0072B5FF")


df4$PFAS

custom_palette_2=c( "PFAS mixture" = ggsci::pal_jama()(n=7)[4],
                    "PFOS" = ggsci::pal_jama()(n=7)[1],
                    "PFOA" = ggsci::pal_jama()(n=7)[2],
                    "PFDA" = ggsci::pal_jama()(n=7)[3],
                    "PFUdA" = ggsci::pal_jama()(n=7)[5],
                    "PFHxS" = ggsci::pal_jama()(n=7)[6],
                    "6:2 Cl-PFESA" = ggsci::pal_jama()(n=7)[7]
)


show_col(ggsci::pal_jama()(n=7))
show_col(ggsci::pal_nejm()(n=7))

show_col(ggsci::pal_d3()(n=7))


# df4$Metabolite[1:33] %in% c(viridis(33))[2]

library(viridis)
df4$Metabolite

omics_color = c(
  "Bisnorcholic acid" = c(viridis(25))[1],
  "tetracosatetraenoate" = c(viridis(25))[2],
  "Tetracosapentaenoic acid (24:5n-3)" = c(viridis(25))[3],
  
  "12-Oxo-20-carboxy-leukotriene B4" = c(viridis(25))[4],
  "(3E,5Z,11Z)-Pentadeca-3,5,11-trienedioylcarnitine" = c(viridis(25))[5],
  "9-Hydroxy-16-oxooctadeca-10,12,14-trienoylcarnitine" = c(viridis(25))[6],
  
  "3-Hydroxyhex-4-enoylcarnitine" = c(viridis(25))[7],
  "13,14-dihydro-15-keto-tetranor Prostaglandin F1beta" = c(viridis(25))[8],
  "Testosterone" = c(viridis(25))[9],
  
  "Linoleic acid" = c(viridis(25))[10],
  "L-Acetylcarnitine" = c(viridis(25))[11],
  "O-Propanoylcarnitine" = c(viridis(25))[12],
  
  "4-Hydroxy-4-methyl-2-oxoadipate" = c(viridis(25))[13],
  "Prostaglandin E1" = c(viridis(25))[14],
  "2-Hydroxyestradiol" = c(viridis(25))[15],
  
  "Tetrahydrocortisol" = c(viridis(25))[16],
  "Eicosapentaenoic acid" = c(viridis(25))[17],
  "3-Hydroxyhexadecanoylcarnitine" = c(viridis(25))[18]
  # 
  # "2-Methoxy-17beta-estradiol" = c(viridis(25))[19],
  # "Tetrahydrocortisol" = c(viridis(25))[20],
  # "Eicosapentaenoic acid" = c(viridis(25))[21],
  # 
  # "Docosahexaenoic acid" = c(viridis(25))[22],
  # "3-Hydroxyhexadecanoylcarnitine" = c(viridis(25))[23],
  # "19-Noraldosterone" = c(viridis(25))[24],
  # 
  # "10-Hydroperoxydocosahexaenoic acid" = c(viridis(25))[25]
  
)



# "L-Acetylcarnitine" = c(viridis(28))[1],omics_color= c(wes_palette("Zissou1", 28, type = "continuous"))

# 获取Zissou1色带，连续的颜色
colors <- c(wes_palette("Zissou1", 18, type = "continuous"))

# 给代谢物变量赋予颜色

viridis(33) ## viridis 主题中提取10个颜色
inferno(28) ## inferno 主题中提取10个颜色


plot(1:28, col=viridis(10),cex=4,pch=20)

barplot(rep(28,28),col=viridis(28),space=0,border = viridis(28))




# "Metabolite" = ggsci::pal_jama()(n=7)[1],
# "Lipidome" = ggsci::pal_aaas()(n=7)[2],

# df5$`Exogenous metals`
# 
# custom_palette_2=c(
#   "Br Intensity"="#00A1D5FF",
#   "Rh Intensity"="#00A1D5FF",
#   "Hg Intensity"="#00A1D5FF",
#   "Rb Intensity"="#00A1D5FF",
#   "Sr Intensity"="#00A1D5FF",
#   "Y Intensity"="#00A1D5FF",
#   "Zr Intensity"="#00A1D5FF",
#   "Sn Intensity"="#00A1D5FF",
#   "Sb Intensity"="#00A1D5FF",
#   "W Intensity"="#00A1D5FF",
#   "Pb Intensity"="#00A1D5FF",
#   "Bi Intensity"="#00A1D5FF",
#   "Rb Percentage"="#FF7F0EFF",
#   "Y Percentage" ="#FF7F0EFF",
#   "Zr Percentage"="#FF7F0EFF",
#   "Mo Percentage"="#FF7F0EFF",
#   "Rh Percentage" ="#FF7F0EFF",
#   "Ag Percentage"="#FF7F0EFF",
#   "Cd Percentage"="#FF7F0EFF",
#   "Sn Percentage"="#FF7F0EFF",
#   "Sb Percentage"="#FF7F0EFF",
#   "W Percentage"="#FF7F0EFF",
#   "Re Percentage"="#FF7F0EFF",
#   "Hg Percentage"="#FF7F0EFF",
#   "Tl Percentage"="#FF7F0EFF",
#   "Pb Percentage"="#FF7F0EFF",
#   "Bi Percentage"="#FF7F0EFF"
# )

show_col(ggsci::pal_jama()(n=7))
# show_col(ggsci::pal_nejm()(n=7))
show_col(ggsci::pal_d3()(n=7))
# omics_color = c(
#   "Metabolite" = ggsci::pal_jama()(n=7)[1],
#   "Lipidome" = ggsci::pal_aaas()(n=7)[2],
#   "Stool microbiome" = ggsci::pal_jama()(n=7)[2],
#   "Skin microbiome" = ggsci::pal_jama()(n=7)[3],
#   "Oral microbiome" = ggsci::pal_jama()(n=7)[4],
#   "Nasal microbiome" = ggsci::pal_jama()(n=7)[5],
#   "Cytokine" = ggsci::pal_jama()(n=7)[6],
#   "Exposome" = ggsci::pal_jama()(n=7)[7],
#   "Proteome" = ggsci::pal_d3()(n=7)[2]
# )

require(wesanderson)
names(wes_palettes)

col = wes_palette("Zissou1", 28, type = "continuous")

show_col(col)

names(col)
# 
# custom_palette_2=c(
#   "Br Intensity"="#3A9AB2",
#   "Rh Intensity"="#4EA3B7",
#   "Hg Intensity"="#62ACBD",
#   "Rb Intensity"="#74B3BF",
#   "Sr Intensity"="#81B6BB",
#   "Y Intensity"="#8EB9B6",
#   "Zr Intensity"="#97BCB0",
#   "Sn Intensity"="#9EBFA8",
#   "Sb Intensity"="#A6C2A0",
#   "W Intensity"="#B0C493",
#   "Pb Intensity"="#B9C786",
#   "Bi Intensity"="#C4C875",
#   "Rb Percentage"="#D0C961",
#   "Y Percentage" ="#DCCB4E",
#   "Zr Percentage"="#DEC336",
#   "Mo Percentage"="#E1BB1E",
#   "Rh Percentage" ="#E3B20E",
#   "Ag Percentage"="#E5A60A",
#   "Cd Percentage"="#E69A05",
#   "Sn Percentage"="#E88E05",
#   "Sb Percentage"="#EA8305",
#   "W Percentage"="#EC7704",
#   "Re Percentage"="#ED6904",
#   "Hg Percentage"="#EE5C03",
#   "Tl Percentage"="#EF4902",
#   "Pb Percentage"="#F03201",
#   "Bi Percentage"="#F11B00"
# )






ggplot(pg2, aes( x = x,
                 next_x = next_x,
                 node = node,
                 next_node = next_node,
                 fill = factor(node),
                 label = node)) +
  geom_sankey(flow.alpha = 1,
              node.color = "gray30") +
  # geom_sankey_label(size = 3.5, color = "black", fill = "white",family = "Arial", face = "bold",position = "identity",nudge_x = 1,nudge_y = 1) +
  scale_fill_manual(values = c(custom_palette_2,omics_color,phenotype_class_color)) +
  # geom_sankey_text(size = 5, vjust = 0.5, hjust = -1,
  #                  color = "black", fontface = "bold") +
  # ggsci::scale_fill_aaas() +
  # theme_sankey(base_size = 18) +
  # labs(x = NULL) +
  # theme(legend.position = "none",
  #       plot.title = element_text(hjust = .5)) +
  # ggtitle("Car features")
  
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(
    panel.grid = element_blank(), # 移除网格线
    axis.title = element_blank(), # 去掉轴标题
    axis.text.y = element_blank(), # 去掉y轴刻度
    axis.ticks.y = element_blank(), # 去掉y轴刻度线
    # legend.title = element_text(size = 14, family = "Arial", face = "bold"),
    legend.text = element_text(size = 14, family = "Arial", face = "bold"),
    axis.text.x = element_text(size = 14, family = "Arial", face = "bold", color = "black"), # 设置x轴文字大小和颜色
    legend.position = "bottom" # 移动图例到底部
    # legend.title = element_blank() # 去掉图例标题
    # plot.background = element_rect(fill = "#F5F5F5", color = NA), # 设置背景颜色
    # panel.background = element_rect(fill = "#F5F5F5", color = NA) # 面板背景色
  ) -> p2

p2

# dev.off()

# 
# ###添加竖直的参考线
# cairo_pdf( "result/liu_nian_202400606/mediation/sankey_picture/pool_cell_sperm_tidy_20241111_merge_lgm_p2_legend.pdf",width = 12, height = 8, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# # par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# p2
# 
# par(opar)
# dev.off()




















ggplot(pg2, aes( x = x,
                 next_x = next_x,
                 node = node,
                 next_node = next_node,
                 fill = factor(node),
                 label = node)) +
  geom_sankey(flow.alpha = 0.8,
              node.color = "gray30") +
  geom_sankey_label(size = 3.5, color = "black", fill = "white",family = "Arial", face = "bold",position = "identity",nudge_x = 1,nudge_y = 1) +
  # geom_sankey_text(size = 5, vjust = 0.5, hjust = -1,
  #                  color = "black", fontface = "bold") +
  scale_fill_manual(values = c(custom_palette_2,omics_color, phenotype_class_color)) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(
    panel.grid = element_blank(), # 移除网格线
    axis.title = element_blank(), # 去掉轴标题
    axis.text.y = element_blank(), # 去掉y轴刻度
    axis.ticks.y = element_blank(), # 去掉y轴刻度线
    # legend.title = element_text(size = 14, family = "Arial", face = "bold"),
    legend.text = element_text(size = 14, family = "Arial", face = "bold"),
    axis.text.x = element_text(size = 14, family = "Arial", face = "bold", color = "black"), # 设置x轴文字大小和颜色
    legend.position = " ", # 移动图例到底部
    legend.title = element_blank() # 去掉图例标题
    # plot.background = element_rect(fill = "#F5F5F5", color = NA), # 设置背景颜色
    # panel.background = element_rect(fill = "#F5F5F5", color = NA) # 面板背景色
  ) -> p3


p3



# 
# 
# ###添加竖直的参考线
# cairo_pdf( "result/liu_nian_202400606/mediation/sankey_picture/pool_cell_sperm_tidy_20241111_merge_lgm.pdf",width = 12, height = 8, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# # par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# p3
# 
# par(opar)
# dev.off()
# 
# 
# 
# 










df4 
# write.xlsx(df4,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/sankey_links.xlsx")




library(sankeywheel)
library(tidyr)

# 假设数据框是 data
links  <- df4 %>%
  pivot_longer(cols = c(PFAS,Metabolite,STL), 
               names_to = "source", 
               values_to = "target") 

###计算权重值值
# links %>%
#   group_by(source, target) %>%
#    summarise(weight = n()) -> links

links %>%
  dplyr::group_by(source, target) %>%
  dplyr::summarise(weight = n()) -> links

# 
# links %>%
#   dplyr::count(source, target) -> links



# write.xlsx(links,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/sankey_links.xlsx")
# 查看转换后的数据
head(links)


nodes<- data.frame(name = c(as.character(links$source),
                            
                            as.character(links$target)) %>% unique)



links$IDsource <- match(links$source, nodes$name) -1

links$IDtarget <- match(links$target, nodes$name) -1

head(nodes)

head(links)

sankeyNetwork(Links = links,
              
              Nodes= nodes,
              
              Source= "IDsource",
              
              Target= "IDtarget",
              
              Value= "weight",
              
              NodeID= "name",
              
              sinksRight= FALSE,
              
              nodeWidth= 20,
              
              fontSize= 13,
              
              nodePadding= 15)








# 创造示例数据-----------------
## 第一层数据
data12 <-  df4 

layer12<- data12 %>%
  select(PFAS,Metabolite) %>% 
  dplyr::rename("source"="PFAS","target"="Metabolite") %>% 
  group_by(source, target) %>%
  dplyr::summarise(weight = n())




## 第二层数据
# data23 <- data.frame(source=data12$target,
#                      target=df4$Metabolite[1:28])

layer23<- data12 %>%
  select(Metabolite,STL) %>% 
  dplyr::rename("source"="Metabolite","target"="STL") %>% 
  group_by(source, target) %>%
  dplyr::summarise(weight = n())

# ## 第三层数据
# data34 <- data.frame(source=data23$target,
#                      target=sample(paste0('layer4_',1:4),100,prob = c(0.4,0.2,0.2,0.2),replace = T))
# layer34<- data34 %>%
#   group_by(source, target) %>%
#   summarise(weight = n())

## 合并三层数据
pdata <- rbind(layer12,layer23)
# pdata <- rbind(pdata,layer34)

# 保存示例数据
# write.csv(pdata,file = 'pdata.csv')



# 修改为符合networkD3包要求的数据格式-------------------
# 创建节点名称数据框
nodes <- data.frame(name = c(as.character(pdata$source), 
                             as.character(pdata$target)) %>% unique())
# 把source、target转换为数字
pdata$IDsource = match(pdata$source, nodes$name)-1
pdata$IDtarget = match(pdata$target, nodes$name)-1
head(pdata)


# 加载包
library(readxl)
library(tidyverse)
library(networkD3) #用于绘图
library(webshot) #用于图片格式转换
# 正式绘图---------------------------------------
p1 <- sankeyNetwork(Links = pdata, # 输入数据1
                    Nodes = nodes, # 输入数据2
                    Source = "IDsource", # 来源变量
                    Target = "IDtarget", # 接受变量
                    Value = "weight", # 关系权重
                    NodeID = "name", #节点名称
                    LinkGroup = 'source', # 颜色分组
                    sinksRight = FALSE, # 设置最后一层标签位置在左/右
                    nodeWidth = 15, #节点格子宽度
                    fontSize = 18, #文本标签字体的大小
                    fontFamily = "Arial",
                    width = 2480,
                    height = 1240,
                    nodePadding = 8) #节点格子间空隙宽度
p1

?sankeyNetwork


# webshot::install_phantomjs()
# 保存   
saveNetwork(p1,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey_revser_20250326.html")

webshot("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey_revser_20250326.html", "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey_revser_20250326.png")

webshot("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey_revser_20250326.html","result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/mediation/sankey_revser_20250326.pdf",zoom = 0.5)


?webshot


?sankeyNetwork














