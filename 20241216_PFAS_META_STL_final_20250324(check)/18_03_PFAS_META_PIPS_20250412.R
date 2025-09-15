# no_source() 

no_source() 
rm(list = ls())
# setwd(r4projects::get_project_wd())
source("R/100-tools.R")

library(tidyverse)
library(tidymass)



pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts,scales,ggsci)
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



str(pfas_meta_stl_all)

names(pfas_meta_stl_all)





pfas_meta_stl_all  %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg     ####取交集之后的最终筛出的代谢物质

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_meeting_q.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))


###########读入每个PFAS化合物的对应代谢物质的pip数值

df_pip <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_final_mix_line.xlsx")

df_pip %>% 
  left_join(pfas_meta_STL_p_kegg,by="metabolite") %>% 
  filter(!BCI_CI95.y=="NA")    -> df_pip_meta


str(df_pip_meta)

df_pip_meta %>% 
  select(name,PFOS.gamma:PF3ONS_9CL.gamma) -> data

# data %>% 
#   as.matrix() -> data


data %>% 
  column_to_rownames(var = "name")    -> data    ####强制设定行变量名称

# 
# # 将 name 列设置为行名
# df_1 <-as.character(data$name)
# 
# rownames(data) <- df_1
# 
# data <- data[, -1]
# 
# data
# 
rownames(data)



df_pip_meta %>% 
  select(name,PFOS.gamma:PF3ONS_9CL.gamma) %>% 
  filter(PFOS.gamma<=0.167)




6/18




df_pip_meta %>% 
  select(name,PFOS.gamma:PF3ONS_9CL.gamma) %>% 
  filter(PFOA.gamma<=0.167)




7/18





df_pip_meta %>% 
  select(name,PFOS.gamma:PF3ONS_9CL.gamma) %>% 
  filter(PFDA.gamma<=0.167)




4/18





df_pip_meta %>% 
  select(name,PFOS.gamma:PF3ONS_9CL.gamma) %>% 
  filter(PFUdA.gamma<=0.167)




15/18




df_pip_meta %>% 
  select(name,PFOS.gamma:PF3ONS_9CL.gamma) %>% 
  filter(L_PFHxS.gamma<=0.167)




6/18






df_pip_meta %>% 
  select(name,PFOS.gamma:PF3ONS_9CL.gamma) %>% 
  filter(PF3ONS_9CL.gamma<=0.167)




14/18













# 绘制热图
# 
# # 绘制热图，并显示每个 gamma 数值
# pheatmap(data,
#          scale = "none",  # 不进行标准化
#          color = colorRampPalette(c("white", "red"))(100),  # 颜色梯度
#          clustering_distance_rows = "euclidean",  # 行聚类距离
#          clustering_distance_cols = "euclidean",  # 列聚类距离
#          display_numbers = TRUE,  # 显示数值
#          number_format = "%.3f",  # 数值格式，保留 3 位小数
#          main = "Gamma Values of Metabolites Across Different PFAS",  # 图表标题
#          fontsize_row = 8,  # 行字体大小
#          fontsize_col = 8)  # 列字体大小
# 
# 
# # 绘制热图，并显示每个 gamma 数值
# pheatmap(data,
#          scale = "none",  # 不进行标准化
#          color = colorRampPalette(c("white", "red"))(100),  # 颜色梯度
#          clustering_distance_rows = "euclidean",  # 行聚类距离
#          clustering_distance_cols = "euclidean",  # 列聚类距离
#          display_numbers = TRUE,  # 显示数值
#          number_format = "%.3f",  # 数值格式，保留 3 位小数
#          main = "Gamma Values of Metabolites Across Different PFAS",  # 图表标题
#          fontsize_row = 8,  # 行字体大小
#          fontsize_col = 8,  # 列字体大小
#          show_rownames = TRUE,  # 显示行名
#          show_colnames = TRUE)  # 显示列名
# 
# 
# 
# library(grid)
# 
# 
# pheatmap(data,
#          scale = "none",  # 不进行标准化
#          color = colorRampPalette(c("white", "red"))(100),  # 颜色梯度
#          clustering_distance_rows = "euclidean",  # 行聚类距离
#          clustering_distance_cols = "euclidean",  # 列聚类距离
#          cluster_rows = FALSE,  # 禁用行聚类
#          cluster_cols = FALSE,  # 禁用列聚类
#          display_numbers = TRUE,  # 显示数值
#          number_format = "%.3f",  # 数值格式，保留 3 位小数
#          main = "Gamma Values of Metabolites Across Different PFAS",  # 图表标题
#          fontsize_row = 8,  # 行字体大小
#          fontsize_col = 8,  # 列字体大小
#          fontface = "bold",
#          fontface_row="bold",
#          fontface_col="bold",
#          # fontsize=20,
#          annotation_colors="black",
#          fontsize_number = 12,   ###热图中字体加粗
#          fontfamily = "Arial",  # 使用支持加粗的字体族
#          show_rownames = TRUE,  # 显示行名
#          show_colnames = TRUE)  # 显示列名
# 
# # 使用 grid 包调整图例位置
# grid.ls()  # 查看所有图形元素的名称
# grid.edit("legend", vp = viewport(x = unit(0.9, "npc"), y = unit(0.5, "npc")))  # 将图例移动到右侧中间
# 
# ?pheatmap
# 
# 
# library(ComplexHeatmap)
# 
# 
# 
# # 创建热图
# Heatmap(data,
#         col = colorRampPalette(c("white", "red"))(100),  # 颜色梯度
#         row_names_gp = gpar(fontsize = 8, fontface = "bold"),  # 行名样式
#         column_names_gp = gpar(fontsize = 8, fontface = "bold"),  # 列名样式
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.3f", data[i, j]), x, y, gp = gpar(fontsize = 12, fontface = "bold"))
#         },  # 显示数值
#         heatmap_legend_param = list(
#           title = "Gamma Values",  # 图例标题
#           legend_height = unit(4, "cm"),  # 图例高度
#           legend_width = unit(1.5, "cm"),  # 图例宽度
#           title_position = "leftcenter-rot",  # 图例标题位置
#           legend_direction = "vertical"  # 图例方向
#         ),
#         name = "heatmap")  # 热图名称
# 
# 
# 



# 创建热图，禁用聚类线并按列变量顺序显示
plot <- Heatmap(data,
        col = colorRampPalette(c("white", "red"))(100),  # 颜色梯度
        row_names_gp = gpar(fontsize = 12, fontface = "bold",fontfamily="Arial"),  # 行名样式
        column_names_gp = gpar(fontsize = 12, fontface = "bold",fontfamily="Arial"),  # 列名样式
        cluster_rows = FALSE,  # 禁用行聚类
        cluster_columns = FALSE,  # 禁用列聚类
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", data[i, j]), x, y, gp = gpar(fontsize = 12, fontface = "bold",fontfamily="Arial"))
        },  # 显示数值
        heatmap_legend_param = list(
          title = "Gamma Values",  # 图例标题
          legend_height = unit(4, "cm"),  # 图例高度
          legend_width = unit(1.5, "cm"),  # 图例宽度
          title_position = "leftcenter-rot",  # 图例标题位置
          legend_direction = "vertical"  # 图例方向
        ),
        name = "heatmap")  # 热图名称





font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/pfas_meta_pip_20250411.pdf",width = 11, height = 9,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

plot  # 使用大写罗马表示

par(opar)
dev.off()












###############缩放瓦片图的长宽高

# 加载 ComplexHeatmap 包
library(ComplexHeatmap)

# 创建热图，禁用聚类线并按列变量顺序显示，缩小方块大小
plot1 <- Heatmap(data,
        col = colorRampPalette(c("white", "red"))(100),  # 颜色梯度
        row_names_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Arial"),  # 行名样式
        column_names_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Arial"),  # 列名样式
        cluster_rows = FALSE,  # 禁用行聚类
        cluster_columns = FALSE,  # 禁用列聚类
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", data[i, j]), x, y, 
                    gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Arial"))
        },  # 显示数值
        width = unit(120, "mm"),  # 设置列宽为 5 毫米
        height = unit(120, "mm"),  # 设置行高为 5 毫米
        heatmap_legend_param = list(
          title = "Gamma Values",  # 图例标题
          legend_height = unit(4, "cm"),  # 图例高度
          legend_width = unit(1.5, "cm"),  # 图例宽度
          title_position = "leftcenter-rot",  # 图例标题位置
          legend_direction = "vertical"  # 图例方向
        ),
        name = "heatmap")  # 热图名称






font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/pfas_meta_pip_20250411_1.pdf",width = 11, height = 9,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

plot1  # 使用大写罗马表示

par(opar)
dev.off()









