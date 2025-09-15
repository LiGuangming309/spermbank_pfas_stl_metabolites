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
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts)
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


####>>>>>>>>读入混合物计算后的P值------------------------------
# write.xlsx(results_df_mix_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/results_df_mix_p_value.xlsx")
pfas_meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q_PFAS.xlsx")
names(pfas_meta_STL)

pfas_meta_STL$meta

pfas_meta_STL%>% 
  dplyr::rename(metabolite="meta") ->pfas_meta_STL_1  


###过滤出PFAS-META-STL<0.2的Q值
pfas_meta_STL_1 %>% 
  filter(q_value<0.2) -> pfas_meta_STL_2



pfas_meta_STL_2 %>% 
  left_join(df_unique_2,by="metabolite") -> pfas_meta_STL_p_kegg

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_q.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))

# 
# df <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy.xlsx")
# 
# 
# 
# library("magrittr")
# 
# 
# #######所有化合物的火山图
# 
# 
# marker_color <-
#   c(
#     "Up" = ggsci::pal_futurama()(n = 9)[6],
#     "Down" = ggsci::pal_futurama()(n = 9)[3],
#     "No" = ggsci::pal_futurama()(n = 9)[9]
#   )
# 
# show_col(ggsci::pal_futurama()(n = 9))
# 
# 
# pfas_meta_STL_p_kegg
# 
# 
# names(pfas_meta_STL_p_kegg)
# 
# library(showtext)
# library(sysfonts)
# 
# font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
# 
# pfas_meta_STL_1 %>% 
#   left_join(df_unique_2,by="metabolite") -> pfas_meta_STL_p_kegg_2
# 
# 
# volcano_plot_all <-
#   pfas_meta_STL_p_kegg_2 %>% 
#   mutate(marker=case_when(
#     q_value < 0.2 & estimate > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
#     q_value< 0.2 & estimate< 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
#     TRUE ~ "No"
#   )) %>%
#   ggplot(aes(x =estimate, y = -log10(q_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
#   # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
#   geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
#   scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
#   geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
#   geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
#   labs(x = "Differential metabolites", y = "-log10(q value)")+
#   theme_bw() +
#   theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         # legend.position = c(0.25,0.15),
#         legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#         # axis.line = element_line(linewidth = 0.9),
#         # axis.ticks = element_line(linewidth = 0.9),
#         axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
#         axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
#         panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
#         panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
#         panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
#   )+
#   geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) 
# 
# # ggplot(aes(x =estimate, y = -log10(q_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
# # # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
# # geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
# # # ylim(c(0,10))+
# # # xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围
# # scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
# # geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
# # # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线
# # geom_hline(yintercept = 0, linetype = 2,size=1) +
# # # facet_wrap(~var.names, scales = "free") +
# # labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
# # # ggtitle("单组火山图") + #标题
# # labs(x = "Estimate", y = "-log10(q_value)")+
# # theme_bw() +
# # theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
# #       text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
# #       # legend.position = c(0.25,0.15),
# #       legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
# #       axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
# #       # axis.line = element_line(linewidth = 0.9),
# #       # axis.ticks = element_line(linewidth = 0.9),
# #       axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
# #       axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
# #       panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
# #       panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
# #       panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
# # )+
# # geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) 
# 
# 
# volcano_plot_all
# 
# # dev.off()
# 
# # 
# # 
# # 
# font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
# cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/volcano_pfas_meta_STL_merger_line.pdf",width = 6, height = 5, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# volcano_plot_all #出图
# 
# par(opar)
# dev.off()
# 







#######2.绘制热图
###读入原始的数据

library(ComplexHeatmap)
df_rawdata <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy_623.xlsx")


str(df_rawdata)
names(df_rawdata)



# 加载必要的R包
library(Mfuzz)  # 用于模糊C均值聚类
library(dplyr)  # 用于数据操作

# 1. 数据预处理
# # 假设 omics_data 是一个矩阵或数据框，行为样本，列为分子（如基因或代谢物）
# omics_data <- read.csv("omics_data.csv", row.names = 1)  # 读取数据
# 
# # log2转换
# omics_data_log2 <- log2(omics_data + 1)  # 加1以避免log2(0)的情况
# 
# # 自动标准化（均值为0，标准差为1）

omics_data_scaled <- df_rawdata[, c(1,25:381)] %>% dplyr::rename(time="number") %>% t()

str(omics_data_scaled)

colnames(omics_data_scaled) <- omics_data_scaled[1,]  # 将第一列设置为行名

omics_data_scaled


write.table(
  omics_data_scaled,
  file = "temp_data.txt",
  sep = '\t',
  quote = FALSE,
  col.names = NA
)

#read it back in as an expression set
data <- table2eset(filename = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/cluster/temp_data.txt")
data.s <- standardise(data)

m1 <- mestimate(data.s)
m1


plot <-
  Dmin(
    data.s,
    m = m1,
    crange = seq(2, 40, 1),
    repeats = 3,
    visu = TRUE
  )

plot <-
  plot %>%
  data.frame(distance = plot,
             k = seq(2, 40, 1)) %>%
  ggplot(aes(k, distance)) +
  geom_point(shape = 21, size = 4, fill = "black") +
  # geom_smooth() +
  geom_segment(aes(
    x = k,
    y = 0,
    xend = k,
    yend = distance
  )) +
  theme_bw() +
  theme(
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(x = "Cluster number",
       y = "Min. centroid distance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

plot

ggsave(plot,
       filename = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/cluster/distance_k_number.pdf",
       width = 7,
       height = 7)





cluster_number <- 5

data.s

c <- mfuzz(data.s, c = cluster_number, m = m1)


# ####any two clusters with correlation > 0.8 should be considered as one
library(corrplot)
layout(1)
# center <- c$centers

membership_cutoff <- 0.5

center <-
  get_mfuzz_center(data = data.s,
                   c = c,
                   membership_cutoff = 0.5)

rownames(center) <- paste("Cluster", rownames(center), sep = ' ')


corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  # addrect = 5,
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
) 

mfuzz.plot(
  eset = data.s,
  min.mem = 0.5,
  cl = c,
  mfrow = c(4, 4),
  # time.labels = time,
  new.window = FALSE
)


# 自动计算绘图布局
nrow <- ceiling(sqrt(cluster_number))
ncol <- ceiling(cluster_number/nrow)

mfuzz.plot(
  eset = data.s,        # 标准化后的表达数据
  min.mem = 0.5,        # 最小隶属度阈值
  cl = c,               # 聚类结果
  mfrow = c(nrow, ncol), # 根据聚类数自动调整布局
  new.window = FALSE    # 不在新窗口中绘图
)

idx <- 1


library(ComplexHeatmap)

###
cluster_color <-
  ggsci::pal_jama()(n = 7)[1:cluster_number]

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

####plot for each cluster
idx <- 1

temp_data <-
  data.s %>% as.data.frame() %>%
  t() %>% as.data.frame()

for (idx in 1:cluster_number) {
  cat(idx, " ")
  
  cluster_data <-
    cluster_info %>%
    # dplyr::filter(cluster == idx) %>%
    dplyr::select(1, 1 + idx, cluster)
  
  colnames(cluster_data)[2] <- c("membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > membership_cutoff)
  # 定义目标路径
  base_path <- "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/cluster"
  
  # 创建子目录并写入 Excel 文件
  path <- file.path(base_path, paste("cluster", idx, sep = "_"))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  
  openxlsx::write.xlsx(
    cluster_data,
    file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  temp_center <-
    center[idx, , drop = TRUE] %>%
    unlist() %>%
    data.frame(time = names(.),
               value = .,
               stringsAsFactors = FALSE) %>%
    dplyr::mutate(time = factor(time, levels = time)) %>%
    dplyr::mutate(time_point = time)
  
  temp_center$time <-
    temp_center$time %>%
    # %>%
    as.numeric()
    # as.character() %>%
    # stringr::str_split(pattern = "_") %>%
    # lapply(function(x) {
    #   mean(as.numeric(x))
    # }) %>%
    # unlist()
  
  temp <-
    temp_data[cluster_data$variable_id, ] %>%
    data.frame(
      membership = cluster_data$membership,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id, membership),
      names_to = "time",
      values_to = "value"
    ) %>%
    dplyr::mutate(time = factor(time, levels = unique(time))) %>%
    dplyr::mutate(time_point = time)
  
  temp$time <-
    temp$time %>%
    as.numeric()
    # as.character() %>%
    # stringr::str_split(pattern = "_") %>%
    # lapply(function(x) {
    #   mean(as.numeric(x))
    # }) %>%
    # unlist()
  
  plot <-
    temp %>%
    dplyr::arrange(membership, variable_id) %>%
    dplyr::mutate(variable_id = factor(variable_id, levels = unique(variable_id))) %>%
    ggplot(aes(time, value, group = variable_id)) +
    geom_line(aes(color = membership), alpha = 0.7) +
    theme_bw() +
    theme(
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 12
      ),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "",
      y = "Z-score",
      title = paste(
        "Cluster ",
        idx,
        " (",
        nrow(cluster_data),
        " metabolic features)",
        sep = ""
      )
    ) +
    geom_line(
      mapping = aes(time, value, group = 1),
      data = temp_center,
      size = 2
    ) +
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis() +
    scale_x_continuous(breaks = c(temp_center$time),
                       labels = temp_center$time_point)
  
  plot
  
  ggsave(
    plot,
    filename = file.path(path, paste("cluster", idx, ".pdf", sep = "")),
    width = 8,
    height = 7
  )
  
}




####导出所有的分类变量
dim(cluster_data)

table(cluster_info$cluster)

cluster_info <-
  unique(cluster_info$cluster) %>%
  purrr::map(function(x) {
    temp <-
      cluster_info %>%
      # dplyr::filter(cluster == x) %>%
      dplyr::select(variable_id, paste0("X", x), cluster)
    colnames(temp)[2] <- "membership"
    temp <-
      temp %>%
      dplyr::filter(membership >= membership_cutoff)
    temp <-
      temp %>%
      dplyr::mutate(cluster_raw = cluster) %>%
      dplyr::mutate(cluster = x)
    temp
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

cluster_info %>%
  dplyr::count(cluster)

cluster_info %>%
  dplyr::filter(membership > 0.5) %>%
  dplyr::count(cluster)


####分类

final_cluster_info <-
  cluster_info

save(final_cluster_info, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/cluster/final_cluster_info")

openxlsx::write.xlsx(
  final_cluster_info,
  file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/cluster/final_cluster_info.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)


####将代谢通路数据与模糊分类数据进行合并
pfas_meta_STL_p_kegg

names(pfas_meta_STL_p_kegg)
names(final_cluster_info)




pfas_meta_STL_p_kegg %>% 
  left_join(final_cluster_info,by=c("metabolite"="variable_id")) -> pfas_meta_STL_p_kegg_clustr


write.xlsx(pfas_meta_STL_p_kegg_clustr,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/cluster/pfas_meta_STL_p_kegg_clustr.xlsx")






names(df_meta_final) %>% as.data.frame() %>%
  dplyr::rename(value=".") %>%
  slice(-1)

final_cluster_info %>%
  filter(variable_id==df_meta_final[,-1])





###筛选出q<0.2的结果

df_rawdata %>% 
  select(number:STL,pfas_meta_STL_p_kegg$metabolite) -> df_rawdata_q

df_rawdata_q

names(df_rawdata_q)



df_rawdata_q %>% 
  select(1,M642:M609) -> df_meta_final




df_meta_final %>% 
  t() -> df_meta_final_4


###将第一行作为列名
colnames(df_meta_final_4) <- df_meta_final_4[1,]

df_meta_final_4[-1,] -> df_meta_final_5


###对行重新生成行名
rownames(df_meta_final_5) <- pfas_meta_STL_p_kegg$name


df_meta_final_5


# write.csv(df_meta_final_5,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/result_pfas_meta_p_value.csv")



Heatmap(df_meta_final_5, cluster_columns =TRUE, cluster_rows = TRUE)

library("pheatmap")

# library()
# 


library(ComplexHeatmap)

ComplexHeatmap::Heatmap(df_meta_final_5,
                        cluster_columns = TRUE,     
                        cluster_rows = TRUE,        
                        show_row_names = TRUE,     
                        show_column_names = FALSE,  
                        col = colorRampPalette(c("blue", "white", "red"))(100),
                        row_names_gp = gpar(fontsize = 10, fontface = "bold",family="Arial"),  # 行标签字体样式
                        column_names_gp = gpar(fontsize = 10, fontface = "bold",family="Arial"),  # 列标签字体样式
                        clustering_distance_rows = "euclidean",  
                        clustering_distance_columns = "euclidean",  
                        clustering_method_rows = "complete",  
                        clustering_method_columns = "complete"
) -> p2

p2


library(pheatmap)

?pheatmap

pheatmap(df_meta_final_5, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         cutree_rows = 5, 
         cutree_cols = NA, 
         col = colorRampPalette(c("blue", "white", "red"))(100), 
         row_names_gp = gpar(fontsize = 10, fontface = "bold", family = "Arial"),  # 行标签字体样式
         column_names_gp = gpar(fontsize = 10, fontface = "bold", family = "Arial"),  # 列标签字体样式
         show_colnames = FALSE,  # 隐藏列名称
         treeheight_row = 30, 
         treeheight_col = 30,
         width = 10) -> ht_list

ht_list


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/meta_STL_1203_merger_heatmap_line_cluster.pdf",width = 6, height = 5, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

ht_list #出图

par(opar)
dev.off()



library("magrittr")






####>>>>>1.利用重构的kegg_pathway信号通路进行富集分析<<<<<<<<################


###>>>>>>1.1利用massdatabase下载、读取和转换kegg数据库
library(massdatabase)
library(metpath)

###下载
# download_kegg_pathway(path = "kegg_human_pathway",
#                       sleep = 1,
#                       organism = "hsa")
#                       

# download_h

# download_smpdb_pathway(path = ".")
# download_smpdb_pathway(path = ".")
# ?download_smpdb_pathway

# read_kegg_pathway(path = "kegg_human_pathway")

# data <- read_smpdb_pathway(path = ".")


# data
# data
# 
# convert_hmdb2metid <-
#   function(data,
#            path = ".",
#            threads = 5,
#            ms1_or_ms2 = c("ms1", "ms2")) {
#     dir.create(path, showWarnings = FALSE, recursive = TRUE)
#     ms1_or_ms2 <- match.arg(ms1_or_ms2)
#     
#     if (ms1_or_ms2 == "ms1") {
#       data <-
#         data %>%
#         dplyr::filter(!is.na(monisotopic_molecular_weight))
#       
#       data[which(data == "", arr.ind = TRUE)] <- NA
#       
#       data <-
#         data %>%
#         dplyr::rename(
#           mz = average_molecular_weight,
#           Create_date = creation_date,
#           Updated_date = update_date,
#           Lab.ID = accession,
#           Compound.name = name,
#           Description = description,
#           Synonyms = synonyms,
#           Formula = chemical_formula,
#           IUPAC_name = iupac_name,
#           Traditional_IUPAC_name = traditional_iupac,
#           CAS.ID = cas_registry_number,
#           SMILES.ID = smiles,
#           INCHI.ID = inchi,
#           INCHIKEY.ID = inchikey,
#           Kingdom = kingdom,
#           Super_class = super_class,
#           Class = class,
#           Sub_class = sub_class,
#           State = state,
#           Biospecimen_locations = biospecimen_locations,
#           Cellular_locations = cellular_locations,
#           Tissue_locations = tissue_locations,
#           CHEMSPIDER.ID = chemspider_id,
#           DRUGBANK.ID = drugbank_id,
#           FOODB.ID = foodb_id,
#           PUBCHEM.ID = pubchem_compound_id,
#           CHEBI.ID = chebi_id,
#           KEGG.ID = kegg_id,
#           BIOCYC.ID = biocyc_id,
#           BIGG.ID = bigg_id,
#           WIKIPEDIA.ID = wikipedia_id,
#           METLIN.ID = metlin_id
#         ) %>%
#         dplyr::mutate(
#           HMDB.ID = Lab.ID,
#           RT = NA,
#           mz.pos = NA,
#           mz.neg = NA,
#           Submitter = "HMDB",
#           From_human = "Yes"
#         ) %>%
#         dplyr::select(
#           Lab.ID,
#           Compound.name,
#           mz,
#           RT,
#           CAS.ID,
#           HMDB.ID,
#           KEGG.ID,
#           Formula,
#           mz.pos,
#           mz.neg,
#           Submitter,
#           everything()
#         )
#       
#       temp_file <- tempfile()
#       dir.create(temp_file, showWarnings = FALSE)
#       
#       readr::write_csv(x = data,
#                        file = file.path(temp_file, "data.csv"))
#       
#       hmdb_ms1 <-
#         metid::construct_database(
#           path = temp_file,
#           version =  as.character(Sys.Date()),
#           metabolite.info.name = "data.csv",
#           source = "HMDB",
#           link = "https://hmdb.ca/",
#           creater = "Xiaotao Shen",
#           email = "shenxt@stanford.edu",
#           rt = FALSE,
#           threads = threads
#         )
#       
#       save(hmdb_ms1, file = file.path(path, "hmdb_ms1"))
#       invisible(hmdb_ms1)
#     } else{
#       remove_idx <-
#         data %>%
#         lapply(function(x) {
#           nrow(x$ms2)
#         }) %>%
#         unlist() %>%
#         `==`(0) %>%
#         which()
#       
#       if (length(remove_idx) > 0) {
#         data <-
#           data[-remove_idx]
#       }
#       
#       spectra_info <-
#         data %>%
#         purrr::map(function(x) {
#           x$ms1_info
#         }) %>%
#         dplyr::bind_rows() %>%
#         as.data.frame()
#       
#       spectra_data <-
#         data %>%
#         purrr::map(function(x) {
#           x$ms2
#         })
#       
#       spectra_info[which(spectra_info == "NA", arr.ind = TRUE)] <-
#         NA
#       spectra_info[which(spectra_info == "n/a", arr.ind = TRUE)] <-
#         NA
#       spectra_info[which(spectra_info == "N/A", arr.ind = TRUE)] <-
#         NA
#       
#       spectra_info <-
#         spectra_info %>%
#         dplyr::select(HMDB.ID,
#                       Instrument_type,
#                       Polarity,
#                       collision_energy_voltage,
#                       adduct)
#       
#       remove_idx <-
#         which(is.na(spectra_info$Polarity))
#       
#       if (length(remove_idx) > 0) {
#         spectra_info <-
#           spectra_info[-remove_idx, ]
#         
#         spectra_data <-
#           spectra_data[-remove_idx]
#       }
#       
#       spectra_info <-
#         spectra_info %>%
#         dplyr::mutate(
#           Polarity = case_when(
#             Polarity == "positive" ~ "Positive",
#             Polarity == "negative" ~ "Negative",
#             TRUE ~ Polarity
#           )
#         )
#       
#       spectra_info$Lab.ID <-
#         masstools::name_duplicated(spectra_info$HMDB.ID) %>%
#         paste("shen", sep = "_")
#       
#       spectra_info2 <-
#         spectra_info %>%
#         plyr::dlply(.variables = .(HMDB.ID)) %>%
#         purrr::map(function(y) {
#           if (sum(is.na(y$collision_energy_voltage)) > 0) {
#             y$collision_energy_voltage[is.na(y$collision_energy_voltage)] <-
#               paste("Unknown",
#                     1:length(y$collision_energy_voltage[is.na(y$collision_energy_voltage)]),
#                     sep = "_")
#           }
#           y
#         }) %>%
#         dplyr::bind_rows() %>%
#         as.data.frame()
#       
#       spectra_info2 <-
#         spectra_info2[match(spectra_info$Lab.ID, spectra_info2$Lab.ID), ]
#       
#       spectra_data2 <-
#         1:length(spectra_data) %>%
#         purrr::map(function(i) {
#           x <- spectra_data[[i]]
#           x <- list(x)
#           names(x) <-
#             spectra_info2$collision_energy_voltage[i]
#           x
#         })
#       
#       names(spectra_data2) <- spectra_info2$Lab.ID
#       
#       ######positive mode
#       spectra_info2$Lab.ID == names(spectra_data2)
#       
#       index_pos <- which(spectra_info2$Polarity == "Positive")
#       index_neg <- which(spectra_info2$Polarity == "Negative")
#       
#       spectra_info_pos <- spectra_info2[index_pos, ]
#       spectra_data_pos <- spectra_data2[index_pos]
#       
#       spectra_info_neg <- spectra_info2[index_neg, ]
#       spectra_data_neg <- spectra_data2[index_neg]
#       
#       colnames(spectra_info2)
#       colnames(hmdb_ms1@spectra.info)
#       
#       spectra_info2 <-
#         spectra_info2 %>%
#         dplyr::rename(CE = "collision_energy_voltage")
#       
#       spectra_info2 <-
#         spectra_info2 %>%
#         dplyr::left_join(hmdb_ms1@spectra.info %>% dplyr::select(-Lab.ID),
#                          by = c("HMDB.ID"))
#       
#       temp_file <- tempfile()
#       dir.create(temp_file, showWarnings = FALSE)
#       
#       readr::write_csv(x = spectra_info2,
#                        file = file.path(temp_file, "spectra_info2.csv"))
#       
#       hmdb_ms2 <-
#         metid::construct_database(
#           path = temp_file,
#           version =  as.character(Sys.Date()),
#           metabolite.info.name = "spectra_info2.csv",
#           source = "HMDB",
#           link = "https://hmdb.ca/",
#           creater = "Xiaotao Shen",
#           email = "shenxt@stanford.edu",
#           rt = FALSE,
#           threads = threads
#         )
#       
#       hmdb_ms2@spectra.data$Spectra.positive <- spectra_data_pos
#       hmdb_ms2@spectra.data$Spectra.negative <- spectra_data_neg
#       
#       save(hmdb_ms2, file = file.path(path, "hmdb_ms2"))
#       message("Done.")
#       invisible(hmdb_ms2)
#     }
#     
#   }
# 
# 
# 
# 
# ###下载了HMDB的pathway文件
# # download_smpdb_pathway(path = ".")
# 
# # convert_smpdb2metpath()
# 
# # data <- read_smpdb_pathway(path = ".", only_primarity_pathway = TRUE) 
# 
# # smpdb_pathway_database <- convert_smpdb2metpath(data = data, path = ".")
# 
# 
# # hmdb_human_pathway <- convert_hmdb2metid(data = data, path = ".")
# 

# hmdb_hsa_pathway = 
#   get_hmdb_pathway(threads = 12)
# 
# hmdb_hsa_pathway

# 
# ####HMDB pathway
# data("hmdb_pathway", package = "metpath")
# hmdb_pathway
# 
# get_pathway_class(hmdb_pathway)
# #get the class of pathways
# pathway_class =
#   metpath::pathway_class(hmdb_pathway)
# 
# remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")
# 
# hmdb_pathway =
#   hmdb_pathway[remain_idx]
# 
# hmdb_pathway
# 
# ###HMDB enrichment
# query_id <-
#   temp_data_metabolomics_crest1$HMDB.ID[!is.na(temp_data_metabolomics_crest1$HMDB.ID)] %>%
#   stringr::str_split("\\{\\}") %>%
#   unlist() %>%
#   unique()

# 
# 
# ###>>>>>>1.2导入待分析数据合并之后的hmdb——id
# 
# str(pfas_meta_STL_p_kegg)
# 
# names(pfas_meta_STL_p_kegg)
# 
# hmdb_id_1<- pfas_meta_STL_p_kegg%>% 
#   select(id_kegg)%>%
#   dplyr::filter(!id_kegg=="NA") %>% 
#   as.data.frame()
# colnames(hmdb_id_1) <- NULL
# 
# 
# hmdb_id_2 <- hmdb_id_1[!is.na(hmdb_id_1)]
# hmdb_id_2 
# 
# 
# # ?get_hmdb_pathway 
# 
# # get_kegg_pathway(local = FALSE, organism = "hsa", threads = 3)    ####联网下载最新的kegg通路
# 
# 
# # hmdb_compound = 
#   # get_hmdb_compound(threads = 12)
# # 
# # hmdb_compound
# # This database is downloaded in 2021-03-02
# 
# 
# 
# #get the class of pathways
# pathway_class =
#   metpath::pathway_class(hmdb_pathway)
# 
# remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")
# 
# remain_idx
# 
# hmdb_pathway =
#   hmdb_pathway[remain_idx]
# 
# hmdb_pathway    ###筛选后的pathway
# 
# result <- 
#   enrich_hmdb(
#     query_id = hmdb_id_2,
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = hmdb_pathway,
#     only_primary_pathway = TRUE,
#     p_cutoff = 0.05,
#     p_adjust_method = "BH",
#     threads = 8
#   )
# 
# # Check the result:
# result
# result@result 
# colnames(result@result)  ###包含每个
# result@result$mapped_id[1]
# result@result$mapped_number[1]
# # write.csv(result@result, file = "pathway_result.csv", row.names = FALSE)
# ## Plot to show pathway enrichment  ###显示通路的条形
# 
# 
# result @result %>%
#   dplyr::filter(p_value < 0.05)
# 
# enrich_bar_plot(
#   object = result,
#   x_axis = "p_value_adjust",
#   cutoff = 0.05,   ####设置curoff
#   top = 10   ###
# )
# 
# enrich_scatter_plot(object = result)    ####查看通路的overlape
# 
# enrich_network(object = result)  ###查看某几个pathway之间的关系 
# 
# 
# result_case_control@result %>% 
#   as_tibble() %>% 
#   write.csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_kegg_meta_fdr_merger_liner_HMDB.csv")
# 
# 
# 
# 
# load("hmdb_metabolites/hmdbMS1Database.rda")
# 
# 
# hmdbMS1Database
# 
# write.xlsx(hmdbMS1Database,"hmdb_metabolites/hmdbMS1Database.csv")
# 
# # download_kegg_pathway(path = "kegg_human_pathway",
# #                       sleep = 1,
# #                       organism = "hsa")
# # 


















                      

###读取
data <- 
  read_kegg_pathway(path = "kegg_human_pathway")







####这个函数不起作用的主要原因是pathway_id和pathway_name这里是list,不是向量

convert_kegg2metpath <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    pathway <-
      new(
        Class = "pathway_database",
        database_info = list(source = "KEGG",
                             version = as.character(Sys.Date())),
        pathway_id = unlist(
          lapply(data, function(x) {
            if (is.null(x$pathway_id) || length(x$pathway_id) == 0 || x$pathway_id == "") {
              return(NA)  # 如果 pathway_id 缺失，返回 NA
            } else {
              return(x$pathway_id)  # 否则返回实际值
            }
          }), use.names = FALSE),
        pathway_name =  unlist(
          lapply(data, function(x) {
            if (is.null(x$pathway_name) || length(x$pathway_name) == 0 || x$pathway_name == "") {
              return(NA)  # 如果 pathway_name 缺失，返回 NA
            } else {
              return(x$pathway_name)  # 否则返回实际值
            }
          }), use.names = FALSE),
        describtion = lapply(data, function(x)
          x$describtion),
        pathway_class = lapply(data, function(x)
          x$pathway_class),
        gene_list = lapply(data, function(x)
          x$gene_list),
        compound_list = lapply(data, function(x)
          x$compound_list),
        protein_list = list(),
        reference_list = list(),
        related_disease = lapply(data, function(x)
          x$related_disease),
        related_module = lapply(data, function(x)
          x$related_module)
      )
    
    if (length(pathway@gene_list) == 0) {
      pathway@gene_list <-
        vector(mode = "list",
               length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }
    
    if (length(pathway@compound_list) == 0) {
      pathway@compound_list <-
        vector(mode = "list",
               length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }
    
    if (length(pathway@protein_list) == 0) {
      pathway@protein_list <-
        vector(mode = "list",
               length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }
    kegg_pathway <- pathway
    rm(list = c("pathway"))
    save(kegg_pathway, file = file.path(path, "kegg_pathway"))
    invisible(kegg_pathway)
  }




kegg_human_pathway <-
  convert_kegg2metpath(data = data,
                       path = "kegg_human_pathway",
                       threads = 14)
kegg_human_pathway



pathway_class <-
  metpath::pathway_class(kegg_human_pathway)
# pathway_class

# get_pathway_class(result_case_control)


remain_idx <-
  pathway_class %>%
  unlist() %>% 
  stringr::str_detect("Disease") %>%
  `!`() %>%
  which()

pathway_database_lgm <-
  kegg_human_pathway[remain_idx]


###>>>>>>1.2导入待分析数据合并之后的kegg——id

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

kegg_id_1<- pfas_meta_STL_p_kegg%>% 
  select(id_kegg)%>%
  dplyr::filter(!id_kegg=="NA") %>% 
  as.data.frame()
colnames(kegg_id_1) <- NULL


kegg_id_2 <- kegg_id_1[!is.na(kegg_id_1)]
kegg_id_2 

result_case_control <- 
  enrich_kegg(query_id = kegg_id_2, 
              query_type = "compound", 
              id_type = "KEGG",
              pathway_database = pathway_database_lgm, 
              p_cutoff = 0.05, 
              p_adjust_method = "BH", 
              threads = 32)


result_case_control


library(clusterProfiler)

# ?enrichKEGG()

# dev.off()

# ?enrich_kegg()
# result_case_control
# dev.off()
# cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/picture/Metabolomics/meta_p.pdf",width = 16, height = 9, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
enrich_bar_plot(object = result_case_control,
                x_axis = "p_value",
                cutoff = 0.05) +
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )-> p1

p1





cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_bar_merger_line.pdf",width = 6, height = 5, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

p1   # 使用大写罗马表示

par(opar)
dev.off()






enrich_bar_plot(object = result_case_control,
                x_axis = "p_value_adjust",
                cutoff = 0.2)+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )-> p1_q 

p1_q





cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_fdr_bar_merger_line.pdf",width = 6, height = 5, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

p1_q  # 使用大写罗马表示

par(opar)
dev.off()



# enrich_network(object = result_case_control,
#                x_axis = "p_value_adjust",
#                cutoff = 0.30)
enrich_scatter_plot(result_case_control)





enrich_scatter_plot(object = result_case_control,
                    y_axis = "p_value",
                    y_axis_cutoff = 0.05)+
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color = "red", size = 1) +
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  ) -> p2_scatter

p2_scatter
# result_case_control@parameter


cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_scater_merger_liner.pdf",width = 6, height = 5, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

p2_scatter  # 使用大写罗马表示

par(opar)
dev.off()





enrich_scatter_plot(object = result_case_control,
                    y_axis = "p_value_adjust",
                    y_axis_cutoff = 0.2)+
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color = "red", size = 1) +
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  ) -> p2_q_scatter




cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_fdr_scater_merger_liner.pdf",width = 6, height = 5, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

p2_q_scatter  # 使用大写罗马表示

par(opar)
dev.off()




#ggplot2画气泡图，修改横坐标的名称为“Enrichment Score”，scale_color_gradient设置蓝红配色
library(ggplot2)
# ggplot(data_kegg,aes(x=EnrichmentScore,y=Term))+geom_point(aes(size=Count,color=pValue))+theme_bw()+labs(y="",x="Enrichment Score")+ scale_color_gradient(low="blue",high="red")
# 
# #修改气泡图横向缩放度：coord_fixed(ratio=1/200)
# ggplot(data_kegg,aes(x=EnrichmentScore,y=Term))+
#   geom_point(aes(size=Count,color=pValue))+
#   theme_bw()+
#   labs(y="",x="Enrichment Score")+
#   scale_color_gradient(low="blue",high="red")+
#    coord_fixed(ratio=1/200)



###>>>>>>>>>>>>>>>>>>>>>>>>>绘制气泡图的KEGG通路<<<<<<<<<<<<<<<<<<<<<<<<<<<####
result_case_control
result_case_control
result_case_control@result 

# conclumn(result_case_control@result) 
# write.csv(result_case_control@result, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/picture/Metabolomics/STL/meta_p_pathway_result_tidymass_STL.csv", row.names =T)


result_case_control@result ->data_kegg

data_kegg %>% 
  dplyr::filter(p_value<0.05) -> data_kegg_1
str(data_kegg_1)

data_kegg_1 %>% 
  dplyr::rename("Count"="mapped_number") -> data_kegg_2


# str()

#查看数据列名
colnames(data_kegg_2)
dim(data_kegg_2)

library(ggsci)
#ggplot2画气泡图，修改横坐标的名称为“Enrichment Score”，scale_color_gradient设置蓝红配色
library(ggplot2)
# ggplot(data_kegg,aes(x=EnrichmentScore,y=Term))+
#   geom_point(aes(size=Count,color=pValue))+
#   theme_bw()+
#   labs(y="",x="Enrichment Score")+ scale_color_gradient(low="blue",high="red")

#修改气泡图横向缩放度：coord_fixed(ratio=1/200)  

ggplot(data_kegg_2,aes(x=mapped_percentage,y=pathway_name))+
  geom_point(aes(size=Count,color=p_value),fill = "lightgray")+
  theme_bw()+
  labs(y="",x="Enrichment Score")+
  # scale_fill_material("red")
  # scale_color_gradient(low="blue",high="red")+
  scale_color_gradientn(colors = pal_lancet()(8))+
  # geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color = "red", size = 1) +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 1.5, colour = "black")  # 设置边框线的粗细和颜色
  )+
  coord_fixed(ratio=12/1)  -> p3_bubble

p3_bubble




cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_bubble_merger_liner.pdf",width = 7, height = 6, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

p3_bubble  # 使用大写罗马表示

par(opar)
dev.off()








###>>>>>>>>>>>>>>>>>>>>>>>>>绘制气泡图的KEGG通路<<<<<<<<<<<<<<<<<<<<<<<<<<<####
result_case_control
result_case_control
result_case_control@result 

# conclumn(result_case_control@result) 
# write.csv(result_case_control@result, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/picture/Metabolomics/STL/meta_p_pathway_result_tidymass_STL.csv", row.names =T)


result_case_control@result ->data_kegg

data_kegg %>% 
  dplyr::filter(p_value_adjust<0.20) -> data_kegg_1
str(data_kegg_1)

data_kegg_1 %>% 
  dplyr::rename("Count"="mapped_number") -> data_kegg_2


# str()

#查看数据列名
colnames(data_kegg_2)
dim(data_kegg_2)

library(ggsci)
#ggplot2画气泡图，修改横坐标的名称为“Enrichment Score”，scale_color_gradient设置蓝红配色
library(ggplot2)
# ggplot(data_kegg,aes(x=EnrichmentScore,y=Term))+
#   geom_point(aes(size=Count,color=pValue))+
#   theme_bw()+
#   labs(y="",x="Enrichment Score")+ scale_color_gradient(low="blue",high="red")

#修改气泡图横向缩放度：coord_fixed(ratio=1/200)  

ggplot(data_kegg_2,aes(x=mapped_percentage,y=pathway_name))+
  geom_point(aes(size=Count,color=p_value),fill = "lightgray")+
  theme_bw()+
  labs(y="",x="Enrichment Score")+
  # scale_fill_material("red")
  # scale_color_gradient(low="blue",high="red")+
  scale_color_gradientn(colors = pal_lancet()(8))+
  # geom_hline(yintercept = -log(0.05,10), linetype = "dashed",color = "red", size = 1) +
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.6),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.6, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 1.5, colour = "black")  # 设置边框线的粗细和颜色
  )+
  coord_fixed(ratio=1/3) -> p3_q_bubble

p3_q_bubble


cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_fdr_bubble_merger_liner.pdf",width = 12, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

p3_q_bubble  # 使用大写罗马表示

par(opar)
dev.off()











# 
# 
# 
# library(patchwork)
# p1+p1_q+plot_annotation(tag_levels = "A") # 使用大写罗马表示
# 
# 
# cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/pfas_STL_merge/pfas_STL_meta_mix20241019.pdf",width = 20, height = 9, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# 
# p1+p1_q+plot_annotation(tag_levels = "A") # 使用大写罗马表示
# 
# par(opar)
# dev.off()
# 
# 




# ?enrich_bar_plot
# par(opar)
# dev.off() 
result_case_control

result_case_control@result %>% 
  as_tibble() %>% 
  write.csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_kegg_meta_fdr_merger_liner.csv")








####>>>>>2.根据BHRMA_G计算的系数进行通路分析<<<<<<<<#################

# ####>>>>>2.1找出代谢通路富集出来的每条通路kegg_id对应的化合物<<<<<<##################


###>>>>>>>>>>>>>>>>>>>>>>>>>>>>画环形图<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##########
###>1.读入代谢数据
temp_data_0 <- result_case_control@result %>% 
  as_tibble()

names(temp_data_0)

temp_data_0 %>% 
  slice(1:5) -> temp_data


edge_data <-
  seq_len(nrow(temp_data)) %>%
  purrr::map(function(i) {
    data.frame(
      from = temp_data$pathway_id[i],
      to = stringr::str_split(temp_data$mapped_id[i], ";")[[1]],
      p_value=temp_data$p_value[i],    ###在这里加入绘制边的变量，设定分组和粗细，这是关键步骤
      group=temp_data$pathway_id[i]
    )
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

temp_data
node_data1 <-
  data.frame(
    node = temp_data$pathway_id,
    node_name = temp_data$pathway_name,
    p_value_adjust = -log(temp_data$p_value, 10),
    mapped_percentage = temp_data$mapped_percentage,
    class = "pathway"
  )

names(pfas_meta_STL_p_kegg)

pfas_meta_STL_p_kegg

pfas_meta_STL_p_kegg

node_data2 <-
  unique(unlist(stringr::str_split(temp_data$mapped_id, ";"))) %>%
  data.frame(node = .) %>%
  left_join(pfas_meta_STL_p_kegg[, c("id_kegg", "name")],
            by = c("node" = "id_kegg")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  left_join(pfas_meta_STL_p_kegg[, c("id_kegg", "name")],
            by = c("node" = "id_kegg")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::mutate(name = case_when(
    !is.na(name.x) ~ name.x,
    !is.na(name.y) ~ name.y
  )) %>%
  dplyr::select(-c(name.x,name.y)) %>%
  dplyr::rename(node_name = name) %>%
  dplyr::mutate(
    p_value_adjust = 1,
    mapped_percentage = NA,
    class = "metabolite"
  )

node_data <-
  rbind(node_data1,
        node_data2)


library(tidygraph)
library(ggraph)
library(igraph)

graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE)
# 

# colors <- colorRampPalette(c("red", "orange", "blue"),space = "rgb")(3)

ggraph(graph_data,
       layout = 'linear',
       circular = TRUE) +
  geom_node_point(
    aes(
      color = class,
      size = p_value_adjust
    ),
    shape = 16,
    alpha = 1,
    # size=5,
    show.legend = TRUE
  )  +
  geom_node_text(
    aes(x = x + 0.05 * (x / sqrt(x^2 + y^2)),  # 向外移动
        y = y + 0.05* (y / sqrt(x^2 + y^2)),   # 向外移动
        angle = node_angle(x, y),  ###必需有这个函数，否者不对
        label = node_name
    ),
    size =8,colour = "black",family = "Arial",face = "bold",
    hjust = 'outward',  ###沿着节点的外侧方向移动
    # hjust = 1,  # 文本水平居中
    # vjust = 1,   # 文本在节点上方
    # size = 4,
    # angle = 0,    # 不需要旋转
    nudge_x = 0,  # 不需水平调整
    nudge_y = 0  # 根据需要微调垂直位置
  )  + # 设置点的注释
  scale_size_continuous(range = c(3, 20)) +  #设置点大小范围，可以设置值越小，点越大
  geom_node_point(size = 2,aes(colour = class))+
  scale_color_manual(values = c("#BC3C29FF",  "#0072B5FF"))+
  geom_edge_arc(
    aes(colour = group),  # 将颜色映射到 class
    strength =1.5,
    alpha = 1,
    show.legend = FALSE,
    # color = "blue",
    # lineheight=3,
    # end_cap=circle(1,'cm'),
    # linemitre = 2,
    width=1,
    nudge =0.1  # 缩短连接线的长度，调整此值以获得最佳效果
  )+
  scale_edge_colour_manual(values = c( "#BC3C29FF",  "#0072B5FF",  "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF","#EE4C97FF")) +
  # scale_edge_width_continuous(range = c(3,5)) +
  ggraph::theme_graph()+
  theme(legend.position = c(0.02, 0.5),  # 图例位置（x, y坐标），调整为合适的值
        legend.justification = c(0, 0.5),  # 图例锚点，左侧中心
        # legend.margin = margin(10, 10, 10, 10),  # 调整图例的边距
        legend.background = element_rect(fill = "white", color = NA),  # 图例背景
        legend.box.background = element_rect(color = "white"),  # 图例框线
        # strip.text.y = element_text(angle = 0, hjust = 0,size =14,colour = "black",family = "Arial",face = "bold"),  ###设置每个分面标题的字体及格式
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold")
  ) +  # 调整图例的边距
  scale_x_continuous(expand = expansion(mult = 0.6)) +  # 控制 x 轴的扩展
  scale_y_continuous(expand = expansion(mult = 0.5)) +  # 控制 y 轴的扩展# 调整坐标范围
  # annotate("rect", xmin = -0.1, xmax = 0.1, ymin = -0.1, ymax = 0.1, fill = "black", alpha = 0.2) +
  coord_fixed(ratio = 1) +  # 固定坐标轴比例
  guides(
    color = guide_legend(title = "Class"),  # 修改节点颜色图例标题
    size = guide_legend(title = "-log(p_value, 10)"),  # 修改大小图例标题
    group= guide_legend(title = "Group")
  ) -> p2_circle


p2_circle

?geom_edge_bend


# 查绘图颜色

library(ggsci)
library(tidyverse)
library(cowplot)
library("scales") 

show_col(pal_nejm("default", alpha = 1)(8))  ####查看NEJM杂志的配色，并进行着色


# scale_color_manual(values = c("Steroid hormone biosynthesis" = "#BC3C29FF", "Ovarian steroidogenesis" = "#0072B5FF", "ABC transporters" = "#E18727FF", "Citrate cycle (TCA cycle)" = "#20854EFF", "Ascorbate and aldarate metabolism" = "#7876B1FF", "Valine, leucine and isoleucine biosynthesis" = "#6F99ADFF","Arginine biosynthesis"="#EE4C97FF"))+


?geom_node_text
# dev.off()


library(patchwork)
# p2+plot_annotation(tag_levels = "A") # 使用大写罗马表示



cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_q_circle_merger_line.pdf",width = 20, height = 20, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
p2_circle # 使用大写罗马表示
par(opar)
dev.off()






####>>>>>2.3每条通路kegg_id对应的化合物的多组图（PFAS）<<<<<<##################

temp_data_0 <-result_case_control@result %>% 
  as_tibble()
names(temp_data_0)

temp_data_0 %>% 
  slice(1:5) -> temp_data


edge_data <-
  seq_len(nrow(temp_data)) %>%
  purrr::map(function(i) {
    data.frame(
      from = temp_data$pathway_id[i],
      to = stringr::str_split(temp_data$mapped_id[i], ";")[[1]],
      p_value=temp_data$p_value[i],    ###在这里加入绘制边的变量，设定分组和粗细，这是关键步骤
      group=temp_data$pathway_id[i]
    )
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



node_data1 <-
  data.frame(
    node = temp_data$pathway_id,
    node_name = temp_data$pathway_name,
    p_value_adjust = -log(temp_data$p_value, 10),
    mapped_percentage = temp_data$mapped_percentage,
    class = "pathway"
  )

names(pfas_meta_STL_p_kegg)

pfas_meta_STL_p_kegg

node_data2 <-
  unique(unlist(stringr::str_split(temp_data$mapped_id, ";"))) %>%
  data.frame(node = .) %>%
  left_join(pfas_meta_STL_p_kegg[, c("id_kegg", "name")],
            by = c("node" = "id_kegg")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  left_join(pfas_meta_STL_p_kegg[, c("id_kegg", "name")],
            by = c("node" = "id_kegg")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::mutate(name = case_when(
    !is.na(name.x) ~ name.x,
    !is.na(name.y) ~ name.y
  )) %>%
  dplyr::select(-c(name.x,name.y)) %>%
  dplyr::rename(node_name = name) %>%
  dplyr::mutate(
    p_value_adjust = 1,
    mapped_percentage = NA,
    class = "metabolite"
  )

node_data <-
  rbind(node_data1,
        node_data2)


node_data1
node_data2

edge_data %>% 
  left_join(node_data1,by=c("group"="node")) -> edge_data_2


edge_data_2 %>% 
  left_join(pfas_meta_STL_p_kegg,by=c("to"="id_kegg")) -> meta_STL_mix_kegg_meta_final_0

colSums(is.na(meta_STL_mix_kegg_meta_final_0))



pfas_meta_STL_p_kegg

names(meta_STL_mix_kegg_meta_final_0)



meta_STL_mix_kegg_meta_final_0 ->meta_STL_mix_kegg_meta_final_3

##将调整后的列名称数据进行合并

names(meta_STL_mix_kegg_meta_final_3)

str(meta_STL_mix_kegg_meta_final_3)

meta_STL_mix_kegg_meta_final_3$name

meta_STL_mix_kegg_meta_final_3$estimate



###（4）取出代谢物

####>>>>>在这里将数据进行了log转换，如果数据中estimate.x大于100，直接取对数；如果数据小于0且绝对值大于100，将绝对值除以100然后添加负号取对数；如果数据小于0，且绝对值小于100，将绝对值直接取对数，并添加负号；如果数据大于0，且小于100，直接取对数 ifelse


# write.xlsx(pfas_STL_mix_kegg_meta_final_3,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/pfas_STL_merge/pfas_mix_kegg_meta_fdr_20241105.csv")

# data.frame(estimate=  pfas_STL_mix_kegg_meta_final_3$estimate,
#            estimate_lcl.x=  pfas_STL_mix_kegg_meta_final_3$estimate_lcl,
#            estimate_ucl.x=   pfas_STL_mix_kegg_meta_final_3$estimate_ucl)

meta_STL_mix_kegg_meta_final_3

names(meta_STL_mix_kegg_meta_final_3)


meta_STL_mix_kegg_meta_final_3 %>% 
  dplyr::rename("estimate.x"="Estimate",
                "estimate_lcl.x"="CI2.5",
                "estimate_ucl.x"="CI97.5") -> meta_STL_mix_kegg_meta_final_4


library(tidyverse)
library(cowplot)
library("scales") 

show_col(pal_nejm("default", alpha = 1)(8))  ####查看NEJM杂志的配色，并进行着色

node


meta_STL_mix_kegg_meta_final_4
str(meta_STL_mix_kegg_meta_final_4)
names(meta_STL_mix_kegg_meta_final_4)

names(meta_STL_mix_kegg_meta_final_4)
# 绘图 
meta_STL_mix_kegg_meta_final_4 %>% 
  mutate(
    marker = case_when(
      estimate.x > 0 ~ "Up",
      estimate.x < 0 ~ "Down",
      TRUE ~ "No"
    )
  ) %>%
  ggplot(aes(x = name, y = estimate.x,color=marker)) + 
  geom_point(size = 1.5, position = position_dodge(width = 0.4),shape=15) +           # 绘制β值的点
  geom_errorbar(aes(ymin = estimate_lcl.x, ymax = estimate_ucl.x), 
                width = 0.2, size = 1, position = position_dodge(width = 0.4)) +  # 绘制置信区间
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",size=1) +  # 基线
  labs(x = "Variables", y = "Estimated β (95% CI)") + 
  theme_classic() + 
  # geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
  scale_color_manual(values=c("#377eb8","#e41a1c","#525252")) + 
  # geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  # geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
  # theme(axis.title.x = element_text(size = 12), 
  #       axis.title.y = element_text(size = 12), 
  #       plot.title = element_text(size = 14, face = "bold")) +
  # ggtitle("Estimates of β with 95% Confidence Intervals by Group") +
  # scale_color_manual(values = c("PFAS" = "blue", "STL" = "orange")) + # 自定义颜色
  scale_fill_manual(values=c("#377eb8","#e41a1c","#525252"))+
  facet_grid(node_name ~ ., scales = "free_y", space = "free_y") + ###分面
  # scale_color_manual(values = c("PFAS" = "#BC3C29FF")) +  # 替换为你的组名和颜色
  # scale_fill_manual(values = c("PFAS" = "blue", "STL" = "red"))+
  # scale_x_discrete(labels = function(name.x) substring(name.x, 1, 15))+   ###只取x轴变量名称的前20个字母
  # scale_fill_manual(values = c( "black", "grey50", "white","green"), name = NULL) +
  scale_shape_manual(values = c(2,2), name = NULL) +
  # scale_linetype_manual(values = c(1,1,1), name = NULL) +
  # scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
  # scale_y_continuous(limits = c(-2.4,2.4)) +
  coord_flip(clip = "off") +    
  ylab("Estimate (95%CI)") +
  theme(axis.title.y = element_blank(),  
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"), 
        panel.background = element_rect(fill="grey95"),
        # legend.position = "none",    ###去除图列
        legend.direction = "vertical",
        # legend.direction = "horizontal",  # 设置图例为水平
        legend.justification = c(0.5, 0.5),  # 图例居中
        legend.title = element_text(size = 14, family = "Arial", face = "bold"),
        legend.text = element_text(size = 14, family = "Arial", face = "bold"),
        # strip.background = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",  # 将分面标题放置在外部
        strip.text.y = element_text(angle = 0, hjust = 0,size =14,colour = "black",family = "Arial",face = "bold"),  ###设置每个分面标题的字体及格式
        text = element_text(size = 15),
        axis.line = element_line(linewidth = 0.7),
        axis.ticks = element_line(linewidth = 0.7),
        panel.grid = element_blank(),
        # axis.line.x = element_line(color = "black"),
        # axis.line.y = element_line(color = "black"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold")) -> p3
p3





cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_path_merger_line.pdf",width = 16, height = 12, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

p3 # 使用大写罗马表示
par(opar)
dev.off()





###>>>>>>>>>>将PFAS与代谢和代谢与端粒的相关性均作在同一个图中######### 



###>>>>>>>>>>>>>>>>>>>>>>三、读入经过BHRMA筛选后的pfas混合物---混合线性模型的<<<<<<<<<<<<<<<<<<

df_pfas_mix_meta <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value_mix_liner.xlsx")


###1.提取经过BH校正之后q_value<0.2的代谢物 
###
names(df_pfas_mix_meta)

df_pfas_mix_meta %>% 
  filter(q_value<0.2) -> df_pfas_mix_meta_q   ###总共获得95个代谢物质   


df_pfas_mix_meta_q 

str(df_pfas_mix_meta_q)
names(df_pfas_mix_meta_q)


###2.筛选出代谢与端粒有意义的，对应的pfas化合物的系数


meta_STL_mix_kegg_meta_final_4 %>%
  left_join(df_pfas_mix_meta_q,by="metabolite") -> meta_STL_mix_kegg_meta_final_4_pfas




###3.将横向数据整理为纵向数据、进行绘图

###（1）取出STL数据 
meta_STL_mix_kegg_meta_final_4_pfas %>% 
  select(metabolite,node_name,name,estimate.x,estimate_lcl.x,estimate_ucl.x,p.value,q_value.x) %>% 
  mutate(group=rep("STL",8))-> meta_STL_mix_kegg_meta_final_4_pfas_1
names(meta_STL_mix_kegg_meta_final_4_pfas_1)



###（2）取出PFAS数据
meta_STL_mix_kegg_meta_final_4_pfas %>% 
  select(metabolite,node_name,name,estimate.y,estimate_lcl,estimate_ucl,p.val,q_value.y) %>% 
  mutate(group=rep("PFAS",8))-> meta_STL_mix_kegg_meta_final_4_pfas_2

names(meta_STL_mix_kegg_meta_final_4_pfas_2) 






colnames(meta_STL_mix_kegg_meta_final_4_pfas_2) <- colnames(meta_STL_mix_kegg_meta_final_4_pfas_1)    ###将pfas_STL_mix_kegg_meta_final_2的名称调整为colnamespfas_STL_mix_kegg_meta_final_1

meta_STL_mix_kegg_meta_final_4_pfas_2

###（3）将PFAS与STL进行列合并
rbind(meta_STL_mix_kegg_meta_final_4_pfas_2,meta_STL_mix_kegg_meta_final_4_pfas_1) -> meta_STL_mix_kegg_meta_final_4_pfas_3 ##将调整后的列名称数据进行合并

names(meta_STL_mix_kegg_meta_final_4_pfas_3)






# names(meta_STL_mix_kegg_meta_final_4)
# 绘图 
meta_STL_mix_kegg_meta_final_4_pfas_3 %>% 
  mutate(
    marker = case_when(
      estimate.x > 0 ~ "Up",
      estimate.x < 0 ~ "Down",
      TRUE ~ "No"
    )
  ) %>%
  ggplot(aes(x = name, y = estimate.x,color = group)) + 
  geom_point(size = 1.5, position = position_dodge(width = 0.4),shape=15) +           # 绘制β值的点
  geom_errorbar(aes(ymin = estimate_lcl.x, ymax = estimate_ucl.x), 
                width = 0.2, size = 1, position = position_dodge(width = 0.4)) +  # 绘制置信区间
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",size=1) +  # 基线
  labs(x = "Variables", y = "Estimated β (95% CI)") + 
  theme_classic() + 
  # geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  # geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
  # theme(axis.title.x = element_text(size = 12), 
  #       axis.title.y = element_text(size = 12), 
  #       plot.title = element_text(size = 14, face = "bold")) +
  # ggtitle("Estimates of β with 95% Confidence Intervals by Group") +
  # scale_color_manual(values = c("PFAS" = "blue", "STL" = "orange")) + # 自定义颜色
  # scale_fill_manual(values = c("PFAS" = "#BC3C29FF", "STL" = "#0072B5FF"))+
  # scale_fill_manual(values=c("#377eb8","#e41a1c","#525252"))+
  facet_grid(node_name ~ ., scales = "free_y", space = "free_y") + ###分面
  # geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
  scale_color_manual(values=c("#377eb8","#e41a1c")) + 
  scale_shape_manual(values = c(16,19)) +  # 修改分组的形状
  # scale_color_manual(values = c("PFAS" = "#BC3C29FF")) +  # 替换为你的组名和颜色
  # scale_fill_manual(values = c("PFAS" = "blue", "STL" = "red"))+
  # scale_x_discrete(labels = function(name.x) substring(name.x, 1, 15))+   ###只取x轴变量名称的前20个字母
  # scale_fill_manual(values = c( "black", "grey50", "white","green"), name = NULL) +
  
  # scale_linetype_manual(values = c(1,1,1), name = NULL) +
  # scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
  # scale_y_continuous(limits = c(-2.4,2.4)) +
  coord_flip(clip = "off") +    
  ylab("Estimate (95%CI)") +
  theme(axis.title.y = element_blank(),  
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"), 
        panel.background = element_rect(fill="grey95"),
        # legend.position = "none",    ###去除图列
        legend.direction = "vertical",
        # legend.direction = "horizontal",  # 设置图例为水平
        legend.justification = c(0.5, 0.5),  # 图例居中
        legend.title = element_text(size = 14, family = "Arial", face = "bold"),
        legend.text = element_text(size = 14, family = "Arial", face = "bold"),
        # strip.background = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",  # 将分面标题放置在外部
        strip.text.y = element_text(angle = 0, hjust = 0,size =14,colour = "black",family = "Arial",face = "bold"),  ###设置每个分面标题的字体及格式
        text = element_text(size = 15),
        axis.line = element_line(linewidth = 0.7),
        axis.ticks = element_line(linewidth = 0.7),
        panel.grid = element_blank(),
        # axis.line.x = element_line(color = "black"),
        # axis.line.y = element_line(color = "black"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold")) -> p4
p4

# dev.off()





cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_path_merger_line_merger.pdf",width = 14, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

p4 # 使用大写罗马表示
par(opar)
dev.off()
























