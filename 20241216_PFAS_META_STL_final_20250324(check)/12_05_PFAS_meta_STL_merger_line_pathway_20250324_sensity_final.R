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
pfas_meta_STL <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_all_meta_STL_zscore_liner_q_PFAS.xlsx")
names(pfas_meta_STL)

pfas_meta_STL$meta

pfas_meta_STL%>% 
  dplyr::rename(metabolite="meta") ->pfas_meta_STL_1  


###过滤出PFAS-META-STL<0.2的Q值
pfas_meta_STL_1 %>% 
  filter(q_value<0.2) -> pfas_meta_STL_2



pfas_meta_STL_2 %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg

# write.xlsx(pfas_meta_STL_p_kegg,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_stl_q.xlsx")

str(pfas_meta_STL_p_kegg)

names(pfas_meta_STL_p_kegg)

colSums(is.na(pfas_meta_STL_p_kegg))

# 
# df <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy.xlsx")
# 


library("magrittr")


#######所有化合物的火山图


marker_color <-
  c(
    "Up" = ggsci::pal_futurama()(n = 9)[6],
    "Down" = ggsci::pal_futurama()(n = 9)[3],
    "No" = ggsci::pal_futurama()(n = 9)[9]
  )

show_col(ggsci::pal_futurama()(n = 9))


pfas_meta_STL_p_kegg


names(pfas_meta_STL_p_kegg)

library(showtext)
library(sysfonts)

font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

pfas_meta_STL_1 %>% 
  left_join(df_unique_3,by="metabolite") -> pfas_meta_STL_p_kegg_2


volcano_plot_all <-
  pfas_meta_STL_p_kegg_2 %>% 
  mutate(marker=case_when(
    q_value < 0.2 & estimate > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.2 & estimate< 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate, y = -log10(q_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, linetype = 2,size=1,color="#e41a1c") +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
  labs(x = "Differential metabolites", y = "-log10(q value)")+
  theme_bw() +
  theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # legend.position = c(0.25,0.15),
        legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
        # axis.line = element_line(linewidth = 0.9),
        # axis.ticks = element_line(linewidth = 0.9),
        axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
        panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
        panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
        panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
  )+
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) 

# ggplot(aes(x =estimate, y = -log10(q_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
# # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
# geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
# # ylim(c(0,10))+
# # xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围
# scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
# geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
# # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线
# geom_hline(yintercept = 0, linetype = 2,size=1) +
# # facet_wrap(~var.names, scales = "free") +
# labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
# # ggtitle("单组火山图") + #标题
# labs(x = "Estimate", y = "-log10(q_value)")+
# theme_bw() +
# theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#       text = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#       # legend.position = c(0.25,0.15),
#       legend.title = element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#       axis.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
#       # axis.line = element_line(linewidth = 0.9),
#       # axis.ticks = element_line(linewidth = 0.9),
#       axis.title.y = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
#       axis.title.x = element_text(size = 14,colour = "black",family = "Arial",face = "bold"),
#       panel.grid.major = element_line(size = 0.8),  # 主网格线粗细 , colour = "red"
#       panel.grid.minor = element_line(size = 0.8, linetype = "dashed"), # 次网格线粗细,
#       panel.border = element_rect(size = 0.8, colour = "black")  # 设置边框线的粗细和颜色
# )+
# geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) 


volcano_plot_all

# dev.off()


# 
# 
# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/volcano_pfas_meta_STL_merger_line.pdf",width = 6, height = 5, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_all #出图

par(opar)
dev.off()







#######2.绘制热图
###读入原始的数据

library(ComplexHeatmap)
df_rawdata <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy_623.xlsx")


names(df_rawdata)

###筛选出q<0.05的结果

df_rawdata %>% 
  select(number:STL,pfas_meta_STL_p_kegg$metabolite) -> df_rawdata_q

df_rawdata_q

names(df_rawdata_q)




df_rawdata_q %>% 
  select(1,M661:M609) -> df_meta_final

df_meta_final %>% 
  t() -> df_meta_final_4


###将第一行作为列名
colnames(df_meta_final_4) <- df_meta_final_4[1,]

df_meta_final_4[-1,] -> df_meta_final_5


###对行重新生成行名
rownames(df_meta_final_5) <- pfas_meta_STL_p_kegg$name





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

pheatmap(df_meta_final_5, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         cutree_rows = 6, 
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

cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/meta_STL_1203_merger_heatmap_line.pdf",width = 6, height = 5, pointsize = 20,family = "Arial")
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


# 
# pathway_class <-
#   metpath::pathway_class(kegg_human_pathway)
# pathway_class

# get_pathway_class(result_case_control)


remain_idx <-
  kegg_human_pathway@pathway_class %>%
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





####>>>>>>>>>>>>>>>>>>>>>>>>>绘制相关网络图<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


##构建object文件 
###1.1首先构建expression_data 
###

df_background_DNA_LT_meta_mean_final <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/result_pfas_rawdata_623.xlsx")  ###此处的pfas和端粒长度均取对数了，并未进行标化或归一化

names(df_background_DNA_LT_meta_mean_final)   ###精液质量和代谢组学数据已经进行取对数和标化
str(df_background_DNA_LT_meta_mean_final)

df_background_DNA_LT_meta_mean_final[,c(1,25:295)] %>% 
  t() -> expression_data_1


##第一行作为列名
colnames(expression_data_1) <- expression_data_1[1,]

##删除第一行

expression_data_1[-1,]  %>% 
  as.data.frame()-> expression_data

names(df_background_DNA_LT_meta_mean_final)

###1.2首先构建sample_info
sample_info <- data.frame(df_background_DNA_LT_meta_mean_final[,c(1:24)] %>% 
                            dplyr::rename(sample_id="number"),class = "Subject",
                          subject_id = "Subject")

###1.3其次构建variable_info
# variable_info <- rownames(expression_data) %>% as.data.frame() %>%
#   dplyr::rename(variable_id=".")
# 
variable_info <- df_unique_3  %>%
  select(metabolite)%>%
  dplyr::rename(variable_id="metabolite") %>% 
  as.data.frame()

# 将第一列设置为行名称
row.names(variable_info) <- variable_info$variable_id


variable_info
# df_unique_3  


object_1 <- 
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

object_1

dim(object_1)

# 
# object %>% 
#   activate_mass_dataset(what = "sample_info") %>% 
#   # filter(sample_id %in% c(control_sample_id)) %>% 
#   activate_mass_dataset(what = "variable_info")  -> a
# 
# a@variable_info


temp_object <- object_1 %>% 
  # activate_mass_dataset(what = "sample_info") %>% 
  # filter(sample_id %in% c(control_sample_id)) %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  filter(variable_id %in% pfas_meta_STL_p_kegg$metabolite) %>% 
  `+`(1) %>% 
  log(2) %>% 
  scale()

str(pfas_meta_STL_p_kegg)

extract_variable_info(temp_object) 

temp_object

library(ggraph)
library(tidygraph)
graph_data <- convert_mass_dataset2graph(
  object = temp_object,
  margin = "variable",
  cor_method = "spearman",
  p_adjust_cutoff = 0.05,
  p_value_cutoff = 0.05,
  pos_cor_cutoff = 0.7,
  neg_cor_cutoff = -0.7
) %>% 
  mutate(Degree = centrality_degree(mode = 'all'))


graph_data




# 提取 graph_data 中的 node 列
graph_nodes <- graph_data %>%
  activate(nodes) %>%
  as_tibble()


###重新提取出来temp_object包涵的数据中的代谢物变量和name

df_unique_3  %>%
  select(metabolite,name)%>%
  dplyr::rename(node="metabolite") %>% 
  dplyr::filter(node  %in% pfas_meta_STL_p_kegg$metabolite) -> new_data  


# 将 new_data 与 graph_nodes 进行合并
merged_data <- left_join(graph_nodes, new_data, by = "node") %>% select(-2)

merged_data    ###生成node对应的name变量


# 将 merged_data 与 graph_data 的节点数据进行合并
graph_data %>%
  activate(nodes) %>%
  left_join(merged_data, by = "node") -> graph_data_2


# 查看结果
# print(merged_data)






# 查看结果
# print(node_list)



library(extrafont)
loadfonts()
library(showtext)
showtext_auto()



library(ggraph)
library(ggplot2)
library(shadowtext)

# ggraph(graph = graph_data_2, layout = "kk") +
#   geom_edge_fan(
#     aes(color = correlation),
#     width = 2,  # 将线加粗 2 倍
#     show.legend = TRUE
#   ) +
#   geom_node_point(
#     aes(size = Degree * 10)  # 将点放大 2 倍
#   ) +
#   shadowtext::geom_shadowtext(
#     aes(x = x, y = y, label = name),
#     bg.colour = "white",
#     colour = "black",
#     size = 12,  # 设置字号为 12
#     family = "Arial"  # 设置字体为 Arial
#   ) +
#   theme_graph(base_family = "Arial") +  # 设置全局字体为 Arial
#   scale_edge_color_gradient2(low = "darkblue", mid = "white", high = "red")



library(ggraph)
library(ggplot2)
library(shadowtext)


ggraph(graph = graph_data_2, layout = "kk") +
  geom_edge_fan(
    aes(color = correlation),
    width = 2,  # 将线加粗 2 倍
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(size = Degree), # 将点放大 10 倍
    shape=21,
    fill = "#377eb8" ,# 设置填充颜色为淡蓝色
    color="#377eb8"
  ) +
  scale_size(range = c(5, 18)) + # 设置点的最小和最大大小
  ggrepel::geom_text_repel(
    aes(x = x, y = y, label = name),  # 设置标签的位置和内容
    color = "black",  # 设置文字颜色
    size = 7,  # 设置字号为 12
    family = "Arial",  # 设置字体为 Arial
    box.padding = 0.5,  # 标签之间的间距
    point.padding = 0.5,  # 标签与点之间的间距
    max.overlaps = Inf,  # 允许最大重叠次数（Inf 表示无限制）
    force = 1,  # 调整标签之间的排斥力
    force_pull = 1,  # 调整标签与点之间的吸引力
    segment.color = "grey50",  # 连接线的颜色
    segment.size = 0.5,  # 连接线的粗细
    min.segment.length = 0.1  # 连接线的最小长度
  ) +
  # guides(
  #   size = guide_legend(override.aes = list(color = "red", fill = "red"))  # 将横线颜色改为与圆点颜色一致
  # )+
  theme_graph(base_family = "Arial") +  # 设置全局字体为 Arial
  scale_edge_color_gradient2(low = "darkblue", mid = "white", high = "#e41a1c")+
  theme(
    legend.key = element_blank(),  # 隐藏图例中点的背景和横线
    legend.text = element_text(family = "Arial", size = 20),  # 设置图例字体为 Arial，12 号
    legend.title = element_text(family = "Arial", size = 20)  # 设置图例标题字体为 Arial，12 号
  ) -> plot


plot



# dev.off()


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pathway/pfas_meta_STL_20241218_fdr_network_liner.pdf",width = 16, height = 9,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)

plot  # 使用大写罗马表示

par(opar)
dev.off()








####>>>>>2.根据BHRMA_G计算的系数进行通路分析<<<<<<<<#################

# ####>>>>>2.1找出代谢通路富集出来的每条通路kegg_id对应的化合物<<<<<<##################


###>>>>>>>>>>>>>>>>>>>>>>>>>>>>画环形图<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##########
###>1.读入代谢数据
temp_data_0 <- result_case_control@result %>% 
  as_tibble()

names(temp_data_0)

temp_data_0 %>% 
  arrange(p_value) %>% 
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
  arrange(p_value) %>% 
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
  mutate(group=rep("STL",10))-> meta_STL_mix_kegg_meta_final_4_pfas_1
names(meta_STL_mix_kegg_meta_final_4_pfas_1)



###（2）取出PFAS数据
meta_STL_mix_kegg_meta_final_4_pfas %>% 
  select(metabolite,node_name,name,estimate.y,estimate_lcl,estimate_ucl,p.val,q_value.y) %>% 
  mutate(group=rep("PFAS",10))-> meta_STL_mix_kegg_meta_final_4_pfas_2

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






























