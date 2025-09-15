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

getwd()
write.table(
  omics_data_scaled,
  file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/cluster/temp_data.txt",
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

names(pfas_meta_STL_p_kegg_clustr)

write.xlsx(pfas_meta_STL_p_kegg_clustr,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/cluster/pfas_meta_STL_p_kegg_clustr.xlsx")


########进一步绘制分类图形  






####>>>>>1.利用重构的kegg_pathway信号通路进行富集分析<<<<<<<<################


###>>>>>>1.1利用massdatabase下载、读取和转换kegg数据库
library(massdatabase)
library(metpath)


###读取
data <- 
  read_kegg_pathway(path = "kegg_human_pathway")



getwd()



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

str(pfas_meta_STL_p_kegg_clustr)

names(pfas_meta_STL_p_kegg_clustr)

kegg_id_1<- pfas_meta_STL_p_kegg_clustr %>% 
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



result_case_control@result %>% 
  as_tibble() 



pfas_meta_STL_p_kegg_clustr



















cor_data <- 
  read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_meta_stl_cor.xlsx")

library(org.Hs.eg.db)
library(clusterProfiler)

###Cluster 1




result_all_cluster1 <-
  tryCatch(
    expr = {
      load("cluster_1/cluster1/result1")
      result_all
    },
    error = function(e) {
      NULL
    }
  )







