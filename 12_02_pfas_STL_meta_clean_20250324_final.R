# no_source() 
rm(list = ls())
pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel)
# get_project_wd()
# masstools::setwd_project()

getwd()


###-----------------------一、读入线性模型的pfas_meta数据---------------------------------

# D:\OneDrive\Postdoctoral_SCI\Sperm_bank\result\20240407_PFOA_sperm\Multi_linear_regression_20240407\metabolomics\pfas_meta_stl_20241216\pfas_meta_result
load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/coefs_PFAS_meta_mean.Rdata")

g_coefs_PFAS_meta_mean

lapply(seq_along(g_coefs_PFAS_meta_mean), function(i) {
  g_coefs_PFAS_meta_mean[[i]]  %>% 
    as.data.frame() 
  # res1 %>% 
  #   dplyr::filter(var.names=="pfos.beta"|var.names=="pfhxs.beta"|var.names=="pfhps.beta"|var.names=="pfoa.beta"|var.names=="pfna.beta"|var.names=="pfda.beta") -> res2
  # BCI_CI95 <- paste0(round(res1$Mean,2)," (",round(res1$X2.5.,2),', ',round(res1$X97.5.,2),")")
  # Uni_glm_model <- data.frame(res1,
  #                             "BCI_CI95"=BCI_CI95
  #                           )
  # return(Uni_glm_model)
}) %>% 
  bind_rows() %>%
  write.csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/coefs_pfas_meta.csv")






####>>>>>>>>>>>>>>>>>>>>>>二、读入混合线性模型的pfas_meta数据<<<<<<<<<<<<<<<#####


load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/coefs_PFAS_mix_liner.Rdata")

g_coefs_PFAS_meta_mix_liner


lapply(seq_along(g_coefs_PFAS_meta_mix_liner), function(i) {
  g_coefs_PFAS_meta_mix_liner[[i]]  %>% 
    as.data.frame() 
  # res1 %>% 
  #   dplyr::filter(var.names=="pfos.beta"|var.names=="pfhxs.beta"|var.names=="pfhps.beta"|var.names=="pfoa.beta"|var.names=="pfna.beta"|var.names=="pfda.beta") -> res2
  # BCI_CI95 <- paste0(round(res1$Mean,2)," (",round(res1$X2.5.,2),', ',round(res1$X97.5.,2),")")
  # Uni_glm_model <- data.frame(res1,
  #                             "BCI_CI95"=BCI_CI95
  #                           )
  # return(Uni_glm_model)
}) %>% 
  bind_rows() %>%
  write.csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/g_coefs_PFAS_meta_mix_liner.csv")








# 
# 
# ###-----------------------读入pfas_meta数据---------------------------------
# 
# load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_pfas_meta.Rdata")
# 
# lapply(seq_along(coefs_pfas_meta), function(i) {
#     coefs_pfas_meta[[i]]  %>% 
#       as.data.frame() 
#     # res1 %>% 
#     #   dplyr::filter(var.names=="pfos.beta"|var.names=="pfhxs.beta"|var.names=="pfhps.beta"|var.names=="pfoa.beta"|var.names=="pfna.beta"|var.names=="pfda.beta") -> res2
#     # BCI_CI95 <- paste0(round(res1$Mean,2)," (",round(res1$X2.5.,2),', ',round(res1$X97.5.,2),")")
#     # Uni_glm_model <- data.frame(res1,
#     #                             "BCI_CI95"=BCI_CI95
#     #                           )
#     # return(Uni_glm_model)
#   }) %>% 
#   bind_rows() %>%
#   write.csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_pfas_meta.csv")
# 
# 
# 
# 
# ###-----------------------读入STL_meta数据---------------------------------
# 
# load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_STL_meta.Rdata")
# 
# STL_coefs
# lapply(seq_along(STL_coefs), function(i) {
#   STL_coefs[[i]]  %>% 
#     as.data.frame() 
#   # res1 %>% 
#   #   dplyr::filter(var.names=="pfos.beta"|var.names=="pfhxs.beta"|var.names=="pfhps.beta"|var.names=="pfoa.beta"|var.names=="pfna.beta"|var.names=="pfda.beta") -> res2
#   # BCI_CI95 <- paste0(round(res1$Mean,2)," (",round(res1$X2.5.,2),', ',round(res1$X97.5.,2),")")
#   # Uni_glm_model <- data.frame(res1,
#   #                             "BCI_CI95"=BCI_CI95
#   #                           )
#   # return(Uni_glm_model)
# }) %>% 
#   bind_rows() %>%
#   write.csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_STL_meta.csv")
# 
# 
# 

