# no_source() 
pacman::p_load(BiocManager,ComplexHeatmap,ggraph,tidygraph,extrafont,shadowtext,
               tidyverse,openxlsx,readr,readxl,tidymass,ggplot2,labelled,
               rstatix,ggpubr,GGally,car,Epi,lme4,lmerTest,emmeans,geepack,
               ggeffects,VIM,mice,gratia,ggrepel,showtext,sysfonts)
# get_project_wd()
# masstools::setwd_project()
rm(list = ls())
getwd()
###-------将交我算-PHC跑的数据加载进来

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





#将导出的数据重新读入数据  


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
  as_tibble() -> coefs_tran
# coefs_tran <- read_csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/coefs_pfas_meta.csv")





####将生成的代谢物传递给新的数据
coefs_tran ->coefs_tran_tidy_0

str(coefs_tran_tidy_0)

length(coefs_tran_tidy_0$metabolite)
table(coefs_tran_tidy_0$metabolite) 
coefs_tran_tidy_0$metabolite

coefs_tran_tidy_0[coefs_tran_tidy_0$metabolite %in% coefs_tran_tidy_0$metabolite[duplicated(coefs_tran_tidy_0$metabolite)], ]


coefs_tran_tidy_0$metabolite


###取出代谢物质单独的变量名称
coefs_tran_tidy_0 %>%
  group_by(metabolite) %>%
  dplyr::summarise(count=n())  -> coefs_tran_tidy_unique

coefs_tran_tidy_unique$metabolite

names(coefs_tran_tidy_unique)

str(coefs_tran_tidy_unique)

###将代谢物重新传递给coefs_tran_tidy

coefs_tran_tidy_0 ->coefs_tran_tidy 

coefs_tran_tidy


results_df <- coefs_tran_tidy  %>% 
  # select(-1) %>% 
  dplyr::mutate(p_value_estimate = na_if(coefs_tran_tidy$p.val, NA) %>% as.numeric)   ###去除p值中含有缺失的值

results_df



## Subset eta ------------------------------------
results_df_2 <- results_df %>%
  # dplyr::filter(str_detect(var.names,"eta")) %>%
  dplyr::rename(estimate = Mean,sd=SD,estimate_lcl=X2.5., estimate_ucl=X97.5.)
# mixtures_eta

# estimate, sd, lcl, ucl, p_value

results_df_2$var.names

## Pivot results wider ---------------------------------


# ## Pivot results wider ---------------------------------
# results_df_w <- results_df_2 %>% 
#   filter(var.names != "eta") %>%
#   tidylog::pivot_wider(id_cols = c(cohort, mixture, metabolite, exposure), 
#                        names_from = term, 
#                        values_from = c(estimate, sd, lcl, ucl, p_value)) %>% 
#   dplyr::select(cohort, mixture, metabolite, exposure, 
#                 contains("beta"), contains("pip")) %>% 
#   dplyr::select(-p_value_pip) %>% 
#   dplyr::rename(p_value = p_value_beta)
# 
# 


###--------转换BCI CI95-----------------------


results_df_3 <- results_df_2 %>% 
  dplyr::filter(var.names %in% c("eta.high","eta.low","psi"))  ###提取出eta.high、eta.low、psi三个重要变量

results_df_3

results_df_3 %>%  
  filter(var.names=="psi") -> results_df_4

paste0(round(results_df_4$estimate,2)," (",round(results_df_4$estimate_lcl,2),', ',round(results_df_4$estimate_ucl,2),")") %>% 
  as_tibble() %>% 
  dplyr::rename(BCI_CI95=value)  %>% 
  cbind(results_df_4)-> results_df_BCI_95CI  ###根据需要将想要的变量提取处理 

str(results_df_BCI_95CI)

str(results_df_2) 
results_df_2

results_df_2 %>% 
  dplyr::filter(!var.names %in% c("eta.high","eta.low","psi")) %>% 
  dplyr::select(metabolite,everything())->results_df_5   ###过滤掉不是eta.high、eta.low、psi三个重要变量


###将转换后的BCI_95与不含eta.high、eta.low、psi三个变量进行合并
results_df_5

results_df_5 %>% 
  left_join(results_df_BCI_95CI,by="metabolite") %>% 
  dplyr::filter(var.names.x %in% c("PFOS.beta","PFOA.beta","PFDA.beta","PFUdA.beta","L_PFHxS.beta","PF3ONS_9CL.beta"))-> results_df_6   ###合并后的数据  ###只保留beta变量,每个化合物的beta和混合物beta

results_df_6 
names(results_df_6)


###-------------只保留后验概率PIP数值-------------------------
results_df_5 %>% 
  dplyr::filter(var.names %in% c("PFOS.gamma","PFOA.gamma","PFDA.gamma","PFUdA.gamma","L_PFHxS.gamma","PF3ONS_9CL.gamma"))  %>% 
  tidylog::pivot_wider( id_cols = c(metabolite), 
                        names_from = var.names,
                        values_from = c(estimate)) -> results_df_7  ###保存每个化合物的后验概率

results_df_7   




# Calculate new p values


# Calculate new p values  计算单一物质的wald检验
names(results_df_6)
results_df_6 %>% 
  dplyr::mutate(wald = (estimate.x/sd.x), 
                # p = 2*(1-pnorm(abs(wald),0,1)),
                p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
                p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
                neg_log_p = -log10(p_value))  %>% 
  group_by(var.names.x) %>% 
  dplyr::mutate(q_value=p.adjust(p_value,method = "fdr"),
                significance = if_else(p_value < 0.05,
                                       "p < 0.05", 
                                       "Not Sig."), 
                significancefdr = if_else(q_value < 0.05,
                                          "q < 0.05",
                                          "Not Sig.")) %>% 
  ungroup() -> results_df_only_p_value

results_df_only_p_value

names(results_df_only_p_value)

write.xlsx(results_df_only_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_only_p_value.xlsx") 




# result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_pfas_meta





########---------------计算混合物质的wald检验---------------------
# Calculate new p values   
# 
results_df_BCI_95CI
names(results_df_BCI_95CI)

results_df_BCI_95CI %>% 
  dplyr::mutate(wald = (estimate/sd), 
                # p = 2*(1-pnorm(abs(wald),0,1)),
                p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
                p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
                neg_log_p = -log10(p_value)) %>% 
  # group_by(metabolite) %>% 
  dplyr::mutate(q_value = p.adjust(p_value,
                                   method = "fdr"), 
                significance=if_else(p_value < 0.05,
                                     "p < 0.05", 
                                     "Not Sig."), 
                significancefdr=if_else(q_value < 0.05,
                                        "q < 0.05",
                                        "Not Sig.")) %>% 
  ungroup() -> results_df_mix_p_value

results_df_mix_p_value

names(results_df_mix_p_value)

write.xlsx(results_df_mix_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value.xlsx")


# 
# q_value = p.adjust(c(0.896,0.928,0.868),
#                    method = "fdr")



###--------------将BCI_CI95与后验概率PIP进行合并-------------


results_df_mix_p_value %>% 
  left_join(results_df_7 ,by="metabolite") %>% 
  dplyr::select(metabolite,BCI_CI95,p_value,q_value,PFOS.gamma:PF3ONS_9CL.gamma,everything())-> results_df_final   ###得到复现文章的Supplemental Excel file文件的结果

write.xlsx(results_df_final,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_final.xlsx")



# install.packages("devtools")
# library("devtools")
# install_github("JAGoodrich/jag2")

# install.packages("jag2")

library("jag2")




names(results_df_only_p_value)


results_df_only_p_value %>% 
  mutate(associated_bci = if_else(estimate_lcl.x > 0 | estimate_ucl.x < 0,
                                  "Associated", 
                                  "Not Associated")) %>% 
  group_by(metabolite, var.names.x) %>% 
  dplyr::summarise(n_features = length(metabolite), 
                   percent_significant_p05 = jag2::npct(significance, 
                                                        "p < 0.05", 
                                                        n.digits = 2), 
                   percent_significant_q05 = jag2::npct(significancefdr, 
                                                        "q < 0.05", 
                                                        n.digits = 2), 
                   pct_associated = jag2::npct(associated_bci, 
                                               "Associated", 
                                               n.digits = 2))  -> mwas_summary


mwas_summary

# pivot summary wider   
mwas_summary_w <- pivot_wider(mwas_summary, 
                              id_cols = metabolite,
                              names_from = var.names.x, 
                              values_from = c(percent_significant_p05, 
                                              percent_significant_q05, 
                                              pct_associated)) 
mwas_summary_w 
# select(exposure, mixture, contains("SOL"), contains("CHS"), contains("Pooled"))  

# write.xlsx(mwas_summary_w,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/results_p_summary_final.xlsx")

###============================画火山图=======================================
####1.###可视化每个化合物有统计学意义的代谢物的火山图
# "Summary of Mixtures MWAS Results_hyper_g_v3.csv")

# -----------SOLAR Volcano Plots ----------------------------------------------
# mwas_results <- read_csv(
#   file = fs::path(dir_results_mixtures, 
#                   "all_pfas_mixtures_results_hyper_g.csv")) %>% 
#   as_tibble()

# Get reduced dataset   ####筛选限定值范围的差异话代谢物
###1.方式画火山图 

# results_df_only_p_value


####>>>>>>>>>>>>>读入代谢物匹配的代谢通路变量名称>>>>>>>>>>>>###################


colSums(is.na(results_df_only_p_value))
range(results_df_only_p_value$q_value)

mwas_reduced <- results_df_only_p_value %>% 
  dplyr::filter(q_value  < 1)
names(mwas_reduced)


####>>>>>>>将代谢物的名称进行合并
####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2

str(df_unique_2)   ###一共有378个确定的代谢物


df_unique_2$metabolite

mwas_reduced$metabolite



######>>>>>>删除ICC<0.3的代谢物  

df_icc <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc_20250222_number.xlsx")
####重复测量的32个同样样本的代谢物特征的一致性

df_icc %>% 
  filter(icc>=0.3) -> df_icc_1   ###筛选出icc>0.3的代谢物特征


df_icc_1$metabolites


df_unique_2 %>% 
  filter(metabolite %in% df_icc_1$metabolites) -> df_unique_3




mwas_reduced  %>% 
  left_join(df_unique_3,by="metabolite")  -> mwas_reduced_2

colSums(is.na(mwas_reduced_2))




# dev.off()

###2.方式画火山图

library(showtext)
library(sysfonts)
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

names(results_df_only_p_value)
results_df_only_p_value$var.names




###由于有些化合物可能没有统计学意义，会导致有些化合物标示
mwas_reduced_2 %>% 
  mutate(marker=case_when(
    q_value< 0.05 & estimate.x > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.05 & estimate.x< 0 ~ "Down" ,   ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  ))  -> res1

mwas_reduced_2 %>% 
  mutate(marker=case_when(
    q_value< 0.05 & estimate.x > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.05 & estimate.x< 0 ~ "Down" ,   ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>% 
  mutate(var.names.x=case_when(
    var.names.x=="PFOS.beta"~"PFOS",
    var.names.x=="PFOA.beta"~"PFOA",
    var.names.x=="PFDA.beta"~"PFDA",
    var.names.x=="L_PFHxS.beta"~"PFHxS",
    var.names.x=="PFUdA.beta"~"PFUdA",
    var.names.x=="PF3ONS_9CL.beta"~"6:2 Cl-PFESA"
  )) %>% 
  ggplot(aes(x =estimate.x, y = -log2(q_value),color= marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
  # xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线 
  geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
  # facet_wrap(~var.names, scales = "free") +
  # labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  # ggtitle("单组火山图") + #标题
  labs(x = "Estimate", y = "-log2(q_value)")+
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
  facet_wrap(~var.names.x, scales = "free") +
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) -> volcano_plot_only

volcano_plot_only




# font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/volcano_only_PFAS.pdf",width = 16, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_only #出图
par(opar)
dev.off()




####2.可视化混合化合物有统计学意义的代谢物的火山图

results_df_mix_p_value 

names(results_df_mix_p_value)

results_df_mix_p_value$p_value





mwas_reduced_mix <-results_df_mix_p_value %>% 
  dplyr::filter(q_value  < 1)
names(mwas_reduced_mix)
# Volcano Plot     ####画火山图
# # 
# ggplot(mwas_reduced,aes(x = estimate, y = -log10(p_value),color=significance)) +
#   geom_point(size = 1, alpha = 0.5) + 
#   geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
#   # geom_hline(yintercept = -log10(q_value), lty=4,col="black",lwd=0.8)+
#   facet_wrap(~var.names, scales = "free") +
#   ylab("-log10 p") +
#   xlab("Effect Estimate")

# dev.off() 
# 


####>>>>>>>匹配变量名称<<<<<<<<<<<<<<<<<<###########################
####>>>>>>>读入单一代谢化合物

df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2

str(df_unique_2)   ###一共有378个确定的代谢物


df_unique_2$metabolite

mwas_reduced$metabolite



######>>>>>>删除ICC<0.3的代谢物  

df_icc <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc_20250222_number.xlsx")
####重复测量的32个同样样本的代谢物特征的一致性

df_icc %>% 
  filter(icc>=0.3) -> df_icc_1   ###筛选出icc>0.3的代谢物特征


df_icc_1$metabolites



df_unique_2 %>% 
  filter(metabolite %in% df_icc_1$metabolites) -> df_unique_3





mwas_reduced_mix  %>% 
  left_join(df_unique_3,by="metabolite")  -> mwas_reduced_mix_2

# 
library(showtext)
library(sysfonts)
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置


names(mwas_reduced_mix_2)


volcano_plot_mix <-
  mwas_reduced_mix_2 %>% 
  mutate(marker=case_when(
    q_value < 0.05 & estimate > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.05 & estimate< 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate, y = -log(q_value,10),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, s llllllllize=2) +  #点的透明度、大小
  geom_point(aes(size = -log(q_value,10),color = marker),alpha = 0.7) +
  # ylim(c(0,75))+
  # xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线 
  # geom_hline(yintercept = 0, linetype = 2,size=1) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
  # facet_wrap(~var.names, scales = "free") +
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  # ggtitle("单组火山图") + #标题
  labs(x = "PFAS and Metabolism correlation", y = "-log10(q_value)")+
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


volcano_plot_mix





font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/volcano_all_PFAS.pdf",width = 12, height = 10, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_mix #出图

par(opar)
dev.off()





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



####>>>>>>>>>>>>>>>>>>>>二、读入混合线性模型的pfas_meta数据<<<<<<<<<<<<<<<#####



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






#将导出的数据重新读入数据  

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
  as_tibble() -> coefs_tran
# coefs_tran <- read_csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/coefs_pfas_meta.csv")





####将生成的代谢物传递给新的数据
coefs_tran ->coefs_tran_tidy_0

str(coefs_tran_tidy_0)

length(coefs_tran_tidy_0$metabolite)
table(coefs_tran_tidy_0$metabolite) 
coefs_tran_tidy_0$metabolite

coefs_tran_tidy_0[coefs_tran_tidy_0$metabolite %in% coefs_tran_tidy_0$metabolite[duplicated(coefs_tran_tidy_0$metabolite)], ]


coefs_tran_tidy_0$metabolite


###取出代谢物质单独的变量名称
coefs_tran_tidy_0 %>%
  group_by(metabolite) %>%
  dplyr::summarise(count=n())  -> coefs_tran_tidy_unique

coefs_tran_tidy_unique$metabolite

names(coefs_tran_tidy_unique)

str(coefs_tran_tidy_unique)

###将代谢物重新传递给coefs_tran_tidy

coefs_tran_tidy_0 ->coefs_tran_tidy 

coefs_tran_tidy

results_df <- coefs_tran_tidy  %>% 
  # select(-1) %>% 
  dplyr::mutate(p_value_estimate = na_if(coefs_tran_tidy$p.val, NA) %>% as.numeric)   ###去除p值中含有缺失的值

results_df



## Subset eta ------------------------------------
results_df_2 <-results_df %>%
  # dplyr::filter(str_detect(var.names,"eta")) %>%
  dplyr::rename(estimate = Mean,sd=SD,estimate_lcl=X2.5., estimate_ucl=X97.5.)
# mixtures_eta

# estimate, sd, lcl, ucl, p_value

results_df_2$var.names

## Pivot results wider ---------------------------------


# ## Pivot results wider ---------------------------------
# results_df_w <- results_df_2 %>% 
#   filter(var.names != "eta") %>%
#   tidylog::pivot_wider(id_cols = c(cohort, mixture, metabolite, exposure), 
#                        names_from = term, 
#                        values_from = c(estimate, sd, lcl, ucl, p_value)) %>% 
#   dplyr::select(cohort, mixture, metabolite, exposure, 
#                 contains("beta"), contains("pip")) %>% 
#   dplyr::select(-p_value_pip) %>% 
#   dplyr::rename(p_value = p_value_beta)
# 
# 


###--------转换BCI CI95-----------------------


results_df_3 <- results_df_2 %>% 
  dplyr::filter(var.names %in% c("eta.high","eta.low","psi"))  ###提取出eta.high、eta.low、psi三个重要变量

results_df_3

results_df_3 %>%  
  filter(var.names=="psi") -> results_df_4

paste0(round(results_df_4$estimate,2)," (",round(results_df_4$estimate_lcl,2),', ',round(results_df_4$estimate_ucl,2),")") %>% 
  as_tibble() %>% 
  dplyr::rename(BCI_CI95=value)  %>% 
  cbind(results_df_4)-> results_df_BCI_95CI  ###根据需要将想要的变量提取处理 

str(results_df_BCI_95CI)

str(results_df_2) 
results_df_2

results_df_2 %>% 
  dplyr::filter(!var.names %in% c("eta.high","eta.low","psi")) %>% 
  dplyr::select(metabolite,everything())->results_df_5   ###过滤掉不是eta.high、eta.low、psi三个重要变量


###将转换后的BCI_95与不含eta.high、eta.low、psi三个变量进行合并
results_df_5

results_df_5 %>% 
  left_join(results_df_BCI_95CI,by="metabolite") %>% 
  dplyr::filter(var.names.x %in% c("PFOS.beta","PFOA.beta","PFDA.beta","PFUdA.beta","L_PFHxS.beta","PF3ONS_9CL.beta"))-> results_df_6   ###合并后的数据  ###只保留beta变量,每个化合物的beta和混合物beta

results_df_6 
names(results_df_6)


###-------------只保留后验概率PIP数值-------------------------
results_df_5 %>% 
  dplyr::filter(var.names %in% c("PFOS.gamma","PFOA.gamma","PFDA.gamma","PFUdA.gamma","L_PFHxS.gamma","PF3ONS_9CL.gamma"))  %>% 
  tidylog::pivot_wider( id_cols = c(metabolite), 
                        names_from = var.names,
                        values_from = c(estimate)) -> results_df_7  ###保存每个化合物的后验概率

results_df_7   




# Calculate new p values


# Calculate new p values  计算单一物质的wald检验
names(results_df_6)

results_df_6 %>% 
  dplyr::mutate(wald = (estimate.x/sd.x), 
                # p = 2*(1-pnorm(abs(wald),0,1)),
                p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
                p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
                neg_log_p = -log10(p_value))  %>% 
  group_by(var.names.x) %>% 
  dplyr::mutate(q_value=p.adjust(p_value,method = "fdr"),
                significance = if_else(p_value < 0.05,
                                       "p < 0.05", 
                                       "Not Sig."), 
                significancefdr = if_else(q_value < 0.05,
                                          "q < 0.05",
                                          "Not Sig.")) %>% 
  ungroup() -> results_df_only_p_value

results_df_only_p_value

names(results_df_only_p_value)

write.xlsx(results_df_only_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_only_p_value_mix_liner.xlsx") 




# result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_pfas_meta





########---------------计算混合物质的wald检验---------------------
# Calculate new p values   
# 
results_df_BCI_95CI
names(results_df_BCI_95CI)

results_df_BCI_95CI %>% 
  dplyr::mutate(wald = (estimate/sd), 
                # p = 2*(1-pnorm(abs(wald),0,1)),
                p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
                p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
                neg_log_p = -log10(p_value)) %>% 
  # group_by(metabolite) %>% 
  dplyr::mutate(q_value = p.adjust(p_value,
                                   method = "fdr"), 
                significance=if_else(p_value < 0.05,
                                     "p < 0.05", 
                                     "Not Sig."), 
                significancefdr=if_else(q_value < 0.05,
                                        "q < 0.05",
                                        "Not Sig.")) %>% 
  ungroup() -> results_df_mix_p_value

results_df_mix_p_value

names(results_df_mix_p_value)

write.xlsx(results_df_mix_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_all_p_value_mix_liner.xlsx")


# 
# q_value = p.adjust(c(0.896,0.928,0.868),
#                    method = "fdr")



###--------------将BCI_CI95与后验概率PIP进行合并-------------


results_df_mix_p_value %>% 
  left_join(results_df_7 ,by="metabolite") %>% 
  dplyr::select(metabolite,BCI_CI95,p_value,q_value,PFOS.gamma:PF3ONS_9CL.gamma,everything())-> results_df_final   ###得到复现文章的Supplemental Excel file文件的结果

write.xlsx(results_df_final,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/results_df_final_mix_line.xlsx")


# install.packages("remotes") 

# remotes::install_github("JAGoodrich/jag2")


names(results_df_only_p_value)


results_df_only_p_value %>% 
  mutate(associated_bci = if_else(estimate_lcl.x > 0 | estimate_ucl.x < 0,
                                  "Associated", 
                                  "Not Associated")) %>% 
  dplyr::group_by(metabolite, var.names.x) %>% 
  dplyr::summarise(n_features = length(metabolite), 
                   percent_significant_p05 = jag2::npct(significance, 
                                                        "p < 0.05", 
                                                        n.digits = 2), 
                   percent_significant_q05 = jag2::npct(significancefdr, 
                                                        "q < 0.05", 
                                                        n.digits = 2), 
                   pct_associated = jag2::npct(associated_bci, 
                                               "Associated", 
                                               n.digits = 2))  -> mwas_summary

mwas_summary


# pivot summary wider   
mwas_summary_w <- pivot_wider(mwas_summary, 
                              id_cols = metabolite,
                              names_from = var.names.x, 
                              values_from = c(percent_significant_p05, 
                                              percent_significant_q05, 
                                              pct_associated)) 

mwas_summary_w 
# select(exposure, mixture, contains("SOL"), contains("CHS"), contains("Pooled"))  

# write.xlsx(mwas_summary_w,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/PFAS_meta_mix/results_p_summary_final.xlsx")

###============================画火山图=======================================
####1.###可视化每个化合物有统计学意义的代谢物的火山图
# "Summary of Mixtures MWAS Results_hyper_g_v3.csv")

# -----------SOLAR Volcano Plots ----------------------------------------------
# mwas_results <- read_csv(
#   file = fs::path(dir_results_mixtures, 
#                   "all_pfas_mixtures_results_hyper_g.csv")) %>% 
#   as_tibble()

# Get reduced dataset   ####筛选限定值范围的差异话代谢物
###1.方式画火山图 

# results_df_only_p_value


####>>>>>>>>>>>>>读入代谢物匹配的代谢通路变量名称>>>>>>>>>>>>###################


colSums(is.na(results_df_only_p_value))
range(results_df_only_p_value$q_value)

mwas_reduced <- results_df_only_p_value %>% 
  dplyr::filter(q_value  < 1)
names(mwas_reduced)


####>>>>>>>将代谢物的名称进行合并
####>>>>>>>读入单一代谢化合物
df_unique <- read_excel("raw_data/Matebolomics_sperm/20241213_haixia/汇总_LGM_20241231.xlsx",sheet = "Sheet1")

names(df_unique)
df_unique %>% 
  dplyr::rename(metabolite="META")  -> df_unique_2

str(df_unique_2)   ###一共有378个确定的代谢物


df_unique_2$metabolite

mwas_reduced$metabolite



######>>>>>>删除ICC<0.3的代谢物  

df_icc <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/ICC_lappy/data_icc_20250222_number.xlsx")
####重复测量的32个同样样本的代谢物特征的一致性

df_icc %>% 
  filter(icc>=0.3) -> df_icc_1   ###筛选出icc>0.3的代谢物特征


df_icc_1$metabolites



df_unique_2 %>% 
  filter(metabolite %in% df_icc_1$metabolites) -> df_unique_3



mwas_reduced  %>% 
  left_join(df_unique_3,by="metabolite")  -> mwas_reduced_2

colSums(is.na(mwas_reduced_2))




# dev.off()

###2.方式画火山图

library(showtext)
library(sysfonts)
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

names(results_df_only_p_value)
results_df_only_p_value$var.names




###由于有些化合物可能没有统计学意义，会导致有些化合物标示
mwas_reduced_2 %>% 
  mutate(marker=case_when(
    q_value< 0.05 & estimate.x > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.05 & estimate.x< 0 ~ "Down" ,   ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  ))  -> res1

mwas_reduced_2 %>% 
  mutate(marker=case_when(
    q_value< 0.05 & estimate.x > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.05 & estimate.x< 0 ~ "Down" ,   ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>% 
  mutate(var.names.x=case_when(
    var.names.x=="PFOS.beta"~"PFOS",
    var.names.x=="PFOA.beta"~"PFOA",
    var.names.x=="PFDA.beta"~"PFDA",
    var.names.x=="L_PFHxS.beta"~"PFHxS",
    var.names.x=="PFUdA.beta"~"PFUdA",
    var.names.x=="PF3ONS_9CL.beta"~"6:2 Cl-PFESA"
  )) %>% 
  ggplot(aes(x =estimate.x, y = -log2(q_value),color= marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
  # xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线 
  geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
  # facet_wrap(~var.names, scales = "free") +
  # labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  # ggtitle("单组火山图") + #标题
  labs(x = "Estimate", y = "-log2(q_value)")+
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
  facet_wrap(~var.names.x, scales = "free") +
  geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) -> volcano_plot_only

volcano_plot_only




# font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/volcano_only_PFAS_mix_liner.pdf",width = 16, height = 9, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_only #出图
par(opar)
dev.off()



####2.可视化混合化合物有统计学意义的代谢物的火山图

results_df_mix_p_value

names(results_df_mix_p_value)

results_df_mix_p_value$p_value





mwas_reduced_mix <-results_df_mix_p_value %>% 
  dplyr::filter(q_value  < 1)
names(mwas_reduced_mix)
# Volcano Plot     ####画火山图
# # 
# ggplot(mwas_reduced,aes(x = estimate, y = -log10(p_value),color=significance)) +
#   geom_point(size = 1, alpha = 0.5) + 
#   geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
#   # geom_hline(yintercept = -log10(q_value), lty=4,col="black",lwd=0.8)+
#   facet_wrap(~var.names, scales = "free") +
#   ylab("-log10 p") +
#   xlab("Effect Estimate")

# dev.off() 
# 


####>>>>>>>匹配变量名称<<<<<<<<<<<<<<<<<<###########################
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

mwas_reduced$metabolite


mwas_reduced_mix  %>% 
  left_join(df_unique_3,by="metabolite")  -> mwas_reduced_mix_2

# 
library(showtext)
library(sysfonts)
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置


names(mwas_reduced_mix_2)


volcano_plot_mix <-
  mwas_reduced_mix_2 %>% 
  mutate(marker=case_when(
    q_value < 0.05 & estimate > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
    q_value< 0.05 & estimate< 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
    TRUE ~ "No"
  )) %>%
  ggplot(aes(x =estimate, y = -log(q_value,10),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
  # geom_point(alpha=0.65, s llllllllize=2) +  #点的透明度、大小
  geom_point(aes(size = -log(q_value,10),color = marker),alpha = 0.7) +
  # ylim(c(0,75))+
  # xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围
  scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
  geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线 
  # geom_hline(yintercept = 0, linetype = 2,size=1) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2,size = 1) +
  # facet_wrap(~var.names, scales = "free") +
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  # ggtitle("单组火山图") + #标题
  labs(x = "PFAS and Metabolism correlation", y = "-log10(q_value)")+
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


volcano_plot_mix



font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/pfas_meta_result/tidy_result/volcano_all_PFAS_mix_liner.pdf",width = 12, height = 10, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
volcano_plot_mix #出图

par(opar)
dev.off()







 
# #################################STL中的差异代谢物############################## 
# 
# load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_STL_meta.Rdata")
# STL_coefs
# 
# library(tidyverse)
# 
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
# #将导出的数据重新读入数据  
# 
# coefs_tran <- read_csv("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_STL_meta.csv")
# 
# coefs_tran
# 
# #
# coefs_tran %>% 
#   dplyr::select(-1) %>% 
#   as_tibble()->coefs_tran_tidy 
# coefs_tran_tidy 
# 
# 
# # 
# # # 定义名称
# # names_first <- c("var.names", "Mean", "SD", "X2.5.", "X97.5.", "p.val", "metabolite")
# # 
# # # 创建一个长度为794的向量，字母重复循环
# # names_new <- rep(names_first, length.out = 14)
# # 
# # # 使用 setNames 函数批量重命名
# # setNames(coefs_tran_tidy, names_new)-> coefs_tran_tidy_newnames
# # 
# # coefs_tran_tidy_newnames
# # 
# # # 将指定列转换为长格式
# # 
# # 
# # # 确定数据框列的总数
# # num_cols <- ncol(coefs_tran_tidy_newnames)
# # num_cols
# # # 每7列拆分
# # num_splits <- ceiling(num_cols/ 7)
# # 
# # # 将数据框按每7列拆分为多个数据框
# # lapply(seq_len(num_splits), function(i) {
# #   start_col <- (i - 1) * 7 + 1
# #   end_col <- min(i * 7, num_cols)
# #   coefs_tran_tidy_newnames %>%
# #     dplyr::select(start_col:end_col)
# # })  -> split_dfs 
# # 
# # do.call(rbind,split_dfs) -> combind_df_1   ###将list值转换为数据框
# # 
# # str(combind_df_1)
# # 
# # combind_df_1 %>% 
# #   dplyr::filter(!var.names=="NA")->temp_final 
# # 
# # str(temp_final)
# # write.csv(temp_final,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/temp_final.csv")
# 
# 
# 
# # gsub("^var.names", "terms", colnames(coefs_tran_tidy)) -> colnames(coefs_tran_tidy)
# 
# # coefs_tran_tidy$var.names[coefs_tran_tidy$var.names=="^var.names"]  <- "term"
# # 
# # coefs_tran_tidy
# # 
# # names(coefs)
# # # 保存修改后的 Excel 文件
# # saveWorkbook(coefs,"D:/Postdoctoral_SCI/Sperm_bank/R/20240407_PFOA_DNA_three_article/PFAS_metabolomics_EHP_2022-master/PFAS_metabolomics_EHP_2022-master/5_reproducible_example/result_pfas_2.xlsx", overwrite = TRUE)
# # 
# # 
# # ###然后再将生成的list重新保存为独立的excel文件
# # lapply(seq_along(res1), function(i) {
# #   write.xlsx(x = res1[[i]], file = paste0("result/liu_nian_202400606/20240623_result/20240720_Metal_follow/result_Metal_(cell_sperm)_follow_adj_20240720_",names(res1) %>% as_tibble() %>% slice(i), ".xlsx"), row.names = FALSE)
# # })
# # 
# # # write.xlsx(simulated_data,"D:/Postdoctoral_SCI/Sperm_bank/R/20240407_PFOA_DNA_three_article/PFAS_metabolomics_EHP_2022-master/PFAS_metabolomics_EHP_2022-master/5_reproducible_example/raw_data.xlsx")
# # # message(paste("SUCCESS from rank", comm.rank()))
# # # 
# # # finalize()
# # # 
# # # slurmR::Slurm_clean(coefs)
# # 
# 
# 
# 
# #######################读入数据整理#########################
# library(purrr)
# library(tidyverse)
# 
# # temp <- read.csv("D:/Postdoctoral_SCI/Sperm_bank/R/20240407_PFOA_DNA_three_article/PFAS_metabolomics_EHP_2022-master/PFAS_metabolomics_EHP_2022-master/5_reproducible_example/result_pfas_rows.csv")
# # 
# # temp
# # purrr::map(read_data_from_hpc(temp_1,n_col_in_df = 7))
# # # Read in all results --------------------------------------
# library(tidyverse)
# library(kableExtra)
# library(broom)
# library(cowplot)
# library(here)
# library(fs)
# 
# # results <- maps::map(temp,
# #                ~read_data_from_hpc(., n_col_in_df = 7))
# 
# ?map
# 
# ?read_data_from_hpc
# ## Subset eta ------------------------------------
# # mixtures_eta <- filter(temp_1,
# #                        term == "eta") %>% 
# #   rename(feature = metabolite)
# 
# # 
# # cohort_name <- data.frame("CHS", "SOL", "Pooled")
# # 
# # mixture_name <- data.frame("pfsas",  "pfcas",  "all pfas")
# # 
# # cohort_mixture = str_c(cohort_name, mixture_name, sep = "_")
# # cohort_mixture
# # 
# # str_split_fixed(cohort_mixture, "_", n = 2)[,1]
# # str_split_fixed(cohort_mixture, "_", n = 2)[,2]
# 
# # results_df <- bind_rows(temp_1, .id = "cohort_mixture") %>% 
# #   as_tibble() %>% 
# #   mutate(cohort = str_split_fixed(cohort_mixture, "_", n = 2)[,1], 
# #          mixture = str_split_fixed(cohort_mixture, "_", n = 2)[,2]) %>% 
# #   select(cohort, mixture, everything(), -cohort_mixture)
# 
# 
# 
# 
# 
# results_df <- coefs_tran_tidy  %>% 
#   # select(-1) %>% 
#   dplyr::mutate(p_value_estimate = na_if(coefs_tran_tidy$p.val, NA) %>% as.numeric)
# 
# results_df
# 
# 
# 
# results_df_2 <-results_df %>%
#   # dplyr::filter(str_detect(var.names,"eta")) %>%
#   dplyr::rename(estimate = Mean,sd=SD,estimate_lcl=X2.5., estimate_ucl=X97.5.)
# # mixtures_eta
# 
# # estimate, sd, lcl, ucl, p_value
# 
# results_df_2$var.names
# 
# ## Pivot results wider ---------------------------------
# 
# ###--------转换BCI CI95-----------------------
# 
# 
# results_df_3 <- results_df_2 %>% 
#   dplyr::filter(var.names %in% c("eta.high","eta.low","psi"))  ###提取出eta.high、eta.low、psi三个重要变量
# 
# results_df_3
# 
# results_df_3 %>%  
#   filter(var.names=="psi") -> results_df_4
# 
# paste0(round(results_df_4$estimate,2)," (",round(results_df_4$estimate_lcl,2),', ',round(results_df_4$estimate_ucl,2),")") %>% 
#   as_tibble() %>% 
#   dplyr::rename(BCI_CI95=value)  %>% 
#   cbind(results_df_4)-> results_df_BCI_95CI  ###根据需要将想要的变量提取处理 
# 
# str(results_df_BCI_95CI)
# 
# str(results_df_2) 
# results_df_2
# 
# results_df_2 %>% 
#   dplyr::filter(!var.names %in% c("eta.high","eta.low","psi")) %>% 
#   dplyr::select(metabolite,everything())->results_df_5   ###过滤掉不是eta.high、eta.low、psi三个重要变量
# 
# 
# ###将转换后的BCI_95与不含eta.high、eta.low、psi三个变量进行合并
# results_df_5
# 
# results_df_5 %>% 
#   left_join(results_df_BCI_95CI,by="metabolite") %>% 
#   dplyr::filter(var.names.x %in% c("STL.beta"))-> results_df_6   ###合并后的数据  ###只保留beta变量
# 
# results_df_6 
# 
# ###-------------只保留后验概率PIP数值-------------------------
# results_df_5 %>% 
#   dplyr::filter(var.names %in% c("STL.gamma"))  %>% 
#   tidylog::pivot_wider( id_cols = c( metabolite), 
#                         names_from = var.names,
#                         values_from = c(estimate)) -> results_df_7  ###保存每个化合物的后验概率
# results_df_7   
# 
# names(results_df_7)
# 
# 
# 
# 
# 
# 
# # Calculate new p values
# 
# names(results_df_6)
# # Calculate new p values  计算单一物质的wald检验
# results_df_6 %>% 
#   dplyr::mutate(wald = (estimate.x/sd.x), 
#                 # p = 2*(1-pnorm(abs(wald),0,1)),
#                 p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
#                 p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
#                 neg_log_p = -log10(p_value),q_value=p.adjust(p_value,
#                                                              method = "fdr")) %>% 
#   group_by(metabolite,var.names.x) %>% 
#   dplyr::mutate(significance = if_else(p_value < 0.05,
#                                        "p < 0.05", 
#                                        "Not Sig."), 
#                 significancefdr = if_else(q_value < 0.05,
#                                           "q < 0.05",
#                                           "Not Sig.")) %>% 
#   ungroup() -> results_df_only_p_value
# 
# results_df_only_p_value
# 
# write.xlsx(results_df_only_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/STL_meta_result/results_df_only_p_value.xlsx") 
# 
# 
# 
# 
# # result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/coefs_pfas_meta
# 
# 
# 
# 
# 
# ########---------------计算混合物质的wald检验---------------------
# # Calculate new p values   
# # 
# results_df_BCI_95CI
# 
# 
# results_df_BCI_95CI %>% 
#   dplyr::mutate(wald = (estimate/sd), 
#                 # p = 2*(1-pnorm(abs(wald),0,1)),
#                 p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
#                 p_value = if_else(wald^2 > 2000, 4.6*(10^(-256)), p_value),
#                 neg_log_p = -log10(p_value),
#                 q_value=p.adjust(p_value,method = "fdr")) %>% 
#   group_by(metabolite) %>% 
#   dplyr::mutate( 
#     significance=if_else(p_value < 0.05,
#                          "p < 0.05", 
#                          "Not Sig."), 
#     significancefdr=if_else(q_value < 0.05,
#                             "q < 0.05",
#                             "Not Sig.")) %>% 
#   ungroup() -> results_df_mix_p_value
# 
# results_df_mix_p_value
# names(results_df_mix_p_value)
# 
# write.xlsx(results_df_mix_p_value,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/STL_meta_result/results_df_mix_p_value.xlsx")
# 
# 
# # 
# # q_value = p.adjust(c(0.896,0.928,0.868),
# #                    method = "fdr")
# 
# 
# 
# ###--------------将BCI_CI95与后验概率PIP进行合并-------------
# 
# 
# results_df_mix_p_value %>% 
#   left_join(results_df_7 ,by="metabolite") %>% 
#   dplyr::select(metabolite,BCI_CI95,p_value,q_value,STL.gamma,everything())-> results_df_final   ###得到复现文章的Supplemental Excel file文件的结果
# 
# write.xlsx(results_df_final,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/STL_meta_result/results_df_final.xlsx")
# 
# 
# 
# 
# 
# names(results_df_only_p_value)
# 
# 
# 
# 
# results_df_only_p_value %>% 
#   mutate(associated_bci = if_else(estimate_lcl.x > 0 | estimate_ucl.x < 0,
#                                   "Associated", 
#                                   "Not Associated")) %>% 
#   group_by(metabolite, var.names.x) %>% 
#   dplyr::summarise(n_features = length(metabolite), 
#                    percent_significant_p05 = jag2::npct(significance, 
#                                                         "p < 0.05", 
#                                                         n.digits = 2), 
#                    percent_significant_q05 = jag2::npct(significancefdr, 
#                                                         "q < 0.05", 
#                                                         n.digits = 2), 
#                    pct_associated = jag2::npct(associated_bci, 
#                                                "Associated", 
#                                                n.digits = 2))  -> mwas_summary
# 
# 
# mwas_summary
# 
# # pivot summary wider   
# mwas_summary_w <- pivot_wider(mwas_summary, 
#                               id_cols = metabolite,
#                               names_from = var.names.x, 
#                               values_from = c(percent_significant_p05, 
#                                               percent_significant_q05, 
#                                               pct_associated)) 
# mwas_summary_w 
# # select(exposure, mixture, contains("SOL"), contains("CHS"), contains("Pooled"))  
# 
# write.xlsx(mwas_summary_w,"result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/STL_meta_result/results_p_summary_final.xlsx")
# 
# ###============================画火山图=======================================
# ####1.###可视化每个化合物有统计学意义的代谢物的火山图
# # "Summary of Mixtures MWAS Results_hyper_g_v3.csv")
# 
# # -----------SOLAR Volcano Plots ----------------------------------------------
# # mwas_results <- read_csv(
# #   file = fs::path(dir_results_mixtures, 
# #                   "all_pfas_mixtures_results_hyper_g.csv")) %>% 
# #   as_tibble()
# 
# # Get reduced dataset   ####筛选限定值范围的差异话代谢物
# ###1.方式画火山图
# colSums(is.na(results_df_only_p_value))
# range(results_df_only_p_value$q_value)
# 
# mwas_reduced <- results_df_only_p_value %>% 
#   dplyr::filter(q_value  < 1)
# names(mwas_reduced)
# 
# 
# 
# ####>>>>>>>
# 
# df_kegg<- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/240823_最终确定版data_lgm.xlsx",sheet = "合并")
# str(df_kegg) 
# 
# 
# 
# df_kegg %>% 
#   dplyr::rename(metabolite="META") -> df_kegg_1
# 
# str(df_kegg_1)
# 
# colSums(is.na(df_kegg_1))
# 
# 
# mwas_reduced  %>% 
#   left_join(df_kegg_1,by="metabolite")  -> mwas_reduced_2
# 
# 
# mwas_reduced_2
# 
# names(mwas_reduced_2)
# 
# 
# 
# 
# 
# # Volcano Plot     ####画火山图
# # # 
# # ggplot(mwas_reduced,aes(x = estimate, y = -log10(p_value),color=significance)) +
# #   geom_point(size = 1, alpha = 0.5) + 
# #   geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
# #   # geom_hline(yintercept = -log10(q_value), lty=4,col="black",lwd=0.8)+
# #   facet_wrap(~var.names, scales = "free") +
# #   ylab("-log10 p") +
# #   xlab("Effect Estimate")
# 
# # dev.off() 
# mwas_reduced_2 %>% 
#   mutate(sig=case_when(
#     q_value < 0.05 & estimate.x > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
#     q_value < 0.05 & estimate.x< 0 ~ "Down" ,   ###p小于0.05，且相关系数小于0代表下调
#     TRUE ~ "No"
#   )) -> mwas_reduced_sig 
# 
# mwas_reduced_sig$sig
# ###绘图——基础火山------fdr校正模型
# ###
# 
# p1 <- ggplot(mwas_reduced_sig, aes(x =estimate.x, y = -log10(q_value),color=sig)) + #x、y轴取值限制，颜色根据"Sig"
#   # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
#   geom_point(aes(size = -log(q_value, 10),color = sig),alpha = 0.7) +
#   scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
#   # xlim(-log(max(abs(mwas_reduced_sig$estimate)), 10),log(max(abs(mwas_reduced_sig$estimate), 10))) + 
#   xlim(c(-30, 30)) +  #调整点的颜色和x轴的取值范围
#   #  #调整点的颜色和x轴的取值范,使图像更好看
#   geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
#   # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线 
#   geom_hline(yintercept = 0, linetype = 2,size=1) +
#   # facet_wrap(~var.names, scales = "free") +
#   labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
#   # ggtitle("单组火山图") + #标题
#   labs(x = "Estimate", y = "-log10(q_value)")+
#   theme_bw() +
#   theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
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
#   facet_wrap(~var.names.x, scales = "free") +
#   geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) 
#   
# p1
# 
# 
# ###2.方式画火山图
# 
# library(showtext)
# library(sysfonts)
# font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
# 
# names(results_df_only_p_value)
# 
# 
# mwas_reduced_2 %>% 
#   mutate(marker=case_when(
#     q_value< 0.05 & estimate.x > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
#     q_value< 0.05 & estimate.x< 0 ~ "Down" ,   ###p小于0.05，且相关系数小于0代表下调
#     TRUE ~ "No"
#   ))  %>% 
#   mutate(var.names.x=case_when(
#     var.names.x=="STL.beta"~"STL"
#   )) %>% 
#   ggplot(aes(x =estimate.x, y = -log10(q_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
#   # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
#   geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
#   xlim(c(-30, 30)) +  #调整点的颜色和x轴的取值范围
#   scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
#   geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
#   # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线 
#   geom_hline(yintercept = 0, linetype = 2,size=1) +
#   # facet_wrap(~var.names, scales = "free") +
#   labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
#   # ggtitle("单组火山图") + #标题
#   labs(x = "Estimate", y = "-log10(q_value)")+
#   theme_bw() +
#   theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
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
#   facet_wrap(~var.names.x, scales = "free")+
#   geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE)   ->volcano_plot_only
# 
# volcano_plot_only
# 
# 
# # font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
# cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/STL_meta_result/volcano_only_STL.pdf",width = 16, height = 9, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# volcano_plot_only #出图
# 
# par(opar)
# dev.off()
# 
# 
# 
# ####2.可视化混合化合物有统计学意义的代谢物的火山图
# 
# results_df_mix_p_value
# 
# names(results_df_mix_p_value)
# 
# results_df_mix_p_value$p_value
# 
# 
# mwas_reduced_mix <-results_df_mix_p_value %>% 
#   dplyr::filter(q_value  < 1)
# names(mwas_reduced_mix)
# 
# 
# ####>>>>>>>
# 
# df_kegg<- read_excel("raw_data/Matebolomics_sperm/20240823_haixia/240823_最终确定版data_lgm.xlsx",sheet = "合并")
# str(df_kegg) 
# 
# 
# 
# df_kegg %>% 
#   dplyr::rename(metabolite="META") -> df_kegg_1
# 
# str(df_kegg_1)
# 
# colSums(is.na(df_kegg_1))
# 
# 
# mwas_reduced_mix %>% 
#   left_join(df_kegg_1,by="metabolite")  -> mwas_reduced_mix_2
# 
# 
# mwas_reduced_2
# 
# names(mwas_reduced_2)
# 
# 
# 
# 
# 
# # Volcano Plot     ####画火山图
# # # 
# # ggplot(mwas_reduced,aes(x = estimate, y = -log10(p_value),color=significance)) +
# #   geom_point(size = 1, alpha = 0.5) + 
# #   geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
# #   # geom_hline(yintercept = -log10(q_value), lty=4,col="black",lwd=0.8)+
# #   facet_wrap(~var.names, scales = "free") +
# #   ylab("-log10 p") +
# #   xlab("Effect Estimate")
# 
# # dev.off() 
# mwas_reduced_mix_2 %>% 
#   mutate(sig=case_when(
#     estimate<0~0,
#     estimate>0~1
#   )) -> mwas_reduced_mix_sig 
# 
# names(mwas_reduced_mix_sig)
# 
# mwas_reduced_mix_sig$sig
# 
# 
# 
# ggplot(mwas_reduced_mix_sig, aes(x =estimate, y = -log10(q_value),color=significancefdr)) + #x、y轴取值限制，颜色根据"Sig"
#   geom_point(alpha=0.65, size=2) +  #点的透明度、大小
#   scale_color_manual(values=c("#546de5", "#ff4757")) +
#   xlim(c(-30, 30)) +  #调整点的颜色和x轴的取值范围
#   geom_vline(xintercept = 0, color = "#ff4757", linetype = 2)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
#   geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线
#   # facet_wrap(~var.names, scales = "free") +
#   labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
#   # ggtitle("单组火山图") + #标题
#   theme_bw() + # 主题，help(theme)查找其他个性化设置
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.position="right", 
#         legend.title = element_blank()
#   ) +
#   geom_text(aes(label=name), hjust = 0.5, vjust = -1.5, size = 3, check_overlap = TRUE,show.legend = FALSE) 
# 
# results_df_mix_p_value$p_value
# 
# source("R/100-tools.R") 
# 
# ###查看调色板颜色
# library("scales")
# show_col(ggsci::pal_futurama("planetexpress")(12))  ###查看调色板颜色
# 
# 
# marker_color <-
#   c("Up" = ggsci::pal_futurama()(n = 9)[6],
#     "Down" = ggsci::pal_futurama()(n = 9)[3],
#     "No" = ggsci::pal_futurama()(n = 9)[9]
#   )
# 
# 
# 
# # scale_color_manual(values=c("#377eb8", "#e41a1c"))
# 
# # marker_color <-
# #   c(
# #     "Up" = c("#e41a1c"),
# #     "Down" = c("#0055aa"),
# #     "No" = c("#999999")
# #   )
# 
# ?ggsci::pal_futurama()
# 
# # install.packages("scales")
# 
# 
# # install.packages("ggsci")
# # pal_futurama(palette = c("planetexpress"), alpha = 1)
# 
# # dev.off()
# # 
# library(showtext)
# library(sysfonts)
# font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
# 
# volcano_plot_mix <-
#   mwas_reduced_mix_2  %>% 
#   mutate(marker=case_when(
#     q_value < 0.05 & estimate > 0 ~ "Up",   ###p小于0.05，且相关系数大于0代表上调
#     q_value< 0.05 & estimate< 0 ~ "Down",    ###p小于0.05，且相关系数小于0代表下调
#     TRUE ~ "No"
#   )) %>%
#   ggplot(aes(x =estimate, y = -log10(q_value),color=marker)) + #x、y轴取值限制，颜色根据"Sig"
#   # geom_point(alpha=0.65, size=2) +  #点的透明度、大小
#   geom_point(aes(size = -log(q_value, 10),color = marker),alpha = 0.7) +
#   ylim(c(0,75))+
#   xlim(c(-30, 30)) +  #调整点的颜色和x轴的取值范围
#   scale_color_manual(values=c("#377eb8","#525252","#e41a1c")) + 
#   geom_vline(xintercept = 0, color = "#e41a1c", linetype = 2,size=1)+#添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
#   # geom_hline(yintercept = 3, lty=4,col="black",lwd=0.8)+  #添加y轴辅助线 
#   geom_hline(yintercept = 0, linetype = 2,size=1) +
#   # facet_wrap(~var.names, scales = "free") +
#   labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
#   # ggtitle("单组火山图") + #标题
#   labs(x = "Estimate", y = "-log10(q_value)")+
#   theme_bw() +
#   theme(legend.text=element_text(size =14,colour = "black",family = "Arial",face = "bold"),
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
# 
# volcano_plot_mix
# 
# 
# 
# 
# 
# font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置
# cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_STL_meta/STL_meta_result/volcano_mix_STL.pdf",width = 16, height = 9, pointsize = 20,family = "Arial")
# opar <- par(no.readonly = TRUE)
# par(pin = c(9,16), mar = c(7,7.5,1,1.5))  ##,mgp=c(3,2,0)
# volcano_plot_mix #出图
# 
# par(opar)
# dev.off()
# 
