
suppressMessages(library(readxl, quietly = TRUE))
suppressMessages(library(R2jags, quietly = TRUE))
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(rjags, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(pbdMPI, quietly = TRUE))
suppressMessages(library(lme4, quietly = TRUE))
suppressMessages(library(brms, quietly = TRUE))




###====================读入一般资料和THM的信息===============================


########>>>>>>>>>>>>一、混合线性回归分析中<<<<<<<<<<<<<<<<<<<<<<<##############


df_background_pfas_stl_zscore_final_tidy <- read_excel("/dssg/home/acct-yxwang0203/yxwang0203/work_data/pfas_STL_meta/df_background_pfas_stl_zscore_final_tidy.xlsx")

# 
# df_background_pfas_stl_zscore_final_tidy <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/df_background_pfas_stl_zscore_final_tidy.xlsx")
###利用取对数和归一化之后的端粒和代谢物，进行混合线性回归分析
df_background_pfas_stl_zscore_final_tidy    ###混合线性回归分析

df_background_pfas_stl_zscore_final_tidy$mtDNAcn

df_background_pfas_stl_zscore_final_tidy$STL

names(df_background_pfas_stl_zscore_final_tidy)




###3.进行贝叶斯和蒙特卡洛回归分析
library(brms)

# 使用 lapply 对第 25 到 26 列进行循环

lapply(names(df_background_pfas_stl_zscore_final_tidy [,c(25:26)]), function(col_name) {
  
  # 定义动态模型公式，插入当前列名
  formula_str <- paste0("STL ~ ", col_name, " + age_class_2 + BMI_class_2 + education + marriage_class_2 + income + smk_class_2 + drk_class_2 + season + abstinence_class_2+(1|number)")
  
  # 加载 brms 包
  library(brms)
  
  # 运行贝叶斯混合效应模型
  fit01 <- brms::brm(
    formula = as.formula(formula_str),
    data = df_background_pfas_stl_zscore_final_tidy,
    family = gaussian,
    chains = 4,
    iter = 500,
    warmup = 100,
    control = list(adapt_delta = 0.95),
    backend = "rstan"
  )
  
  # 获取模型摘要
  summary(fit01) -> fit1
  fit1$fixed -> result
  
  result %>% 
    mutate(p_value=case_when(
      `l-95% CI` < 0 & `u-95% CI` > 0 ~ ">0.05",  # 置信区间跨越零，p值大于0.05
      TRUE ~ "<0.05"  # 否则 p值小于0.05
    )) -> result
  # 返回模型结果
  return(result)
  
}) -> res1  ###将循环中每个y生成的list结果保存为res1



res1 -> res5



###然后再将生成的list重新保存为独立的excel文件
###给每个list添加行变量名称
lapply(seq_along(res5), function(i) {
  res5[[i]] %>% cbind(rownames(res5[[i]]))
}) -> res6


save(res6, file = 'coefs_meta_stl_mix_liner_bayes.Rdata')
