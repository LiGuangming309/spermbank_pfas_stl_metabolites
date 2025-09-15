##构建object文件 
###1.1首先构建expression_data 
###
names(df_background_DNA_LT_meta_mean_final)   ###精液质量和代谢组学数据已经进行取对数和标化
str(df_background_DNA_LT_meta_mean_final)

df_background_DNA_LT_meta_mean_final[,c(1,25:381)] %>% 
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
variable_info <- rownames(expression_data) %>% as.data.frame() %>% 
  dplyr::rename(variable_id=".") 


object <- 
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

object


##find marker which are change according to aging

###linear regression
library(tidyverse)
library(ggpubr)
library(rstatix)


expression_data <-
  extract_expression_data(object)  %>% 
  # `+`(1) %>%
  # log(2) %>%
  # apply(1, function(x) {
  #   (x - mean(x)) / sd(x)
  # }) %>%
  # t() %>%
  as.data.frame()

library(plyr)

sample_info <-
  object@sample_info

expression_data


colSums(is.na(expression_data))


lm(M33 ~ age_class_2 + BMI_class_2 + education + 
     marriage_class_2 + income + smk_class_2 + drk_class_2 +
     season + abstinence_class_2 , data = df_background_DNA_LT_meta_mean_final) %>%
  residuals()


class(expression_data)

library(tidyverse)
library(tidymass)
str(sample_info)

expression_data %>% as.data.frame() %>% rownames()

library(plyr)
# as.numeric(expression_data[M33,])

# expression_data[1,]
# length(rownames(expression_data))
###探究协变量对代谢物质的影响

lm_adjust <- function(expression_data,
                      sample_info,
                      threads = 5) {
  library(future)
  library(furrr)
  # plan(strategy = multisession(workers = threads))
  lapply(seq_along(rownames(expression_data)), function(i){
    # cat(name, " ")
    x = as.numeric(expression_data[i,])
    temp_data =
      data.frame(x = x, sample_info)
    adjusted_x <-
      lm(x ~ age_class_2+BMI_class_2+education+marriage_class_2+income+smk_class_2+drk_class_2+season+abstinence_class_2, data = temp_data) %>%
      residuals()
    adjusted_x
  }) %>%
    bind_rows() %>%
    as.data.frame() ->  new_expression_data
  
  colnames(new_expression_data) <-
    colnames(expression_data)
  
  rownames(new_expression_data) <-
    rownames(expression_data)
  new_expression_data
}


#######adjust BMI, sex, and IRIS, ethnicity
expression_data_2 <-
  lm_adjust(expression_data = expression_data,
            sample_info = object@sample_info,
            threads = 16)

temp_object <- object

temp_object@expression_data <- expression_data_2

temp_object@expression_data

save(temp_object, file = "result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/meta_cor_residuals.RData")
# save(p_fc, file = "p_fc")

temp_object@sample_info

names(temp_object@sample_info)

load("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/meta_cor_residuals.RData")

temp_object
# ##correlation
# ##cor
# 


as.numeric(expression_data[1,])

length(as.numeric(expression_data[1,]))

length(temp_object@sample_info[,c(24)])

cor.test(as.numeric(expression_data[1,]), temp_object@sample_info[,c(24)], method = "spearman")


####>>>>>>>>>>>>>>2.对精子端粒进行spearman回归分析<<<<<<<<<<<<<<<<<<<<<<<
cor_data <-
  lapply(seq_along(rownames(temp_object)), function(i){
    value <- as.numeric(expression_data[i,]) 
    cor_result <- cor.test(value,temp_object@sample_info[,c(24)], method = "spearman")
    
    data.frame(variable_id = rownames(temp_object)[i],
               cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()



names(temp_object@sample_info)



cor_data %>% 
  write.xlsx("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/meta_stl_result/result_meta_stl_cor.xlsx")



