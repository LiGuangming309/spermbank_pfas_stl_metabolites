################### 2.计算MWAS-BHRMA-G结果 #####################################  

suppressMessages(library(readxl, quietly = TRUE))
suppressMessages(library(R2jags, quietly = TRUE))
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(rjags, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(pbdMPI, quietly = TRUE))
suppressMessages(library(lme4, quietly = TRUE))
.comm.size <- comm.size()
.comm.rank <- comm.rank()

# Create Function ---------------------------------------------
# library(lme4)

BHRMA.g <- function(X=NULL, Y=NULL, U=NULL, LOD=NULL, profiles=NULL) {
  # JAGS model
  ridge.BDL.model <- 
    "model {
  for(i in 1:N) {
    Y[i] ~ dnorm(mu[i], prec.sigma.Y)
    mu[i] <- alpha + inprod(beta[1:P], X.s[i,1:P]) + inprod(delta[1:Q], U[i,1:Q])
    
    # imputation BDL
    for(p in 1:P) {
      X[i,p] ~ dnorm(X.true[i,p],prec.X[p]) 
      X.true[i,p] <- X.notmiss[i,p]*(1-R[i,p]) + X.miss[i,p]*R[i,p]
      X.notmiss[i,p] ~ dnorm(mu.X[p], tau.X[p])T(LOD[p], )
      X.miss[i,p] ~ dnorm(mu.X[p], tau.X[p])T( , LOD[p])
      X.s[i,p] <- (X.true[i,p] - mu.X[p])/sigma.X[p]
    }
  }
  # prior on outcome variance
  prec.sigma.Y <- 1/(sigma.Y*sigma.Y)
  sigma.Y ~ dunif(0,3)
  
  # prior on covariate effects
  for(q in 1:Q) { delta[q] ~ dnorm(0, 1.0E-06) }
  
  # prior on intercept
  alpha ~ dnorm(0, 1.0E-06)
  
  # prior on exposure effects
  beta[1:P] ~ dmnorm(mu.beta[1:P], T[1:P, 1:P])
  for(j in 1:P) {
    mu.beta[j] <- (1-gamma[j])*prop.mu.beta[j]
    b[j] <- beta[j]*gamma[j]
    gamma[j] ~ dbern(pi)
    for(k in 1:P) {
      T[j,k] <- gamma[j]*gamma[k]*XtX[j,k]/(G) + (1-gamma[j]*gamma[k])*equals(j,k)*pow(prop.sd.beta[j],-2)
    }
    tau.X[j] <- 1/(sigma.X[j]*sigma.X[j])
    sigma.X[j] ~ dunif(0,5)
    mu.X[j] ~ dnorm(0, 1.0E-06)
    prec.X[j] <- 10000
  }
  pi ~ dbeta(1,P) 

  # Hyper-g prior (following Perrakis 2018, note that this is on the G^-1 so the Beta distribution is switchd in terms of a and b from Li and Clyde 2019 equation 34)
  a <- 3
  bw <- a/2 - 1
  w~dbeta(1,bw) T(0.00001,0.99999) # JG: added T(0.0001,0.9999) 
  G <- w/(1-w)

  # g-estimation
  eta.low <- inprod(b[1:P], profiles[1,1:P])
  eta.high <- inprod(b[1:P], profiles[2,1:P])
  psi <- eta.high-eta.low
  
}"
  
  # Other Stuff
  N <- length(Y)
  P <- ncol(X)
  Q <- ncol(U)
  R <- ifelse(is.na(X), 1,0)
  exposure.Names <- colnames(X)
  
  ### get the univariate result
  univariate.results <- t(sapply(1:P, FUN=function(p) {  # using index p facilitate write
    x <- as.matrix(X[,p])
    reg <- glm(Y~x, family=gaussian)    # perform logistic regression 
    s.reg <- summary(reg)     # get the summary for the regression
    c.reg <- s.reg$coef["x",] # select the coefficients for the exposure
    return(c.reg)             # to avoid potential memory issues only return coefficients if small number of exposures
  }, simplify=T))
  univariate.results <- data.frame(exposure.Names,univariate.results)
  
  ### g prior model result
  prop.mu.beta <- rep(0, P)
  prop.sd.beta <- univariate.results$Std..Error
  XtX <- t(as.matrix(X))%*%as.matrix(X) 
  
  # run jags
  jdata <- list(N=N, Y=Y, X=X, R=R, U=U, P=P, Q=Q, profiles=profiles, LOD=LOD,XtX=XtX, prop.mu.beta=prop.mu.beta, prop.sd.beta=prop.sd.beta)
  var.s <- c("beta", "gamma", "eta.low", "eta.high",  "psi")
  model.fit <- jags.model(file=textConnection(ridge.BDL.model), data=jdata, n.chains=1, n.adapt=10000, quiet=T)
  update(model.fit, n.iter=10000, progress.bar="none")
  model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=10000, thin=1, progress.bar="none")
  
  # summarize results
  r <- summary(model.fit)
  var.names <- c(paste(exposure.Names, "beta", sep="."),
                 "eta.high",
                 "eta.low",
                 paste(exposure.Names, "gamma", sep="."),
                 "psi")
  ridge.BDL.results <- data.frame(var.names, 
                                  r$statistics[,1:2], 
                                  r$quantiles[,c(1,5)])
  wald = abs(ridge.BDL.results[,"Mean"]/ridge.BDL.results[,"SD"])
  ridge.BDL.results$p.val = (2*(1-pnorm(wald,0,1)))
  return(ridge.BDL.results)
}
# example function call:
# BHRMA(X=X.obs, Y=Y, U=U, LOD=LOD, profiles=profiles)

############################读入数据#################################
df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA <- read_excel("/dssg/home/acct-yxwang0203/yxwang0203/work_data/pfas_STL_meta/pfas_20250324/result_pfas_rawdata_623.xlsx")


# df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA <- read_excel("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/metabolomics/pfas_meta_stl_20241216/data_preparation/result_pfas_rawdata_623.xlsx")
# X: A NxP matrix of exposures for mixture analysis (on the original scale with NA's for individuals with BLD)
X.obs = df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA%>% 
  select(PFOS,PFOA,PFDA,PFUdA,L_PFHxS, PF3ONS_9CL) %>% 
  mutate(across(everything(), ~scale(.)))


# X.obs

# P: the number of exposures in the mixture
P = ncol(X.obs)
# P

# Y: A N-length vector for a continuous outcome
Y = df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA %>% 
  select(M661:M625) %>% 
  scale(center = F,scale = F) ####代谢化合物已经进行了对数转换和归一化,这里不再进行对数，只是将这里的数据格式进行转化
# Y
# U: A NxQ matrix of covariates (variables included in the regression model but not included in the g-estimation)
# Covariates

# str(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA)

U = cbind.data.frame(age_class_2.num = as.numeric(as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$age_class_2)),
                     BMI_class_2.num = as.numeric(as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$BMI_class_2)),
                     education.num = as.numeric(as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$education)),
                     marriage_class_2.num = as.numeric(as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$marriage_class_2)),
                     income.num = as.numeric(as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$income)),  ###基线变量与代谢组的变量重复了
                     smk_class_2.num = as.numeric(as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$smk_class_2)),
                     drk_class_2.num = as.numeric(as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$drk_class_2)),
                     season.num = as.numeric(as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$season)),  ###基线变量与代谢组的变量重复了
                     abstinence_class_2.num = as.numeric(as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$abstinence_class_2)))

# U = cbind.data.frame(age_class_2.num = as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$age_class_2),
#                      BMI_class_2.num = as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$BMI_class_2),
#                      education.num = as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$education),
#                      marriage_class_2.num = as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$marriage_class_2),
#                      income.num = as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$income),  ###基线变量与代谢组的变量重复了
#                      smk_class_2.num = as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$smk_class_2),
#                      drk_class_2.num = as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$drk_class_2),
#                      season.num = as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$season),  ###基线变量与代谢组的变量重复了
#                      abstinence_class_2.num = as.factor(df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$abstinence_class_2))





# str(U)
id <- df_background_DNA_LT_sperm_PFOA_meta_tidy_noNA$number
# LOD: A P-length vector of LODs for each exposure. Individuals with missing data will have data imputed below this level of detection  
LOD = c(0.0516440006885867,0.00675113081441142,0.00397719740156436, 0.00351658656663932, 0.00338104361546264, 0.000503330369276714)
# LOD
# profiles: A 2xP matrix of two counterfactual profiles of exposures for which a potential outcomes risk difference is calculated (as the exposures are standardized within the function, these profiles should be on the standard normal scale)
profiles = c(-1,1)*matrix(.5, nrow=2, ncol=P)
# profiles
# Names of exposures
exposure.Names = colnames(X.obs)
# exposure.Names

# Get number of metabolites
n_met <- ncol(Y)

# Model

# Model
model <- function(i){
  output <- BHRMA.g(X=X.obs,
                    Y=Y[,i],
                    U=U,
                    LOD=LOD,
                    profiles=profiles)
  # Sys.sleep(1) # 模拟耗时操作
  # Name Metabolite
  output$metabolite = colnames(Y)[i]
  return(output)
}
pbdLapply(1:n_met, model, pbd.mode = "spmd") -> coefs



g_coefs_PFAS_meta_mean <- allgather(coefs)

if(comm.rank() == 0){
  save(g_coefs_PFAS_meta_mean, file = 'coefs_PFAS_meta_mean.Rdata')
}

# 结束 MPI 环境
finalize()
