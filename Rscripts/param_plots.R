## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(loo)
library(tidyverse)
####################################################################################


## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "Rscripts")
modelDir <- file.path(projectDir, "models")
dataDir <- file.path(projectDir, "datafiles")
toolsDir <- file.path(scriptDir, "tools")
outputDir <- file.path(projectDir, "output_fit")
saveDir <- file.path(projectDir, 'save_csv')
LooDir <- file.path('loo_fit') 

## model specific details that needs to be change for every run
modelName1 <- "totFOB_MZB_timeinflux"
modelName2 <- "carGC_MZB_timeinflux"

# compiling multiple stan objects together that ran on different nodes
stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName1, "_1", ".csv")))
stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName1, "_2",".csv")))
stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName1, "_3", ".csv")))
stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName1, "_4",".csv")))
stanfit5 <- read_stan_csv(file.path(saveDir, paste0(modelName1, "_5", ".csv")))
stanfit6 <- read_stan_csv(file.path(saveDir, paste0(modelName1, "_6",".csv")))

fit1 <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4, stanfit5, stanfit6))


#compiling multiple stan objects together that ran on different nodes
stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName2, "_1", ".csv")))
stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName2, "_2",".csv")))
stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName2, "_3", ".csv")))
stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName2, "_4",".csv")))
stanfit5 <- read_stan_csv(file.path(saveDir, paste0(modelName2, "_5", ".csv")))
stanfit6 <- read_stan_csv(file.path(saveDir, paste0(modelName2, "_6",".csv")))

fit2 <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4, stanfit5, stanfit6))

#### parameter plots
ts_pred <- seq(4, 30, length.out = 500)
car_fob_time <- function(t){
  F0 = exp(11.763019); B0 = exp(4.717021); n = 5.092933;
  X = 7.121932 ;  q = 6.475884;
  F0 + (B0 * t^n) * (1 - ((t^q)/((X^q) + (t^q))))
}
CAR_FOB <- sapply(ts_pred, car_fob_time)

matrix_of_draws1 <- as.data.frame(fit1)   #matrix of parameter draws

alpha1_pred1 <- quantile(matrix_of_draws1$alpha1* (CAR_FOB + exp(17.15)), probs = c(0.5, 0.025, 0.975))
lambda_WT_pred1 <- quantile(1/matrix_of_draws1$lambda_WT, probs = c(0.5, 0.025, 0.975))
lambda_N2KO_pred1 <- quantile(1/matrix_of_draws1$lambda_N2KO, probs = c(0.5, 0.025, 0.975))
delta_pred1 <- quantile(1/matrix_of_draws1$delta, probs = c(0.5, 0.025, 0.975))
mu_pred1 <- quantile(1/matrix_of_draws1$mu, probs = c(0.5, 0.025, 0.975))

parnames <- c('lambda_WT', 'lambda_N2KO', 'delta', "mu")
df_pars1 <- data.frame(t(data.frame(lambda_WT_pred1, lambda_N2KO_pred1, delta_pred1, mu_pred1))) %>%
  mutate(parname = parnames, 
         model = "M1")

matrix_of_draws2 <- as.data.frame(fit2)   #matrix of parameter draws

alpha1_pred2 <- quantile(matrix_of_draws2$alpha1* (CAR_FOB + exp(17.15)), probs = c(0.5, 0.025, 0.975))
lambda_WT_pred2 <- quantile(1/matrix_of_draws2$lambda_WT, probs = c(0.5, 0.025, 0.975))
lambda_N2KO_pred2 <- quantile(1/matrix_of_draws2$lambda_N2KO, probs = c(0.5, 0.025, 0.975))
delta_pred2 <- quantile(1/matrix_of_draws2$delta, probs = c(0.5, 0.025, 0.975))
mu_pred2 <- quantile(1/matrix_of_draws2$mu, probs = c(0.5, 0.025, 0.975))

parnames <- c('lambda_WT', 'lambda_N2KO', 'delta', "mu")
df_pars2 <- data.frame(t(data.frame(lambda_WT_pred2, lambda_N2KO_pred2, delta_pred2, mu_pred2))) %>%
  mutate(parname = parnames, 
         model = "M2")

all_pars_df <- rbind(df_pars1, df_pars2)
names(all_pars_df) <- c('Estimates', 'par_lb', 'par_ub', 'param', "model")

#### plotting style
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                 axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())
ggplot(all_pars_df, aes(y=Estimates, x=param, col= as.factor(model)))+
  geom_point(position=position_dodge(width=0.4), size=3) +
  geom_errorbar(aes(y=Estimates, ymin=par_lb, ymax=par_ub, x=param, col=model),
                 width=0.2, linetype=1,  position=position_dodge(0.4)) +
  scale_y_log10() +
  #scale_color_manual(values = c(4, 3, 2, 2, 2), label = real_modelname)+
  #facet_wrap(~ param) + scale_y_log10() + 
  myTheme 
