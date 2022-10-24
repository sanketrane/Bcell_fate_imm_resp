## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(loo)
library(tidyverse)
####################################################################################

## model specific details that needs to be change for every run
modelName <- "Null_timeinflux1"
data_der <- "Bcell_imm_data.csv"    
data_der2 <- "N2KO_imm_data.csv"    

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "Rscripts")
modelDir <- file.path(projectDir, "models")
dataDir <- file.path(projectDir, "datafiles")
toolsDir <- file.path(scriptDir, "tools")
outputDir <- file.path(projectDir, "output_fit")
saveDir <- file.path(projectDir, 'save_csv')
LooDir <- file.path('loo_fit') 

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))                # save results in new folder

# compiling multiple stan objects together that ran on different nodes
stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_1", ".csv")))
stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_2",".csv")))
stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_3", ".csv")))
stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_4",".csv")))
stanfit5 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_5", ".csv")))
stanfit6 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_6",".csv")))

fit <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4, stanfit5, stanfit6))

# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars <- which(fit@model_pars %in% "sigma3")      # the variable "sigma4" will change depdending on the data used
parametersToPlot <- fit@model_pars[1:num_pars]

# number of post-burnin samples that are used for plotting 
nPost <- nrow(fit)

################################################################################################
################################################################################################

## loading required datasets for plotting
imm_data <- read_csv(file.path(dataDir, data_der))
imm_N2ko_data <- read_csv(file.path(dataDir, data_der2))

# ################################################################################################
# calculating PSIS-L00-CV for the fit
MZ_fractions_loglik <- extract_log_lik(fit, parameter_name = "log_lik2", merge_chains = TRUE)
GC_fractions_loglik <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = TRUE)
MZN2_fractions_loglik <- extract_log_lik(fit, parameter_name = "log_lik4", merge_chains = TRUE)
GCN2_fractions_loglik <- extract_log_lik(fit, parameter_name = "log_lik3", merge_chains = TRUE)

log_lik_comb <- cbind(MZ_fractions_loglik, GC_fractions_loglik,
                      MZN2_fractions_loglik, GCN2_fractions_loglik)
# optional but recommended
ll_array <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = FALSE)
r_eff <- relative_eff(exp(ll_array))

# loo-ic values
loo_loglik <- loo(log_lik_comb, save_psis = FALSE, cores = 8)
loofilename <- paste0("loosave_", modelName, "_", data_der, ".rds")
write_rds(loo_loglik, file  = file.path(LooDir, loofilename))


# Widely applicable AIC
AICw_lok <- waic(MZ_fractions_loglik, GC_counts_loglik, GC_fractions_loglik)
loo_loglik

ploocv <- data.frame("Model" = modelName,
                     "LooIC" = loo_loglik$estimates[3],
                     "SE" = loo_loglik$estimates[6], 
                     "PLoo" = loo_loglik$estimates[2])

write.table(ploocv, file = file.path(outputDir, 'timeinfluxfit', "stat_table_MZB1.csv"),
            sep = ",", append = T, quote = FALSE,
            col.names = F, row.names = FALSE)


################################################################################################
################################################################################################
## posterior predictive distributions
# time sequence for predictions 
ts_pred <- seq(4, 30, length.out = 500)
numPred <- length(ts_pred)


#### plotting style
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                 axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank())

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())


####### plotting
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(3:8))


Y1pred <- as.data.frame(fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)


Y2pred <- as.data.frame(fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)

Y3pred <- as.data.frame(fit, pars = "y3_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)


Y4pred <- as.data.frame(fit, pars = "y4_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

GCcounts_pred <- as.data.frame(fit, pars = "CAR_GCcounts_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)

MZcounts_pred <- as.data.frame(fit, pars = "CAR_MZcounts_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)


GCN2counts_pred <- as.data.frame(fit, pars = "CAR_GCN2counts_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)


MZN2counts_pred <- as.data.frame(fit, pars = "CAR_MZN2counts_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)



#### plots
p1 <- ggplot() +
  #geom_hline(yintercept = exp(10.8))+
  geom_line(data = Y2pred, aes(x = timeseries, y = median), col =2) +
  geom_ribbon(data = Y2pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=2, alpha = 0.25)+
  geom_line(data = Y4pred, aes(x = timeseries, y = median), col =4) +
  geom_ribbon(data = Y4pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=4, alpha = 0.25)+
  #geom_ribbon(data = MZfractions_pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=2, alpha = 0.25)+
  geom_point(data = imm_data, aes(x = days_post_imm, y = CARpos_MZB), col=2) +
  geom_point(data = imm_N2ko_data, aes(x = days_post_imm, y = CARpos_MZB), col=4) +
  labs(title=paste("CAR positive MZ B cells"),  y=NULL, x="Days post immunization") + 
  xlim(0, 30) +
  scale_y_continuous(limits = c(2e3, 3e5), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e3, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p2 <- ggplot() +
  geom_line(data = Y1pred, aes(x = timeseries, y = median), col =2) +
  geom_ribbon(data = Y1pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=2, alpha = 0.15)+
  #geom_ribbon(data = MZfractions_pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=2, alpha = 0.25)+
  geom_point(data = imm_data, aes(x = days_post_imm, y = CARpos_GCB), col=2) +
  geom_line(data = Y3pred, aes(x = timeseries, y = median), col =2) +
  geom_ribbon(data = Y3pred, aes(x = timeseries, ymin = lb, ymax = ub), fill="#ba6dd1", alpha = 0.15)+
  #geom_ribbon(data = MZfractions_pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=2, alpha = 0.25)+
  geom_point(data = imm_N2ko_data, aes(x = days_post_imm, y = CARpos_GCB), col=4) +
  labs(title=paste("CAR positive GC B cells"),  y=NULL, x="Days post immunization") + 
  xlim(0, 30) +
  scale_y_continuous(limits = c(5e3, 1e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")


### Residual plots
resid_d1  <- t(as.data.frame(fit, pars = "resid_d1"))[,1]
resid_d2  <- t(as.data.frame(fit, pars = "resid_d2"))[,1]
resid_d3  <- t(as.data.frame(fit, pars = "resid_d3"))[,1]
resid_d4  <- t(as.data.frame(fit, pars = "resid_d4"))[,1]


p1.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = imm_data, aes(x = days_post_imm, y = resid_d1), col=6) +
  labs(title=paste("Residuals CAR GC counts WT"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p2.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = (imm_data), aes(x = days_post_imm, y = resid_d2), col=2) +
  labs(title=paste("Residuals CAR MZ counts WT"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p3.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = imm_N2ko_data, aes(x = days_post_imm, y = resid_d3), col=4) +
  labs(title=paste("Residuals GC counts N2KO"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p4.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = imm_N2ko_data, aes(x = days_post_imm, y = resid_d4), col=4) +
  labs(title=paste("Residuals MZ counts N2KO"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

cowplot::plot_grid(p1, p2, p1.1, p2.1, p3.1, p4.1, nrow  = 3)


## saving  plots for quality control 
pdf(file = file.path(outputDir, 'timeinfluxfit', paste(modelName,"StanPlots%03d.pdf", sep = "")),
    width = 9, height = 4, onefile = FALSE, useDingbats = FALSE)
cowplot::plot_grid(p1, p2, ncol  = 2)
dev.off()

################################################################################################
### parameters table
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
out_table

write.csv(out_table, file = file.path(outputDir, 'timeinfluxfit', paste0('params_', modelName, ".csv")))

#time_shape <- function(Time, delta, nu){
#  delta/(1 + exp(-nu *(Time-4)^2))
#}
#
#Time_pred <- time_shape(seq(4, 30, length.out=100), out_table$mean[4], out_table$mean[6])
delta_pred <- as.data.frame(fit, pars = "delta_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.055),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.945)) %>%
  bind_cols("timeseries" = ts_pred)

p4 <- ggplot() +
  #\geom_line(aes(x = seq(4, 30, length.out=100), y = Time_pred), col =4, size=1.52) +
  geom_line(data = delta_pred, aes(x = timeseries, y = median), col ="#e63590") +
  geom_ribbon(data = delta_pred, aes(x = timeseries, ymin = lb, ymax = ub), fill="#ba6dd1", alpha = 0.25)+
  labs(title=paste("Time dependent Loss rate of GC B cells"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal") +
  scale_y_log10(limits=c(0.2, 1.5), breaks=c(0.25, 0.5,1))

p4
## saving  plots for quality control 
pdf(file = file.path(outputDir, paste(modelName,"ExtraPlots%03d.pdf", sep = "")),
    width = 5, height = 4, onefile = FALSE, useDingbats = FALSE)
p4
dev.off()


#alpha_shape <-  function(Time, alpha, nu){
#  alpha/(1 + exp(nu * (Time-4)^2))
#}

#alpha_pred <- alpha_shape(seq(4, 30, length.out=500), out_table$mean[2], out_table$mean[7])
alpha_pred <- as.data.frame(fit, pars = "alpha_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.16),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.84))%>%
  bind_cols("timeseries" = ts_pred)

p4 <- ggplot() +
  #geom_line(aes(x = seq(4, 30, length.out=500), y = alpha_pred), col =4, size=1.52) +
  geom_line(data = alpha_pred, aes(x = timeseries, y = median), col ="#e63590") +
  geom_ribbon(data = alpha_pred, aes(x = timeseries, ymin = lb, ymax = ub), fill="#ba6dd1", alpha = 0.25)+
  labs(title=paste("Time dependent influx of FO B cells into GCB"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal") +
  scale_x_log10() + scale_y_log10(limits=c(1e-3, 0.1))

p4


mu_pred <- as.data.frame(fit, pars = "mu_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.16),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.84))%>%
  bind_cols("timeseries" = ts_pred)

p5 <- ggplot() +
  #geom_line(aes(x = seq(4, 30, length.out=500), y = alpha_pred), col =4, size=1.52) +
  geom_line(data = mu_pred, aes(x = timeseries, y = median), col ="#e63590") +
  geom_ribbon(data = mu_pred, aes(x = timeseries, ymin = lb, ymax = ub), fill="#ba6dd1", alpha = 0.25)+
  geom_line(data = alpha_pred, aes(x = timeseries, y = median), col ="#e63590") +
  geom_ribbon(data = alpha_pred, aes(x = timeseries, ymin = lb, ymax = ub), fill="#ba6dd1", alpha = 0.25)+
  labs(title=paste("Time dependent influx of FO B cells into MZB"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal") +
  scale_x_log10() + scale_y_log10(limits=c(1e-5, 1e-1))

p5
## saving  plots for parameter estimates  
pdf(file = file.path(outputDir, paste(modelName,"ExtraPlots%03d.pdf", sep = "")),
    width = 5, height = 4, onefile = FALSE, useDingbats = FALSE)
p4;p5
dev.off()




## open graphics device 
## saving  plots for quality control 
pdf(file = file.path(outputDir, paste(modelName,"StanPlots%03d.pdf", sep = "")),
    width = 8, height = 5, onefile = FALSE, useDingbats = FALSE)
pairs(fit, pars = parametersToPlot)
options(bayesplot.base_size = 15,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "viridis")

myTheme <- theme(text = element_text(size = 10), axis.text = element_text(size = 10), axis.title =  element_text(size = 10, face = "bold"),
               plot.title = element_text(size=10,  hjust = 0.5, face = "bold"),
               legend.background = element_blank(), legend.key = element_blank(),
               legend.text = element_text(size=9), legend.title = element_text(9))

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

posterior <- as.array(fit)
mcmc_acf(posterior, pars = parametersToPlot) + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

mcmc_dens_overlay(posterior, parametersToPlot)
mcmc_dens(posterior, parametersToPlot) + myTheme

dev.off()


