## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(loo)
library(tidyverse)
####################################################################################

## model specific details that needs to be change for every run
modelName <- "GC_desc_dens_timeLoss"
data_der <- "Bcell_imm_data.csv"    

## Setting all the directories for opeartions
projectDir <- getwd()
scriptDir <- file.path(projectDir, "Rscripts")
modelDir <- file.path(projectDir, "models")
dataDir <- file.path(projectDir, "datafiles")
toolsDir <- file.path(scriptDir, "tools")
outputDir <- file.path(projectDir, "output_fit")
saveDir <- file.path(projectDir, 'save_csv')

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))                # save results in new folder

# compiling multiple stan objects together that ran on different nodes
stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_1", ".csv")))
stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_2",".csv")))
stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_3", ".csv")))
stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_4",".csv")))

fit <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4))

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

# ################################################################################################
# calculating PSIS-L00-CV for the fit
MZ_fractions_loglik <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = TRUE)
GC_counts_loglik <- extract_log_lik(fit, parameter_name = "log_lik3", merge_chains = TRUE)
GC_fractions_loglik <- extract_log_lik(fit, parameter_name = "log_lik2", merge_chains = TRUE)

log_lik_comb <- cbind(MZ_fractions_loglik, GC_counts_loglik, GC_fractions_loglik)
# optional but recommended
ll_array <- extract_log_lik(fit, parameter_name = "log_lik1", merge_chains = FALSE)
r_eff <- relative_eff(exp(ll_array))

# loo-ic values
loo_loglik <- loo(log_lik_comb, save_psis = FALSE, cores = 8)

# Widely applicable AIC
AICw_lok <- waic(MZ_fractions_loglik, GC_counts_loglik, GC_fractions_loglik)
loo_loglik

ploocv <- data.frame("Model" = modelName,
                     "LooIC" = loo_loglik$estimates[3],
                     "SE" = loo_loglik$estimates[6], 
                     "PLoo" = loo_loglik$estimates[2])

write.table(ploocv, file = file.path(outputDir, "stat_table.csv"),
            sep = ",", append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE)


################################################################################################
################################################################################################
## posterior predictive distributions
# time sequence for predictions 
ts_pred <- seq(0, 30, length.out = 300)
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

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))


Y1pred <- as.data.frame(fit, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)


MZfractions_pred <- as.data.frame(fit, pars = "MZfractions_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)


Y2pred <- as.data.frame(fit, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)

GCfractions_pred <- as.data.frame(fit, pars = "GCfractions_pred") %>%
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

GCcounts_pred <- as.data.frame(fit, pars = "GCcounts_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)


#### plots
p1 <- ggplot() +
  geom_line(data = Y1pred, aes(x = timeseries, y = median), col =2) +
  geom_ribbon(data = Y1pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=2, alpha = 0.25)+
  #geom_ribbon(data = MZfractions_pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=2, alpha = 0.25)+
  geom_point(data = imm_data, aes(x = days.post.imm, y = fraction_CAR_MZ), col=2) +
  labs(title=paste("Fraction CAR in MZ B cells"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p2 <- ggplot() +
  geom_line(data = Y2pred, aes(x = timeseries, y = median), col =6) +
  geom_ribbon(data = Y2pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=6, alpha = 0.25)+
  geom_point(data = imm_data, aes(x = days.post.imm, y = fraction_CAR_GC), col=6) +
  labs(title=paste("Fraction CAR in GC B cells"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p3 <- ggplot() +
  geom_line(data = Y3pred, aes(x = timeseries, y = median), col =4) +
  geom_ribbon(data = Y3pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=4, alpha = 0.25)+
  #geom_ribbon(data = GCcounts_pred, aes(x = timeseries, ymin = lb, ymax = ub), fill=4, alpha = 0.25)+
  geom_point(data = imm_data, aes(x = days.post.imm, y = GCB_cell_numbers), col=4) +
  labs(title=paste("Total numbers of GC B cells"),  y=NULL, x="Days post immunization") + 
  scale_y_continuous(limits = c(1e4, 2e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")


### Residual plots
resid_d1  <- t(as.data.frame(fit, pars = "resid_d1"))[,1]
resid_d2  <- t(as.data.frame(fit, pars = "resid_d2"))[,1]
resid_d3  <- t(as.data.frame(fit, pars = "resid_d3"))[,1]

p1.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = (imm_data), aes(x = days.post.imm, y = resid_d1), col=2) +
  labs(title=paste("Residuals MZ Fractions fit"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p2.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = imm_data, aes(x = days.post.imm, y = resid_d2), col=6) +
  labs(title=paste("Residuals GC Fractions fit"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p3.1 <- ggplot() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(data = imm_data, aes(x = days.post.imm, y = resid_d3), col=4) +
  labs(title=paste("Residuals GC counts fit"),  y=NULL, x="Time") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

cowplot::plot_grid(p1, p2,  p3, p1.1, p2.1, p3.1, nrow  = 2)


## saving  plots for quality control 
pdf(file = file.path(outputDir, paste(modelName,"StanPlots%03d.pdf", sep = "")),
    width = 11, height = 6.5, onefile = FALSE, useDingbats = FALSE)
cowplot::plot_grid(p1, p2,  p3, p1.1, p2.1, p3.1, nrow  = 2)
dev.off()

################################################################################################
### parameters table
ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
out_table <- ptable[1:num_pars, c(1, 3, 4, 8)]
out_table

write.csv(out_table, file = file.path(outputDir, paste0('params_', modelName, ".csv")))

time_shape <- function(Time, delta, nu, b0){
  tau=14
  #delta * exp(nu1 * Time)/(1 + exp(nu2 * (Time - tau)))
  delta*(1 + exp(-nu * (Time - b0)^2))#/(1 + exp(nu2 * (Time - tau)))
  #(delta*(1 - exp(-nu * Time)^1)) + b0
  #(delta*(1 + (Time/nu)^2)) + b0
}

Time_pred <- time_shape(seq(0, 30, length.out=100), out_table$mean[4], out_table$mean[3], out_table$mean[5])

p4 <- ggplot() +
  geom_line(aes(x = seq(0, 30, length.out=100), y = Time_pred), col =4, size=1.52) +
  labs(title=paste("Time dependent Loss rate"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p4
## saving  plots for quality control 
pdf(file = file.path(outputDir, paste(modelName,"ExtraPlots%03d.pdf", sep = "")),
    width = 6, height = 4.5, onefile = FALSE, useDingbats = FALSE)
p4
dev.off()


alpha_shape <- function(Time, alpha, nu, b0){
  tau=14
  #delta * exp(nu1 * Time)/(1 + exp(nu2 * (Time - tau)))
  alpha*(1 + exp(-nu * (Time - b0)^2))#/(1 + exp(nu2 * (Time - tau)))
  #alpha/(1 + exp(-nu * Time)^1) + b0
}

alpha_pred <- alpha_shape(seq(0, 50), out_table$mean[2], out_table$mean[3], out_table$mean[5])

p4 <- ggplot() +
  geom_line(aes(x = seq(0, 50), y = alpha_pred), col =4, size=1.52) +
  labs(title=paste("Time dependent activation rate"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p4

## saving  plots for quality control 
pdf(file = file.path(outputDir, paste(modelName,"ExtraPlots%03d.pdf", sep = "")),
    width = 6, height = 4.5, onefile = FALSE, useDingbats = FALSE)
p4
dev.off()


dens_shape <- function(counts, delta, size_dens_Log){
  size_dens = exp(size_dens_Log)
  delta/(1 + (counts/size_dens)^2);
}

dens_pred <- dens_shape(Y3pred$median, out_table$mean[7], out_table$mean[2])

ggplot() +
  geom_line(aes(x = ts_pred, y = dens_pred), col =4) +
  labs(title=paste("Total numbers of GC B cells"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")


write.csv(ploocv, file = file.path(outputDir, paste0('stats_', modelName, ".csv")))

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


pdf(file = file.path(outputDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 7, height = 5, onefile = FALSE, useDingbats = FALSE )

lay1 <- rbind(c(1,1,2,2),
              c(NA,3,3,NA))


gridExtra::grid.arrange(counts_facet, nfd_comb, ki67_comb, layout_matrix = lay1)

toprow <- cowplot::plot_grid(counts_facet, nfd_comb,  labels = c("A", "B"), ncol = 2)
cowplot::plot_grid(toprow, ki67_plot,  labels = c("", "C"), nrow = 2)


dev.off()

