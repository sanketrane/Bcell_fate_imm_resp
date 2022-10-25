## clearing the environment
rm(list = ls())  
gc()    

library(rstan)
library(loo)
library(tidyverse)
####################################################################################

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
modelName1 <- "Branched_timeinflux3"
modelName2 <- "Null_timeinflux1"

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

fit2 <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4, stanfit5))

## Time seq for predictions
ts_pred <- seq(4, 30, length.out = 500)

## combined plots
Y1pred1 <- as.data.frame(fit1, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)


Y2pred1 <- as.data.frame(fit1, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)


Y3pred1 <- as.data.frame(fit1, pars = "y3_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Y4pred1 <- as.data.frame(fit1, pars = "y4_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)




Y1pred2 <- as.data.frame(fit2, pars = "y1_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)


Y2pred2 <- as.data.frame(fit2, pars = "y2_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955))%>%
  bind_cols("timeseries" = ts_pred)


Y3pred2 <- as.data.frame(fit2, pars = "y3_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)

Y4pred2 <- as.data.frame(fit2, pars = "y4_mean_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.045),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.955)) %>%
  bind_cols("timeseries" = ts_pred)



## loading required datasets for plotting
imm_data <- read_csv(file.path(dataDir, "Bcell_imm_data.csv"))
imm_N2ko_data <- read_csv(file.path(dataDir, "N2KO_imm_data.csv"))


#### plots
p1 <- ggplot() +
  geom_line(data = Y2pred1, aes(x = timeseries, y = median), col =2) +
  #geom_line(data = Y2pred2, aes(x = timeseries, y = median), linetype=2, col ="#923347") +
  geom_ribbon(data = Y2pred1, aes(x = timeseries, ymin = lb, ymax = ub), fill=2, alpha = 0.25)+
  geom_line(data = Y4pred1, aes(x = timeseries, y = median), col =4) +
  geom_ribbon(data = Y4pred1, aes(x = timeseries, ymin = lb, ymax = ub), fill=4, alpha = 0.25)+
  #geom_line(data = Y4pred2, aes(x = timeseries, y = median), linetype=2, col ="darkblue") +
  geom_point(data = imm_data, aes(x = days_post_imm, y = CARpos_MZB), col=2) +
  geom_point(data = imm_N2ko_data, aes(x = days_post_imm, y = CARpos_MZB), col=4) +
  labs(title=paste("CAR positive MZ B cells"),  y=NULL, x="Days post immunization") + 
  xlim(3, 30) +
  scale_y_continuous(limits = c(2e3, 3e5), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e3, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")

p2 <- ggplot() +
  geom_line(data = Y1pred1, aes(x = timeseries, y = median), col =2) +
  #geom_line(data = Y1pred2, aes(x = timeseries, y = median), linetype=2, col ="#923347") +
  geom_ribbon(data = Y1pred1, aes(x = timeseries, ymin = lb, ymax = ub), fill=2, alpha = 0.15)+
  geom_point(data = imm_data, aes(x = days_post_imm, y = CARpos_GCB), col=2) +
  geom_point(data = imm_N2ko_data, aes(x = days_post_imm, y = CARpos_GCB), col=4) +
  labs(title=paste("CAR positive GC B cells"),  y=NULL, x="Days post immunization") + 
  xlim(3, 30) +
  scale_y_continuous(limits = c(5e3, 1e7), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal")


## saving  plots for quality control 
pdf(file = file.path(outputDir, paste("Combinedfit.pdf", sep = "")),
    width = 10, height = 4.5, onefile = FALSE, useDingbats = FALSE)
cowplot::plot_grid(p1, p2, ncol  = 2)
dev.off()



#### parameter plots
matrix_of_draws1 <- as.data.frame(fit1)   #matrix of parameter draws
alpha_pred1 <- quantile(0.5 * matrix_of_draws1$alpha, probs = c(0.5, 0.025, 0.975))
beta_pred1 <- quantile(0.5 * matrix_of_draws1$beta, probs = c(0.5, 0.025, 0.975))
lambda_WT_pred1 <- quantile(matrix_of_draws1$lambda_WT, probs = c(0.5, 0.025, 0.975))
lambda_N2KO_pred1 <- quantile(matrix_of_draws1$lambda_N2KO, probs = c(0.5, 0.025, 0.975))
delta_pred1 <- quantile(matrix_of_draws1$delta, probs = c(0.5, 0.025, 0.975))
mu_pred1 <- quantile(0.5 * matrix_of_draws1$mu, probs = c(0.5, 0.025, 0.975))
nu_pred1 <- quantile(sqrt(log(3)/matrix_of_draws1$nu), probs = c(0.5, 0.025, 0.975))

params_table <- t(round(data.frame(1/alpha_pred1,
                           1/beta_pred1,
                           lambda_WT_pred1,
                           lambda_N2KO_pred1,
                           delta_pred1,
                           1/mu_pred1,
                           nu_pred1), 4))

write.csv(params_table, 'params_table.csv', row.names = T)

pars_plot <- c("lambda_WT", "lambda_N2KO",  "delta", "beta")
parnames <- c('Clonal half-life of CAR+ MZ in WT mice', 'Clonal half-life of CAR+ MZ in N2KO mice', 'GC clonal half-life', "Propensity to gain CAR expression (%) for CAR- MZ B cells")
df_pars1 <- data.frame(t(data.frame(lambda_WT_pred1, lambda_N2KO_pred1, delta_pred1,  beta_pred1))) %>%
  mutate(parname = parnames, 
         pars_plot = pars_plot,
         Model = "Branched")

parnames_lambda <- c('Control/CAR', 'N2KO//CAR', 'Control/CAR', 'N2KO//CAR')
df_pars_lambda <- data.frame(t(data.frame(lambda_WT_pred1, lambda_N2KO_pred1, delta_pred1, delta_pred1))) %>%
  mutate(parname = parnames_lambda,
         Param = c("Clonal half-life of CAR positive MZB (days)", "Clonal half-life of CAR positive MZB (days)",
                   "Clonal half-life of GCB (days)", "Clonal half-life of GCB (days)"),
         cell_subset = c("MZ", "MZ", "GC", "GC"))
names(df_pars_lambda) <- c('Estimates', 'par_lb', 'par_ub', 'parname', 'Param', "Subset")
blank_data <- data.frame( Param = c("Clonal half-life of CAR positive MZB (days)", "Clonal half-life of CAR positive MZB (days)",
                                    "Clonal half-life of GCB (days)", "Clonal half-life of GCB (days)"),
                         Estimates = c(5,5,5,5),
                         parname = parnames_lambda,
                         Subset =  c("MZ", "MZ", "GC", "GC"))

ggplot(df_pars_lambda, aes(y=Estimates, x=factor(Subset), col=parname))+
  labs(y=NULL) +
  geom_errorbar(aes(y=Estimates, ymin=par_lb, ymax=par_ub, x=Subset),
                width=0.2, linetype=1,  position=position_dodge(0.4)) +
  geom_blank(data = blank_data)+
  geom_point(position=position_dodge(width=0.4), stat = "identity", size=4) + 
  facet_wrap(~ factor(Param), scales = "free") + 
  expand_limits(y = 0) + scale_y_continuous(expand = c(0.15, 0.1))+
  scale_color_manual(values=c(2, 4), name="Mouse strain")+
  myTheme + theme(axis.text.x=element_blank(),
                  axis.title.x=element_blank())+ theme(legend.background = element_blank(), legend.position = c(0.88, 0.85))


matrix_of_draws2 <- as.data.frame(fit2)   #matrix of parameter draws

alpha_pred2 <- quantile(0.5 * matrix_of_draws2$alpha, probs = c(0.5, 0.025, 0.975))
beta_pred2 <- quantile(0.5 * matrix_of_draws2$beta, probs = c(0.5, 0.025, 0.975))
lambda_WT_pred2 <- quantile(log(2)/matrix_of_draws2$lambda_WT, probs = c(0.5, 0.025, 0.975))
lambda_N2KO_pred2 <- quantile(log(2)/matrix_of_draws2$lambda_N2KO, probs = c(0.5, 0.025, 0.975))
delta_pred2 <- quantile(log(2)/matrix_of_draws2$delta, probs = c(0.5, 0.025, 0.975))
mu_pred2 <- quantile(matrix_of_draws2$mu, probs = c(0.5, 0.025, 0.975))
nu_pred2 <- quantile(sqrt(log(3)/matrix_of_draws2$nu), probs = c(0.5, 0.025, 0.975))


df_pars2 <- data.frame(t(data.frame(lambda_WT_pred2, lambda_N2KO_pred2, delta_pred2,  beta_pred2))) %>%
  mutate(parname = parnames, 
         pars_plot = pars_plot,
         Model = "Linear")

#formattable::formattable(all_pars_df)
all_pars_df <- rbind(df_pars1, df_pars2)
names(all_pars_df) <- c('Estimates', 'par_lb', 'par_ub', 'param', "pars_plot",  "Model")


blank_data <- data.frame(param = parnames,
                        Estimates = c(8, 8, 8, 8, 1.5, 1.5, 2, 2),
                         Model = rep(c('Branched', "Linear"), 4))


ggplot(all_pars_df, aes(y=Estimates, x=pars_plot, col=Model))+
  labs(y=NULL)+
  geom_point(position=position_dodge(width=0.4), size=3) +
  #geom_blank(data = blank_data)+
  geom_errorbar(aes(y=Estimates, ymin=par_lb, ymax=par_ub, x=pars_plot, col=Model),
                width=0.2, linetype=1,  position=position_dodge(0.4)) +
  facet_wrap(~ param, scales = "free", labeller = label_wrap_gen(width=42)) + 
  expand_limits(y = 0) + scale_y_continuous(expand = c(0, 0))+
  guides(col='none', fill='none')+
  myTheme + theme(axis.text.x=element_blank(),
                axis.title.x=element_blank())


FOtoCARMZ_pred1 <- as.data.frame(fit1, pars = "FOtoCARMZ_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.16),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.84))%>%
  bind_cols("timeseries" = ts_pred,
            "param" = "FOB to CAR+ MZB",
            "Model" = "CARMZ")

MZtoCARMZ_pred1 <- as.data.frame(fit1, pars = "MZtoCARMZ_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.16),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.84))%>%
  bind_cols("timeseries" = ts_pred,
            "param" = "CAR- MZB to CAR+ MZB",
            "Model" = "Branched")


## plots
FOtoBranch_pred <- rbind(FOtoCARMZ_pred1, MZtoCARMZ_pred1)

ggplot() +
  geom_line(data = FOtoBranch_pred, aes(x = timeseries, y = median*100, col = param)) +
  #geom_line(data = MZinflux_pred2, aes(x = timeseries, y = median*100), col =4) +
  geom_ribbon(data = FOtoBranch_pred, aes(x = timeseries, ymin = lb*100, ymax = ub*100, fill=param), alpha = 0.25) +
  #geom_ribbon(data = MZinflux_pred2, aes(x = timeseries, ymin = lb*100, ymax = ub*100), fill=4, alpha = 0.25)+
  labs(title=paste("Influx into CAR+ MZ (as % of CAR+ MZ)"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal") +
  #scale_x_log10(limits=c(3, 30)) +  scale_y_log10(limits=c(3, 125), breaks=c(3, 10, 30, 100)) +
  facet_wrap(~ param) + guides(col="none", fill="none")


FOtoCARGC_pred1 <- as.data.frame(fit1, pars = "FOtoCARGC_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.45),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.975))%>%
  bind_cols("timeseries" = ts_pred,
            "param" = "FO_to_GC",
            "Model" = "FOB to CAR+ GCB")


FOtoCARGC_pred2 <- as.data.frame(fit1, pars = "FOtoCARGC_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.45),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.975))%>%
  bind_cols("timeseries" = ts_pred,
            "param" = "FO_to_GC",
            "Model" = "Linear")


ggplot() +
  geom_line(data = FOtoCARGC_pred1, aes(x = timeseries, y = median, col = Model), size=1.2) +
  geom_ribbon(data = FOtoCARGC_pred1, aes(x = timeseries, ymin = lb, ymax = ub, fill=Model), alpha = 0.25) +
  #geom_line(data = FOtoCARGC_pred2, aes(x = timeseries, y = median, col = Model), size=1.2) +
  #geom_ribbon(data = FOtoCARGC_pred2, aes(x = timeseries, ymin = lb, ymax = ub, fill=Model), alpha = 0.25) +
  labs(title=paste("Influx into GC B cells (as % of GC)"),  y=NULL, x="Days post immunization") + 
  myTheme + theme(legend.position = c(0.5, 0.85), legend.direction = "horizontal") + 
  ylim(0, 8) + guides(col="none", fill="none") +
  scale_x_log10(limits=c(4, 30), breaks=c(5, 10, 20)) + #scale_y_log10(limits=c(1e-0, 10)) 
  facet_wrap(~ Model, scales = "free") 

cowplot::plot_grid(p1, p1.2, align = "v", nrow = 2)
cowplot::plot_grid(p2, p2.2, align = "v", nrow = 2)


muvec <- function(t, nu){
  2/(1 + exp(nu * t^2));
}
ts_pred <- seq(0, 30)
nuva = 0.0054
mus <- sapply(ts_pred, muvec, nu=nuva)
qplot(ts_pred, mus)
sqrt(log(3)/nuva)
muvec(14.26, nuva)


