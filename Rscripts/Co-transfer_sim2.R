library(deSolve)
library(tidyverse)
library(readxl)
library(viridis)
library(hrbrthemes)
library(rstan)

## model specific details that needs to be change for every run
modelName <- "Linear_timeinflux1"
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

myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12), axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank(),
                 legend.text = element_text(size=11), legend.title = element_text(12))


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

# loadiong the scr# loadiong the script that contains functions for plotting stan parameters
source(file.path(toolsDir, "stanTools.R"))                # save results in new folder

# compiling multiple stan objects together that ran on different nodes
stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_1", ".csv")))
stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_2",".csv")))
stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_3", ".csv")))
stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_4",".csv")))
#stanfit5 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_5", ".csv")))
#stanfit6 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_6",".csv")))

fit <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4))
# finding the parameters used in the model 
# using the last parameter("sigma4") in the array to get the total number of parameters set in the model
num_pars <- which(fit@model_pars %in% "sigma1") -2      # the variable "sigma4" will change depdending on the data used
parametersToPlot <- fit@model_pars[1:num_pars] 

## extract posterior distribution  of parameters as a matrix
fit_ss <- as.matrix(fit, pars= parametersToPlot) # matrix of posterior samples
write.csv(fit_ss, paste0("PostDF_", modelName, ".csv"), row.names = FALSE)

#### donor CAR+ FO and GC B cells in recipients
## ode function for Branched model
ode_func <-  function(t, state, parms){
  with(as.list(c(state, parms)),{
    ## parameters estimated from modeling immunization data
    alpha = parms[1]
    beta = parms[2]
    mu = parms[3]
    delta = parms[4]
    lambda_WT = parms[5]
    lambda_N2KO = parms[6]
    nu = parms[7]
    
    ## new parameter
    epsilon = parms[8]
    
    t0 = 4.0;
    ## form for influx of activated FOB into GC varying with time
    alpha_tau = alpha/(1 + exp(nu * (t-t0)^2))
    
    #Fo b cells
    ## epsilon denotes the loss of Activated FO B cells in addition to their conversion to GC and MZ
    dY1 <- - (alpha_tau + mu + epsilon) * Y1
    
    # CAR positive GCB cells in WT
    dY2 = alpha_tau * Y1 - delta * Y2;
    
    # CAR positive MZB cells in WT
    ##influx into CAR+MZ is constant -- fitted better to immunization data!
    dY3 = mu * Y1 - lambda_WT * Y3;
    
    #return the rate of change
    list(c(dY1, dY2, dY3))
    
  })  # end with(as.list ...
}

#initial conditions
state <- c(Y1 = 1e6, Y2 = 0, Y3=0)

# time points for which conc is reported
# include the points where data is available
ts_pred <- seq(4, 90, length.out=100)


boot_func <- function(eps){
  ## Sampling random numbers to select the row from posterior dist matrix
  set.seed(1456)
  rowvec <- sample(1:nrow(fit_ss), 1000, replace = FALSE) %>% sort()
  
  ## initiating empty matrices to store the output
  Fo_df <- matrix(nrow = length(rowvec), ncol = length(ts_pred))
  GC_df <- matrix(nrow = length(rowvec), ncol = length(ts_pred))
  MZ_df <- matrix(nrow = length(rowvec), ncol = length(ts_pred))
  
  ## run the model with different sets of parameters sampled from the posterior dist matrix
  for (i in 1:length(rowvec)) {
    ## ith row from posterior dist matrix
    parms_vec <- c(fit_ss[rowvec[i], ], "epsilon" = eps)  ## epsilon denotes the loss of Activated FO B cells in addition to their conversion to GC and MZ
    ### simulation using the selected params
    sol_ode <- data.frame(ode(y= state, times=ts_pred, func = ode_func, parms = parms_vec))
    
    ### outpur
    Fo_df[i, ] <- sol_ode$Y1
    GC_df[i, ] <- sol_ode$Y2
    MZ_df[i, ] <- sol_ode$Y3
  }
  
  Fo_summary <- data.frame(Fo_df) %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.045),
              median = quantile(value, probs = 0.5),
              ub = quantile(value, probs = 0.955)) %>%
    bind_cols("timeseries" = ts_pred) %>%
    mutate(key = "FoB")
  
  GC_summary <- data.frame(GC_df) %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.045),
              median = quantile(value, probs = 0.5),
              ub = quantile(value, probs = 0.955)) %>%
    bind_cols("timeseries" = ts_pred) %>%
    mutate(key = "GCB")
  
  MZ_summary <- data.frame(MZ_df) %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.045),
              median = quantile(value, probs = 0.5),
              ub = quantile(value, probs = 0.955)) %>%
    bind_cols("timeseries" = ts_pred) %>%
    mutate(key = "MZB")
  
  
  rbind(Fo_summary, GC_summary, MZ_summary)%>%
    filter(timeseries >= 5)
}

pooled_df <- boot_func(1.35)
  
#ggplot(pooled_df) +
#  geom_line(aes(x=timeseries, y=median, col=key))+
#  geom_ribbon(aes(x = timeseries, ymin = lb, ymax = ub, fill=key), alpha = 0.25) +
#  labs(x="Days since immunization", y="Cell counts", title = NULL)+
#  scale_x_log10() + scale_y_log10(limits=c(1, 1e5))


sol_df <- pooled_df %>%
  select(-lb, -ub) %>%
  spread(key = key, value = median)
  

lab_vc <- c("MZB", "FoB", "GCB")

ggplot(data = sol_df, aes(x= timeseries)) +
  geom_area(aes(y= FoB, fill="FO_counts" ), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= GCB, fill="GC_counts"), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= MZB, fill="CARMZ_WT"), alpha=0.6 , linewidth=.5, colour="white") +
  #scale_fill_viridis(discrete = T) +
  labs(x="Days since immunization", title="Cell counts", y = NULL)+
  scale_fill_discrete(labels=lab_vc, name=NULL)+
  #theme_ipsum() + 
  scale_x_log10(limits = c(5, 100), breaks=c(7, 14, 35, 70)) + scale_y_log10(limits = c(1, 1e6)) + 
  labs(x='Days since immunization', y=NULL)  + 
  #guides(fill="none") + 
  myTheme

plot_df <- pooled_df %>%
  filter(key != "FoB") %>%
  filter(timeseries >= 5)

aresf <- plot_df %>%
  select(-lb, -ub) %>%
  group_by(timeseries) %>%
  mutate(n_tot = sum(median), 
         percentage = median/n_tot)

# Give a specific order:
aresf$popln <- factor(aresf$key , levels=c("CARMZ_WT", "GC_counts", "FO_counts") )

ggplot(data = aresf, aes(x= (timeseries), y=percentage*100, fill=key)) +
  geom_area(alpha=0.6 , linewidth=.5, colour="white") +
  scale_fill_viridis(discrete = T) +
  #theme_ipsum() + 
  xlim(5, 35) + #scale_y_log10() + 
  labs(x='Days since immunization', title ="% of differentiated cells in rep+ B cells", y=NULL)  + 
  #facet_wrap(.~popln) + 
  guides(col="none")
  

