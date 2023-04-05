library(deSolve)
library(tidyverse)
library(hrbrthemes)
library(rstan)

## model specific details that needs to be change for every run
modelName <- "Branched_timeinflux1"
modelName <- "Linear_timeinflux1"

### extract posterior distribution  of parameters as a df
fit_ss <- read.csv(file = paste0("PostDF_", modelName, ".csv"))
mu_b <- median(fit_ss$mu)
lambda_b <- median(fit_ss$lambda_WT)
alpha_b <- median(fit_ss$alpha)
nu_b <- median(fit_ss$nu)


fit_ss2 <- read.csv(file = paste0("PostDF_", modelName, ".csv"))
delta_l <- median(fit_ss2$delta)
lambda_l <- median(fit_ss2$lambda_WT)
mu_l <- median(fit_ss2$mu)

#### donor CAR+ FO and GC B cells in recipients
## ode function for Branched model
ode_func <-  function(t, state, parms){
  with(as.list(c(state, parms)),{
    ## parameters estimated from modeling immunization data
    alpha = parms[1]
    mu_branched = parms[2]
    mu_linear = parms[3]
    delta = parms[4]
    lambda_Branched = parms[5]
    lambda_linear = parms[6]
    nu = parms[7]
    
    ## new parameter
    epsilon = parms[8]
    
    t0 = 4.0;
    ## form for influx of activated FOB into GC varying with time
    alpha_tau = alpha/(1 + exp(nu * (t-t0)^2))
    
    #Fo b cells cd45.1
    ## epsilon denotes the loss of Activated FO B cells in addition to their conversion to GC and MZ
    dY1 <- - (alpha_tau + mu_branched + epsilon) * Y1
    
    # CAR positive MZB CD45.1
    ##influx into CAR+MZ is constant -- fitted better to immunization data!
    dY2 = mu_branched * Y1 - lambda_Branched * Y2;
    
    # CAR positive GCB CD45.2
    dY3 = - (delta + mu_linear) * Y3;
    
    # CAR positive MZB CD45.2
    dY4 = mu_linear * Y3 - lambda_linear * Y4;
    
    #return the rate of change
    list(c(dY1, dY2, dY3, dY4))
    
  })  # end with(as.list ...
}

#initial conditions
state <- c(Y1 = 1e5, Y2 = 0, Y3=1e5, Y4=0)

# time points for which conc is reported
# include the points where data is available
ts_pred <- seq(7, 60, length.out=100)


out_df <- data.frame()
eps_vec <- c(1/10, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5)

for(i in 1:length(eps_vec)) {
  parms_vec <- c("alpha" = alpha_b, "mu_branched" = mu_b , "mu_linear" = mu_l, 
                 "delta" = delta_l, "lambda_branched" = lambda_b, "lambda_linear" = lambda_l,
                 "nu" = nu_b, "epsilon" = eps_vec[i] * alpha_b)  ## epsilon denotes the loss of Activated FO B cells in addition to their conversion to GC and MZ
  ### simulation using the selected params
  sol_ode <- data.frame(ode(y= state, times=ts_pred, func = ode_func, parms = parms_vec))
  
  sol_df <- sol_ode %>%
    rename(timeseries = time,
           FO.1 = Y1,
           MZ.1 = Y2,
           GC.2 = Y3,
           MZ.2 = Y4) %>%
    mutate(ratio_MZ = replace_na(MZ.2/MZ.1, 0),
           EPS = round(eps_vec[i], 1))
  out_df <- rbind(out_df, sol_df)
  out_df
}



#ggplot(data = sol_df, aes(x= timeseries)) +
#  geom_area(aes(y= FO.1 + MZ.1 + GC.2 + MZ.2, fill="FO CD45.1" ), alpha=0.6, colour="white") +
#  geom_area(aes(y= MZ.1 + GC.2 + MZ.2, fill="MZ CD45.1" ), alpha=0.6, colour="white") +
#  geom_area(aes(y= GC.2 + MZ.2, fill="GC CD45.2" ), alpha=0.5, colour="white") +
#  geom_area(aes(y= MZ.2, fill="MZ CD45.2"), alpha=0.5, colour="white") +
#  #scale_fill_viridis(discrete = T) +
#  #theme_ipsum() + 
#  scale_x_log10(limits = c(7, 60), breaks=c(7, 14, 28, 56)) + scale_y_log10() + 
#  labs(x='Days since immunization', y=NULL)  +  theme_ipsum() +
#  theme(legend.title = element_blank())


plot_labl <- c(`0.1` = "eps=0.1", `0.2` = "eps=0.25", `0.3` = "eps=0.33", `0.5` = "eps=0.5", 
               `1` = "eps=1", `2` = "eps=2", `3` = "eps=3", `4` = "eps=4", 
               `5` = "eps=5" )

ggplot(data = out_df, aes(x= timeseries)) +
  geom_area(aes(y= FO.1 + MZ.1 + GC.2 + MZ.2, fill="FO CD45.1" ), alpha=0.6, colour="white") +
  geom_area(aes(y= MZ.1 + GC.2 + MZ.2, fill="MZ CD45.1" ), alpha=0.6, colour="white") +
  geom_area(aes(y= GC.2 + MZ.2, fill="GC CD45.2" ), alpha=0.5, colour="white") +
  geom_area(aes(y= MZ.2, fill="MZ CD45.2"), alpha=0.5, colour="white") +
  #scale_fill_viridis(discrete = T) +
  #theme_ipsum() + 
  scale_x_log10(limits = c(7, 60), breaks=c(7, 14, 28, 56)) + scale_y_log10() + 
  labs(x='Days since immunization', y=NULL)  +  theme_ipsum() +
  theme(legend.title = element_blank()) +
  guides(fill='none')+
  facet_wrap(~EPS, nrow = 3, labeller = as_labeller(plot_labl)) 



ggplot(out_df, aes(x=timeseries, y=ratio_MZ))+
  geom_line(col=4, size=1.25)+ ylim(0, 1)+
  scale_x_log10(limits = c(7, 60), breaks=c(7, 14, 28, 56)) +
  labs(x='Days since immunization', y=NULL)  + theme_grey()+
  facet_wrap(~EPS)
  


plot_df <- out_df %>%
  select(timeseries, MZ.1, MZ.2, EPS) %>%
  gather(-c(timeseries, EPS), key = "popl", value = "counts") 

aresf <- plot_df %>%
  group_by(timeseries, EPS) %>%
  mutate(n_tot = sum(counts), 
         percentage = counts/n_tot)

# Give a specific order:
aresf$popln <- factor(aresf$popl , levels=c("MZ 45.1", "MZ 45.2") )

ggplot(data = aresf, aes(x= (timeseries), y=percentage*100, fill=popl)) +
  geom_area(alpha=0.6, colour="white") +
  scale_fill_manual(values = c(4,2))+
  theme_ipsum() +
  scale_x_log10(limits = c(7, 60), breaks=c(7, 14, 28, 56)) +
  #xlim(7, 56) + #scale_y_log10() + 
  labs(x='Days since immunization', y=NULL)  + 
  facet_wrap(.~EPS, labeller = as_labeller(plot_labl)) +
  theme(legend.title = element_blank())







#############
stop()

boot_func <- function(eps){
  ## Sampling random numbers to select the row from posterior dist matrix
  set.seed(1456)
  rowvec <- sample(1:nrow(fit_ss), 300, replace = FALSE) %>% sort()
  rowvec <- 
  
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

pooled_df <- boot_func(0.01)
  
ggplot(pooled_df) +
  geom_line(aes(x=timeseries, y=median, col=key))+
  geom_ribbon(aes(x = timeseries, ymin = lb, ymax = ub, fill=key), alpha = 0.25) +
  scale_x_log10() + scale_y_log10(limits=c(1, 1e5))


sol_df <- pooled_df %>%
  select(-lb, -ub) %>%
  spread(key = key, value = median)
  


ggplot(data = sol_df, aes(x= timeseries)) +
  geom_area(aes(y= FoB+ GCB + MZB, fill="FO_counts" ), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= GCB + MZB, fill="GC_counts"), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= MZB, fill="CARMZ_WT"), alpha=0.6 , linewidth=.5, colour="white") +
  #scale_fill_viridis(discrete = T) +
  #theme_ipsum() + 
  scale_x_log10(limits = c(5, 35), breaks=c(7, 14, 35)) + scale_y_log10() + 
  labs(x='Days since immunization', y=NULL)  + 
  guides(fill="none")

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
  theme_ipsum() + 
  xlim(5, 35) + #scale_y_log10() + 
  labs(x='Days since immunization', y=NULL)  + 
  #facet_wrap(.~popln) + 
  guides(col="none")
  

