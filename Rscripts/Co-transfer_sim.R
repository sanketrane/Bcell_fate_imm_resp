library(deSolve)
library(tidyverse)
library(hrbrthemes)
library(rstan)

## model specific details that needs to be change for every run
modelName1 <- "Branched_timeinflux1"
modelName2 <- "Linear_timeinflux1"

theme_set(theme_light())
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12), axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(), legend.key = element_blank(),
                 legend.text = element_text(size=12), legend.title = element_text(12),
                 strip.text.x = element_text(size = 14))


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

log10minorbreaks=as.numeric(1:10 %o% 10^(1:8))


### extract posterior distribution  of parameters as a df
fit_ss <- read.csv(file = paste0("PostDF_", modelName1, ".csv"))
mu_b <- median(fit_ss$mu)
lambda_b <- median(fit_ss$lambda_WT)
alpha_b <- median(fit_ss$alpha)
nu_b <- median(fit_ss$nu)
delta_b <- median(fit_ss$delta)


fit_ss2 <- read.csv(file = paste0("PostDF_", modelName2, ".csv"))
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
    
    #Fo B cells cd45.1
    ## epsilon denotes the loss of Activated FO B cells in addition to their conversion to GC and MZ
    dY1 <- - (alpha_tau + mu_branched + epsilon) * Y1
    
    # CAR positive GCB cells cd45.1
    dY2 = alpha_tau * Y1 - (delta + mu_linear)  * Y2;
    
    # CAR positive MZB CD45.1
    ##influx into CAR+MZ is constant -- fitted better to immunization data!
    dY3 = mu_branched * Y1 + mu_linear * Y2 - lambda_linear * Y3;
    
    # CAR positive GCB CD45.2
    dY4 = - (delta + mu_linear) * Y4;
    
    # CAR positive MZB CD45.2
    dY5 = mu_linear * Y4 - lambda_linear * Y5;
    
    #return the rate of change
    list(c(dY1, dY2, dY3, dY4, dY5))
    
  })  # end with(as.list ...
}

#initial conditions
state <- c(Y1 = 1e5, Y2 = 0,  Y3=0, Y4=1e5, Y5=0)

# time points for which conc is reported
# include the points where data is available
ts_pred <- seq(7, 60, length.out=100)


out_df <- data.frame()
eps_vec <- c(0.01, 0.1, 0.7, 1)

for(i in 1:length(eps_vec)) {
  parms_vec <- c("alpha" = alpha_b, "mu_branched" = mu_b , "mu_linear" = mu_l, 
                 "delta" = delta_l, "lambda_branched" = lambda_b, "lambda_linear" = lambda_l,
                 "nu" = nu_b, "epsilon" = eps_vec[i])  ## epsilon denotes the loss of Activated FO B cells in addition to their conversion to GC and MZ
  ### simulation using the selected params
  sol_ode <- data.frame(ode(y= state, times=ts_pred, func = ode_func, parms = parms_vec))
  
  sol_df <- sol_ode %>%
    rename(timeseries = time,
           FO.1 = Y1,
           GC.1 = Y2,
           MZ.1 = Y3,
           GC.2 = Y4,
           MZ.2 = Y5) %>%
    mutate(ratio_MZ = replace_na(MZ.2/MZ.1, 0),
           EPS = round(eps_vec[i], 4))
  out_df <- rbind(out_df, sol_df)
  out_df
}

plot_df <- out_df %>% 
  mutate(FO.1 = replace(FO.1, FO.1 <= 1, 1),
         GC.1 = replace(GC.1, GC.1 <= 1, 1),
         MZ.1 = replace(MZ.1, MZ.1 <= 1, 1),
         GC.2 = replace(GC.2, GC.2 <= 1, 1),
         MZ.2 = replace(MZ.2, MZ.2 <= 1, 1))

plot_df$EPS <- factor(plot_df$EPS, levels = c("0.01", "0.1", "0.7", "1"), ordered = TRUE,
                      labels = c(expression(paste(phi, "=0.01")), expression(paste(phi, "=0.1")),
                                 expression(paste(phi, "=0.7")), expression(paste(phi, "=1"))))


#ggplot(data = out_df, aes(x= timeseries)) +
#  geom_area(aes(y= GC.1 + MZ.1 + GC.2 + MZ.2 , fill="GC CD45.1" ), alpha=0.6, colour="white") +
#  geom_area(aes(y= GC.2 + MZ.2 , fill="GC CD45.2" ), alpha=0.5, colour="white") +
#  geom_area(aes(y= MZ.2, fill="MZ CD45.2"), alpha=0.5, colour="white") +
#  geom_area(aes(y= MZ.1 + GC.2 + MZ.2 , fill="MZ CD45.1" ), alpha=0.6, colour="white") +
#  scale_x_log10(limits = c(7, 60), breaks=c(7, 14, 28, 56)) + 
#  scale_y_continuous(limits = c(1, 1e6), trans = "log10", labels = fancy_scientific, minor_breaks = log10minorbreaks) + 
#  labs(x='Days since immunization', y=NULL)  +  #theme_ipsum() +
#  guides(fill='none')+
#  facet_wrap(~EPS, nrow = 3, labeller = as_labeller(plot_labl)) + myTheme

ggplot(data = plot_df, aes(x= timeseries)) +
  geom_area(aes(y= GC.1, fill="GC CD45.1" ), alpha=0.6, colour="white") +
  geom_area(aes(y= GC.2, fill="GC CD45.2" ), alpha=0.5, colour="white") +
  geom_area(aes(y= MZ.2, fill="MZ CD45.2"), alpha=0.5, colour="white") +
  geom_area(aes(y= MZ.1, fill="MZ CD45.1" ), alpha=0.6, colour="white") +
  scale_x_log10(limits = c(7, 60), breaks=c(7, 14, 28, 56)) + 
  scale_y_continuous(limits = c(1, 2e5), trans = "log10", labels = fancy_scientific, minor_breaks = log10minorbreaks) + 
  labs(x='Days since immunization', y=NULL)  +  #theme_ipsum() +
  guides(fill='none') +
  facet_wrap(~EPS, nrow = 3, labeller = label_parsed) + myTheme
  


plot_df2 <- out_df %>%
  select(timeseries, MZ.1, MZ.2, EPS) %>%
  gather(-c(timeseries, EPS), key = "popl", value = "counts") 

aresf <- plot_df2 %>%
  group_by(timeseries, EPS) %>%
  mutate(n_tot = sum(counts), 
         percentage = counts/n_tot)

aresf$EPS <- factor(aresf$EPS, levels = c("0.01", "0.1", "0.7", "1"), ordered = TRUE,
                      labels = c(expression(paste(phi, "=0.01")), expression(paste(phi, "=0.1")),
                                 expression(paste(phi, "=0.7")), expression(paste(phi, "=1"))))

# Give a specific order:
aresf$popln <- factor(aresf$popl , levels=c("MZ 45.1", "MZ 45.2") )

ggplot(data = aresf, aes(x= (timeseries), y=percentage*100, fill=popl)) +
  geom_area(alpha=0.6, colour="white") +
  scale_fill_manual(values = c(4,2))+
  scale_x_log10(limits = c(7, 60), breaks=c(7, 14, 28, 56)) +
  #xlim(7, 56) + #scale_y_log10() + 
  labs(x='Days since immunization', y=NULL)  + 
  facet_wrap(.~EPS, labeller = label_parsed) +
  guides(fill='none') + 
  theme(legend.title = element_blank()) + myTheme






  

