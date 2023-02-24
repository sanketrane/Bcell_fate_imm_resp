library(deSolve)
library(tidyverse)
library(readxl)
library(viridis)
library(hrbrthemes)


#### donor CAR+ MZ and GC B cells in recipients
## ode function for Branched model

CAR_negative_MZB <- function(Time){
  M0 = exp(14.06);  nu = 0.0033;  b0 = 20.58;
  value = M0 * (1 + exp(-nu * (Time - b0)^2));
  return(value)
}

ode_func <-  function(t, state, parms){
  with(as.list(c(state, parms)),{
    alpha = parms[1]
    beta = parms[2]
    mu = parms[3]
    delta = parms[4]
    lambda_WT = parms[5]
    lambda_N2KO = parms[6]
    nu = parms[7]
    epsilon = parms[8]
    
    t0 = 4.0;
    alpha_tau = alpha/(1 + exp(nu * (t-t0)^2))
    beta_tau = beta/(1 + exp(nu * (t-t0)^2))
    
    #Fo b cells
    dY1 <- - (alpha_tau + mu + epsilon) * Y1
    
    # CAR positive GCB cells in WT
    dY2 = alpha_tau * Y1 - delta * Y2;
    
    # CAR positive MZB cells in WT
    dY3 = mu * Y1 - lambda_WT * Y3;
    
    #return the rate of change
    list(c(dY1, dY2, dY3))
    
  })  # end with(as.list ...
}

#initial conditions
state <- c(Y1 = 1e6, Y2 = 0, Y3=0)

# time points for which conc is reported
# include the points where data is available
times <- seq(4, 35, length.out=150)
parms_vec <- c(alpha=0.049, beta=0.0048, mu=0.0025, delta=1.05, lambda_WT=0.81, lambda_N2KO=0.73, nu=-0.68, epsilon=0.2)

sol_ode <- ode(y= state, times=times, func = ode_func, parms = parms_vec)

sol_df <-  data.frame(sol_ode) %>%
  rename(FO_counts = Y1,
         CARMZ_WT = Y3,
         GC_counts = Y2,
         time_since_imm = time) %>%
  filter(time_since_imm >= 5)


ggplot(data = sol_df, aes(x= time_since_imm)) +
  geom_area(aes(y= FO_counts+ GC_counts + CARMZ_WT, fill="FO_counts" ), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= GC_counts + CARMZ_WT, fill="GC_counts"), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= CARMZ_WT, fill="CARMZ_WT"), alpha=0.6 , linewidth=.5, colour="white") +
  #scale_fill_viridis(discrete = T) +
  #theme_ipsum() + 
  scale_x_log10(limits = c(5, 35), breaks=c(7, 14, 35)) + scale_y_log10() + 
  labs(x='Days since immunization', y=NULL)  + 
  guides(fill="none")

plot_df <- sol_df %>% 
  gather(-time_since_imm, key = 'popln', value = 'counts') %>%
  filter(popln != "FO_counts") %>%
  filter(time_since_imm >= 5)

aresf <- plot_df %>%
  group_by(time_since_imm) %>%
  mutate(n_tot = sum(counts), 
         percentage = counts/n_tot)

# Give a specific order:
aresf$popln <- factor(aresf$popln , levels=c("CARMZ_WT", "GC_counts", "FO_counts") )

ggplot(data = aresf, aes(x= (time_since_imm), y=percentage, fill=popln)) +
  geom_area(alpha=0.6 , linewidth=.5, colour="white") +
  scale_fill_viridis(discrete = T) +
  theme_ipsum() + 
  xlim(5, 95) + #scale_y_log10() + 
  labs(x='Days since immunization', y=NULL)  + 
  #facet_wrap(.~popln) + 
  guides(col="none")
  



#### data transformation functions
logit_trans <- function(x){
  
  #log(x/(1-x))
  asin(sqrt(x))
}

expit_trans <- function(x){
  
  exp(x)/(1 + exp(x))
}


