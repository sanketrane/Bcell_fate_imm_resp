library(deSolve)
library(tidyverse)
library(readxl)

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
    
    t0 = 4.0;
    alpha_tau = alpha/(1 + exp(nu * (t-t0)^2))
    beta_tau = beta/(1 + exp(nu * (t-t0)^2))
    
    #Fo b cells
    dY1 <- - (alpha_tau + mu + epsilon) * Y1
    
    # CAR positive GCB cells in WT
    dY2 = alpha_tau * Y1 - delta * Y2;
    
    # CAR positive MZB cells in WT
    dY3 = mu * Y1 + beta_tau * CAR_negative_MZB(t) - lambda_WT * Y3;
    
    # CAR positive MZB cells in N2KO
    dY4 = beta_tau * CAR_negative_MZB(t) - lambda_N2KO * Y4;
    
    #return the rate of change
    list(c(dY1, dY2, dY3, dY4))
    
  })  # end with(as.list ...
}

#initial conditions
state <- c(Y1 = 1e6, Y2 = 0, Y3=0, Y4=0)

# time points for which conc is reported
# include the points where data is available
times <- c(0, 1, 4, 9, 14)
parms_vec <- c(alpha=0.049, beta=0.0048, mu=0.0025, delta=1.05, lambda_WT=0.81, lambda_N2KO=0.73, nu=-0.68)

sol_ode <- ode(y= state, times=times, func = ode_func, parms = parms_vec)

sol_df <-  data.frame(sol_ode) %>%
  mutate(FO_counts = Y1,
         CARMZ_counts = Y3,
         GC_counts = Y2) #%>% na.omit

ode_solver(solve_t, par_start)

## conversion of Adoptively transferred of FO B cells into MZ B cells
adopt_tra_df <- read_excel("Datafiles/TAM_all_samples.xlsx", sheet =4) %>%na.omit() %>%
  select(-contains('host')) %>%
  mutate(donor_prop = donor_in_Bcells/100,
         mzbdonor_prop = (mzb_in_donor/100),
         fobdonor_prop = (fob_in_donor/100)) %>%
  select(days.post.transfer, contains('prop')) 

#### data transformation functions
logit_trans <- function(x){
  
  #log(x/(1-x))
  asin(sqrt(x))
}

expit_trans <- function(x){
  
  exp(x)/(1 + exp(x))
}



#fitting log N_total & logit ki_prop to the transformed NCD4 data
LL_fit <- function(params, boot_data) { 
  
  ## transforming data
  observed_data.df <- boot_data%>%
    mutate(logit_mz = logit_trans(mzbdonor_prop),
           logit_fo = logit_trans(fobdonor_prop)) 
  
  ## logit transformed poroportions from the model
  pred_mz_prop <- logit_trans(ode_solver(observed_data.df$days.post.transfer, params,1))
  pred_fo_prop <- logit_trans(ode_solver(observed_data.df$days.post.transfer, params,2))
  
  ## calculating the sun of squared residuals
  ssqres1 <- sum((pred_mz_prop - observed_data.df$logit_mz)^2)
  ssqres2 <- sum((pred_fo_prop - observed_data.df$logit_fo)^2)
  
  k  <- length(params)                #number of unknown parameters 
  n1 <- nrow(observed_data.df)            #number of observations in dataset1
  
  #cost function
  #log-likelihood ignoring all the terms dependent only on the number of observations n
  #matrix multipltication of residual and transpose of residuals
  logl <-  - (n1/2)* log(ssqres1)   - (n1/2)* log(ssqres2)  
  
  return(-logl)     #since optim minimizes the function by default, ML
} 


### optimising the LL function to maximise the loglikelihood
optim_fit <- optim(par = par_start, fn = LL_fit, boot_data = adopt_tra_df,
                 method=c("Nelder-Mead"),
                 hessian = T,  control = list(trace = 1, maxit = 1000))


## calculating AIC 
AIC_est <- 2 * length(optim_fit$par) + 2 * optim_fit$value


## Parametere estimates from the best fit method
par_est <- c("delta" = optim_fit$par[1], 
             "mu"  = optim_fit$par[2],
             "lambda" = optim_fit$par[3])

pred_time <- seq(0, 15, 0.1)
pred_mz_prop <- ode_solver(pred_time, par_est,1)
pred_fo_prop <- ode_solver(pred_time, par_est,2)


ggplot() +
  geom_point(data = adopt_tra_df, aes(x= (days.post.transfer), y=mzbdonor_prop)) +
  geom_line(aes(x=pred_time, y=pred_mz_prop)) +
  xlim(0, 15) +
  #scale_y_log10() + 
  labs(x='Days post trabsfer', y=NULL) +
  myTheme + guides(col="none")

ggplot() +
  geom_point(data = adopt_tra_df, aes(x= (days.post.transfer), y=fobdonor_prop)) +
  geom_line(aes(x=pred_time, y=pred_fo_prop)) +
  xlim(0, 15) + ylim(0.5,1) + 
  labs(x='Days post trabsfer', y=NULL) +
  myTheme + guides(col="none")

1/exp(par_est)





