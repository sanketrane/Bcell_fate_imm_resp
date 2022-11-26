
## testing models from stan
library(tidyverse)
library(rstan)

stanmodel_file <- file.path("stan_models", "sourcevar_model2.stan")
expose_stan_functions(stanmodel_file)

set.seed(18497)
ts_seq <- round(sample(seq(1, 100, length.out = 100), 45, replace = FALSE)) %>% sort()
init_cond <- c(exp(11.5))
params <- c(0.02, 0.1, 14, 0.01, 18)
ode_sol <- solve_ODE_sys(ts_seq, init_cond, params)

stan_pred_df <- data.frame("time_seq" = ts_seq,
                           "y_pred" = matrix(unlist(ode_sol), nrow = length(ts_seq), byrow = TRUE))%>%
  mutate(total_counts = rnorm(length(ts_seq), y_pred, 3e4),
         theta_obs = rnorm(length(ts_seq), sapply(ts_seq, theta_spline, parms=params), 1.5e5)) %>%
  select(-contains("y_pred")) 


ggplot(stan_pred_df)+
  geom_point(aes(x=time_seq, y=theta_obs), col=2)+
  geom_point(aes(x=time_seq, y=total_counts), col=4)+
  scale_y_log10()


## Data to import in Stan
numObs <- length(ts_seq)
solve_time <- stan_pred_df$time_seq
total_counts <- stan_pred_df$total_counts
theta_obs <-stan_pred_df$theta_obs

# time sequence for predictions specific to age bins within the data
ts_pred <- seq(1, 100, length.out = 300)
numPred <- length(ts_pred)


rstan::stan_rdump(c("numObs",  "solve_time", 
                    "total_counts",   "theta_obs",
                    "ts_pred", "numPred"),
                  file = file.path('datafiles', paste0('sorcevar',".Rdump")))
