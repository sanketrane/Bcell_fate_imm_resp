
## testing models from stan
library(tidyverse)
library(rstan)

stanmodel_file <- file.path("stan_models", "sourcevar_model2.stan")
expose_stan_functions(stanmodel_file)


ts_seq <- c(7, 14, 30)
init_cond <- c(exp(11.5))
params <- c(0.05, 0.1)
ode_sol <- solve_ODE_sys(ts_seq, init_cond, params)

ts_pred <- seq(1.1, 30, length.out=100)
ode_df <- solve_ODE_sys(ts_pred, init_cond, params)
stan_pred_df <- data.frame("time_pred" = ts_pred,
                           "y_pred" = matrix(unlist(ode_df), nrow = length(ode_df), byrow = TRUE))%>%
  mutate(CAR_GC = y_pred) %>%
  select(time_pred, contains("CAR")) %>%
  gather(-time_pred, key="subplop", value="counts")


ggplot(stan_pred_df)+
  geom_point(aes(x=time_pred, y=counts, col=subplop))+
  scale_y_log10()+
  facet_grid(.~subplop)



