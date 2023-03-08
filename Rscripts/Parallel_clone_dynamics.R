## Simulations of clonal dynamics using the model fits to immunization data

## clear env and mem cache
rm(list = ls()); gc();

## loading libraries
library(tidyverse)
library(wesanderson)
library(parallel)
library(rstan)

### model specific details that needs to be change for every run
modelName <- "Branched_timeinflux"


Wes_pallete <- wesanderson::wes_palette('Darjeeling1', 75, type = "continuous")


NITER = 300
NUM_CLONES = 75

### extract posterior distribution  of parameters as a df
fit_ss <- read.csv(file = paste0("PostDF_", modelName, ".csv"))
mu_med <- median(fit_ss$mu)
lambda_med <- median(fit_ss$lambda_WT)


## Phenomenological function describing number of activated (CAR expressing) FoB cells varying with time
CAR_FoB <- function(Time){
  # spline fitted to FoB numbers 
  F0 = exp(0); B0 = exp(5);  n = 4 ;  X = 9.4 ;  q = 6;
  value = (B0 * Time^n) * (1 - ((Time^q)/((X^q) + (Time^q))));
  return(value)
}

### Phenomenological function describing number of activated (CAR expressing) FoB cells varying with time
#CAR_MZB <- function(Time){
#  # spline fitted to FoB numbers 
#  #F0 = exp(11.722278); B0 = exp(4.475064);  n = 4.781548 ;  X = 6.943644 ;  q = 5;
#  F0 = exp(5); B0 = exp(5.5);  n = 3 ;  X = 9.2 ;  q = 4;
#  value =  (B0 * Time^n) * (1 - ((Time^q)/((X^q) + (Time^q)))) #- exp(11.7);
#  return(value)

## Assuming power-law distribution for the antigen-specific B cell clones in the activated pool
## We are simulating for 30 different clones.
## The pool size of activated FoB cells varies with time. 
power_law_prob <- function(x, k, A, m, n){
  value = (A * x^-k)/((n^(-k+1) - m^(-k+1))/(-k+1))
  return(value)
}

power_law_prob(seq(1, 30), k=0.3, A=1, m=1, n=30)

PLdist_Clone_Time <- function(Time, nClones){
  ## power law distribution of clones varying with time
  m_dot=1; n_dot=nClones;
  k0= 0.1; alpha = 0.03
  x_vec <- seq(1, nClones, 1)
  if (Time >= 4) {
    A_dot = CAR_FoB(Time)
    k_dot = k0 + alpha * (Time -4);
    value = as.integer(power_law_prob(x = x_vec, k=k_dot, A=A_dot, m=m_dot, n=n_dot))
  } else {
    value = 0;
  }
  return(value)
}

## Cells that differentiate into MZ and GC phenotype within a short interval time dt
### from the activated clones sample cells with propensity mu * X for MZ and alpha(Time) * X for GC B cells,
### Where X is the pool size of activated FoB cells 
MZ_Clone_recruitment_Time <- function(Time, distrib="powerlaw", nClones){
  mu = mu_med;  # median value of recruitment rate
  
  clone_sample <-  if(distrib == "unif"){
    sample(uniform_dist(Time, nClones), mu * CAR_FoB(Time))
  } else {
    sample(PLdist_Clone_Time(Time, nClones), mu * CAR_FoB(Time))
  }
  return(sort(clone_sample))
}


singleRun <- function(nClones, distrib){
  ### !!!TSTEP needs to be very small otherwise the propensities are not accurate!!!
  TSTEP = 0.04; updated_time = 4; Tmax = 60; current_time = 0; nClones = nClones
  start_clone_dist <- MZ_time(Time = updated_time, distrib = distrib, nClones = nClones)
  
  initial_dist <- data.frame("CloneID" = as.factor(seq(1, nClones, 1)),
                             "T0" = rep(0, nClones))
  start_dist <- data.frame(table(start_clone_dist)) %>%
    rename(CloneID=start_clone_dist,
           !!paste0('TS_', updated_time) := Freq)
  
  out_dist <- initial_dist %>%
    left_join(start_dist, by="CloneID")
  
  time_vec <- c(updated_time)
  
  while(current_time < Tmax){
    
    current_time = round(updated_time + TSTEP, 2);
    lambda = lambda_med; persist = 1 - lambda; mu = mu_med;
    prob_loss = 1 - exp(-lambda * TSTEP)
    
    ## pre-existing clones
    poolsize <- length(start_clone_dist)
    
    ## update the pool
    binomial_persist <- rbinom(1, size = poolsize, prob = 1 - prob_loss)
    clone_persist <- sample(start_clone_dist, binomial_persist)
    
    ## influx
    prob_influx = 1 - exp(-mu * TSTEP)
    number_new_cells <- as.integer(prob_influx * CAR_FoB(current_time))
    parent_dist <- MZ_time(Time = current_time, distrib = distrib, nClones = nClones)
    new_clones <- sample(parent_dist, number_new_cells)
    
    ## update the pool
    final_clone_dist <- c(clone_persist, new_clones)
    current_poolsize <- length(final_clone_dist)
    
    ## output
    updated_dist <- data.frame(table(final_clone_dist)) %>%
      rename(CloneID=final_clone_dist,
             !!paste0('TS_', current_time) := Freq)
    
    out_dist <- out_dist %>%
      left_join(updated_dist, by="CloneID")
    
    ## reset the time an clone dist
    start_clone_dist <- final_clone_dist
    updated_time <- current_time
    
    time_vec <- append(time_vec, values = updated_time)
    
  }
  
  timeseries <- c()
  for (i in 1: length(time_vec)){
    ind_a = 1 + (i-1) * nClones
    ind_b = nClones*i
    timeseries[ind_a:ind_b] <- rep(time_vec[i], nClones)
  }
  
  single_run_plot <- out_dist %>%
    select(- T0) %>% replace(is.na(.), 0) %>%
    gather(-CloneID, key="Timesteps", value = "Clonefreq")  %>%
    bind_cols("Timeseries" = timeseries)%>%
    select(- Timesteps)
  
  return(single_run_plot)
}


BatchRun <- function(nIter, nClones, distrib){
  ## function to iterate
  MultiRun_Func <- function(nIter, nClones, distrib){
    new_run <- singleRun(nClones=nClones, distrib=distrib)
    res1 <- data.frame(new_run$Clonefreq)
    colnames(res1) <- paste0('Clonefreq_', nIter)
    return(res1)
  }
  
  InitRun <- singleRun(nClones=nClones, distrib=distrib) 
  
  ## iterate MultiRun_Func for nIter times using lapply 
  iter_seq <- seq(1, nIter)
  para_run <- data.frame(mclapply(iter_seq, MultiRun_Func, nClones=nClones, distrib=distrib,
                                  mc.cores = detectCores()-2))
  #para_run <- data.frame(lapply(iter_seq, MultiRun_Func, NUM_Clones=NUM_Clones, distrib=distrib))
  MultiRun_df <- InitRun %>%
    bind_cols(para_run) 
  return(MultiRun_df)
}

print("START BATCH RUN!!")

system.time(
  MultiRunPlot <- BatchRun(nIter = NITER, nClones = NUM_CLONES, distrib =  "powerlaw")
)

MultiRunPlot <- MultiRunPlot %>% 
  gather(-c(CloneID, Timeseries), key='key', value = 'value') %>%
  group_by(CloneID, Timeseries) %>% replace(is.na(.), 0) %>%
  summarize(lb = as.integer(quantile(value, probs = 0.045)),
            ub = as.integer(quantile(value, probs = 0.955)),
            median = quantile(value, probs = 0.5)) 

ggplot(MultiRunPlot)+
  geom_line(aes(x=Timeseries, y=median, col=CloneID)) +
  geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill=CloneID), alpha = 0.25) +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d()+ scale_fill_viridis_d()+
  guides(col='none', fill = 'none') +
  facet_wrap(.~ CloneID, nrow = 5)


ggsave("target_dist_facets.pdf", last_plot(), width = 6, height = 4.5, device = 'pdf')

ggplot(MultiRunPlot)+
  geom_line(aes(x=Timeseries, y=median, col=CloneID)) +
  geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill=CloneID), alpha = 0.25) +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d()+ scale_fill_viridis_d()+
  guides(col='none', fill = 'none') 

ggsave("target_dist.pdf", last_plot(), width = 6, height = 4.5, device = 'pdf')


MultiRun_clonefreq <- MultiRunPlot %>%
  group_by(Timeseries) %>%
  mutate(clone_freq = median/sum(median))

ggplot(MultiRun_clonefreq)+
  geom_line(aes(x=Timeseries, y=clone_freq, col=CloneID)) +
  #geom_ribbon(aes(x = Timeseries, ymin = clone_freq_LB, ymax = clone_freq_UB, fill=CloneID), alpha = 0.25) +
  guides(col='none', fill='none') +
  ylim(0, 0.15)

ggsave("target_clonefreq.pdf", last_plot(), width = 6, height = 4.5, device = 'pdf')


print("DONE!")








































