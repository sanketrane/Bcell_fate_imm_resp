### Simulations of clonal dynamics using the model fits to immunization data

## loading libraries
library(tidyverse)
library(wesanderson)

## Phenomenological function describing number of activated (CAR expressing) FoB cells varying with time
CAR_FoB <- function(Time){
  # spline fitted to FoB numbers 
  F0 = exp(11.722278); B0 = exp(4.475064);  n = 4.781548 ;  X = 6.943644 ;  q = 5;
  value = F0 + (B0 * Time^n) * (1 - ((Time^q)/((X^q) + (Time^q))));
  return(value)
}

## Assuming that antigen-specific B cell clones are distributed uniformly in the activated pool
## We are simulating for 40 different clones.
## The pool size of acivated FoB cells varies with time. 
uniform_dist <- function(Time){ 
  rdunif(as.integer(CAR_FoB(Time)), a = 1, b = 40)
}

ggplot()+
  geom_histogram(aes(x=uniform_dist(30)), binwidth = 1, fill=4)


## Cells that differentiate into MZ and GC phenotype within a short interval time dt
### from the activated clones sample cells with propensity mu * X for MZ and alpha(Time) * X for GC B cells,
### Where X is the pool size of acivated FoB cells s
MZ_time <- function(Time){
  mu= 0.000322374412043;
  clone_sample <- sample(uniform_dist(Time), mu * CAR_FoB(Time))
  return(sort(clone_sample))
}


singleRun <- function(nClones){
  time_vec <- c(4, 7, 10, 14, 18, 21, 24, 27, 30)
  
  MZ_bin <- data.frame("CloneID" = as.factor(seq(1, nClones, 1)),
                       "T0" = rep(0, nClones))
  
  for (i in 1:length(time_vec)){
    new_freq_df <- data.frame(table(MZ_time(time_vec[i]))) %>%
      rename(CloneID=Var1,
             !!paste0('TSTEP_', i) := Freq)
    MZ_bin <- MZ_bin %>%
      left_join(new_freq_df, by="CloneID")
  }           
  
  MZ_bin_plot <- MZ_bin %>%
    select(-T0) %>%
    gather(-CloneID, key="Timesteps", value = "Clonefreq") %>%
    mutate(Timeseries = ifelse(Timesteps == "TSTEP_1", 4,
                               ifelse(Timesteps == "TSTEP_2", 7,
                                      ifelse(Timesteps == "TSTEP_3", 14,
                                             ifelse(Timesteps == "TSTEP_4", 21,
                                                    ifelse(Timesteps == "TSTEP_5", 18,
                                                           ifelse(Timesteps == "TSTEP_6", 24,
                                                                  ifelse(Timesteps == "TSTEP_7", 27, 30)))))))) %>%
    select(-Timesteps)
}

#ggplot(single_run_plot)+
#  geom_point(aes(x=Timeseries, y=Clonefreq, col=CloneID)) +
#  guides(col='none') +
#  facet_wrap(.~ CloneID, nrow = 5)


BatchRun <- function(nIter, nClones){
  multiRun <- singleRun(nClones) 
  for (i in 1:nIter){
    new_run <- singleRun(nClones) 
    multiRun <- multiRun %>%
      mutate(!!paste0("Clonefreq_", i) := new_run$Clonefreq)
  }
  multiRun_plot <- multiRun %>%
    gather(-c(CloneID, Timeseries), key='key', value = 'value') %>%
    group_by(CloneID, Timeseries) %>% replace(is.na(.), 0) %>%
    summarize(lb = quantile(value, probs = 0.045),
              median = quantile(value, probs = 0.5),
              ub = quantile(value, probs = 0.955)) 
}

  
multiRun_plot <- BatchRun(30, 30)

Wes_pallete <- wesanderson::wes_palette('Darjeeling1', 30, type = "continuous")

ggplot(multiRun_plot)+
  geom_line(aes(x=Timeseries, y=median, col=CloneID)) +
  geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill=CloneID), alpha = 0.25) +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d()+ scale_fill_viridis_d()+
  guides(col='none', fill = 'none') +
  facet_wrap(.~ CloneID, nrow = 5)


GC_time <- function(Time){
  t0=4
  alpha = 0.1175329567; nu = 0.00536442741333333;
  alpha_t= alpha/(1 + exp(nu * (Time-t0)^2));
  
  sample(uniform_dist(Time), alpha_t * CAR_FoB(Time))
}


ggplot()+
  geom_histogram(aes(x=MZ_time(10)), binwidth = 1, fill=2)+
  geom_histogram(aes(x=GC_time(10)), binwidth = 1, fill=4)


