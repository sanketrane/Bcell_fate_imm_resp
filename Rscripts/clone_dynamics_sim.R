## Simulations of clonal dynamics using the model fits to immunization data

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
## We are simulating for 30 different clones.
## The pool size of activated FoB cells varies with time. 
uniform_dist <- function(Time, nClones){ 
  rdunif(as.integer(CAR_FoB(Time)), a = 1, b = nClones)
}

ggplot()+
  geom_histogram(aes(x=uniform_dist(10, 30)), binwidth = 1, fill=4)

## Assuming power-law distribution for the antigen-specific B cell clones in the activated pool
## We are simulating for 30 different clones.
## The pool size of activated FoB cells varies with time. 
Norm_power_law <- function(Time, nClones){
  k= 0.5; a=1; b=nClones;
  
  x_vec <- seq(1, nClones, 1)
 as.integer(CAR_FoB(Time) * (x_vec^-k)/((b^(-k+1) - a^(-k+1))/(-k+1)))
}

power_law_dist <- function(Time, nClones){
  ## generate a distribution
  dist <-c()  ## start an empty vector
  j=0
  for (i in 1:nClones){
    npl <- Norm_power_law(Time, nClones)
    ### indices to distribute frequencies to each clone ID
    ind_a = j + 1
    j = j + npl[i]
    ## assortment
    dist[ind_a:j] <- rep(i, npl[i]) 
  }
  return(dist)
}

pl_dist <- power_law_dist(10, 30)

ggplot()+
  geom_histogram(aes(x=pl_dist), binwidth = 1, fill=2)



## Cells that differentiate into MZ and GC phenotype within a short interval time dt
### from the activated clones sample cells with propensity mu * X for MZ and alpha(Time) * X for GC B cells,
### Where X is the pool size of activated FoB cells 
MZ_time <- function(Time, distrib="powerlaw", nClones){
  mu= 0.000322374412043;
  
  clone_sample <-  if(distrib == "unif"){
    sample(uniform_dist(Time, nClones), mu * CAR_FoB(Time))
  } else {
    sample(power_law_dist(Time, nClones), mu * CAR_FoB(Time))
  }
  return(sort(clone_sample))
}



singleRun <- function(nClones, distrib="powerlaw"){
  time_vec <- c(4, 7, 10, 14, 18, 21, 24, 27, 30)
  
  MZ_bin <- data.frame("CloneID" = as.factor(seq(1, nClones, 1)),
                       "T0" = rep(0, nClones))
  
  for (i in 1:length(time_vec)){
    new_freq_df <- data.frame(table(MZ_time(Time = time_vec[i], distrib = distrib, nClones = nClones))) %>%
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

single_run_plot <- singleRun(30)

ggplot(single_run_plot)+
  geom_point(aes(x=Timeseries, y=Clonefreq, col=CloneID)) +
  guides(col='none') +
  facet_wrap(.~ CloneID, nrow = 5)


BatchRun <- function(nIter, distrib="powerlaw", nClones){
  multiRun <- singleRun(nClones) 
  for (i in 1:nIter){
    new_run <- singleRun(nClones = nClones, distrib = distrib) 
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

  
multiRun_plot <- BatchRun(30, distrib="unif", 30)

Wes_pallete <- wesanderson::wes_palette('Darjeeling1', 30, type = "continuous")

ggplot(multiRun_plot)+
  geom_line(aes(x=Timeseries, y=median, col=CloneID)) +
  geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill=CloneID), alpha = 0.25) +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d()+ scale_fill_viridis_d()+
  guides(col='none', fill = 'none') +
  facet_wrap(.~ CloneID, nrow = 5)


TSTEP = 0.1; updated_time = 4; Tmax = 10; current_time = 0; NUM_Clones = 30
start_clone_dist <- MZ_time(Time = updated_time, distrib = 'unif', nClones = NUM_Clones)

initial_dist <- data.frame("CloneID" = as.factor(seq(1, NUM_Clones, 1)),
                          "T0" = rep(0, NUM_Clones))
start_dist <- data.frame(table(start_clone_dist)) %>%
  rename(CloneID=start_clone_dist,
         !!paste0('TS_', updated_time) := Freq)

out_dist <- initial_dist %>%
  left_join(start_dist, by="CloneID")

time_vec <- c(updated_time)

while(current_time <= Tmax){
  
  current_time = round(updated_time + TSTEP, 2);
  lambda = 0.3; persist = 1 - lambda; mu = 0.00032;
  prob_loss = 1- exp(-lambda * TSTEP)
 
  ## pre-existing clones
  poolsize <- length(start_clone_dist)
  
  ## update the pool
  binomial_loss <- rbinom(1, size = poolsize, prob = prob_loss)
  clone_persist <- sample(start_clone_dist, poolsize - binomial_loss)
  
  ## influx
  number_new_cells <- as.integer(mu * CAR_FoB(current_time) * TSTEP)
  parent_dist <- MZ_time(Time = current_time, distrib = 'unif', nClones = NUM_Clones)
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

single_run_plot <- out_dist %>%
  select(- T0) %>%
  gather(-CloneID, key="Timesteps", value = "Clonefreq")  %>%
  bind_cols("Timeseries" = rep(time_vec, NUM_Clones))


ggplot(single_run_plot)+
  geom_line(aes(x=Timeseries, y=Clonefreq, col=CloneID)) +
  guides(col='none') +
  facet_wrap(.~ CloneID, nrow = 5)
  
  
  
  
  
  
  
  
  
  
  



























