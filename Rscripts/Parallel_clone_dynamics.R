## Simulations of clonal dynamics using the model fits to immunization data

## clear env and mem cache
rm(list = ls()); gc();

## loading libraries
## Installing r pachages on the go:
if(!("tidyverse" %in% rownames(installed.packages())) ){
  install.packages("tidyverse", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}
if(!("wesanderson" %in% rownames(installed.packages())) ){
  install.packages("wesanderson", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}
if(!("parallel" %in% rownames(installed.packages())) ){
  install.packages("parallel", repos = "https://cloud.r-project.org/", dependencies=TRUE)
}

library(tidyverse)
library(wesanderson)
library(parallel)

### model specific details that needs to be change for every run
modelName <- "Branched_timeinflux1"

NITER = 300
NUM_CLONES = 75
T_MAX = 60


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

rUN_seed <- paste0("Run_", args[1])
#rUN_seed <- paste0("Run_", 3)
print(rUN_seed)

Wes_pallete <- wesanderson::wes_palette('Darjeeling1', NUM_CLONES, type = "continuous")

### extract posterior distribution  of parameters as a df
fit_ss <- read.csv(file = paste0("PostDF_", modelName, ".csv"))
mu_med <- median(fit_ss$mu)
lambda_med <- median(fit_ss$lambda_WT)
lambda_sd <- sd(fit_ss$lambda_WT)


## Phenomenological function describing number of activated (CAR expressing) FoB cells varying with time
CAR_FoB <- function(Time){
  # spline fitted to FoB numbers 
  #F0 = exp(0); B0 = exp(5);  n = 4 ;  X = 9.4 ;  q = 6;
  #value = (B0 * Time^n) * (1 - ((Time^q)/((X^q) + (Time^q))));
  F0 = exp(11.722278);  B0 = exp(4.475064);  n = 4.781548 ;  X = 6.943644 ;  q = 5;
  value = F0 + (B0 * Time^n) * (1 - ((Time^q)/((X^q) + (Time^q))))
  return(value)
}

## Assuming power-law distribution for the antigen-specific B cell clones in the activated pool
## We are simulating for 30 different clones.
## The pool size of activated FoB cells varies with time. 
power_law_prob <- function(x, k, A, m, n){
  value = (A * x^-k)/((n^(-k+1) - m^(-k+1))/(-k+1))
  return(value)
}

K_Time <- function(Time){
  ## power law distribution of clones varying with time
  if(rUN_seed == "Run_1"){
    k0= 0.3; alpha = 0.025
    k_dot = k0 - alpha * (Time - 0)
  } else if (rUN_seed== "Run_2"){
    k0= 0.3; alpha = 0.3
    k_dot = k0 - alpha * (Time - 9)
  } else if (rUN_seed == "Run_3"){
    k0= 0.15; alpha = 0.01
    k_dot = k0+ alpha * (Time - 0)
  }
  return(k_dot)
}

PLdist_Clone_Time <- function(Time, nClones){
  ## power law distribution of clones varying with time
  m_dot=1; n_dot=nClones;
  x_vec <- seq(1, nClones, 1)
  if (Time >= 0) {
    A_dot = CAR_FoB(Time)
    k_dot = K_Time(Time);
    value = as.integer(power_law_prob(x = x_vec, k=k_dot, A=A_dot, m=m_dot, n=n_dot))
  } else {
    value = 0;
  }
  return(value)
}


pl_df <- data.frame("d07" = PLdist_Clone_Time(7, NUM_CLONES),
                    "d14" = PLdist_Clone_Time(14, NUM_CLONES),
                    "d21" = PLdist_Clone_Time(21, NUM_CLONES),
                    "d28" = PLdist_Clone_Time(28, NUM_CLONES),
                    "d35" = PLdist_Clone_Time(35, NUM_CLONES),
                    "d56" = PLdist_Clone_Time(56, NUM_CLONES), 
                    "Clone.ID" = seq(1, 75)) %>%
  gather(-Clone.ID,  key='Time.Days', value='Clone.Size')

ggplot(pl_df)+
  geom_line(aes(x=log(Clone.ID), y=log(Clone.Size), col = Clone.ID)) +
  scale_color_viridis_c() + 
  labs(x = "log(Clone Rank)", y = "log(Clone Size)") +
  guides(col='none') + #scale_x_log10()+
  facet_wrap(~Time.Days)

ggplot(pl_df)+
  geom_line(aes(x=log(Clone.ID), y=log(Clone.Size), col = Time.Days), size=0.75) +
  #scale_color_viridis_d(option = "D") + 
  scale_color_discrete(name = "Time since \n immunization")+
  labs(x = "log(Clone Rank)", y = "log(Clone Size)") +
  theme(legend.position = c(0.9, 0.75), legend.background = element_blank(), legend.title.align=0.5)

ggsave(paste0(rUN_seed, "_parent_Hist.pdf"), last_plot(), width = 6, height = 4.5, device = 'pdf')


cloneid <-c()  ## start an empty vector
j=0
for (i in 1:NUM_CLONES){
  ind_a = j + 1
  j = j + length(seq(0, T_MAX,0.3))
  ## assortment\
  cloneid[ind_a:j] <- rep(i, length(seq(0, T_MAX,0.3))) 
}


pwvec <- data.frame(t(sapply(seq(0, T_MAX,0.3), PLdist_Clone_Time, nClones=NUM_CLONES))) %>%
  bind_cols("Timeseries" = seq(0, T_MAX,0.3)) %>%
  gather(-Timeseries, key='Clone_num', value="CloneSize") %>%
  bind_cols("CloneID" = cloneid)


Parent_clonefreq <- pwvec %>%
  group_by(Timeseries) %>%
  mutate(clone_freq = CloneSize/sum(CloneSize))

p1 <- ggplot(Parent_clonefreq)+
  geom_line(aes(x=Timeseries, y=clone_freq, col=Clone_num)) +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  labs(x = "Time since immunization (days)", y = "Clone Frequnecy") +
  guides(col='none') + scale_y_log10(limits=c(0.0001, 1)) #+ scale_x_log10(limits=c(3, 60))

ggsave(paste0(rUN_seed, "_parent_freq_dist.pdf"), p1, width = 6, height = 4.5, device = 'pdf')

#sapply(seq(1, 60), K_Time)

p2 <- ggplot(pwvec)+
  geom_line(aes(x=Timeseries, y=CloneSize, col=as.factor(CloneID))) +
  #scale_y_log10() + scale_x_log10() +
  labs(x = "Time since immunization (days)", y = "Clone size") +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d() + +
  #facet_wrap(.~ CloneID, nrow = 5)
  guides(col='none')

ggsave(paste0(rUN_seed, "_parent_count_dist.pdf"), p2, width = 6, height = 4.5, device = 'pdf')

### ditribution of antigen-specific clones in FoB population
PL_parent_dist <- function(Time, nClones){
  ## generate a distribution
  dist <-c()  ## start an empty vector
  j=0
  for (i in 1:nClones){
    npl <- PLdist_Clone_Time(Time, nClones)
    ### indices to distribute frequencies to each clone ID
    ind_a = j + 1
    j = j + npl[i]
    ## assortment
    dist[ind_a:j] <- rep(i, npl[i]) 
  }
  return(dist)
}


singleRun <- function(nClones){
  ### !!!TSTEP needs to be very small otherwise the propensities are not accurate!!!
  TSTEP = 0.04; updated_time = 1; Tmax = T_MAX; current_time = 0; nClones = nClones
  errp <- try(MZ_Clone_dist_Time(Time = updated_time, nClones = nClones), silent = T)
  while("try-error" %in% class(errp)) {
    updated_time = updated_time +TSTEP
    errp <- try(MZ_Clone_dist_Time(Time = updated_time, nClones = nClones), silent = T)
  }
  start_clone_dist <- MZ_Clone_dist_Time(Time = updated_time, nClones = nClones) %>% sort()
  
  while(length(start_clone_dist)<2){
    updated_time = updated_time +TSTEP
    start_clone_dist <- MZ_Clone_dist_Time(Time = updated_time, nClones = nClones) %>% sort()
  } 
  
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
    lambda_Gauss = rnorm(nClones, lambda_med, lambda_sd)
    prob_loss_Gauss = 1 - exp(-lambda_Gauss * TSTEP)
    prob_loss = 1 - exp(-lambda * TSTEP)
    prob_influx = 1 - exp(-mu * TSTEP)
    
    ## pre-existing clones
    poolsize <- length(start_clone_dist)
    
    ## update the pool
    binomial_persist <- rbinom(1, size = poolsize, prob = 1 - prob_loss)
    #clone_persist <- sample(start_clone_dist, binomial_persist, prob = prob_loss_Gauss[start_clone_dist]) %>% sort()
    clone_persist <- sample(start_clone_dist, binomial_persist) %>% sort()
    
    ## influx
    size_of_recruit <- rbinom(1, size = as.integer(CAR_FoB(current_time)), prob = prob_influx)
    new_clones <- sample(PL_parent_dist(current_time, nClones), size_of_recruit)
    
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
    start_clone_dist <- final_clone_dist %>% sort()
    updated_time <- current_time
    if (updated_time %% 5 == 0) print(updated_time)
    
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

#single_run_plot <- singleRun(NUM_CLONES)
#
#ggplot(single_run_plot)+
#  geom_line(aes(x=Timeseries, y=Clonefreq, col=as.factor(CloneID))) +
#  #scale_y_log10() + scale_x_log10() +
#  labs(x = "Time since immunization (days)", y = "Clone size") +
#  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete) +
#  #geom_text(aes(x=Timeseries, y=Clonefreq, label=CloneID, col=as.factor(CloneID)))+
#  #scale_color_viridis_d() + +
#  #facet_wrap(.~ CloneID, nrow = 5)
#  guides(col='none')
#
#ggsave("MZ_count_dist.pdf", last_plot(), width = 6, height = 4.5, device = 'pdf')
#
#ggplot(single_run_plot)+
#  geom_line(aes(x=Timeseries, y=Clonefreq, col=CloneID)) +
#  labs(x = "Time since immunization (days)", y = "Clone size") +
#  guides(col='none') +  
#  facet_wrap(.~ CloneID, nrow = 5)
#
#
#ggsave("MZ_clone_dist.pdf", last_plot(), width = 10, height = 9, device = 'pdf')
#
#singleRun_clonefreq <- single_run_plot %>%
#  group_by(Timeseries) %>%
#  mutate(clone_freq = Clonefreq/sum(Clonefreq))
#
#ggplot(singleRun_clonefreq)+
#  geom_line(aes(x=Timeseries, y=clone_freq, col=CloneID)) +
#  labs(x = "Time since immunization (days)", y = "Clone Frequency") +
#  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
#  guides(col='none') + scale_y_log10(limits=c(0.0001, 1))
#
#
#ggsave("MZ_freq_dist.pdf", last_plot(), width = 6, height = 4.5, device = 'pdf')

BatchRun <- function(nIter, nClones){
  ## function to iterate
  MultiRun_Func <- function(nIter, nClones){
    new_run <- singleRun(nClones=nClones)
    res1 <- data.frame(new_run$Clonefreq)
    colnames(res1) <- paste0('Clonefreq_', nIter)
    return(res1)
  }
  
  InitRun <- singleRun(nClones=nClones) 
  
  ## iterate MultiRun_Func for nIter times using lapply 
  iter_seq <- seq(1, nIter)
  para_run <- data.frame(mclapply(iter_seq, MultiRun_Func, nClones=nClones,
                                  mc.cores = detectCores()-2))
  #para_run <- data.frame(lapply(iter_seq, MultiRun_Func, NUM_Clones=NUM_Clones, distrib=distrib))
  MultiRun_df <- InitRun %>%
    bind_cols(para_run) 
  return(MultiRun_df)
}

print("START BATCH RUN!!")

system.time(
  MultiRunPlot <- BatchRun(nIter = NITER, nClones = NUM_CLONES)
)

MultiRunPlot <- MultiRunPlot %>% 
  gather(-c(CloneID, Timeseries), key='key', value = 'value') %>%
  group_by(CloneID, Timeseries) %>% replace(is.na(.), 0) %>%
  summarize(lb = as.integer(quantile(value, probs = 0.045)),
            ub = as.integer(quantile(value, probs = 0.955)),
            median = quantile(value, probs = 0.5)) 

p4 <- ggplot(MultiRunPlot)+
  geom_line(aes(x=Timeseries, y=median, col=CloneID)) +
  geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill=CloneID), alpha = 0.25) +
  labs(x = "Time since immunization (days)", y = "Clone size") +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d()+ scale_fill_viridis_d()+
  guides(col='none', fill = 'none') +
  facet_wrap(.~ CloneID, nrow = 5)


ggsave(paste0(rUN_seed, "_target_dist_facets.pdf"), p4, width = 10, height = 9, device = 'pdf')

p5 <- ggplot(MultiRunPlot)+
  geom_line(aes(x=Timeseries, y=median, col=CloneID), alpha = 0.75) +
  geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill=CloneID), alpha = 0.25) +
  labs(x = "Time since immunization (days)", y = "Clone size") +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d()+ scale_fill_viridis_d()+
  guides(col='none', fill = 'none') 

ggsave(paste0(rUN_seed, "_target_dist.pdf"), p5, width = 6, height = 4.5, device = 'pdf')


MultiRun_clonefreq <- MultiRunPlot %>%
  group_by(Timeseries) %>%
  mutate(clone_freq = median/sum(median),
         clone_freq_LB = lb/sum(lb),
         clone_freq_UB = ub/sum(ub))

p6 <- ggplot(MultiRun_clonefreq)+
  geom_line(aes(x=Timeseries, y=clone_freq, col=CloneID)) +
  geom_ribbon(aes(x = Timeseries, ymin = clone_freq_LB, ymax = clone_freq_UB, fill=CloneID), alpha = 0.25) +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  labs(x = "Time since immunization (days)", y = "Clone Frequency") +
  guides(col='none', fill='none') +scale_y_log10(limits=c(0.0001, 1)) #+ scale_x_log10(limits=c(1, 60))

ggsave(paste0(rUN_seed, "_target_clonefreq.pdf"), p6, width = 6, height = 4.5, device = 'pdf')


print("DONE!")








































