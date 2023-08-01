## Simulations of clonal dynamics using the model fits to immunization data

## clear env and mem cache
rm(list = ls()); gc();

# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly = TRUE)
# 
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).\n", call.=FALSE)
# } else if (length(args)==1) {
#   # default output file
#   args[3] = "out.txt"
# }

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

NITER = 100
NUM_CLONES = 75
T_MAX = 60



#rUN_seed <- paste0("Run_", args[1])
rUN_seed <- paste0("Run_", 2)
print(rUN_seed)
#simnum <- args[2]

Wes_pallete <- wesanderson::wes_palette('Darjeeling1', NUM_CLONES, type = "continuous")

### extract posterior distribution  of parameters as a df
fit_ss <- read.csv(file = paste0("PostDF_", modelName, ".csv"))
mu_med <- median(fit_ss$mu)
lambda_med <- median(fit_ss$lambda_WT)
lambda_sd <- sd(fit_ss$lambda_WT)

## Phenomenological function describing number of activated (CAR expressing) FoB cells varying with time
CAR_FoB <- function(Time){
  # spline fitted to FoB numbers
  F0 = exp(11.722279);  B0 = exp(4.475058);  n = 4.781548 ;  X = 6.943657 ;  q = 5.838882;
  #F0 = exp(11.2);  B0 = exp(7.5);  n = 3 ;  X = 10 ;  q = 5;
  value = F0 + (B0 * Time^n) * (1 - ((Time^q)/((X^q) + (Time^q))))
  # F0 = exp(16.7);  nu = 0.004;  b0 = 20;
  # value = F0 * (1 + exp(-nu * (Time - b0)^2));
  return(value)
}
ts_pred <- seq(0, T_MAX,0.3)
fobvec <- sapply(ts_pred, CAR_FoB)
ggplot()+
  geom_line(aes(x=ts_pred, y=fobvec)) + scale_y_log10()

## Assuming power-law distribution for the antigen-specific B cell clones in the activated pool
## We are simulating for 30 different clones.
## The pool size of activated FoB cells varies with time. 
power_law_prob <- function(x, k, nclones){
  value = (x^-k * (1-k))/nclones^(1-k)
  return(value)
}


K_Time <- function(Time){
  ## power law distribution of clones varying with time
  if(rUN_seed == "Run_1"){
    k0= 0.2; alpha = 0.005
    k_dot = k0 + alpha * (Time - 0)
  } else if (rUN_seed== "Run_2"){
    k0= 0.15; alpha = 0.013
    k_dot = k0 - alpha * (Time - 0)
  } 
  return(k_dot)
}

PLdist_Clone_Time <- function(Time, nClones){
  ## power law distribution of clones varying with time
  n_dot=nClones;
  x_vec <- seq(1, nClones, 1)
  if (Time >= 0) {
    A_dot = CAR_FoB(Time)
    k_dot = K_Time(Time);
    value = as.integer(A_dot * power_law_prob(x = x_vec, k=k_dot, n=n_dot))
  } else {
    value = 0;
  }
  return(value)
}

PLdist_Clone_Time(40, 20)


pl_df <- data.frame("d00" = PLdist_Clone_Time(0, NUM_CLONES),
                    "d07" = PLdist_Clone_Time(7, NUM_CLONES),
                    "d14" = PLdist_Clone_Time(14, NUM_CLONES),
                    "d21" = PLdist_Clone_Time(21, NUM_CLONES),
                    "d28" = PLdist_Clone_Time(28, NUM_CLONES),
                    "d35" = PLdist_Clone_Time(35, NUM_CLONES),
                    "d56" = PLdist_Clone_Time(56, NUM_CLONES), 
                    "Clone.ID" = seq(1, NUM_CLONES)) %>%
  gather(-Clone.ID,  key='Time.Days', value='Clone.Size')

cc <- scales::seq_gradient_pal("#2eaafa", "#af2f98", "Lab")(seq(0, 1, length.out=7))
ggplot(pl_df)+
  geom_line(aes(x=(Clone.ID), y=log(Clone.Size), col = Time.Days), linewidth=1.0) +
  #scale_color_viridis_d(option = "D") + 
  scale_color_manual(values = cc, name = "Time since\nimmunization")+ ylim(3, 11)+
  labs(x = "log(Clone Rank)", y = "log(Clone Size)") + scale_x_log10(breaks=c(1, 3, 10, 30, 100))+
  theme(legend.position = c(0.9, 0.75), legend.background = element_blank(), legend.title.align=0.5)+ theme_bw()

#ggsave(paste0(rUN_seed, "_parent_Hist.pdf"), last_plot(), width = 6, height = 4.5, device = 'pdf')


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

ggplot(Parent_clonefreq)+
  geom_line(aes(x=Timeseries, y=clone_freq, col=Clone_num)) +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  labs(x = "Time since immunization (days)", y = "Clone Frequnecy") +
  guides(col='none') + scale_y_log10(limits=c(0.0001, 1)) #+ scale_x_log10(limits=c(3, 60))

#ggsave(paste0(rUN_seed, "_parent_freq_dist.pdf"), last_plot(), width = 6, height = 4.5, device = 'pdf')

ggplot(pwvec)+
  geom_line(aes(x=Timeseries, y=CloneSize, col=as.factor(CloneID))) +
  #scale_y_log10() + scale_x_log10() +
  labs(x = "Time since immunization (days)", y = "Clone size") +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d() + +
  #facet_wrap(.~ CloneID, nrow = 5)
  guides(col='none')

#ggsave(paste0(rUN_seed, "_parent_count_dist.pdf"), p2, width = 6, height = 4.5, device = 'pdf')

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
  ## latent variables to start the process
  TSTEP = 0.01; updated_time = 0; Tmax = T_MAX; nClones = nClones
  ## initial dist
  CARMZ_dist <- data.frame("CloneID" = as.factor(seq(1, nClones, 1)),
                             "T0" = rep(0, nClones))
  start_clone_dist <- CARMZ_dist$T0
  ## timesteps to track
  time_vec <- c(updated_time)
  
  ## looping upto Tmax
  while(updated_time < Tmax){
    current_time = round(updated_time + TSTEP, 2);
    lambda = lambda_med; persist = 1 - lambda; mu = mu_med;
    lambda_Gauss = rnorm(nClones, lambda_med, lambda_sd)
    prob_loss_Gauss = 1 - exp(-lambda_Gauss * TSTEP)
    prob_loss = 1 - exp(-lambda * TSTEP)
    prob_influx = 1 - exp(-mu * TSTEP)
    
    ## pre-existing clones
    poolsize <- ifelse(updated_time == 0, 0, length(start_clone_dist))
    
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
    
    CARMZ_dist <- CARMZ_dist %>%
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
  
  single_run_plot <- CARMZ_dist %>% replace(is.na(.), 0) %>%
    gather(-CloneID, key="Timesteps", value = "Clonefreq")  %>%
    bind_cols("Timeseries" = timeseries)%>%
    select(- Timesteps)
  
  return(single_run_plot)
}

# single_run_plot <- singleRun(NUM_CLONES)
# 
# #
# ggplot(single_run_plot)+
#  geom_line(aes(x=Timeseries, y=Clonefreq, col=as.factor(CloneID))) +
#  #scale_y_log10() + scale_x_log10() +
#  labs(x = "Time since immunization (days)", y = "Clone size") +
#  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete) +
#  #geom_text(aes(x=Timeseries, y=Clonefreq, label=CloneID, col=as.factor(CloneID)))+
#  #scale_color_viridis_d() + +
#  #facet_wrap(.~ CloneID, nrow = 5)
#  guides(col='none')
# #
# #ggsave("MZ_count_dist.pdf", last_plot(), width = 6, height = 4.5, device = 'pdf')
# #
# ggplot(single_run_plot)+
#  geom_line(aes(x=Timeseries, y=Clonefreq, col=CloneID)) +
#  labs(x = "Time since immunization (days)", y = "Clone size") +
#  guides(col='none') +
#  facet_wrap(.~ CloneID, nrow = 5)

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
  # para_run <- data.frame(mclapply(iter_seq, MultiRun_Func, nClones=nClones,
  #                                 mc.cores = detectCores()-2))
  para_run <- data.frame(mclapply(iter_seq, MultiRun_Func, nClones=nClones,
                                  mc.cores = 25))
  #para_run <- data.frame(lapply(iter_seq, MultiRun_Func, NUM_Clones=NUM_Clones, distrib=distrib))
  MultiRun_df <- InitRun %>%
    bind_cols(para_run) 
  
  return(MultiRun_df)
}

print("START BATCH RUN!!")

system.time(
  MultiRunPlot <- BatchRun(nIter = NITER, nClones = NUM_CLONES)
)
write.csv(MultiRunPlot, file = file.path("plot_sim", paste0(rUN_seed, simnum, ".csv")))

# MultiRunPlot1 <- read.csv(file.path("plot_sim", paste0(rUN_seed, 1, ".csv")))
# MultiRunPlot2 <- read.csv(file.path("plot_sim", paste0(rUN_seed, 2, ".csv")))
# MultiRunPlot3 <- read.csv(file.path("plot_sim", paste0(rUN_seed, 3, ".csv")))
# MultiRunPlot3 <- read.csv(file.path("plot_sim", paste0(rUN_seed, 3, ".csv")))
# MultiRunPlot3 <- read.csv(file.path("plot_sim", paste0(rUN_seed, 3, ".csv")))
# 
# MultiRunPlot <- rbind(MultiRunPlot1, MultiRunPlot2, MultiRunPlot3)

# MultiRunPlot <- MultiRunPlot %>% 
#   gather(-c(CloneID, Timeseries), key='key', value = 'value') %>%
#   group_by(CloneID, Timeseries) %>% replace(is.na(.), 0) %>%
#   summarize(lb = as.integer(quantile(value, probs = 0.045)),
#             ub = as.integer(quantile(value, probs = 0.955)),
#             median = quantile(value, probs = 0.5)) 


# p4 <- ggplot(MultiRunPlot)+
#   geom_line(aes(x=Timeseries, y=median, col=CloneID)) +
#   geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill=CloneID), alpha = 0.25) +
#   labs(x = "Time since immunization (days)", y = "Clone size") +
#   #scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
#   scale_color_viridis_c()+ scale_fill_viridis_c()+
#   guides(col='none', fill = 'none') +
#   facet_wrap(.~ CloneID, nrow = 5)
# # 
# 
# ggsave(paste0(rUN_seed, "_target_dist_facets.pdf"), p4, width = 10, height = 9, device = 'pdf')
# 
# ggplot(MultiRunPlot)+
#   geom_line(aes(x=Timeseries, y=median, col= as.factor(CloneID)), alpha = 0.75) +
#   geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill= as.factor(CloneID)), alpha = 0.25) +
#   labs(x = "Time since immunization (days)", y = "Clone size") +
#   scale_color_manual(values = Wes_pallete) + scale_fill_manual(values=Wes_pallete)+
#   #scale_color_viridis_c()+ scale_fill_viridis_c()+
#   guides(col='none', fill = 'none')
# 
# ggsave(paste0(rUN_seed, "_target_dist.pdf"), p5, width = 6, height = 4.5, device = 'pdf')
# 
# 
# MultiRun_clonefreq <- MultiRunPlot %>%
#   group_by(Timeseries) %>%
#   mutate(clone_freq = median/sum(median),
#          clone_freq_LB = lb/sum(lb),
#          clone_freq_UB = ub/sum(ub))
# 
# p6 <- ggplot(MultiRun_clonefreq)+
#   geom_line(aes(x=Timeseries, y=clone_freq, col=CloneID)) +
#   geom_ribbon(aes(x = Timeseries, ymin = clone_freq_LB, ymax = clone_freq_UB, fill=CloneID), alpha = 0.25) +
#   scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
#   labs(x = "Time since immunization (days)", y = "Clone Frequency") +
#   guides(col='none', fill='none') +scale_y_log10(limits=c(0.0001, 1)) #+ scale_x_log10(limits=c(1, 60))
# 
# ggsave(paste0(rUN_seed, "_target_clonefreq.pdf"), p6, width = 6, height = 4.5, device = 'pdf')


print("DONE!")








































