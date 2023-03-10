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

## Setting all the directories for opeartions
projectDir <- getwd()
saveDir <- file.path(projectDir, 'save_csv')
## compiling multiple stan objects together that ran on different nodes
#stanfit1 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_1", ".csv")))
#stanfit2 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_2",".csv")))
#stanfit3 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_3", ".csv")))
#stanfit4 <- read_stan_csv(file.path(saveDir, paste0(modelName, "_4",".csv")))
#
#fit <- sflist2stanfit(list(stanfit1, stanfit2, stanfit3, stanfit4))
## finding the parameters used in the model 
## using the last parameter("sigma4") in the array to get the total number of parameters set in the model
#num_pars <- which(fit@model_pars %in% "sigma1") -2      # the variable "sigma4" will change depdending on the data used
#parametersToPlot <- fit@model_pars[1:num_pars] 
#
### extract posterior distribution  of parameters as a matrix
#fit_ss <- as.data.frame(fit, pars= parametersToPlot) # matrix of posterior samples
#write.csv(fit_ss, file = paste0("PostDF_", modelName, ".csv"), row.names = F)

fit_ss <- read.csv(file = paste0("PostDF_", modelName, ".csv"))
mu_med <- median(fit_ss$mu)
lambda_med <- median(fit_ss$lambda_WT)


### Immune response dynamics of B cells -- reporter induction upon b cell activation
#B_cell_data <- read_excel("datafiles/New_dataCAR.xlsx", sheet =2) %>%
#  mutate(mouse.id = paste0('m_', seq(1, 72))) %>%
#  filter(mouse.id != "m_1")
#
#
## generating data for fitting
#imm_data <- B_cell_data %>% filter(genotype == "CAR") %>%
#  select(days_post_imm, contains("CARpos"))  %>%
#  filter(days_post_imm > 0) %>%
#  select(days_post_imm, CARpos_FoB, CARpos_MZB) 
#
### Immune response dynamics of B cells -- reporter induction upon b cell activation
#OLD_data <- read_excel("datafiles/TD_response_OD.xlsx", sheet =1) %>% 
#  filter(genotype == "CAR") %>%
#  mutate(tot_B_count = total_splenic_counts_millions * 1e6 * (total_B_cells/100),
#         tot_CARB_count = tot_B_count * (fraction_CAR_in_Bcells/100),
#         CARFoB = (total_FOBs/100) * tot_B_count * (fraction_CAR_in_FOBs/100),
#         CARMZB = (total_MZBs/100) * tot_B_count * (fraction_CAR_in_MZBs/100)) %>%
#  select(days.post.imm, CARFoB, CARMZB) %>%
#  gather(-days.post.imm, key = 'popln', value = 'cell_counts')
#
#FoB_baseline <- OLD_data %>% filter(popln == "CARFoB") %>% filter(days.post.imm == 0) %>% summarise(mean(cell_counts))
#MZB_baseline <- OLD_data %>% filter(popln == "CARMZB") %>% filter(days.post.imm == 0) %>% summarise(mean(cell_counts))
#
#norm_data <- imm_data %>%
#  mutate(norm_CARFoB = CARpos_FoB - FoB_baseline$`mean(cell_counts)`,
#         norm_CARMZB = CARpos_MZB - MZB_baseline$`mean(cell_counts)`) %>%
#  select(days_post_imm, contains('norm')) %>%
#  gather(-days_post_imm, key = 'popln', value = 'cell_counts')
#
#
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
#}
#
time_vec <- seq(4, 30, 0.5)
carFoB_vec <- sapply(time_vec, CAR_FoB)
#carMZB_vec <- sapply(time_vec, CAR_MZB)
ggplot() +
  geom_line(aes(x=time_vec, y=(carFoB_vec)), col=2, size=0.8) +
#  geom_line(aes(x=time_vec, y=(carMZB_vec)), col=4, size=0.8) +
#  geom_point(data=norm_data,
#             aes(x=days_post_imm, y=cell_counts, col=popln)) + 
  scale_color_manual(values = c(2, 4))+
  scale_x_continuous(limits=c(4, 30))+
  scale_y_log10() +
  labs(x='Days post immunization', y="Cell counts", title='CAR positive FoB Cells') 

ggsave("CARFoB_count.pdf", last_plot(), width = 6, height = 4.5, device = 'pdf')


## Assuming that antigen-specific B cell clones are distributed uniformly in the activated pool
## We are simulating for 30 different clones.
## The pool size of activated FoB cells varies with time. 
uniform_dist <- function(Time, nClones){ 
  rdunif(as.integer(CAR_FoB(Time)), a = 1, b = nClones)
}

#ggplot()+
#  geom_histogram(aes(x=uniform_dist(10, 30)), binwidth = 1, fill=4)

## Assuming power-law distribution for the antigen-specific B cell clones in the activated pool
## We are simulating for 30 different clones.
## The pool size of activated FoB cells varies with time. 
power_law_prob <- function(x, k, A, m, n){
  value = (A * x^-k)/((n^(-k+1) - m^(-k+1))/(-k+1))
  return(value)
}

power_law_prob(seq(1, 30), k=0.3, A=1, m=1, n=30)

Dist_power_law <- function(Time, nClones){
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

## The pool size of activated FoB cells varies with time. 
Norm_power_law <- function(Time, nClones){
  k0= 0.1; a=1; b=nClones;
  alpha = 0.03
  x_vec <- seq(1, nClones, 1)
  if (Time >= 4) {
    k = k0 + alpha * (Time -4);
    value = as.integer(CAR_FoB(Time) * (x_vec^-k)/((b^(-k+1) - a^(-k+1))/(-k+1)))
  } else {
    value = 0;
  }
  
  return(value)
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


Hist_plot_func <- function(Time){
  set.seed(1345)
  Time_week = Time * 7
  k_pow = 0.1 + 0.03 * Time_week 
  #app1 <- data.frame(f0 = sample(power_law_dist(Time, 75), as.integer(CAR_FoB(Time))),
  #                   keycol = "app1")
  app <- data.frame(f0 = sample(seq(1, 75), as.integer(CAR_FoB(Time_week)), 
                                 prob = power_law_prob(seq(1, 75), k_pow, A=1, 1, 75), replace = T),
                     keycol = paste0("Week_", Time))
  return(app)
}
time_obs <- c(1, 2, 3, 4, 8, 25)
dist_global <- data.frame()
for (i in 1:length(time_obs)){
  dist_global <- dist_global %>% rbind(Hist_plot_func(time_obs[i]))
}

ggplot(dist_global, aes(f0)) +
  geom_histogram(fill=4, alpha=0.7, binwidth=1) +
  labs(x="Clone rank", y="Count") +
  facet_wrap(~keycol, )

#NUM_Clones = 75
#pl_dist <- data.frame("Clone_Num" = (seq(1, NUM_Clones)))
#
#for (i in 1:length(time_obs)) {
#  pl_dist[1:NUM_Clones, paste0("Time_", time_obs[i])] <- Norm_power_law(time_obs[i], NUM_Clones)
#}
#
#pl_distPlot <- pl_dist %>%
#  gather(-Clone_Num, key = 'TimeObs', value = 'CloneDist') %>%
#  group_by(TimeObs) %>%
#  mutate(clone_freq = CloneDist/sum(CloneDist),
#         Time_obs = ifelse(TimeObs == "Time_4", 4, 
#                           ifelse(TimeObs == "Time_7", 7, 
#                                  ifelse(TimeObs == "Time_10", 10, 
#                                         ifelse(TimeObs == "Time_14", 14, 
#                                                ifelse(TimeObs == "Time_21", 21, 28))))))
#
#
#ggplot(pl_distPlot)+
#  geom_line(aes(x=Time_obs, y=clone_freq, col=as.factor(Clone_Num)))+ guides(col='none') +
#  ylim(0, 0.3)

cloneid <-c()  ## start an empty vector
j=0
for (i in 1:75){
  ind_a = j + 1
  j = j + length(time_vec)
  ## assortment\
  cloneid[ind_a:j] <- rep(i, length(time_vec)) 
}
pwvec <- data.frame(t(sapply(time_vec, Norm_power_law, nClones=75))) %>%
  bind_cols("Timeseries" = time_vec) %>%
  gather(-Timeseries, key='Clone_num', value="CloneSize") %>%
  bind_cols("CloneID" = cloneid)

#ggplot(pwvec)+
#  geom_point(aes(x=Timeseries, y=CloneSize, col=(CloneID))) +
#  scale_y_log10() + scale_x_log10() +
#  labs(x = "Time since immunization (days)", y = "Clone size") +
#  scale_color_viridis_c() 
#


Wes_pallete <- wesanderson::wes_palette('Darjeeling1', 75, type = "continuous")

ggplot(pwvec)+
  geom_line(aes(x=Timeseries, y=CloneSize, col=Clone_num)) +
  scale_y_log10() + scale_x_log10() +
  labs(x = "Time since immunization (days)", y = "Clone size") +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d() + +
  #facet_wrap(.~ CloneID, nrow = 5)
  guides(col='none')

ggsave("parent_count_dist.pdf", last_plot(), width = 6, height = 4.5, device = 'pdf')
  
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

time_obs <- c(4, 7, 14, 21, 28, 35)
dist_plots <- list()
for (i in 1:length(time_obs)) {
  dist_plots[[i]] <- power_law_dist(time_obs[i], 75)
}

#MP <- wesanderson::wes_palettes$Rushmore1
colfunc <- colorRampPalette(c("lightblue", "navy"))
MP<- colfunc(6)

P1 <- ggplot()+
  geom_histogram(aes(x=dist_plots[[1]]), binwidth = 1, fill=MP[1]) +
  labs(x="Clone rank", y="Count")
P2 <- ggplot()+
  geom_histogram(aes(x=dist_plots[[2]]), binwidth = 1, fill=MP[2]) +
  labs(x="Clone rank", y="Count")
P3 <- ggplot()+
  geom_histogram(aes(x=dist_plots[[3]]), binwidth = 1, fill=MP[3]) +
  labs(x="Clone rank", y="Count")
P4 <- ggplot()+
  geom_histogram(aes(x=dist_plots[[4]]), binwidth = 1, fill=MP[4]) +
  labs(x="Clone rank", y="Count")
P5 <- ggplot()+
  geom_histogram(aes(x=dist_plots[[5]]), binwidth = 1, fill=MP[5]) +
  labs(x="Clone rank", y="Count")
P6 <- ggplot()+
  geom_histogram(aes(x=dist_plots[[6]]), binwidth = 1, fill=MP[6]) +
  labs(x="Clone rank", y="Count")

cowplot::plot_grid(P1, P2, P3, P4, P5, P6)

ggsave("parent_Hist.pdf", last_plot(), width = 10, height = 6, device = 'pdf')

## Cells that differentiate into MZ and GC phenotype within a short interval time dt
### from the activated clones sample cells with propensity mu * X for MZ and alpha(Time) * X for GC B cells,
### Where X is the pool size of activated FoB cells 
MZ_time <- function(Time, distrib="powerlaw", nClones){
  mu = mu_med;
  
  clone_sample <-  if(distrib == "unif"){
    sample(uniform_dist(Time, nClones), mu * CAR_FoB(Time))
  } else {
    sample(power_law_dist(Time, nClones), mu * CAR_FoB(Time))
  }
  return(sort(clone_sample))
}


singleRun <- function(NUM_Clones, distrib){
  ### !!!TSTEP needs to be very small otherwise the propensities are not accurate!!!
  TSTEP = 0.004; updated_time = 4; Tmax = 30; current_time = 0; NUM_Clones = NUM_Clones
  start_clone_dist <- MZ_time(Time = updated_time, distrib = distrib, nClones = NUM_Clones)
  
  initial_dist <- data.frame("CloneID" = as.factor(seq(1, NUM_Clones, 1)),
                             "T0" = rep(0, NUM_Clones))
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
    parent_dist <- MZ_time(Time = current_time, distrib = distrib, nClones = NUM_Clones)
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
    ind_a = 1 + (i-1) * NUM_Clones
    ind_b = NUM_Clones*i
    timeseries[ind_a:ind_b] <- rep(time_vec[i], NUM_Clones)
  }
  
  single_run_plot <- out_dist %>%
    select(- T0) %>% replace(is.na(.), 0) %>%
    gather(-CloneID, key="Timesteps", value = "Clonefreq")  %>%
    bind_cols("Timeseries" = timeseries)%>%
    select(- Timesteps)
  
  return(single_run_plot)
}

#single_run_plot <- singleRun(75, 'powerlaw')
#
#ggplot(single_run_plot)+
#  geom_line(aes(x=Timeseries, y=Clonefreq, col=CloneID)) +
#  guides(col='none') +  
#  facet_wrap(.~ CloneID, nrow = 5)
#
#
#singleRun_clonefreq <- single_run_plot %>%
#  group_by(Timeseries) %>%
#  mutate(clone_freq = Clonefreq/sum(Clonefreq))
#
#ggplot(singleRun_clonefreq)+
#  geom_line(aes(x=Timeseries, y=clone_freq, col=CloneID)) +
#  guides(col='none') +
#  ylim(0, 0.3)

BatchRun <- function(nIter, NUM_Clones, distrib){
  ## function to iterate
  MultiRun_Func <- function(nIter, NUM_Clones, distrib){
    new_run <- singleRun(NUM_Clones=NUM_Clones, distrib=distrib)
    res1 <- data.frame(new_run$Clonefreq)
    colnames(res1) <- paste0('Clonefreq_', nIter)
    return(res1)
  }
  
  InitRun <- singleRun(NUM_Clones=NUM_Clones, distrib=distrib) 
  
  ## iterate MultiRun_Func for nIter times using lapply 
  iter_seq <- seq(1, nIter)
  para_run <- data.frame(mclapply(iter_seq, MultiRun_Func, NUM_Clones=NUM_Clones, distrib=distrib,
                                  mc.cores = detectCores()-2))
  #para_run <- data.frame(lapply(iter_seq, MultiRun_Func, NUM_Clones=NUM_Clones, distrib=distrib))
  MultiRun_df <- InitRun %>%
    bind_cols(para_run) 
  return(MultiRun_df)
}

print("START BATCH RUN!!")

system.time(
MultiRunPlot <- BatchRun(100, NUM_Clones = 75, distrib =  "powerlaw")
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


#clonevec <- c()
#for (i in 1:30) {
#  clonevec[i] <- paste0("clone_", i)
#}
#clonerep <- rep(clonevec, 62)
#  
#area_df <- new_run %>%
#  select(-CloneID, -Timeseries)
#
#area_dff <- t(area_df)  
#
#colnames(area_dff) <- clonevec
#
#area_plot <- area_dff %>%
#  bind_cols("Timeseries" = time_vec)
#
#
#  
#
#
#ggplot(data = area_plot, aes(x= Timeseries)) +
#  geom_area(aes(y= clone_1+ clone_2 + clone_3, clone_4+ clone_5 + clone_6, fill="clone_3" ), alpha=0.6 , linewidth=.5, colour="white") +
#  geom_area(aes(y= clone_1+ clone_2 + clone_3, clone_4+ clone_5 , fill="clone_3" ), alpha=0.6 , linewidth=.5, colour="white") +
#  geom_area(aes(y= clone_1+ clone_2 + clone_3, clone_4, fill="clone_3" ), alpha=0.6 , linewidth=.5, colour="white") +
#  geom_area(aes(y= clone_1+ clone_2 + clone_3, fill="clone_3" ), alpha=0.6 , linewidth=.5, colour="white") +
#  geom_area(aes(y= clone_1 + clone_2, fill="clone_2"), alpha=0.6 , linewidth=.5, colour="white") +
#  geom_area(aes(y= clone_1, fill="clone_1"), alpha=0.6 , linewidth=.5, colour="white") +
#  #scale_fill_viridis(discrete = T) +
#  #theme_ipsum() + 
#  scale_x_log10(limits = c(4, 30), breaks=c(7, 14, 35)) + scale_y_log10() + 
#  labs(x='Days since immunization', y=NULL)  + 
#  guides(fill="none")
#  
#  
#  



























