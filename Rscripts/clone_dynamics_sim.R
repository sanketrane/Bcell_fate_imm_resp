## Simulations of clonal dynamics using the model fits to immunization data

## loading libraries
library(tidyverse)
library(wesanderson)
library(parallel)

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
  k0= 0.001; a=1; b=nClones;
  alpha = 0.02
  k = k0 + alpha * Time
  
  x_vec <- seq(1, nClones, 1)
 as.integer(CAR_FoB(Time) * (x_vec^-k)/((b^(-k+1) - a^(-k+1))/(-k+1)))
}


time_obs <- c(4, 7, 10, 14, 21, 28)
NUM_Clones = 10
pl_dist <- data.frame("Clone_Num" = (seq(1, NUM_Clones)))

for (i in 1:length(time_obs)) {
  pl_dist[1:NUM_Clones, paste0("Time_", time_obs[i])] <- Norm_power_law(time_obs[i], NUM_Clones)
}

pl_distPlot <- pl_dist %>%
  gather(-Clone_Num, key = 'TimeObs', value = 'CloneDist') %>%
  group_by(TimeObs) %>%
  mutate(clone_freq = CloneDist/sum(CloneDist),
         Time_obs = ifelse(TimeObs == "Time_4", 4, 
                           ifelse(TimeObs == "Time_7", 7, 
                                  ifelse(TimeObs == "Time_10", 10, 
                                         ifelse(TimeObs == "Time_14", 14, 
                                                ifelse(TimeObs == "Time_21", 21, 28))))))


ggplot(pl_distPlot)+
  geom_line(aes(x=Time_obs, y=clone_freq, col=as.factor(Clone_Num)))+ guides(col='none') +
  ylim(0, 0.3)


time_vec <- seq(4, 30, 0.5)
pwvec <- data.frame(t(sapply(time_vec, Norm_power_law, nClones=10))) %>%
  bind_cols("Timeseries" = time_vec) %>%
  gather(-Timeseries, key='Clone_num', value="CloneSize")

ggplot(pwvec)+
  geom_line(aes(x=Timeseries, y=CloneSize, col=Clone_num)) +
  scale_y_log10() + scale_x_log10() + guides(col='none')
  


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

pl_histPlot <- power_law_dist(10, 30)

ggplot()+
  geom_histogram(aes(x=pl_histPlot), binwidth = 1, fill=2) +
  labs(x="Clone rank", y="Count")



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



#singleRun <- function(nClones, distrib="powerlaw"){
#  time_vec <- c(4, 7, 10, 14, 18, 21, 24, 27, 30)
#  
#  MZ_bin <- data.frame("CloneID" = as.factor(seq(1, nClones, 1)),
#                       "T0" = rep(0, nClones))
#  
#  for (i in 1:length(time_vec)){
#    new_freq_df <- data.frame(table(MZ_time(Time = time_vec[i], distrib = distrib, nClones = nClones))) %>%
#      rename(CloneID=Var1,
#             !!paste0('TSTEP_', i) := Freq)
#    MZ_bin <- MZ_bin %>%
#      left_join(new_freq_df, by="CloneID")
#  }           
#  
#  MZ_bin_plot <- MZ_bin %>%
#    select(-T0) %>%
#    gather(-CloneID, key="Timesteps", value = "Clonefreq") %>%
#    mutate(Timeseries = ifelse(Timesteps == "TSTEP_1", 4,
#                               ifelse(Timesteps == "TSTEP_2", 7,
#                                      ifelse(Timesteps == "TSTEP_3", 14,
#                                             ifelse(Timesteps == "TSTEP_4", 21,
#                                                    ifelse(Timesteps == "TSTEP_5", 18,
#                                                           ifelse(Timesteps == "TSTEP_6", 24,
#                                                                  ifelse(Timesteps == "TSTEP_7", 27, 30)))))))) %>%
#    select(-Timesteps)
#}
#
#single_run_plot <- singleRun(30)
#
#ggplot(single_run_plot)+
#  geom_point(aes(x=Timeseries, y=Clonefreq, col=CloneID)) +
#  guides(col='none') +
#  facet_wrap(.~ CloneID, nrow = 5)
#
#
#BatchRun <- function(nIter, distrib="powerlaw", nClones){
#  multiRun <- singleRun(nClones) 
#  for (i in 1:nIter){
#    new_run <- singleRun(nClones = nClones, distrib = distrib) 
#    multiRun <- multiRun %>%
#      mutate(!!paste0("Clonefreq_", i) := new_run$Clonefreq)
#  }
#  multiRun_plot <- multiRun %>%
#    gather(-c(CloneID, Timeseries), key='key', value = 'value') %>%
#    group_by(CloneID, Timeseries) %>% replace(is.na(.), 0) %>%
#    summarize(lb = quantile(value, probs = 0.045),
#              median = quantile(value, probs = 0.5),
#              ub = quantile(value, probs = 0.955)) 
#}
#
#  
#multiRun_plot <- BatchRun(30, distrib="unif", 30)
#
#Wes_pallete <- wesanderson::wes_palette('Darjeeling1', 30, type = "continuous")
#
#ggplot(multiRun_plot)+
#  geom_line(aes(x=Timeseries, y=median, col=CloneID)) +
#  geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill=CloneID), alpha = 0.25) +
#  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
#  #scale_color_viridis_d()+ scale_fill_viridis_d()+
#  guides(col='none', fill = 'none') +
#  facet_wrap(.~ CloneID, nrow = 5)

#influx_func <- function(s){
#  mu = 0.0025;
#  CAR_FoB(s) * mu * exp(-mu * s)
#}
#
#integrate(influx_func, lower = 10, upper = 10.04)$value


singleRun <- function(NUM_Clones, distrib){
  ### !!!TSTEP needs to be very small otherwise the propensities are not accurate!!!
  TSTEP = 0.04; updated_time = 4; Tmax = 30; current_time = 0; NUM_Clones = NUM_Clones
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
    lambda = 0.81; persist = 1 - lambda; mu = 0.0025;
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

single_run_plot <- singleRun(10, 'powerlaw')

ggplot(single_run_plot)+
  geom_line(aes(x=Timeseries, y=Clonefreq, col=CloneID)) +
  guides(col='none') +  
  facet_wrap(.~ CloneID, nrow = 5)


singleRun_clonefreq <- single_run_plot %>%
  group_by(Timeseries) %>%
  mutate(clone_freq = Clonefreq/sum(Clonefreq))

ggplot(singleRun_clonefreq)+
  geom_line(aes(x=Timeseries, y=clone_freq, col=CloneID)) +
  guides(col='none') +
  ylim(0, 0.3)

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
    bind_cols(para_run) #%>%
    #gather(-c(CloneID, Timeseries), key='key', value = 'value') %>%
    #group_by(CloneID, Timeseries) %>% replace(is.na(.), 0) %>%
    #summarize(lb = quantile(value, probs = 0.045),
    #          median = quantile(value, probs = 0.5),
    #          ub = quantile(value, probs = 0.955)) 
  
  return(MultiRun_df)
}

system.time(
MultiRunPlot <- BatchRun(10, NUM_Clones = 30, distrib =  "powerlaw")
)

Wes_pallete <- wesanderson::wes_palette('Darjeeling1', 30, type = "continuous")

ggplot(MultiRunPlot)+
  geom_line(aes(x=Timeseries, y=median, col=CloneID)) +
  geom_ribbon(aes(x = Timeseries, ymin = lb, ymax = ub, fill=CloneID), alpha = 0.25) +
  scale_color_manual(values = Wes_pallete)+ scale_fill_manual(values=Wes_pallete)+
  #scale_color_viridis_d()+ scale_fill_viridis_d()+
  guides(col='none', fill = 'none') +
  facet_wrap(.~ CloneID, nrow = 5)


clonevec <- c()
for (i in 1:30) {
  clonevec[i] <- paste0("clone_", i)
}
clonerep <- rep(clonevec, 62)
  
area_df <- new_run %>%
  select(-CloneID, -Timeseries)

area_dff <- t(area_df)  

colnames(area_dff) <- clonevec

area_plot <- area_dff %>%
  bind_cols("Timeseries" = time_vec)


  


ggplot(data = area_plot, aes(x= Timeseries)) +
  geom_area(aes(y= clone_1+ clone_2 + clone_3, clone_4+ clone_5 + clone_6, fill="clone_3" ), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= clone_1+ clone_2 + clone_3, clone_4+ clone_5 , fill="clone_3" ), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= clone_1+ clone_2 + clone_3, clone_4, fill="clone_3" ), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= clone_1+ clone_2 + clone_3, fill="clone_3" ), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= clone_1 + clone_2, fill="clone_2"), alpha=0.6 , linewidth=.5, colour="white") +
  geom_area(aes(y= clone_1, fill="clone_1"), alpha=0.6 , linewidth=.5, colour="white") +
  #scale_fill_viridis(discrete = T) +
  #theme_ipsum() + 
  scale_x_log10(limits = c(4, 30), breaks=c(7, 14, 35)) + scale_y_log10() + 
  labs(x='Days since immunization', y=NULL)  + 
  guides(fill="none")
  
  
  



























