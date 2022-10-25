rm(list = ls()); gc();

############################################################
### Preamble
############################################################

## loading libraries
library(tidyverse)
library(readxl)
library("ggsci")


#### plotting style
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                 axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 #legend.title = element_text(size=12,  hjust = 0.5, face = "bold"),
                 legend.background = element_blank(),
                 legend.box.background =  element_blank(),
                 legend.key = element_blank())

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks = as.numeric(1:10 %o% 10^(-4:8))


############################################################
############################################################

## import data

############################################################
## Immune response dynamics of B cells -- reporter induction upon b cell activation
B_cell_data <- read_excel("datafiles/New_dataCAR.xlsx", sheet =2) %>%
  mutate(mouse.id = paste0('m_', seq(1, 72))) %>%
  filter(mouse.id != "m_1")


B_cell_countsWT_df <- B_cell_data %>% 
  filter(genotype == "CAR", days_post_imm > 0) %>%
  select(days_post_imm, genotype, contains('CARpos')) %>%
  gather(-c(days_post_imm, genotype), key = "subpop", value="cell_counts")%>%
  mutate(#CARexpr = ifelse(grepl("CARpos", subpop), "CARpos", "CARneg"),
         subplop = ifelse(grepl("Fo", subpop), "FoB", ifelse(grepl("GC", subpop), "GCB", ifelse(grepl("MZ", subpop), "MZB", "TotalB")))) 

B_cell_countsN2_df <- B_cell_data %>% 
  filter(genotype == "N2KO", days_post_imm > 0) %>%
  select(days_post_imm, genotype, contains('CARpos')) %>%
  gather(-c(days_post_imm, genotype), key = "subpop", value="cell_counts")%>%
  mutate(#CARexpr = ifelse(grepl("CARpos", subpop), "CARpos", "CARneg"),
    subplop = ifelse(grepl("Fo", subpop), "FoB", ifelse(grepl("GC", subpop), "GCB", ifelse(grepl("MZ", subpop), "MZB", "TotalB")))) 


blank_data <- data.frame(subplop = c('FoB', 'FoB', 'GCB', 'GCB', 'MZB', 'MZB','TotalB', 'TotalB'),
                         cell_counts = c(5e4, 3e6, 1e4, 3e6, 1e4, 5e5, 1e5, 5e6),
                         days_post_imm = rep(c(4, 30), 4),
                         CARexpr = rep("CARpos", 8))

Bcell_all_data <- rbind(B_cell_countsWT_df, B_cell_countsN2_df) %>%
  mutate(mouse_strain = ifelse(grepl("CAR", genotype), "WT", "N2KO"))

ggplot() +
  geom_point(data=Bcell_all_data, aes(x= days_post_imm, y=cell_counts, col= mouse_strain), size=1.7) + 
  scale_color_manual(values = c("#CA0636", "#0B78BA"), name=NULL) +
  #geom_point(data=CAR_cell_counts_df, aes(x= (days_post_imm), y=(cell_counts)), col="navy",size=2, alpha=0.7) +
  geom_blank(data=blank_data, aes(x= (days_post_imm), y=(cell_counts))) +
  scale_y_log10(labels=fancy_scientific) + 
  labs(x='Days post immunization', y='Cell counts')+
  facet_wrap(. ~ subplop, scales = 'free') + myTheme 

ggsave('plots/New_plotsNatComm/dataplot1.pdf', last_plot(), device = 'pdf', width = 12, height = 8)

#### plotting correlations
CAR_cell_counts_df <- B_cell_data %>% 
  filter(genotype == "CAR", days_post_imm >=0) %>% 
  mutate(days_bin = ifelse(days_post_imm < 7, "1Wk",
                           ifelse(days_post_imm <14, "2Wk",
                                  ifelse(days_post_imm <21, "3Wk", "4Wk"))))


corFOMZ <- cor.test(B_cell_data$CARpos_FoB, B_cell_data$CARpos_MZB)
corGCMZ <- cor.test(B_cell_data$CARpos_GCB, B_cell_data$CARpos_MZB)
corFOGC <- cor.test(B_cell_data$CARpos_FoB, B_cell_data$CARpos_GCB)

p11 <-ggplot(data = CAR_cell_counts_df)+
  geom_point(aes(x=log(CARpos_FoB), y=log(CARpos_MZB), col=days_bin), size=2)+
  #scale_color_manual(values = wesanderson::wes_palette("Zissou1", n=4, type = 'continuous'))+
  scale_color_npg()+ guides(col='none') +
  labs(x="Log(Counts of CAR positive FoB)", y="Log(Counts of CAR positive MZB)")+
  geom_text(x=12.3, y=12.4, label=paste0("r=", round(corFOMZ$estimate, 2)), col='darkred', size=5) +myTheme


p22 <-ggplot(data = CAR_cell_counts_df) +
  geom_point(aes(x=log(CARpos_GCB), y=log(CARpos_MZB), col=days_bin), size=2)+
  #scale_color_manual(values = wesanderson::wes_palette("Zissou1", n=4, type = 'continuous'))+
  scale_color_npg()+ guides(col='none') +
  labs(x="Log(Counts of CAR positive GCB)", y="Log(Counts of CAR positive MZB)")+
  geom_text(x=11, y=12.4, label=paste0("r=", round(corGCMZ$estimate, 2)), col='darkred', size=5)+ myTheme


p33 <- ggplot(data = CAR_cell_counts_df)+
  geom_point(aes(x=log(CARpos_FoB), y=log(CARpos_GCB), col=days_bin), size=2) +
  labs(x="Log(Counts of CAR positive FoB)", y="Log(Counts of CAR positive GCB)") +
  geom_text(x=12.3, y=14.8, label=paste0("r=", round(corFOGC$estimate, 2)), col='darkred',  size=5) +
  scale_color_npg(name = "Days post\nimmunization") + myTheme + theme(legend.position = c(0.85, 0.2))


cowplot::plot_grid( p11, p22, p33, nrow  = 1)
ggsave('plots/natcomm/dataplot2.pdf', last_plot(), device = 'pdf', width = 12, height = 4)



###############################
##############################

# Modeling Precursor Dynamics 
phi_func <- function(Time, basl, theta, n, X, q){
  t = Time - 0
  theta_exp = exp(theta)
  basl_exp = exp(basl)
  return( basl_exp + (theta_exp * t^n) * (1 - ((t^q)/((X^q) + (t^q)))))
}

sol_time <- seq(0, 30, length.out=100)
qplot(sol_time, y= log(phi_func(sol_time, basl = 17.2, theta = 7.5, n = 3 , X = 10, q = 5)))+
  geom_line()+
  geom_point(data = B_cell_data, aes(x=days_post_imm, y=log(total_FoB)), col=2)


repFOB_nlm <- nls(log(total_FoB) ~ log(phi_func(days_post_imm, basl, theta, n, X, q)),
                  data = B_cell_data,
                  start = list(basl = 11.2, theta = 7.5, n = 3 , X = 10, q = 5))

par_est <- coef(repFOB_nlm)
phi_vec_m1 <- phi_func(sol_time, basl = par_est[1], theta =  par_est[2], n =  par_est[3] , X =  par_est[4], q =  par_est[5])
phi_vec_m <- phi_func(sol_time, basl = 11.2, theta = 7.5, n = 3 , X = 10, q = 5)

ggplot() +
  #geom_line(aes(x=sol_time, y=(phi_vec_m)), col="#6082B6", size=0.8) +
  geom_line(aes(x=sol_time, y=(phi_vec_m1)), col=2, size=0.8) +
  geom_point(data=filter(B_cell_data, days_post_imm >=4),
             aes(x=days_post_imm, y=(total_FoB)), size=2, col="#6082B6") + 
  scale_x_continuous(limits=c(3.9, 30))+
  scale_y_log10(limits = c(5e6, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y="Cell counts", title='CAR positive FoB Cells') + myTheme

ggsave('plots/New_plotsNatComm/CARposFOB.pdf', last_plot(), device = 'pdf', width = 6, height = 4.5)

### CAR neg MZ B cells as precursors

phi2_func <- function(Time, basl, nu, b0){
  t = Time - 0
  basl_exp = exp(basl)
  return( basl_exp * (1 + exp(-nu * (Time - b0)^2)))
}
qplot(sol_time, y= log(phi2_func(sol_time, basl = 14, nu = 0.01, b0=8)))+
  geom_line()+
  geom_point(data = B_cell_data, aes(x=days_post_imm, y=log(CARneg_MZB)), col=2)


Neg_MZB_nlm <- nls(log(CARneg_MZB) ~ log(phi2_func(days_post_imm,  basl, nu, b0)),
                  data = B_cell_data,
                  start = list(basl = 11, nu = 0.01, b0=8))

par_est <- coef(Neg_MZB_nlm)
sol_time <- seq(4, 30, length.out=100)
phi_vec_m <- phi2_func(sol_time, basl=par_est[1], nu=par_est[2], b0=par_est[3])
phi_vec_mm <- phi2_func(sol_time, basl=14.06, nu=0.0033, b0=20.58)

ggplot() +
  geom_line(aes(x=sol_time, y=(phi_vec_m)), col="blue", size=2) +
  geom_line(aes(x=sol_time, y=(phi_vec_mm)), col="#CD7F32", size=0.8) +
  geom_point(data=filter(B_cell_data, days_post_imm >=4),
             aes(x=days_post_imm, y=(CARneg_MZB)), size=2, col="#CD7F32") +
  scale_x_continuous(limits=c(3.9, 30))+
  scale_y_log10(limits = c(1e5, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y="Cell counts", title='CAR negative MZ B Cells') + myTheme


### tot FOb

phi2_func <- function(Time, basl, nu, b0){
  t = Time - 0
  basl_exp = exp(basl)
  return( basl_exp * (1 + exp(-nu * (Time - b0)^2)))
}

totFOB_nlm <- nls(log(total_FoB) ~ log(phi2_func(days_post_imm,  basl, nu, b0)),
                   data = B_cell_data,
                   start = list(basl = 11, nu = 0.01, b0=8))

par_est <- coef(totFOB_nlm)
sol_time <- seq(4, 30, length.out=100)
phi_vec_m <- phi2_func(sol_time, basl=par_est[1], nu=par_est[2], b0=par_est[3])
phi_vec_mm <- phi2_func(sol_time, basl=16.7, nu=0.004, b0=20)


ggplot() +
  #geom_line(aes(x=sol_time, y=(phi_vec_m)), col="blue", size=2) +
  geom_line(aes(x=sol_time, y=(phi_vec_mm)), col="#CD7F32", size=0.8) +
  geom_hline(yintercept = exp(17.1), col=2, size=0.8) +
  geom_point(data=filter(B_cell_data, days_post_imm >=4),
             aes(x=days_post_imm, y=(total_FoB)), size=2, col="#CD7F32") +
  scale_x_continuous(limits=c(3.9, 30))+
  scale_y_log10(limits = c(1e6, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y="Cell counts", title='total FoB Cells') + myTheme


ggsave('plots/New_plotsNatComm/CARnegMZB.pdf', last_plot(), device = 'pdf', width = 6, height = 4.5)

###############################
##############################

# generating data for fitting
imm_data <- B_cell_data %>% filter(genotype == "CAR") %>%
  select(days_post_imm, contains("CARpos"))  %>%
  filter(days_post_imm >= 0)

write.csv(imm_data, file.path("datafiles", "Bcell_imm_data.csv"), row.names = F)

imm_data %>% 
  filter(days_post_imm == 0) %>%
  summarise("CAR_countsMZ" = log(mean(CARpos_MZB)),
            "CAR_countsGC" = log(mean(CARpos_GCB)))

# generating data for fitting
imm_N2ko_data <- B_cell_data %>% filter(genotype == "N2KO") %>%
  select(days_post_imm, contains("CARpos"))  %>%
  filter(days_post_imm > 0)

imm_N2ko_data %>% 
  filter(days_post_imm == 0) %>%
  summarise("CAR_countsMZ" = log(mean(CARpos_MZB)),
            "CAR_countsGC" = log(mean(CARpos_GCB)))


write.csv(imm_N2ko_data, file.path("datafiles", "N2KO_imm_data.csv"), row.names = F)

### all time points in data
data_times1 <- imm_data$days_post_imm 
data_times2 <- imm_N2ko_data$days_post_imm 

## Unique time points with indices to map
unique_times_df1 <- imm_data %>% distinct(days_post_imm, .keep_all = TRUE)
solve_times <- unique_times_df1$days_post_imm
## Map of the unique time points on all the timepoints
time_index1 <- purrr::map_dbl(data_times1, function(x) which(x == solve_times))    # keeping track of index of time point in relation to solve_time
time_index2 <- purrr::map_dbl(data_times2, function(x) which(x == solve_times))    # keeping track of index of time point in relation to solve_time

## Data to import in Stan
numObs1 <- length(imm_data$days_post_imm)
numObs2 <- length(imm_N2ko_data$days_post_imm)
n_shards <- length(solve_times)
solve_time <- solve_times
time_index1 <- time_index1
time_index2 <- time_index2
CAR_MZ_counts <- imm_data$CARpos_MZB
CAR_MZN2_counts <- imm_N2ko_data$CARpos_MZB
CAR_GC_counts <- imm_data$CARpos_GCB
CAR_GCN2_counts <- imm_N2ko_data$CARpos_GCB

# time sequence for predictions specific to age bins within the data
ts_pred <- seq(4, 30, length.out = 500)
numPred <- length(ts_pred)


rstan::stan_rdump(c("numObs1",  "n_shards", "solve_time", "time_index1", "numObs2", "time_index2",
             "CAR_MZ_counts",   "CAR_GC_counts", "CAR_MZN2_counts", "CAR_GCN2_counts", 
             "ts_pred", "numPred"),
           file = file.path('datafiles', paste0('Bcell_Imm',".Rdump")))

###############################
##############################

## testing models from stan
library(rstan)

stanmodel_file <- file.path("stan_models/New_modelsCAR", "Branched_timeinflux.stan")
expose_stan_functions(stanmodel_file)


ts_seq <- c(7, 14, 30)
init_cond <- c(exp(11.5), exp(10.8), exp(11.5))
params <- c(0.05, 0.02, 0.01, 0.04, 0.2, 0.1, 0.01)
ode_sol <- solve_ODE_sys(ts_seq, init_cond, params)

solve_ODE_init(exp(9.9), params)

ts_pred <- seq(4.1, 30, length.out=100)
ode_df <- solve_ODE_sys(ts_pred, init_cond, params)
stan_pred_df <- data.frame("time_pred" = ts_pred,
                           "y_pred" = matrix(unlist(ode_df), nrow = length(ode_df), byrow = TRUE))%>%
  mutate(CAR_GC = y_pred.1,
         CAR_MZ = y_pred.2,
         CAR_GCN2 = y_pred.3) %>%
  select(time_pred, contains("CAR")) %>%
  gather(-time_pred, key="subplop", value="counts")


ggplot(stan_pred_df)+
  geom_point(aes(x=time_pred, y=counts, col=subplop))+
  scale_y_log10()+
  facet_grid(.~subplop)












