rm(list = ls()); gc();

############################################################
### Preamble
############################################################

## loading libraries
library(tidyverse)
library(readxl)

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
B_cell_data <- read_excel("datafiles/CAR_corrected_data.xlsx", sheet =1) %>%
  mutate(#### Deriving the counts of each B cell subset
    total_spleen_counts = total_splenic_counts_millions * 1e6,
    B_cell_counts = (total_B_cells/100) * total_spleen_counts,
    CARpos_counts = (CARpos_Bcells/100) * B_cell_counts,
    CARpos_nonGC_numbers = (nonGC_in_CARpos/100) * CARpos_counts,
    CARpos_GC_numbers = (GC_in_CARpos/100) * CARpos_counts,
    CARpos_Fo_numbers = (FoB_in_CARpos/100) * CARpos_nonGC_numbers,
    CARpos_MZ_numbers = (MZB_in_CARpos/100) * CARpos_nonGC_numbers,
    CARneg_counts = (CARneg_Bcells/100) * B_cell_counts,
    CARneg_nonGC_numbers = (nonGC_in_CARneg/100) * CARneg_counts,
    CARneg_GC_numbers = (GC_in_CARneg/100) * CARneg_counts,
    CARneg_Fo_numbers = (FoB_in_CARneg/100) * CARneg_nonGC_numbers,
    CARneg_MZ_numbers = (MZB_in_CARneg/100) * CARneg_nonGC_numbers) 


B_cell_countsWT_df <- B_cell_data %>%
  select(days_post_imm, genotype, contains('numbers'), -contains('nonGC')) %>%
  gather(-c(days_post_imm, genotype), key = "subpop", value="cell_counts")%>%
  mutate(CARexpr = ifelse(grepl("CARpos", subpop), "CARpos", "CARneg"),
         subplop = ifelse(grepl("Fo", subpop), "FoB", ifelse(grepl("GC", subpop), "GC", "MZB"))) %>% 
filter(genotype == "CAR", CARexpr == "CARpos") 

B_cell_countsN2_df <- B_cell_data %>%
  select(days_post_imm, genotype, contains('numbers'), -contains('nonGC')) %>%
  gather(-c(days_post_imm, genotype), key = "subpop", value="cell_counts")%>%
  mutate(CARexpr = ifelse(grepl("CARpos", subpop), "CARpos", "CARneg"),
         subplop = ifelse(grepl("Fo", subpop), "FoB", ifelse(grepl("GC", subpop), "GC", "MZB"))) %>% 
  filter(genotype == "N2KO", CARexpr == "CARpos") 


Old_Bcell_data <-  read_excel("datafiles/NEW_CAR_fractions_TD_response.xlsx", sheet =1) %>%
  # filter(genotype == "CAR") 
  mutate(B_cell_numbers = (total_B_cells/100) * total_splenic_counts_millions * 1e6,
         FOB_cell_numbers = (total_FOBs/100) * B_cell_numbers,
         MZB_cell_numbers = (total_MZBs/100) * B_cell_numbers,
         GCB_cell_numbers = (total_Gcs/100) * B_cell_numbers,
         PCB_cell_numbers = (total_PCs/100) * B_cell_numbers,
         CAR_B_numbers = (fraction_CAR_in_Bcells/100) * B_cell_numbers,
         CAR_FOB_numbers = (fraction_CAR_in_FOBs/100) * FOB_cell_numbers,
         CAR_MZB_numbers = (fraction_CAR_in_MZBs/100) * MZB_cell_numbers,
         CAR_GCB_numbers = (fraction_CAR_in_GCs/100)* GCB_cell_numbers,
         CAR_PCB_numbers = (fraction_CAR_in_PCs/100)* PCB_cell_numbers) 


CAR_cell_counts_df <- Old_Bcell_data %>% 
  select(days.post.imm, genotype, contains("CAR"))%>%
  select(days.post.imm, genotype, contains('numbers'), - CAR_B_numbers, -CAR_PCB_numbers) %>%
  gather(-c(days.post.imm, genotype), key = 'subpop', value='cell_counts') %>%
  mutate(subplop = ifelse(grepl("FO", subpop), "FoB", ifelse(grepl("GC", subpop), "GC", "MZB"))) %>% 
  filter(genotype == "CAR", days.post.imm>0) 


blank_data <- data.frame(subplop = c('FoB', 'FoB', 'GC', 'GC', 'MZB', 'MZB'),
                         cell_counts = c(1e5, 2e6, 1e4, 4e6, 1e4, 1e6),
                         days_post_imm = rep(c(4, 30), 3),
                         CARexpr = rep("CARpos", 6))

Bcell_all_data <- rbind(B_cell_countsWT_df, B_cell_countsN2_df) %>%
  mutate(mouse_strain = ifelse(grepl("CAR", genotype), "WT", "N2KO"))

p1 <- ggplot() +
  geom_point(data=Bcell_all_data, aes(x= days_post_imm, y=cell_counts, col= mouse_strain), size=1.7) + 
  scale_color_manual(values = c("#CA0636", "#0B78BA"), name=NULL) +
  #geom_point(data=CAR_cell_counts_df, aes(x= (days.post.imm), y=(cell_counts)), col="navy",size=2, alpha=0.7) +
  geom_blank(data=blank_data, aes(x= (days_post_imm), y=(cell_counts))) +
  scale_y_log10() + 
  labs(x='Days post immunization', y='Cell counts')+
  facet_wrap(. ~ subplop, scales = 'free') + myTheme 


#### plotting correlations
CAR_cell_counts22_df <- B_cell_data %>% 
  filter(genotype == "CAR", days_post_imm > 0) %>% 
  select(days_post_imm, contains('numbers'), -contains('nonGC')) %>%
  mutate(days_bin = ifelse(days_post_imm < 7, "1Wk",
                           ifelse(days_post_imm <14, "2Wk",
                                  ifelse(days_post_imm <21, "3Wk", "4Wk"))))

library("ggsci")


corFOMZ <- cor.test(CAR_cell_counts22_df$CARpos_Fo_numbers, CAR_cell_counts22_df$CARpos_MZ_numbers)
corGCMZ <- cor.test(CAR_cell_counts22_df$CARpos_GC_numbers, CAR_cell_counts22_df$CARpos_MZ_numbers)
corFOGC <- cor.test(CAR_cell_counts22_df$CARpos_Fo_numbers, CAR_cell_counts22_df$CARpos_GC_numbers)

p11 <-ggplot(data = CAR_cell_counts22_df)+
  geom_point(aes(x=log(CARpos_Fo_numbers), y=log(CARpos_MZ_numbers), col=days_bin), size=2)+
  #scale_color_manual(values = wesanderson::wes_palette("Zissou1", n=4, type = 'continuous'))+
  scale_color_npg()+ guides(col='none') +
  labs(x="Log(Counts of CAR positive FoB)", y="Log(Counts of CAR positive MZB)")+
  geom_text(x=12.1, y=12.4, label=paste0("r=", round(corFOMZ$estimate, 2)), size=5) +myTheme


p22 <-ggplot(data = CAR_cell_counts22_df) +
  geom_point(aes(x=log(CARpos_GC_numbers), y=log(CARpos_MZ_numbers), col=days_bin), size=2)+
  #scale_color_manual(values = wesanderson::wes_palette("Zissou1", n=4, type = 'continuous'))+
  scale_color_npg()+ guides(col='none') +
  labs(x="Log(Counts of CAR positive GCB)", y="Log(Counts of CAR positive MZB)")+
  geom_text(x=10.2, y=12.4, label=paste0("r=", round(corGCMZ$estimate, 2)), size=5)+ myTheme


p33 <- ggplot(data = CAR_cell_counts22_df)+
  geom_point(aes(x=log(CARpos_Fo_numbers), y=log(CARpos_GC_numbers), col=days_bin), size=2) +
  labs(x="Log(Counts of CAR positive FoB)", y="Log(Counts of CAR positive GCB)") +
  geom_text(x=12.1, y=14.2, label=paste0("r=", round(corFOGC$estimate, 2)), size=5) +
  scale_color_npg(name = "Days post\nimmunization") + myTheme + theme(legend.position = c(0.875, 0.2))



## saving  plots for quality control 
pdf(file = file.path("plots", paste("New_plotsNatComm/DataPlots%03d.pdf", sep = "")),
    width = 13, height = 4, onefile = FALSE, useDingbats = FALSE)
p1
cowplot::plot_grid( p11, p22, p33, nrow  = 1)
dev.off()


###############################
##############################

# Modeling Precursor Dynamics 
phi_func <- function(Time, basl, theta, n, X, q){
  t = Time - 0
  theta_exp = exp(theta)
  basl_exp = exp(basl)
  return( basl_exp + (theta_exp * t^n) * (1 - ((t^q)/((X^q) + (t^q)))))
}
qplot(sol_time, y= log(phi_func(sol_time, basl = 11.5, theta = 7, n = 3 , X = 10, q = 5)))+
  geom_line()+
  geom_point(data = CAR_cell_counts22_df, aes(x=days_post_imm, y=log(CARpos_Fo_numbers)), col=2)


repFOB_nlm <- nls(log(CARpos_Fo_numbers) ~ log(phi_func(days_post_imm, basl, theta, n, X, q)),
                  data = CAR_cell_counts22_df,
                  start = list( basl = 11.5, theta = 7, n = 3 , X = 10, q = 5))

par_est <- coef(repFOB_nlm)
sol_time <- seq(4, 30, length.out=100)
phi_vec_m <- phi_func(sol_time, basl = 11.5, theta = 7, n = 3 , X = 10, q = 5)

ggplot() +
  geom_line(aes(x=sol_time, y=(phi_vec_m)), col="#6082B6", size=0.8) +
  geom_point(data=filter(CAR_cell_counts22_df, days_post_imm >=4),
             aes(x=days_post_imm, y=(CARpos_Fo_numbers)), size=2, col="#6082B6") + 
  scale_x_continuous(limits=c(4, 30))+
  scale_y_log10(limits = c(5e4, 2e6), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y="Cell counts", title='CAR positive FoB Cells') + myTheme


### CAR neg MZ B cells as precursors

phi_func <- function(Time, basl, nu, b0){
  t = Time - 0
  basl_exp = exp(basl)
  return( basl_exp * (1 + exp(-nu * (Time - b0)^2)))
}
qplot(sol_time, y= log(phi_func(sol_time, basl = 14, nu = 0.01, b0=8)))+
  geom_line()+
  geom_point(data = B_cell_data, aes(x=days_post_imm, y=log(CARneg_MZ_numbers)), col=2)


Neg_MZB_nlm <- nls(log(CARneg_MZ_numbers) ~ log(phi_func(days_post_imm,  basl, nu, b0)),
                  data = B_cell_data,
                  start = list(basl = 11, nu = 0.01, b0=8))

par_est <- coef(Neg_MZB_nlm)
sol_time <- seq(4, 30, length.out=100)
phi_vec_m <- phi_func(sol_time, basl=par_est[1], nu=par_est[2], b0=par_est[3])
phi_vec_mm <- phi_func(sol_time, basl=14, nu=0.01, b0=18)

ggplot() +
  #geom_line(aes(x=sol_time, y=(phi_vec_m)), col=2, size=0.8) +
  geom_line(aes(x=sol_time, y=(phi_vec_mm)), col="#CD7F32", size=0.8) +
  geom_point(data=filter(CAR_PosNeg_df, days.post.imm >=4),
             aes(x=days.post.imm, y=(CAR_Neg_MZB)), size=2, col="#CD7F32") + scale_x_continuous(limits=c(0, 30))+
  scale_y_log10(limits = c(1e5, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y="Cell counts", title='CAR negative MZ B Cells') + myTheme


###############################
##############################

# generating data for fitting
imm_data <- B_cell_data %>% filter(genotype == "CAR") %>%
  select(days.post.imm, contains("MZ"), contains("GC"), -contains("fraction"), -total_MZBs) %>%
  mutate(fraction_CAR_MZ = CAR_MZB_numbers/MZB_cell_numbers,
         fraction_CAR_GC = CAR_GCB_numbers/GCB_cell_numbers)%>%
  select(days.post.imm, CAR_MZB_numbers, GCB_cell_numbers, CAR_GCB_numbers, MZB_cell_numbers) %>%
  filter(days.post.imm > 0)

write.csv(imm_data, file.path("datafiles", "Bcell_imm_data.csv"), row.names = F)

imm_data %>% 
  filter(days.post.imm == 4) %>%
  summarise("CAR_countsMZ" = log(mean(CAR_MZB_numbers)),
            "countsGC" = log(mean(GCB_cell_numbers)),
            "CAR_countsGC" = log(mean(CAR_GCB_numbers)))

# generating data for fitting
imm_N2ko_data <- B_cell_data %>% filter(genotype == "N2KO") %>%
  select(days.post.imm, contains("MZ"), contains("GC"), -contains("fraction"), -total_MZBs) %>%
  mutate(fraction_CAR_MZ = CAR_MZB_numbers/MZB_cell_numbers,
         fraction_CAR_GC = CAR_GCB_numbers/GCB_cell_numbers,
         CAR_Neg_MZB = MZB_cell_numbers - CAR_MZB_numbers)%>%
  select(days.post.imm, CAR_MZB_numbers, GCB_cell_numbers, CAR_GCB_numbers, MZB_cell_numbers)%>%
  filter(days.post.imm > 0)

imm_N2ko_data %>% 
  filter(days.post.imm == 0) %>%
  summarise("CAR_countsMZ" = log(mean(CAR_MZB_numbers)),
            "countsGC" = log(mean(GCB_cell_numbers)),
            "CAR_countsGC" = log(mean(CAR_GCB_numbers)))


write.csv(imm_N2ko_data, file.path("datafiles", "N2KO_imm_data.csv"), row.names = F)

imm_N2ko_plot <- imm_N2ko_data %>%
  gather(- days.post.imm, key="subpop", value="cell_numbers")

imm_WT_plot <- imm_data %>%
  gather(- days.post.imm, key="subpop", value="cell_numbers")

ggplot() +
  geom_point(data = imm_WT_plot,
             aes(x= (days.post.imm), y=cell_numbers), col=2) +
  geom_point(data = imm_N2ko_plot,
             aes(x= (days.post.imm), y=cell_numbers), col=4) +
  scale_y_log10(limits = c(1e3, 1e7)) + xlim(0, 30)+
  labs(x='Days post immunization', y='Cell counts')+
  facet_wrap(.~ subpop) + myTheme + guides(col="none")



## Unique time points with indices to map
unique_times_df1 <- imm_data %>% distinct(days.post.imm, .keep_all = TRUE) 
data_times1 <- imm_data$days.post.imm 
solve_times1 <- unique_times_df1$days.post.imm  ## unique time points in the data
## Map of the unique time points on all the timepoints
time_index1 <- purrr::map_dbl(data_times1, function(x) which(x == solve_times1))    # keeping track of index of time point in relation to solve_time

## Unique time points with indices to map
unique_times_df2 <- imm_N2ko_data %>% distinct(days.post.imm, .keep_all = TRUE) 
data_times2 <- imm_N2ko_data$days.post.imm 
solve_times2 <- unique_times_df2$days.post.imm  ## unique time points in the data
## Map of the unique time points on all the timepoints
time_index2 <- purrr::map_dbl(data_times2, function(x) which(x == solve_times2))    # keeping track of index of time point in relation to solve_time

## Data to import in Stan
numObs1 <- length(imm_data$days.post.imm)
numObs2 <- length(imm_N2ko_data$CAR_MZB_numbers)
n_shards1 <- length(solve_times1)
solve_time1 <- solve_times1
time_index1 <- time_index1
n_shards2 <- length(solve_times2)
solve_time2 <- solve_times2
time_index2 <- time_index2
CAR_MZ_counts <- imm_data$CAR_MZB_numbers
CAR_MZN2_counts <- imm_N2ko_data$CAR_MZB_numbers
GC_counts <- imm_data$GCB_cell_numbers
CAR_GC_counts <- imm_data$CAR_GCB_numbers
CAR_GCN2_counts <- imm_N2ko_data$CAR_GCB_numbers
GCN2_counts <- imm_N2ko_data$GCB_cell_numbers

# time sequence for predictions specific to age bins within the data
ts_pred <- seq(4, 30, length.out = 500)
numPred <- length(ts_pred)


rstan::stan_rdump(c("numObs1",  "n_shards1", "solve_time1", "time_index1",
                    "numObs2",  "n_shards2", "solve_time2", "time_index2",
             "CAR_MZ_counts",  "GC_counts", "CAR_GC_counts", 
             "CAR_MZN2_counts", "CAR_GCN2_counts", "GCN2_counts",
             "ts_pred", "numPred"),
           file = file.path('datafiles', paste0('Bcell_Imm',".Rdump")))

###############################
##############################

## testing models from stan
library(rstan)

stanmodel_file <- file.path("stan_models", "totFOB_MZB_timeinflux.stan")
expose_stan_functions(stanmodel_file)


ts_seq <- c(4, 7, 14, 30)
init_cond <- c(exp(9.82), 0.299, exp(11.8))
params <- c(0.05, 0.02, 0.01, 0.04, 0.2, 0.1)


ts_pred <- seq(5, 30, length.out=100)
ode_df <- solve_ODE(ts_pred, init_cond, params)
stan_pred_df <- data.frame("time" = ts_pred,
                           "y_pred" = matrix(unlist(ode_df), nrow = length(ode_df), byrow = TRUE))%>%
  mutate(GC_counts = y_pred.3 + y_pred.2,
         CAR_GC = y_pred.2/GC_counts,
         CAR_MZ = y_pred.1)


ggplot(stan_pred_df)+
  geom_point(aes(x=time, y=GC_counts))+
  scale_y_log10()


carf_vec <- sapply(ts_pred, CAR_positive_FOB)
carf_vec2 <- sapply(ts_seq, CAR_positive_FOB)

ggplot() +
  geom_line(aes(x=ts_pred, y=carf_vec)) +
  geom_point(aes(x=ts_seq, y=carf_vec2))












