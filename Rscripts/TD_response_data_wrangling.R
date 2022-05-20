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
                 legend.background = element_blank(), legend.key = element_blank())

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
NEW_CAR_prop_df <- read_excel("datafiles/NEW_CAR_fractions_TD_response.xlsx", sheet =1) %>%
  filter(genotype == "CAR") %>%
  select(-contains("counts"), -contains("PCs")) %>%
  gather(-c(days.post.imm, genotype), key = 'subpop', value='fraction_cells')

ggplot() +
  geom_point(data=NEW_CAR_prop_df, 
             aes(x= as.factor(days.post.imm), y=fraction_cells), col = 4) +
  scale_y_log10() + 
  labs(x='Days post immunization', y='% CAR+')+
  facet_wrap(.~ subpop) + myTheme + guides(col="none")



## Immune response dynamics of B cells -- reporter induction upon b cell activation
B_cell_data <- read_excel("datafiles/NEW_CAR_fractions_TD_response.xlsx", sheet =1) %>%
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


B_cell_counts_df <- B_cell_data %>% 
  select(-contains("CAR"))%>%
  select(days.post.imm, genotype, contains('numbers'), - B_cell_numbers) %>%
  gather(-c(days.post.imm, genotype), key = 'subpop', value='cell_numbers') %>% na.omit()

CAR_cell_counts_df <- B_cell_data %>% 
  select(days.post.imm, genotype, contains("CAR"))%>%
  select(days.post.imm, genotype, contains('numbers'), - CAR_B_numbers) %>%
  gather(-c(days.post.imm, genotype), key = 'subpop', value='cell_numbers') %>% na.omit()


CAR_PosNeg_df <- B_cell_data %>% 
  select(days.post.imm, genotype, 
         contains("FO"), contains("GC"), contains("MZ"),
         -contains('total'), -contains('fraction')) %>%
  mutate(CAR_Neg_FOB = FOB_cell_numbers - CAR_FOB_numbers,
         CAR_Neg_MZB = MZB_cell_numbers - CAR_MZB_numbers,
         CAR_Neg_GCB = GCB_cell_numbers - CAR_GCB_numbers) #%>% 
  #select(days.post.imm, genotype, contains("CAR")) %>%
  #gather(-c(days.post.imm, genotype), key = 'subpop', value='cell_numbers') %>% na.omit()


facet_labels1 <- c(`B_cell_numbers` = "Total B cells",
                  `FOB_cell_numbers` = "FO B cells",
                  `MZB_cell_numbers` = "MZ B cells",
                  `GCB_cell_numbers` = "GC B cells",
                  `PCB_cell_numbers` = "PC B cells")

ggplot() +
  geom_point(data = filter(B_cell_counts_df, genotype == "CAR"),
             aes(x= (days.post.imm), y=cell_numbers, col = subpop)) +
  scale_y_log10(limits = c(1e4, 1e8)) + 
  labs(x='Days post immunization', y='Cell counts')+
  facet_wrap(.~ subpop, labeller = as_labeller(facet_labels1)) + myTheme + guides(col="none")


p1 <- ggplot() +
  geom_point(data = filter(B_cell_counts_df, genotype == "CAR", subpop == "FOB_cell_numbers"),
             aes(x= (days.post.imm), y=cell_numbers), col=2, size=2) +
  scale_y_log10(limits = c(5e6, 1e8), breaks=c(1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) + 
  labs(x='Days post immunization', y='Cell counts', title = "FO B cells") + myTheme

p2 <- ggplot() +
  geom_point(data = filter(B_cell_counts_df, genotype == "CAR", subpop == "MZB_cell_numbers"),
             aes(x= (days.post.imm), y=cell_numbers), col=6, size=2) +
  scale_y_log10(limits = c(1e5, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) + 
  labs(x='Days post immunization', y='Cell counts', title = "MZ B cells")+ myTheme


p3 <- ggplot() +
  geom_point(data = filter(B_cell_counts_df, genotype == "CAR", subpop == "GCB_cell_numbers"),
             aes(x= (days.post.imm), y=cell_numbers), col=4, size=2) +
  scale_y_log10(limits = c(1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) + 
  labs(x='Days post immunization', y='Cell counts', title = "GC B cells")+ myTheme



facet_labels2 <- c(`CAR_B_numbers` = "Total B cells",
                   `CAR_FOB_numbers` = "FO B cells",
                   `CAR_MZB_numbers` = "MZ B cells",
                   `CAR_GCB_numbers` = "GC B cells",
                   `CAR_PCB_numbers` = "PC B cells")

ggplot() +
  geom_point(data = filter(CAR_cell_counts_df, genotype == "CAR"),
             aes(x= (days.post.imm), y=cell_numbers, col = subpop)) +
  scale_y_log10(limits = c(1e3, 1e7)) +
  labs(x='Days post immunization', title ='Number of CAR+ Cells')+
  facet_wrap(.~ subpop, labeller = as_labeller(facet_labels2)) + myTheme + guides(col="none")


p1.1 <- ggplot() +
  geom_point(data = filter(CAR_cell_counts_df, genotype == "CAR", subpop == "CAR_FOB_numbers"),
             aes(x= (days.post.imm), y=cell_numbers), col=2, size=2) +
  scale_y_log10(limits = c(1e3, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y='Number of CAR+ Cells') + myTheme


p2.1 <- ggplot() +
  geom_point(data = filter(CAR_cell_counts_df, genotype == "CAR", subpop == "CAR_MZB_numbers"),
             aes(x= (days.post.imm), y=cell_numbers), col=6, size=2) +
  scale_y_log10(limits = c(1e3, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y='Number of CAR+ Cells') + myTheme


p3.1 <- ggplot() +
  geom_point(data = filter(CAR_cell_counts_df, genotype == "CAR", subpop == "CAR_GCB_numbers"),
             aes(x= (days.post.imm), y=cell_numbers), col=4, size=2) +
  scale_y_log10(limits = c(1e3, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y='Number of CAR+ Cells') + myTheme



## saving  plots for quality control 
pdf(file = file.path("plots", paste("DataPlots%03d.pdf", sep = "")),
    width = 12, height = 7.5, onefile = FALSE, useDingbats = FALSE)
cowplot::plot_grid(p1, p3, p2, p1.1, p3.1, p2.1,
                   nrow  = 2)
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

repFOB_nlm <- nls(log(CAR_FOB_numbers) ~ log(phi_func(days.post.imm, basl, theta, n, X, q)),
                  data = CAR_PosNeg_df,
                  start = list(basl = 11.2, theta = 7.5, n = 3.5 , X = 8, q = 5))

par_est <- coef(repFOB_nlm)
sol_time <- seq(0, 30, length.out=100)
phi_vec_m <- phi_func(sol_time, basl=par_est[1], theta=par_est[2], n=par_est[3],
                    X = par_est[4], q=par_est[5])

ggplot() +
  geom_line(aes(x=sol_time, y=(phi_vec_m)), col=2, size=0.8) +
  geom_point(data=CAR_PosNeg_df,
             aes(x=days.post.imm, y=(CAR_FOB_numbers)), size=2, col=2) +
  scale_y_log10(limits = c(1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y=NULL, title='Number of CAR+ FO B Cells') + myTheme


ggplot() +
  #geom_line(aes(x=sol_time, y=exp(phi_vec_m) + exp(17.1527)), col=2, size=0.8) +
  #geom_point(data=CAR_PosNeg_df,aes(x=days.post.imm, y=(CAR_Neg_FOB)), size=2, col=2)+
  geom_point(data=CAR_PosNeg_df,
             aes(x=days.post.imm, y=(FOB_cell_numbers)), size=2, col=4)+
  scale_y_log10(limits = c(1e6, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y=NULL, title='Number of CAR- FO B Cells') + myTheme



### CAR neg MZ B cells as precursors

phi_func <- function(Time, basl, nu, b0){
  t = Time - 0
  basl_exp = exp(basl)
  return( basl_exp * (1 + exp(-nu * (Time - b0)^2)))
}


Neg_MZB_nlm <- nls(log(CAR_Neg_MZB) ~ log(phi_func(days.post.imm,  basl, nu, b0)),
                  data = CAR_PosNeg_df,
                  start = list(basl = 11, nu = 0.01, b0=8))

par_est <- coef(Neg_MZB_nlm)
sol_time <- seq(0, 30, length.out=100)
phi_vec_m <- phi_func(sol_time, basl=par_est[1], nu=par_est[2], b0=par_est[3])
phi_vec_mm <- phi_func(sol_time, basl=14, nu=0.01, b0=18)

ggplot() +
  #geom_line(aes(x=sol_time, y=(phi_vec_m)), col=2, size=0.8) +
  geom_line(aes(x=sol_time, y=(phi_vec_mm)), col=2, size=0.8) +
  geom_point(data=CAR_PosNeg_df,
             aes(x=days.post.imm, y=(CAR_Neg_MZB)), size=2, col=2) +
  scale_y_log10(limits = c(1e5, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y=NULL, title='Number of CAR- MZ B Cells') + myTheme


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
numObs1 <- length(data_times)
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
ts_pred <- seq(4, 30, length.out = 300)
numPred <- length(ts_pred)

ggplot() +
  geom_point(aes(data_times, log(CAR_GC_counts)), col=2)


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












