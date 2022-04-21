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
  geom_point(data=NEW_NEW_CAR_prop_df, 
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
  select(days.post.imm, genotype, contains('numbers')) %>%
  gather(-c(days.post.imm, genotype), key = 'subpop', value='cell_numbers') %>% na.omit()

CAR_cell_counts_df <- B_cell_data %>% 
  select(days.post.imm, genotype, contains("CAR"))%>%
  select(days.post.imm, genotype, contains('numbers')) %>%
  gather(-c(days.post.imm, genotype), key = 'subpop', value='cell_numbers') %>% na.omit()


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

facet_labels2 <- c(`CAR_B_numbers` = "Total B cells",
                   `CAR_FOB_numbers` = "FO B cells",
                   `CAR_MZB_numbers` = "MZ B cells",
                   `CAR_GCB_numbers` = "GC B cells",
                   `CAR_PCB_numbers` = "PC B cells")

ggplot() +
  geom_point(data = filter(CAR_cell_counts_df, genotype == "CAR"),
             aes(x= (days.post.imm), y=cell_numbers, col = subpop)) +
  scale_y_log10(limits = c(1e3, 1e7)) +
  labs(x='Days post immunization', y='Number of CAR+ Cells')+
  facet_wrap(.~ subpop, labeller = as_labeller(facet_labels2)) + myTheme + guides(col="none")


###############################
##############################

# Modeling Precursor Dynamics 

phi_func <- function(t, b0, r, nu){
  tau = 7
  #exp(b0) * (1 +  t * exp(-nu * t))
  exp(b0) * exp(r *t)/(1 +  exp(nu * (t - tau)))
}

sol_time <- seq(0, 30, length.out=200)

repFOB_nlm <- nls((fraction_cells/100) ~ phi_func(days.post.imm, b0,r, nu),
                  data = filter(NEW_CAR_prop_df, subpop == "fraction_CAR_in_FOBs"),
                  start = list(b0=-3, r=0.01, nu=0.5))

par_est <- coef(repFOB_nlm)

phi_vec_m <- sapply(sol_time, phi_func, b0=par_est[1], r=par_est[2], nu=par_est[3])

ggplot() +
  #geom_line(aes(x=sol_time, y=phi_vec))+
  geom_line(aes(x=sol_time, y=phi_vec_m), col=2)+
  geom_point(data=filter(NEW_CAR_prop_df, subpop == "fraction_CAR_in_FOBs"),
             aes(x=days.post.imm, y=fraction_cells/100))



###############################
##############################

# generating data for fitting
imm_data <- B_cell_data %>% filter(genotype == "CAR") %>%
  select(days.post.imm, contains("MZ"), contains("GC"), -contains("fraction"), -total_MZBs) %>%
  mutate(fraction_CAR_MZ = CAR_MZB_numbers/MZB_cell_numbers,
         fraction_CAR_GC = CAR_GCB_numbers/GCB_cell_numbers)%>%
  select(days.post.imm, MZB_cell_numbers, GCB_cell_numbers, fraction_CAR_MZ, fraction_CAR_GC)


imm_data %>% 
  filter(days.post.imm == 0) %>%
  summarise("countsMZ" = log(mean(MZB_cell_numbers)),
            "fractionMZ" = mean(fraction_CAR_MZ),
            "countsGC" = log(mean(GCB_cell_numbers)),
            "fractionGC" = mean(fraction_CAR_GC))

# generating data for fitting
imm_N2ko_data <- B_cell_data %>% filter(genotype == "N2KO") %>%
  select(days.post.imm, contains("MZ"), contains("GC"), -contains("fraction"), -total_MZBs) %>%
  mutate(fraction_CAR_MZ = CAR_MZB_numbers/MZB_cell_numbers,
         fraction_CAR_GC = CAR_GCB_numbers/GCB_cell_numbers)%>%
  select(days.post.imm, MZB_cell_numbers, GCB_cell_numbers, fraction_CAR_MZ, fraction_CAR_GC)


## Unique time points with indices to map
unique_times_df <- imm_data %>% distinct(days.post.imm, .keep_all = TRUE) 
data_times <- imm_data$days.post.imm 
solve_times <- unique_times_df$days.post.imm  ## unique time points in the data
## Map of the unique time points on all the timepoints
time_index <- purrr::map_dbl(data_times, function(x) which(x == solve_times))    # keeping track of index of time point in relation to solve_time

## Data to import in Stan
numObs <- length(data_times)
n_shards <- length(solve_times)
solve_time <- solve_times
time_index <- time_index
MZ_counts <- imm_data$MZB_cell_numbers
MZ_fractions <- imm_data$fraction_CAR_MZ
GC_counts <- imm_data$MZB_cell_numbers
GC_fractions <- imm_data$fraction_CAR_GC

# time sequence for predictions specific to age bins within the data
ts_pred <- seq(0, 30, length.out = 300)
numPred <- length(ts_pred)

logit_trans <- function(x) log(x/(1-x))
asin_trans <- function(x) asin(sqrt(x))

ggplot() +
  geom_point(aes(data_times, logit_trans(GC_fractions)), col=2)+
  #geom_point(aes(data_times, asin_trans(MZ_fractions)), col=4)+
  geom_point(aes(data_times, (MZ_fractions))) + 
  ylim(0, .42)



stan_rdump(c("numObs",  "n_shards", "solve_time", "time_index",
             "MZ_counts",  "MZ_fractions", "GC_counts", "GC_fractions",
             "ts_pred", "numPred"),
           file = file.path('datafiles', paste0('Bcell_Imm',".Rdump")))

###############################
##############################

## testing models from stan
library(rstan)

stanmodel_file <- file.path("stan_models", "direct_desc_time_act.stan")
expose_stan_functions(stanmodel_file)

ts_seq <- seq(0, 4, 7, 14, 30)
init_cond <- c(0.0136, 0.299, exp(11.8))
params <- c(0.05, 0.02, 0.01, 0.04, 0.2, 0.1)
solve_ODE(ts_seq, init_cond, params)












