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
## Notch2 induction using TAM in CD19hom
N2_hom_df <- read_excel("Datafiles/TAM_all_samples.xlsx", sheet =2) %>%
  select(mouse.id, days.post.TAM, contains('hCD2')) %>% na.omit() %>%
  mutate(hcd2_prop = hCD2/100,
         mzb_prop = (MZBs_in_hCD2/100)*hcd2_prop) %>%
  select(mouse.id, days.post.TAM, contains('prop')) %>%
  #filter(days.post.TAM != 21) %>%
  gather(c(hcd2_prop, mzb_prop), key = subpop, value=prop_B)

## Notch2 induction using TAM in CD19het 
N2_het_df <- read_excel("Datafiles/TAM_all_samples.xlsx", sheet =3) %>%
  select(mouse.id, days.post.TAM, contains('hCD2')) %>% na.omit()

legend_labels <- c(`hcd2_prop` = "hCD2+",
                   `mzb_prop` = "MZB")

  
## plots
ggplot(N2_hom_df) +
  geom_boxplot(aes(x= as.factor(days.post.TAM), y=prop_B*100, col = subpop)) +
  scale_color_discrete(name=NULL, labels=legend_labels)+
  scale_y_continuous(trans="log10", breaks=c(0.1, 0.3, 1, 3, 10), 
                     minor_breaks = log10minorbreaks) +
  labs(x='Days post TAM treatment', y='% of total B cells') 


## plots
ggplot(N2_hom_df) +
  geom_jitter(aes(x= (days.post.TAM), y=prop_B*100, col = subpop), width = 0.2) +
  scale_color_discrete(name=NULL, labels=legend_labels)+
  scale_y_continuous(trans="log10", breaks=c(0.1, 0.3, 1, 3, 10), 
                     minor_breaks = log10minorbreaks) +
  labs(x='Days post TAM treatment', y='% of total B cells') 

############################################################
############################################################

## conversion of Adoptively transferred of FO B cells into MZ B cells
adopt_tra_df <- read_excel("Datafiles/TAM_all_samples.xlsx", sheet =4) %>%na.omit() %>%
  select(-contains('host')) %>%
  mutate(donor_prop = donor_in_Bcells/100,
         mzbdonor_prop = (mzb_in_donor/100),
         Nfd = mzbdonor_prop) %>%
  select(days.post.transfer, contains('prop')) %>%
  gather(- days.post.transfer, key=subpop, value= prop_B)


facet_labels <- c(`donor_prop` = "% total donor in B cells",
                  `mzbdonor_prop` = "% MZB cells in donor compartment")

ggplot(adopt_tra_df, aes(x= as.factor(days.post.transfer), y=prop_B*100, col = subpop)) +
  geom_boxplot() + geom_jitter(width = 0.1, alpha=0.5)+
  #scale_y_log10() + 
  labs(x='Days post trabsfer', y=NULL)+
  facet_wrap(.~ subpop, labeller = as_labeller(facet_labels)) +
  myTheme + guides(col="none")

ggplot(adopt_tra_df, aes(x= (days.post.transfer), y=prop_B*100, col = subpop)) +
  geom_point()+ xlim(0, 15)+
  scale_y_log10() + 
  labs(x='Days post trabsfer', y=NULL)+
  facet_wrap(.~ subpop, labeller = as_labeller(facet_labels)) +
  myTheme + guides(col="none")



############################################################
############################################################

## Immune response dynamics of B cells -- reporter induction upon b cell activation
car_prop_df <- read_excel("Datafiles/MZB_TDimmune_response.xlsx", sheet =3)  %>%
  select( - contains('total'), -contains('MZP')) %>% #filter(genotype == "CAR")  %>%
  gather(- c(days.post.imm, genotype), key=subpop, value= prop_B)

facet_labels <- c(`rep_pos_B` = "Total B cells",
                  `rep_pos_MZP` = "MZP cells",
                  `rep_pos_MZB` = "MZB cells")

ggplot(car_prop_df) +
  geom_boxplot(aes(x= as.factor(days.post.imm), y=prop_B, col = subpop)) +
  scale_y_log10() + 
  labs(x='Days post immunization', y='% CAR+')+
  facet_wrap(.~ subpop, labeller = as_labeller(facet_labels)) + myTheme + guides(col="none")

## Immune response dynamics of B cells -- reporter induction upon b cell activation
car_counts_df <- read_excel("Datafiles/MZB_TDimmune_response.xlsx", sheet =3)  %>%
  select(-contains('MZP')) %>%
  mutate(num_totB = (total_B/100) * total_spleen_counts *1e6,
         num_rep_posB= (rep_pos_B/100) * num_totB,
         num_rep_pos_MZB = (rep_pos_MZB/100) * num_rep_posB) %>%
  select(days.post.imm, genotype, contains('num_rep')) %>%
  gather(- c(days.post.imm, genotype), key=subpop, value=cell_num)

legend_labels <- c(`num_rep_posB` = "Total B",
                  `num_rep_pos_MZB` = "MZB")

ggplot(car_counts_df, aes(x= as.factor(days.post.imm), y=cell_num, col=subpop)) +
  geom_boxplot() + geom_jitter(width = 0.2, alpha=0.5)+
  scale_color_discrete(name=NULL, labels=legend_labels)+
  scale_y_continuous(trans="log10", breaks=c(1e4, 1e5, 1e6), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y='# of Rep+ cells') +
  facet_wrap(.~ genotype)



legend_labels <- c(`num_rep_posB` = "Total B",
                   `num_rep_pos_MZB` = "MZB")

ggplot(car_prop_df, aes(x= as.factor(days.post.imm), y=prop_B, col=genotype)) +
  geom_boxplot() + geom_jitter(width = 0.2, alpha=0.5)+
  scale_color_discrete(name=NULL, labels=legend_labels)+
  #scale_y_continuous(trans="log10", breaks=c(1e4, 1e5, 1e6), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x='Days post immunization', y='Proportion of Rep+ cells') +
  facet_wrap(.~ subpop)


tot_B_df <- car_counts_df %>%
  filter(genotype == "CAR") %>% filter(subpop == "num_rep_posB") %>% 
  select(days.post.imm, cell_num)

phi_func <- function(t, b0, nu){
 10^b0 * t * exp(-nu * t)
  # b0 *  exp(-nu * t)
}

sol_time <- seq(0.1, max(tot_B_df$days.post.imm), length.out=200)

phi_vec <- sapply(sol_time, phi_func, b0=5.8,  nu=0.1)

ggplot() +
  geom_line(aes(x=sol_time, y=phi_vec))+
  geom_point(data=tot_B_df, aes(x=days.post.imm, y=cell_num))+
  scale_y_log10()


repB_nlm <- nls(cell_num ~ phi_func(days.post.imm, b0, nu), data = tot_B_df,
                start = list(b0=6, nu=0.1))

par_est <- coef(repB_nlm)


phi_vec <- sapply(sol_time, phi_func, b0=par_est[1], nu=par_est[2])
phi_vec_fid <- sapply(sol_time, phi_func, b0=6, nu=0.131)

ggplot() +
  #geom_line(aes(x=sol_time, y=phi_vec))+
  geom_line(aes(x=sol_time, y=phi_vec_fid), size=1.0, col=4)+
  geom_point(data=tot_B_df, aes(x=days.post.imm, y=cell_num), col=4)+
  scale_y_log10(limits=c(1e4, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  labs(x="Days post immunization", y=NULL, title="Fraction of CAR+ cells in Total B cells") +
  myTheme


phi_funct <- function(t, psi){
  b0=6; nu=0.131
 exp(psi) * 10^b0 * t * exp(-nu * t)
}

integran_func <- function(Time, alpha, psi){
  alpha_transf = exp(alpha)
  phi_funct(Time, psi) * exp(-alpha_transf * Time)
}

integ_func <- function(Time, alpha, psi){
  integrate(integran_func, 0, Time, alpha=alpha, psi=psi)$value
}

ode_func <- function(Time, alpha, psi, M0){
  M0_transf = 10^M0
 #y0 * exp(alpha * Time) + 
    (exp(alpha * Time)/M0_transf) * integ_func(Time, alpha, psi)
}

ode_vec <- Vectorize(ode_func)


MZ_B_df <- car_prop_df %>%
  filter(genotype == "CAR") %>% filter(subpop == "rep_pos_MZB") %>% 
  select(days.post.imm, prop_B) %>%
  mutate(rep_pos_frac = prop_B/100)

asin_transf <- function(x){asin(sqrt(x))}

asin_transf(MZ_B_df$rep_pos_frac)

#model optimisation 
LL_MZB <- function(param, boot_data) {
  alpha  <- param[1]             #parametrs to be estimated as part of a vector
  psi  <- param[2]
  M0  <- param[3]
  
  k <- length(param)        #numner of unknown parameters 
  n <-  nrow(boot_data)
  
  pred <- ode_vec(boot_data$days.post.imm, alpha, psi, M0) 
  dat <-  boot_data$rep_pos_frac
  
  R1 <- sum((asin_transf(pred) - asin_transf(dat))^2)
  
  logl <- -(n/2)*(log(R1)) #-(n/2)*(log(R2))
  
  aiccFP <<- -2*logl + 2*k*n/(n-k-1)
  #n1=n2=n since number of obeservations isnt different
  
  return(-logl)
} 

fit_LL_MZB <- optim(par=c(-3, -8, 6), fn=LL_MZB, boot_data=MZ_B_df,
                    method = "Nelder-Mead",
                    control = list(trace = 6))
fit_LL_MZB
par_est <- fit_LL_MZB$par
aiccFP

ode_sol_df <- ode_vec(sol_time, par_est[1],par_est[2], par_est[3])
ode_sol_m <- ode_vec(sol_time, 0.08, -3, par_est[3])

ggplot() +
  geom_line(aes(x=sol_time, y=ode_sol_df), size=1.0, col=4)+
  #geom_line(aes(x=sol_time, y=ode_sol_m), size=1.0, col=2)+
  geom_point(data=MZ_B_df, aes(x=days.post.imm, y=rep_pos_frac), col=4)+
  labs(x="Days post immunization", y=NULL, title="Fraction of CAR+ cells in MZ compartment") +
  myTheme



