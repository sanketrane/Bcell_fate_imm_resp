rm(list = ls())  
gc()    
library(loo)
library(tidyverse)
library(gridExtra)
library(formattable)

## directories for saving outputs
OutputDir <- file.path('output_fit')

customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customOrange = "#FF8C00"
customRed = "#FF6347"

## function for exporting a data table of delta loo ic values and akaike weights
#takes 2 separate lists of the name of the models and loo-ic values
model_compare <- function(looiclist){
  # delta loo-ic
  deltaloo_list <- looiclist - min(looiclist)
  ## akaike wt
  akaikewt_numerator <- function(x) exp(-0.5 * x)
  akaikewt_list <- sapply(deltaloo_list, akaikewt_numerator) * 100/sum(exp(- 0.5 * deltaloo_list))
  
  export_table <- data.frame('deltaloo' = round(deltaloo_list, 2),
                             'Akaike_wt' = round(akaikewt_list, 2))
  colnames(export_table)[1:2] <-  c(paste0('\u0394', 'LooIC'),  paste0('Akaike weight ', '\u0025'))
  
  return(export_table)
}

LooDir <- file.path("loo_fit")

### reading the loo objects for each model
Branched_timeinflux_filename = paste0('loosave_Branched_timeinflux_Bcell_imm_data.csv.rds')
Branched_timeinflux1_filename = paste0('loosave_Branched_timeinflux1_Bcell_imm_data.csv.rds')
Branched_timeinflux2_filename = paste0('loosave_Branched_timeinflux2_Bcell_imm_data.csv.rds')
Branched_timeinflux3_filename = paste0('loosave_Branched_timeinflux3_Bcell_imm_data.csv.rds')
#Branched_timeinflux4_filename = paste0('loosave_Branched_timeinflux4_Bcell_imm_data.csv.rds')
#Branched_timeinflux5_filename = paste0('loosave_Branched_timeinflux5_Bcell_imm_data.csv.rds')
Linear_timeinflux_filename = paste0('loosave_Linear_timeinflux_Bcell_imm_data.csv.rds')
Linear_timeinflux1_filename = paste0('loosave_Linear_timeinflux1_Bcell_imm_data.csv.rds')
Linear_timeinflux2_filename = paste0('loosave_Linear_timeinflux2_Bcell_imm_data.csv.rds')
Linear_timeinflux3_filename = paste0('loosave_Linear_timeinflux3_Bcell_imm_data.csv.rds')
#Linear_timeinflux4_filename = paste0('loosave_Linear_timeinflux4_Bcell_imm_data.csv.rds')
#Linear_timeinflux5_filename = paste0('loosave_Linear_timeinflux5_Bcell_imm_data.csv.rds')
Null_timeinflux_filename = paste0('loosave_Null_timeinflux_Bcell_imm_data.csv.rds')
Null_timeinflux1_filename = paste0('loosave_Null_timeinflux1_Bcell_imm_data.csv.rds')
#Null_timeinflux2_filename = paste0('loosave_Null_timeinflux2_Bcell_imm_data.csv.rds')
Branched_neutral_filename = paste0('loosave_Branched_neutral_Bcell_imm_data.csv.rds')
Linear_neutral_filename = paste0('loosave_Linear_neutral_Bcell_imm_data.csv.rds')
Null_neutral_filename = paste0('loosave_Null_neutral_Bcell_imm_data.csv.rds')
Branched_timeloss_filename = paste0('loosave_Branched_timeloss_Bcell_imm_data.csv.rds')
Linear_timeloss_filename = paste0('loosave_Linear_timeloss_Bcell_imm_data.csv.rds')
Null_timeloss_filename = paste0('loosave_Null_timeloss_Bcell_imm_data.csv.rds')


model_list <- list(#'Branched_timeinflux' = readRDS(file.path(LooDir, Branched_timeinflux_filename)), 
                  #'Branched_timeinflux1' = readRDS(file.path(LooDir, Branched_timeinflux1_filename)), 
                  #'Branched_timeinflux2' = readRDS(file.path(LooDir, Branched_timeinflux2_filename)), 
                  'Branched_timeinflux3' = readRDS(file.path(LooDir, Branched_timeinflux3_filename)), 
                  #'Linear_timeinflux' = readRDS(file.path(LooDir, Linear_timeinflux_filename)), 
                  #'Linear_timeinflux1' = readRDS(file.path(LooDir, Linear_timeinflux1_filename)), 
                  #'Linear_timeinflux2' = readRDS(file.path(LooDir, Linear_timeinflux2_filename)), 
                  'Linear_timeinflux3' = readRDS(file.path(LooDir, Linear_timeinflux3_filename)), 
                   #"Null_timeinflux" = readRDS(file.path(LooDir, Null_timeinflux_filename)), 
                   "Null_timeinflux1" = readRDS(file.path(LooDir, Null_timeinflux1_filename)), 
                   "Branched_neutral" = readRDS(file.path(LooDir, Branched_neutral_filename)), 
                   "Linear_neutral" = readRDS(file.path(LooDir, Linear_neutral_filename)), 
                   "Null_neutral" = readRDS(file.path(LooDir, Null_neutral_filename)), 
                   "Branched_timeloss" = readRDS(file.path(LooDir, Branched_timeloss_filename)), 
                   "Linear_timeloss" = readRDS(file.path(LooDir, Linear_timeloss_filename)), 
                   "Null_timeloss" = readRDS(file.path(LooDir, Null_timeloss_filename)))
compare_mods <- loo_compare(model_list)
print(compare_mods, simplify = F)

compare_mods
loo_model_weights(model_list, method = 'pseudobma') * 100


