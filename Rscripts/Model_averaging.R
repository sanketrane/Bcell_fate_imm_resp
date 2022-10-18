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
MZB_const_filename = paste0('loosave_Branched_timeinflux_Bcell_imm_data.csv.rds')
totFOB_const_filename = paste0('loosave_Linear_timeinflux_Bcell_imm_data.csv.rds')
carGC_const_filename = paste0('loosave_Null_timeinflux_Bcell_imm_data.csv.rds')
MZB_timeloss_filename = paste0('loosave_Branched_neutral_Bcell_imm_data.csv.rds')
totFOB_timeloss_filename = paste0('loosave_Linear_neutral_Bcell_imm_data.csv.rds')
carGC_timeloss_filename = paste0('loosave_Null_neutral_Bcell_imm_data.csv.rds')
MZB_timeinflux_filename = paste0('loosave_Branched_timeloss_Bcell_imm_data.csv.rds')
totFOB_timeinflux_filename = paste0('loosave_Linear_timeloss_Bcell_imm_data.csv.rds')
carGC_timeinflux_filename = paste0('loosave_Null_timeloss_Bcell_imm_data.csv.rds')

MZB_const_loo <- readRDS(file.path(LooDir, MZB_const_filename))
totFOB_const_loo <- readRDS(file.path(LooDir, totFOB_const_filename))
carGC_const_loo <- readRDS(file.path(LooDir, carGC_const_filename))
MZB_timeloss_loo <- readRDS(file.path(LooDir, MZB_timeloss_filename))
totFOB_timeloss_loo <- readRDS(file.path(LooDir, totFOB_timeloss_filename))
carGC_timeloss_loo <- readRDS(file.path(LooDir, carGC_timeloss_filename))
MZB_timeinflux_loo <- readRDS(file.path(LooDir, MZB_timeinflux_filename))
totFOB_timeinflux_loo <- readRDS(file.path(LooDir, totFOB_timeinflux_filename))
carGC_timeinflux_loo <- readRDS(file.path(LooDir, carGC_timeinflux_filename))

model_list <- list('Branched_timeinflux' = MZB_const_loo, 
                   'Linear_timeinflux' = totFOB_const_loo, 
                   "Null_timeinflux" = carGC_const_loo, 
                   "Branched_neutral" = MZB_timeloss_loo,
                   "Linear_neutral" = totFOB_timeloss_loo,
                   "Null_neutral" = carGC_timeloss_loo,
                   "Branched_timeloss" = MZB_timeinflux_loo,
                   "Linear_timeloss" = totFOB_timeinflux_loo,
                   "Null_timeloss" = carGC_timeinflux_loo)
compare_mods <- loo_compare(model_list)
print(compare_mods, simplify = F)

compare_mods
formattable::formattable(
mw_vec <- loo_model_weights(model_list, method = 'pseudobma') * 100
)
data.frame(mw_vec)

