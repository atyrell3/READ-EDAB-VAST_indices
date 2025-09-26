
# Script Header -----------------------------------------------------------

# Title: Create VAST Inputs for Zooplankton Index Indicators
# Author: Adelle Molina
# Date: September 2025
# Description: This script processes ECOMON plankton data 
# For the 2025 SOE, steps to add OISST excluded (as is done for forage) since not used in the models
# The primary output is a dataset aggregated by station with columns for species groupings at various depths of interest 
# It is suitable for VAST modeling or other analyses of
# zooplankton index indicators for the State of the Ecosystem report.

# see decisions in https://noaa-edab.github.io/zooplanktonindex/ZoopVASTworkflow.html
# and detailed workflow description here https://noaa-edab.github.io/zooplanktonindex/ZoopSOEworkflow.html

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(here)
library(dplyr)
library(rlang)

# Data --------------------------------------------------------------------

# read in ecomon plankton data, this is what needs to be updated 
ecomonall <- read_csv(here::here("zooplankton/inputs/EcoMon_Plankton_Data_v3_8.csv"))
#ecomonall <- read_csv(here::here("data/EcoMon Data_NCEI/EcoMon_Plankton_Data_v3_09_Do_Not_Distribute.csv"))

# SST data - from forage index script output --> 
# this must be run afterwards & not even used here plus, the stations will be different 
# OISST.dat <- readRDS(here("data-raw/dietstn_OISST_1982_2023.rds"))

# Summarize -----------------------------------------------------------------

#new dataset names all caps
names(ecomonall) <- stringr::str_to_lower(names(ecomonall))

# Define column groups 
copepods_10m2 <- c("ctyp_10m2", "calfin_10m2", "pseudo_10m2", "tlong_10m2", "cham_10m2", "para_10m2", 
                   "acarspp_10m2", "mlucens_10m2", "calminor_10m2", "clauso_10m2", "acarlong_10m2", 
                   "euc_10m2", "fur_10m2", "calspp_10m2", "ost_10m2", "temspp_10m2", "tort_10m2", "paraspp_10m2")
copepods_100m3 <- c("ctyp_100m3", "calfin_100m3", "pseudo_100m3", "tlong_100m3", "cham_100m3", "para_100m3",
                    "acarspp_100m3", "mlucens_100m3", "calminor_100m3", "clauso", "acarlong_100m3", "euc_100m3",
                    "fur_100m3", "calspp_100m3", "ost_100m3", "temspp_100m3", "tort_100m3", "paraspp_100m3")

smallcope_all_10m2 <- c("ctyp_10m2", "pseudo_10m2", "tlong_10m2", "cham_10m2", "para_10m2", 
                        "acarspp_10m2", "clauso_10m2", "acarlong_10m2", "fur_10m2", "ost_10m2",
                        "temspp_10m2", "tort_10m2", "paraspp_10m2")
smallcope_all_100m3 <- c("ctyp_100m3", "pseudo_100m3", "tlong_100m3", "cham_100m3", "para_100m3", 
                         "acarspp_100m3", "clauso", "acarlong_100m3", "fur_100m3", "ost_100m3", 
                         "temspp_100m3", "tort_100m3", "paraspp_100m3")

lgcope_all_10m2 <- c("calfin_10m2", "mlucens_10m2", "calminor_10m2", "euc_10m2", "calspp_10m2")
lgcope_all_100m3 <- c("calfin_100m3", "mlucens_100m3", "calminor_100m3", "euc_100m3", "calspp_100m3")

smallcope_soe_10m2 <- c("ctyp_10m2", "pseudo_10m2", "tlong_10m2", "cham_10m2")
smallcope_soe_100m3 <- c("ctyp_100m3", "pseudo_100m3", "tlong_100m3", "cham_100m3")

# sum across all copepods
sumcopepods <- ecomonall |>
  filter(!is.na(cruise_name)) |>
  rowwise() |>
  mutate(
    allcopepods_10m2 = sum(c_across(all_of(copepods_10m2)), na.rm = TRUE),
    allcopepods_100m3 = sum(c_across(all_of(copepods_100m3)), na.rm = TRUE),
    smallcopeALL_10m2 = sum(c_across(all_of(smallcope_all_10m2)), na.rm = TRUE),
    smallcopeALL_100m3 = sum(c_across(all_of(smallcope_all_100m3)), na.rm = TRUE),
    smallcopeSOE_10m2 = sum(c_across(all_of(smallcope_soe_10m2)), na.rm = TRUE),
    smallcopeSOE_100m3 = sum(c_across(all_of(smallcope_soe_100m3)), na.rm = TRUE),
    lgcopeALL_10m2 = sum(c_across(all_of(lgcope_all_10m2)), na.rm = TRUE),
    lgcopeALL_100m3 = sum(c_across(all_of(lgcope_all_100m3)), na.rm = TRUE),
    lgcopeSOE_10m2 = calfin_10m2,
    lgcopeSOE_100m3 = calfin_100m3
  ) |>
  ungroup() |>
  select(cruise_name, station, lat, lon, date, time, depth, sfc_temp, sfc_salt,
         btm_temp, btm_salt, volume_1m2, volume_100m3, starts_with("allcopepods_"),
         starts_with("smallcopeALL_"), starts_with("lgcopeALL_"),
         starts_with("smallcopeSOE_"), starts_with("lgcopeSOE_"))

# add summed copepod columns, fix dates, establish seasons
copepod_dat <- ecomonall |>
  dplyr::filter(!is.na(cruise_name)) |>
  dplyr::left_join(sumcopepods) |>
  dplyr::select(cruise_name, station, lat, lon, date, time, depth, 
                sfc_temp, sfc_salt, btm_temp, btm_salt, volume_1m2, volume_100m3,
                allcopepods_10m2, allcopepods_100m3, euph_10m2, euph_100m3,
                hyper_10m2, hyper_100m3, calfin_10m2, calfin_100m3,
                smallcopeALL_10m2, smallcopeALL_100m3,
                lgcopeALL_10m2, lgcopeALL_100m3,
                smallcopeSOE_10m2, smallcopeSOE_100m3,) |>
                #cluhar_10m2, cluhar_100m3) |>
  dplyr::mutate(date = lubridate::dmy(date),
                year = lubridate::year(date),
                month = lubridate::month(date),
                day = lubridate::day(date),
                stationid = paste0(cruise_name, "_", station),
                season_ng = case_when(month <= 6 ~ "SPRING",
                                      month >= 7 ~ "FALL",
                                      TRUE ~ as.character(NA)),
                season_4 = case_when(month %in% c(12,1,2) ~ "winter",
                                     month %in% c(3,4,5) ~ "spring",
                                     month %in% c(6,7,8) ~ "summer",
                                     month %in% c(9,10,11) ~ "fall",
                                     TRUE ~ as.character(NA)))

# Save out -----------------------------------------------------------------

#saveRDS(copepod_dat, here::here("zooplankton/outputs/test_VAST_input.rds"))
saveRDS(copepod_dat, here::here("zooplankton/outputs/copepod_VAST_input.rds"))

