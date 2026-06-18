# Script Header -----------------------------------------------------------

# Title: Run VAST models to create Zooplankton Index Indicators
# Author: AST
# Date: November 2025
# Description: This script runs many VAST models
# For the 2025 SOE, no day of year covariates or temperature have been included
# The primary output are various seasonal zooplankton index indicators for the State of the Ecosystem report.
# Spring and fall Euphausiids
# Spring and fall Zooplankton volume
#

# This is modified from original VAST code from S.Gaichas
# Use documentation here to debug code:
# https://noaa-edab.github.io/zooplanktonindex/ZoopSOEworkflow.html

# VAST attempt 2 univariate model as a script
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/Index-standardization

# Debugging a common error
# this error is caused by memory issues
# if you get this error, it can be solved by closing out other applications
# and re-running the model

# Error in TransformADFunObject(ADFun, method = "laplace", config = cfg,  :
# Caught exception 'std::bad_alloc' in function 'TransformADFunObjectTemplate'

## trying to run this to recreate the 2025 index and set up a workflow to update indicators when new data is available
## 2025-09-024 AIM
## 2025-09-30 AST - added functions to create SOE formatted outputs. Zooplankton volume and copepod models still need to be run.
## 2025-11-20 AST - ran fall models for everything except zooplankton volume, which is erroring out and not a priority to fix right now
## euphausiid, small and large copepod, calfin spring & fall models need to be re-run with updated data from Harvey

# Libraries and functions -----------------------------------------------------------
# Load packages
library(here)
library(dplyr)
library(VAST)
# rlang

# Functions
source(here::here("zooplankton/functions.R"))

## Load main dataset ----

#Read in main data created in zooindex_process_input_data.R
copepod_dat <- readRDS(here::here(
  "zooplankton/outputs/zooplankton_VAST_input_2025.rds"
))

## General Settings  ---- ----##################################################
# Make settings (turning off bias.correct to save time for example)
# NEFSC strata limits https://github.com/James-Thorson-NOAA/VAST/issues/302

MAB <- c(1010:1080, 1100:1120, 1600:1750, 3010:3450, 3470, 3500, 3510)
GB <- c(1090, 1130:1210, 1230, 1250, 3460, 3480, 3490, 3520:3550)
GOM <- c(1220, 1240, 1260:1290, 1360:1400, 3560:3830)
SS <- c(1300:1352, 3840:3990)

# Mid Atlantic Bight EPU
MAB2 <- FishStatsUtils::northwest_atlantic_grid %>%
  dplyr::filter(stratum_number %in% MAB) %>%
  dplyr::select(stratum_number) %>%
  dplyr::distinct()

# Georges Bank EPU
GB2 <- FishStatsUtils::northwest_atlantic_grid %>%
  dplyr::filter(stratum_number %in% GB) %>%
  dplyr::select(stratum_number) %>%
  dplyr::distinct()

# gulf of maine EPU (for SOE)
GOM2 <- FishStatsUtils::northwest_atlantic_grid %>%
  dplyr::filter(stratum_number %in% GOM) %>%
  dplyr::select(stratum_number) %>%
  dplyr::distinct()

# scotian shelf EPU (for SOE)
SS2 <- FishStatsUtils::northwest_atlantic_grid %>%
  dplyr::filter(stratum_number %in% SS) %>%
  dplyr::select(stratum_number) %>%
  dplyr::distinct()

# needed to cover the whole northwest atlantic grid--lets try without
# allother2 <- coast3nmbuffst %>%
#   dplyr::filter(!stratum_number %in% c(MAB, GB, GOM, SS)) %>%
#   dplyr::select(stratum_number2) %>%
#   dplyr::distinct()

# all epus
allEPU2 <- FishStatsUtils::northwest_atlantic_grid %>%
  dplyr::filter(stratum_number %in% c(MAB, GB, GOM, SS)) %>%
  dplyr::select(stratum_number) %>%
  dplyr::distinct()

## Model Configurations  ----############################################################

### create model parameters ----

# default for index2: works for calfin, lgcope, small cop, euph
ObsModel1 <- c(2, 1) # this is "Index2", Gamma distribution for positive catches and Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number

# alternative for models where we encounter stuff everywhere (zooplankton volume, maybe small copepods)
ObsModel2 <- c(2, 4) # should be Gamma distribution for positive catches and Poisson-link fixing encounter probability=1 for any year where all samples encounter the species

# default configs
FieldConfig = matrix(
  "IID",
  ncol = 2,
  nrow = 3,
  dimnames = list(
    c("Omega", "Epsilon", "Beta"),
    c("Component_1", "Component_2")
  )
)

RhoConfig <- c(
  "Beta1" = 0, # temporal structure on years (intercepts)
  "Beta2" = 0,
  "Epsilon1" = 0, # temporal structure on spatio-temporal variation
  "Epsilon2" = 0
)
# not testing alternative RhoConfigs here just noted for completeness
# 0 off (fixed effects)
# 1 independent
# 2 random walk
# 3 constant among years (fixed effect)
# 4 AR1

use_anisotropy <- TRUE

OverdispersionConfig <- c("eta1" = 0, "eta2" = 0)
# eta1 = vessel effects on prey encounter rate
# eta2 = vessel effects on prey weight

strata.limits <- as.list(c(
  "AllEPU" = allEPU2,
  "MAB" = MAB2,
  "GB" = GB2,
  "GOM" = GOM2,
  "SS" = SS2
))

### create model settings ----
# first set of settings
# used for euphausiids
settings <- make_settings(
  n_x = 500,
  Region = "northwest_atlantic",
  Version = "VAST_v14_0_1", #needed to prevent error from newer dev version number
  #strata.limits = list('All_areas' = 1:1e5), full area
  strata.limits = strata.limits,
  purpose = "index2",
  ObsModel = c(2, 1), # this is "Index2", Gamma distribution for positive catches and Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number
  bias.correct = TRUE,
  use_anisotropy = use_anisotropy,
  FieldConfig = FieldConfig,
  RhoConfig = RhoConfig, #always default
  OverdispersionConfig = OverdispersionConfig
)

# second set of settings
# used for zooplankton volume, small copepod spring
settings2 <- make_settings(
  n_x = 500,
  Region = "northwest_atlantic",
  Version = "VAST_v14_0_1", #needed to prevent error from newer dev version number
  #strata.limits = list('All_areas' = 1:1e5), full area
  strata.limits = strata.limits,
  purpose = "index2",
  ObsModel = c(2, 4), # should be Gamma distribution for positive catches and Poisson-link fixing encounter probability=1 for any year where all samples encounter the species
  bias.correct = TRUE,
  use_anisotropy = use_anisotropy,
  FieldConfig = FieldConfig,
  RhoConfig = RhoConfig, #always default
  OverdispersionConfig = OverdispersionConfig
)

# Euphausiids ----

## Euphasiid data prep ---- ----##################################################################

#Create three separate data frames using a single loop

# Create a list to hold the configurations for each variable (annual vs seasonal)
# annual models aren't run at this point in time, but keeping for future use
# This can be used for euph and zoopvol
configs <- list(
  ann = list(season_filter = FALSE),
  fall = list(season_filter = "FALL"),
  spring = list(season_filter = "SPRING")
)

# Loop through the configurations to generate each data frame
for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]

  # Start with the base data and filter for the year
  df <- copepod_dat %>%
    filter(year > 1981)

  # Apply the season filter if it is specified in the configuration
  if (cfg$season_filter != FALSE) {
    df <- df %>% filter(season_ng == cfg$season_filter)
  }

  df <- df %>%
    mutate(
      AreaSwept_km2 = 1,
      Vessel = 1,
      Dayofyear = lubridate::yday(date)
    ) %>%
    dplyr::select(
      Catch_g = euph_100m3,
      Year = year,
      Vessel,
      AreaSwept_km2,
      Lat = lat,
      Lon = lon,
      Dayofyear
    ) %>%
    na.omit() %>%
    as.data.frame()

  # Assign the resulting data frame to the looping variable name
  assign(paste0("euph_stn_", cfg_name), df)
}

## Euphausiid model runs (separately by season, no annual) ---------------------------------------

### Fall ----
# no DOY or temp covariate
run_vast_model(
  data = euph_stn_fall,
  out_dir = here::here(sprintf("zooplankton/outputs/fall_euph_model")),
  season = "fall",
  vast_settings = settings
)

### spring ----
# run spring model
## TODO: run this model with a DOY covariate
## Sarah's model selection showed that the DOY model converged and had lower AIC than the default model
## if model fails, remove DOY covariate
run_vast_model(
  data = euph_stn_spring,
  out_dir = here::here(sprintf("zooplankton/outputs/spring_euph_model")),
  season = "spring",
  vast_settings = settings,
  doy = TRUE
)

# Copepods ----

## Copepod data prep ---- ----##################################################################

# --- New Configuration for Taxa and Season ---

# Define the list of copepod taxa short names and their corresponding column names
# should probably add argument for euph and zoopvol to do this all at once in one loop at the top of the script
taxa_cols <- list(
  "smcope" = "smallcopeSOE_100m3",
  "calfin" = "calfin_100m3",
  "lgcope" = "lgcopeALL_100m3",
  "smcopeall" = "smallcopeALL_100m3"
)

# Define the season configurations (No annual)
season_configs <- list(
  fall = "FALL",
  spring = "SPRING"
)

# --- Data Processing Loop ---

# Outer loop: Iterate over each taxa
for (taxa_name in names(taxa_cols)) {
  taxa_col_name <- taxa_cols[[taxa_name]]

  # Inner loop: Iterate over each season
  for (season_name in names(season_configs)) {
    season_filter_value <- season_configs[[season_name]]

    # Build the data frame for the current taxa and season
    df <- copepod_dat %>%
      # Filter for the specific season and year
      filter(
        season_ng == season_filter_value, # Only Fall or Spring seasons
        year > 1981
      ) %>%
      # Apply common mutations
      mutate(
        AreaSwept_km2 = 1,
        Vessel = 1,
        Dayofyear = lubridate::yday(date)
      ) %>%
      # Dynamically select the Catch_g column based on the current taxa
      dplyr::select(
        # Use !!rlang::sym() for dynamic column selection
        Catch_g = !!rlang::sym(taxa_col_name),
        Year = year,
        Vessel,
        AreaSwept_km2,
        Lat = lat,
        Lon = lon,
        Dayofyear
      ) %>%
      # Apply final steps: remove NA values and convert to a data frame
      na.omit() %>%
      as.data.frame()

    # Assign the resulting data frame using the taxa name and season name
    # e.g., 'euph_fall', 'calfin_spring'
    assign(paste0(taxa_name, "_", season_name), df)
  }
}

## Copepod model runs (separately by season, no annual) ---------------------------------------

## small SOE Copepod models ----

### fall ----
# Fall first, no DOY or temp covariate
run_vast_model(
  data = smcope_fall,
  out_dir = here::here(sprintf("zooplankton/outputs/fall_smcope_model")),
  season = "fall",
  vast_settings = settings
)

### spring ----
# run spring models, also with no DOY or temp covariates
run_vast_model(
  data = smcope_spring,
  out_dir = here::here(sprintf("zooplankton/outputs/spring_smcope_model")),
  season = "spring",
  vast_settings = settings2 # different settings because 100% encounters
)

## small ALL Copepod models ----

### fall ----
# Fall first, no DOY or temp covariate
# did not run in 2025 SOE, but including here for completeness
run_vast_model(
  data = smcopeall_fall,
  out_dir = here::here(sprintf("zooplankton/outputs/fall_smcopeall_model")),
  season = "fall",
  vast_settings = settings2 # different settings because 100% encounters
)

### spring ----
# run spring models, also with no DOY or temp covariates
run_vast_model(
  data = smcopeall_spring,
  out_dir = here::here(sprintf("zooplankton/outputs/spring_smcopeall_model")),
  season = "spring",
  vast_settings = settings2 # different settings because 100% encounters
)

## large Copepod models ----

### fall ----
# Fall first, no DOY or temp covariate
run_vast_model(
  data = lgcope_fall,
  out_dir = here::here(sprintf("zooplankton/outputs/fall_lgcope_model")),
  season = "fall",
  vast_settings = settings
)

### spring ----
# run spring model
## TODO: run this model with a DOY covariate
## Sarah's model selection showed that the DOY model converged and had lower AIC than the default model
## if model fails, remove DOY covariate
run_vast_model(
  data = lgcope_spring,
  out_dir = here::here(sprintf("zooplankton/outputs/spring_lgcope_model")),
  season = "spring",
  vast_settings = settings,
  doy = TRUE
)

## calfin models ----

### fall ----
# Fall first, no DOY or temp covariate
run_vast_model(
  data = calfin_fall,
  out_dir = here::here(sprintf("zooplankton/outputs/fall_calfin_model")),
  season = "fall",
  vast_settings = settings
)

### spring ----
# run spring models, also with no DOY or temp covariates
run_vast_model(
  data = calfin_spring,
  out_dir = here::here(sprintf("zooplankton/outputs/spring_calfin_model")),
  season = "spring",
  vast_settings = settings
)

# Create SOE indices ----

# set up areas
stratlook_EPUonly <- data.frame(
  Stratum = c(
    "Stratum_1",
    "Stratum_2",
    "Stratum_3",
    "Stratum_4",
    "Stratum_5"
  ),
  Region = c("AllEPU", "MAB", "GB", "GOM", "SS")
)

## map SOE index code ----
## map code over all models with list of directory names
## this list_dirs needs to be updated to include spring and zoopvol models (when results exist)
list_dirs <- list.dirs(here::here("zooplankton/outputs")) |>
  stringr::str_subset(pattern = "model") |>
  stringr::str_subset(pattern = "zoopvol", negate = TRUE) # exclude zoopvol for now since those models are erroring out
# create_soe_indices(dir = list_dirs[1])

purrr::map(
  unique(list_dirs),
  ~ create_soe_indices(dir = .x, strat = stratlook_EPUonly)
)

######

# Zooplankton volume ----
## DO NOT RUN -- ERRORS

## Zooplankton Volume data prep ---- ----##################################################################

#Create three separate data frames using a single loop
# use same seasonal configurations as for euphausiids

# Loop through the configurations to generate each data frame
for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]

  # Start with the base data and filter for the year
  df <- copepod_dat %>%
    filter(year > 1981)

  # Apply the season filter if it is specified in the configuration
  if (cfg$season_filter != FALSE) {
    df <- df %>% filter(season_ng == cfg$season_filter)
  }

  df <- df %>%
    mutate(
      AreaSwept_km2 = 1,
      Vessel = 1,
      Dayofyear = lubridate::yday(date)
    ) %>%
    dplyr::select(
      Catch_g = volume_100m3,
      Year = year,
      Vessel,
      AreaSwept_km2,
      Lat = lat,
      Lon = lon,
      Dayofyear
    ) %>%
    na.omit() %>%
    as.data.frame()

  # Assign the resulting data frame to the looping variable name
  assign(paste0("zoopvol_stn_", cfg_name), df)
}

## Zooplankton Volume model runs (separately by season, no annual) ---------------------------------------

### fall ----
# Fall first, no DOY or temp covariate
working_dir <- here::here(sprintf("zooplankton/outputs/fall_zoopvol_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings2,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = zoopvol_stn_fall$Lat,
  Lon_i = zoopvol_stn_fall$Lon,
  t_i = zoopvol_stn_fall$Year,
  b_i = as_units(zoopvol_stn_fall[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(zoopvol_stn_fall)),
  v_i = zoopvol_stn_fall$Vessel,
  #Q_ik = as.matrix(zoopvol_stn_fall[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/")
)

# ## ERRORS WITH MESSAGE:
# The following parameters appear to be approaching zero:
#   Param starting_value Lower          MLE Upper final_gradient
# 45 L_omega2_z              1  -Inf 6.259014e-06   Inf    0.001255733
# Please turn off factor-model variance parameters `L_` that are approaching zero and re-run the model
#
#
# Error: Please change model structure to avoid problems with parameter estimates and then re-try; see details in `?check_fit`

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot(fit, working_dir = paste0(working_dir, "/"))

# extract center of gravity
# had to run this function separately....issue with the utils script that wont' run right
cog <- extract_cog(fit)
saveRDS(cog, here::here(working_dir, paste0("fall_cog.rds")))

### spring ----
# run spring models, also with no DOY or temp covariates
working_dir <- here::here(sprintf("zooplankton/outputs/spring_zoopvol_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = zoopvol_stn_spring$Lat,
  Lon_i = zoopvol_stn_spring$Lon,
  t_i = zoopvol_stn_spring$Year,
  b_i = as_units(zoopvol_stn_spring[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(zoopvol_stn_spring)),
  v_i = zoopvol_stn_spring$Vessel,
  #Q_ik = as.matrix(zoopvol_stn_spring[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/")
)

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot(fit, working_dir = paste0(working_dir, "/"))

# extract center of gravity
cog <- extract_cog(fit)
saveRDS(cog, here::here(working_dir, paste0("spring_cog.rds")))
