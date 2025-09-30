# Script Header -----------------------------------------------------------

# Title: Run VAST models to create Zooplankton Index Indicators
# Author: Adelle Molina
# Date: September 2025
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

## trying to run this to recreate the 2025 index and set up a workflow to update indicators when new data is available
## 2025-09-024 AIM
## 2025-09-30 AST - added functions to create SOE formatted outputs. Zooplankton volume and copepod models still need to be run.

# Libraries and functions -----------------------------------------------------------
# Load packages
library(here)
library(dplyr)
library(VAST)
# rlang

# Functions

SOEinputs <- function(infile, season, stratlook, taxa, outfile) {
  splitoutput <- read.csv(infile)
  zoopindex <- splitoutput %>%
    left_join(stratlook) %>%
    dplyr::select(
      Time,
      EPU = Region,
      "Abundance Index Estimate" = Estimate,
      "Abundance Index Estimate SE" = Std..Error.for.Estimate
    ) %>%
    tidyr::pivot_longer(
      c("Abundance Index Estimate", "Abundance Index Estimate SE"),
      names_to = "Var",
      values_to = "Value"
    ) %>%
    dplyr::filter(EPU %in% c("MAB", "GB", "GOM", "AllEPU")) %>%
    dplyr::mutate(Units = "numbers per 100 cu m volume") %>%
    dplyr::select(Time, Var, Value, EPU, Units)

  zoopindex$Var <- stringr::str_c(season, taxa, zoopindex$Var, sep = " ")

  saveRDS(zoopindex, outfile)
}

extract_cog <- function(model_fit) {
  dir.create(here::here("temp"))
  cog <- FishStatsUtils::plot_range_index(
    Sdreport = model_fit$parameter_estimates$SD,
    Report = model_fit$Report,
    TmbData = model_fit$data_list,
    year_labels = as.numeric(model_fit$year_labels),
    years_to_plot = model_fit$years_to_plot,
    Znames = colnames(model_fit$data_list$Z_xm),
    PlotDir = here::here("temp")
  ) #already have plots, will delete this directory
  unlink(here::here("temp"), recursive = TRUE)
  return(cog)
}

SOEinputsCOG <- function(infile, season, taxa, outfile) {
  cogout <- readRDS(infile)
  zoocog <- as.data.frame(cogout$COG_Table) |>
    dplyr::mutate(direction = ifelse(m == 1, "Eastward", "Northward")) |>
    dplyr::select(
      "Time" = Year,
      "Center of Gravity" = COG_hat,
      "Center of Gravity SE" = SE,
      direction
    ) |>
    tidyr::pivot_longer(
      c("Center of Gravity", "Center of Gravity SE"),
      names_to = "Var",
      values_to = "Value"
    ) |>
    #direction into Var
    tidyr::unite(Var, direction:Var, sep = " ") |>
    dplyr::mutate(Units = "km", EPU = "ALLEPU") |>
    dplyr::select(Time, Var, Value, EPU, Units)

  zoocog$Var <- stringr::str_c(
    stringr::str_to_title(season),
    stringr::str_to_title(taxa),
    zoocog$Var,
    sep = " "
  )
  saveRDS(zoocog, outfile)
}

create_soe_indices <- function(
  dir, # directory names where model results are saved
  strat = stratlook_EPUonly
) {
  simple_file <- basename(dir)

  split <- strsplit(simple_file, split = "_")

  season <- split[[1]][1]
  taxa <- split[[1]][2]

  ## center of gravity
  cog_file <- paste0(dir, "/", season, "_cog.rds")

  SOEinputsCOG(
    infile = cog_file,
    season = season,
    taxa = taxa,
    outfile = paste0(
      here::here("zooplankton/outputs/Indices"),
      "/",
      season,
      "_",
      taxa,
      "_",
      "cog",
      "_",
      Sys.Date(),
      ".rds"
    )
  )

  ## index
  index_file <- paste0(dir, "/Index.csv")
  SOEinputs(
    infile = index_file,
    season = season,
    taxa = taxa,
    stratlook = strat,
    outfile = paste0(
      here::here("zooplankton/outputs/Indices"),
      "/",
      season,
      "_",
      taxa,
      "_",
      "index",
      "_",
      Sys.Date(),
      ".rds"
    )
  )
}
## Load main dataset ----

#Read in main data created in zooindex_process_input_data.R
copepod_dat <- readRDS(here::here(
  "zooplankton/outputs/zooplankton_VAST_input.rds"
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

# default for index2: works for calfin and lgcopeALL
ObsModel1 <- c(2, 1) # this is "Index2", Gamma distribution for positive catches and Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number

# alternative for models where we encounter stuff everywhere (zooplankton volume, maybe small copepods)
ObsModel2 <- c(2, 4) # should be Gamma distribution for positive catches and Poisson-link fixing encounter probability=1 for any year where all samples encounter the species

# Model selection 1 (spatial, spatio-temporal effects, no covariates) options and naming:
# selected _alleffectson for all seasons and copeopod groups
# Season_knots + suffix below
# _alleffectson             FieldConfig default (all IID)

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
  #"her_sp" = her_spr,
  #"her_fa" = her_fall,
  "MAB" = MAB2,
  "GB" = GB2,
  "GOM" = GOM2,
  "SS" = SS2
))


# this is run within the loop if you run the models all at once for each taxa, here run it first and use for all models
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
# used for zooplankton volume
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
# This can be used for euph and zoopvol
configs <- list(
  ann = list(
    season_filter = FALSE
  ),
  fall = list(
    season_filter = "FALL"
  ),
  spring = list(
    season_filter = "SPRING"
  )
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
      AreaSwept_km2 = 1, #Elizabeth's code
      #declon = -declon already done before neamap merge
      Vessel = 1,
      Dayofyear = lubridate::yday(date) #as.numeric(as.factor(vessel))-1
    ) %>%
    dplyr::select(
      Catch_g = euph_100m3, #use megabenwt for individuals input in example
      Year = year,
      Vessel,
      AreaSwept_km2,
      Lat = lat,
      Lon = lon,
      # sstfill,
      Dayofyear
    ) %>%
    na.omit() %>%
    as.data.frame()

  # Assign the resulting data frame to the looping variable name
  assign(paste0("euph_stn_", cfg_name), df)
}


## Euphausiid model runs (simultaneous, slower, do not run)######################
# Run bias corrected models for annual, spring, and fall simultaneously

# list of data, settings, and directory for output for each option

mod.season <- c("euph_fall_500", "euph_spring_500", "euph_ann_500") #includes n knots

mod.dat <- list(euph_stn_fall, euph_stn_spring, euph_stn_ann)

names(mod.dat) <- mod.season

mod.obsmod <- list(ObsModel1, ObsModel1, ObsModel1)

names(mod.obsmod) <- mod.season

# Define covariate combinations

mod.covar <- c("biascorrect", "biascorrect_doy")

for (season in mod.season) {
  season <- season

  dat <- mod.dat[[season]]

  Q_ikbase <- NULL
  Q_ikdoy <- as.matrix(dat[, c("Dayofyear")])

  mod.Qik <- list(Q_ikbase, Q_ikdoy)

  names(mod.Qik) <- mod.covar

  for (covar in mod.covar) {
    name <- paste0(season, "_", covar)

    working_dir <- here::here(sprintf("pyindex/%s/", name))

    if (!dir.exists(working_dir)) {
      dir.create(working_dir)
    }

    ObsModel <- mod.obsmod[[season]]
    Q_ik <- mod.Qik[[covar]]

    settings <- make_settings(
      n_x = 500,
      Region = "northwest_atlantic",
      Version = "VAST_v14_0_1", #needed to prevent error from newer dev version number
      #strata.limits = list('All_areas' = 1:1e5), full area
      strata.limits = strata.limits,
      purpose = "index2",
      ObsModel = ObsModel,
      bias.correct = TRUE,
      use_anisotropy = use_anisotropy,
      FieldConfig = FieldConfig,
      RhoConfig = RhoConfig, #always default
      OverdispersionConfig = OverdispersionConfig
    )

    fit <- try(fit_model(
      settings = settings,
      #extrapolation_list = New_Extrapolation_List,
      Lat_i = dat$Lat,
      Lon_i = dat$Lon,
      t_i = dat$Year,
      b_i = as_units(dat[, 'Catch_g'], 'g'),
      a_i = rep(1, nrow(dat)),
      v_i = dat$Vessel,
      Q_ik = Q_ik,
      #Use_REML = TRUE,
      #test_fit = FALSE,
      working_dir = paste0(working_dir, "/")
    ))

    saveRDS(fit, file = paste0(working_dir, "/fit.rds"))

    # Plot results
    if (!class(fit) == "try-error") {
      plot(fit, working_dir = paste0(working_dir, "/"))
    }
  } # end config loop
} # end season loop


## Euphausiid model runs (separately by season, no annual) ---------------------------------------

# Fall first, no DOY or temp covariate
working_dir <- here::here(sprintf("zooplankton/outputs/fall_euph_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = euph_stn_fall$Lat,
  Lon_i = euph_stn_fall$Lon,
  t_i = euph_stn_fall$Year,
  b_i = as_units(euph_stn_fall[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(euph_stn_fall)),
  v_i = euph_stn_fall$Vessel,
  #Q_ik = as.matrix(euph_stn_fall[,c("Catch_g")]),
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
saveRDS(cog, here::here(working_dir, paste0("fall_cog.rds")))

###

# run spring models, also with no DOY or temp covariates
working_dir <- here::here(sprintf("zooplankton/outputs/spring_euph_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = euph_stn_spring$Lat,
  Lon_i = euph_stn_spring$Lon,
  t_i = euph_stn_spring$Year,
  b_i = as_units(euph_stn_spring[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(euph_stn_spring)),
  v_i = euph_stn_spring$Vessel,
  #Q_ik = as.matrix(euph_stn_spring[,c("Catch_g")]),
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

# Zooplankton volume ----

## TODO: zoop models are erroring out with message:
# Error in (function (b_i, a_i, t_i, c_iz = rep(0, length(b_i)), e_i = c_iz[,  :
# Some years and/or categories have 100% encounters, and this requires either temporal structure of a different link-function

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
      AreaSwept_km2 = 1, #Elizabeth's code
      #declon = -declon already done before neamap merge
      Vessel = 1,
      Dayofyear = lubridate::yday(date) #as.numeric(as.factor(vessel))-1
    ) %>%
    dplyr::select(
      Catch_g = volume_100m3,
      Year = year,
      Vessel,
      AreaSwept_km2,
      Lat = lat,
      Lon = lon,
      # sstfill,
      Dayofyear
    ) %>%
    na.omit() %>%
    as.data.frame()

  # Assign the resulting data frame to the looping variable name
  assign(paste0("zoopvol_stn_", cfg_name), df)
}

## Zooplankton Volume model runs (separately by season, no annual) ---------------------------------------

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

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot(fit, working_dir = paste0(working_dir, "/"))

# extract center of gravity
# had to run this function separately....issue with the utils script that wont' run right
cog <- extract_cog(fit)
saveRDS(cog, here::here(working_dir, paste0("fall_cog.rds")))

###

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

# Copepods ----

## Copepod data prep ---- ----##################################################################

# --- New Configuration for Taxa and Season ---

# Define the list of copepod taxa short names and their corresponding column names
# should probably add argument for euph and zoopvol to do this all at once in one loop at the top of the script
taxa_cols <- list(
  "smcope" = "smallcopeSOE_100m3",
  "calfin" = "calfin_100m3",
  "lgcope" = "lgcopeALL_100m3"
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
        AreaSwept_km2 = 1, #Elizabeth's code
        Vessel = 1,
        Dayofyear = lubridate::yday(date) #as.numeric(as.factor(vessel))-1
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
        # sstfill,
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
## small Copepod models ----

# Fall first, no DOY or temp covariate
working_dir <- here::here(sprintf("zooplankton/outputs/fall_smcope_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = smcope_fall$Lat,
  Lon_i = smcope_fall$Lon,
  t_i = smcope_fall$Year,
  b_i = as_units(smcope_fall[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(smcope_fall)),
  v_i = smcope_fall$Vessel,
  #Q_ik = as.matrix(smcope_fall[,c("Catch_g")]),
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
saveRDS(cog, here::here(working_dir, paste0("fall_cog.rds")))

###

# run spring models, also with no DOY or temp covariates
working_dir <- here::here(sprintf("zooplankton/outputs/spring_smcope_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = smcope_spring$Lat,
  Lon_i = smcope_spring$Lon,
  t_i = smcope_spring$Year,
  b_i = as_units(smcope_spring[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(smcope_spring)),
  v_i = smcope_spring$Vessel,
  #Q_ik = as.matrix(smcope_spring[,c("Catch_g")]),
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

## large Copepod models ----

# Fall first, no DOY or temp covariate
working_dir <- here::here(sprintf("zooplankton/outputs/fall_lgcope_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = lgcope_fall$Lat,
  Lon_i = lgcope_fall$Lon,
  t_i = lgcope_fall$Year,
  b_i = as_units(lgcope_fall[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(lgcope_fall)),
  v_i = lgcope_fall$Vessel,
  #Q_ik = as.matrix(lgcope_fall[,c("Catch_g")]),
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
saveRDS(cog, here::here(working_dir, paste0("fall_cog.rds")))

###

# run spring models, also with no DOY or temp covariates
working_dir <- here::here(sprintf("zooplankton/outputs/spring_lgcope_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = lgcope_spring$Lat,
  Lon_i = lgcope_spring$Lon,
  t_i = lgcope_spring$Year,
  b_i = as_units(lgcope_spring[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(lgcope_spring)),
  v_i = lgcope_spring$Vessel,
  #Q_ik = as.matrix(lgcope_spring[,c("Catch_g")]),
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

## calfin models ----

# Fall first, no DOY or temp covariate
working_dir <- here::here(sprintf("zooplankton/outputs/fall_calfin_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = calfin_fall$Lat,
  Lon_i = calfin_fall$Lon,
  t_i = calfin_fall$Year,
  b_i = as_units(calfin_fall[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(calfin_fall)),
  v_i = calfin_fall$Vessel,
  #Q_ik = as.matrix(calfin_fall[,c("Catch_g")]),
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
saveRDS(cog, here::here(working_dir, paste0("fall_cog.rds")))

###

# run spring models, also with no DOY or temp covariates
working_dir <- here::here(sprintf("zooplankton/outputs/spring_calfin_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = calfin_spring$Lat,
  Lon_i = calfin_spring$Lon,
  t_i = calfin_spring$Year,
  b_i = as_units(calfin_spring[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(calfin_spring)),
  v_i = calfin_spring$Vessel,
  #Q_ik = as.matrix(calfin_spring[,c("Catch_g")]),
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
## map COG code over all combinations
list_dirs <- paste0(
  here::here("zooplankton/outputs"),
  "/",
  rep(c("spring", "fall"), each = 5),
  "_",
  rep(
    c(
      # commented these out for testing
      # all models have not run yet
      "euph" #, "zoopvol", "smcope", "lgcope", "calfin"
    ),
    2
  ),
  "_model"
)

# create_soe_indices(dir = list_dirs[1])

purrr::map(unique(list_dirs), ~ create_soe_indices(dir = .x))
