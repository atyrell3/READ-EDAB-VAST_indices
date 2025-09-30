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

# VAST attempt 2 univariate model as a script
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/Index-standardization

## trying to run this to recreate the 2025 index and set up a workflow to update indicators when new data is available
## 2025-09-024 AIM

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
    dplyr::select(Time, 
                  EPU = Region, 
                  "Abundance Index Estimate" = Estimate, 
                  "Abundance Index Estimate SE" = Std..Error.for.Estimate) %>%
    tidyr::pivot_longer(c("Abundance Index Estimate", "Abundance Index Estimate SE"), 
                        names_to = "Var", values_to = "Value") %>%
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
    dplyr::mutate(direction = ifelse(m==1, "Eastward", "Northward")) |>
    dplyr::select("Time" = Year,
                  "Center of Gravity" = COG_hat, 
                  "Center of Gravity SE" = SE,
                  direction) |>
    tidyr::pivot_longer(c("Center of Gravity", "Center of Gravity SE"), 
                        names_to = "Var", values_to = "Value") |>
    #direction into Var
    tidyr::unite(Var, direction:Var, sep = " ") |>
    dplyr::mutate(Units = "km",
                  EPU = "ALLEPU") |>
    dplyr::select(Time, Var, Value, EPU, Units)
  
  zoocog$Var <- stringr::str_c(stringr::str_to_title(season), stringr::str_to_title(taxa), zoocog$Var, sep = " ")
  saveRDS(zoocog, outfile)
}
## Load main dataset ----

#Read in main data created in zooindex_process_input_data.R
copepod_dat <- readRDS(here::here("zooplankton/outputs/zooplankton_VAST_input.rds"))

## General Settings  ---- ----##################################################
# Make settings (turning off bias.correct to save time for example)
# NEFSC strata limits https://github.com/James-Thorson-NOAA/VAST/issues/302

MAB <- c(1010:1080, 1100:1120, 1600:1750, 3010:3450, 3470, 3500, 3510)
GB  <- c(1090, 1130:1210, 1230, 1250, 3460, 3480, 3490, 3520:3550)
GOM <- c(1220, 1240, 1260:1290, 1360:1400, 3560:3830)
SS  <- c(1300:1352, 3840:3990)

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
ObsModel1  <-  c(2,1) # this is "Index2", Gamma distribution for positive catches and Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number

# alternative for models where we encounter stuff everywhere (zooplankton volume, maybe small copepods)
ObsModel2  <-  c(2,4) # should be Gamma distribution for positive catches and Poisson-link fixing encounter probability=1 for any year where all samples encounter the species

# Model selection 1 (spatial, spatio-temporal effects, no covariates) options and naming:
# selected _alleffectson for all seasons and copeopod groups
# Season_knots + suffix below
# _alleffectson             FieldConfig default (all IID)

# default configs
FieldConfig = matrix( "IID", ncol=2, nrow=3, 
                      dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")))

RhoConfig <- c(
  "Beta1" = 0,      # temporal structure on years (intercepts) 
  "Beta2" = 0, 
  "Epsilon1" = 0,   # temporal structure on spatio-temporal variation
  "Epsilon2" = 0
) 
# not testing alternative RhoConfigs here just noted for completeness
# 0 off (fixed effects)
# 1 independent
# 2 random walk
# 3 constant among years (fixed effect)
# 4 AR1

use_anisotropy <- TRUE

OverdispersionConfig	<- c("eta1"=0, "eta2"=0)
# eta1 = vessel effects on prey encounter rate
# eta2 = vessel effects on prey weight

strata.limits <- as.list(c("AllEPU" = allEPU2,
                           #"her_sp" = her_spr,
                           #"her_fa" = her_fall,
                           "MAB" = MAB2,
                           "GB" = GB2,
                           "GOM" = GOM2,
                           "SS" = SS2
))


# this is run within the loop if you run the models all at once for each taxa, here run it first and use for all models
settings <- make_settings( n_x = 500, 
                           Region = "northwest_atlantic",
                           Version = "VAST_v14_0_1", #needed to prevent error from newer dev version number
                           #strata.limits = list('All_areas' = 1:1e5), full area
                           strata.limits = strata.limits,
                           purpose = "index2", 
                           bias.correct = TRUE,
                           use_anisotropy = use_anisotropy,
                           FieldConfig = FieldConfig,
                           RhoConfig = RhoConfig, #always default
                           OverdispersionConfig = OverdispersionConfig
)

## Euphasiid data prep ---- ----##################################################################

#Create three separate data frames from the input data using a single loop

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

for(season in mod.season){
  
  season <- season 
  
  dat <- mod.dat[[season]]
  
  Q_ikbase  <-  NULL
  Q_ikdoy <- as.matrix(dat[,c("Dayofyear")])
  
  mod.Qik <- list(Q_ikbase, Q_ikdoy)
  
  names(mod.Qik) <- mod.covar
  
  for(covar in mod.covar) {
    
    name <- paste0(season,"_", covar)
    
    working_dir <- here::here(sprintf("pyindex/%s/", name))
    
    if(!dir.exists(working_dir)) {
      dir.create(working_dir)
    }
    
    ObsModel <- mod.obsmod[[season]]
    Q_ik <- mod.Qik[[covar]]
    
    settings <- make_settings( n_x = 500, 
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
      b_i = as_units(dat[,'Catch_g'], 'g'),
      a_i = rep(1, nrow(dat)),
      v_i = dat$Vessel,
      Q_ik = Q_ik,
      #Use_REML = TRUE,
      #test_fit = FALSE,
      working_dir = paste0(working_dir, "/"))
    )
    
    saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
    
    # Plot results
    if(!class(fit)=="try-error"){
      plot( fit,
            working_dir = paste0(working_dir, "/"))
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
  b_i = as_units(euph_stn_fall[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(euph_stn_fall)),
  v_i = euph_stn_fall$Vessel,
  #Q_ik = as.matrix(euph_stn_fall[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

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
  b_i = as_units(euph_stn_spring[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(euph_stn_spring)),
  v_i = euph_stn_spring$Vessel,
  #Q_ik = as.matrix(euph_stn_spring[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

# extract center of gravity
cog <- extract_cog(fit)
saveRDS(cog, here::here(working_dir, paste0("spring_cog.rds"))) 

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
  settings = settings, 
  #extrapolation_list = New_Extrapolation_List,
  Lat_i = zoopvol_stn_fall$Lat, 
  Lon_i = zoopvol_stn_fall$Lon, 
  t_i = zoopvol_stn_fall$Year, 
  b_i = as_units(zoopvol_stn_fall[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(zoopvol_stn_fall)),
  v_i = zoopvol_stn_fall$Vessel,
  #Q_ik = as.matrix(zoopvol_stn_fall[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

# extract center of gravity
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
  b_i = as_units(zoopvol_stn_spring[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(zoopvol_stn_spring)),
  v_i = zoopvol_stn_spring$Vessel,
  #Q_ik = as.matrix(zoopvol_stn_spring[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

# extract center of gravity
cog <- extract_cog(fit)
saveRDS(cog, here::here(working_dir, paste0("spring_cog.rds"))) 

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
  b_i = as_units(smcope_fall[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(smcope_fall)),
  v_i = smcope_fall$Vessel,
  #Q_ik = as.matrix(smcope_fall[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

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
  b_i = as_units(smcope_spring[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(smcope_spring)),
  v_i = smcope_spring$Vessel,
  #Q_ik = as.matrix(smcope_spring[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

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
  b_i = as_units(lgcope_fall[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(lgcope_fall)),
  v_i = lgcope_fall$Vessel,
  #Q_ik = as.matrix(lgcope_fall[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

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
  b_i = as_units(lgcope_spring[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(lgcope_spring)),
  v_i = lgcope_spring$Vessel,
  #Q_ik = as.matrix(lgcope_spring[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

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
  b_i = as_units(calfin_fall[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(calfin_fall)),
  v_i = calfin_fall$Vessel,
  #Q_ik = as.matrix(calfin_fall[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

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
  b_i = as_units(calfin_spring[,'Catch_g'], 'g'),
  a_i = rep(1, nrow(calfin_spring)),
  v_i = calfin_spring$Vessel,
  #Q_ik = as.matrix(calfin_spring[,c("Catch_g")]),
  #Use_REML = TRUE,
  #test_fit = FALSE,
  working_dir = paste0(working_dir, "/"))

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
#fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot( fit, working_dir = paste0(working_dir, "/"))

# extract center of gravity
cog <- extract_cog(fit)
saveRDS(cog, here::here(working_dir, paste0("spring_cog.rds"))) 

## Create SOE indices - hardcoded ----
# this should probably be a loop

# set up areas
stratlook_EPUonly <- data.frame(
  Stratum = c(
    "Stratum_1",
    "Stratum_2",
    "Stratum_3",
    "Stratum_4",
    "Stratum_5"),
  Region = c("AllEPU", "MAB", "GB", "GOM", "SS")
)
  ## Fall indices ----
# Euphausiids
SOEinputs(
  infile = "zooplankton/outputs/fall_euph_model/Index.csv",
  season = "Fall",
  taxa = "Euph",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/fallEuphindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/fall_euph_model/fall_cog.rds",
  season = "Fall",
  taxa = "Euph",
  outfile = paste0("zooplankton/outputs/Indices/fallEuphcog", Sys.Date(), ".rds"))

# Zooplankton volume
SOEinputs(
  infile = "zooplankton/outputs/fall_zoopvol_model/Index.csv",
  season = "Fall",
  taxa = "Zoopvol",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/fallZoopvolindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/fall_zoopvol_model/fall_cog.rds",
  season = "Fall",
  taxa = "Zoopvol",
  outfile = paste0("zooplankton/outputs/Indices/fallZoopvolcog", Sys.Date(), ".rds"))

## small Copepod
SOEinputs(
  infile = "zooplankton/outputs/fall_smcope_model/Index.csv",
  season = "Fall",
  taxa = "Smcope",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/fallSmcopeindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/fall_smcope_model/fall_cog.rds",
  season = "Fall",
  taxa = "Smcope",
  outfile = paste0("zooplankton/outputs/Indices/fallSmcopecog", Sys.Date(), ".rds"))

## large Copepod
SOEinputs(
  infile = "zooplankton/outputs/fall_lgcope_model/Index.csv",
  season = "Fall",
  taxa = "Lgcope",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/fallLgcopeindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/fall_lgcope_model/fall_cog.rds",
  season = "Fall",
  taxa = "Lgcope",
  outfile = paste0("zooplankton/outputs/Indices/fallLgcopecog", Sys.Date(), ".rds"))

## calanus finmarchicus
SOEinputs(
  infile = "zooplankton/outputs/fall_calfin_model/Index.csv",
  season = "Fall",
  taxa = "Calfin",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/fallCalfinindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/fall_calfin_model/fall_cog.rds",
  season = "Fall",
  taxa = "Calfin",
  outfile = paste0("zooplankton/outputs/Indices/fallCalfincog", Sys.Date(), ".rds"))

  ## Spring indices ----
# Euphausiids
SOEinputs(
  infile = "zooplankton/outputs/spring_euph_model/Index.csv",
  season = "Spring",
  taxa = "Euph",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/springEuphindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/spring_euph_model/spring_cog.rds",
  season = "Spring",
  taxa = "Euph",
  outfile = paste0("zooplankton/outputs/Indices/springEuphcog", Sys.Date(), ".rds"))

# Zooplankton volume
SOEinputs(
  infile = "zooplankton/outputs/spring_zoopvol_model/Index.csv",
  season = "Spring",
  taxa = "Zoopvol",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/springZoopvolindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/spring_zoopvol_model/spring_cog.rds",
  season = "Spring",
  taxa = "Zoopvol",
  outfile = paste0("zooplankton/outputs/Indices/springZoopvolcog", Sys.Date(), ".rds"))

## small Copepod
SOEinputs(
  infile = "zooplankton/outputs/spring_smcope_model/Index.csv",
  season = "Spring",
  taxa = "Smcope",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/springSmcopeindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/spring_smcope_model/spring_cog.rds",
  season = "Spring",
  taxa = "Smcope",
  outfile = paste0("zooplankton/outputs/Indices/springSmcopecog", Sys.Date(), ".rds"))

## large Copepod
SOEinputs(
  infile = "zooplankton/outputs/spring_lgcope_model/Index.csv",
  season = "Spring",
  taxa = "Lgcope",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/springLgcopeindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/spring_lgcope_model/spring_cog.rds",
  season = "Spring",
  taxa = "Lgcope",
  outfile = paste0("zooplankton/outputs/Indices/springLgcopecog", Sys.Date(), ".rds"))

## calanus finmarchicus
SOEinputs(
  infile = "zooplankton/outputs/spring_calfin_model/Index.csv",
  season = "Spring",
  taxa = "Calfin",
  stratlook = stratlook_EPUonly,
  outfile = paste0("zooplankton/outputs/Indices/springCalfinindex_",Sys.Date(),".rds"))

SOEinputsCOG(
  infile = "zooplankton/outputs/spring_calfin_model/spring_cog.rds",
  season = "Spring",
  taxa = "Calfin",
  outfile = paste0("zooplankton/outputs/Indices/springCalfincog", Sys.Date(), ".rds"))

## Create SOE indices (Loop makes all simultaneously ----
# generated by gemini, did not test since not all models run

# Define the taxa and their file path mapping for SOE inputs
# Key: Taxa argument for SOEinputs/SOEinputsCOG (e.g., "Euph")
# Value: Model directory component (e.g., "euph")
soe_taxa_map <- list(
  Euph = "euph",
  Zoopvol = "zoopvol",
  Smcope = "smcope",
  Lgcope = "lgcope",
  Calfin = "calfin"
)

# Define the seasons to loop over for the SOE reports
soe_seasons <- list(
  Spring = "Spring",
  Fall = "Fall"
)

# Outer loop: Loop through seasons
for (SOE_season in names(soe_seasons)) {
  # Convert season to lowercase prefix (e.g., 'spring' or 'fall')
  SOE_prefix <- tolower(SOE_season)
  
  # Inner loop: Loop through the SOE taxa
  for (taxa_name in names(soe_taxa_map)) {
    # Get the path component (e.g., 'euph' for 'Euph')
    path_component <- soe_taxa_map[[taxa_name]]
    
    # Construct file paths and names dynamically
    # E.g., zooplankton/outputs/spring_euph_model
    model_dir <- paste0("zooplankton/outputs/", SOE_prefix, "_", path_component, "_model")
    
    # Index inputs
    infile_index <- paste0(model_dir, "/Index.csv")
    # E.g., zooplankton/outputs/Indices/springEuphindex_YYYY-MM-DD.rds
    outfile_index <- paste0("zooplankton/outputs/Indices/", SOE_prefix, taxa_name, "index_", Sys.Date(), ".rds")
    
    # COG inputs
    # E.g., zooplankton/outputs/spring_euph_model/spring_cog.rds
    infile_cog <- paste0(model_dir, "/", SOE_prefix, "_cog.rds")
    # E.g., zooplankton/outputs/Indices/springEuphcogYYYY-MM-DD.rds
    outfile_cog <- paste0("zooplankton/outputs/Indices/", SOE_prefix, taxa_name, "cog", Sys.Date(), ".rds")
    
    # Execute SOEinputs and SOEinputsCOG
    
    # Index
    SOEinputs(
      infile = infile_index,
      season = SOE_season, # Capitalized Season name for the function argument
      taxa = taxa_name, 
      stratlook = stratlook_EPUonly,
      outfile = outfile_index
    )
    
    # COG
    SOEinputsCOG(
      infile = infile_cog,
      season = SOE_season, # Capitalized Season name for the function argument
      taxa = taxa_name,
      outfile = outfile_cog
    )
  }
}