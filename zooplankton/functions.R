#' @title Format VAST Index Output for SOE Reporting
#'
#' @description
#' Reads a CSV file containing VAST index estimates and standard errors (SE),
#' processes the data using \code{dplyr} and \code{tidyr} to conform to a
#' specific long format, and saves the resulting data frame as an RDS file.
#' The processing includes filtering by specified EPUs (Ecological Production Units)
#' and creating a combined variable name from the season, taxa, and variable type.
#' *This function is not used directly, but is wrapped in `create_soe_indices`.*
#'
#' @param infile A character string specifying the path to the input CSV file
#'   containing VAST index output (e.g., "Index.csv"). This file is expected
#'   to have columns like 'Time', 'Region', 'Estimate', and 'Std.Error.for.Estimate'.
#' @param season A character string specifying the season (e.g., "Spring").
#'   Used to create the final \code{Var} column names.
#' @param stratlook A data frame used for joining to map stratification units
#'   to Ecological Production Units (EPU). It is expected to contain columns
#'   that link to the input data's stratification and a 'Region' column for the EPU.
#' @param taxa A character string specifying the taxa/species (e.g., "Calanus finmarchicus").
#'   Used to create the final \code{Var} column names.
#' @param outfile A character string specifying the path and name for the output
#'   RDS file (e.g., "Spring_Calanus_index_2023-10-27.rds").
#'
#' @return The function saves an RDS file to the path specified by \code{outfile}.
#'   The saved object is a data frame with columns: \code{Time}, \code{Var}
#'   (e.g., "Spring Calanus finmarchicus Abundance Index Estimate"), \code{Value},
#'   \code{EPU}, and \code{Units}. Returns \code{NULL} invisibly.
#'
#' @importFrom dplyr %>% left_join select filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_c
#' @importFrom here here

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

  zoopindex$Var <- stringr::str_c(
    stringr::str_c(season, taxa, sep = " ") |>
      stringr::str_to_title(),
    zoopindex$Var,
    sep = " "
  )

  saveRDS(zoopindex, outfile)
}


#' @title Extract Center of Gravity from a VAST Model Fit
#'
#' @description
#' Calculates the Center of Gravity (COG) from a fitted VAST model object
#' using the \code{FishStatsUtils::plot_range_index} function. It temporarily
#' creates a directory for diagnostic plots (which is immediately deleted),
#' and returns the resulting COG table.
#' *This function is not used directly, but is wrapped in `run_vast_model` to create center of gravity VAST output.*
#'
#' @param model_fit A list object representing a fitted VAST model, typically
#'   the output of \code{VAST::fit_model()}. This object must contain the
#'   necessary components like \code{parameter_estimates}, \code{Report},
#'   \code{data_list}, and \code{year_labels}.
#'
#' @return A list object containing the COG table and related information,
#'   which is the direct output of \code{FishStatsUtils::plot_range_index}.
#'   The COG table is typically a data frame containing columns for Year,
#'   COG estimate, and SE.
#'
#' @importFrom here here
#' @importFrom FishStatsUtils plot_range_index

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

#' @title Format VAST Center of Gravity (COG) Output for SOE Reporting
#'
#' @description
#' Reads an RDS file containing the raw VAST Center of Gravity (COG) output,
#' processes it using \code{dplyr} and \code{tidyr} to conform to a
#' specific long format, and saves the resulting data frame as an RDS file.
#' The processing includes pivoting the COG estimate and SE to a long format,
#' combining the Eastward/Northward direction with the variable name, and
#' assigning common units and EPU.
#' *This function is not used directly, but is wrapped in `create_soe_indices`.*
#'
#' @param infile A character string specifying the path to the input RDS file
#'   containing the raw COG output (a list with a \code{COG_Table} element).
#' @param season A character string specifying the season (e.g., "spring").
#'   Used to create the final \code{Var} column names.
#' @param taxa A character string specifying the taxa/species (e.g., "calanus").
#'   Used to create the final \code{Var} column names.
#' @param outfile A character string specifying the path and name for the output
#'   RDS file (e.g., "Spring_Calanus_cog_2023-10-27.rds").
#'
#' @return The function saves an RDS file to the path specified by \code{outfile}.
#'   The saved object is a data frame with columns: \code{Time}, \code{Var}
#'   (e.g., "Spring Calanus Eastward Center of Gravity"), \code{Value},
#'   \code{EPU}, and \code{Units}. Returns \code{NULL} invisibly.
#'
#' @importFrom dplyr %>% mutate select
#' @importFrom tidyr pivot_longer unite
#' @importFrom stringr str_c str_to_title

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

#' @title Create State of the Ecosystem (SOE) Index and COG Files from VAST Output Directory
#'
#' @description
#' Takes a directory path containing VAST model results (e.g., Index.csv and a raw COG RDS file)
#' and processes both the abundance index and the Center of Gravity (COG) into the
#' standardized SOE reporting format using \code{SOEinputs} and \code{SOEinputsCOG}.
#' It extracts the season and taxa from the directory name.
#'
#' @param dir A character string specifying the directory path where the VAST model
#'   results are saved. The directory name is expected to be in the format
#'   "SEASON_TAXA" (e.g., "Spring_Calanus").
#' @param strat A data frame used for joining in the \code{SOEinputs} function to map
#'   stratification units to Ecological Production Units (EPU). Defaults to
#'   \code{stratlook_EPUonly}.
#'
#' @return The function does not return a value but saves two RDS files:
#'   one for the COG data and one for the abundance index data, to the
#'   \code{zooplankton/outputs/Indices} directory (using \code{here::here}).
#'
#' @importFrom here here
#' @importFrom base strsplit basename Sys.Date paste0
#'
#' @seealso \code{\link{SOEinputs}}, \code{\link{SOEinputsCOG}}

create_soe_indices <- function(
  dir, # directory names where model results are saved
  strat
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
#' @title Run a VAST Model and Process Outputs
#'
#' @description
#' Executes the main VAST model fitting routine using \code{VAST::fit_model()}.
#' It handles the setup of the output directory, fits the model with or without
#' a Day-of-Year (DOY) covariate based on a flag, saves the fitted object,
#' generates standard diagnostic plots, and extracts the Center of Gravity (COG).
#'
#' @param data A data frame containing the input data for VAST, including columns
#'   like 'Lat', 'Lon', 'Year', 'Catch_g', and 'Vessel'. Must also contain
#'   'Dayofyear' if \code{doy = TRUE}.
#' @param out_dir A character string specifying the path to the directory where
#'   all model outputs (fit object, plots, COG RDS) will be saved.
#' @param season A character string specifying the season (e.g., "Spring").
#'   Used for naming the saved COG file.
#' @param vast_settings A list of settings for the VAST model, typically created
#'   by \code{VAST::make_settings()}.
#' @param doy A logical value (\code{TRUE} or \code{FALSE}) indicating whether
#'   the 'Dayofyear' column should be included as a covariate (\code{Q_ik}) in the VAST model.
#'   Defaults to \code{FALSE}.
#'
#' @return The function returns the fitted VAST model object (a list). It also has
#'   multiple side effects: creates \code{out_dir} if it doesn't exist, saves the
#'   fitted model object (\code{fit.rds}), saves plots, and saves the COG results
#'   (\code{SEASON_cog.rds}).
#'
#' @importFrom VAST fit_model
#' @importFrom units as_units
#' @importFrom here here
#'
#' @seealso \code{\link{extract_cog}}

run_vast_model <- function(data, out_dir, season, vast_settings, doy = FALSE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  if (doy == TRUE) {
    fit <- fit_model(
      settings = vast_settings,
      Lat_i = data$Lat,
      Lon_i = data$Lon,
      t_i = data$Year,
      b_i = as_units(data[, 'Catch_g'], 'g'),
      a_i = rep(1, nrow(data)),
      v_i = data$Vessel,
      Q_ik = as.matrix(data[, c("Dayofyear")]),
      working_dir = paste0(out_dir, "/")
    )
  } else {
    fit <- fit_model(
      settings = vast_settings,
      Lat_i = data$Lat,
      Lon_i = data$Lon,
      t_i = data$Year,
      b_i = as_units(data[, 'Catch_g'], 'g'),
      a_i = rep(1, nrow(data)),
      v_i = data$Vessel,
      working_dir = paste0(out_dir, "/")
    )
  }

  saveRDS(fit, file = paste0(out_dir, "/fit.rds"))

  # Plot results
  plot(fit, working_dir = paste0(out_dir, "/"))

  # extract center of gravity
  cog <- extract_cog(fit)
  saveRDS(cog, here::here(out_dir, paste0(season, "_cog.rds")))
}

# run_vast_model(data = euph_stn_fall,
#                out_dir = here::here(sprintf("zooplankton/outputs/fall_euph_model")),
#                season = "fall",
#                vast_settings = settings)
