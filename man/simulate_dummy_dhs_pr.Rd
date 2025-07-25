% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_dummy_data.R
\name{simulate_dummy_dhs_pr}
\alias{simulate_dummy_dhs_pr}
\title{Simulate Processed DHS Survey Personal Records Data}
\usage{
simulate_dummy_dhs_pr(
  country = "Gambia",
  country_code_iso3 = "GMB",
  country_code_dhs = "GM",
  year_of_survey = 2024,
  total_population = 266,
  urban_proportion = 0.602,
  lon_range = c(-16.802, -13.849),
  lat_range = c(13.149, 13.801),
  mean_web_x = -1764351,
  mean_web_y = 1510868,
  log_scale_mean = 2.82,
  log_scale_sd = 0.2,
  log_shape_mean = 0.331,
  log_shape_sd = 0.1,
  b1_mean = 0.0142,
  b1_sd = 0.002,
  c_mean = -0.00997,
  c_sd = 0.001,
  b2_mean = 0.00997,
  b2_sd = 0.002,
  nsampled_range = c(180, 220),
  output_path = here::here("01_data", "1a_survey_data", "processed",
    "dhs_pr_records_combined.rds"),
  seed = 123
)
}
\arguments{
\item{country}{Character. Country name. Default is "Gambia".}

\item{country_code_iso3}{Character. ISO3 country code. Default is "GMB".}

\item{country_code_dhs}{Character. DHS country code. Default is "GM".}

\item{year_of_survey}{Integer. Year of the survey. Default is 2024.}

\item{total_population}{Integer. Total number of individuals to simulate.
Default is 266.}

\item{urban_proportion}{Numeric. Proportion of the population in urban areas.
Default is 0.602.}

\item{lon_range}{Numeric vector of length 2. Longitude bounds for simulation.
Default is c(-16.802, -13.849).}

\item{lat_range}{Numeric vector of length 2. Latitude bounds for simulation.
Default is c(13.149, 13.801).}

\item{mean_web_x}{Numeric. Mean Web Mercator X coordinate. Default is
-1764351.}

\item{mean_web_y}{Numeric. Mean Web Mercator Y coordinate. Default is
1510868.}

\item{log_scale_mean}{Numeric. Mean of log-scale parameter. Default is 2.82.}

\item{log_scale_sd}{Numeric. SD of log-scale parameter. Default is 0.2.}

\item{log_shape_mean}{Numeric. Mean of log-shape parameter. Default is 0.331.}

\item{log_shape_sd}{Numeric. SD of log-shape parameter. Default is 0.1.}

\item{b1_mean}{Numeric. Mean of b1. Default is 0.0142.}

\item{b1_sd}{Numeric. SD of b1. Default is 0.002.}

\item{c_mean}{Numeric. Mean of c. Default is -0.00997.}

\item{c_sd}{Numeric. SD of c. Default is 0.001.}

\item{b2_mean}{Numeric. Mean of b2. Default is 0.00997.}

\item{b2_sd}{Numeric. SD of b2. Default is 0.002.}

\item{nsampled_range}{Integer vector of length 2. Range of sample sizes per
cluster. Default is c(180, 220).}

\item{output_path}{Character. Path to save the resulting RDS file.}

\item{seed}{Integer. Random seed for reproducibility. Default is 123.}
}
\value{
Saves a list containing a tibble (`age_param_data`) to the specified
   path.
}
\description{
This function simulates a synthetic DHS survey dataset, generating age
distribution parameters, coordinates, and metadata. It mimics typical
DHS/MICS output and saves the result to disk. The metadata includes:
- `country`: Name of the country ("Gambia" by default)
- `country_code_iso3`: ISO3 country code ("GMB" by default)
- `country_code_dhs`: DHS country code ("GM" by default)
- `year_of_survey`: Year of the simulated survey (2024 by default)
}
