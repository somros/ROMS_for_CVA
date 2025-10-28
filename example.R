library(tidyverse)
library(tidync)
library(ncdf4)
library(lubridate)
library(here)
library(sf)
library(rnaturalearth)
library(arrow)

rm(list = ls())
gc()

# Source your functions
source("functions.R")

this_file <- "data/hindcast_annual_data.parquet"

# open without loading:
all_data <- open_dataset(this_file)

# filter to the variables and layer of interest
data_subset <- all_data %>%
  filter(layer == "surface") %>%
  filter(variable == "temp") %>%
  collect()  # Only now loads to RAM

# 4. Create spatial maps
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62)

plot_spatial_map(data = data_subset, 
                 variable = "temp", 
                 year = 2014, 
                 month = 8, 
                 coastline = coast,
                 psize = 1)

# 4. Create time series plots
plot_time_series(all_data, 
                 variable = "temp",
                 start_year = 2010,
                 end_year = 2020)

