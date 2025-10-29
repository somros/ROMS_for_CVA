# this script allows you to read in a ROMS run and make some plots / explore the data
# this script leverages the parquet format to load only the data you need into the session
# NOTE: the full data in parquet format are cropped to the GOA and to 1000 m, to minimize their size
# If you want plots for broader areas you need to read the nc files again, which you can do with the function process_annual_file()

rm(list = ls())
gc()

# load packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow)
pacman::p_load_gh("ropensci/rnaturalearthhires")

# what run do you want to look at?
this_run <- "hindcast"

# Source your functions
source("functions.R")

this_file <- paste0("data/", this_run, "_annual_data.parquet")

# open without loading
all_data <- open_dataset(this_file) # leverage arrow

# filter to the variables and layer of interest
data_subset <- all_data %>%
  filter(layer == "surface") %>%
  filter(variable == "temp") %>%
  collect()  # Only now loads to RAM

# Create spatial maps
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%
  st_shift_longitude()

plot_spatial_map(data = data_subset, 
                 variable = "temp", 
                 year = 2015, 
                 month = 8, 
                 coastline = coast,
                 psize = 1)

# Create time series plots
plot_time_series(data_subset, 
                 variable = "temp",
                 start_year = 2025,
                 end_year = 2090)
