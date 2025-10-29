# ============================================================================
# ROMS Data Exploration Script
# ============================================================================
# Purpose: This script demonstrates how to efficiently load and visualize
# processed ROMS ocean model data. It uses the Apache Parquet file format
# for fast, selective data loading without reading entire datasets into memory.
#
# What this script does:
#   1. Opens a processed ROMS dataset (stored as Parquet)
#   2. Filters to specific variables, depths, and time periods
#   3. Creates spatial maps showing conditions at a specific time
#   4. Generates time series plots showing trends over multiple years
#
# IMPORTANT NOTES:
# - The Parquet files are pre-filtered to the Gulf of Alaska region and 
#   depths ≤1000m to minimize file size
# - If you need broader spatial coverage or deeper depths, you must reprocess
#   the original netCDF files using the process_annual_file() function
# - Projection runs (ssp126, ssp585) are NOT bias-corrected yet
# - This script assumes you have already processed raw ROMS outputs into
#   Parquet format using the main processing workflow
# ============================================================================

# Clear workspace and free up memory
rm(list = ls())
gc()

# Load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow)
pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data

# Load custom analysis functions
source("functions.R")

# ============================================================================
# SECTION 1: Load and Filter Data
# ============================================================================
# The Arrow package allows "lazy evaluation" - you can filter the data before
# loading it into RAM, which is much faster than loading everything first.
# ============================================================================

# Specify which model run to analyze
# Options  include: "hindcast", "historical", "ssp126", "ssp245", "ssp585"
# WARNING: Projection runs (ssp*) are not yet bias-corrected
this_run <- "ssp585"

# Construct the file path to the processed Parquet file
this_file <- paste0("data/processed/", this_run, "_annual_data.parquet")

# Open the dataset WITHOUT loading it into memory
# This creates a connection that allows you to query the data structure
all_data <- open_dataset(this_file)  # Uses Apache Arrow for efficient access

# Filter to specific data subset before loading
# This is much more memory-efficient than loading everything first
# You don't have to apply filters but it will occupy less RAM if you do
# If you want to read many data stream at once, you should definitely filter
data_subset <- all_data %>%
  filter(layer == "surface") %>%  # Select only surface layer (other options: "bottom", "30m", etc.)
  filter(variable == "temp") %>%  # Choose variable of interest (temp, oxygen, salinity, etc.)
  filter(date > as.Date("2030-01-01")) %>%  # Restrict to dates after Jan 1, 2030
  collect()  # NOW load the filtered data into RAM

# ============================================================================
# SECTION 2: Create Spatial Maps
# ============================================================================
# Visualize ocean conditions at a specific point in time across the
# Gulf of Alaska region.
# ============================================================================

# Load coastline data for map visualization
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude to match ROMS data

# Generate a spatial map for a specific month and year
# This uses a custom plotting function defined in functions.R
plot_spatial_map(data = data_subset,  # The filtered dataset
                 variable = "temp",  # Variable to plot (must match filter above)
                 year = 2039,  # Year to visualize
                 month = 8,  # Month to visualize (1-12)
                 coastline = coast,  # Coastline overlay for context
                 psize = 1)  # Point size for grid cells

# ============================================================================
# SECTION 3: Create Time Series Plots
# ============================================================================
# Visualize how conditions change over time, aggregated across the study area
# or for specific regions/depth layers.
# ============================================================================

# Generate a time series showing trends from 2025 to 2090
# This uses a custom plotting function defined in functions.R
plot_time_series(data_subset,  # The filtered dataset
                 variable = "temp",  # Variable to plot (must match filter above)
                 start_year = 2025,  # Beginning of time series
                 end_year = 2090)  # End of time series

# ============================================================================
# Next Steps and Customization
# ============================================================================
# To explore different scenarios:
# - Change 'this_run' to compare historical vs. future projections
# - Modify the filter() calls to examine different variables, depths, or time periods
# - Adjust plot parameters (year, month, date ranges) to focus on specific events
# - Use the plot functions with different subsets to compare conditions
#
# To access data outside the pre-processed bounds:
# - Use process_annual_file() to reprocess netCDF files with custom spatial
#   or depth filters
# - See the main processing scripts for examples of how to do this
# ============================================================================