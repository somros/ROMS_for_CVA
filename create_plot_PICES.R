# ============================================================================
# PICES Conference Figure Generation Script
# ============================================================================
# Purpose: This script creates a spatial map of sea surface
# temperature in the Gulf of Alaska for presentation at PICES 
#
# What this script does:
#   1. Reads the ROMS grid and spatial mask
#   2. Processes multiple annual netCDF files for a specific year
#   3. Combines data from all files
#   4. Generates a high-resolution spatial map
#   5. Exports the figure as a PNG file
#
# Key Features:
# - Processes raw netCDF files rather than pre-processed Parquet files
# - Allows custom depth filtering (here set to include all depths)
# - Focuses on a single year for a snapshot visualization
#
# Customize This Script:
# - Change which run you want to plot
# - Change min_year/max_year to visualize different time periods
# - Modify variables parameter to plot other ocean properties
# - Adjust month parameter in plot_spatial_map() for seasonal comparisons
# - Change layer filter to show bottom temperatures or other depths
# ============================================================================

# Load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow)
pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data

# Load custom analysis functions
source("functions.R")

# ============================================================================
# SECTION 1: Load Grid and Spatial Mask
# ============================================================================

# Read the ROMS grid structure
# This contains the spatial coordinates (lat/lon) for each grid cell
grid <- read_roms_grid()

# Load the spatial mask identifying cells to exclude
# This mask removes cells outside NMFS management areas (e.g., Bering Sea)
goa_mask <- scan("idx_to_drop.txt", "character", sep = " ")

# ============================================================================
# SECTION 2: Process Annual NetCDF Files
# ============================================================================
# The data are stored as multiple annual netCDF files with monthly time steps. This section
# reads all files, extracts the relevant data, and combines them into a
# single dataset. Processing from raw netCDF allows custom filtering not
# available in the pre-processed Parquet files.
# ============================================================================

# Which run do you want?
this_run <- "ssp585"

# Find all annual netCDF files in the hindcast directory
# Pattern matches files like "annual_2010.nc", "annual_2011.nc", etc.
nc_files <- list.files(paste0("data/annual_files/", this_run), 
                       pattern = "annual_.*\\.nc$", 
                       full.names = FALSE)

# Report how many files were found
message(paste("Found", length(nc_files), "files to process"))

# Initialize an empty list to store data from each file
all_data_list <- list()

# Loop through each netCDF file and process it
for (i in seq_along(nc_files)) {
  file <- nc_files[i]
  
  # Process the current file using custom function
  # This extracts variables, applies spatial mask, and filters by depth/time
  data <- process_annual_file(
    ncfile = file,  # Current file name
    data_dir = paste0("data/annual_files/", this_run),  # Directory containing files
    roms_grid = grid,  # Grid structure for coordinate mapping
    variables = "temp",  # Extract only temperature (faster than all variables)
    min_year = 2039,  # Start of time period to extract
    max_year = 2039,  # End of time period
    maxdepth = 9999  # Include all depths (no depth filtering)
  )
  
  # Add processed data to the list
  all_data_list[[i]] <- data
  
  # Print progress message to track processing
  message(paste("Completed", i, "of", length(nc_files)))
  
  # Periodically free up memory to prevent crashes with large datasets
  # Runs garbage collection every 5 files
  if (i %% 5 == 0) {
    gc()
  }
}

# Combine all processed data into a single dataframe
all_data <- bind_rows(all_data_list)

# Sort the data for easier subsetting and analysis
all_data <- all_data %>%
  arrange(date, variable, layer)

# ============================================================================
# SECTION 3: Create and Export Spatial Map
# ============================================================================
# Generate a map showing sea surface temperature across the Gulf of Alaska.
# ============================================================================

# Load high-resolution coastline data for the map
coast <- ne_coastline(scale = "large", returnclass = "sf") %>%  # Use detailed coastline
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to Gulf of Alaska
  st_shift_longitude()  # Convert to 0-360° longitude to match ROMS coordinates

# Generate the spatial map
# This uses a custom plotting function that handles ROMS grid data
p <- plot_spatial_map(
  all_data %>% filter(layer == "surface"),  # Show only surface layer (SST)
  variable = "temp",  # Plot temperature
  year = 2039,  # Year to visualize
  month = 8,  # August (summer conditions)
  coastline = coast,  # Add coastline overlay
  psize = 0.75  # Point size for grid cells - play with this depending on the size you want the image
)
p

# Export the figure as a high-resolution PNG
ggsave("sst.png",  # Output filename
       p,  # Plot object to save
       dpi = 600,  # Resolution
       width = 8,  # Width in inches
       height = 8)  # Height in inches (square format)

# ============================================================================
# Output Location and Next Steps
# ============================================================================
# The figure has been saved to: sst.png (in your working directory)
#
# To modify this figure:
# - Change the year/month to show different time periods
# - Adjust layer filter to "bottom" for seafloor temperatures
# - Modify variables parameter to plot salinity, oxygen, etc.
# - Change psize for different visual detail levels
# - Adjust figure dimensions (width/height) for different aspect ratios
#
# For multiple figures:
# - Wrap the plotting code in a loop over years or months
# - Use paste0() to create unique filenames for each figure
# ============================================================================
