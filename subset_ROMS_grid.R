# ============================================================================
# ROMS Grid Subsetting Script
# ============================================================================
# Purpose: This script helps identify spatial coordinates for subsetting the 
# NEP 10K ROMS grid to the Gulf of Alaska (GOA) region. The script has two
# main sections:
#   1. Define the spatial subset in ROMS grid coordinates (xi and eta)
#   2. Create a mask to exclude cells outside NMFS management areas
#
# IMPORTANT NOTES:
# - This is a reference script - you only need to rerun it if you want to 
#   change the spatial subset boundaries
# - ROMS uses a curvilinear grid indexed by xi (columns) and eta (rows)
# - The grid.nc file maps between lat/lon coordinates and xi/eta indices
# - Subsetting must be done using xi/eta coordinates, not lat/lon
# - This approach works for variables on the rho grid (temperature, salinity)
#   but NOT for velocity components (u and v), which use staggered grids
# - If you identify new subset coordinates here, you must rerun the bash 
#   extraction script server-side and repeat the entire data processing workflow
# ============================================================================

rm(list = ls())
gc()

# Load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow)

# ============================================================================
# SECTION 1: Identify GOA Bounding Box in ROMS Grid Space
# ============================================================================
# Goal: Convert a lat/lon bounding box into xi/eta grid indices that can be
# used to subset the raw ROMS netCDF files. This reduces file size to ~1/3
# of the original by trimming excess regions outside the Gulf of Alaska.
# ============================================================================

# Read the ROMS grid file (contains coordinate mapping information)
romsfile <- here::here('data','NEP_grid_5a.nc')
roms <- tidync(romsfile)

# Extract all available grids and their associated variables from the netCDF
roms_vars <- hyper_grids(roms) %>%  # Get list of all grids in the file
  pluck("grid") %>%  # Extract grid identifiers
  purrr::map_df(function(x){  # For each grid, create a reference table
    roms %>% 
      activate(x) %>%  # Switch to this grid
      hyper_vars() %>%  # Get variables on this grid
      mutate(grd=x)  # Add grid identifier column
  })

# Identify which grid contains the rho point lat/lon coordinates
latlon_rhogrd <- roms_vars %>% 
  filter(name=="lat_rho") %>%  # Find the lat_rho variable
  pluck('grd')  # Extract its grid identifier

# Extract rho grid coordinates and convert to spatial features
roms_rho <- roms %>% 
  activate(latlon_rhogrd) %>%  # Activate the rho grid
  hyper_tibble() %>%  # Read as a tibble
  dplyr::select(lon_rho, lat_rho, xi_rho, eta_rho) %>%  # Keep only coordinate columns
  st_as_sf(coords=c('lon_rho','lat_rho'), crs=4326)  # Convert to spatial object (WGS84)

# Get overall extent of the ROMS grid
rho_bbox <- roms_rho %>% st_bbox() 

# Define the Gulf of Alaska region in lat/lon coordinates
# NOTE: Longitude uses positive values (190-230° instead of -170 to -130°)
# Adjust these limits if you need to change the spatial extent
min_lon <- 190  # Western boundary
max_lon <- 230  # Eastern boundary
min_lat <- 52   # Southern boundary
max_lat <- 60   # Northern boundary

# Create a bounding box from these limits
mask <- st_bbox(c(xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat), 
                crs = st_crs(roms_rho))

# Crop the ROMS grid to the GOA bounding box
subset_rho <- st_crop(roms_rho, mask)

# Load coastline data for visualization
wc <- ne_countries(continent="North America", returnclass = 'sf') %>% 
  filter(name %in% c('Canada','United States of America')) %>%  # Keep only US and Canada
  st_set_crs(4326)  # Set coordinate system

# Visualize the subset region
p <- ggplot()+
  geom_sf(data=subset_rho)+  # Plot ROMS grid cells in subset
  geom_sf(data=wc)+  # Overlay coastline
  labs(x='',y='')
p

# Extract the xi and eta index ranges that define this spatial box
# These are the indices you'll use to subset the raw netCDF files
xi_subset <- subset_rho %>% 
  pull(xi_rho) %>%  # Extract xi indices
  unique()  # Get unique values
range(xi_subset)  # Returns: 67 226 (your xi bounds for subsetting)

eta_subset <- subset_rho %>% 
  pull(eta_rho) %>%  # Extract eta indices
  unique()  # Get unique values
range(eta_subset)  # Returns: 237 450 (your eta bounds for subsetting)

# Calculate the proportion of the original grid retained after subsetting
xi_all <- roms_rho %>% pull(xi_rho) %>% unique()  # All xi indices
eta_all <- roms_rho %>% pull(eta_rho) %>% unique()  # All eta indices
length(xi_subset)/length(xi_all)  # Keeps 70.8% of columns
length(eta_subset)/length(eta_all)  # Keeps 33.3% of rows
# Overall, this reduces file size to approximately 1/3 of the original

# Final visualization: show the retained grid cells
rho_final <- roms_rho %>%
  filter(xi_rho %in% xi_subset, eta_rho %in% eta_subset)  # Keep only cells in subset

ggplot()+
  geom_sf(data=rho_final %>% st_shift_longitude())+  # Shift longitude for better display
  geom_sf(data=wc %>% st_shift_longitude())+  # Shift coastline to match
  labs(x='',y='')

# ============================================================================
# SECTION 2: Create a Mask Based on NMFS Management Areas
# ============================================================================
# Goal: Identify which ROMS grid cells fall outside the NMFS groundfish 
# management areas in the Gulf of Alaska. These cells will be excluded when
# reading data to reduce processing time and memory usage.
#
# CRITICAL NOTE: The xi_rho and eta_rho coordinates here are DIFFERENT from
# Section 1. After the netCDF subsetting in Section 1, the CDO commands 
# renumber the grid starting from 1. Therefore, these are called 
# xi_rho_new and eta_rho_new to distinguish them.
#
# This step is applied to reprocessed netCDF files AFTER data extraction,
# not to the raw ROMS outputs.
# ============================================================================

# Clear workspace and start fresh for this section
rm(list = ls())
gc()

# Load custom functions (assumes you have a functions.R file)
source("functions.R")

# Read the reprocessed ROMS grid (already subset from Section 1)
# Note: This grid has been renumbered, so indices start from 1
grid <- read_roms_grid() %>%
  st_as_sf(coords = c(x = "lon_rho", y = "lat_rho"), crs = 4326)  # Convert to spatial object

# Load NMFS groundfish management area shapefile
nmfs <- st_read("data/NMFS management area shapefiles/gf95_nmfs.shp")
nmfs <- nmfs %>% 
  st_transform(crs = 4326) %>%  # Ensure same coordinate system as grid
  st_shift_longitude() %>%  # Shift to 0-360° longitude for consistency
  filter(NMFS_AREA %in% c(610,620,630,640,650))  # Keep only GOA management areas

# Visualize the grid and management areas together
grid %>% 
  ggplot()+
  geom_sf()+  # Plot all grid cells
  geom_sf(data = nmfs)  # Overlay NMFS boundaries

# Identify grid cells that do NOT intersect with any NMFS management area
# These cells (e.g., in the Bering Sea) will be excluded from analysis
non_overlapping <- grid %>% 
  filter(lengths(st_intersects(., nmfs)) == 0)  # Keep cells with no intersection

# Visualize which cells will be dropped
non_overlapping %>% 
  ggplot()+
  geom_sf()

# Create a vector of cell identifiers to exclude
# Format: "xi_eta" (e.g., "15_23" represents the cell at xi=15, eta=23)
to_drop <- non_overlapping %>%
  mutate(idx_drop = paste(xi_rho_new, eta_rho_new, sep = "_")) %>%  # Create unique ID
  pull(idx_drop) %>%  # Extract the ID column
  unique()  # Ensure no duplicates

# Write the exclusion list to a text file
# This file will be read when loading data to skip these cells
# Commented out to prevent accidental overwriting - uncomment to update
# cat(file = "idx_to_drop.txt", to_drop)