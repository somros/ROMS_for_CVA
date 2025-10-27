# script to make a figure for PICES
library(tidyverse)
library(tidync)
library(ncdf4)
library(lubridate)
library(here)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(stars)

# Source your functions
source("functions.R")

# 1. Read ROMS grid
grid <- read_roms_grid()

# 2. Read all annual files and bind them together
nc_files <- list.files("data/annual_files/hindcast/", pattern = "annual_.*\\.nc$", full.names = FALSE)

message(paste("Found", length(nc_files), "files to process"))

# Initialize empty list to store results
all_data_list <- list()

# Loop through files
for (i in seq_along(nc_files)) {
  file <- nc_files[i]
  
  # Process file
  data <- process_annual_file(
    ncfile = file,
    data_dir = "data/annual_files/hindcast",
    roms_grid = grid,
    variables = "temp",
    min_year = 2015,
    max_year = 2015,
    maxdepth = 9999
  )
  
  # Add to list
  all_data_list[[i]] <- data
  
  # Progress message
  message(paste("Completed", i, "of", length(nc_files)))
  
  # Optional: garbage collection every few files to free up memory
  if (i %% 5 == 0) {
    gc()
  }
}

# Bind all data together
all_data <- bind_rows(all_data_list)

# Sort by date
all_data <- all_data %>%
  arrange(date, variable, layer)

# 4. Create spatial maps
coast <- ne_coastline(scale = "large", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62)

p <- plot_spatial_map(all_data %>% filter(layer == "surface"), 
                 variable = "temp", 
                 year = 2015, 
                 month = 8, 
                 coastline = coast,
                 psize = 0.8) +
  coord_sf(xlim = c(-170,-133), ylim = c(50, 62))
p

ggsave("sst.png", p, dpi = 300, width = 8, height = 8)

