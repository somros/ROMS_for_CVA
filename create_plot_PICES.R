# script to make a figure for PICES
# load packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow)
pacman::p_load_gh("ropensci/rnaturalearthhires")

# Source functions
source("functions.R")

# Read ROMS grid
grid <- read_roms_grid()
goa_mask <- scan("idx_to_drop.txt", "character", sep = " ")

# Read all annual files and bind them together
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

# Create spatial maps
coast <- ne_coastline(scale = "large", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%
  st_shift_longitude()

p <- plot_spatial_map(all_data %>% filter(layer == "surface"), 
                 variable = "temp", 
                 year = 2015, 
                 month = 8, 
                 coastline = coast,
                 psize = 0.8)
p

ggsave("sst.png", p, dpi = 300, width = 8, height = 8)

