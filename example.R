library(tidyverse)
library(tidync)
library(ncdf4)
library(lubridate)
library(here)
library(sf)
library(rnaturalearth)

# Source your functions
source("functions.R")

hind_file <- "data/combined_annual_data.RDS"

if(file.exists(hind_file)){
  
  all_data <- readRDS(hind_file)
  
} else {
  
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
      roms_grid = grid
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
  
  # Summary
  message(paste("Total rows:", nrow(all_data)))
  message(paste("Date range:", min(all_data$date), "to", max(all_data$date)))
  message(paste("Variables:", paste(unique(all_data$variable), collapse = ", ")))
  
  # 3. Save the combined dataset
  saveRDS(all_data, hind_file, compress = "xz")
  message("Data saved to data/combined_annual_data.RDS")
  
}

# 4. Create spatial maps
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62)

plot_spatial_map(all_data, 
                 variable = "temp", 
                 year = 2010, 
                 month = 6, 
                 coastline = coast)

# 4. Create time series plots
plot_time_series(all_data, 
                 variable = "temp",
                 start_year = 1995,
                 end_year = 2000)

# # make gif for talk
# library(magick)
# 
# # Filter data for 2010 and calculate proper range
# temp_range <- all_data %>% 
#   mutate(yr = year(date)) %>%
#   filter(yr == 2010, variable == "temp") %>%
#   pull(value) %>%
#   range()
# 
# print(temp_range)
# 
# temp_dir <- tempdir()
# 
# # Generate plots with fixed COLOR scale (not fill)
# for(i in 1:12) {
#   p <- plot_spatial_map(all_data, 
#                         variable = "temp", 
#                         year = 2010, 
#                         month = i, 
#                         coastline = coast)
#   
#   # Override with fixed color scale limits (note: scale_COLOR not scale_FILL)
#   p <- p + scale_color_viridis_c(option = "plasma", limits = temp_range)
#   
#   ggsave(filename = file.path(temp_dir, paste0("month_", sprintf("%02d", i), ".png")),
#          plot = p,
#          width = 8, 
#          height = 6, 
#          dpi = 150)
# }
# 
# # Create GIF
# img_files <- list.files(temp_dir, pattern = "month_.*\\.png", full.names = TRUE)
# img_files <- img_files[order(img_files)]
# 
# img_list <- lapply(img_files, image_read)
# img_joined <- image_join(img_list)
# img_animated <- image_animate(img_joined, fps = 2)
# 
# image_write(img_animated, path = "temperature_2010.gif")
# file.remove(img_files)
