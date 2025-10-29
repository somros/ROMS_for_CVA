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

# TODO: wrap this as a function where the argument is the run and years

# Source your functions
source("functions.R")

# which run do you want to create files for? options are:
# hindcast
# historical
# spp585
# ssp245
# ssp126

run <- "hindcast"

outfile <- paste0("data/", run, "_annual_data.parquet")

# Read ROMS grid
grid <- read_roms_grid()

# Read in mask for which ROMS point to drop
goa_mask <- scan("idx_to_drop.txt", "character", sep = " ")

# Read all annual files and bind them together
nc_files <- list.files(paste0("data/annual_files/", run), pattern = "annual_.*\\.nc$", full.names = FALSE)

message(paste("Found", length(nc_files), "files to process"))

# Initialize empty list to store results
all_data_list <- list()

# Loop through files
for (i in seq_along(nc_files)) {
  file <- nc_files[i]
  
  # Process file
  data <- process_annual_file(
    ncfile = file,
    data_dir = paste0("data/annual_files/", run),
    roms_grid = grid,
    min_year = 1990,
    max_year = 2020
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

rm(all_data_list)
gc()

# Sort by date
all_data <- all_data %>%
  arrange(date, variable, layer) %>%
  mutate(run = factor(run))

# Summary
message(paste("Total rows:", nrow(all_data)))
message(paste("Date range:", min(all_data$date), "to", max(all_data$date)))
message(paste("Variables:", paste(unique(all_data$variable), collapse = ", ")))

# 3. Save the combined dataset
write_parquet(all_data, outfile)
message("Data saved to parquet")
