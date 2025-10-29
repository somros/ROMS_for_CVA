# identify which ROMS cells are in the EBS and are therefore to be dropped at the stage of reading the nc files
# doing this now will allow us to ditch a good amount of cells
# work on the grid object, match it with the NMFS areas, write out the voctor of xi and eta to ditch
# do not need to run this more than once

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

# 1. Read ROMS grid
grid <- read_roms_grid() %>%
  st_as_sf(coords = c(x = "lon_rho", y = "lat_rho"), crs = 4326)

# read in NMFS areas, NOT trimmed to depth
nmfs <- st_read("data/NMFS management area shapefiles/gf95_nmfs.shp")
nmfs <- nmfs %>% 
  st_transform(crs = 4326) %>%
  filter(NMFS_AREA %in% c(610,620,630,640,650)) # subset to GOA areas

# view
grid %>% 
  ggplot()+
  geom_sf()+
  geom_sf(data = nmfs)

# write out a blacklist of xi and eta pairs
# Find points that DO NOT intersect with any nmfs polygon
non_overlapping <- grid %>%
  filter(lengths(st_intersects(., nmfs)) == 0)

# view
non_overlapping %>% 
  ggplot()+
  geom_sf()

to_drop <- non_overlapping %>%
  mutate(idx_drop = paste(xi_rho, eta_rho, sep = "_")) %>%
  pull() %>%
  unique()

# write out
cat(file = "idx_to_drop.txt", to_drop)
