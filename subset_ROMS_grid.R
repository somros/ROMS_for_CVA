# Use this script to view the ROMS grid and identify where to subset it based on coordinates
# Note: you only need to rerun this script if you wan tto change the spatial subset of the ROMS data

# The NEP 10K grid is over a very large area, for CVA we ony need GOA
# This will be a crude subsetting just to trim off most of the excess
# later on we'll filter by depth and NMFS area too
# NOTE: the goal is to subset based on the rho grid (xi and eta coordinates)
# The ROMS output files do not have lat-lon, the variables are on the rho grid
# the grid.nc file is the key to use between lat-lon and xi-eta
# if you rerun this script and identify new coordinates, you'll need to rerun the bash extraction script server-side and repeat the whole process
# note that this won't work for velocities (u and v) as those are not on the rho grid

# load packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow)
pacman::p_load_gh("ropensci/rnaturalearthhires")

# cell lat/lons from ROMS
romsfile <- here::here('data','NEP_grid_5a.nc')
roms <- tidync(romsfile)
roms_vars <- hyper_grids(roms) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    roms %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

latlon_rhogrd <- roms_vars %>% filter(name=="lat_rho") %>% pluck('grd')
roms_rho <- roms %>% activate(latlon_rhogrd) %>% hyper_tibble() %>%
  dplyr::select(lon_rho,lat_rho,xi_rho,eta_rho) %>% 
  st_as_sf(coords=c('lon_rho','lat_rho'),crs=4326)

# bbox
rho_bbox <- roms_rho %>% st_bbox() 

# limits: these you'll have to eyeball yourself if you want to change them
# note that these use positive coordinates for longitude
min_lon <- 190
max_lon <- 230 
min_lat <- 52
max_lat <- 60
mask <- st_bbox(c(xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat), 
                crs = st_crs(roms_rho))

subset_rho <- st_crop(roms_rho, mask)

# continent
wc <- ne_countries(continent="North America",returnclass = 'sf') %>% 
  filter(name %in% c('Canada','United States of America')) %>%
  st_set_crs(4326)

# plot
p <- ggplot()+
  geom_sf(data=subset_rho)+
  geom_sf(data=wc)+
  labs(x='',y='')
p

# identify the ranges of xi and eta that define this spatial box
# those indices will be what you need to use to subset the netcdf files
# you cannot subset the netcdf files by lat and lon
xi_subset <- subset_rho %>% pull(xi_rho) %>% unique()
range(xi_subset) # 67 226 are your xi bounds
eta_subset <- subset_rho %>% pull(eta_rho) %>% unique()
range(eta_subset) # 237 450 are your eta bounds

# compare to original
xi_all <- roms_rho %>% pull(xi_rho) %>% unique()
eta_all <- roms_rho %>% pull(eta_rho) %>% unique()
length(xi_subset)/length(xi_all) # 0.7079646
length(eta_subset)/length(eta_all) # 0.3333333

# view
rho_final <- roms_rho %>%
  filter(xi_rho %in% xi_subset, eta_rho %in% eta_subset)

ggplot()+
  geom_sf(data=rho_final %>% st_shift_longitude())+
  geom_sf(data=wc %>% st_shift_longitude())+
  labs(x='',y='')

# this should reduce the size of the ROMS files to about 1/3 of their original size

# Create a mask based on the NMFS areas -----------------------------------
# Very important note before you begin
# xi_rho and eta_rho have a different meaning here than above
# this step is applied to the reprocessed netcdf files
# the grid has changed in those compared to the raw ROMS outputs

# identify which ROMS cells are in the EBS and are therefore to be dropped at the stage of reading the nc files
# doing this now will allow us to ditch a good amount of cells
# work on the grid object, match it with the NMFS areas, write out the voctor of xi and eta to ditch
# do not need to run this more than once

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
