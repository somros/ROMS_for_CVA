# Use this script to view the ROMS grid and identify where to subset it based on coordinates
# The NEP 10K grid is over a very large area, for CVA we ony need GOA
# This will be a crude subsetting just to trim off most of the excess
# later on we'll filter by depth too
# NOTE: the goal is to subset based on the rho grid (xi and eta coordinates)
# The ROMS output files do not have lat-lon, the variables are on the rho grid
# the grid.nc file is the key to use between lat-lon and xi-eta
# if you rerun this script and identify new coordinates, you'll need to rerun the bash extraction script server-side and repeat the whole process

library(sf)
library(here)
library(tidyverse)
library(tidync)
library(rnaturalearth)

# cell lat/lons from ROMS
romsfile <- here::here('data','NEP_grid_5a.nc')
roms <- tidync(romsfile)
roms_vars <- hyper_grids(roms) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables asssociated with that grid and make a reference table
  purrr::map_df(function(x){
    roms %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

latlon_rhogrd <- roms_vars %>% filter(name=="lat_rho") %>% pluck('grd')
roms_rho <- roms %>% activate(latlon_rhogrd) %>% hyper_tibble() %>%
  dplyr::select(lon_rho,lat_rho,xi_rho,eta_rho) %>% 
  st_as_sf(coords=c('lon_rho','lat_rho'),crs=4326)

# bbox
rho_bbox <- roms_rho %>%
  st_bbox() 

# limits: these you'll have to eyeball yourself if you want to change them
min_lon <- 190 # until end of 610 and then 2.5 degree buffer
max_lon <- 230 # catch the easternmost edge of Atlantis 
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

# now find if this can translate into a simple xi / eta subsetting

xi_subset <- subset_rho %>% pull(xi_rho) %>% unique()
eta_subset <- subset_rho %>% pull(eta_rho) %>% unique()

# compare to original
xi_all <- roms_rho %>% pull(xi_rho) %>% unique()
eta_all <- roms_rho %>% pull(eta_rho) %>% unique()

length(xi_subset)/length(xi_all)
length(eta_subset)/length(eta_all)

# view
rho_final <- roms_rho %>%
  filter(xi_rho %in% xi_subset, eta_rho %in% eta_subset)

p <- ggplot()+
  geom_sf(data=rho_final)+
  geom_sf(data=wc)+
  labs(x='',y='')
p

rho_final %>%
  st_bbox() 

# can reduce the size of the ROMS files to about 1/3
# better than nothing 