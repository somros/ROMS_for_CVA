rm(list = ls())
gc()

# load packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, arrow)

# Source your functions
source("functions.R")

# which run do you want to create files for? options are:
# hindcast (1990-2020)
# historical (1980-2014)
# spp585 (2015-2099)
# ssp245 (2015-2099)
# ssp126 (2015-2099)

# Do not need to rerun these, the data already exists
# create_parquet_file(run = "hindcast")
# create_parquet_file(run = "historical")
# create_parquet_file(run = "ssp126")
# create_parquet_file(run = "ssp245")
# create_parquet_file(run = "ssp585")


# # example using custom filters
# create_parquet_file(
#   run = "historical",
#   min_year = 1985,      
#   max_year = 2010,     
#   variables = c("temp"),
#   maxdepth = 300
# )
