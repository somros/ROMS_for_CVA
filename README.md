# ROMS Data Processing and Visualization

This repository contains code for processing and visualizing Regional Ocean Modeling System (ROMS) ocean model output data for the Gulf of Alaska region.

*NOTE*: Most GOACLIM models use the approach outlines [here](https://github.com/GOA-CLIM/ROMS_to_Index/tree/main), which uses ROMS monthly output to produce area- and depth-aggregated indices at the level of NMFS areas or aroader GOA. *You should only be looking at the repo you are currently in if you need spatially-explicit ROMS indices*.

This code only pulls the surface and bottom layers from ROMS nc files. Depth-integrated indices are not currently being produced.

## Workflow

### 1. Server-side Pre-processing (on loon)

Run `pre-processing.sh` on the loon server where NetCDF files are hosted:

```bash
./pre-processing.sh
```

This script:
- Extracts surface and bottom layers for 10 variables (temp, salt, PhS, PhL, MZS, MZL, Cop, NCa, Eup, Det)
- Spatially subsets data to the Gulf of Alaska region (xi_rho: 67-226, eta_rho: 237-450)
- Combines monthly files into annual files (1990-2020)
- Outputs processed files to `./annual_files/`

**Requirements:** CDO (Climate Data Operators)

### 2. Local Analysis (R)

Transfer the processed annual files to your local machine under `data/annual_files/hindcast/`, then:

```r
source("example.R")
```

This script:
- Reads the ROMS grid and processed annual files
- Combines all years into a single dataset
- Creates spatial maps and time series visualizations

## Files

- `functions.R` - Core functions for reading, processing, and plotting ROMS data
- `example.R` - Example workflow demonstrating data loading and visualization
- `pre-processing.sh` - Server-side bash script for data extraction and subsetting
- `subset_ROMS_grid.R` - Utility script to identify grid coordinates (for setup only, not needed for routine use)

## R Dependencies

```r
tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth
```

## Data Structure

Processed files contain:
- **Variables:** temp, salt, PhS, PhL, MZS, MZL, Cop, NCa, Eup, Det
- **Vertical layers:** Surface and bottom only
- **Temporal resolution:** Monthly (12 timesteps per year)
- **Spatial domain:** Gulf of Alaska (filtered to depths < 1000m)