# ROMS Data Processing and Visualization

This repository contains code for processing and visualizing NEP 10K ROMS ocean model output data for the Gulf of Alaska region. The workflow consists of server-side preprocessing to extract and subset large NetCDF files, followed by local R-based analysis using either raw NetCDF files or preprocessed Parquet files.

**NOTE**: Most GOACLIM models use the approach outlined [here](https://github.com/GOA-CLIM/ROMS_to_Index/tree/main), which produces area- and depth-aggregated indices at the level of NMFS areas or broader GOA regions. You should only be using this repository if you need spatially-explicit ROMS indices.

This code extracts only surface and bottom layers from ROMS NetCDF files. Depth-integrated indices are not currently produced.

This repository uses Git LFS for large data files. Install it before cloning:
- Download from https://git-lfs.github.com/
- Run `git lfs install`

## Workflow Overview

### 1. Server-side Preprocessing (loon)

Run `pre-processing.sh` on the loon server where NetCDF files are hosted. This bash script:
- Extracts surface (layer 42) and bottom (layer 1) data for 10 variables: temp, salt, PhS, PhL, MZS, MZL, Cop, NCa, Eup, Det
- Spatially subsets data to the Gulf of Alaska region (xi_rho: 67-226, eta_rho: 237-450)
- Combines monthly files into annual files
- Outputs processed files to `./annual_files/`

Requirements: CDO (Climate Data Operators)

```bash
./pre-processing.sh
```

The script processes files according to patterns defined for different model runs:
- hindcast: 1990-2020
- historical: 1980-2014
- ssp585, ssp245, ssp126: 2015-2099

#### Location of NetCDF files with monthyl time steps from NEP 10K

*The major subdirectories of /ahr0/hermann/goa-output/ are as follows:*

- cgoa 3km model hindcast: `monthly_aves`

- nep 10km model hindcast: `monthly_aves_nep_revised_hind`

- nep 10km model GFDL "historical" run: `monthly_aves_nep_wb_hist`

- nep 10km model GFDL ssp126 projection: `monthly_aves_nep_wb_ssp126`

- nep 10km model GFDL ssp245 projection: `/ahr2/hermann/monthly_aves_nep_wb_ssp245` (note different path)

- nep 10km model GFDL ssp585 projection: `monthly_aves_nep_wb_ssp585`

### 2. Local R Processing

Transfer the processed annual NetCDF files from loon to your local machine under `data/annual_files/{run}/`, where `{run}` is one of: hindcast, historical, ssp585, ssp245, or ssp126.

## Files

### Core Files

**functions.R** - Contains all core functions for data processing and visualization:
- `read_roms_grid()` - Reads and subsets the ROMS grid file
- `process_annual_file()` - Processes individual annual NetCDF files with filtering options
- `create_parquet_file()` - Converts NetCDF files to Parquet format
- `plot_spatial_map()` - Creates spatial maps for specific months/years
- `plot_time_series()` - Generates time series plots with spatial variability ribbons

### Workflow Scripts

**create_parquet_from_nc.R** - Converts processed NetCDF files to Parquet format. Run this script to create Parquet files for different model runs. Parquet files enable faster data loading and more efficient memory usage for large datasets.

**example.R** - Demonstrates how to work with Parquet files. Shows how to:
- Open Parquet datasets with lazy evaluation
- Filter data before loading into memory
- Create spatial maps and time series visualizations

**create_plot_PICES.R** - Example script that works directly with NetCDF files (bypassing Parquet). Useful when you need:
- Custom filtering not available in preprocessed Parquet files
- To process specific years or variables on-the-fly
- High-resolution figures for publications/presentations

### Reference Scripts

**subset_ROMS_grid.R** - Utility script for identifying grid coordinates (setup only). Contains two sections:
1. Identifies spatial coordinates (xi_rho, eta_rho) for subsetting the NEP 10K ROMS grid to GOA
2. Creates a mask to exclude cells outside NMFS management areas (saves to `idx_to_drop.txt`)

This script only needs to be run if you want to change the spatial subset boundaries. If you modify the spatial subset, you must rerun the server-side extraction and repeat the entire processing workflow.

**pre-processing.sh** - Server-side bash script for initial data extraction (runs on loon).

## Understanding Parquet Files

Parquet is a columnar storage file format optimized for analytical queries. I am testing using this format here because loading the whole ROMS series into your R session will choke your RAM, especially if you load multiple projections at once.

Parquet files support:

- **Lazy evaluation**: Query the data structure without loading it
- **Selective loading**: Read only the columns and rows you need
- **Efficient compression**: Smaller file sizes than RDS files
- **Fast filtering**: Apply filters before loading data into memory

### Working with Parquet Files

```r
# Open dataset without loading into memory
library(arrow)
data <- open_dataset("data/processed/hindcast_annual_data.parquet")

# Apply filters BEFORE loading (much faster and less memory)
filtered_data <- data %>%
  filter(variable == "temp") %>%
  filter(layer == "surface") %>%
  filter(date > as.Date("2010-01-01")) %>%
  collect()  # Now load filtered data into memory
```

The preprocessed Parquet files are filtered to:
- Gulf of Alaska region (based on spatial subset in pre-processing.sh)
- Depths less than or equal to 1000m
- Cells within NMFS management areas (excluding Bering Sea, etc.)

If you need data outside these bounds, you must reprocess the raw NetCDF files using `process_annual_file()` with custom parameters.

## Data Structure

### Processed NetCDF Files
Annual files contain:
- **Variables**: temp, salt, PhS, PhL, MZS, MZL, Cop, NCa, Eup, Det
- **Vertical layers**: Surface (layer 42) and bottom (layer 1) only
- **Temporal resolution**: Monthly (12 timesteps per year, except 1990 with 11)
- **Spatial domain**: Gulf of Alaska subset (xi_rho: 67-226, eta_rho: 237-450)

### Parquet Files
Contain the same data as NetCDF files but with additional filtering:
- **Spatial filtering**: Gulf of Alaska region, depths less than or equal to 1000m
- **Area masking**: Cells outside NMFS management areas excluded
- **Format**: Long format with columns: date, variable, layer, lon_rho, lat_rho, value, run

## Model Runs Available

- **hindcast**: 1990-2020 (11 months in 1990, 12 months for other years)
- **historical**: 1980-2014
- **ssp585**: 2015-2099 (SSP5-8.5 scenario, high emissions)
- **ssp245**: 2015-2099 (SSP2-4.5 scenario, intermediate emissions)
- **ssp126**: 2015-2099 (SSP1-2.6 scenario, low emissions)

NOTE: Projection runs (ssp scenarios) are not yet bias-corrected.

## Important Notes

- This workflow is designed for the rho grid (temperature, salinity, biogeochemical variables). It does not work for velocity components (u, v) which use staggered grids.
- The spatial subset coordinates (xi_rho: 67-226, eta_rho: 237-450) were identified using `subset_ROMS_grid.R`. Changing these requires rerunning the entire preprocessing pipeline.
- Grid indices are renumbered after CDO subsetting. Functions in `functions.R` handle this automatically by using `xi_rho_new` and `eta_rho_new`.
- All Parquet files include a mask that excludes cells outside NMFS management areas. This mask is stored in `idx_to_drop.txt` and applied during processing.

## File Paths

The repository expects the following directory structure:

```
.
├── data/
│   ├── NEP_grid_5a.nc                    # ROMS grid file
│   ├── annual_files/                     # Processed NetCDF files from loon
│   │   ├── hindcast/
│   │   │   └── annual_*.nc
│   │   ├── historical/
│   │   └── ssp585/
│   └── processed/                        # Parquet files
│       ├── hindcast_annual_data.parquet
│       ├── historical_annual_data.parquet
│       └── ssp585_annual_data.parquet
├── functions.R
├── create_parquet_from_nc.R
├── example.R
├── create_plot_PICES.R
├── subset_ROMS_grid.R
├── idx_to_drop.txt                       # Spatial mask
└── pre-processing.sh                     # Run on loon server
```
