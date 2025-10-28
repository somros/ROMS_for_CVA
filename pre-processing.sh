#!/bin/bash
# -------------------------------------------------------------------------
# Run this on loon
# Combined script to:
# 1. Extract surface and bottom temperature layers with spatial subsetting
# 2. Combine monthly files into annual files
# -------------------------------------------------------------------------

echo "=== ROMS Data Processing Pipeline ==="
echo "Step 1: Extracting surface and bottom layers..."

# Create output directory for subsetted files
mkdir -p ./bottom_surface_subset

# Process each file from the list to extract surface and bottom layers
while read file_path; do
    # Extract filename without path
    filename=$(basename "$file_path")

    # Create output file name
    output_file="./bottom_surface_subset/b_s_roms_${filename}"

    # Extract temp variable at surface (42) and bottom (1) levels
    # With spatial subsetting using your specific indices
    # I got these manually from looking at the NC grid file and working out which xi and eta rho correspond to the area of interest
    cdo -f nc4 -z zip_9 selindexbox,67,226,237,450 -sellevidx,1,42 -selname,temp,salt,Det,Eup,NCa,MZL,PhL,Cop,MZS,PhS "$file_path" "$output_file"

    echo "Processed $filename"
done < ssp585_files.txt

echo "Step 1 complete: Surface and bottom layer extraction finished"
echo ""
echo "Step 2: Combining monthly files into annual files..."

# -------------------------------------------------------------------------
# Annual file creation section
# -------------------------------------------------------------------------

# Directory containing the subsetted NetCDF files
INPUT_DIR="./bottom_surface_subset"
# Directory for annual output files
OUTPUT_DIR="./annual_files"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Years to process: 
# hindcast 1990-2020
# historical run: 1980-2014
# projections: 2015-2099
START_YEAR=2015
END_YEAR=2099

# Process each year
for YEAR in $(seq $START_YEAR $END_YEAR); do
    echo "Processing year $YEAR..."
    
    # Create a temporary directory for this year's files
    TEMP_DIR="./temp_${YEAR}"
    mkdir -p $TEMP_DIR
    
    # Find all monthly files for this year and copy them to temp directory
    # Files follow pattern: surf_bottom_roms_nep_revised_hind_moave_YYYY_MM.nc
    FOUND_FILES=0
    for MONTH in $(seq -w 1 12); do
        MONTHLY_FILE="$INPUT_DIR/b_s_roms_nep_wb_ssp585_moave_${YEAR}_${MONTH}.nc"
        if [ -f "$MONTHLY_FILE" ]; then
            cp "$MONTHLY_FILE" "$TEMP_DIR/"
            FOUND_FILES=$((FOUND_FILES + 1))
            echo "  Found: $(basename $MONTHLY_FILE)"
        else
            echo "  Missing: b_s_roms_nep_wb_ssp585_moave_${YEAR}_${MONTH}.nc"
        fi
    done
    
    # Check if we found any files for this year
    if [ $FOUND_FILES -eq 0 ]; then
        echo "  No files found for year $YEAR."
        rm -rf $TEMP_DIR
        continue
    fi
    
    echo "  Found $FOUND_FILES monthly files for $YEAR"
    
    # Merge all monthly files for this year
    echo "  Merging monthly files for $YEAR..."
    cdo mergetime $TEMP_DIR/b_s_roms_nep_wb_ssp585_moave_${YEAR}_*.nc $OUTPUT_DIR/annual_${YEAR}.nc
    
    # Check how many timesteps we got
    NUM_TIMESTEPS=$(cdo ntime $OUTPUT_DIR/annual_${YEAR}.nc)
    echo "  Created annual file with $NUM_TIMESTEPS timesteps"
    
    if [ $NUM_TIMESTEPS -ne 12 ]; then
        echo "  WARNING: Expected 12 monthly timesteps, but got $NUM_TIMESTEPS"
        # Show which months are present
        echo "  Timestamps in file:"
        cdo showtimestamp $OUTPUT_DIR/annual_${YEAR}.nc
    else
        echo "  Complete: All 12 months present for year $YEAR"
    fi
    
    # Clean up temporary directory
    rm -rf $TEMP_DIR
done

echo ""
echo "=== Processing Pipeline Complete ==="
echo "- Subsetted files (surface + bottom layers): ./surface_subset/"
echo "- Annual files: ./annual_files/"
echo "- Each annual file contains all variables with 2 vertical layers and 12 monthly timesteps (except 1990 that has 11)"