#' Read ROMS Grid Data with Spatial Subset
#' 
#' @description
#' Reads the ROMS grid file and subsets it to match the indices used in the 
#' pre-processed NetCDF files (xi_rho: 67-226, eta_rho: 237-450)
#' 
#' @param grid_file Path to the ROMS grid NetCDF file (default: 'data/NEP_grid_5a.nc')
#' 
#' @return A tibble with columns: xi_rho, eta_rho, lon_rho, lat_rho, h (depth)
#' 
#' @examples
#' roms_grid <- read_roms_grid()
#' 
#' @import tidync
#' @import dplyr
#' @importFrom purrr pluck map_df
read_roms_grid <- function(grid_file = here::here('data', 'NEP_grid_5a.nc')) {
  
  # Open grid file
  roms <- tidync(grid_file)
  
  # Find available grids
  roms_vars <- hyper_grids(roms) %>%
    pluck("grid") %>%
    purrr::map_df(function(x) {
      roms %>% activate(x) %>% hyper_vars() %>% 
        mutate(grd = x)
    })
  
  # Get rho grid
  latlon_rhogrd <- roms_vars %>% filter(name == "lat_rho") %>% pluck('grd')
  roms_rho <- roms %>% 
    activate(latlon_rhogrd) %>% 
    hyper_tibble() %>%
    dplyr::select(lon_rho, lat_rho, xi_rho, eta_rho, h) %>% 
    mutate(lon_rho = lon_rho - 360) # Convert to -180 to 180
  
  # Subset to match the pre-processed data indices
  # Original indices: xi_rho 67-226, eta_rho 237-450
  roms_rho <- roms_rho %>%
    filter(between(xi_rho, 67, 226), between(eta_rho, 237, 450)) %>%
    mutate(xi_rho = xi_rho - 66,  # Renumber starting from 1
           eta_rho = eta_rho - 236)
  
  return(roms_rho)
}


#' Process Single Annual NetCDF File
#' 
#' @description
#' Reads one annual NetCDF file containing surface and bottom layer data for
#' multiple variables, with optional filtering by year range and depth.
#' Files outside the year range are skipped before reading.
#' 
#' @param ncfile Name of the NetCDF file to process (e.g., "annual_1990.nc")
#' @param data_dir Directory containing the NetCDF files
#' @param roms_grid ROMS grid data from read_roms_grid()
#' @param variables Character vector of variable names to extract. 
#'   Default is all 10 variables: temp, salt, PhS, PhL, MZS, MZL, Cop, NCa, Eup, Det
#' @param min_year Minimum year to include in output (required)
#' @param max_year Maximum year to include in output (required)
#' @param maxdepth Maximum depth (h) to include in meters (default: 1000)
#' 
#' @return A tibble with columns: date, variable, layer, lon_rho, lat_rho, value
#'   Returns NULL if file year is outside the specified range
#' 
#' @examples
#' grid <- read_roms_grid()
#' data <- process_annual_file("annual_1990.nc", "data/annual_files", grid,
#'                            min_year = 1990, max_year = 2020)
#' data <- process_annual_file("annual_2010.nc", "data/annual_files", grid, 
#'                            min_year = 2005, max_year = 2015, maxdepth = 500)
#' 
#' @import tidync
#' @import dplyr
#' @import tidyr
#' @import ncdf4
#' @import lubridate
process_annual_file <- function(ncfile, data_dir, roms_grid, 
                                variables = c("temp", "salt", "PhS", "PhL", 
                                              "MZS", "MZL", "Cop", "NCa", 
                                              "Eup", "Det"),
                                min_year = NA,
                                max_year = NA,
                                maxdepth = 1000,
                                mask = goa_mask) {
  
  # Check if min_year and max_year are provided
  if (is.na(min_year) || is.na(max_year)) {
    stop("Both min_year and max_year must be provided.")
  }
  
  # Extract year from filename (e.g., "annual_1990.nc" -> 1990)
  file_year <- as.numeric(gsub("annual_(\\d{4})\\.nc", "\\1", ncfile))
  
  # Check if file year is within range - skip if not
  if (file_year < min_year || file_year > max_year) {
    print(paste("Skipping file:", ncfile, "(outside year range)"))
    return(NULL)
  }
  
  print(paste("Processing file:", ncfile))
  
  # Full path to file
  filepath <- here::here(data_dir, ncfile)
  
  # Open with tidync
  nc <- tidync(filepath)
  
  # Get all variables in long format
  nc_data <- nc %>% 
    hyper_tibble(na.rm = FALSE) %>%
    dplyr::select(xi_rho, eta_rho, ocean_time, s_rho, all_of(variables))
  
  # Convert to long format for easier plotting
  nc_data_long <- nc_data %>%
    pivot_longer(cols = all_of(variables), 
                 names_to = "variable", 
                 values_to = "value") %>%
    drop_na(value)  # Remove NA values (land cells)
  
  # Join with grid to get coordinates
  nc_data_long <- nc_data_long %>%
    left_join(roms_grid, by = c("xi_rho", "eta_rho"))
  
  # Filter by depth using maxdepth parameter
  nc_data_long <- nc_data_long %>% filter(h < maxdepth)
  
  # eliminate xi and eta points outside the ROMS mask
  # this slows down the function but it will produce much smaller masks
  nc_data_long <- nc_data_long %>%
    mutate(idx_drop = paste(xi_rho, eta_rho, sep = "_")) %>%
    filter(!idx_drop %in% mask) %>%
    select(-idx_drop)
  
  # Convert ocean_time to dates
  nc_file <- nc_open(filepath)
  time_data <- ncvar_get(nc_file, "ocean_time")
  time_units <- ncatt_get(nc_file, "ocean_time", "units")$value
  time_parts <- strsplit(time_units, " ")[[1]]
  ref_date_str <- paste(time_parts[3:length(time_parts)], collapse = " ")
  nc_close(nc_file)
  
  # Create date lookup
  dates <- data.frame(
    ocean_time = time_data,
    date = as.POSIXct(time_data, origin = ref_date_str, tz = "UTC")
  ) %>%
    mutate(date = as.Date(date))
  
  # Add dates to data
  nc_data_long <- nc_data_long %>%
    left_join(dates, by = "ocean_time")
  
  # Add layer labels
  nc_data_long <- nc_data_long %>%
    mutate(layer = case_when(
      s_rho == max(s_rho) ~ "surface",
      s_rho == min(s_rho) ~ "bottom",
      TRUE ~ "other"
    ))
  
  # drop unneeded cols and transform to factor where possible
  nc_data_long <- nc_data_long %>%
    select(date, variable, layer, lon_rho, lat_rho, value) %>%
    mutate(
      variable = factor(variable),
      layer = factor(layer)
    )
  
  return(nc_data_long)
}


#' Create Spatial Maps for Surface and Bottom Layers
#' 
#' @description
#' Creates faceted spatial maps showing surface and bottom layers for a 
#' specific variable, year, and month with coastline overlay
#' 
#' @param data Data frame from reading annual files
#' @param variable Variable name to plot (e.g., "temp", "salt")
#' @param year Year to plot (e.g., 2010)
#' @param month Month to plot (1-12)
#' @param coastline sf object containing coastline data from rnaturalearth
#' @param title Plot title (optional, auto-generated if NULL)
#' 
#' @return A ggplot object with faceted maps for surface and bottom layers
#' 
#' @examples
#' # Load coastline data once
#' library(rnaturalearth)
#' coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
#'   st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62)
#' 
#' # Create map for June 2010
#' plot_spatial_map(all_data, "temp", 2010, 6, coast)
#' 
#' # Create map for different variable
#' plot_spatial_map(all_data, "salt", 2015, 8, coast)
#' 
#' @import ggplot2
#' @import dplyr
#' @import sf
#' @import rnaturalearth
plot_spatial_map <- function(data, variable, year, month, coastline, title = NULL, psize = 0.5) {
  
  title <- paste0(variable, " in ", year, "-", month)
  
  dat <- data %>%
    mutate(yr = year(date),
           mo = month(date)) %>%
    filter(variable == !!variable, yr == year, mo == month)
  
  # Create plot
  p <- dat %>%
    st_as_sf(coords = c("lon_rho", "lat_rho"), crs = 4326) %>%
    ggplot() +
    geom_sf(aes(color = value), size = psize, alpha = 0.8) +
    geom_sf(data = coastline, color = "black", linewidth = 0.3) +
    facet_grid(~layer) +
    scale_color_viridis_c(option = "plasma") +
    labs(title = title,
         x = "Longitude",
         y = "Latitude",
         color = variable) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(p)
}

#' Create Time Series Plot with Ribbon
#' 
#' @description
#' Creates a faceted time series plot showing spatial mean with ribbons 
#' representing spatial variability (5th-95th percentiles) for surface 
#' and bottom layers
#' 
#' @param data Data frame from reading annual files
#' @param variable Variable name to plot (e.g., "temp", "salt")
#' @param start_year Starting year for the plot (optional, e.g., 2000)
#' @param end_year Ending year for the plot (optional, e.g., 2010)
#' @param title Plot title (optional, auto-generated if NULL)
#' 
#' @return A ggplot object with faceted panels for surface and bottom layers
#' 
#' @examples
#' # Full time series
#' plot_time_series(all_data, "temp")
#' 
#' # Specific year range
#' plot_time_series(all_data, "temp", start_year = 2000, end_year = 2010)
#' 
#' # Just years after 2005
#' plot_time_series(all_data, "salt", start_year = 2005)
#' 
#' @import ggplot2
#' @import dplyr
#' @import lubridate
plot_time_series <- function(data, variable, start_year = NULL, end_year = NULL, title = NULL) {
  
  # Filter by variable
  plot_data <- data %>%
    filter(variable == !!variable)
  
  # Filter by year range if provided
  if (!is.null(start_year)) {
    start_date <- as.Date(paste0(start_year, "-01-01"))
    plot_data <- plot_data %>%
      filter(date >= start_date)
  }
  
  if (!is.null(end_year)) {
    end_date <- as.Date(paste0(end_year, "-12-31"))
    plot_data <- plot_data %>%
      filter(date <= end_date)
  }
  
  # Calculate spatial statistics with 5th and 95th percentiles
  plot_data <- plot_data %>%
    group_by(date, layer) %>%
    summarise(mean_value = mean(value, na.rm = TRUE),
              lower = quantile(value, 0.05, na.rm = TRUE),
              upper = quantile(value, 0.95, na.rm = TRUE),
              .groups = "drop")
  
  # Set title
  if (is.null(title)) {
    year_range <- if (!is.null(start_year) || !is.null(end_year)) {
      paste0(" (", 
             ifelse(is.null(start_year), "", start_year),
             "-",
             ifelse(is.null(end_year), "", end_year),
             ")")
    } else {
      ""
    }
    title <- paste0("Spatial mean ", variable, " over time", year_range, 
                    "\n(shaded area: 5th-95th percentiles)")
  }
  
  # Create plot with facets
  p <- ggplot(plot_data, aes(x = date, y = mean_value)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "steelblue") +
    geom_line(color = "steelblue", linewidth = 0.8) +
    facet_wrap(~layer, ncol = 1) +
    labs(title = title,
         x = "Date",
         y = paste("Mean", variable)) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"))
  
  return(p)
}
