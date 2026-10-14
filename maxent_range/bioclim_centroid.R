library(geodata)
library(raster)
library(ncdf4)
library(data.table)
library(dplyr)

setwd("~/Desktop/Plants Trait Project")

dataset <- "basic"
format <- "default"

project_path <- "~/Desktop/Plants Trait Project"
basic_path   <- file.path(project_path, "basic")
species_csv <- file.path(project_path, "plant_genome/plant_genome.csv")
metadata_file <- file.path(basic_path, "metadata_default.csv")

metadata <- read.csv(metadata_file)
metadata$scientificname <- tolower(trimws(metadata$scientificname))
name_to_idx <- setNames(seq_len(nrow(metadata)), metadata$scientificname)

plant_genome <- read.csv(species_csv)
species_list <- tolower(trimws(plant_genome$scientific_name))

bio <- worldclim_global(
  var = "bio",
  res = 10,
  path = file.path(project_path, "climate_data")
)

# take one species name and pull out the two raster layers (native and maxent) for that species from range data
get_species_layers <- function(species_name) {
  species_name <- tolower(trimws(species_name))
  i <- name_to_idx[species_name]
  
  if (is.na(i)) return(NULL)
  
  speciesID <- metadata$speciesID[i]
  nc_file <- file.path(project_path, dataset, paste0("range_data_", format, ".nc"))
  
  native_layer <- raster(nc_file, varname = "Native region", band = speciesID)
  maxent_layer <- raster(nc_file, varname = "Maxent prediction", band = speciesID)
  
  e <- extent(
    metadata$extent.xmin[i], metadata$extent.xmax[i],
    metadata$extent.ymin[i], metadata$extent.ymax[i]
  )
  
  native_layer <- crop(native_layer, e)
  maxent_layer <- crop(maxent_layer, e)
  
  list(
    native_layer = native_layer,
    maxent_layer = maxent_layer
  )
}

# Compute one weighted lon and lat for a species based on centroids_maxent.R and its results
# Run centroids_maxent.R first to get the final_extracted_areas_centroids.csv
centroid_df <- read.csv("maxent_range/final_extracted_areas_centroids.csv", stringsAsFactors = FALSE)
centroid_df$scientific_name <- tolower(trimws(centroid_df$scientific_name))

get_weighted_coords <- function(species_name) {
  row <- centroid_df[centroid_df$scientific_name == species_name, ]
  
  if (nrow(row) == 0) {
    return(data.frame(
      weighted_lon = NA_real_,
      weighted_lat = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    weighted_lon = row$centroid_x[1],
    weighted_lat = row$centroid_y[1],
    stringsAsFactors = FALSE
  )
}

# Compute weighted mean of one BIO climate variable for one species
get_bioclim_value <- function(species_name, bio_index = 19) {
  layers <- get_species_layers(species_name)
  
  if (is.null(layers)) return(NA_real_)
  
  native_layer <- layers$native_layer
  maxent_layer <- layers$maxent_layer
  
  bio_layer <- raster(bio[[bio_index]])
  
  if (!compareCRS(bio_layer, maxent_layer)) {
    bio_layer <- projectRaster(bio_layer, crs = crs(maxent_layer), method = "bilinear")
  }
  
  bio_crop  <- crop(bio_layer, extent(maxent_layer))
  bio_align <- resample(bio_crop, maxent_layer, method = "bilinear")
  
  maxent_use <- mask(maxent_layer, native_layer, maskvalue = 0)
  bio_use    <- mask(bio_align, native_layer, maskvalue = 0)
  
  suit_vals <- getValues(maxent_use)
  bio_vals  <- getValues(bio_use)
  
  keep <- !is.na(suit_vals) & !is.na(bio_vals) & suit_vals > 0
  
  if (!any(keep)) return(NA_real_)
  
  weighted.mean(bio_vals[keep], suit_vals[keep])
}

# __main__
# Start process to collect all 19 BIO variable
all_results_list <- vector("list", length(species_list))

cat("Starting run for", length(species_list), "species and 19 BIO variables\n")

for (s in seq_along(species_list)) {
  sp <- species_list[s]
  cat("Species", s, "of", length(species_list), ":", sp, "\n")
  
  coords <- tryCatch(
    get_weighted_coords(sp),
    error = function(e) {
      cat("   COORD ERROR:", conditionMessage(e), "\n")
      data.frame(
        weighted_lon = NA_real_,
        weighted_lat = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
  
  row_out <- data.frame(
    scientific_name = sp,
    weighted_lon = coords$weighted_lon[1],
    weighted_lat = coords$weighted_lat[1],
    stringsAsFactors = FALSE
  )
  
  for (j in 1:19) {
    cat("   -> BIO", j, "\n")
    
    bio_val <- tryCatch(
      get_bioclim_value(sp, bio_index = j),
      error = function(e) {
        cat("      ERROR:", conditionMessage(e), "\n")
        NA_real_
      }
    )
    
    row_out[[paste0("BIO", j, "_weighted_bioclim")]] <- bio_val
  }
  
  all_results_list[[s]] <- row_out
  cat("Finished species:", sp, "\n\n")
}

all_results <- do.call(rbind, all_results_list)

write.csv(
  all_results,
  file.path(project_path, "all_19_bioclim_weighted_summary_all_species.csv"),
  row.names = FALSE
)

