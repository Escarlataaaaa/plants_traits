library(raster)
library(ncdf4)
library(data.table)

setwd("~/Desktop/Plants Trait Project")
dataset <- "basic"
species_csv <- "plant_genome/plant_genome.csv"
format <- "default"
path <- "~/Desktop/Plants Trait Project"
output_dir <- "maxent_range/extracted_ranges"
iucn_range <- "maxent_range/IUCN_extracted_range.csv"

output_csv <- file.path(output_dir, "final_extracted_areas_no_omission.csv")

if (!dir.exists(output_dir)) {dir.create(output_dir)}

# metadata of the suggested (basic) dataset (only the best prediction for each species)
metadata<-read.csv(paste(path, "/", dataset,"/","metadata_",format,".csv",sep=""))

# set species names to lowercase
metadata$scientificname <- tolower(metadata$scientificname)

name_to_idx <- setNames(seq_len(nrow(metadata)), metadata$scientificname)

# get species list from species_csv
plant_genome <- read.csv(species_csv)
species_list <- tolower(plant_genome$scientific_name)

# IUCN file
iucn_df <- fread(iucn_range)
iucn_df[, scientificName := tolower(scientificName)]

read_data<-function(species_name, variables, path, dataset, format){
  i <- unname(name_to_idx[species_name])
  if (is.na(i)) return(NULL)
  
  speciesID <- metadata$speciesID[i]
  
  ras<-stack()
  
  for(v in variables){
    ras<-stack(ras,raster(paste(path, "/",dataset,"/","range_data_",format,".nc",sep=""), varname = v, band = speciesID))
  }
  
  names(ras)<-variables
  
  e <- extent(metadata$extent.xmin[i], metadata$extent.xmax[i],
              metadata$extent.ymin[i], metadata$extent.ymax[i])
  
  return(crop(ras, e))
}

results_list <- list()
variables <- c("Native region", "Maxent prediction")
total_species <- length(species_list)

for (i in seq_along(species_list)) {
  sp <- species_list[i]
  pct <- (i / total_species) * 100
  cat(sprintf("Processing %d of %d (%.1f%%): %s...\n", i, total_species, pct, sp))
  
  results_list[[sp]] <- tryCatch({
    ras_cropped <- read_data(
      species_name = sp,
      variables = variables,
      path = path,
      dataset = dataset,
      format = format
    )
    
    eoo_value <- iucn_df[scientificName == sp, IUCN_ExtentOfOccurrence_km2][1]
    if (length(eoo_value) == 0) eoo_value <- NA_real_
    
    if (!is.null(ras_cropped)) {
      idx <- unname(name_to_idx[sp])
      threshold <- metadata$cutoff.no.omission[idx]
      if (!is.finite(threshold)) threshold <- 0
      
      native_layer <- ras_cropped[[1]]
      maxent_layer <- ras_cropped[[2]]
      
      # identify where maxent meets the threshold
      is_suitable <- maxent_layer >= threshold
      
      # identify native regions
      is_native <- native_layer > 0
      
      # 1 only where both are true
      binary_map <- is_suitable * is_native
      binary_map[is.na(binary_map)] <- 0
      
      # calculate cell areas
      cell_areas_km2 <- area(binary_map)
      
      # sum the area of all cells where the species is predicted to be present
      presence_area_raster <- cell_areas_km2 * binary_map
      maxent_area_km2 <- cellStats(presence_area_raster, stat = "sum", na.rm = TRUE)
      
    } else {
      threshold <- NA_real_
      maxent_area_km2 <- NA_real_
    }
    final_area_km2 <- ifelse(!is.na(eoo_value), eoo_value, maxent_area_km2)
    area_source <- ifelse(!is.na(eoo_value), "EOO", "Maxent")
    
    data.frame(
      scientific_name = sp,
      eoo_area_km2 = eoo_value,
      maxent_area_km2 = maxent_area_km2,
      final_area_km2 = final_area_km2,
      area_source = area_source,
      threshold_used = threshold,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("  -> ERROR processing", sp, ":", conditionMessage(e), "\n")
    
    eoo_value <- iucn_df[scientificName == sp, IUCN_ExtentOfOccurrence_km2][1]
    if (length(eoo_value) == 0) eoo_value <- NA_real_
    
    data.frame(
      scientific_name = sp,
      eoo_area_km2 = eoo_value,
      maxent_area_km2 = NA_real_,
      final_area_km2 = eoo_value,
      area_source = ifelse(!is.na(eoo_value), "EOO", NA),
      threshold_used = NA_real_,
      stringsAsFactors = FALSE
    )
  })
}


final_results <- do.call(rbind, results_list)
write.csv(final_results, output_csv, row.names = FALSE)

sum(!is.na(final_results$maxent_area_km2))
