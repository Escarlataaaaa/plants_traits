library(raster)
library(ncdf4)
library(data.table)

setwd("~/Desktop/Plants Trait Project")
dataset <- "basic"
species_csv <- "plant_genome/plant_genome.csv"
format <- "default"
path <- "~/Desktop/Plants Trait Project"
output_dir <- "maxent_range/extracted_ranges"

output_csv <- file.path(output_dir, "final_extracted_areas_equal_sens_spec.csv")

if (!dir.exists(output_dir)) {dir.create(output_dir)}

# metadata of the suggested (basic) dataset (only the best prediction for each species)
metadata<-read.csv(paste(path, "/", dataset,"/","metadata_",format,".csv",sep=""))

# set species names to lowercase
metadata$scientificname <- tolower(metadata$scientificname)

name_to_idx <- setNames(seq_len(nrow(metadata)), metadata$scientificname)

# get species list from species_csv
plant_genome <- read.csv(species_csv)
species_list <- plant_genome$scientific_name

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
  
  return (crop(ras, e))
}
results_list <- list()
variables <- c("Native region", "Maxent prediction")
total_species <- length(species_list)

for (i in seq_along(species_list)) {
  sp <- species_list[i]
  pct <- (i / total_species) * 100
  cat(sprintf("Processing %d of %d (%.1f%%): %s...\n", i, total_species, pct, sp))
  
  results_list[[sp]] <- tryCatch({
    ras_cropped <- read_data(species_name = sp, variables = variables, path = path, dataset = dataset, format = format)
    
    if (!is.null(ras_cropped)) {
      # get the threshold from metadata to turn probabilities into binary presence/absence
      # using the threshold where sensitivity and specificity are equal 
      idx <- unname(name_to_idx[sp])
      threshold <- metadata$cutoff.equal.sens.spec[idx] 
      if (!is.finite(threshold)) threshold <- 0
      
      native_layer <- ras_cropped[[1]]
      maxent_layer <- ras_cropped[[2]]
      
      # identify where maxent meets the threshold
      is_suitable <- maxent_layer >= threshold
      
      # identify native regions (assuming 1 = native, adjust if the dataset uses >0)
      is_native <- native_layer > 0
      
      # multiply to combine masks: 1 only where BOTH are true, 0 otherwise
      binary_map <- is_suitable * is_native
      binary_map[is.na(binary_map)] <- 0
      
      # calculate cell areas
      cell_areas_km2 <- area(binary_map)
      
      # multiply area by the final combined mask
      presence_area_raster <- cell_areas_km2 * binary_map
      # sum the area of all cells where the species is predicted to be present
      total_area_km2 <- cellStats(presence_area_raster, stat = "sum", na.rm = TRUE)
      
      data.frame(
        scientific_name = sp,
        maxent_area_km2 = total_area_km2,
        threshold_used = threshold,
        stringsAsFactors = FALSE
      )
      
    } else {
      data.frame(scientific_name = sp, maxent_area_km2 = NA, threshold_used = NA, stringsAsFactors = FALSE)    
      }
  
  }, error = function(e) {
    cat("  -> ERROR processing", sp, ":", conditionMessage(e), "\n")
    data.frame(scientific_name = sp, maxent_area_km2 = NA, threshold_used = NA, stringsAsFactors = FALSE)  
  })
}

final_results <- do.call(rbind, results_list)
write.csv(final_results, output_csv, row.names = FALSE)



#debug
sp <- "ficus carica"  # change to one of the 0-area species
sp <- tolower(sp)

ras <- read_data(sp, c("Native region","Maxent prediction"), path, dataset, format)
native_layer <- ras[[1]]
maxent_layer <- ras[[2]]

idx <- unname(name_to_idx[sp])
threshold <- metadata$cutoff.no.omission[idx]

cat("threshold:", threshold, "\n")
cat("maxent min/max:", cellStats(maxent_layer, "min", na.rm=TRUE), cellStats(maxent_layer, "max", na.rm=TRUE), "\n")
cat("native cells:", cellStats(native_layer > 0, "sum", na.rm=TRUE), "\n")
cat("cells >= threshold & native:", cellStats((maxent_layer >= threshold) * (native_layer > 0), "sum", na.rm=TRUE), "\n")
