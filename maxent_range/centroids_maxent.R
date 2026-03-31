library(raster)
library(ncdf4)
library(data.table)

setwd("~/Desktop/Plants Trait Project")
dataset <- "basic"
species_csv <- "plant_genome/plant_genome.csv"
format <- "default"
path <- "~/Desktop/Plants Trait Project"
output_dir <- "maxent_range/extracted_ranges"

output_csv <- file.path(output_dir, "final_extracted_areas_centroids.csv")
if (!dir.exists(output_dir)) dir.create(output_dir)

metadata <- read.csv(paste0(path, "/", dataset, "/", "metadata_", format, ".csv"))
metadata$scientificname <- tolower(metadata$scientificname)
name_to_idx <- setNames(seq_len(nrow(metadata)), metadata$scientificname)

plant_genome <- read.csv(species_csv)
species_list <- tolower(plant_genome$scientific_name)

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

calc_weighted_centroid <- function(native_ras, pred_ras = NULL) {
  native_mask <- native_ras > 0
  native_mask[is.na(native_mask)] <- 0
  
  # if maxent exists
  if (!is.null(pred_ras) && !all(is.na(values(pred_ras)))) {
    p_native <- pred_ras * native_mask
    p_native[is.na(p_native)] <- 0
    
    w_ras <- area(p_native) * p_native
    df <- as.data.frame(w_ras, xy = TRUE, na.rm = TRUE)
    colnames(df)[3] <- "weight"
  
    if (nrow(df) == 0 || sum(df$weight, na.rm = TRUE) == 0) {
      return(data.frame(
        centroid_x = NA_real_,
        centroid_y = NA_real_,
        centroid_method = "native_plus_maxent"
      ))
    }
    
    total_weight <- sum(df$weight, na.rm = TRUE)
    weighted_x <- sum(df$x * df$weight, na.rm = TRUE) / total_weight
    weighted_y <- sum(df$y * df$weight, na.rm = TRUE) / total_weight
    
    return(data.frame(
      centroid_x = weighted_x,
      centroid_y = weighted_y,
      centroid_method = "native_plus_maxent"
    ))
  }
      
  # otherwise use unweighted centroid of native region only
  native_cells <- rasterToPoints(native_mask, spatial = FALSE)
  native_cells <- native_cells[native_cells[, 3] > 0, , drop = FALSE]
  
  if (nrow(native_cells) == 0) {
    return(data.frame(
      centroid_x = NA_real_,
      centroid_y = NA_real_,
      centroid_method = "native_only"
    ))
  }
  
  return(data.frame(
    centroid_x = mean(native_cells[, 1], na.rm = TRUE),
    centroid_y = mean(native_cells[, 2], na.rm = TRUE),
    centroid_method = "native_only"
  ))
}

target_variable <- c("Native region", "Maxent prediction")

results_list <- list()
for (species in species_list) {
  ras <- read_data(species_name = species,
                   variables = c(target_variable),
                   path = path,
                   dataset = dataset,
                   format = format)
  if (!is.null(ras)) {
    native_layer <- ras[[1]]
    pred_layer   <- ras[[2]]
    centroid <- calc_weighted_centroid(native_layer, pred_layer)
    results_list[[species]] <- data.frame(scientific_name = species,
                                          centroid_x = centroid$centroid_x,
                                          centroid_y = centroid$centroid_y,
                                          centroid_method = centroid$centroid_method,
                                          stringsAsFactors = FALSE
                                          )
  } else {
    results_list[[species]] <- data.frame(scientific_name = species,
                                          centroid_x = NA,
                                          centroid_y = NA,
                                          centroid_method = NA_character_,
                                          stringsAsFactors = FALSE)
  }
}

final_results <- rbindlist(results_list)
write.csv(final_results, output_csv, row.names = FALSE)

sum(!is.na(final_results$centroid_x))
          


