library(raster)
library(ncdf4)
library(data.table)
library(ggplot2)

setwd("~/Desktop/Plants Trait Project")

dataset <- "basic"
species_csv <- "plant_genome/plant_genome.csv"
format <- "default"
path <- "~/Desktop/Plants Trait Project"
iucn_range <- "maxent_range/IUCN_extracted_range.csv"
output_dir <- "maxent_range/extracted_ranges"
output_csv <- file.path(output_dir, "final_extracted_maxent_area_with_eoo.csv")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# metadata
metadata <- read.csv(
  file.path(path, dataset, paste0("metadata_", format, ".csv")),
  stringsAsFactors = FALSE
)
metadata$scientificname <- tolower(metadata$scientificname)
name_to_idx <- setNames(seq_len(nrow(metadata)), metadata$scientificname)

# species list
plant_genome <- read.csv(species_csv, stringsAsFactors = FALSE)
species_list <- tolower(plant_genome$scientific_name)

# IUCN
iucn_df <- fread(iucn_range)
iucn_df[, scientificName := tolower(scientificName)]

# thresholds to evaluate
threshold_candidates <- c(
  "cutoff.kappa",
  "cutoff.spec.sens",
  "cutoff.no.omission",
  "cutoff.prevalence",
  "cutoff.equal.sens.spec",
  "cutoff.sensitivity"
)
threshold_candidates <- threshold_candidates[threshold_candidates %in% names(metadata)]

read_data <- function(species_name, variables, path, dataset, format) {
  i <- unname(name_to_idx[species_name])
  if (is.na(i)) return(NULL)
  
  speciesID <- metadata$speciesID[i]
  ras <- stack()
  
  for (v in variables) {
    ras <- stack(
      ras,
      raster(
        file.path(path, dataset, paste0("range_data_", format, ".nc")),
        varname = v,
        band = speciesID
      )
    )
  }
  
  names(ras) <- variables
  
  e <- extent(
    metadata$extent.xmin[i], metadata$extent.xmax[i],
    metadata$extent.ymin[i], metadata$extent.ymax[i]
  )
  
  crop(ras, e)
}

results_list <- vector("list", length(species_list))
variables <- c("Native region", "Maxent prediction")
total_species <- length(species_list)

for (i in seq_along(species_list)) {
  sp <- species_list[i]
  pct <- 100 * i / total_species
  cat(sprintf("Processing %d of %d (%.1f%%): %s\n", i, total_species, pct, sp))
  
  results_list[[i]] <- tryCatch({
    idx <- unname(name_to_idx[sp])
    
    eoo_value <- iucn_df[scientificName == sp, IUCN_ExtentOfOccurrence_km2][1]
    if (length(eoo_value) == 0) eoo_value <- NA_real_
    
    sp_data <- list(
      scientific_name = sp,
      eoo_area_km2 = eoo_value
    )
    
    # initialize all threshold columns as NA first
    for (thr_type in threshold_candidates) {
      sp_data[[thr_type]] <- NA_real_
      sp_data[[paste0(thr_type, "_val")]] <- NA_real_
    }
    
    if (!is.na(idx)) {
      ras_cropped <- read_data(
        species_name = sp,
        variables = variables,
        path = path,
        dataset = dataset,
        format = format
      )
      
      if (!is.null(ras_cropped)) {
        native_layer <- ras_cropped[[1]]
        maxent_layer <- ras_cropped[[2]]
        
        is_native <- native_layer > 0
        cell_areas_km2 <- area(native_layer)
        
        for (thr_type in threshold_candidates) {
          thr_val <- metadata[[thr_type]][idx]
          
          if (is.finite(thr_val)) {
            is_suitable <- maxent_layer >= thr_val
            binary_map <- is_suitable * is_native
            binary_map[is.na(binary_map)] <- 0
            
            presence_area_raster <- cell_areas_km2 * binary_map
            area_km2 <- cellStats(presence_area_raster, stat = "sum", na.rm = TRUE)
            
            sp_data[[thr_type]] <- area_km2
            sp_data[[paste0(thr_type, "_val")]] <- thr_val
          }
        }
      }
    }
    
    as.data.frame(sp_data, stringsAsFactors = FALSE)
    
  }, error = function(e) {
    cat("  -> ERROR processing", sp, ":", conditionMessage(e), "\n")
    
    eoo_value <- iucn_df[scientificName == sp, IUCN_ExtentOfOccurrence_km2][1]
    if (length(eoo_value) == 0) eoo_value <- NA_real_
    
    sp_data <- list(
      scientific_name = sp,
      eoo_area_km2 = eoo_value
    )
    
    for (thr_type in threshold_candidates) {
      sp_data[[thr_type]] <- NA_real_
      sp_data[[paste0(thr_type, "_val")]] <- NA_real_
    }
    
    as.data.frame(sp_data, stringsAsFactors = FALSE)
  })
}

all_areas_df <- rbindlist(results_list, fill = TRUE)

write.csv(all_areas_df, output_csv, row.names = FALSE)


# choose threshold that matched EOO the best

eval_df <- copy(all_areas_df)
eval_df[, eoo_area_km2 := as.numeric(gsub(",", "", as.character(eoo_area_km2)))]

score_table <- eval_df[, {
  log_eoo <- log(eoo_area_km2 + 1)
  
  rbindlist(lapply(threshold_candidates, function(thr) {
    valid <- is.finite(eoo_area_km2) & is.finite(get(thr))
    log_maxent <- log(get(thr) + 1)
    
    data.table(
      threshold_candidates = thr,
      n_species = sum(valid),
      correlation_log = cor(log_eoo[valid], log_maxent[valid]),
      mean_abs_log_error = mean(abs(log_eoo[valid] - log_maxent[valid]))
    )
  }))
}]

score_table <- score_table[order(-correlation_log, mean_abs_log_error)]
score_table

best_threshold <- score_table$threshold_candidates[1]
best_r <- score_table$correlation_log[1]

# make plotting dataframe
plot_df <- copy(eval_df)[
  is.finite(eoo_area_km2) & is.finite(get(best_threshold)),
  .(
    scientific_name,
    eoo_area_km2,
    maxent_area_km2 = get(best_threshold)
  )
]

# log columns
plot_df[, `:=`(
  log_eoo = log(eoo_area_km2 + 1),
  log_maxent = log(maxent_area_km2 + 1)
)]

# plot
corr_plot <- ggplot(plot_df, aes(x = log_eoo, y = log_maxent)) +
  geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    title = paste0("EOO vs Maxent Area (", best_threshold, ")"),
    subtitle = paste0("Pearson r = ", round(best_r, 3), 
                      " | n = ", nrow(plot_df)),
    x = "log(EOO area + 1)",
    y = paste0("log(Maxent area + 1) [", best_threshold, "]")
  ) +
  theme_minimal()

print(corr_plot)

ggsave(
  file.path(output_dir, "best_threshold_correlation.png"),
  corr_plot,
  width = 7,
  height = 7,
  dpi = 300,
  bg = "white"
)

# compare EOO of other species not studied and compare with threshold to show no bias

