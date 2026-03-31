library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

setwd("~/Desktop/Plants Trait Project")

status_colors <- c(
  "LC" = "#0c7345",
  "NT" = "#6cca0e",
  "VU" = "#F9E814",
  "EN" = "#FC7F3F",
  "CR" = "#D81E05",
  "EW" = "#4A0E1E"
)

centroids <- read.csv("maxent_range/extracted_ranges/final_extracted_areas_centroids.csv")
genome <- read.csv("plant_genome/plant_genome.csv")

genome <- genome %>%
  mutate(iucn_status = case_when(
    redlistCategory == "Least Concern" ~ "LC",
    redlistCategory == "Lower Risk/least concern" ~ "LC",
    redlistCategory == "Near Threatened" ~ "NT",
    redlistCategory == "Lower Risk/near threatened" ~ "NT",
    redlistCategory == "Vulnerable" ~ "VU",
    redlistCategory == "Endangered" ~ "EN",
    redlistCategory == "Critically Endangered" ~ "CR",
    redlistCategory == "Extinct in the Wild" ~ "EW",
    TRUE ~ NA_character_
  ))

centroids <- centroids %>%
  filter(!is.na(centroid_x), !is.na(centroid_y))

centroids <- centroids %>%
  left_join(
    genome %>% select(scientific_name, iucn_status),
    by = "scientific_name"
  )

centroids_sf <- st_as_sf(centroids, coords = c("centroid_x", "centroid_y"), crs = 4326)
world <- ne_countries(scale = "medium", returnclass = "sf")

map <- ggplot() +
  
  geom_sf(data = world,
          fill = "#cfe8d6",
          color = "#9dbfa7",
          size = 0.2) +
  
  geom_sf(data = centroids_sf,
          aes(color = iucn_status),
          size = 3,
          alpha = 0.9) +
  
  scale_color_manual(
    values = status_colors,
    name = "IUCN Status"
  ) +
  
  coord_sf(expand = FALSE) +
  
  theme_minimal() +
  
  theme(
    panel.background = element_rect(fill = "#e6f2ff"),
    plot.background = element_rect(fill = "#e6f2ff"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "white")
  )

print(map)

ggsave(
  "maxent_range/figures/weighted_centroids_global_map.pdf",
  map,
  width = 16,
  height = 8,
  dpi = 600
)
