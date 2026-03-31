# plants_traits
## Maxent range folder
- **[choose_best_threshold](maxent_range/choose_best_threshold.R)**
  Test different threshold and compute a predicted native-range area for each species. Then compare with EOO and choose the threshold with strongest log-correlation.

- **[maxent_range](maxent_range/maxent_range.R)**
  Builds final species range-area using the best threshold. It uses EOO if available, else compute Maxent-derived area.
  
- **[centroids_maxent](maxent_range/centroids_maxent.R)**
  Computes one centroid per species. If maxent predictions exist, computes a centroid weighted by prediction intensity within native region. If not, compute unweighted centroid from native-region cells.
  
- **[map_maxent](maxent_range/map_maxent.R)**
  Visualization of all species centroids on world map with colors for threat status. 
  
