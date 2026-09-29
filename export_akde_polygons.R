# export_akde_polygons.R
# Einmalig laufen lassen — exportiert aKDE-Polygone aus Block4b_Ergebnisse.rds
# als GeoJSON damit Python die Populationskarte zeichnen kann.

library(sf)

rds_pfad  <- "/Users/MaintenantPret/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF/output/Block4b_Einzeltier/Block4b_Ergebnisse.rds"
out_pfad  <- "/Users/MaintenantPret/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF/output/Block4b_Einzeltier/akde_polygons.geojson"

res <- readRDS(rds_pfad)
ergebnisse <- res$akde_ergebnisse   # Liste pro Tier

polys <- lapply(ergebnisse, function(x) {
  p95 <- x$poly_95
  p50 <- x$poly_50
  ig  <- x$igel
  a95 <- x$a95ha
  a50 <- x$a50ha

  rows <- list()
  if (!is.null(p95) && inherits(p95, "sf")) {
    p95$igel  <- ig
    p95$level <- "95"
    p95$area_ha <- round(a95, 3)
    rows[[1]] <- p95[, c("igel", "level", "area_ha", "geometry")]
  }
  if (!is.null(p50) && inherits(p50, "sf")) {
    p50$igel  <- ig
    p50$level <- "50"
    p50$area_ha <- round(a50, 3)
    rows[[2]] <- p50[, c("igel", "level", "area_ha", "geometry")]
  }
  if (length(rows) > 0) do.call(rbind, rows) else NULL
})

polys <- Filter(Negate(is.null), polys)
all_polys <- do.call(rbind, polys)

# In WGS84 umrechnen (für einfacheres Plotten in Python)
all_polys_wgs <- st_transform(all_polys, 4326)

st_write(all_polys_wgs, out_pfad, delete_dsn = TRUE)
cat(sprintf("Exportiert: %d Polygone nach %s\n", nrow(all_polys_wgs), out_pfad))
