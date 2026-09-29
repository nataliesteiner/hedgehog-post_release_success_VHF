# ==============================================================================
#  Block7_DetektionsValidierung.R
#  Hedgehog VHF telemetry, Wildtierstation Sachsenhagen
#
#  Question: Are all hedgehog locations actually detected, or are there
#            systematic gaps due to limited transmitter range?
#
#  Relevance: Critical prerequisite for interpreting the kernel home-range
#             results (Block4). Reviewers will ask this question.
#
#  Analyses:
#    1. Detection distance:     distance of each fix to the nearest antenna
#    2. Detection gaps:         nights without signal during active period
#    3. Station count:          share of fixes with only 2 receiving stations
#    4. Detection bias enclosure: blind spot near the release site?
#    5. Word methods report:    methods + limitations (officer/flextable)
#    6. Word results report:    quantitative results + interpretation (officer)
#
#  Input:
#    data/kernel_files/multilateration_April2026/*.multilaterations.gpkg

# ── 0. Pakete ──────────────────────────────────────────────────────────────────

pakete <- c("sf", "data.table", "ggplot2", "patchwork", "scales",
            "readxl", "openxlsx", "dplyr", "lubridate",
            "officer", "flextable")
neu <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) {
  cat("Installiere fehlende Pakete:", paste(neu, collapse = ", "), "\n")
  install.packages(neu)
}

suppressPackageStartupMessages({
  library(sf)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(lubridate)
})
cat("✓ Pakete geladen\n\n")

# ── Hilfsfunktion: sichere Spaltensuche ───────────────────────────────────────
# Gibt NA_character_ zurück (nie character(0)) — verhindert is.na()-Crashes
get_col <- function(df, patterns) {
  found <- unlist(lapply(patterns, function(p)
    grep(p, names(df), ignore.case = TRUE, value = TRUE)))
  if (length(found) > 0) found[1] else NA_character_
}

# ── 1. Pfade ───────────────────────────────────────────────────────────────────
projekt_root <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"

gpkg_ordner  <- file.path(projekt_root, "data", "kernel_files")
meta_datei   <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
out_ordner   <- file.path(projekt_root, "output", "Block7_DetektionsValidierung")
dir.create(out_ordner, showWarnings = FALSE, recursive = TRUE)

cat("📁 GPKG-Ordner:", gpkg_ordner, "\n")
cat("📁 Output:     ", out_ordner, "\n\n")

# ── 2. Antennen-Koordinaten ───────────────────────────────────────────────────
# tRackIT Stationen Wildtierstation Sachsenhagen (WGS84)
ant_wgs84 <- data.frame(
  station_id = c("A1",       "A2",       "A3",       "A4",       "A5",       "A6"),
  lat        = c(52.398057,  52.396738,  52.397434,  52.399223,  52.398802,  52.400136),
  lon        = c( 9.217069,   9.216513,   9.215485,   9.217175,   9.214216,   9.215782)
)

ant_sf    <- st_as_sf(ant_wgs84, coords = c("lon", "lat"), crs = 4326)
ant_utm   <- st_transform(ant_sf, crs = 25832)
ant_coords <- st_coordinates(ant_utm)

ANTENNEN <- data.frame(
  station_id = ant_wgs84$station_id,
  x_utm      = ant_coords[, 1],
  y_utm      = ant_coords[, 2]
)

cat(sprintf("✓ %d Antennen geladen.\n", nrow(ANTENNEN)))
cat("  Stationen:", paste(ANTENNEN$station_id, collapse = ", "), "\n\n")

COL_IGEL <- "#4C9A52"

# ── 3. Metadaten laden ─────────────────────────────────────────────────────────
meta_raw <- as.data.table(read_excel(meta_datei))
names(meta_raw) <- tolower(gsub("[^a-zA-Z0-9]", "_", names(meta_raw)))

# Bug-Fix: get_col gibt NA_character_ zurück, nie character(0)
col_indiv <- get_col(meta_raw, c("individual", "igel", "tier", "animal"))
col_rel   <- get_col(meta_raw, c("date_release", "release_date", "auswilderung"))
col_last  <- get_col(meta_raw, c("date_last_signal", "last_signal", "letztes_signal"))

if (is.na(col_indiv) || is.na(col_rel))
  stop("[FEHLER] Spalten 'individual' oder 'date_release' nicht in Metadaten gefunden.")

meta <- meta_raw %>%
  transmute(
    tier_label       = trimws(.data[[col_indiv]]),
    date_release     = as.Date(.data[[col_rel]]),
    date_last_signal = if (!is.na(col_last)) as.Date(.data[[col_last]]) else as.Date(NA)
  ) %>%
  filter(!is.na(tier_label), !is.na(date_release))

cat(sprintf("✓ Metadaten: %d Tiere geladen.\n\n", nrow(meta)))

# ── 4. GeoPackages laden ───────────────────────────────────────────────────────
cat("── Lade alle Multilateralisierungs-Fixes ────────────────────\n")

# Bug-Fix: recursive = TRUE damit Unterordner (z. B. multilateration_April2026/) durchsucht werden
gpkg_dateien <- list.files(gpkg_ordner,
                            pattern    = "\\.gpkg$",
                            full.names = TRUE,
                            recursive  = TRUE)
cat(sprintf("  %d GeoPackage-Dateien gefunden.\n", length(gpkg_dateien)))

alle_fixes <- rbindlist(lapply(gpkg_dateien, function(fp) {
  tryCatch({
    sf_obj  <- st_read(fp, quiet = TRUE)

    bn      <- basename(fp)
    lbl_m   <- regmatches(bn, regexpr("Igel\\d+", bn, perl = TRUE))
    lbl     <- if (length(lbl_m) > 0) lbl_m[1] else sub("_.*", "", bn)

    if (is.na(st_crs(sf_obj))) sf_obj <- st_set_crs(sf_obj, 4326)
    if (st_crs(sf_obj)$epsg != 25832)
      sf_obj <- st_transform(sf_obj, crs = 25832)
    sf_obj <- sf_obj[!st_is_empty(sf_obj), ]
    if (nrow(sf_obj) == 0) return(NULL)
    coords <- st_coordinates(sf_obj)

    # Bug-Fix: get_col statt direktem grep → kein character(0)-Crash
    time_col <- get_col(sf_obj, c("X_time", "_time", "time", "timestamp", "Timestamp"))
    if (is.na(time_col)) {
      if (all(c("Date", "Time") %in% names(sf_obj))) {
        sf_obj$ts <- as.POSIXct(paste(sf_obj$Date, sf_obj$Time), tz = "UTC")
        time_col <- "ts"
      } else {
        cat(sprintf("  [WARN] Keine Zeitspalte in %s\n", bn))
        return(NULL)
      }
    }
    ts <- as.POSIXct(sf_obj[[time_col]], tz = "UTC")

    # Bug-Fix: get_col für Station Count
    sc_col  <- get_col(sf_obj, c("Station Count", "station_count", "StationCount",
                                  "station.count", "stationcount"))
    sc_vals <- if (!is.na(sc_col)) as.integer(sf_obj[[sc_col]]) else NA_integer_

    dt <- data.table(
      tier_label    = lbl,
      timestamp     = ts,
      date          = as.Date(ts, tz = "Europe/Berlin"),
      x_utm         = coords[, 1],
      y_utm         = coords[, 2],
      station_count = sc_vals
    )
    dt <- dt[!is.na(x_utm) & !is.na(y_utm) & !is.na(timestamp)]
    cat(sprintf("  ✓ %s: %d Fixes\n", lbl, nrow(dt)))
    dt

  }, error = function(e) {
    cat(sprintf("  [WARN] Fehler in %s: %s\n", basename(fp), conditionMessage(e)))
    NULL
  })
}))

if (is.null(alle_fixes) || nrow(alle_fixes) == 0)
  stop("[FEHLER] Keine Fixes geladen. Bitte gpkg_ordner prüfen.")

alle_fixes <- merge(alle_fixes,
                    meta[, .(tier_label, date_release, date_last_signal)],
                    by = "tier_label", all.x = TRUE)

alle_fixes <- alle_fixes[
  !is.na(date_release) &
    date >= date_release &
    (is.na(date_last_signal) | date <= date_last_signal)
]

cat(sprintf("\n  Gesamt (aktive Periode): %d Fixes von %d Tieren\n",
            nrow(alle_fixes), n_distinct(alle_fixes$tier_label)))

# ══════════════════════════════════════════════════════════════════════════════
#  ANALYSE 1 — Detektionsdistanz
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Analyse 1: Detektionsdistanz ─────────────────────────────\n")

ant_mat <- as.matrix(ANTENNEN[, c("x_utm", "y_utm")])

alle_fixes[, dist_nearest_ant_m := {
  fix_mat <- matrix(c(x_utm, y_utm), ncol = 2)
  apply(fix_mat, 1, function(pt) {
    min(sqrt((ant_mat[, 1] - pt[1])^2 + (ant_mat[, 2] - pt[2])^2))
  })
}]

dist_stats <- alle_fixes[, .(
  N_fixes  = .N,
  Median_m = round(median(dist_nearest_ant_m, na.rm = TRUE)),
  Mean_m   = round(mean(dist_nearest_ant_m, na.rm = TRUE)),
  P75_m    = round(quantile(dist_nearest_ant_m, 0.75, na.rm = TRUE)),
  P90_m    = round(quantile(dist_nearest_ant_m, 0.90, na.rm = TRUE)),
  P95_m    = round(quantile(dist_nearest_ant_m, 0.95, na.rm = TRUE)),
  Max_m    = round(max(dist_nearest_ant_m, na.rm = TRUE))
)]

cat("\n  Distanz zur nächsten Antenne [m]:\n")
print(dist_stats)
cat(sprintf("\n  ► Median: %d m | 95. Pz.: %d m | Max: %d m\n",
            dist_stats$Median_m, dist_stats$P95_m, dist_stats$Max_m))

if (dist_stats$P95_m < 300) {
  cat("  [OK] 95% aller Fixes liegen innerhalb von 300 m zur nächsten Antenne.\n")
} else {
  cat("  [!] 5% der Fixes liegen > 300 m von der nächsten Antenne entfernt.\n")
  cat("      → In Publikation als Limitation diskutieren.\n")
}

dist_pro_tier <- alle_fixes[, .(
  N_fixes  = .N,
  Median_m = round(median(dist_nearest_ant_m, na.rm = TRUE)),
  P90_m    = round(quantile(dist_nearest_ant_m, 0.90, na.rm = TRUE)),
  Max_m    = round(max(dist_nearest_ant_m, na.rm = TRUE))
), by = tier_label][order(tier_label)]

# Plot 1a: Dichteverteilung
p1a <- ggplot(alle_fixes, aes(x = dist_nearest_ant_m)) +
  geom_density(fill = COL_IGEL, color = COL_IGEL, alpha = 0.35, linewidth = 1.1) +
  geom_vline(xintercept = dist_stats$P95_m,
             linetype = "dashed", color = "firebrick", linewidth = 0.9) +
  annotate("label", x = dist_stats$P95_m, y = Inf,
           label = paste0("95. Pz. = ", dist_stats$P95_m, " m"),
           vjust = 1.5, hjust = -0.05, size = 3.5,
           fill = "white", color = "firebrick", label.size = 0.3) +
  geom_vline(xintercept = dist_stats$Median_m,
             linetype = "dotted", color = "grey40", linewidth = 0.9) +
  annotate("label", x = dist_stats$Median_m, y = Inf,
           label = paste0("Median = ", dist_stats$Median_m, " m"),
           vjust = 3.0, hjust = 1.05, size = 3.2,
           fill = "white", color = "grey40", label.size = 0) +
  scale_x_continuous(labels = label_number(suffix = " m"),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title    = "Detektionsdistanz — Abstand jedes Fixes zur nächsten Antenne",
       subtitle = sprintf("n = %d Fixes von %d Tieren  |  Dashed = 95. Perzentile",
                          nrow(alle_fixes), n_distinct(alle_fixes$tier_label)),
       x = "Abstand zur nächsten Antenne (m)", y = "Dichte") +
  theme_classic(base_size = 13) +
  theme(panel.grid.major.y = element_line(color = "grey93", linewidth = 0.4))

# Plot 1b: ECDF
p1b <- ggplot(alle_fixes, aes(x = dist_nearest_ant_m)) +
  stat_ecdf(color = COL_IGEL, linewidth = 1.3) +
  geom_vline(xintercept = dist_stats$P95_m,
             linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept = 0.95,
             linetype = "dashed", color = "firebrick", linewidth = 0.5, alpha = 0.6) +
  scale_x_continuous(labels = label_number(suffix = " m")) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(title    = "Kumulative Verteilung der Detektionsdistanz",
       subtitle = "Zeigt bei welchem Abstand welcher Anteil der Fixes liegt",
       x = "Abstand zur nächsten Antenne (m)", y = "Kumulativer Anteil der Fixes") +
  theme_classic(base_size = 13)

p_range <- p1a / p1b
ggsave(file.path(out_ordner, "07a_detection_distance.png"), p_range,
       width = 10, height = 10, dpi = 200, bg = "white")
cat("  ✓ 07a_detection_distance.png\n")

# ══════════════════════════════════════════════════════════════════════════════
#  ANALYSE 2 — Detektionslücken
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Analyse 2: Detektionslücken ──────────────────────────────\n")

nacht_grid <- rbindlist(lapply(unique(alle_fixes$tier_label), function(lbl) {
  sub_meta <- meta[tier_label == lbl]
  if (nrow(sub_meta) == 0 || is.na(sub_meta$date_release)) return(NULL)
  last_dt <- if (!is.na(sub_meta$date_last_signal))
               sub_meta$date_last_signal
             else
               max(alle_fixes[tier_label == lbl]$date)
  all_dates <- seq.Date(sub_meta$date_release, last_dt, by = "day")
  data.table(tier_label = lbl, date = all_dates)
}))

detektiert <- alle_fixes[, .(detected = TRUE), by = .(tier_label, date)]
nacht_grid <- merge(nacht_grid, detektiert, by = c("tier_label", "date"), all.x = TRUE)
nacht_grid[is.na(detected), detected := FALSE]

luecken_stats <- nacht_grid[, .(
  total_days       = .N,
  days_detected    = sum(detected),
  days_missing     = sum(!detected),
  detection_rate   = round(mean(detected) * 100, 1)
), by = tier_label][order(tier_label)]

cat("\n  Detektionsrate pro Tier:\n")
print(luecken_stats)
cat(sprintf("\n  Median Detektionsrate: %.1f%%\n",
            median(luecken_stats$detection_rate, na.rm = TRUE)))
cat(sprintf("  Min: %.1f%%  |  Max: %.1f%%\n",
            min(luecken_stats$detection_rate),
            max(luecken_stats$detection_rate)))

schlechte_det <- luecken_stats[detection_rate < 50]
if (nrow(schlechte_det) > 0) {
  cat(sprintf("  [!] %d Tier(e) mit < 50%% Detektionsrate: %s\n",
              nrow(schlechte_det),
              paste(schlechte_det$tier_label, collapse = ", ")))
  cat("      → Möglicherweise außerhalb Detektionsbereich oder Sender ausgefallen.\n")
}

p2 <- ggplot(luecken_stats,
             aes(x = detection_rate, y = reorder(tier_label, detection_rate))) +
  geom_segment(aes(x = 0, xend = detection_rate,
                   yend = reorder(tier_label, detection_rate)),
               color = "grey80", linewidth = 0.5) +
  geom_point(aes(fill = detection_rate), size = 4, shape = 21, color = "grey30") +
  geom_vline(xintercept = median(luecken_stats$detection_rate),
             linetype = "dashed", color = COL_IGEL, linewidth = 0.9) +
  annotate("label",
           x    = median(luecken_stats$detection_rate),
           y    = 0.5,
           label = paste0("Median\n", round(median(luecken_stats$detection_rate), 1), "%"),
           vjust = 0, hjust = -0.1, size = 3.2, color = COL_IGEL,
           fill = "white", label.size = 0) +
  scale_fill_gradient(low = "#D73027", high = COL_IGEL, guide = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 105),
                     expand = expansion(add = c(0, 5))) +
  labs(title    = "Detektionsrate pro Tier (aktive Periode)",
       subtitle = sprintf(
         "Anteil der Tage mit mindestens 1 Fix (Auswilderung → letztes Signal)\nn = %d Tiere | Median = %.1f%%",
         nrow(luecken_stats), median(luecken_stats$detection_rate)),
       x = "Detektionsrate (% der aktiven Tage)", y = NULL) +
  theme_classic(base_size = 13) +
  theme(panel.grid.major.x = element_line(color = "grey93", linewidth = 0.4))

ggsave(file.path(out_ordner, "07b_detection_gaps.png"), p2,
       width = 10, height = max(6, nrow(luecken_stats) * 0.35 + 2),
       dpi = 200, bg = "white")
cat("  ✓ 07b_detection_gaps.png\n")

# ══════════════════════════════════════════════════════════════════════════════
#  ANALYSE 3 — Station Count
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Analyse 3: Station Count ─────────────────────────────────\n")

if (!all(is.na(alle_fixes$station_count))) {

  sc_data <- alle_fixes[!is.na(station_count)]

  sc_stats <- sc_data[, .(
    N_fixes     = .N,
    pct_count2  = round(mean(station_count == 2) * 100, 1),
    pct_count3  = round(mean(station_count == 3) * 100, 1),
    pct_count4p = round(mean(station_count >= 4) * 100, 1),
    median_sc   = median(station_count)
  )]

  cat("\n  Station-Count-Verteilung:\n")
  print(sc_stats)
  cat(sprintf("\n  ► %.1f%% der Fixes mit Station Count = 2 (Mindestqualität)\n",
              sc_stats$pct_count2))
  cat(sprintf("  ► %.1f%% der Fixes mit Station Count ≥ 4 (hohe Qualität)\n",
              sc_stats$pct_count4p))

  if (sc_stats$pct_count2 > 50)
    cat("  [!] Mehr als die Hälfte der Fixes mit nur 2 Stationen.\n",
        "      → Homerange-Schätzungen konservativ interpretieren.\n")

  sc_freq <- sc_data[, .N, by = .(station_count = factor(station_count))]

  p3a <- ggplot(sc_freq, aes(x = station_count, y = N, fill = station_count)) +
    geom_col(alpha = 0.80, color = "grey30", linewidth = 0.5, width = 0.65) +
    geom_text(aes(label = paste0(round(N / sum(N) * 100, 1), "%")),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("2" = "#D73027", "3" = "#FDB863",
                                  "4" = "#74C476", "5" = "#2B8CBE",
                                  "6" = "#084594"),
                      guide = "none") +
    scale_y_continuous(labels = label_number(big.mark = "."),
                       expand = expansion(mult = c(0, 0.12))) +
    labs(title    = "Anzahl empfangender Stationen pro Fix",
         subtitle = sprintf("n = %d Fixes  |  Station Count = 2: minimale Multilateralisierung",
                            nrow(sc_data)),
         x = "Anzahl Stationen (Station Count)", y = "Anzahl Fixes") +
    theme_classic(base_size = 13) +
    theme(panel.grid.major.y = element_line(color = "grey93", linewidth = 0.4))

  p3b <- ggplot(sc_data, aes(x = dist_nearest_ant_m, y = factor(station_count))) +
    geom_jitter(height = 0.30, alpha = 0.25, size = 1.2, color = COL_IGEL) +
    stat_summary(aes(group = station_count),
                 fun = median, geom = "point",
                 shape = 18, size = 5, color = "firebrick") +
    scale_x_continuous(labels = label_number(suffix = " m")) +
    labs(title    = "Station Count vs. Abstand zur nächsten Antenne",
         subtitle = "Erwartung: weniger Stationen bei größerem Abstand  |  ◆ = Median",
         x = "Abstand zur nächsten Antenne (m)", y = "Station Count") +
    theme_classic(base_size = 13)

  p_sc <- p3a / p3b
  ggsave(file.path(out_ordner, "07c_station_count.png"), p_sc,
         width = 10, height = 10, dpi = 200, bg = "white")
  cat("  ✓ 07c_station_count.png\n")

} else {
  cat("  [INFO] Keine Station-Count-Daten in GeoPackages gefunden.\n")
  sc_stats <- NULL
}

# ══════════════════════════════════════════════════════════════════════════════
#  ANALYSE 4 — Detektionsbias nahe der Voliere
#
#  Hintergrund: In Block6_Dispersal starten alle Tiere bei ~150–200 m Distanz
#  zur Voliere. Das kann biologisch sein — oder ein Artefakt der Detektions-
#  reichweite (blinder Fleck im Zentrum des Antennennetzes).
#
#  Mechanismus: Tier nahe Voliere = symmetrische Signalausbreitung =
#  schlechte Hyperbel-Geometrie = kaum Multilateralisierungen möglich.
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Analyse 4: Detektionsbias nahe der Voliere ───────────────\n")

VOLIERE_LAT <- 52.398608
VOLIERE_LON <-  9.216834

vol_wgs <- st_sfc(st_point(c(VOLIERE_LON, VOLIERE_LAT)), crs = 4326)
vol_utm <- st_transform(vol_wgs, crs = 25832)
vol_xy  <- st_coordinates(vol_utm)

cat(sprintf("  Voliere UTM: X = %.1f, Y = %.1f\n", vol_xy[1,1], vol_xy[1,2]))

alle_fixes[, dist_voliere_m := sqrt((x_utm - vol_xy[1,1])^2 +
                                      (y_utm - vol_xy[1,2])^2)]

ANTENNEN$dist_vol_m <- round(sqrt((ANTENNEN$x_utm - vol_xy[1,1])^2 +
                                    (ANTENNEN$y_utm - vol_xy[1,2])^2))

cat("\n  Distanz der Antennen zur Voliere:\n")
print(ANTENNEN[, c("station_id", "dist_vol_m")])
naechste_ant_dist <- min(ANTENNEN$dist_vol_m)
cat(sprintf("\n  ► Nächste Antenne: %d m von der Voliere entfernt\n", naechste_ant_dist))
cat(sprintf("  ► Erwarteter blinder Fleck: ca. %d–%d m (30–50%% der Antennen-Distanz)\n",
            round(naechste_ant_dist * 0.3), round(naechste_ant_dist * 0.5)))

vol_stats <- alle_fixes[, .(
  N_fixes     = .N,
  P05_m       = round(quantile(dist_voliere_m, 0.05, na.rm = TRUE)),
  P10_m       = round(quantile(dist_voliere_m, 0.10, na.rm = TRUE)),
  Median_m    = round(median(dist_voliere_m, na.rm = TRUE)),
  Mean_m      = round(mean(dist_voliere_m, na.rm = TRUE)),
  P90_m       = round(quantile(dist_voliere_m, 0.90, na.rm = TRUE)),
  Max_m       = round(max(dist_voliere_m, na.rm = TRUE)),
  N_unter50m  = sum(dist_voliere_m <  50, na.rm = TRUE),
  N_unter100m = sum(dist_voliere_m < 100, na.rm = TRUE),
  N_unter150m = sum(dist_voliere_m < 150, na.rm = TRUE)
)]

cat("\n  Distanz zur Voliere (alle Fixes, aktive Periode):\n")
print(vol_stats)
cat(sprintf("\n  ► Fixes < 50 m:  %d (%.1f%%)\n",
            vol_stats$N_unter50m,  vol_stats$N_unter50m  / vol_stats$N_fixes * 100))
cat(sprintf("  ► Fixes < 100 m: %d (%.1f%%)\n",
            vol_stats$N_unter100m, vol_stats$N_unter100m / vol_stats$N_fixes * 100))
cat(sprintf("  ► Fixes < 150 m: %d (%.1f%%)\n",
            vol_stats$N_unter150m, vol_stats$N_unter150m / vol_stats$N_fixes * 100))

bias_bestaetigt <- (vol_stats$N_unter100m / vol_stats$N_fixes) < 0.02
bias_moeglich   <- (vol_stats$N_unter100m / vol_stats$N_fixes) < 0.10

if (bias_bestaetigt) {
  cat("\n  [!] DETEKTIONSBIAS BESTÄTIGT: < 2% der Fixes < 100 m → blinder Fleck.\n")
  cat("      → Block6-Distanzen spiegeln teilweise den Detektionshorizont wider.\n")
} else if (bias_moeglich) {
  cat(sprintf("\n  [~] LEICHTER DETEKTIONSBIAS möglich: %.1f%% der Fixes < 100 m.\n",
              vol_stats$N_unter100m / vol_stats$N_fixes * 100))
} else {
  cat(sprintf("\n  [OK] Kein offensichtlicher Detektionsbias nahe der Voliere.\n"))
}

breaks_vol <- c(0, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500)

p4a <- ggplot(alle_fixes, aes(x = dist_voliere_m)) +
  annotate("rect", xmin = 0, xmax = naechste_ant_dist * 0.4,
           ymin = 0, ymax = Inf, fill = "#D73027", alpha = 0.08) +
  annotate("text", x = naechste_ant_dist * 0.2, y = Inf,
           vjust = 1.6, size = 3.0, color = "#D73027", fontface = "italic",
           label = "Mögl.\nBlindzone") +
  geom_histogram(breaks = breaks_vol,
                 fill = COL_IGEL, color = "white", alpha = 0.80, linewidth = 0.4) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  annotate("text", x = 50, y = Inf, vjust = 1.5, hjust = -0.1,
           label = "Volieren-\nRadius (50 m)", size = 3.0, color = "grey40") +
  geom_vline(xintercept = naechste_ant_dist, linetype = "dotted",
             color = "firebrick", linewidth = 0.9) +
  annotate("text", x = naechste_ant_dist, y = Inf, vjust = 1.5, hjust = 1.05,
           label = sprintf("Nächste Ant.\n(%d m)", naechste_ant_dist),
           size = 3.0, color = "firebrick") +
  geom_vline(xintercept = vol_stats$Median_m, linetype = "solid",
             color = COL_IGEL, linewidth = 1.0, alpha = 0.7) +
  annotate("text", x = vol_stats$Median_m, y = Inf, vjust = 3.0, hjust = -0.1,
           label = sprintf("Median\n%d m", vol_stats$Median_m),
           size = 3.0, color = COL_IGEL) +
  scale_x_continuous(labels = label_number(suffix = " m"),
                     breaks = c(0, 50, 100, 150, 200, 250, 300, 400),
                     expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Verteilung der Fixes nach Distanz zur Auswilderungsvoliere",
       subtitle = sprintf(
         "n = %d Fixes von %d Tieren  |  Fixes < 100 m: %d (%.1f%%)  |  Rote Zone = mögl. Blindbereich",
         nrow(alle_fixes), n_distinct(alle_fixes$tier_label),
         vol_stats$N_unter100m, vol_stats$N_unter100m / vol_stats$N_fixes * 100),
       x = "Distanz zur Voliere (m)", y = "Anzahl Fixes") +
  theme_classic(base_size = 13) +
  theme(panel.grid.major.y = element_line(color = "grey93", linewidth = 0.4))

p4b <- ggplot(alle_fixes, aes(x = dist_voliere_m)) +
  annotate("rect", xmin = 0, xmax = naechste_ant_dist * 0.4,
           ymin = 0, ymax = 1, fill = "#D73027", alpha = 0.07) +
  stat_ecdf(color = COL_IGEL, linewidth = 1.3) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_vline(xintercept = vol_stats$P05_m, linetype = "dotted",
             color = "#9E3030", linewidth = 0.9) +
  annotate("text", x = vol_stats$P05_m, y = 0.05,
           hjust = -0.1, size = 3.0, color = "#9E3030",
           label = sprintf("5. Pz. = %d m\n(Mindest-Detektionsdistanz)", vol_stats$P05_m)) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "#9E3030",
             linewidth = 0.5, alpha = 0.6) +
  scale_x_continuous(labels = label_number(suffix = " m"),
                     breaks = c(0, 50, 100, 150, 200, 250, 300, 400)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(title    = "Kumulativverteilung: Distanz der Fixes zur Voliere",
       subtitle = sprintf(
         "5. Pz. = %d m → Mindestens %d m von der Voliere detektierbar",
         vol_stats$P05_m, vol_stats$P05_m),
       x = "Distanz zur Voliere (m)", y = "Kumulativer Anteil der Fixes") +
  theme_classic(base_size = 13)

p_vol <- p4a / p4b
ggsave(file.path(out_ordner, "07d_volieren_detektionsbias.png"), p_vol,
       width = 11, height = 11, dpi = 200, bg = "white")
cat("  ✓ 07d_volieren_detektionsbias.png\n")

# ══════════════════════════════════════════════════════════════════════════════
#  EXCEL-EXPORT
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Excel-Export ─────────────────────────────────────────────\n")

wb <- createWorkbook()

s_hdr  <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white",
                       fgFill = "#2E5E34", textDecoration = "bold",
                       halign = "center", border = "Bottom",
                       borderColour = "#BFBFBF")
s_body <- createStyle(fontName = "Arial", fontSize = 10)
s_alt  <- createStyle(fontName = "Arial", fontSize = 10, fgFill = "#F2F9F2")
s_warn <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "#991B1B",
                       textDecoration = "bold")

add_sheet_styled <- function(wb, name, data, warn_rows = NULL) {
  addWorksheet(wb, name)
  writeData(wb, name, as.data.frame(data), startRow = 1)
  nr <- nrow(data); nc <- ncol(data)
  addStyle(wb, name, s_hdr,  rows = 1,         cols = seq_len(nc), gridExpand = TRUE)
  addStyle(wb, name, s_body, rows = 2:(nr + 1), cols = seq_len(nc), gridExpand = TRUE)
  for (r in seq(2, nr + 1, by = 2))
    addStyle(wb, name, s_alt, rows = r, cols = seq_len(nc), gridExpand = TRUE, stack = TRUE)
  if (!is.null(warn_rows))
    addStyle(wb, name, s_warn, rows = warn_rows + 1, cols = seq_len(nc),
             gridExpand = TRUE, stack = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(nc), widths = "auto")
}

summary_rows <- data.frame(
  Analyse = c(
    "1. Detektionsdistanz — Median",
    "1. Detektionsdistanz — 90. Pz.",
    "1. Detektionsdistanz — 95. Pz.",
    "1. Detektionsdistanz — Maximum",
    "2. Detektionsrate — Median (%)",
    "2. Detektionsrate — Minimum (%)",
    "2. Detektionsrate — Maximum (%)",
    "3. Station Count = 2 (%)",
    "3. Station Count >= 4 (%)",
    "4. Nächste Antenne zur Voliere (m)",
    "4. 5. Pz. Fix-Distanz zur Voliere (m)",
    "4. Fixes < 100 m zur Voliere (%)"
  ),
  Wert = c(
    paste0(dist_stats$Median_m, " m"),
    paste0(dist_stats$P90_m,    " m"),
    paste0(dist_stats$P95_m,    " m"),
    paste0(dist_stats$Max_m,    " m"),
    paste0(median(luecken_stats$detection_rate, na.rm = TRUE), "%"),
    paste0(min(luecken_stats$detection_rate,    na.rm = TRUE), "%"),
    paste0(max(luecken_stats$detection_rate,    na.rm = TRUE), "%"),
    if (!is.null(sc_stats)) paste0(sc_stats$pct_count2, "%") else "n. v.",
    if (!is.null(sc_stats)) paste0(sc_stats$pct_count4p, "%") else "n. v.",
    paste0(naechste_ant_dist, " m"),
    paste0(vol_stats$P05_m, " m"),
    sprintf("%.1f%%", vol_stats$N_unter100m / vol_stats$N_fixes * 100)
  ),
  Interpretation = c(
    "Typischer Abstand eines Fixes zur nächsten Antenne",
    "90% der Fixes liegen innerhalb dieses Radius",
    "95% der Fixes liegen innerhalb dieses Radius",
    "Maximale detektierte Entfernung",
    "Anteil der aktiven Tage mit mindestens 1 Fix (Median aller Tiere)",
    "Tier mit den wenigsten Detektionen",
    "Tier mit den meisten Detektionen",
    "Fixes mit minimalem Stationsnachweis — geringste Positionsgenauigkeit",
    "Fixes mit hoher Multilateralisierungsqualität",
    "Geometrische Basis für Blindzonenabschätzung",
    "Unterste 5% der Distanzen = Mindest-Detektionsdistanz zur Voliere",
    "Anteil der Fixes im potentiellen Blindbereich"
  ),
  stringsAsFactors = FALSE
)
add_sheet_styled(wb, "00_Zusammenfassung", summary_rows)
add_sheet_styled(wb, "01_Detektionsdist_Gesamt",  dist_stats)
add_sheet_styled(wb, "01_Detektionsdist_ProTier", dist_pro_tier,
                 warn_rows = which(dist_pro_tier$P90_m > 300))

fixes_summary <- alle_fixes[, .(
  N_fixes       = .N,
  dist_median_m = round(median(dist_nearest_ant_m, na.rm = TRUE), 1),
  dist_p90_m    = round(quantile(dist_nearest_ant_m, 0.90, na.rm = TRUE), 1),
  dist_max_m    = round(max(dist_nearest_ant_m, na.rm = TRUE), 1),
  sc_median     = if (!all(is.na(station_count)))
                    round(median(station_count, na.rm = TRUE), 1) else NA_real_
), by = .(tier_label, date)][order(tier_label, date)]
add_sheet_styled(wb, "01_Tageszusammenfassung", fixes_summary)
add_sheet_styled(wb, "02_Detektionsluecken_ProTier", luecken_stats,
                 warn_rows = which(luecken_stats$detection_rate < 50))

if (!is.null(sc_stats)) add_sheet_styled(wb, "03_StationCount_Gesamt", sc_stats)

vol_export <- cbind(as.data.frame(vol_stats),
                    Naechste_Ant_dist_m   = naechste_ant_dist,
                    Blindzone_Abschaetz_m = round(naechste_ant_dist * 0.4))
add_sheet_styled(wb, "04_Volieren_Bias", vol_export)
add_sheet_styled(wb, "04_Antennen_Abstand_Voliere",
                 ANTENNEN[, c("station_id", "x_utm", "y_utm", "dist_vol_m")])

ant_export <- merge(ant_wgs84, ANTENNEN[, c("station_id", "x_utm", "y_utm")], by = "station_id")
add_sheet_styled(wb, "Antennen_Koordinaten", ant_export)

saveWorkbook(wb, file.path(out_ordner, "07_validierung_rohdaten.xlsx"), overwrite = TRUE)
cat("  ✓ 07_validierung_rohdaten.xlsx\n")

# ══════════════════════════════════════════════════════════════════════════════
#  5. WORD METHODENBERICHT
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Word Methodenbericht ─────────────────────────────────────\n")

if (!requireNamespace("officer",   quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")

if (requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)) {
  library(officer)
  library(flextable)

  make_ft <- function(df, hl_rows = NULL, warn_rows = NULL) {
    brd <- fp_border(color = "#BFBFBF", width = 0.5)
    ft  <- flextable(df) |>
      font(fontname = "Arial", part = "all") |>
      fontsize(size = 10, part = "all") |>
      bold(part = "header") |>
      color(part = "header", color = "white") |>
      bg(part = "header", bg = "#2E5E34") |>
      hline(border = brd, part = "all") |>
      vline(border = brd, part = "all") |>
      hline_top(border = brd, part = "header") |>
      bg(part = "body", bg = "white") |>
      set_table_properties(layout = "autofit")
    if (!is.null(hl_rows) && length(hl_rows) > 0)
      ft <- bg(ft, i = hl_rows, bg = "#E8F4EA", part = "body")
    if (!is.null(warn_rows) && length(warn_rows) > 0)
      ft <- bg(ft, i = warn_rows, bg = "#FEE2E2", part = "body") |>
              color(i = warn_rows, color = "#991B1B", part = "body")
    ft
  }

  doc_m <- read_docx()

  doc_m <- doc_m |>
    body_add_par("Block 7: Detektionsreichweiten-Validierung — Methodenbericht",
                 style = "heading 1") |>
    body_add_par("Igelbesenderung Wildtierstation Sachsenhagen — TiHo Hannover | Natalie Steiner",
                 style = "Normal") |>
    body_add_par(paste0("Erstellt: ", format(Sys.Date(), "%d. %B %Y")), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("1. Fragestellung und Relevanz", style = "heading 2") |>
    body_add_par(paste0(
      "Block 7 beantwortet die methodische Kernfrage: Werden alle räumlichen Aufenthalte ",
      "der Igel durch das tRackIT-System zuverlässig detektiert, oder gibt es systematische ",
      "Lücken aufgrund begrenzter Senderreichweite, schlechter Antennengeometrie oder ",
      "eines blinden Flecks nahe der Auswilderungsvoliere? ",
      "Diese Validierung ist Voraussetzung für die korrekte Interpretation von ",
      "Block 4 (Kernel Home Range) und Block 6 (Dispersal-Analyse). ",
      "Insbesondere der Befund aus Block 6, dass Tiere bereits am Auswilderungstag ",
      "in ~150–200 m Entfernung detektiert werden, könnte biologisch (Meidungsverhalten) ",
      "oder technisch (blinder Fleck im Antennennetz) bedingt sein."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("2. Datenbasis", style = "heading 2") |>
    body_add_par(paste0(
      "Grundlage sind alle VHF-Multilateralisierungen (GeoPackages, .gpkg) aus dem ",
      "Ordner data/kernel_files/ (inkl. Unterordner, recursive = TRUE). ",
      "Jeder Fix repräsentiert eine Positionsschätzung des tRackIT-Systems basierend auf ",
      "Signallaufzeitdifferenzen zwischen mindestens 2 Antennen (Station Count >= 2). ",
      "Analysiert werden ausschließlich Fixes aus der aktiven Tracking-Periode ",
      "(Auswilderungsdatum bis letztes registriertes Signal gemäß Metadatentabelle)."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("3. Antenneninfrastruktur (tRackIT)", style = "heading 2") |>
    body_add_par(paste0(
      "Das tRackIT-System umfasst 6 stationäre Empfangsantennen (A1–A6), ",
      "die das Untersuchungsgelände umschließen. Die Koordinaten wurden im ",
      "Siebenschläfer-Projekt (Steiner N., gleiche Infrastruktur) vermessen. ",
      "Koordinatensystem: WGS84 für Eingabe, UTM Zone 32N (EPSG 25832) für Distanzberechnungen."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("Tabelle 1: Antennenabstand zur Auswilderungsvoliere", style = "Normal") |>
    body_add_flextable(make_ft(
      setNames(ANTENNEN[, c("station_id", "dist_vol_m")],
               c("Antenne", "Distanz zur Voliere [m]"))
    )) |>
    body_add_par("", style = "Normal") |>

    body_add_par("4. Methoden der vier Analysen", style = "heading 2") |>

    body_add_par("4.1 Detektionsdistanz (Analyse 1)", style = "heading 3") |>
    body_add_par(paste0(
      "Für jeden Fix wird die euklidische Distanz zur nächstgelegenen der 6 Antennen berechnet (UTM 32N). ",
      "Die Verteilung dieser Distanzen zeigt den effektiven Abdeckungsbereich des Systems. ",
      "Beurteilungskriterium: Wenn das 95. Perzentil unter dem mittleren Antennen-Antennen-Abstand liegt, ",
      "war das Netz für das Untersuchungsgebiet ausreichend dicht. ",
      "Kennzahlen: Median, 75., 90., 95. Perzentile und Maximum [m], gesamt und pro Tier."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("4.2 Detektionslücken (Analyse 2)", style = "heading 3") |>
    body_add_par(paste0(
      "Pro Tier wird ein vollständiges Datumsraster von Auswilderung bis letztem Signal erzeugt. ",
      "Für jeden Tag wird geprüft, ob mindestens 1 Fix vorliegt. ",
      "Detektionsrate = Anteil detektierter Tage x 100 %. ",
      "Tiere mit < 50 % Detektionsrate werden als kritisch markiert (Signalverlust oder Randaufenthalt). ",
      "Visualisierung: Lollipop-Chart sortiert nach Detektionsrate, Medianlinie eingezeichnet."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("4.3 Station Count — Multilateralisierungsqualität (Analyse 3)", style = "heading 3") |>
    body_add_par(paste0(
      "Der Station Count gibt an, wie viele Antennen das Signal gleichzeitig empfangen haben. ",
      "Station Count 2: Minimalfall — Hyperbel-Schnittpunkt ergibt 2 mögliche Positionen, ",
      "Software wählt die wahrscheinlichere aus. Positionsgenauigkeit ca. 20–50 m. ",
      "Station Count >= 4: Überbestimmtes System, Positionsgenauigkeit ca. 5–15 m. ",
      "Plots: (a) Häufigkeitsverteilung nach Station Count, ",
      "(b) Station Count vs. Abstand zur nächsten Antenne (Erwartung: Inverse Korrelation)."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("4.4 Detektionsbias nahe der Voliere (Analyse 4)", style = "heading 3") |>
    body_add_par(paste0(
      "Hypothese: Das tRackIT-System hat einen blinden Fleck direkt bei der Auswilderungsvoliere. ",
      "Mechanismus: Ein Tier im geometrischen Zentrum des Antennennetzes sendet gleichmaessig ",
      "in alle Richtungen — die Signallaufzeitdifferenzen zwischen den Antennen sind minimal, ",
      "die Hyperbel-Geometrie ist schlecht konditioniert, und es entstehen kaum Multilateralisierungen. ",
      "Test: Die Verteilung der Fix-Distanzen zur Voliere wird auf einen scharfen cutoff bei < 100 m geprueft. ",
      "Schwellenwert: < 2 % der Fixes innerhalb 100 m → Detektionsbias bestätigt. ",
      "Erwarteter Blindfleck-Radius: 30–50 % der Distanz zur nächsten Antenne."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("5. Limitationen", style = "heading 2") |>
    body_add_par(paste0(
      "Die Detektionsreichweite ist abhängig von Vegetation (Laubdichte, saisonale Variation), ",
      "Topographie, Wetterbedingungen (Regen dämpft VHF-Signale) und Senderzustand (Batteriestand). ",
      "Detektionslücken können technisch (Senderabfall, Batterieerschöpfung) oder biologisch ",
      "(Tier außerhalb Reichweite) bedingt sein — eine Trennung ist anhand der VHF-Daten nicht möglich. ",
      "Der Detektionsbias nahe der Voliere kann nicht vollständig von echtem Meidungsverhalten ",
      "unterschieden werden. Für eine vollständige Validierung wäre ein Testsender mit bekanntem ",
      "Aufenthaltsort optimal (Testdaten liegen für dieses Projekt nicht vor)."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("6. Referenz", style = "heading 2") |>
    body_add_par(paste0(
      "Gottwald J, Wild T, Zeidler S, Friess N, Meinecke I, Warth P, Nauss T, Farwig N (2021). ",
      "tRackIT OS: An open-source software system for fine-grained analysis of animal movement ",
      "from radio-telemetry data. Methods in Ecology and Evolution 12(10):1820–1832."
    ), style = "Normal") |>
    body_add_par(paste0(
      "Steiner N (2025). Igelbesenderung Wildtierstation Sachsenhagen — ",
      "Methodendokumentation Block 7. TiHo Hannover."
    ), style = "Normal")

  out_m <- file.path(out_ordner, "Block7_Detektionsvalidierung_Methodenbericht.docx")
  print(doc_m, target = out_m)
  cat("  ✓ Block7_Detektionsvalidierung_Methodenbericht.docx\n")

} else {
  cat("  [WARN] officer/flextable nicht verfügbar — kein Word-Methodenbericht.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
#  6. WORD ERGEBNISBERICHT (dynamisch — nutzt berechnete Variablen)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Word Ergebnisbericht ─────────────────────────────────────\n")

if (requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)) {

  make_ft2 <- function(df, hl_rows = NULL, warn_rows = NULL) {
    brd <- fp_border(color = "#BFBFBF", width = 0.5)
    ft  <- flextable(df) |>
      font(fontname = "Arial", part = "all") |>
      fontsize(size = 10, part = "all") |>
      bold(part = "header") |>
      color(part = "header", color = "white") |>
      bg(part = "header", bg = "#2E5E34") |>
      hline(border = brd, part = "all") |>
      vline(border = brd, part = "all") |>
      hline_top(border = brd, part = "header") |>
      bg(part = "body", bg = "white") |>
      set_table_properties(layout = "autofit")
    if (!is.null(hl_rows) && length(hl_rows) > 0)
      ft <- bg(ft, i = hl_rows, bg = "#E8F4EA", part = "body")
    if (!is.null(warn_rows) && length(warn_rows) > 0)
      ft <- bg(ft, i = warn_rows, bg = "#FEE2E2", part = "body") |>
              color(i = warn_rows, color = "#991B1B", part = "body")
    ft
  }

  doc_e <- read_docx()

  doc_e <- doc_e |>
    body_add_par("Block 7: Detektionsreichweiten-Validierung — Ergebnisbericht",
                 style = "heading 1") |>
    body_add_par("Igelbesenderung Wildtierstation Sachsenhagen — TiHo Hannover | Natalie Steiner",
                 style = "Normal") |>
    body_add_par(paste0("Erstellt: ", format(Sys.Date(), "%d. %B %Y")), style = "Normal") |>
    body_add_par("", style = "Normal")

  # Überblicktabelle
  uebersicht_df <- data.frame(
    Analyse  = c("Tiere analysiert", "Fixes analysiert (aktive Periode)",
                 "Detektionsdistanz — Median", "Detektionsdistanz — 95. Pz.",
                 "Detektionsrate — Median", "Detektionsrate — Minimum",
                 "Nächste Antenne zur Voliere", "5. Pz. Fix-Distanz zur Voliere",
                 "Fixes < 100 m zur Voliere"),
    Ergebnis = c(
      sprintf("%d Tiere", n_distinct(alle_fixes$tier_label)),
      sprintf("%d", nrow(alle_fixes)),
      sprintf("%d m", dist_stats$Median_m),
      sprintf("%d m", dist_stats$P95_m),
      sprintf("%.1f%%", median(luecken_stats$detection_rate)),
      sprintf("%.1f%%", min(luecken_stats$detection_rate)),
      sprintf("%d m", naechste_ant_dist),
      sprintf("%d m", vol_stats$P05_m),
      sprintf("%d (%.1f%%)", vol_stats$N_unter100m,
              vol_stats$N_unter100m / vol_stats$N_fixes * 100)
    ), stringsAsFactors = FALSE
  )
  doc_e <- doc_e |>
    body_add_par("Tabelle 1: Überblick Kennzahlen", style = "Normal") |>
    body_add_flextable(make_ft2(uebersicht_df)) |>
    body_add_par("", style = "Normal")

  # 1. Detektionsdistanz
  ok_p95 <- dist_stats$P95_m < 300
  doc_e <- doc_e |>
    body_add_par("1. Detektionsdistanz zur nächsten Antenne", style = "heading 2") |>
    body_add_par(sprintf(paste0(
      "Die Distanz jedes Fixes zur nächstgelegenen Antenne betrug im Median %d m ",
      "(90. Pz.: %d m, 95. Pz.: %d m, Maximum: %d m). "),
      dist_stats$Median_m, dist_stats$P90_m, dist_stats$P95_m, dist_stats$Max_m),
      style = "Normal") |>
    body_add_par(
      if (ok_p95)
        sprintf(paste0(
          "95%% aller Fixes lagen innerhalb von %d m zur nächsten Antenne — ",
          "das Antennennetz war für das Untersuchungsgebiet ausreichend dicht. ",
          "Die Detektionsreichweite stellte keine wesentliche Einschränkung für ",
          "die Homerange-Schätzungen (Block 4) dar."), dist_stats$P95_m)
      else
        sprintf(paste0(
          "Das 95. Perzentil (%d m) überschreitet den 300-m-Richtwert. ",
          "5%% der Fixes liegen möglicherweise in Bereichen mit reduzierter Genauigkeit. ",
          "Dies sollte als Limitation im Manuskript genannt werden."), dist_stats$P95_m),
      style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("Tabelle 2: Detektionsdistanz pro Tier", style = "Normal") |>
    body_add_flextable(make_ft2(
      setNames(dist_pro_tier, c("Tier", "N Fixes", "Median [m]", "90. Pz. [m]", "Max [m]")),
      warn_rows = which(dist_pro_tier$P90_m > 300)
    )) |>
    body_add_par("Rot markiert: 90. Pz. > 300 m.", style = "Normal") |>
    body_add_par("", style = "Normal")

  # 2. Detektionslücken
  median_rate <- median(luecken_stats$detection_rate)
  min_rate    <- min(luecken_stats$detection_rate)
  n_schlecht  <- sum(luecken_stats$detection_rate < 50)

  doc_e <- doc_e |>
    body_add_par("2. Detektionslücken — Tage ohne Signal", style = "heading 2") |>
    body_add_par(sprintf(paste0(
      "Die Detektionsrate betrug im Median %.1f%% (Min: %.1f%%, Max: %.1f%%). ",
      "%d Tier(e) erreichten eine Detektionsrate < 50%%. "),
      median_rate, min_rate, max(luecken_stats$detection_rate), n_schlecht),
      style = "Normal") |>
    body_add_par(
      if (n_schlecht == 0)
        "Alle Tiere wurden an der Mehrheit ihrer aktiven Tage detektiert — die Datenlage ist gut."
      else
        sprintf(paste0(
          "%d Tier(e) mit < 50%% Detektionsrate sind in Tabelle 3 rot markiert. ",
          "Ob die Lücken biologisch (Tier außerhalb Reichweite) oder technisch ",
          "(Senderabfall) bedingt sind, lässt sich nicht eindeutig klären."), n_schlecht),
      style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("Tabelle 3: Detektionsrate pro Tier", style = "Normal") |>
    body_add_flextable(make_ft2(
      setNames(luecken_stats,
               c("Tier", "Tage gesamt", "Tage detektiert", "Tage fehlend", "Detektionsrate [%]")),
      warn_rows = which(luecken_stats$detection_rate < 50)
    )) |>
    body_add_par("Rot markiert: Detektionsrate < 50%.", style = "Normal") |>
    body_add_par("", style = "Normal")

  # 3. Station Count
  if (!is.null(sc_stats)) {
    doc_e <- doc_e |>
      body_add_par("3. Station Count — Multilateralisierungsqualität", style = "heading 2") |>
      body_add_par(sprintf(paste0(
        "Von %d Fixes: %.1f%% mit Station Count 2 (Mindestqualität), ",
        "%.1f%% mit Count 3, %.1f%% mit Count >= 4 (hohe Qualität). Median Count: %.0f."),
        sc_stats$N_fixes, sc_stats$pct_count2,
        sc_stats$pct_count3, sc_stats$pct_count4p, sc_stats$median_sc),
        style = "Normal") |>
      body_add_par(
        if (sc_stats$pct_count2 > 50)
          paste0(
            "Der hohe Anteil von Station Count 2 (> 50%) bedeutet eingeschränkte Positionsgenauigkeit ",
            "für den Großteil der Fixes. Homerange-Schätzungen (Block 4) sind als ",
            "konservative Mindestabschätzungen zu interpretieren.")
        else
          "Die Mehrheit der Fixes entstand mit Station Count >= 3 — zufriedenstellende Qualität.",
        style = "Normal") |>
      body_add_par("", style = "Normal")
  }

  # 4. Volieren-Detektionsbias
  pct_unter100 <- vol_stats$N_unter100m / vol_stats$N_fixes * 100

  doc_e <- doc_e |>
    body_add_par("4. Detektionsbias nahe der Auswilderungsvoliere", style = "heading 2") |>
    body_add_par(sprintf(paste0(
      "Die nächste Empfangsantenne befand sich %d m von der Auswilderungsvoliere entfernt. ",
      "Erwarteter Blindfleck-Radius: ca. %d–%d m. ",
      "Tatsächlich lagen nur %d Fixes (%.1f%%) innerhalb von 100 m zur Voliere. ",
      "Das 5. Perzentil der Distanzverteilung betrug %d m."),
      naechste_ant_dist,
      round(naechste_ant_dist * 0.3), round(naechste_ant_dist * 0.5),
      vol_stats$N_unter100m, pct_unter100,
      vol_stats$P05_m),
      style = "Normal") |>
    body_add_par(
      if (bias_bestaetigt)
        sprintf(paste0(
          "DETEKTIONSBIAS BESTÄTIGT: < 2%% der Fixes liegen innerhalb 100 m zur Voliere. ",
          "Das tRackIT-System hat einen blinden Fleck von ca. %d m um die Voliere. ",
          "Dies erklärt, warum in Block 6 die mittlere Tages-Distanz bereits am Auswilderungstag ",
          "bei ~150–200 m liegt: Aufenthalte nahe der Voliere sind systematisch unterrepräsentiert. ",
          "Block-6-Distanzwerte sind als Mindestabschätzungen zu verstehen. ",
          "Diese Limitation ist im Manuskript explizit zu benennen."), vol_stats$P05_m)
      else if (bias_moeglich)
        sprintf(paste0(
          "Leichter Detektionsbias möglich: %.1f%% der Fixes < 100 m. ",
          "Im Manuskript als mögliche Limitation erwähnen."), pct_unter100)
      else
        "Kein offensichtlicher Detektionsbias nahe der Voliere.",
      style = "Normal") |>
    body_add_par("", style = "Normal")

  # 5. Empfehlungen
  doc_e <- doc_e |>
    body_add_par("5. Empfehlungen für Publikation und Folgestudien", style = "heading 2") |>
    body_add_par(paste0(
      "(1) Detektionsreichweite (95. Pz.) und Station-Count-Verteilung als Supplementary-Table einbinden. ",
      "(2) Den Detektionsbias nahe der Voliere (blinder Fleck ~ 5. Pz. der Fix-Distanz) ",
      "als explizite Limitation in der Diskussion adressieren — ",
      "Block-6-Distanzwerte repräsentieren Mindestwerte des tatsächlichen Aufenthalts. ",
      "(3) Tiere mit Detektionsrate < 50% in Homerange-Analysen mit Hinweis interpretieren. ",
      "(4) Für Folgestudien: Eine zusätzliche Antenne nahe der Voliere (~50–80 m) ",
      "würde den blinden Fleck deutlich reduzieren."
    ), style = "Normal")

  out_e <- file.path(out_ordner, "Block7_Detektionsvalidierung_Ergebnisbericht.docx")
  print(doc_e, target = out_e)
  cat("  ✓ Block7_Detektionsvalidierung_Ergebnisbericht.docx\n")

} else {
  cat("  [WARN] officer/flextable nicht verfügbar — kein Word-Ergebnisbericht.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
#  ZUSAMMENFASSUNG
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("  ZUSAMMENFASSUNG — Block 7 Detektionsreichweiten-Validierung\n")
cat("══════════════════════════════════════════════════════════════\n\n")
cat(sprintf("  Tiere:              %d\n", n_distinct(alle_fixes$tier_label)))
cat(sprintf("  Fixes analysiert:   %d\n", nrow(alle_fixes)))
cat(sprintf("  Antennen (tRackIT): %d Stationen\n", nrow(ANTENNEN)))
cat(sprintf("\n  Analyse 1 — Detektionsdistanz:\n"))
cat(sprintf("    Median:  %d m  |  90. Pz.: %d m  |  95. Pz.: %d m  |  Max: %d m\n",
            dist_stats$Median_m, dist_stats$P90_m, dist_stats$P95_m, dist_stats$Max_m))
cat(sprintf("\n  Analyse 2 — Detektionsrate:\n"))
cat(sprintf("    Median: %.1f%%  |  Min: %.1f%%  |  Max: %.1f%%\n",
            median(luecken_stats$detection_rate),
            min(luecken_stats$detection_rate),
            max(luecken_stats$detection_rate)))
if (!is.null(sc_stats)) {
  cat(sprintf("\n  Analyse 3 — Station Count 2: %.1f%%  |  Count >= 4: %.1f%%\n",
              sc_stats$pct_count2, sc_stats$pct_count4p))
}
cat(sprintf("\n  Analyse 4 — Volieren-Bias:\n"))
cat(sprintf("    Nächste Antenne: %d m  |  5. Pz. zur Voliere: %d m  |  Fixes < 100 m: %.1f%%\n",
            naechste_ant_dist, vol_stats$P05_m,
            vol_stats$N_unter100m / vol_stats$N_fixes * 100))
cat(sprintf("    Status: %s\n",
            if (bias_bestaetigt) "BIAS BESTÄTIGT" else if (bias_moeglich) "BIAS MÖGLICH" else "KEIN BIAS"))
cat("\n══════════════════════════════════════════════════════════════\n")
cat(sprintf("✓ Block7 fertig. Output: %s\n", out_ordner))
