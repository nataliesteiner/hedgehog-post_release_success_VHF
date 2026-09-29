# ==============================================================
# Block 4 — space use: fix loading & spatial metrics
# ==============================================================
# Project:  Hedgehog VHF telemetry, Sachsenhagen, Lower Saxony
# Author:   Natalie Steiner, TiHo Hannover
#
# This script:
#   1. Loads multilateration fixes (*.gpkg) and metadata
#   2. Computes spatial metrics per night and per animal:
#        - maximum night distance (m)
#        - radius of gyration (RoG, m) — fix-count-independent
#        - night-centroid distance to the release enclosure (m)
#        - return rate (<= 100 m to the enclosure)
#   3. Saves kernel_ergebnisse.rds for Block4b
#
# IMPORTANT: no standard KDE (adehabitatHR) — home-range sizes
#            come exclusively from Block4b (aKDE, ctmm package).
#            Pairwise overlap is also in Block4b (aKDE polygons).
#
# Order:

# ── 0. Pakete ──────────────────────────────────────────────────
pakete <- c("sf", "data.table", "ggplot2", "patchwork", "scales",
            "openxlsx", "readxl", "lubridate", "viridis")
fehlend <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(fehlend) > 0) install.packages(fehlend)

suppressPackageStartupMessages({
  library(sf); library(data.table); library(ggplot2)
  library(patchwork); library(scales); library(openxlsx)
  library(readxl); library(lubridate); library(viridis)
})
cat("Pakete geladen\n\n")

# ── 1. Einstellungen ───────────────────────────────────────────
projekt_root <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
gpkg_ordner  <- file.path(projekt_root, "data", "kernel_files")
meta_datei   <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
out_ordner   <- file.path(projekt_root, "output", "Block4_Kernel")
dir.create(out_ordner, showWarnings = FALSE, recursive = TRUE)

# Auswilderungsgehege (UTM 32N, EPSG:25832)
release_x  <- 514743.7
release_y  <- 5805363.8
ziel_epsg  <- 25832

# Rueckkehr-Schwellenwert
rueckkehr_m <- 100  # Distanz zum Gehege fuer Rueckkehr-Klassifikation

cat("Output:", out_ordner, "\n\n")

# ── 2. Metadaten laden ─────────────────────────────────────────
meta_raw <- as.data.table(read_excel(meta_datei))

meta <- data.table(
  igel         = trimws(meta_raw$individual),
  sex          = meta_raw$sex,
  release_date = as.Date(meta_raw$date_release),
  time_reha    = as.numeric(meta_raw$time_reha),
  diagnose     = meta_raw$diagnosis_main,
  tagging_date = as.Date(meta_raw$tagging_date),
  # Altersklasse: Einlieferungsgewicht < 300 g = Jungtier (Paper-Kriterium)
  alter        = ifelse(as.numeric(meta_raw$weigth_entry) < 300,
                        "Jungtier", "Adult")
)
meta <- meta[!is.na(igel) & igel != "" & !is.na(release_date)]

# Saison
meta[, saison_auswild := {
  m <- month(release_date)
  fcase(m %in% 3:5,  "Frühling",
        m %in% 6:8,  "Sommer",
        m %in% 9:11, "Herbst",
        default      = "Winter")
}]

# Diagnosegruppe (ohne Altersklasse — kommt aus Block1)
meta[, diag_gruppe := fcase(
  grepl("(?i)parasit",               diagnose), "Parasiten",
  grepl("(?i)trauma|injury|verletz", diagnose), "Trauma/Verletzung",
  grepl("(?i)blind|eye|auge|ocular", diagnose), "Augen-Pathologie",
  default = "Andere/Unbekannt"
)]

cat(sprintf("Metadaten: %d Igel geladen\n\n", nrow(meta)))

# ── 3. GPKG-Fixes laden ────────────────────────────────────────
# Identische Ladefunktion wie Block4d_TemporaleAnalyse.R
lade_fixes <- function(pfad) {
  tryCatch({
    layer    <- st_layers(pfad)$name[1]
    sf_all   <- st_read(pfad, layer = layer, quiet = TRUE)
    crs_pfad <- st_crs(sf_all)
    dat      <- as.data.table(sf_all)

    if (!"Night" %in% names(dat)) {
      message(sprintf("  [SKIP] %s: Keine Night-Spalte", basename(pfad)))
      return(NULL)
    }
    nacht <- dat[!is.na(Night) & Night > 0]
    if (nrow(nacht) < 15) {
      message(sprintf("  [SKIP] %s: nur %d Nacht-Fixes", basename(pfad), nrow(nacht)))
      return(NULL)
    }

    # Koordinaten extrahieren
    geom_col <- attr(sf_all, "sf_column")
    if (is.null(geom_col) || !geom_col %in% names(nacht))
      geom_col <- intersect(c("geometry","geom","the_geom"), names(nacht))[1]

    xy_ok <- FALSE
    if (!is.na(geom_col)) {
      tryCatch({
        sf_n <- st_as_sf(nacht, sf_column_name = geom_col, crs = crs_pfad)
        xy   <- as.data.table(st_coordinates(sf_n))
        if (!anyNA(xy$X)) {
          nacht[, x := xy$X]
          nacht[, y := xy$Y]
          xy_ok <- TRUE
        }
      }, error = function(e) NULL)
    }
    if (!xy_ok && geom_col %in% names(nacht)) {
      gl <- nacht[[geom_col]]
      nacht[, x := vapply(gl, function(g) {
        v <- tryCatch(as.numeric(g), error=function(e) c(NA,NA))
        if (length(v)>=1) v[1] else NA_real_}, numeric(1))]
      nacht[, y := vapply(gl, function(g) {
        v <- tryCatch(as.numeric(g), error=function(e) c(NA,NA))
        if (length(v)>=2) v[2] else NA_real_}, numeric(1))]
    }

    igel_col <- grep("^Individual[ .]Name$", names(nacht), value=TRUE)[1]
    date_col <- grep("^Date$",  names(nacht), value=TRUE)[1]
    time_col <- grep("^Time$",  names(nacht), value=TRUE)[1]

    if (is.na(igel_col)) {
      nacht[, igel := sub("^(Igel\\d+)_.*$", "\\1", basename(pfad))]
    } else {
      nacht[, igel := trimws(nacht[[igel_col]])]
    }
    nacht[, datum    := as.Date(if (!is.na(date_col)) nacht[[date_col]] else NA_character_)]
    nacht[, datetime := with_tz(
      as.POSIXct(paste(if (!is.na(date_col)) nacht[[date_col]] else "1970-01-01",
                       if (!is.na(time_col)) nacht[[time_col]] else "00:00:00"),
                 format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      "Europe/Berlin"
    )]
    nacht <- nacht[!is.na(x) & !is.na(y) & x != 0 & y != 0]
    if (nrow(nacht) < 15) return(NULL)
    nacht[, .(igel, datum, datetime, x, y, Night)]

  }, error = function(e) {
    message(sprintf("  [FEHLER] %s: %s", basename(pfad), conditionMessage(e)))
    NULL
  })
}

gpkg_dateien <- list.files(gpkg_ordner, pattern="\\.gpkg$",
                            full.names=TRUE, recursive=TRUE)
gpkg_dateien <- gpkg_dateien[!grepl("\\(\\d+\\)\\.gpkg$", gpkg_dateien)]

cat(sprintf("Lade GPKG-Fixes aus %d Dateien...\n", length(gpkg_dateien)))
alle_fixes_liste <- lapply(gpkg_dateien, lade_fixes)
alle_fixes <- rbindlist(Filter(Negate(is.null), alle_fixes_liste), fill=TRUE)

if (nrow(alle_fixes) == 0)
  stop("Keine GPKG-Fixes geladen! Pfad pruefen: ", gpkg_ordner)

alle_fixes <- merge(
  alle_fixes,
  meta[, .(igel, release_date, sex, time_reha, diagnose, diag_gruppe, saison_auswild)],
  by = "igel", all.x = TRUE
)
alle_fixes[, tage_seit := as.numeric(datum - release_date)]
alle_fixes <- alle_fixes[!is.na(tage_seit) & tage_seit >= 0]
setorder(alle_fixes, igel, datetime)

cat(sprintf("  %d Fixes, %d Igel geladen\n\n",
            nrow(alle_fixes), uniqueN(alle_fixes$igel)))

# ── 4. Raeumliche Metriken (KEIN KDE benoetigt) ────────────────
cat("=========================================\n")
cat("4. Raeumliche Metriken pro Nacht\n")
cat("=========================================\n\n")

dist_to_release <- function(x, y)
  sqrt((x - release_x)^2 + (y - release_y)^2)

radius_of_gyration <- function(x, y) {
  cx <- mean(x, na.rm=TRUE); cy <- mean(y, na.rm=TRUE)
  sqrt(mean((x-cx)^2 + (y-cy)^2, na.rm=TRUE))
}

nacht_metriken <- alle_fixes[, {
  xi <- x; yi <- y; n_f <- length(xi)
  cx  <- mean(xi, na.rm=TRUE); cy <- mean(yi, na.rm=TRUE)
  dist_zentroid     <- sqrt((xi-cx)^2 + (yi-cy)^2)
  dist_gehege_alle  <- dist_to_release(xi, yi)
  # Maximale Nachtdistanz: groesste Entfernung vom ersten Fix der Nacht
  # = wie weit entfernte sich das Tier maximal vom Startpunkt?
  # Robuster als erste-zu-letzte (erfasst auch Ausfluege mit Rueckkehr)
  max_dist <- if (n_f >= 2) {
    x0 <- xi[1]; y0 <- yi[1]
    max(sqrt((xi - x0)^2 + (yi - y0)^2), na.rm = TRUE)
  } else NA_real_
  list(
    n_fixes           = n_f,
    zentroid_x        = cx,
    zentroid_y        = cy,
    max_dist_m        = max_dist,
    dist_gehege_m     = dist_to_release(cx, cy),
    min_dist_gehege_m = min(dist_gehege_alle, na.rm=TRUE),
    rog_m             = radius_of_gyration(xi, yi),
    akt_radius_m      = as.numeric(quantile(dist_zentroid, 0.95, na.rm=TRUE)),
    rueckkehr         = (min(dist_gehege_alle, na.rm=TRUE) <= rueckkehr_m)
  )
}, by = .(igel, datum, tage_seit)]

# Zusammenfassung pro Tier
metriken_summary <- nacht_metriken[, .(
  n_naechte          = .N,
  n_fixes_gesamt     = sum(n_fixes),
  max_dist_median_m  = round(median(max_dist_m,        na.rm=TRUE), 0),
  max_dist_max_m     = round(max(max_dist_m,           na.rm=TRUE), 0),
  max_dist_mean_m    = round(mean(max_dist_m,          na.rm=TRUE), 0),
  rog_m              = round(median(rog_m,             na.rm=TRUE), 1),
  dist_gehege_med_m  = round(median(dist_gehege_m,     na.rm=TRUE), 0),
  dist_gehege_mean_m = round(mean(dist_gehege_m,       na.rm=TRUE), 0),
  min_dist_gehege_m  = round(min(min_dist_gehege_m,    na.rm=TRUE), 0),
  n_rueckkehr        = sum(rueckkehr, na.rm=TRUE),
  pct_rueckkehr      = round(100 * mean(rueckkehr, na.rm=TRUE), 1)
), by = igel]

metriken_summary <- merge(
  metriken_summary,
  meta[, .(igel, alter, sex, saison_auswild, diag_gruppe, time_reha)],
  by = "igel", all.x = TRUE
)

n_rueck <- sum(metriken_summary$n_rueckkehr > 0, na.rm=TRUE)

cat("Populationsmediane:\n")
cat(sprintf("  Max. Nachtdistanz:  %d m  (Range %d - %d m)\n",
    round(median(metriken_summary$max_dist_median_m, na.rm=TRUE)),
    round(min(metriken_summary$max_dist_median_m, na.rm=TRUE)),
    round(max(metriken_summary$max_dist_max_m, na.rm=TRUE))))
cat(sprintf("  Radius of Gyration: %d m  (Range %d - %d m)\n",
    round(median(metriken_summary$rog_m, na.rm=TRUE)),
    round(min(metriken_summary$rog_m, na.rm=TRUE)),
    round(max(metriken_summary$rog_m, na.rm=TRUE))))
cat(sprintf("  Zentroid-Distanz Gehege: %d m (Median)\n",
    round(median(metriken_summary$dist_gehege_med_m, na.rm=TRUE))))
cat(sprintf("  Min. Zentroid-Distanz Gehege: %d m (Min. ueber alle Tiere)\n",
    round(min(metriken_summary$min_dist_gehege_m, na.rm=TRUE))))
cat(sprintf("  Rueckkehr <= %dm: %d/%d Tiere (%.0f%%)\n\n",
    rueckkehr_m, n_rueck, nrow(metriken_summary),
    100*n_rueck/nrow(metriken_summary)))

# ── 5. Tieruebersicht (Metadaten-Zusammenfassung) ──────────────
igel_uebersicht <- nacht_metriken[, .(n_naechte=.N, n_fixes=sum(n_fixes)),
                                   by=igel]
igel_uebersicht <- merge(igel_uebersicht,
  meta[, .(igel, sex, saison_auswild, diag_gruppe, time_reha, release_date)],
  by="igel", all.x=TRUE)

cat(sprintf("Tiere mit Nacht-Fixes: %d\n", nrow(igel_uebersicht)))
print(igel_uebersicht[order(release_date), .(igel, n_naechte, n_fixes, sex,
                                              saison_auswild, diag_gruppe)])
cat("\n")

# ── 6. kernel_ergebnisse.rds speichern ─────────────────────────
# WIRD VON Block4b_Einzeltier_aKDE.R BENOETIGT
# Enthaelt: alle_fixes (GPS-Daten mit Metadaten)
#           kde_gesamt = Metadaten-Tabelle (OHNE KDE-Flaechen —
#                        diese kommen aus Block4b!)
cat("Speichere kernel_ergebnisse.rds...\n")

# kde_gesamt: nur Metadaten + Nacht/Fix-Zaehler (keine KDE-Flaechen)
# Block4b ergaenzt diese Tabelle um akde_95ha und akde_50ha
kde_gesamt_meta <- igel_uebersicht[, .(
  igel, n_naechte, n_fixes,
  sex, saison_auswild, diag_gruppe, time_reha,
  kde_95ha = NA_real_,   # Platzhalter — wird von Block4b befuellt
  kde_50ha = NA_real_    # Platzhalter — wird von Block4b befuellt
)]

saveRDS(
  list(
    alle_fixes       = alle_fixes,        # GPS-Fixes (fuer Block4b)
    kde_gesamt       = kde_gesamt_meta,   # Metadaten (fuer Block4b-Merge)
    nacht_metriken   = nacht_metriken,    # Pro-Nacht-Metriken
    metriken_summary = metriken_summary,  # Zusammenfassung pro Tier
    release_punkt    = data.frame(x=release_x, y=release_y)
  ),
  file.path(out_ordner, "kernel_ergebnisse.rds")
)
cat("  kernel_ergebnisse.rds gespeichert\n\n")

# ── 7. Excel-Uebersicht ────────────────────────────────────────
# Enthaelt nur Distanzmetriken (KEINE Homerange-Groessen —
# diese kommen nach Block4b aus akde_statistik.xlsx)
cat("Erstelle Block4_Kernel_Uebersicht.xlsx...\n")

wb <- createWorkbook()
hs   <- createStyle(fgFill="#2C3E6B", fontColour="white",
                    textDecoration="bold", halign="center", border="Bottom")
grau <- createStyle(fgFill="#F2F2F2")

schreibe_sheet <- function(wb, name, titel, daten) {
  addWorksheet(wb, name)
  writeData(wb, name, data.frame(V1=titel), startRow=1, colNames=FALSE)
  writeData(wb, name, as.data.frame(daten), startRow=3, colNames=TRUE)
  n <- ncol(daten)
  addStyle(wb, name, hs, rows=3, cols=seq_len(n), gridExpand=TRUE)
  for (i in seq_len(nrow(daten)))
    if (i%%2==0) addStyle(wb, name, grau, rows=i+3,
                          cols=seq_len(n), gridExpand=TRUE, stack=TRUE)
}

# Sheet 1: Kennzahlen
addWorksheet(wb, "Kennzahlen")
kenn <- data.frame(
  Kennzahl = c(
    "Block 4: Raeumliche Metriken  |  Igelbesenderung Niedersachsen",
    "HINWEIS: Homerange-Groessen (aKDE) kommen aus Block4b_Einzeltier_aKDE.R",
    "Kennzahl",
    "Tiere mit Nacht-Fixes (N)",
    "Nacht-Fixes gesamt",
    "Median Beobachtungsnaechte",
    paste0("Rueckkehr <= ", rueckkehr_m, " m: N Tiere (%)"),
    "Median Rueckkehrrate (%)",
    "Max. Nachtdistanz: Populationsmedian (m)",
    "Max. Nachtdistanz: Range (m)",
    "Radius of Gyration: Populationsmedian (m)",
    "Radius of Gyration: Range (m)",
    "Zentroid-Distanz Gehege: Populationsmedian (m)",
    "Min. Distanz Gehege: Populationsmin (m)"
  ),
  Wert = c("", "", "Wert",
    as.character(nrow(metriken_summary)),
    format(nrow(alle_fixes), big.mark="'"),
    as.character(round(median(metriken_summary$n_naechte))),
    paste0(n_rueck, " / ", round(100*n_rueck/nrow(metriken_summary)), "%"),
    as.character(round(median(metriken_summary$pct_rueckkehr, na.rm=TRUE), 1)),
    as.character(round(median(metriken_summary$max_dist_median_m, na.rm=TRUE))),
    paste0(round(min(metriken_summary$max_dist_median_m, na.rm=TRUE)), " - ",
           round(max(metriken_summary$max_dist_max_m, na.rm=TRUE))),
    as.character(round(median(metriken_summary$rog_m, na.rm=TRUE))),
    paste0(round(min(metriken_summary$rog_m, na.rm=TRUE)), " - ",
           round(max(metriken_summary$rog_m, na.rm=TRUE))),
    as.character(round(median(metriken_summary$dist_gehege_med_m, na.rm=TRUE))),
    as.character(round(min(metriken_summary$min_dist_gehege_m, na.rm=TRUE)))
  ), stringsAsFactors=FALSE
)
writeData(wb, "Kennzahlen", kenn, startRow=1, colNames=FALSE)
addStyle(wb, "Kennzahlen", hs, rows=3, cols=1:2, gridExpand=TRUE)
setColWidths(wb, "Kennzahlen", cols=1:2, widths=c(55, 25))

# Sheet 2: Distanzmetriken
schreibe_sheet(wb, "Distanzmetriken",
  paste0("Distanzmetriken pro Tier (Populationsmediane)\n",
         "Radius of Gyration = sqrt(mittl. quadr. Abstand vom Nacht-Schwerpunkt)"),
  metriken_summary[order(rog_m), .(
    Igel             = igel,
    `Naechte`        = n_naechte,
    `N Fixes`        = n_fixes_gesamt,
    `Max. Dist. (m)` = max_dist_max_m,
    `Mittl. Max. Dist. (m)` = max_dist_mean_m,
    `Median Max. Dist. (m)` = max_dist_median_m,
    `Radius of Gyration (m)` = rog_m,
    `Zentroid-Dist. Gehege (m)` = dist_gehege_med_m,
    `Min. Dist. Gehege (m)` = min_dist_gehege_m,
    Geschlecht = sex, Saison = saison_auswild
  )])
setColWidths(wb, "Distanzmetriken", cols=1:11,
  widths=c(10,9,10,14,18,18,20,22,20,11,14))

# Sheet 3: Pro-Nacht-Distanzen
schreibe_sheet(wb, "Dist_Naechte",
  "Naechtliche Raumnutzungsmetriken pro Igel x Nacht",
  nacht_metriken[order(igel, datum), .(
    Igel              = igel,
    Datum             = as.character(datum),
    `Tage seit Ausw.` = tage_seit,
    `N Fixes`         = n_fixes,
    `Max. Dist. (m)`  = round(max_dist_m, 0),
    `Zentroid-Dist. Gehege (m)` = round(dist_gehege_m, 0),
    `Min. Dist. Gehege (m)` = round(min_dist_gehege_m, 0),
    `RoG (m)`         = round(rog_m, 0),
    Rueckkehr         = ifelse(rueckkehr, "ja", "nein")
  )])
setColWidths(wb, "Dist_Naechte", cols=1:9, widths=c(10,13,16,9,14,22,20,10,11))

# Sheet 4: Rueckkehr
schreibe_sheet(wb, "Rueckkehr",
  paste0("Rueckkehr-Analyse: Naechte mit Annaeherung <= ", rueckkehr_m,
         " m ans Auswilderungsgehege\nSortiert nach Rueckkehrrate"),
  metriken_summary[order(-pct_rueckkehr), .(
    Igel                 = igel,
    `Naechte (N)`        = n_naechte,
    `Rueckkehren (N)`    = n_rueckkehr,
    `Rueckkehrrate (%)`  = pct_rueckkehr,
    `Med. Dist. (m)`     = dist_gehege_med_m,
    `Min. Dist. (m)`     = min_dist_gehege_m,
    Bewertung = fcase(
      pct_rueckkehr >= 75, "Haeufige Rueckkehr",
      pct_rueckkehr >= 40, "Gelegentliche Rueckkehr",
      default             = "Selten / nie"
    )
  )])
setColWidths(wb, "Rueckkehr", cols=1:7, widths=c(10,12,15,16,15,15,22))

# Sheet 5: Rueckkehr pro Nacht
schreibe_sheet(wb, "Rueckkehr_Naechte",
  "Naechte x Minimaldistanz Gehege (pro Igel x Nacht)",
  nacht_metriken[order(igel, datum), .(
    Igel              = igel,
    Datum             = as.character(datum),
    `Tage seit Ausw.` = tage_seit,
    `Min. Dist. (m)`  = round(min_dist_gehege_m, 0),
    Rueckkehr         = ifelse(rueckkehr, "ja", "nein")
  )])
setColWidths(wb, "Rueckkehr_Naechte", cols=1:5, widths=c(10,13,16,14,12))

saveWorkbook(wb, file.path(out_ordner, "Block4_Kernel_Uebersicht.xlsx"), overwrite=TRUE)
cat("  Block4_Kernel_Uebersicht.xlsx gespeichert\n\n")

# ── 8. Plots ───────────────────────────────────────────────────
cat("Erstelle Plots...\n")

# Plot 1: Rueckkehr und Zentroid-Distanz
metriken_ord <- metriken_summary[order(-pct_rueckkehr)]
metriken_ord[, igel_f := factor(igel, levels=igel)]

p_rueck <- ggplot(metriken_ord, aes(x=igel_f, y=dist_gehege_med_m)) +
  geom_hline(yintercept=rueckkehr_m, linetype="dashed",
             color="grey40", linewidth=0.8) +
  geom_col(aes(fill=pct_rueckkehr), alpha=0.85) +
  geom_point(aes(y=min_dist_gehege_m), shape=18, size=3.5, color="#C0392B") +
  scale_fill_gradient(low="#f7fbff", high="#2171b5", name="Return rate (%)") +
  annotate("text", x=Inf, y=rueckkehr_m+8,
           label=paste0(rueckkehr_m, " m threshold"),
           hjust=1.1, size=3, color="grey40") +
  coord_flip() +
  labs(title="Nightly centroid distance to release aviary",
       subtitle="Bar = median centroid | Diamond (red) = minimum | Dashed = return threshold",
       x=NULL, y="Distance to aviary (m)") +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")

ggsave(file.path(out_ordner, "metriken_rueckkehr.png"),
       p_rueck, width=9, height=8, dpi=150)
cat("  metriken_rueckkehr.png\n")

# Plot 2: Radius of Gyration
p_rog <- ggplot(metriken_ord, aes(x=igel_f, y=rog_m, fill=saison_auswild)) +
  geom_col(alpha=0.85) +
  scale_fill_viridis_d(option="plasma", end=0.85, na.value="grey70",
                       name="Release season") +
  coord_flip() +
  labs(title="Radius of Gyration per hedgehog",
       subtitle="RoG = sqrt(mean squared distance from nightly centroid) — fix-count independent",
       x=NULL, y="Radius of Gyration (m)") +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")

ggsave(file.path(out_ordner, "metriken_rog.png"),
       p_rog, width=9, height=8, dpi=150)
cat("  metriken_rog.png\n\n")

# ── 9. Zusammenfassung ─────────────────────────────────────────
cat("=========================================\n")
cat("ERGEBNISZUSAMMENFASSUNG Block 4\n")
cat("=========================================\n")
cat(sprintf("N = %d Tiere mit Nacht-Fixes\n", nrow(metriken_summary)))
cat(sprintf("Max. Nachtdistanz: Median=%d m (Range %d-%d m)\n",
    round(median(metriken_summary$max_dist_median_m, na.rm=TRUE)),
    round(min(metriken_summary$max_dist_median_m, na.rm=TRUE)),
    round(max(metriken_summary$max_dist_max_m, na.rm=TRUE))))
cat(sprintf("Radius of Gyration: Median=%d m (Range %d-%d m)\n",
    round(median(metriken_summary$rog_m, na.rm=TRUE)),
    round(min(metriken_summary$rog_m, na.rm=TRUE)),
    round(max(metriken_summary$rog_m, na.rm=TRUE))))
cat(sprintf("Zentroid-Distanz Gehege: Median=%d m\n",
    round(median(metriken_summary$dist_gehege_med_m, na.rm=TRUE))))
cat(sprintf("Rueckkehr <=%dm: %d/%d Tiere\n",
    rueckkehr_m, n_rueck, nrow(metriken_summary)))

cat(sprintf("\nOutputs in: %s\n", out_ordner))
cat("  kernel_ergebnisse.rds           --> Block4b_Einzeltier_aKDE.R\n")
cat("  Block4_Kernel_Uebersicht.xlsx   --> Distanzmetriken (OHNE aKDE-Homeranges)\n")
cat("  metriken_rueckkehr.png | metriken_rog.png\n\n")
cat("Reihenfolge:\n")
cat("  Schritt 1 (fertig):    Block4_Kernel_HomeRange.R  -> kernel_ergebnisse.rds\n")
cat("  Schritt 2 (naechster): Block4b_Einzeltier_aKDE.R  -> aKDE-Homeranges\n")
cat("  Schritt 3 (danach):    Block4_Kernel_Statistik.R  -> Gruppenvergleiche\n")
