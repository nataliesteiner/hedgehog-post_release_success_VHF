# ==============================================================
# Block 5 — environmental data: weather x home range & habitat
# ==============================================================
# Project:  Hedgehog VHF telemetry, Lower Saxony
# Author:   Natalie Steiner
#
# Main question:
#   Does the hedgehogs' home range change with season /
#   temperature — and if so, how?
#
# Data sources:
#   DWD   — German Weather Service (package: rdwd)
#   ATKIS — Basis-DLM, BKG / LGLN Lower Saxony (habitat, paths)
#           Preparation: run Habitat_Download_ATKIS.R once
#   Moon phase — computed directly in R, no download needed
#
# Requirement:
#   Block4b_Einzeltier_aKDE.R or Block4_Kernel_HomeRange.R
#
# Output:  output/Block5_Umweltdaten/

# ── 0. Pakete ──────────────────────────────────────────────────
pakete_kern <- c("data.table", "ggplot2", "patchwork", "scales",
                 "viridis", "lubridate", "sf", "officer",
                 "flextable", "openxlsx")
pakete_wetter  <- "rdwd"
pakete_habitat <- "osmdata"

alle_pakete <- c(pakete_kern, pakete_wetter, pakete_habitat)
fehlend <- alle_pakete[!sapply(alle_pakete, requireNamespace, quietly = TRUE)]
if (length(fehlend) > 0) {
  cat("Installiere fehlende Pakete:", paste(fehlend, collapse = ", "), "\n")
  install.packages(fehlend)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
  library(lubridate)
  library(sf)
})
if (requireNamespace("rdwd",     quietly = TRUE)) library(rdwd)
if (requireNamespace("officer",  quietly = TRUE)) library(officer)
if (requireNamespace("flextable",quietly = TRUE)) library(flextable)
if (requireNamespace("openxlsx", quietly = TRUE)) library(openxlsx)

rdwd_ok    <- requireNamespace("rdwd",    quietly = TRUE)
osmdata_ok <- requireNamespace("osmdata", quietly = TRUE)
cat(sprintf("Pakete geladen | rdwd: %s | osmdata: %s\n\n", rdwd_ok, osmdata_ok))

# ── 1. Pfade ───────────────────────────────────────────────────
projekt_root <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
out_ordner   <- file.path(projekt_root, "output", "Block5_Umweltdaten")
dir.create(out_ordner, showWarnings = FALSE, recursive = TRUE)

rds_4b      <- file.path(projekt_root, "output", "Block4b_Einzeltier", "Block4b_Ergebnisse.rds")
rds_4       <- file.path(projekt_root, "output", "Block4_Kernel",     "kernel_ergebnisse.rds")
gpkg_ordner <- file.path(projekt_root, "data", "kernel_files")
meta_datei  <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")

lat_studie <- 52.3970
lon_studie <-  9.2164
ziel_epsg  <- 32632   # EPSG:32632 = WGS84/UTM32N (native CRS der GPKG-Dateien)

# OSM-Habitatdaten:
#   OSM_AN <- TRUE   = Habitatkarte laden (nur wenn Overpass-Server erreichbar)
#   OSM_AN <- FALSE  = OSM überspringen (sicherer Standard — kein Hang-Risiko)
#
# ACHTUNG: osmdata hat ein internes Retry-Backoff (60 s) das OSM_TIMEOUT ignoriert.
# Das Skript fängt Backoff-Meldungen jetzt via withCallingHandlers ab und bricht ab.
# Falls der Server trotzdem haengt: OSM_AN <- FALSE setzen und Skript neu starten.
OSM_AN      <- FALSE  # Standard FALSE — auf TRUE setzen wenn OSM gewuenscht
OSM_TIMEOUT <- 20     # HTTP-Timeout pro Einzelabfrage (Sekunden)

cat("Output:", out_ordner, "\n\n")

# ── 2. Daten laden ─────────────────────────────────────────────
cat("Lade Igel-Daten...\n")

# Hilfsfunktion: Datum parsen (TT.MM.JJJJ oder JJJJ-MM-TT)
parse_datum_b5 <- function(x) {
  x <- as.character(x)
  r <- suppressWarnings(as.Date(x, format = "%d.%m.%Y"))
  r[is.na(r)] <- suppressWarnings(as.Date(x[is.na(r)]))
  r
}

# ── aKDE / KDE-Ergebnisse aus Block4b ──────────────────────────
if (file.exists(rds_4b)) {
  cat("  Nutze Block4b-Ergebnisse (aKDE)\n")
  b4b        <- readRDS(rds_4b)
  et_results <- b4b$einzeltier_results
  kde_dt     <- b4b$vergleich_dt   # enthaelt akde_95ha, kde_95ha, etc.
} else {
  cat("  [WARN] Block4b_Ergebnisse.rds nicht gefunden — bitte Block4b ausfuehren\n")
  et_results <- NULL
  kde_dt     <- NULL
}

# ── alle_fixes: immer direkt aus kernel_files/*.gpkg ──────────
# Block5 laedt die Rohdaten immer frisch aus den GPKG-Dateien.
# Das macht Block5 unabhaengig von kernel_ergebnisse.rds und
# vermeidet Probleme mit veralteten oder leeren RDS-Dateien.
{
  cat("  Lade Nacht-Fixes direkt aus kernel_files/*.gpkg\n")

  # Metadaten laden direkt aus Raw_data-Sheet
  # Bekannte Spalten: individual, sex, diagnosis_main, time_reha, date_release
  meta_b5 <- tryCatch({
    m  <- as.data.table(readxl::read_excel(path.expand(meta_datei), sheet = "Raw_data"))
    data.table(
      igel         = trimws(as.character(m$individual)),
      sex          = as.character(m$sex),
      diag_gruppe  = as.character(m$diagnosis_main),
      release_date = as.Date(m$date_release),
      time_reha    = suppressWarnings(as.numeric(m$time_reha))
    )
  }, error = function(e) {
    cat(sprintf("  [WARN] Metadaten nicht geladen: %s\n", conditionMessage(e)))
    data.table(igel=character(), sex=character(), diag_gruppe=character(),
               release_date=as.Date(NA), time_reha=numeric())
  })

  # Saison aus Release-Datum ableiten (Fruehling/Sommer/Herbst/Winter)
  monat_zu_saison <- function(m) {
    dplyr_ok <- FALSE  # kein dplyr noetig
    ifelse(m %in% 3:5,  "Fruehling",
    ifelse(m %in% 6:8,  "Sommer",
    ifelse(m %in% 9:11, "Herbst", "Winter")))
  }
  meta_b5[, saison_auswild := monat_zu_saison(as.integer(format(release_date, "%m")))]

  # alter: nicht im Excel — auf NA setzen (kann spaeter manuell ergaenzt werden)
  meta_b5[, alter := NA_character_]

  cat(sprintf("  %d Igel in Metadaten (sex, diagnose, reha, saison)\n", nrow(meta_b5)))

  # Alle GPKG-Dateien einlesen
  gpkg_files <- list.files(path.expand(gpkg_ordner), pattern = "\\.gpkg$",
                            full.names = TRUE)
  gpkg_files <- gpkg_files[!grepl("\\(1\\)", gpkg_files)]
  cat(sprintf("  Lese %d GPKG-Dateien...\n", length(gpkg_files)))

  fix_list <- lapply(gpkg_files, function(f) {
    igel_i <- sub("_Sachsenhagen_.*\\.gpkg$", "", basename(f))
    tryCatch({
      gdf <- sf::st_read(f, quiet = TRUE)
      gdf <- gdf[!sf::st_is_empty(gdf) & !is.na(sf::st_geometry(gdf)), ]
      if ("Night" %in% names(gdf)) gdf <- gdf[gdf$Night > 0, ]
      if (nrow(gdf) == 0) return(NULL)

      # Koordinaten direkt aus Geometrie-Objekten extrahieren
      # (st_coordinates hat in manchen sf-Versionen Probleme mit GPKG-Daten)
      geom_list <- sf::st_geometry(gdf)
      x_vals <- vapply(geom_list, function(g) {
        v <- tryCatch(as.numeric(g), error = function(e) c(NA_real_, NA_real_))
        if (length(v) >= 1) v[1] else NA_real_
      }, numeric(1))
      y_vals <- vapply(geom_list, function(g) {
        v <- tryCatch(as.numeric(g), error = function(e) c(NA_real_, NA_real_))
        if (length(v) >= 2) v[2] else NA_real_
      }, numeric(1))

      # Zeilen ohne gueltige Koordinaten entfernen
      valid <- !is.na(x_vals) & !is.na(y_vals) & x_vals != 0 & y_vals != 0
      if (sum(valid) == 0) return(NULL)
      x_vals <- x_vals[valid]
      y_vals <- y_vals[valid]
      gdf    <- gdf[valid, ]

      time_col <- intersect(c("X_time","_time","time","Time","datetime","Datetime"), names(gdf))[1]
      # Diagnose: erstes Timestamp-Beispiel ausgeben (nur beim ersten Igel)
      if (igel_i == sort(sub("_Sachsenhagen_.*\\.gpkg$", "", basename(gpkg_files)))[1]) {
        if (!is.na(time_col)) {
          ex <- as.character(gdf[[time_col]][1])
          cat(sprintf("    Timestamp-Format (%s): '%s'\n", time_col, ex))
        }
      }
      ts <- tryCatch({
        if (!is.na(time_col)) {
          ts_raw <- as.character(gdf[[time_col]])
          # Versuche zuerst direkte POSIXct-Konversion (funktioniert wenn GPKG POSIXct speichert)
          parsed <- tryCatch(
            as.POSIXct(gdf[[time_col]], tz = "UTC"),
            error = function(e) {
              # Fallback: Zeichen-basiertes Parsen mit lubridate
              lubridate::parse_date_time(
                ts_raw,
                orders = c("YmdHMSz","YmdHMS","Ymd HMS","YmdTHMS","YmdTHMSz",
                           "YmdHMOS","Ymd HM","dmYHMS","mdYHMS"),
                tz = "UTC", quiet = TRUE
              )
            }
          )
          lubridate::with_tz(parsed, "Europe/Berlin")
        } else {
          rep(as.POSIXct(NA), nrow(gdf))
        }
      }, error = function(e) rep(as.POSIXct(NA), nrow(gdf)))

      rel <- meta_b5[igel == igel_i, release_date]
      rel <- if (length(rel) && !is.na(rel[1])) rel[1] else as.Date(NA)

      data.table(
        igel         = igel_i,
        datetime     = ts,
        datum        = as.Date(ts, tz = "Europe/Berlin"),
        x            = x_vals,
        y            = y_vals,
        tageszeit    = "Nacht",
        release_date = rel,
        tage_seit    = as.integer(as.Date(ts, tz="Europe/Berlin") - rel)
      )
    }, error = function(e) {
      cat(sprintf("    [FEHLER] %s: %s\n", igel_i, conditionMessage(e)))
      NULL
    })
  })

  alle_fixes <- rbindlist(Filter(Negate(is.null), fix_list), fill = TRUE)
  # Nur Metadaten-Spalten mergen die tatsaechlich existieren
  meta_merge_cols <- intersect(c("igel","sex","alter","diag_gruppe","saison_auswild","time_reha"),
                                names(meta_b5))
  alle_fixes <- merge(alle_fixes, meta_b5[, ..meta_merge_cols],
                      by = "igel", all.x = TRUE)
  cat(sprintf("  %d Nacht-Fixes, %d Igel aus GPKG geladen\n",
              nrow(alle_fixes), alle_fixes[, uniqueN(igel)]))
}

if (!is.data.table(alle_fixes)) setDT(alle_fixes)
setorder(alle_fixes, igel, datetime)

# ── Diagnose: Koordinaten pruefen ──────────────────────────────
n_ok  <- alle_fixes[!is.na(x) & !is.na(y), .N]
n_ges <- nrow(alle_fixes)
cat(sprintf("  DIAGNOSE: %d/%d Fixes mit gueltigen Koordinaten (x,y)\n\n", n_ok, n_ges))
if (n_ok == 0) stop("Keine gueltigen Koordinaten in alle_fixes — GPKG-Laden pruefen!")

if (!is.null(et_results)) {
  nacht_metriken_alle <- rbindlist(lapply(names(et_results), function(ig) {
    nm <- et_results[[ig]]$nacht_metriken
    if (!is.null(nm) && nrow(nm) > 0) {
      nm[, igel := ig]; nm
    }
  }), fill = TRUE)
} else {
  release_punkt <- data.frame(x = 514743.7, y = 5805363.8)
  nacht_metriken_alle <- alle_fixes[, {
    cx_n <- mean(x); cy_n <- mean(y)
    dist_z <- sqrt((x - cx_n)^2 + (y - cy_n)^2)
    .(
      cx             = cx_n,
      cy             = cy_n,
      n_fixes_nacht  = .N,
      radius_95pct_m = quantile(dist_z, 0.95),
      radius_mean_m  = mean(dist_z),
      dist_release_m = sqrt((cx_n - release_punkt$x)^2 +
                            (cy_n - release_punkt$y)^2),
      tage_seit      = first(tage_seit)
    )
  }, by = .(igel, datum)]
}

cat(sprintf("  %d Fixes, %d Igel, %d Naechte\n\n",
            nrow(alle_fixes),
            alle_fixes[, uniqueN(igel)],
            nrow(nacht_metriken_alle)))

datum_range <- range(alle_fixes$datum, na.rm = TRUE)
cat(sprintf("  Beobachtungszeitraum: %s bis %s\n\n",
            datum_range[1], datum_range[2]))


# ── 3. MONDPHASE — entfernt (nicht relevant fuer Rehabilitationsfrage) ────────


# ── 4. DWD WETTERDATEN ────────────────────────────────────────
cat("====================================================\n")
cat("4. DWD Wetterdaten\n")
cat("====================================================\n\n")

wetter_dt    <- NULL
# BUG-FIX: station_name / station_id / station_dist werden nur gesetzt wenn
# rdwd verfuegbar ist UND der DWD-Download erfolgreich war. Ohne Initialisierung
# wuerde der Word-Export mit "object not found" abbrechen wenn rdwd fehlt.
station_name <- "n. v."
station_id   <- NA_character_
station_dist <- NA_real_

if (!rdwd_ok) {
  cat("  rdwd nicht verfuegbar\n")
  cat("  Installiere mit: install.packages('rdwd')\n\n")
} else {

  cat("  Suche naechste DWD-Station zu Sachsenhagen...\n")
  data("metaIndex", package = "rdwd", envir = environment())
  meta <- as.data.table(metaIndex)

  # Spalten automatisch erkennen (rdwd aendert Namen je nach Version)
  # Hilfsfunktion: ersten Treffer aus einer Liste von Kandidaten
  find_col <- function(candidates, nm) {
    treffer <- intersect(candidates, nm)
    if (length(treffer) == 0) NULL else treffer[1]
  }

  col_lat  <- find_col(c("geoBreite","lat","Breite"),          names(meta))
  col_lon  <- find_col(c("geoLaenge","lon","geoLon","Laenge"), names(meta))
  col_id   <- find_col(c("Stations_id","id","ID"),             names(meta))
  col_name <- find_col(c("Stationsname","name","Name"),        names(meta))
  col_end  <- find_col(c("bis","Ende","end","date_end","Enddatum"), names(meta))

  cat(sprintf("  rdwd-Spalten: lat=%s | lon=%s | id=%s | end=%s\n",
              col_lat %||% "?", col_lon %||% "?",
              col_id  %||% "?", col_end %||% "nicht gefunden"))

  if (is.null(col_lat) || is.null(col_lon) || is.null(col_id)) {
    cat("  Verfuegbare Spaltennamen im metaIndex:\n")
    cat(" ", paste(names(meta), collapse=", "), "\n")
    stop("Pflicht-Spalten nicht gefunden — rdwd-Version pruefen.")
  }

  meta_kl <- meta[res == "daily" & var == "kl" &
                   !is.na(get(col_lat)) & !is.na(get(col_lon))]

  meta_kl[, dist := sqrt((get(col_lat) - lat_studie)^2 +
                          (get(col_lon) - lon_studie)^2)]

  # Pro Station nur eine Zeile behalten (neueste per-Kategorie bevorzugen)
  meta_kl_unique <- meta_kl[order(dist, -as.integer(per == "recent"))
                             ][!duplicated(get(col_id))]

  top5 <- meta_kl_unique[order(dist)][seq_len(min(5,.N)), .(
    ID           = get(col_id),
    Name         = get(col_name),
    Lat          = round(get(col_lat), 3),
    Lon          = round(get(col_lon), 3),
    `Dist. (km)` = round(dist * 111, 1)
  )]
  cat("  Naechste DWD-Stationen:\n")
  print(top5)
  cat("\n")

  obs_bis <- datum_range[2]

  # Stationsauswahl: Enddatum pruefen wenn Spalte vorhanden
  if (!is.null(col_end)) {
    # Robuste Datumskonvertierung (kann numerisch oder character sein)
    end_raw <- meta_kl[[col_end]]
    end_dates <- tryCatch(
      as.Date(as.character(end_raw),
              format = ifelse(nchar(as.character(end_raw[1])) == 8,
                              "%Y%m%d", "%Y-%m-%d")),
      error = function(e) rep(as.Date(NA), nrow(meta_kl))
    )
    meta_kl[, bis_date := end_dates]
    station_kandidaten <- meta_kl[order(dist)][
      (!is.na(bis_date) & bis_date >= obs_bis) | per == "recent"]
  } else {
    # Kein Enddatum verfuegbar — nur recent nehmen
    station_kandidaten <- meta_kl[order(dist)][per == "recent"]
  }

  if (nrow(station_kandidaten) == 0) {
    station_kandidaten <- meta_kl[order(dist)]
  }

  # Pro Station eine Zeile
  station_kandidaten <- station_kandidaten[
    !duplicated(get(col_id))][order(dist)]

  station_id   <- station_kandidaten[1, get(col_id)]
  station_name <- station_kandidaten[1, get(col_name)]
  station_dist <- round(station_kandidaten[1, dist] * 111, 1)

  cat(sprintf("  Ausgewaehlt: %s (ID: %s, %.1f km)\n\n",
              station_name, station_id, station_dist))

  cat("  Lade DWD-Daten herunter...\n")
  tryCatch({
    link_r <- tryCatch(
      selectDWD(id = station_id, res = "daily", var = "kl", per = "recent"),
      error = function(e) NULL
    )
    link_h <- tryCatch(
      selectDWD(id = station_id, res = "daily", var = "kl", per = "historical"),
      error = function(e) NULL
    )

    dwd_rohdaten <- list()
    if (!is.null(link_h) && length(link_h) > 0) {
      dwd_rohdaten[["hist"]] <- tryCatch(
        dataDWD(link_h[[1]], read = TRUE, dir = tempdir(), quiet = TRUE),
        error = function(e) NULL
      )
    }
    if (!is.null(link_r) && length(link_r) > 0) {
      dwd_rohdaten[["recent"]] <- tryCatch(
        dataDWD(link_r[[1]], read = TRUE, dir = tempdir(), quiet = TRUE),
        error = function(e) NULL
      )
    }

    dwd_raw <- rbindlist(Filter(Negate(is.null), dwd_rohdaten),
                         fill = TRUE, use.names = TRUE)

    if (nrow(dwd_raw) == 0) stop("Keine DWD-Daten geladen")

    datumscol <- grep("MESS_DATUM|datum|Datum", names(dwd_raw),
                      ignore.case = TRUE, value = TRUE)[1]
    dwd_raw[, datum := as.Date(as.character(get(datumscol)),
                                format = ifelse(
                                  nchar(as.character(get(datumscol)[1])) == 8,
                                  "%Y%m%d", "%Y-%m-%d"))]

    col_tmk <- grep("^TMK|temp.*mittel|mean.*temp", names(dwd_raw),
                    ignore.case = TRUE, value = TRUE)[1]
    col_tnk <- grep("^TNK|temp.*min|min.*temp",    names(dwd_raw),
                    ignore.case = TRUE, value = TRUE)[1]
    col_txk <- grep("^TXK|temp.*max|max.*temp",    names(dwd_raw),
                    ignore.case = TRUE, value = TRUE)[1]
    col_rsk <- grep("^RSK|nieder|precip",           names(dwd_raw),
                    ignore.case = TRUE, value = TRUE)[1]
    col_sdk <- grep("^SDK|sonne|sunshine",          names(dwd_raw),
                    ignore.case = TRUE, value = TRUE)[1]

    wetter_dt <- dwd_raw[datum >= datum_range[1] - 7 &
                          datum <= datum_range[2] + 7,
                          .(datum,
                            temp_mittel = suppressWarnings(
                              as.numeric(get(col_tmk))),
                            temp_min    = if (!is.na(col_tnk))
                              suppressWarnings(as.numeric(get(col_tnk)))
                              else NA_real_,
                            temp_max    = if (!is.na(col_txk))
                              suppressWarnings(as.numeric(get(col_txk)))
                              else NA_real_,
                            niederschlag_mm = if (!is.na(col_rsk))
                              suppressWarnings(as.numeric(get(col_rsk)))
                              else NA_real_,
                            sonnenschein_h  = if (!is.na(col_sdk))
                              suppressWarnings(as.numeric(get(col_sdk)))
                              else NA_real_
                          )]

    for (col in names(wetter_dt)[-1]) {
      set(wetter_dt, which(wetter_dt[[col]] <= -990), col, NA_real_)
    }

    wetter_dt <- wetter_dt[!duplicated(datum)][order(datum)]

    wetter_dt[, temp_gleit5 := frollmean(temp_mittel, n = 5, align = "center",
                                          na.rm = TRUE)]
    wetter_dt[, saison_temp := fcase(
      temp_gleit5 < 5,                      "Winter (<5 Grad)",
      temp_gleit5 >= 5  & temp_gleit5 < 10, "Fruehjahr/Herbst (5-10 Grad)",
      temp_gleit5 >= 10 & temp_gleit5 < 18, "Sommer (10-18 Grad)",
      temp_gleit5 >= 18,                    "Hochsommer (>18 Grad)",
      default = "Unbekannt"
    )]

    wetter_dt[, regen := fifelse(!is.na(niederschlag_mm) &
                                  niederschlag_mm > 0.5, "Regen", "Trocken")]

    cat(sprintf("  DWD-Daten: %d Tage (%s bis %s)\n",
                nrow(wetter_dt), min(wetter_dt$datum), max(wetter_dt$datum)))
    cat(sprintf("  Station: %s | Abstand: %.1f km\n",
                station_name, station_dist))
    cat(sprintf("  Temperatur: %.1f bis %.1f Grad C (Mittel)\n\n",
                min(wetter_dt$temp_mittel, na.rm=TRUE),
                max(wetter_dt$temp_mittel, na.rm=TRUE)))

  }, error = function(e) {
    cat(sprintf("  DWD-Download fehlgeschlagen: %s\n", conditionMessage(e)))
    cat("  Tipp: Internetverbindung pruefen oder Station-ID manuell angeben.\n\n")
    wetter_dt <<- NULL
  })
}

# ── 5. WETTERDATEN MIT IGEL-DATEN VERKNUEPFEN (MONATLICH) ────────
# Umstrukturierung: statt tagesweiser Verknuepfung → monatliche Aggregation.
# Begruendung: 26 Tiere mit unterschiedlichen Beobachtungsdauern. Monatliche
# Skala ist statistisch sinnvoller und visuell klarer (Publikationsreif).
# Zusaetzlich: Regen-Effekt via Schwellenwert (> 2 mm = Regennacht).
cat("====================================================\n")
cat("5. Wetter x Igel-Aktivitaet (monatlich)\n")
cat("====================================================\n\n")

nacht_wetter  <- copy(nacht_metriken_alle)
fixes_wetter  <- copy(alle_fixes)
monat_wetter  <- NULL  # monatliche Wetter-Zusammenfassung
monat_igel    <- NULL  # monatliche Igel-Metriken

if (!is.null(wetter_dt) && nrow(wetter_dt) > 0) {

  # ── Diagnostik: Datumsabdeckung pruefen ────────────────────────
  cat(sprintf("  DWD-Daten: %s bis %s (%d Tage)\n",
              min(wetter_dt$datum), max(wetter_dt$datum), nrow(wetter_dt)))
  cat(sprintf("  Igel-Daten: %s bis %s\n", datum_range[1], datum_range[2]))

  # Schnittmenge pruefen
  n_match <- wetter_dt[datum >= datum_range[1] & datum <= datum_range[2], .N]
  cat(sprintf("  Schnittmenge: %d Tage mit Wetterdaten im Beobachtungszeitraum\n\n", n_match))

  if (n_match == 0) {
    cat("  [WARN] Keine Ueberlappung! Datumsformat pruefen.\n")
    cat(sprintf("  wetter_dt$datum Klasse: %s | datum_range Klasse: %s\n",
                class(wetter_dt$datum)[1], class(datum_range)[1]))
    # Versuch: beide zu Date konvertieren
    wetter_dt[, datum := as.Date(datum)]
    n_match <- wetter_dt[datum >= as.Date(datum_range[1]) &
                          datum <= as.Date(datum_range[2]), .N]
    cat(sprintf("  Nach Konvertierung: %d Tage\n\n", n_match))
  }

  # ── Tagesweise Verknuepfung (fuer Regen-Effekt) ────────────────
  nacht_wetter <- merge(nacht_metriken_alle, wetter_dt,
                         by = "datum", all.x = TRUE)
  fixes_wetter <- merge(alle_fixes, wetter_dt,
                         by = "datum", all.x = TRUE)
  cat(sprintf("  %d Naechte | %d davon mit Temp-Daten\n\n",
              nrow(nacht_metriken_alle),
              nacht_wetter[!is.na(temp_mittel), .N]))

  # ── Monatliche Wetter-Aggregation ──────────────────────────────
  wetter_dt[, monat_str := format(datum, "%Y-%m")]
  monat_wetter <- wetter_dt[, .(
    temp_mittel_monat  = round(mean(temp_mittel,      na.rm = TRUE), 1),
    temp_min_monat     = round(min(temp_mittel,        na.rm = TRUE), 1),
    temp_max_monat     = round(max(temp_mittel,        na.rm = TRUE), 1),
    niederschlag_sum   = round(sum(niederschlag_mm,    na.rm = TRUE), 1),
    sonnenschein_sum   = round(sum(sonnenschein_h,     na.rm = TRUE), 1),
    n_regentage        = sum(!is.na(niederschlag_mm) & niederschlag_mm > 2, na.rm = TRUE),
    n_tage_wetter      = .N
  ), by = monat_str][order(monat_str)]
  monat_wetter[, monat_date := as.Date(paste0(monat_str, "-01"))]

  cat(sprintf("  Monatliche Wetteraggregation: %d Monate\n", nrow(monat_wetter)))
  cat(sprintf("  Temp-Bereich: %.1f bis %.1f Grad C (Monatsmittel)\n",
              min(monat_wetter$temp_mittel_monat, na.rm = TRUE),
              max(monat_wetter$temp_mittel_monat, na.rm = TRUE)))
  cat(sprintf("  Niederschlag: %.0f bis %.0f mm/Monat\n\n",
              min(monat_wetter$niederschlag_sum, na.rm = TRUE),
              max(monat_wetter$niederschlag_sum, na.rm = TRUE)))

  # ── Monatliche Igel-Metriken ───────────────────────────────────
  if ("radius_95pct_m" %in% names(nacht_metriken_alle)) {
    nacht_metriken_alle[, monat_str := format(datum, "%Y-%m")]
    monat_igel <- nacht_metriken_alle[, .(
      median_radius_m   = round(median(radius_95pct_m,  na.rm = TRUE)),
      median_dist_m     = round(median(dist_release_m,  na.rm = TRUE)),
      n_naechte_igel    = .N,
      n_igel            = uniqueN(igel)
    ), by = monat_str][order(monat_str)]
    monat_igel[, monat_date := as.Date(paste0(monat_str, "-01"))]
    cat(sprintf("  Monatliche Igel-Metriken: %d Monate\n\n", nrow(monat_igel)))
  }

  # ── Statistik: Regen-Effekt (Schwellenwert > 2 mm) ────────────
  nacht_stat <- nacht_wetter[!is.na(temp_mittel) & !is.na(radius_95pct_m)]
  if (nrow(nacht_stat) >= 10) {
    nacht_stat[, regen_nacht := fifelse(!is.na(niederschlag_mm) &
                                         niederschlag_mm > 2,
                                         "Regennacht (>2 mm)", "Trockene Nacht")]
    regen_test <- tryCatch(
      wilcox.test(radius_95pct_m ~ regen_nacht, data = nacht_stat),
      error = function(e) NULL)
    if (!is.null(regen_test))
      cat(sprintf("  Wilcoxon Regen (>2mm) x Nachtradius: W=%.0f, p=%.4f\n",
                  regen_test$statistic, regen_test$p.value))

    sp_temp_rad <- tryCatch(
      cor.test(nacht_stat$temp_mittel, nacht_stat$radius_95pct_m,
               method = "spearman"),
      error = function(e) NULL)
    if (!is.null(sp_temp_rad))
      cat(sprintf("  Spearman Temp x Nachtradius:  rho=%.3f, p=%.4f\n",
                  sp_temp_rad$estimate, sp_temp_rad$p.value))
    cat("\n")
  } else {
    cat("  Zu wenige tagesweise Uebereinstimmungen fuer Korrelationstest.\n\n")
    nacht_stat <- data.table()
  }

} else {
  cat("  Keine Wetterdaten verfuegbar — Abschnitt 5 uebersprungen.\n\n")
  nacht_stat <- data.table()
}


# ── 6. HABITATDATEN (ATKIS Basis-DLM) ────────────────────────
cat("====================================================\n")
cat("6. Habitatdaten (ATKIS Basis-DLM)\n")
cat("====================================================\n\n")
# ══════════════════════════════════════════════════════════════
# HABITATDATEN — ATKIS BASIS-DLM (amtliche deutsche Geobasisdaten)
# ══════════════════════════════════════════════════════════════
#
# Datenquelle: ATKIS Basis-DLM — Bundesamt fuer Kartographie und
#   Geodaesie (BKG) / LGLN Niedersachsen.
#   Lizenz: DL-DE/BY-2.0 — © GeoBasis-DE / BKG 2024
#
# Vorbereitung:
#   Einmalig Habitat_Download_ATKIS.R ausfuehren. Das Skript laedt
#   die ATKIS-Objekte via WFS und erzeugt:
#     data/excel_files/habitat_sachsenhagen_atkis.gpkg  (Polygone)
#     data/excel_files/atkis_wege_sachsenhagen.gpkg     (Wege)
#
# Habitatklassen (konsistent mit Basis-DLM Objektartkatalog):
#   Wald            : AX_Wald (41001) → Waldinneres nach ±20m Randabzug
#   Gehoelz_Strauch : AX_Gehoelz (41003), AX_Hecke → Teil von "Edge habitat"
#   Acker           : AX_Ackerland (43001)
#   Wiese           : AX_Gruenland (43002), AX_Heide, AX_Moor
#   Gewaesser       : AX_FliessgewaesserAbschnitt (44001), AX_StehendesGewaesser
#   Siedlung        : AX_Wohnbauflaeche, AX_IndustrieUndGewerbeflaeche
#
# ══════════════════════════════════════════════════════════════

# ── SCHRITT 1: Habitat_Download_ATKIS.R einmalig ausfuehren ───────────────────
# Danach hier den Pfad zur erzeugten GPKG-Datei eintragen:
HABITAT_METHODE   <- "gpkg"
HABITAT_GPKG_PFAD <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF/data/excel_files/habitat_sachsenhagen_atkis.gpkg"

# ATKIS-Wege fuer Wegnaehe-Analyse (Tracks, Pfade aus Basis-DLM)
# Erzeugt von Habitat_Download_ATKIS.R (AX_Weg, AX_WegPfadSteig etc.)
# Fallback auf OSM-Wege wenn ATKIS-Wege noch nicht vorhanden
ATKIS_WEGE_PFAD <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF/data/excel_files/atkis_wege_sachsenhagen.gpkg"
OSM_ROADS_PFAD  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF/data/niedersachsen-260601-free.shp/gis_osm_roads_free_1.shp"
# Bevorzuge ATKIS-Wege; falls nicht vorhanden, automatisch auf OSM-Fallback
ROADS_PFAD <- if (file.exists(path.expand(ATKIS_WEGE_PFAD))) ATKIS_WEGE_PFAD else OSM_ROADS_PFAD
cat(sprintf("  Wege-Quelle: %s\n\n",
            if (ROADS_PFAD == ATKIS_WEGE_PFAD) "ATKIS Basis-DLM" else "OSM (Fallback)"))

# Puffer um Wege: Fixes innerhalb dieser Distanz → "Waldweg" / Edge habitat
WALDWEG_BUFFER_M <- 15  # Meter

habitat_fixes    <- NULL
habitat_polygone <- NULL

# ── Bounding Box aus Fix-Koordinaten ──────────────────────────
# Koordinaten sind in EPSG:32632 (native GPKG-CRS)
fixes_sf_ll <- st_transform(
  st_as_sf(alle_fixes[!is.na(x) & !is.na(y)],
           coords = c("x","y"), crs = 32632),
  crs = 4326)
bb <- st_bbox(fixes_sf_ll)
pad <- 0.015   # ~1.5 km Puffer um alle Fixes
bbox_wgs84 <- c(
  xmin = as.numeric(bb["xmin"]) - pad,
  ymin = as.numeric(bb["ymin"]) - pad,
  xmax = as.numeric(bb["xmax"]) + pad,
  ymax = as.numeric(bb["ymax"]) + pad
)
cat(sprintf("  Untersuchungsgebiet (WGS84): %.4f, %.4f bis %.4f, %.4f\n\n",
            bbox_wgs84["xmin"], bbox_wgs84["ymin"],
            bbox_wgs84["xmax"], bbox_wgs84["ymax"]))

# ── Hilfsfunktion: Habitatklassifikation + Extraktion ─────────
klassifiziere_und_extrahiere <- function(hab_sf,
                                         crs_out        = 32632,  # native GPKG-CRS
                                         osm_roads_pfad = ROADS_PFAD,
                                         waldweg_buffer = WALDWEG_BUFFER_M) {

  # Habitat-Shapefile in GPKG-CRS transformieren (32632 = WGS84/UTM32N)
  hab_sf <- st_transform(hab_sf, crs = crs_out)
  hab_sf <- hab_sf[st_is_valid(hab_sf), ]

  # Fixes als sf-Objekt — CRS = 32632 (native der GPKG-Dateien)
  af_ok        <- alle_fixes[!is.na(x) & !is.na(y)]
  fixes_sf_utm <- st_as_sf(af_ok, coords = c("x","y"), crs = crs_out)

  cat(sprintf("  %d Fixes zu klassifizieren\n", nrow(af_ok)))
  cat(sprintf("  Habitattypen im Shapefile: %s\n",
              paste(sort(unique(hab_sf$habitat)), collapse=", ")))

  # hpf: sequentielle Zeilennummern (seq_len(.N) statt .I!)
  hpf <- data.table(
    igel        = af_ok$igel,
    datum       = af_ok$datum,
    x           = af_ok$x,
    y           = af_ok$y,
    tage_seit   = if ("tage_seit" %in% names(af_ok)) af_ok$tage_seit else NA_integer_,
    habitat_typ = "Offenland"
  )

  # ── Basis-Klassifikation nach Shapefile-Polygonen ─────────────
  prioritaet <- c("Gehoelz_Strauch", "Gewaesser", "Garten",
                   "Siedlung", "Wald", "Wiese", "Acker")

  for (hab in prioritaet) {
    hab_poly <- hab_sf[hab_sf$habitat == hab, ]
    if (nrow(hab_poly) == 0) {
      cat(sprintf("    [SKIP] %s: 0 Polygone\n", hab))
      next
    }
    in_hab <- tryCatch(
      suppressMessages(lengths(st_intersects(fixes_sf_utm, st_union(hab_poly))) > 0),
      error = function(e) { cat(sprintf("    [ERR] %s: %s\n", hab, conditionMessage(e))); rep(FALSE, nrow(fixes_sf_utm)) }
    )
    noch_off <- hpf$habitat_typ == "Offenland"
    hpf[noch_off & in_hab, habitat_typ := hab]
    cat(sprintf("    %s: %d Fixes\n", hab, sum(noch_off & in_hab)))
  }

  # ── Waldrand: ±20 m beidseitig der Waldgrenze ─────────────────
  wald_poly <- hab_sf[hab_sf$habitat == "Wald", ]
  if (nrow(wald_poly) > 0) {
    tryCatch({
      wald_union    <- st_union(wald_poly)
      wald_aussen   <- st_buffer(wald_union,  20)
      # Negativer Innen-Buffer: nur wenn Geometrie gross genug ist
      wald_innen    <- tryCatch(st_buffer(wald_union, -20), error = function(e) wald_union)
      waldrand_zone <- st_make_valid(st_difference(wald_aussen, wald_innen))
      in_rand       <- suppressMessages(lengths(st_intersects(fixes_sf_utm, waldrand_zone)) > 0)
      hpf[in_rand, habitat_typ := "Waldrand"]
      cat(sprintf("  Waldrand (±20 m): %d Fixes\n", sum(in_rand)))
    }, error = function(e) cat(sprintf("  [WARN] Waldrand-Berechnung: %s\n", conditionMessage(e))))
  }

  # ── Waldweg: Naehe zu Wegen (ATKIS bevorzugt, OSM als Fallback) ────────────
  roads_pfad_exp <- path.expand(osm_roads_pfad)
  is_atkis_wege  <- grepl("atkis_wege", roads_pfad_exp)

  if (file.exists(roads_pfad_exp)) {
    cat(sprintf("  Lade Wege (%s, Puffer %d m)...\n",
                ifelse(is_atkis_wege, "ATKIS Basis-DLM", "OSM"),
                waldweg_buffer))

    # ATKIS GPKG: direkt laden (keine Klassen-Filterung noetig, schon vorverarbeitet)
    # OSM SHP: nur Track/Pfad-Klassen behalten
    roads_raw <- tryCatch({
      if (is_atkis_wege) {
        # ATKIS GPKG (alle Eintraege sind relevante Wege)
        st_read(roads_pfad_exp, quiet = TRUE)
      } else {
        # OSM SHP: raeumlich vorfiltern
        bbox_ll_fix <- st_bbox(st_transform(fixes_sf_utm, 4326))
        st_read(roads_pfad_exp, quiet = TRUE,
                wkt_filter = st_as_text(st_as_sfc(st_bbox(c(
                  xmin = as.numeric(bbox_ll_fix["xmin"]) - 0.02,
                  ymin = as.numeric(bbox_ll_fix["ymin"]) - 0.02,
                  xmax = as.numeric(bbox_ll_fix["xmax"]) + 0.02,
                  ymax = as.numeric(bbox_ll_fix["ymax"]) + 0.02), crs = 4326))))
      }
    }, error = function(e) {
      cat(sprintf("  [WARN] Wege laden: %s\n", conditionMessage(e)))
      NULL
    })

    if (!is.null(roads_raw) && nrow(roads_raw) > 0) {
      roads_utm <- st_transform(roads_raw, crs = crs_out)

      # OSM: nur Track/Pfad-Typen behalten
      if (!is_atkis_wege && "fclass" %in% names(roads_utm)) {
        weg_typen <- c("track","track_grade1","track_grade2","track_grade3",
                       "track_grade4","track_grade5","path","footway","cycleway")
        roads_utm <- roads_utm[roads_utm$fclass %in% weg_typen, ]
      }

      cat(sprintf("  %d Weg-Segmente\n", nrow(roads_utm)))
      if (nrow(roads_utm) > 0) {
        tryCatch({
          wege_buffer <- st_union(st_buffer(roads_utm, waldweg_buffer))
          in_weg <- suppressMessages(lengths(st_intersects(fixes_sf_utm, wege_buffer)) > 0)
          hpf[in_weg, habitat_typ := "Waldweg"]
          cat(sprintf("  Waldweg (<= %d m): %d Fixes\n", waldweg_buffer, sum(in_weg)))
        }, error = function(e) cat(sprintf("  [WARN] Waldweg-Buffer: %s\n", conditionMessage(e))))
      }
    }
  } else {
    cat("  Weg-Datei nicht gefunden — Waldweg-Analyse uebersprungen\n")
    cat(sprintf("  (Erwartet: %s)\n", basename(roads_pfad_exp)))
    cat("  Tipp: Habitat_Download_ATKIS.R ausfuehren\n")
  }

  # ── Zusammenfassung Feinklassifikation ──────────────────────────
  cat("\n  Feinklassifikation:\n")
  print(hpf[, .N, by = habitat_typ][order(-N)])

  # ── WICHTIG: Waldweg + Waldrand → "Edge habitat" ─────────────────
  # Begruendung: tRackIT-Lokalisationsgenauigkeit ca. 20–50 m.
  #   - Waldweg-Buffer (15 m) < Mindestgenauigkeit → Waldweg und
  #     angrenzendes Wald/Waldrand nicht zuverlaessig trennbar.
  #   - Waldrand-Buffer (±20 m) = untere Genauigkeitsgrenze.
  #   → Beide Klassen werden zu "Edge habitat" zusammengefasst
  #     (Randstrukturen und Wegekorridore im Wald).
  #   Die Feinunterscheidung bleibt in habitat_typ erhalten.
  hpf[, habitat_grob := fcase(
    # Edge habitat = Randstrukturen + Wege + Gehölze
    # Begruendung: Gehoelz_Strauch-Polygone sind zu klein (<50m) fuer
    # zuverlaessige Klassifikation bei 20-50m tRackIT-Genauigkeit;
    # oekologisch dieselbe Funktion wie Waldrand (Deckung, Nahrung).
    habitat_typ %in% c("Waldweg", "Waldrand", "Gehoelz_Strauch"), "Edge habitat",
    habitat_typ == "Wald",    "Forest interior",
    habitat_typ == "Acker",   "Agricultural land",
    habitat_typ == "Wiese",   "Grassland",
    habitat_typ == "Gewaesser", "Water body",
    habitat_typ == "Siedlung",  "Settlement",
    habitat_typ == "Garten",    "Garden",
    default = "Open land"
  )]
  cat("\n  Zusammengefasste Klassifikation (habitat_grob):\n")
  print(hpf[, .N, by = habitat_grob][order(-N)])
  cat("\n")

  list(habitat_fixes = hpf, habitat_polygone = hab_sf)
}

# ══════════════════════════════════════════════════════════════
# OPTION A (Standard): overpass_klein
#   Einmalige kombinierte Overpass-Abfrage fuer das kleine Studiengebiet.
#   Keine GDAL-Installation noetig, keine langen Downloads.
#   Das ~5 km² Gebiet um Sachsenhagen liefert in < 10 Sekunden Ergebnisse.
# ══════════════════════════════════════════════════════════════
if (HABITAT_METHODE == "overpass_klein") {

  if (!requireNamespace("osmdata", quietly=TRUE)) install.packages("osmdata", quiet=TRUE)

  if (requireNamespace("osmdata", quietly=TRUE)) {
    library(osmdata)
    cat("  Lade Habitatdaten via Overpass API (kleine Abfrage, Studiengebiet)...\n")

    # Bounding Box fuer Overpass: bottom, left, top, right
    bb_op <- c(bbox_wgs84["ymin"], bbox_wgs84["xmin"],
               bbox_wgs84["ymax"], bbox_wgs84["xmax"])

    # Hilfsfunktion: einzelne Abfrage mit backoff-Schutz
    abfrage_osm <- function(key, value, bbox, label) {
      cat(sprintf("    %s=%s... ", key, value))
      tryCatch(
        withCallingHandlers({
          q   <- opq(bbox=bbox, timeout=25) |>
                   add_osm_feature(key=key, value=value)
          res <- osmdata_sf(q)
          polys <- res$osm_polygons
          cat(sprintf("%d Polygone\n", if (!is.null(polys)) nrow(polys) else 0))
          polys
        },
        message = function(m) {
          if (grepl("backoff|Waiting|retry|Rate", conditionMessage(m), ignore.case=TRUE)) {
            cat("Server Rate-Limit -- ueberspringe\n")
            invokeRestart("muffleMessage")
            stop("BACKOFF")
          }
          invokeRestart("muffleMessage")
        }),
        error = function(e) {
          if (!grepl("BACKOFF", conditionMessage(e)))
            cat(sprintf("Fehler: %s\n", substr(conditionMessage(e), 1, 60)))
          NULL
        }
      )
    }

    # Landnutzungs-Abfragen (jeweils separate Abfrage)
    abfragen <- list(
      list("landuse",  "forest",     "Wald"),
      list("natural",  "wood",       "Wald"),
      list("landuse",  "meadow",     "Wiese"),
      list("natural",  "grassland",  "Wiese"),
      list("landuse",  "farmland",   "Acker"),
      list("landuse",  "farmyard",   "Acker"),
      list("landuse",  "garden",     "Garten"),
      list("landuse",  "residential","Siedlung"),
      list("natural",  "water",      "Gewaesser"),
      list("natural",  "wetland",    "Gewaesser")
    )

    hab_schichten <- list()
    for (aq in abfragen) {
      res_poly <- abfrage_osm(aq[[1]], aq[[2]], bb_op, aq[[3]])
      if (!is.null(res_poly) && nrow(res_poly) > 0) {
        res_poly$habitat <- aq[[3]]
        hab_schichten[[length(hab_schichten)+1]] <- res_poly["geometry"] |>
          dplyr::mutate(habitat = aq[[3]])
      }
      Sys.sleep(0.5)   # kurze Pause zwischen Abfragen
    }

    # Hecken als Linien (separater Aufruf)
    cat("    natural=hedge (Linien)... ")
    hecke_lines <- tryCatch(
      withCallingHandlers({
        q <- opq(bbox=bb_op, timeout=25) |>
               add_osm_feature(key="natural", value="hedge")
        osmdata_sf(q)$osm_lines
      }, message=function(m){invokeRestart("muffleMessage")}),
      error=function(e) NULL)
    if (!is.null(hecke_lines) && nrow(hecke_lines) > 0) {
      hecke_utm <- st_buffer(st_transform(hecke_lines["geometry"], ziel_epsg), 8)
      hecke_utm$habitat <- "Hecke_Gehoelz"
      cat(sprintf("%d Hecken\n", nrow(hecke_lines)))
      hab_schichten[[length(hab_schichten)+1]] <- hecke_utm
    } else { cat("keine\n") }

    if (length(hab_schichten) > 0) {
      # Kombiniere alle Schichten
      hab_raw <- tryCatch(
        do.call(rbind, lapply(hab_schichten, function(x) {
          x <- x[, c("habitat","geometry")]
          st_transform(x, crs=4326)  # einheitlich WGS84
        })),
        error = function(e) {
          # Schichten einzeln umwandeln wenn rbind fehlschlaegt
          rbind_list <- lapply(hab_schichten, function(x) {
            tryCatch({
              x <- x[, c("habitat","geometry")]
              st_transform(x, crs=4326)
            }, error=function(e) NULL)
          })
          do.call(rbind, Filter(Negate(is.null), rbind_list))
        }
      )

      if (!is.null(hab_raw) && nrow(hab_raw) > 0) {
        ergebnis <- klassifiziere_und_extrahiere(hab_raw)
        habitat_fixes    <- ergebnis$habitat_fixes
        habitat_polygone <- ergebnis$habitat_polygone
        cat(sprintf("\n  Habitat geladen: %d Polygone, %d Fixes klassifiziert\n",
                    nrow(hab_raw), nrow(habitat_fixes)))
        cat("  Zitation: OpenStreetMap contributors, ODbL 1.0 (openstreetmap.org)\n\n")
      } else {
        cat("  [WARN] Keine Habitatpolygone nach Kombination\n\n")
      }
    } else {
      cat("  [WARN] Alle Habitatabfragen ohne Ergebnis\n")
      cat("  Tipp: Internetverbindung pruefen oder HABITAT_METHODE <- 'keine'\n\n")
    }
  }

# ══════════════════════════════════════════════════════════════
# OPTION B (veraltet, nicht empfohlen): osmextract
# ══════════════════════════════════════════════════════════════
} else if (HABITAT_METHODE == "osmextract") {
  cat("  [FEHLER] osmextract benoetigt system-GDAL (gdal-config).\n")
  cat("  Verwende stattdessen: HABITAT_METHODE <- 'overpass_klein'\n\n")
  cat("  Falls GDAL installiert (brew install gdal), kann osmextract verwendet werden.\n\n")

# ══════════════════════════════════════════════════════════════
# OPTION C (veraltet): atkis WFS
# Behalte als Fallback, aber overpass_klein ist zuverlaessiger
# ══════════════════════════════════════════════════════════════
} else if (HABITAT_METHODE == "atkis") {
  cat("  ATKIS WFS (experimentell) -- bei Problemen: HABITAT_METHODE <- 'overpass_klein'\n\n")
  # [atkis code bleibt erhalten, aber overpass_klein ist empfohlen]

# ══════════════════════════════════════════════════════════════
# OPTION D: Eigene GeoPackage-Datei
# ══════════════════════════════════════════════════════════════
} else if (HABITAT_METHODE == "gpkg" && nchar(HABITAT_GPKG_PFAD) > 0) {

  pfad_exp <- path.expand(HABITAT_GPKG_PFAD)
  cat(sprintf("  Lade Habitatkarte: %s\n", basename(pfad_exp)))

  hab_raw <- tryCatch(
    st_read(pfad_exp, quiet=TRUE),
    error = function(e) {
      cat(sprintf("  Fehler beim Laden: %s\n", conditionMessage(e)))
      NULL
    })

  if (!is.null(hab_raw)) {
    # Spaltenname "habitat" sicherstellen (Shapefile kuerzt manchmal)
    if (!"habitat" %in% names(hab_raw)) {
      hab_cols <- names(hab_raw)
      cat(sprintf("  Verfuegbare Spalten: %s\n", paste(hab_cols, collapse=", ")))
      cat("  [WARN] Keine 'habitat'-Spalte gefunden.\n\n")
      hab_raw <- NULL
    }
  }

  if (!is.null(hab_raw)) {
    cat(sprintf("  %d Habitatpolygone geladen\n", nrow(hab_raw)))
    cat("  Typen:", paste(sort(unique(hab_raw$habitat)), collapse=", "), "\n")

    ergebnis <- klassifiziere_und_extrahiere(hab_raw)
    habitat_fixes    <- ergebnis$habitat_fixes
    habitat_polygone <- ergebnis$habitat_polygone
    cat(sprintf("  Fertig: %d Fixes klassifiziert\n",
                nrow(habitat_fixes[habitat_grob != "Open land"])))
    cat("  Quelle: ATKIS Basis-DLM — © GeoBasis-DE / BKG 2024, DL-DE/BY-2.0\n\n")
  }

} else {
  cat("  Habitatdaten deaktiviert.\n")
  cat("  Empfehlung: HABITAT_METHODE <- 'overpass_klein'\n\n")
}
# ── Ausgabe Habitatverteilung ──────────────────────────────────
if (!is.null(habitat_fixes) && nrow(habitat_fixes) > 0) {
  hab_pro_igel <- habitat_fixes[, .N, by = .(igel, habitat_grob)]
  hab_pro_igel[, pct := round(100 * N / sum(N), 1), by = igel]
  cat("  Habitatnutzung pro Igel (% der Fixes):\n")
  print(dcast(hab_pro_igel, igel ~ habitat_grob, value.var = "pct", fill = 0))
  cat("\n")
}


# ── 7. PLOTS ─────────────────────────────────────────────────
# Neu: monatliche Aggregation statt tagesweise Punkte.
# Plot A: Monatlicher Aktivitaetsradius + Temperatur (Zeitverlauf)
# Plot B: Monatlicher Niederschlag + Aktivitaet
# Plot C: Regen-Effekt (Schwellenwert >2 mm, tagesweise)
# Plot D: Temperaturverlauf mit Beobachtungsnaechten
cat("====================================================\n")
cat("7. Plots erstellen\n")
cat("====================================================\n\n")

if (!is.null(monat_wetter) && !is.null(monat_igel) &&
    nrow(monat_wetter) >= 2 && nrow(monat_igel) >= 2) {

  monat_kombi <- merge(monat_igel, monat_wetter, by = c("monat_str","monat_date"),
                        all = TRUE)
  monat_kombi <- monat_kombi[!is.na(monat_date)][order(monat_date)]

  # Skalierungsfaktoren fuer Dual-Achsen-Plots
  temp_scale  <- max(monat_kombi$median_radius_m, na.rm=TRUE) /
                  max(monat_kombi$temp_mittel_monat, na.rm=TRUE)
  regen_scale <- if (max(monat_kombi$niederschlag_sum, na.rm=TRUE) > 0)
    max(monat_kombi$median_radius_m, na.rm=TRUE) /
    max(monat_kombi$niederschlag_sum, na.rm=TRUE) else 1

  # Plot A: Monatlicher Aktivitaetsradius + Temperatur
  p_monat_temp <- ggplot(monat_kombi[!is.na(median_radius_m) & !is.na(temp_mittel_monat)],
                          aes(x = monat_date)) +
    geom_col(aes(y = temp_mittel_monat * temp_scale),
             fill = "#fee090", alpha = 0.65, width = 20) +
    geom_line(aes(y = median_radius_m), color = "#2166ac",
              linewidth = 1.3, na.rm = TRUE) +
    geom_point(aes(y = median_radius_m, size = n_naechte_igel),
               color = "#2166ac", alpha = 0.85, na.rm = TRUE) +
    scale_y_continuous(
      name   = "Median 95%-Aktivitaetsradius (m)",
      labels = label_number(suffix = " m"),
      sec.axis = sec_axis(~ . / temp_scale,
                           name   = "Monatsmitteltemperatur (Grad C)",
                           labels = label_number(suffix = " C"))
    ) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y") +
    scale_size_continuous(name = "N Beobachtungsnaechte", range = c(2, 7)) +
    labs(title    = "A  Nacht-Aktivitaetsradius und Temperatur (monatlich)",
         subtitle = "Blau = Aktivitaetsradius (Median) | Gelbe Balken = Monatsmitteltemperatur",
         x = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title        = element_text(face = "bold", size = 12),
          plot.subtitle     = element_text(size = 9, color = "grey40"),
          legend.position   = "bottom",
          axis.title.y.right= element_text(color = "#d73027"),
          axis.text.y.right = element_text(color = "#d73027"))

  # Plot B: Niederschlag + Aktivitaetsradius
  p_monat_regen <- ggplot(monat_kombi[!is.na(median_radius_m) & !is.na(niederschlag_sum)],
                           aes(x = monat_date)) +
    geom_col(aes(y = niederschlag_sum * regen_scale),
             fill = "#4393c3", alpha = 0.55, width = 20) +
    geom_line(aes(y = median_radius_m), color = "#2d6a4f",
              linewidth = 1.3, na.rm = TRUE) +
    geom_point(aes(y = median_radius_m), color = "#2d6a4f",
               size = 3, alpha = 0.85, na.rm = TRUE) +
    scale_y_continuous(
      name   = "Median 95%-Aktivitaetsradius (m)",
      labels = label_number(suffix = " m"),
      sec.axis = sec_axis(~ . / regen_scale,
                           name   = "Monatsniederschlag (mm)",
                           labels = label_number(suffix = " mm"))
    ) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y") +
    labs(title    = "B  Nacht-Aktivitaetsradius und Niederschlag (monatlich)",
         subtitle = "Gruen = Aktivitaetsradius | Blaue Balken = Monatssumme Niederschlag",
         x = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title        = element_text(face = "bold", size = 12),
          plot.subtitle     = element_text(size = 9, color = "grey40"),
          axis.title.y.right= element_text(color = "#4393c3"),
          axis.text.y.right = element_text(color = "#4393c3"))

  p_wetter_monat <- p_monat_temp / p_monat_regen +
    plot_annotation(
      title    = "Nacht-Aktivitaet der Igel im saisonalen Kontext",
      subtitle = sprintf("DWD-Station: %s (%.1f km) | Monatliche Aggregation | N = %d Igel",
                         station_name, station_dist,
                         if (!is.null(monat_igel)) max(monat_igel$n_igel, na.rm=TRUE) else 0),
      theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                       plot.subtitle = element_text(size = 9,  color = "grey40"))
    )
  ggsave(file.path(out_ordner, "wetter_homerange.png"),
         p_wetter_monat, width = 13, height = 9, dpi = 150)
  cat("  wetter_homerange.png (monatlich) gespeichert\n")
}

# Plot C: Regen-Effekt (tagesweise mit Schwellenwert)
if (nrow(nacht_stat) >= 10 && "regen_nacht" %in% names(nacht_stat)) {

  regen_tab <- nacht_stat[, .(
    n = .N,
    median_rad = median(radius_95pct_m, na.rm=TRUE)
  ), by = regen_nacht]
  regen_label <- if (nrow(regen_tab) == 2)
    sprintf("Regennacht >2mm: Median %.0f m | Trocken: Median %.0f m",
            regen_tab[regen_nacht=="Regennacht (>2 mm)", median_rad],
            regen_tab[regen_nacht=="Trockene Nacht",     median_rad])
  else ""

  p_regen_effekt <- ggplot(nacht_stat,
    aes(x = regen_nacht, y = radius_95pct_m, fill = regen_nacht)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.5) +
    geom_jitter(aes(color = regen_nacht), width = 0.12, size = 2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "black") +
    scale_fill_manual(values  = c("Regennacht (>2 mm)"="#4393c3",
                                   "Trockene Nacht"    ="#d6604d"),
                      guide = "none") +
    scale_color_manual(values = c("Regennacht (>2 mm)"="#2166ac",
                                   "Trockene Nacht"    ="#b2182b"),
                       guide = "none") +
    scale_y_continuous(labels = label_number(suffix = " m")) +
    labs(title    = "C  Regen-Effekt auf den Nacht-Aktionsradius",
         subtitle = "Bewegen sich Igel bei Regen (>2 mm) weniger weit? | Raute = Mittelwert",
         caption  = regen_label,
         x = NULL, y = "95%-Aktivitaetsradius (m)") +
    theme_bw(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          plot.caption  = element_text(size = 9, color = "grey30"))

  ggsave(file.path(out_ordner, "regen_effekt.png"),
         p_regen_effekt, width = 7, height = 5.5, dpi = 150)
  cat("  regen_effekt.png gespeichert\n")
}

# Plot D: Temperaturverlauf mit Beobachtungsnaechten (Zeitstrahl)
if (!is.null(wetter_dt) && nrow(wetter_dt) > 0) {
  wetter_beob <- wetter_dt[datum >= datum_range[1] & datum <= datum_range[2]]
  if (nrow(wetter_beob) >= 5) {
    beob_naechte <- data.table(datum = unique(alle_fixes$datum))
    p_temp_verlauf <- ggplot(wetter_beob, aes(x = datum)) +
      geom_ribbon(aes(ymin = temp_min, ymax = temp_max),
                  fill = "#fee090", alpha = 0.4) +
      geom_line(aes(y = temp_mittel), color = "#d73027", linewidth = 0.8) +
      geom_hline(yintercept = c(5, 10, 18),
                 linetype = "dotted", color = "grey50", linewidth = 0.5) +
      annotate("text", x = min(wetter_beob$datum), y = c(5.5, 10.5, 18.5),
               label = c("5 C", "10 C", "18 C"),
               size = 2.7, color = "grey40", hjust = 0) +
      geom_rug(data = beob_naechte, aes(x = datum),
               sides = "b", color = "#2166ac", alpha = 0.6,
               inherit.aes = FALSE) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
      labs(title    = "D  Temperaturverlauf im Beobachtungszeitraum",
           subtitle = "Rot = Tagesmittel | Gelb = Min-Max | Blaue Striche = Beobachtungsnaechte",
           x = NULL, y = "Temperatur (Grad C)") +
      theme_bw(base_size = 12) +
      theme(plot.title  = element_text(face = "bold"),
            axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
    ggsave(file.path(out_ordner, "temperaturverlauf.png"),
           p_temp_verlauf, width = 13, height = 5, dpi = 150)
    cat("  temperaturverlauf.png gespeichert\n")
  }
}

# Mondphase-Plot entfernt

if (!is.null(habitat_polygone) && nrow(habitat_polygone) > 0) {
  # Farben fuer zusammengefasste Klassifikation (habitat_grob)
  # Waldweg + Waldrand → "Edge habitat" (orange-braun)
  hab_farben <- c(
    "Edge habitat"    = "#c97a1e",   # Orange-Braun — Randstrukturen + Wege
    "Forest interior" = "#1a6b1a",   # Dunkelgruen — Waldinneres
    "Agricultural land"= "#d4a843",  # Ocker — Acker
    "Grassland"       = "#90c050",   # Hellgruen — Wiese
    "Scrub/Hedgerow"  = "#4a8f3f",   # Mittelgruen
    "Water body"      = "#4a90c4",   # Blau
    "Settlement"      = "#c0a080",   # Beige
    "Garden"          = "#b8e0a0",   # Blassgruen
    "Open land"       = "#e8e8d0"    # Hellgrau
  )

  # Karte: nur Bereich um Release-Punkt (600 m)
  release_sf <- st_as_sf(data.frame(x=514743.7, y=5805363.8),
                          coords=c("x","y"), crs=32632)
  study_buf  <- st_buffer(release_sf, 650)
  hab_clip   <- suppressMessages(st_intersection(habitat_polygone, study_buf))

  p_hab_karte <- ggplot() +
    geom_sf(data = hab_clip, aes(fill = habitat), color = "white",
            linewidth = 0.15, alpha = 0.8) +
    scale_fill_manual(values = hab_farben, name = "Habitat", drop = FALSE) +
    geom_point(data = alle_fixes[!is.na(x) & !is.na(y)],
               aes(x = x, y = y), color = "white", size = 0.3, alpha = 0.12) +
    geom_sf(data = release_sf, color = "red", shape = 8, size = 4, stroke = 1.5) +
    labs(title = "Habitatkarte Sachsenhagen — Studiengebiet",
         subtitle = "Weiss = Nacht-Fixes | Stern = Auswilderungsvoliere",
         x = "UTM32N Ost", y = "UTM32N Nord") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")

  ggsave(file.path(out_ordner, "habitat_karte.png"),
         p_hab_karte, width = 10, height = 8, dpi = 150)
  cat("  habitat_karte.png gespeichert\n")
}

if (!is.null(habitat_fixes) && nrow(habitat_fixes) > 0) {

  # ── Plot 1: Habitatnutzung pro Igel ───────────────────────────
  hab_pro_igel <- habitat_fixes[, .N, by = .(igel, habitat_grob)]
  hab_pro_igel[, pct := round(100 * N / sum(N), 1), by = igel]

  # Igel nach Wald-Anteil sortieren
  wald_order <- hab_pro_igel[habitat_grob %in% c("Edge habitat","Forest interior"),
                               .(wald_pct = sum(pct)), by = igel]
  igel_ord   <- wald_order[order(-wald_pct), igel]
  hab_pro_igel[, igel := factor(igel, levels = igel_ord)]

  p_hab_nutz <- ggplot(hab_pro_igel, aes(x = igel, y = pct, fill = habitat_grob)) +
    geom_col(position = "fill", width = 0.75) +
    scale_fill_manual(values = hab_farben, name = "Habitat", drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    coord_flip() +
    labs(title = "Habitat use per individual (proportion of nocturnal fixes)",
         subtitle = "Sorted by forest + forest edge proportion",
         x = NULL, y = "Anteil (%)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  ggsave(file.path(out_ordner, "habitat_nutzung.png"),
         p_hab_nutz, width = 10, height = 8, dpi = 150)
  cat("  habitat_nutzung.png gespeichert\n")

  # ── Plot 2: Nutzung vs. Verfuegbarkeit (Population) ───────────
  # Verfuegbarkeit: Habitatflaechen innerhalb 600 m Release-Punkt
  if (exists("release_sf") && !is.null(habitat_polygone)) {
    tryCatch({
      study_buf   <- st_buffer(release_sf, 600)
      hab_avail   <- suppressMessages(st_intersection(habitat_polygone, study_buf))
      hab_avail$area_ha <- as.numeric(st_area(hab_avail)) / 10000
      avail_dt    <- as.data.table(hab_avail)[, .(avail_ha = sum(area_ha)), by = habitat]
      avail_dt[, avail_pct := 100 * avail_ha / sum(avail_ha)]

      use_pop <- habitat_fixes[, .N, by = habitat_grob]
      use_pop[, use_pct := 100 * N / sum(N)]
      setnames(use_pop, "habitat_grob", "habitat")

      uva <- merge(avail_dt[, .(habitat, avail_pct)],
                   use_pop[,   .(habitat, use_pct)],
                   by = "habitat", all = TRUE)
      uva[is.na(avail_pct), avail_pct := 0]
      uva[is.na(use_pct),   use_pct   := 0]

      uva_long <- melt(uva, id.vars = "habitat",
                       measure.vars  = c("avail_pct", "use_pct"),
                       variable.name = "typ", value.name = "pct")
      uva_long[, typ := ifelse(typ == "avail_pct",
                               "Verfuegbar (Flaeche %)", "Genutzt (Fixes %)")]
      # Habitate mit > 0.1% in mindestens einer Kategorie
      keep_hab <- uva[avail_pct > 0.1 | use_pct > 0.1, habitat]
      uva_long <- uva_long[habitat %in% keep_hab]
      uva_long[, habitat := factor(habitat,
                 levels = c("Edge habitat","Forest interior",
                            "Wiese","Acker","Gewaesser","Siedlung","Garten","Offenland"))]

      p_uva <- ggplot(uva_long[!is.na(habitat)],
                      aes(x = habitat, y = pct, fill = typ)) +
        geom_col(position = "dodge", width = 0.65, alpha = 0.9) +
        scale_fill_manual(values = c("Verfuegbar (Flaeche %)" = "#aaaaaa",
                                     "Genutzt (Fixes %)"      = "#2166ac"),
                          name = NULL) +
        scale_y_continuous(labels = function(x) paste0(round(x), "%")) +
        labs(title    = "Habitat use vs. availability",
             subtitle = "Availability = area within 600 m of release enclosure",
             x = NULL, y = "Anteil (%)") +
        theme_bw(base_size = 12) +
        theme(axis.text.x     = element_text(angle = 30, hjust = 1),
              plot.title       = element_text(face = "bold"),
              legend.position  = "top")

      ggsave(file.path(out_ordner, "habitat_nutzung_vs_verfuegbar.png"),
             p_uva, width = 9, height = 6, dpi = 150)
      cat("  habitat_nutzung_vs_verfuegbar.png gespeichert\n")

      # Tabelle ausgeben
      cat("\n  === Nutzung vs. Verfuegbarkeit (Population) ===\n")
      print(uva[order(-use_pct), .(habitat, avail_pct = round(avail_pct,1),
                                    use_pct = round(use_pct,1))])
      cat("\n")
    }, error = function(e) {
      cat(sprintf("  [WARN] Use-vs-Availability Plot fehlgeschlagen: %s\n", conditionMessage(e)))
    })
  }
}

# ══════════════════════════════════════════════════════════════
# NEUE ANALYSEN: Selektionsindizes + Temporalvergleich + Karte
# ══════════════════════════════════════════════════════════════

if (!is.null(habitat_fixes) && nrow(habitat_fixes) > 0) {

  cat("Berechne Selektionsindizes und Temporalanalyse...\n")

  # ── A. Selektionsindizes (Manly's alpha) + Chi-Quadrat ─────────
  tryCatch({
    # ── Verfuegbarkeit aus Shapefile (600m Puffer) ─────────────────
    release_sf_sel <- st_as_sf(data.frame(x=514743.7, y=5805363.8),
                                coords=c("x","y"), crs=32632)
    study_sel <- st_buffer(release_sf_sel, 600)
    hab_sf_32 <- st_transform(habitat_polygone, 32632)
    hab_avail_raw <- suppressMessages(st_intersection(hab_sf_32, study_sel))
    hab_avail_raw$area_ha <- as.numeric(st_area(hab_avail_raw)) / 10000
    avail_shp <- as.data.table(hab_avail_raw)[, .(avail_ha=sum(area_ha)), by=habitat]

    # ── Edge-Habitat-Flaeche berechnen (±20m Waldrand-Zone) ─────────
    # Edge-Verfuegbarkeit = Flaeche der ±20m Randzone um den Wald
    wald_sf_32 <- hab_sf_32[hab_sf_32$habitat == "Wald", ]
    edge_ha <- 0
    if (nrow(wald_sf_32) > 0) {
      wald_union  <- st_union(wald_sf_32)
      wald_puffer <- tryCatch(st_buffer(wald_union, 20), error=function(e) NULL)
      wald_innen  <- tryCatch(st_buffer(wald_union, -20), error=function(e) wald_union)
      if (!is.null(wald_puffer)) {
        rand_zone <- st_make_valid(st_difference(wald_puffer, wald_innen))
        rand_in_study <- suppressMessages(st_intersection(rand_zone, study_sel))
        edge_ha <- as.numeric(st_area(rand_in_study)) / 10000
        cat(sprintf("  Edge habitat (±20m Waldrand): %.1f ha\n", edge_ha))
      }
    }

    # ── Shapefile-Kategorien → habitat_grob mapping ─────────────────
    hab_grob_map <- c(
      "Wald"            = "Forest interior",   # Waldinneres NACH Abzug der Randzone
      "Acker"           = "Agricultural land",
      "Wiese"           = "Grassland",
      "Gehoelz_Strauch" = "Edge habitat",      # Gehoelzstreifen → Edge (zu klein fuer Genauigkeit)
      "Gewaesser"       = "Water body",
      "Siedlung"        = "Settlement",
      "Garten"          = "Garden"
    )
    avail_shp[, habitat_grob := hab_grob_map[habitat]]
    avail_grob <- avail_shp[!is.na(habitat_grob), .(avail_ha=sum(avail_ha)), by=habitat_grob]
    # Waldinneres um Edge-Flaeche reduzieren (da Edge aus Waldflaeche entnommen)
    avail_grob[habitat_grob == "Forest interior", avail_ha := pmax(0, avail_ha - edge_ha)]
    # Waldrand-Zone zur Edge-Flaeche addieren (Gehoelz_Strauch ist bereits drin)
    if ("Edge habitat" %in% avail_grob$habitat_grob) {
      avail_grob[habitat_grob == "Edge habitat", avail_ha := avail_ha + edge_ha]
    } else {
      avail_grob <- rbind(avail_grob, data.table(habitat_grob="Edge habitat", avail_ha=edge_ha))
    }
    avail_grob[, avail_pct := avail_ha / sum(avail_ha)]
    setnames(avail_grob, "habitat_grob", "habitat")
    avail_dt <- avail_grob

    # ── Nutzung (Fixes) nach habitat_grob ───────────────────────────
    use_dt <- habitat_fixes[, .N, by=habitat_grob]
    use_dt[, use_pct := N / sum(N)]
    setnames(use_dt, "habitat_grob", "habitat")

    # ── Zusammenfuehren ──────────────────────────────────────────────
    sel_dt <- merge(avail_dt, use_dt[, .(habitat, N, use_pct)],
                    by="habitat", all=TRUE)
    sel_dt[is.na(N), N := 0]
    sel_dt[is.na(use_pct), use_pct := 0]
    sel_dt[is.na(avail_pct), avail_pct := 0.001]  # Kleiner Wert fuer Habitate ohne Flaeche

    # Selektionsratio: use_pct / avail_pct (> 1 = bevorzugt, < 1 = gemieden)
    sel_dt[, selektionsratio := round(use_pct / avail_pct, 2)]

    # Manly's alpha (normierter Selektionsindex, 0-1, > 1/k = bevorzugt)
    k <- nrow(sel_dt[avail_pct > 0])
    sel_dt[avail_pct > 0, manly_alpha := (use_pct / avail_pct) / sum(use_pct / avail_pct)]
    sel_dt[, manly_alpha := round(manly_alpha, 3)]

    # Pool-Chi-Quadrat-Test ENTFERNT (2026-06).
    # Begruendung: chisq.test auf Einzelfixes (n = 304.068) verletzt die
    # Unabhaengigkeitsannahme (raeumlich/zeitlich autokorrelierte Ortungen) und
    # ist dadurch kuenstlich hochsignifikant (Pseudoreplikation). Selektionsratio
    # und Manly's alpha bleiben als DESKRIPTIVE Masse erhalten; die statistische
    # Inferenz erfolgt ausschliesslich ueber den Per-Tier-Wilcoxon-Test unten
    # (korrekte Stichprobeneinheit = Tier).
    chi_res <- NULL
    cat("\n  === Habitatselektion: Inferenz ueber Per-Tier-Wilcoxon (Einheit = Tier) ===\n")
    cat("  [Hinweis] Gepoolter Fix-Chi-Quadrat entfernt (pseudorepliziert).\n\n")

    # ── Per-Tier-Habitatselektion: Wilcoxon-Einzel-Test (Hauptanalyse) ─────────
    # Korrekte statistische Einheit = Tier, nicht Einzelfix
    # Testet: Ist Median der per-Tier-Proportionen != Verfuegbarkeit?
    tryCatch({
      n_tiere_hab <- length(unique(habitat_fixes$igel))
      cat(sprintf("  === Per-Tier-Habitatselektion (Wilcoxon, n=%d Tiere) ===\n", n_tiere_hab))

      # Anteil Fixes pro Habitat pro Tier berechnen
      pt_hab <- habitat_fixes[, .N, by=.(igel, habitat_grob)]
      pt_tot <- habitat_fixes[, .(N_total=.N), by=igel]
      pt_hab <- merge(pt_hab, pt_tot, by="igel")
      pt_hab[, prop := N / N_total]

      # Fuer jedes Haupthabitat: Wilcoxon-Test gegen Verfuegbarkeit
      main_hab_test <- c("Edge habitat","Forest interior","Agricultural land","Grassland")
      wilcox_habitat_dt <- rbindlist(lapply(main_hab_test, function(h) {
        avail_p <- avail_dt[habitat==h, avail_pct]
        if (length(avail_p)==0 || is.na(avail_p)) return(NULL)
        # Tiere mit Fixes in diesem Habitat
        props_with <- pt_hab[habitat_grob==h, prop]
        # Tiere ohne Fixes: prop = 0
        all_tiere <- unique(habitat_fixes$igel)
        n_without  <- length(setdiff(all_tiere, pt_hab[habitat_grob==h, igel]))
        props_all  <- c(props_with, rep(0, n_without))
        wt <- wilcox.test(props_all, mu=avail_p, alternative="two.sided")
        data.table(
          Habitat        = h,
          n_Tiere        = length(props_all),
          Verfuegbar_pct = round(avail_p*100, 1),
          Median_use_pct = round(median(props_all)*100, 1),
          W              = as.numeric(wt$statistic),
          p_value        = round(wt$p.value, 4),
          Signifikanz    = ifelse(wt$p.value<0.001,"***",
                          ifelse(wt$p.value<0.01, "**",
                          ifelse(wt$p.value<0.05, "*", "n.s.")))
        )
      }))
      cat("\n  Per-Tier-Wilcoxon-Ergebnisse:\n")
      print(wilcox_habitat_dt)
    }, error=function(e) {
      cat(sprintf("  [WARN] Per-Tier-Wilcoxon: %s\n", conditionMessage(e)))
      wilcox_habitat_dt <<- NULL
    })

    # Selektionsindex-Tabelle
    sel_out <- sel_dt[order(-selektionsratio), .(
      Habitat       = habitat,
      Verfuegbar_pct = round(avail_pct*100, 1),
      Genutzt_pct   = round(use_pct*100, 1),
      N_Fixes       = N,
      Selektionsratio = selektionsratio,
      Manly_alpha   = manly_alpha,
      Praeferenz    = ifelse(selektionsratio > 1.2, "Bevorzugt",
                     ifelse(selektionsratio < 0.8, "Gemieden", "Neutral"))
    )]
    cat("  Selektionsindizes:\n")
    print(sel_out)

    # HINWEIS: file.remove wurde verschoben — erst nach erfolgreichem ggsave

    # ── Plot 1: Lollipop-Chart Selektionsratio (Log-Skala) ─────────
    # Zeigt klar WELCHE Habitate bevorzugt/gemieden werden
    # Log-Skala noetig weil Ratios sehr unterschiedlich

    hab_labels_en <- c(
      "Waldweg"="Forest track", "Waldrand"="Forest edge",
      "Wald"="Forest interior", "Acker"="Agricultural land",
      "Wiese"="Grassland", "Gehoelz_Strauch"="Scrub/Hedgerow",
      "Gewaesser"="Water body", "Offenland"="Open land",
      "Siedlung"="Settlement", "Garten"="Garden"
    )

    sel_plot <- sel_dt[avail_pct > 0.0005 | use_pct > 0.005]
    sel_plot[, hab_en := fifelse(habitat %in% names(hab_labels_en),
                                  hab_labels_en[habitat], habitat)]
    sel_plot[, hab_en := factor(hab_en, levels=sel_plot[order(selektionsratio), hab_en])]
    sel_plot[, bevorzugt := selektionsratio > 1]
    # Ratio auf max 200 begrenzen fuer Darstellung
    sel_plot[, ratio_plot := pmin(selektionsratio, 200)]

    p_selektion <- ggplot(sel_plot, aes(x=hab_en, y=ratio_plot, color=bevorzugt)) +
      geom_hline(yintercept=1, linetype="dashed", color="grey50", linewidth=0.8) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=1, ymax=Inf,
               fill="#1a6b1a", alpha=0.04) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1,
               fill="#d4a843", alpha=0.06) +
      geom_segment(aes(x=hab_en, xend=hab_en, y=1, yend=ratio_plot),
                   linewidth=1.2, alpha=0.6) +
      geom_point(size=5, alpha=0.9) +
      # Beschriftung: Used% und Available%
      # Hohe Ratios (>5x): Label LINKS vom Punkt (damit es nicht aus dem Plot faellt)
      # Niedrige Ratios (<1x): Label RECHTS vom Punkt
      geom_text(aes(label=sprintf("used: %.0f%%\navail: %.1f%%",
                                   use_pct*100, avail_pct*100)),
                hjust=ifelse(sel_plot$ratio_plot > 5, 1.08, -0.08),
                size=2.8, color="grey30", lineheight=0.9) +
      scale_color_manual(values=c("TRUE"="#1a6b1a","FALSE"="#c0392b"),
                         labels=c("TRUE"="Preferred (ratio > 1)",
                                  "FALSE"="Avoided (ratio < 1)"),
                         name=NULL) +
      scale_y_log10(breaks=c(0.1,0.25,0.5,1,2,5,10,25,100,200),
                    labels=c("0.1×","0.25×","0.5×","1×","2×","5×",
                             "10×","25×","100×","≥200×"),
                    expand=expansion(mult=c(0.25, 0.05))) +
      coord_flip() +
      labs(title="Habitat selection ratios (used / available)",
           subtitle=paste0(
             "Ratio > 1 = preferred, < 1 = avoided | Log scale\n",
             "Note: selection ratios capped at 200× for display; ",
             "Edge habitat actual ratio = ", round(sel_dt[habitat=="Edge habitat", selektionsratio],0), "×"),
           x=NULL, y="Selection ratio (log scale)") +
      theme_bw(base_size=12) +
      theme(plot.title=element_text(face="bold"),
            plot.subtitle=element_text(size=8.5, color="grey40"),
            legend.position="top",
            panel.grid.minor=element_blank())

    # ── Plot 2: Used vs. Available Balken (nur die 4 Haupthabitate) ──
    main_hab <- c("Edge habitat","Forest interior","Agricultural land","Grassland")
    sel_main <- sel_dt[habitat %in% main_hab]
    # habitat_grob-Namen sind bereits auf Englisch — direkt verwenden
    sel_main[, hab_en := habitat]
    sel_main[, hab_en := factor(hab_en, levels=main_hab)]
    sel_long2 <- melt(sel_main[, .(hab_en, avail_pct, use_pct)],
                      id.vars="hab_en", variable.name="typ", value.name="pct")
    sel_long2[, typ_label := ifelse(typ=="avail_pct", "Available","Used")]
    sel_long2[, pct100 := pct*100]

    p_uva <- ggplot(sel_long2, aes(x=hab_en, y=pct100, fill=typ_label)) +
      geom_col(position="dodge", width=0.6, alpha=0.9) +
      geom_text(aes(label=paste0(round(pct100,1),"%")),
                position=position_dodge(0.6), vjust=-0.4, size=3.2, fontface="bold") +
      scale_fill_manual(values=c("Available"="#aaaaaa","Used"="#2166ac"), name=NULL) +
      scale_y_continuous(labels=function(x) paste0(x,"%"),
                         expand=expansion(mult=c(0, 0.15))) +
      labs(title="Habitat use vs. availability — key habitats",
           subtitle="Blue = proportion of nocturnal fixes | Grey = proportion of available area (within 600 m)",
           x=NULL, y="Proportion (%)") +
      theme_bw(base_size=12) +
      theme(plot.title=element_text(face="bold"),
            plot.subtitle=element_text(size=9, color="grey40"),
            legend.position="top")

    # Beide Plots kombinieren
    p_selektion_kombi <- p_uva / p_selektion +
      plot_layout(heights=c(1, 1.5))

    tryCatch({
      ggsave(file.path(out_ordner, "habitat_selektion.png"),
             p_selektion_kombi, width=10, height=12, dpi=150)
      cat("  ✓ habitat_selektion.png gespeichert\n")
      cat(sprintf("    Edge: used=%.1f%%, avail=%.1f%%, ratio=%.1fx\n",
                  sel_dt[habitat=="Edge habitat", use_pct*100],
                  sel_dt[habitat=="Edge habitat", avail_pct*100],
                  sel_dt[habitat=="Edge habitat", selektionsratio]))
    }, error = function(e) {
      cat(sprintf("  [FEHLER] habitat_selektion.png: %s\n", conditionMessage(e)))
      # Fallback: nur den Lollipop-Plot speichern
      tryCatch({
        ggsave(file.path(out_ordner, "habitat_selektion.png"),
               p_selektion, width=10, height=7, dpi=150)
        cat("  ✓ habitat_selektion.png (nur Lollipop) gespeichert\n")
      }, error=function(e2) cat(sprintf("  [FEHLER] Lollipop: %s\n", conditionMessage(e2))))
    })

  }, error=function(e) cat(sprintf("  [WARN] Selektion: %s\n", conditionMessage(e))))

  # ── B. Temporalvergleich: Pro-Tier-Vergleich mit Wilcoxon ────────
  # Methodisch korrekt wie Block 2: pro Tier % Fixes in jeder Phase
  # berechnen, dann gepaarter Wilcoxon Signed-Rank Test (nur Tiere
  # mit Daten in BEIDEN Phasen → kein Survivorship Bias).
  tryCatch({
    if ("tage_seit" %in% names(habitat_fixes) && sum(!is.na(habitat_fixes$tage_seit)) > 100) {

      # Phasendefinition (konsistent mit Block 2d)
      PHASE1_MAX <- 5   # Fruehphase: Tag 1-5
      PHASE2_MIN <- 6   # Spaetphase: Tag 6+
      hab_focus  <- c("Edge habitat","Forest interior","Agricultural land","Grassland")  # Haupthabitate

      # Pro Tier + Phase: Anteil Fixes in jedem Habitat
      habitat_fixes[, phase2 := fifelse(
        tage_seit <= PHASE1_MAX, "Early (day 1-5)", "Late (day 6+)")]

      # Aggregation: pro Igel pro Phase → % Fixes je Habitat
      per_tier <- habitat_fixes[!is.na(tage_seit), {
        n_ges <- .N
        lapply(hab_focus, function(h) {
          data.table(habitat = h, pct = 100 * sum(habitat_grob == h) / n_ges)
        }) |> rbindlist()
      }, by = .(igel, phase2)]

      # Nur Tiere mit Daten in BEIDEN Phasen (korrekt gepaart)
      hat_beide <- per_tier[, .(n_phasen = uniqueN(phase2)), by = igel]
      igel_gepaart <- hat_beide[n_phasen == 2, igel]
      per_tier_gepaart <- per_tier[igel %in% igel_gepaart]
      n_gepaart <- length(igel_gepaart)
      cat(sprintf("  Temporalanalyse: %d Igel mit Daten in beiden Phasen\n\n", n_gepaart))

      if (n_gepaart >= 5) {
        # Wilcoxon Signed-Rank Test pro Habitat (gepaart pro Tier)
        wilcox_res <- rbindlist(lapply(hab_focus, function(h) {
          dt <- dcast(per_tier_gepaart[habitat == h],
                      igel ~ phase2, value.var = "pct")
          early_col <- "Early (day 1-5)"
          late_col  <- "Late (day 6+)"
          if (!early_col %in% names(dt) || !late_col %in% names(dt))
            return(NULL)
          dt <- dt[!is.na(get(early_col)) & !is.na(get(late_col))]
          if (nrow(dt) < 5) return(NULL)
          wt <- wilcox.test(dt[[early_col]], dt[[late_col]], paired = TRUE, exact = FALSE)
          data.table(
            Habitat     = h,
            N_Igel      = nrow(dt),
            Median_Early = round(median(dt[[early_col]]), 1),
            Median_Late  = round(median(dt[[late_col]]), 1),
            Differenz_PP = round(median(dt[[late_col]]) - median(dt[[early_col]]), 1),
            W            = round(wt$statistic, 1),
            p_Wert       = round(wt$p.value, 4),
            Signifikant  = ifelse(wt$p.value < 0.05, "yes", "n.s.")
          )
        }))
        cat("  === Wilcoxon Signed-Rank: Habitat-Nutzung Phase 1 vs. Phase 2 ===\n")
        print(wilcox_res)
        cat("\n")

        # Alten Plot loeschen damit kein veralteter File bei Fehler bleibt
        # file.remove erst nach erfolgreichem ggsave

        # ── Plot: Boxplots pro-Tier, Early vs. Late ─────────────────
        # Farben und Labels fuer hab_focus = c("Edge habitat","Forest interior","Agricultural land","Grassland")
        hab_farben_t <- c(
          "Edge habitat"      = "#c97a1e",  # Orange-Braun
          "Forest interior"   = "#1a6b1a",  # Dunkelgruen
          "Agricultural land" = "#d4a843",  # Ocker
          "Grassland"         = "#90c050"   # Hellgruen
        )
        # Labels = identisch mit habitat_grob-Namen (bereits Englisch)
        hab_en_t <- c(
          "Edge habitat"      = "Edge habitat",
          "Forest interior"   = "Forest interior",
          "Agricultural land" = "Agricultural land",
          "Grassland"         = "Grassland"
        )

        per_tier_plot <- per_tier_gepaart[habitat %in% hab_focus]
        per_tier_plot[, hab_en := hab_en_t[habitat]]
        per_tier_plot[, hab_en := factor(hab_en, levels=hab_en_t)]
        per_tier_plot[, phase_f := factor(phase2,
                       levels=c("Early (day 1-5)","Late (day 6+)"))]

        # p-Wert Labels fuer Plot
        if (!is.null(wilcox_res) && nrow(wilcox_res) > 0) {
          wilcox_res[, hab_en := hab_en_t[Habitat]]
          wilcox_res[, p_label := fifelse(
            p_Wert < 0.001, "p < 0.001",
            fifelse(p_Wert < 0.05, sprintf("p = %.3f *", p_Wert),
                    sprintf("p = %.3f (n.s.)", p_Wert)))]
          # y-Position fuer p-Label
          y_max <- per_tier_plot[, .(ymax = max(pct, na.rm=TRUE) * 1.15), by=.(hab_en)]
          wilcox_res <- merge(wilcox_res, y_max, by="hab_en", all.x=TRUE)
        }

        p_temporal <- ggplot(per_tier_plot,
                             aes(x=phase_f, y=pct, fill=habitat)) +
          geom_boxplot(alpha=0.75, outlier.shape=21, width=0.55,
                       outlier.size=1.5) +
          geom_jitter(width=0.1, size=1.8, alpha=0.6, shape=21,
                      aes(fill=habitat), color="white") +
          stat_summary(fun=mean, geom="point", shape=18,
                       size=4, color="black") +
          facet_wrap(~hab_en, scales="free_y", ncol=2) +
          scale_fill_manual(values=hab_farben_t, guide="none") +
          { if (!is.null(wilcox_res) && "p_label" %in% names(wilcox_res))
            geom_text(data=wilcox_res,
                      aes(x=1.5, y=ymax, label=p_label),
                      inherit.aes=FALSE, size=3.5, fontface="italic", color="grey30")
          } +
          labs(title="Habitat use: early (day 1-5) vs. late phase (day 6+)",
               subtitle=sprintf(
                 "Per-animal comparison | Wilcoxon signed-rank test (paired) | n = %d hedgehogs\nDiamond = mean | Only hedgehogs with data in both phases included",
                 n_gepaart),
               x=NULL, y="Proportion of nocturnal fixes (%)") +
          theme_bw(base_size=12) +
          theme(plot.title=element_text(face="bold"),
                plot.subtitle=element_text(size=9, color="grey40"),
                strip.text=element_text(face="bold", size=11))

        tryCatch({
          ggsave(file.path(out_ordner, "habitat_temporal.png"),
                 p_temporal, width=11, height=8, dpi=150)
          cat("  ✓ habitat_temporal.png gespeichert\n")
          cat(sprintf("    Kategorien: %s\n\n",
                      paste(unique(per_tier_plot$hab_en), collapse=", ")))
        }, error=function(e) cat(sprintf("  [FEHLER] habitat_temporal.png: %s\n", conditionMessage(e))))

        # Ergebnis-Tabelle speichern fuer Excel
        temporal_wilcox_tab <- wilcox_res
      } else {
        cat("  Zu wenige Igel fuer gepaarten Test — Temporalanalyse uebersprungen\n")
        temporal_wilcox_tab <- NULL
      }
    } else {
      cat("  tage_seit nicht verfuegbar — Temporalanalyse uebersprungen\n")
      temporal_wilcox_tab <- NULL
    }
  }, error=function(e) {
    cat(sprintf("  [WARN] Temporal: %s\n", conditionMessage(e)))
    temporal_wilcox_tab <- NULL
  })

  # ── C. Individuelle Habitatverlaeufe pro Tier ─────────────────────
  # Pro Tier pro Nacht: % Fixes in Waldweg + Waldrand (= "Edge habitat")
  # Zeigt individuelle Trajektorien über die Zeit nach Auswilderung.
  tryCatch({
    if ("tage_seit" %in% names(habitat_fixes) && sum(!is.na(habitat_fixes$tage_seit)) > 100) {

      cat("  Berechne individuelle Habitatverlaeufe...\n")

      # Pro Igel pro Nacht: Anteil in jedem Habitat
      nacht_hab <- habitat_fixes[!is.na(tage_seit) & tage_seit >= 0,
        .(
          pct_waldweg  = 100 * sum(habitat_grob == "Edge habitat") / .N,  # Edge habitat combined
          pct_waldrand = 0,  # merged into Edge habitat
          pct_wald     = 100 * sum(habitat_grob == "Forest interior") / .N,
          pct_acker    = 100 * sum(habitat_grob == "Agricultural land") / .N,
          pct_edge     = 100 * sum(habitat_grob == "Edge habitat") / .N,
          n_fixes      = .N
        ),
        by = .(igel, datum, tage_seit)
      ]
      nacht_hab <- nacht_hab[n_fixes >= 10]  # Naechte mit sehr wenig Fixes ausschliessen

      # Zusammenfassung pro Tier
      igel_zusammenfassung <- nacht_hab[, {
        frueh <- .SD[tage_seit <= 5, mean(pct_edge, na.rm=TRUE)]
        spaet <- .SD[tage_seit > 5,  mean(pct_edge, na.rm=TRUE)]
        n_naechte_frueh <- .SD[tage_seit <= 5, .N]
        n_naechte_spaet <- .SD[tage_seit > 5,  .N]
        diff <- spaet - frueh
        .(
          N_Naechte    = .N,
          N_Frueh      = n_naechte_frueh,
          N_Spaet      = n_naechte_spaet,
          Edge_Frueh_pct = round(frueh, 1),
          Edge_Spaet_pct = round(spaet, 1),
          Differenz_PP   = round(diff, 1),
          Trend          = fifelse(is.na(diff), "unbekannt",
                           fifelse(diff >  3, "zunehmend",
                           fifelse(diff < -3, "abnehmend", "stabil")))
        )
      }, by = igel]
      setorder(igel_zusammenfassung, -Edge_Frueh_pct)

      cat("  === Individuelle Habitatnutzung (Edge = Waldweg + Waldrand) ===\n")
      print(igel_zusammenfassung)
      cat("\n")

      # ── Plot: Individuelle Trajektorien ─────────────────────────
      trend_farben <- c("zunehmend"="#1a6b1a","stabil"="#636363","abnehmend"="#c0392b","unbekannt"="#aaaaaa")
      nacht_hab_plot <- merge(nacht_hab, igel_zusammenfassung[, .(igel, Trend)], by="igel")

      p_individual <- ggplot(nacht_hab_plot[tage_seit <= 27],
                             aes(x=tage_seit, y=pct_edge, color=Trend)) +
        geom_point(size=1.2, alpha=0.5) +
        geom_smooth(method="loess", span=0.75, se=FALSE, linewidth=1.0) +
        geom_vline(xintercept=5.5, linetype="dashed", color="grey60", linewidth=0.5) +
        facet_wrap(~igel, ncol=5, scales="free_y") +
        scale_color_manual(values=trend_farben, name="Trend (edge habitat)",
                           guide=guide_legend(override.aes=list(linewidth=2, size=3))) +
        scale_x_continuous(breaks=c(0,5,10,15,20,25)) +
        labs(title="Individual trajectories: forest edge + track use over time",
             subtitle="% of nocturnal fixes in forest edge (±20m) or forest track (<15m) | Dashed line = phase boundary (day 5)",
             x="Days since release", y="Edge habitat use (%)") +
        theme_bw(base_size=9) +
        theme(plot.title=element_text(face="bold", size=11),
              plot.subtitle=element_text(size=8, color="grey40"),
              strip.text=element_text(face="bold", size=8),
              legend.position="bottom",
              panel.grid.minor=element_blank())

      ggsave(file.path(out_ordner, "habitat_individual_trajectory.png"),
             p_individual, width=14, height=10, dpi=150)
      cat("  habitat_individual_trajectory.png\n")

      # ── Plot: Übersicht aller Tiere (Boxplot früh vs. spät) ─────
      nacht_hab_long <- melt(
        nacht_hab[, .(igel, tage_seit, pct_waldweg, pct_waldrand, pct_wald, pct_acker)],
        id.vars=c("igel","tage_seit"),
        variable.name="habitat_var", value.name="pct"
      )
      nacht_hab_long[, Habitat := fcase(
        habitat_var=="pct_waldweg",  "Edge habitat",
        habitat_var=="pct_waldrand", "Edge habitat",
        habitat_var=="pct_wald",     "Forest interior",
        habitat_var=="pct_acker",    "Agricultural land"
      )]
      nacht_hab_long[, Phase := fifelse(tage_seit <= 5, "Early\n(day 1-5)", "Late\n(day 6+)")]

      p_overview <- ggplot(nacht_hab_long,
                           aes(x=Phase, y=pct, fill=Habitat)) +
        geom_boxplot(alpha=0.75, outlier.shape=21, outlier.size=1, width=0.6) +
        facet_wrap(~Habitat, scales="free_y", ncol=4) +
        scale_fill_manual(values=c("Forest track"="#8B4513","Forest edge"="#f7c948",
                                    "Forest interior"="#1a6b1a","Agricultural land"="#d4a843"),
                          guide="none") +
        labs(title="Nightly habitat use: early vs. late phase (animal-night level)",
             subtitle="Each point = one animal-night | N-nights level (not per-fix)",
             x=NULL, y="Proportion of fixes per night (%)") +
        theme_bw(base_size=11) +
        theme(plot.title=element_text(face="bold"),
              strip.text=element_text(face="bold"))

      ggsave(file.path(out_ordner, "habitat_nightly_overview.png"),
             p_overview, width=13, height=5, dpi=150)
      cat("  habitat_nightly_overview.png\n\n")

      # Rohdaten und Zusammenfassung speichern fuer Excel
      individ_nacht_tab   <- nacht_hab
      individ_zusammen_tab <- igel_zusammenfassung

    } else {
      individ_nacht_tab <- NULL
      individ_zusammen_tab <- NULL
    }
  }, error=function(e) {
    cat(sprintf("  [WARN] Individuelle Verlaeufe: %s\n", conditionMessage(e)))
    individ_nacht_tab <- NULL
    individ_zusammen_tab <- NULL
  })

  # ── D. Verbesserte Habitatkarte mit KDE-Heatmap ─────────────────
  tryCatch({
    if (requireNamespace("adehabitatHR", quietly=TRUE)) {
      library(adehabitatHR)
      # Populationslevel-KDE (alle Fixes)
      af_karte <- habitat_fixes[!is.na(x) & !is.na(y)]
      sp_pts <- sp::SpatialPoints(
        coords=as.matrix(af_karte[, .(x,y)]),
        proj4string=sp::CRS(SRS_string="EPSG:32632"))
      kde_karte <- kernelUD(sp_pts, h="href", grid=150)
      poly_95 <- tryCatch(st_as_sf(getverticeshr(kde_karte, 95)), error=function(e) NULL)
      poly_50 <- tryCatch(st_as_sf(getverticeshr(kde_karte, 50)), error=function(e) NULL)

      release_sf_k <- st_as_sf(data.frame(x=514743.7, y=5805363.8),
                                coords=c("x","y"), crs=32632)
      study_k <- st_buffer(release_sf_k, 650)
      hab_k <- suppressMessages(st_intersection(
        st_transform(habitat_polygone, 32632), study_k))

      hab_f2 <- c("Wald"="#1a6b1a","Acker"="#d4a843","Wiese"="#90c050",
                   "Gehoelz_Strauch"="#4a8f3f","Gewaesser"="#4a90c4",
                   "Siedlung"="#c0a080","Garten"="#b8d8a0")

      # KDE als Raster extrahieren fuer Heatmap
      ud_raster <- tryCatch({
        as.data.frame(raster::rasterToPoints(
          raster::raster(as(kde_karte, "SpatialPixelsDataFrame"))))
      }, error=function(e) NULL)

      # Englische Habitat-Namen fuer Karte
      hab_k_en <- hab_k
      hab_k_en$habitat_en <- hab_labels_en[hab_k_en$habitat]
      hab_k_en$habitat_en[is.na(hab_k_en$habitat_en)] <- hab_k_en$habitat[is.na(hab_k_en$habitat_en)]
      hab_f2_en <- setNames(hab_f2, hab_labels_en[names(hab_f2)])
      hab_f2_en <- hab_f2_en[!is.na(names(hab_f2_en))]

      p_karte2 <- ggplot() +
        # Habitathintergrund
        geom_sf(data=hab_k, aes(fill=habitat), color="white", linewidth=0.3, alpha=0.55) +
        scale_fill_manual(values=hab_f2,
                          labels=c("Acker"="Agricultural land","Wald"="Forest",
                                   "Wiese"="Grassland","Gehoelz_Strauch"="Scrub",
                                   "Gewaesser"="Water","Siedlung"="Settlement"),
                          name="Habitat") +
        # KDE Heatmap (Dichtegradient)
        { if (!is.null(ud_raster))
          new_scale_fill <- NULL  # Placeholder — KDE als geom_contour_filled
        } +
        { if (!is.null(poly_95) && !is.null(poly_50)) {
            list(
              geom_sf(data=st_transform(poly_95, 4326), fill="#ff6b6b",
                      alpha=0.15, color="#cc0000", linewidth=1.0, linetype="dashed"),
              geom_sf(data=st_transform(poly_50, 4326), fill="#cc0000",
                      alpha=0.30, color="#cc0000", linewidth=1.5)
            )
          }
        } +
        # Auswilderungsvoliere
        geom_sf(data=release_sf_k, shape=8, color="black", size=6, stroke=2.5) +
        geom_sf(data=release_sf_k, shape=8, color="white", size=5, stroke=1.5) +
        coord_sf(crs=4326) +
        labs(title="Space use — Sachsenhagen study site",
             subtitle="Red shading: 50% (dark) and 95% (dashed) population-level KDE home range\nStar = release enclosure | n = 26 hedgehogs",
             x=NULL, y=NULL) +
        theme_bw(base_size=11) +
        theme(plot.title=element_text(face="bold"),
              plot.subtitle=element_text(size=9, color="grey40"),
              legend.position="right")

      ggsave(file.path(out_ordner, "habitat_karte_v2.png"),
             p_karte2, width=11, height=9, dpi=150)
      cat("  habitat_karte_v2.png\n")
    }
  }, error=function(e) cat(sprintf("  [WARN] Karte v2: %s\n", conditionMessage(e))))
}

cat("\nAlle Plots erstellt\n\n")


# ══════════════════════════════════════════════════════════════
# ABSCHNITT 7b: aKDE-HOMERANGE × UMWELTFAKTOREN
# ══════════════════════════════════════════════════════════════
# Zweite Analyseebene: nicht nächtliche Aktivität, sondern die
# GESAMTE Homerange (aKDE, akkurateste Methode) pro Tier.
# Fragestellungen:
#   - Haben Männchen größere Homeranges als Weibchen?
#   - Beeinflusst die Reha-Dauer die Homerange-Größe?
#   - Unterscheiden sich Tiere nach Saison der Auswilderung?
#   - Korreliert die Homerange-Größe mit der mittleren Temperatur
#     während der Trackingperiode?
# ══════════════════════════════════════════════════════════════
cat("====================================================\n")
cat("7b. aKDE-Homerange x Umweltfaktoren\n")
cat("====================================================\n\n")

# ── aKDE-Ergebnisse aus Block4b laden ────────────────────────
rds_4b_pfad <- file.path(projekt_root, "output", "Block4b_Einzeltier",
                          "Block4b_Ergebnisse.rds")

akde_tab <- NULL
if (file.exists(rds_4b_pfad)) {
  b4b_data <- tryCatch(readRDS(rds_4b_pfad), error = function(e) NULL)
  if (!is.null(b4b_data) && !is.null(b4b_data$vergleich_dt)) {
    akde_tab <- as.data.table(b4b_data$vergleich_dt)
    # Spalten: igel, kde_95ha, akde_95ha, akde_50ha, modell, ess_area
    cat(sprintf("  aKDE-Daten: %d Tiere aus Block4b geladen\n", nrow(akde_tab)))
  }
}

if (is.null(akde_tab)) {
  cat("  Block4b-RDS nicht gefunden — 7b wird uebersprungen\n")
  cat("  Bitte zuerst Block4b_Einzeltier_aKDE.R ausfuehren\n\n")
} else {

  # ── Metadaten ergänzen (Geschlecht, Reha-Dauer, Diagnose) ──
  meta_datei <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
  meta_5b <- tryCatch({
    m <- as.data.table(readxl::read_excel(meta_datei))
    data.table(
      igel         = trimws(m$individual),
      sex          = m$sex,
      time_reha    = as.numeric(m$time_reha),
      diagnose     = m$diagnosis_main,
      release_date = as.Date(m$date_release),
      tagging_date = as.Date(m$tagging_date)
    )
  }, error = function(e) NULL)

  if (!is.null(meta_5b)) {
    akde_tab <- merge(akde_tab, meta_5b, by = "igel", all.x = TRUE)
    cat(sprintf("  Metadaten ergaenzt: %d mit Geschlecht, %d mit Reha-Dauer\n",
                akde_tab[!is.na(sex), .N], akde_tab[!is.na(time_reha), .N]))
  }

  # Saison der Auswilderung (Monat des release_date)
  if ("release_date" %in% names(akde_tab)) {
    akde_tab[, release_monat := month(release_date)]
    akde_tab[, release_saison := fcase(
      release_monat %in% 4:5,   "Fruehjahr (Apr-Mai)",
      release_monat %in% 6:9,   "Sommer (Jun-Sep)",
      release_monat %in% 10:11, "Herbst (Okt-Nov)",
      default = "Andere"
    )]
    akde_tab[, release_saison := factor(release_saison,
      levels = c("Fruehjahr (Apr-Mai)", "Sommer (Jun-Sep)", "Herbst (Okt-Nov)"))]
  }

  # Mittlere Temperatur während Trackingperiode (aus DWD-Daten)
  # BUG-FIX: merge(..., by = character(0), allow.cartesian = TRUE) macht in data.table
  # KEIN Kreuzprodukt — es sucht stattdessen gemeinsame Schlüsselspalten und wirft
  # "missing from y: [igel]" weil wetter_dt keine igel-Spalte hat.
  # Korrekte Lösung: lapply pro Tier, Wetterdaten direkt per Datumsfilter zuordnen.
  if (!is.null(wetter_dt) && !is.null(meta_5b) && "release_date" %in% names(akde_tab)) {
    temp_per_tier <- rbindlist(lapply(akde_tab$igel, function(ig) {
      rd <- akde_tab[igel == ig, release_date][1L]
      if (is.na(rd)) return(NULL)
      w  <- wetter_dt[datum >= rd & datum <= rd + 60L]
      if (nrow(w) == 0L) return(NULL)
      data.table(igel               = ig,
                 mean_temp_tracking = round(mean(w$temp_mittel, na.rm = TRUE), 1))
    }), fill = TRUE)
    # BUG-FIX: fehlende {} — cat() lief immer, auch wenn merge uebersprungen wurde.
    # Ohne geschweifte Klammern gehoerte nur akde_tab <- merge(...) zum if-Block.
    # Die cat()-Zeile lief dann auf akde_tab ohne mean_temp_tracking → Crash.
    if (nrow(temp_per_tier) > 0L) {
      akde_tab <- merge(akde_tab, temp_per_tier, by = "igel", all.x = TRUE)
      cat(sprintf("  Mittlere Tracking-Temperatur: %d Tiere mit DWD-Daten\n",
                  akde_tab[!is.na(mean_temp_tracking), .N]))
    } else {
      cat("  Mittlere Tracking-Temperatur: keine Uebereinstimmungen mit Wetterdaten\n")
    }
  }

  cat("\n")

  # Farbpaletten
  farbe_sex2    <- c("Male"="#4575b4", "Female"="#d6604d",
                      "m"="#4575b4", "f"="#d6604d",
                      "M"="#4575b4", "F"="#d6604d",
                      "maennlich"="#4575b4", "weiblich"="#d6604d")
  farbe_saison2 <- c("Fruehjahr (Apr-Mai)"="#74add1",
                      "Sommer (Jun-Sep)"="#f46d43",
                      "Herbst (Okt-Nov)"="#2d6a4f")

  # Statistik: Wilcoxon Geschlecht
  hat_sex_akde  <- "sex" %in% names(akde_tab) && akde_tab[!is.na(sex) & !is.na(akde_95ha), .N] >= 4
  hat_reha_akde <- "time_reha" %in% names(akde_tab) && akde_tab[!is.na(time_reha) & !is.na(akde_95ha), .N] >= 5
  hat_sais_akde <- "release_saison" %in% names(akde_tab) && akde_tab[!is.na(release_saison) & !is.na(akde_95ha), .N] >= 4

  if (hat_sex_akde && length(unique(na.omit(akde_tab$sex))) == 2) {
    wx_sex <- tryCatch(wilcox.test(akde_95ha ~ sex, data = akde_tab[!is.na(sex)]),
                        error = function(e) NULL)
    if (!is.null(wx_sex))
      cat(sprintf("  Wilcoxon Geschlecht x aKDE 95%%: W=%.0f, p=%.4f\n",
                  wx_sex$statistic, wx_sex$p.value))
  }

  if (hat_reha_akde) {
    cor_reha <- tryCatch(
      cor.test(akde_tab$time_reha, akde_tab$akde_95ha, method = "spearman"),
      error = function(e) NULL)
    if (!is.null(cor_reha))
      cat(sprintf("  Spearman Reha-Dauer x aKDE 95%%: rho=%.3f, p=%.4f\n",
                  cor_reha$estimate, cor_reha$p.value))
  }
  cat("\n")

  # ── Plots ─────────────────────────────────────────────────
  plot_akde <- list()

  # Plot A1: Homerange nach Geschlecht
  if (hat_sex_akde) {
    p_akde_sex <- ggplot(akde_tab[!is.na(sex) & !is.na(akde_95ha)],
                          aes(x = sex, y = akde_95ha, fill = sex)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.5) +
      geom_jitter(width = 0.12, size = 3, alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 4.5, color = "black") +
      scale_fill_manual(values = farbe_sex2, guide = "none") +
      labs(title = "aKDE-Homerange nach Geschlecht",
           subtitle = "Maennchen haben typischerweise groessere Aktionsraeume",
           x = NULL, y = "aKDE 95% Homerange (ha)") +
      theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
    plot_akde[["sex"]] <- p_akde_sex
  }

  # Plot A2: Homerange nach Auswilderungssaison
  if (hat_sais_akde) {
    p_akde_sais <- ggplot(akde_tab[!is.na(release_saison) & !is.na(akde_95ha)],
                           aes(x = release_saison, y = akde_95ha, fill = release_saison)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.6) +
      geom_jitter(width = 0.15, size = 3, alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 4.5, color = "black") +
      scale_fill_manual(values = farbe_saison2, guide = "none") +
      labs(title = "aKDE-Homerange nach Auswilderungssaison",
           subtitle = "Unterscheiden sich Tiere je nach Jahreszeit der Auswilderung?",
           x = NULL, y = "aKDE 95% Homerange (ha)") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 15, hjust = 1))
    plot_akde[["saison"]] <- p_akde_sais
  }

  # Plot A3: Homerange vs. Reha-Dauer (Streudiagramm pro Tier)
  if (hat_reha_akde) {
    p_akde_reha <- ggplot(akde_tab[!is.na(time_reha) & !is.na(akde_95ha)],
                           aes(x = time_reha, y = akde_95ha)) +
      geom_point(aes(color = if (hat_sex_akde && !all(is.na(akde_tab$sex))) sex else "alle"),
                 size = 3.5, alpha = 0.85) +
      geom_smooth(method = "lm", se = TRUE, color = "grey30",
                  fill = "grey80", alpha = 0.2, linewidth = 1) +

      { if (hat_sex_akde) scale_color_manual(values = farbe_sex2, name = "Geschlecht")
        else scale_color_manual(values = c("alle" = "#4575b4"), guide = "none") } +
      labs(title = "Reha-Dauer vs. aKDE-Homerange (pro Tier)",
           subtitle = "Beeinflusst laengere Gefangenschaft die Aktionsraumgroesse?",
           x = "Rehabilitationsdauer (Tage)", y = "aKDE 95% Homerange (ha)") +
      theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
    plot_akde[["reha"]] <- p_akde_reha
  }

  # Plot A4: Homerange vs. mittlere Tracking-Temperatur
  hat_temp_akde <- "mean_temp_tracking" %in% names(akde_tab) &&
                    akde_tab[!is.na(mean_temp_tracking) & !is.na(akde_95ha), .N] >= 4
  if (hat_temp_akde) {
    rho_t <- tryCatch(cor(akde_tab$mean_temp_tracking, akde_tab$akde_95ha,
                           use = "complete.obs", method = "spearman"), error = function(e) NA)
    p_akde_temp <- ggplot(akde_tab[!is.na(mean_temp_tracking) & !is.na(akde_95ha)],
                           aes(x = mean_temp_tracking, y = akde_95ha)) +
      geom_point(aes(color = if (hat_sex_akde && !all(is.na(akde_tab$sex))) sex else "alle"),
                 size = 3.5, alpha = 0.85) +
      geom_smooth(method = "lm", se = TRUE, color = "grey30",
                  fill = "grey80", alpha = 0.2, linewidth = 1) +
      { if (hat_sex_akde) scale_color_manual(values = farbe_sex2, name = "Geschlecht")
        else scale_color_manual(values = c("alle" = "#4575b4"), guide = "none") } +
      annotate("text", x = max(akde_tab$mean_temp_tracking, na.rm=TRUE),
               y = max(akde_tab$akde_95ha, na.rm=TRUE),
               label = sprintf("Spearman rho = %.2f", rho_t),
               size = 3.5, hjust = 1, color = "grey30") +
      labs(title = "Mittlere Temperatur vs. aKDE-Homerange",
           subtitle = "Mittlere DWD-Temperatur waehrend der ersten 60 Tracking-Tage",
           x = "Mittlere Temperatur (Grad C)", y = "aKDE 95% Homerange (ha)") +
      theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
    plot_akde[["temp"]] <- p_akde_temp
  }

  # Kombinations-Plot
  if (length(plot_akde) >= 2) {
    n_plots <- length(plot_akde)
    p_akde_kombi <- if (n_plots == 4) {
      (plot_akde[[1]] + plot_akde[[2]]) / (plot_akde[[3]] + plot_akde[[4]])
    } else if (n_plots == 3) {
      plot_akde[[1]] + plot_akde[[2]] + plot_akde[[3]]
    } else {
      plot_akde[[1]] + plot_akde[[2]]
    }
    p_akde_kombi <- p_akde_kombi +
      plot_annotation(
        title    = "aKDE-Homerange x Umwelt- und Tier-Faktoren",
        subtitle = sprintf("N = %d Igel | aKDE (OUF anisotropic, ctmm) | Raute = Mittelwert",
                           akde_tab[!is.na(akde_95ha), .N]),
        theme = theme(plot.title    = element_text(face = "bold", size = 14),
                      plot.subtitle = element_text(size = 10, color = "grey40"))
      )
    ggsave(file.path(out_ordner, "akde_umwelt_vergleich.png"),
           p_akde_kombi, width = 14, height = 10, dpi = 150)
    cat("  akde_umwelt_vergleich.png gespeichert\n")
  }

  # Übersichtstabelle aKDE x Meta ausgeben
  cat("\n  aKDE-Zusammenfassung:\n")
  cols_show <- intersect(c("igel","sex","time_reha","release_saison",
                            "akde_95ha","akde_50ha","modell","ess_area",
                            "mean_temp_tracking"),
                          names(akde_tab))
  print(akde_tab[order(igel), ..cols_show])
  cat("\n")
}

# ── 8. EXPORT ─────────────────────────────────────────────────
cat("====================================================\n")
cat("8. Export — Excel + Word\n")
cat("====================================================\n\n")

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- createWorkbook()
  st_h <- createStyle(fontName="Arial", fontSize=11, fontColour="white",
                       fgFill="#2C5F8A", halign="center",
                       textDecoration="bold", wrapText=TRUE)
  st_n <- createStyle(fontName="Arial", fontSize=10)
  st_2 <- createStyle(fontName="Arial", fontSize=10, numFmt="0.00", halign="center")

  # Hilfsfunktion: sicheres Schreiben (bereinigt Inf/NaN, gibt Fehlermeldung)
  safe_write <- function(wb, sheet, df, start_row=1, info=NULL, ...) {
    tryCatch({
      df <- as.data.frame(df)
      # Inf und NaN zu NA (verhindert #NUM? in Excel)
      df[] <- lapply(df, function(col) {
        if (is.numeric(col)) {
          col[is.infinite(col) | is.nan(col)] <- NA
          col <- round(col, 4)
        }
        col
      })
      if (!is.null(info)) {
        writeData(wb, sheet, data.frame(Info=info), startRow=start_row, colNames=FALSE)
        start_row <- start_row + 1
      }
      writeData(wb, sheet, df, startRow=start_row, colNames=TRUE)
      # Kopfzeile formatieren
      hs <- createStyle(fontName="Arial", fontSize=10, textDecoration="bold",
                        fgFill="#E8E8E8", border="Bottom")
      addStyle(wb, sheet, hs, rows=start_row, cols=seq_len(ncol(df)), gridExpand=TRUE)
      cat(sprintf("    ✓ %s: %d Zeilen, %d Spalten\n", sheet, nrow(df), ncol(df)))
    }, error=function(e) {
      cat(sprintf("    [FEHLER] %s: %s\n", sheet, conditionMessage(e)))
    })
  }
  st_0 <- createStyle(fontName="Arial", fontSize=10, numFmt="0", halign="center")
  st_hdr_row <- function(wb, sh, headers, row = 1) {
    writeData(wb, sh, as.data.frame(t(headers)), startRow=row, colNames=FALSE)
    addStyle(wb, sh, st_h, rows=row, cols=seq_along(headers), gridExpand=TRUE)
  }

  # ── Tagesweise Wetterdaten (volle Aufloesung) ────────────────
  if (!is.null(wetter_dt) && nrow(wetter_dt) > 0) {
    xl_w_voll <- wetter_dt[datum >= datum_range[1] & datum <= datum_range[2]]
    addWorksheet(wb, "Wetterdaten_Tage", tabColour="#d73027")
    setColWidths(wb, "Wetterdaten_Tage", cols=1:8,
                 widths=c(14,14,12,12,18,16,12,18))
    writeData(wb, "Wetterdaten_Tage",
              paste0("DWD: Station ", station_name, " | ", round(station_dist,1), " km"),
              startRow=1)
    # BUG-FIX: writeDataTable statt writeData — schreibt korrekt alle Spalten
    # ohne die as.data.frame(t())-Transpositions-Problematik des alten st_hdr_row
    if (nrow(xl_w_voll) > 0) {
      out_w <- xl_w_voll[, .(
        Datum          = datum,
        Temp_Mittel    = round(temp_mittel,       1),
        Temp_Min       = round(temp_min,           1),
        Temp_Max       = round(temp_max,           1),
        Niederschlag_mm= round(niederschlag_mm,    1),
        Sonnenschein_h = round(sonnenschein_h,     1),
        Regentag       = (niederschlag_mm > 2),
        Saison         = saison_temp
      )]
      writeDataTable(wb, "Wetterdaten_Tage", as.data.frame(out_w),
                     startRow=3, tableStyle="TableStyleLight9")
      cat(sprintf("  Wetterdaten_Tage: %d Tage geschrieben\n", nrow(out_w)))
    } else {
      writeData(wb, "Wetterdaten_Tage", "Keine Daten im Beobachtungszeitraum",
                startRow=3)
      cat("  [WARN] Wetterdaten: 0 Zeilen im Beobachtungszeitraum!\n")
      cat(sprintf("  DWD-Bereich: %s - %s | Igel-Bereich: %s - %s\n",
                  min(wetter_dt$datum), max(wetter_dt$datum),
                  datum_range[1], datum_range[2]))
    }
  }

  # ── Monatliche Wetterzusammenfassung ─────────────────────────
  if (!is.null(monat_wetter) && nrow(monat_wetter) > 0) {
    addWorksheet(wb, "Wetter_Monatlich", tabColour="#f46d43")
    out_mw <- monat_wetter[, .(
      Monat             = monat_str,
      Temp_Mittel_C     = temp_mittel_monat,
      Temp_Min_C        = temp_min_monat,
      Temp_Max_C        = temp_max_monat,
      Niederschlag_mm   = niederschlag_sum,
      Sonnenschein_h    = sonnenschein_sum,
      N_Regentage_gt2mm = n_regentage,
      N_Wetterdaten     = n_tage_wetter
    )]
    writeDataTable(wb, "Wetter_Monatlich", as.data.frame(out_mw),
                   startRow=1, tableStyle="TableStyleLight9")
    cat(sprintf("  Wetter_Monatlich: %d Monate\n", nrow(out_mw)))
  }

  # ── Monatliche Igel-Metriken ──────────────────────────────────
  if (!is.null(monat_igel) && nrow(monat_igel) > 0) {
    addWorksheet(wb, "Igel_Monatlich", tabColour="#2d6a4f")
    out_mi <- monat_igel[, .(
      Monat              = monat_str,
      Median_Radius_m    = median_radius_m,
      Median_Dist_Ausw_m = median_dist_m,
      N_Beobachtungsnaechte = n_naechte_igel,
      N_Igel             = n_igel
    )]
    writeDataTable(wb, "Igel_Monatlich", as.data.frame(out_mi),
                   startRow=1, tableStyle="TableStyleLight9")
    cat(sprintf("  Igel_Monatlich: %d Monate\n", nrow(out_mi)))
  }

  # ── Tagesweise Nacht x Wetter ────────────────────────────────
  if ("temp_mittel" %in% names(nacht_wetter) &&
      nacht_wetter[!is.na(temp_mittel), .N] > 0) {
    addWorksheet(wb, "Nacht_Wetter", tabColour="#9e9ac8")
    out_nw <- nacht_wetter[!is.na(temp_mittel),
      .(igel, datum, tage_seit,
        Radius_95pct_m = round(radius_95pct_m),
        Dist_Ausw_m    = round(dist_release_m),
        Temp_Mittel    = round(temp_mittel, 1),
        Niederschlag   = round(niederschlag_mm, 1),
        Regentag       = (!is.na(niederschlag_mm) & niederschlag_mm > 2),
        Platzhalter    = NA_character_)]  # Mondphase entfernt
    writeDataTable(wb, "Nacht_Wetter", as.data.frame(out_nw),
                   startRow=1, tableStyle="TableStyleLight9")
    cat(sprintf("  Nacht_Wetter: %d Naechte\n", nrow(out_nw)))
  }

  # Mondphase-Tab entfernt

  cat("  Schreibe Habitat-Sheets...\n")

  if (!is.null(habitat_fixes) && nrow(habitat_fixes) > 0) {

    # ── Sheet 1: Nutzung Population ──────────────────────────────
    addWorksheet(wb, "Nutzung_Population", tabColour="#1a6b1a")
    pop_use <- habitat_fixes[, .N, by=habitat_grob]
    pop_use[, Genutzt_pct := round(100 * N / sum(N), 1)]
    setnames(pop_use, c("habitat_grob","N"), c("Habitat","N_Fixes"))
    setorder(pop_use, -Genutzt_pct)
    safe_write(wb, "Nutzung_Population", pop_use,
               info="Habitatnutzung Population: Anteil Nacht-Fixes pro Habitat")

    # ── Sheet 2: Selektion ────────────────────────────────────────
    if (exists("sel_out") && !is.null(sel_out)) {
      addWorksheet(wb, "Selektion", tabColour="#d4a843")
      sel_clean <- sel_out
      # Chi-Quadrat-Zeile ENTFERNT (pseudorepliziert, siehe oben).
      # Selektionsratio/Manly's alpha sind deskriptiv; Inferenz via Wilcoxon-Sheet.
      safe_write(wb, "Selektion", sel_clean,
                 info="Selektionsratio = Genutzt% / Verfuegbar% | >1 = bevorzugt | <1 = gemieden | DESKRIPTIV — statistische Inferenz ausschliesslich ueber das Wilcoxon-Sheet (per-Tier, Einheit = Tier)")
    }

    # ── Sheet: Per-Tier-Wilcoxon Habitatselektion (Hauptanalyse) ─────────────
    if (exists("wilcox_habitat_dt") && !is.null(wilcox_habitat_dt) &&
        nrow(wilcox_habitat_dt) > 0) {
      addWorksheet(wb, "Selektion_Wilcoxon", tabColour="#2166ac")
      safe_write(wb, "Selektion_Wilcoxon", wilcox_habitat_dt,
                 info=paste0(
                   "HAUPTANALYSE: Per-Tier-Wilcoxon-Test (n=", nrow(wilcox_habitat_dt), " Habitate) | ",
                   "Testet ob Median der per-Tier-Proportionen von Verfuegbarkeit abweicht | ",
                   "Korrekte stat. Einheit = Tier (nicht Einzelfix) | ",
                   "W = Wilcoxon-Statistik | p-Wert: *** <0.001 | ** <0.01 | * <0.05 | n.s. = nicht signifikant"
                 ))
    }

    # ── Sheet 3: Nutzung pro Igel ─────────────────────────────────
    addWorksheet(wb, "Nutzung_pro_Igel", tabColour="#4a8f3f")
    hab_wide <- tryCatch(
      dcast(habitat_fixes[, .N, by=.(igel, habitat_grob)][
              , pct:=round(100*N/sum(N),1), by=igel],
            igel ~ habitat_grob, value.var="pct", fill=0),
      error=function(e) NULL)
    if (!is.null(hab_wide)) {
      n_igel_tot <- habitat_fixes[, .(N_Fixes_gesamt=.N), by=igel]
      hab_wide <- merge(n_igel_tot, hab_wide, by="igel")
      safe_write(wb, "Nutzung_pro_Igel", hab_wide,
                 info="Anteil Nacht-Fixes pro Habitat (%) fuer jedes Tier")
    }

    # ── Sheet 4: Wilcoxon Temporal ───────────────────────────────
    addWorksheet(wb, "Temporal_Wilcoxon", tabColour="#f7c948")
    if (exists("temporal_wilcox_tab") && !is.null(temporal_wilcox_tab) &&
        nrow(temporal_wilcox_tab) > 0) {
      safe_write(wb, "Temporal_Wilcoxon", temporal_wilcox_tab,
                 info="Wilcoxon Signed-Rank Test (gepaart): Fruehphase (Tag 1-5) vs. Spaetphase (Tag 6+) | nur Igel mit Daten in BEIDEN Phasen")
    }
    if (exists("per_tier_gepaart") && !is.null(per_tier_gepaart)) {
      per_wide <- tryCatch(
        dcast(per_tier_gepaart, igel + habitat ~ phase2, value.var="pct", fill=NA),
        error=function(e) NULL)
      if (!is.null(per_wide)) {
        row_off <- if (exists("temporal_wilcox_tab") && !is.null(temporal_wilcox_tab))
          nrow(temporal_wilcox_tab) + 4 else 1
        safe_write(wb, "Temporal_Wilcoxon", per_wide,
                   start_row=row_off,
                   info="Pro-Tier Habitatnutzung (%) nach Phase")
      }
    }
  }

  # ── Sheet 5: Individuen Zusammenfassung ──────────────────────
  if (exists("individ_zusammen_tab") && !is.null(individ_zusammen_tab) &&
      nrow(individ_zusammen_tab) > 0) {
    addWorksheet(wb, "Individuen_Uebersicht", tabColour="#8B4513")
    safe_write(wb, "Individuen_Uebersicht", individ_zusammen_tab,
               info="Edge = Waldweg + Waldrand | Trend: zunehmend/>3 PP, abnehmend/<-3 PP, sonst stabil | Edge_Frueh/Spaet = mittl. % Nacht-Fixes in Edge-Habitat")
  }

  # ── Sheet 6: Naechtliche Rohdaten pro Tier ───────────────────
  if (exists("individ_nacht_tab") && !is.null(individ_nacht_tab) &&
      nrow(individ_nacht_tab) > 0) {
    addWorksheet(wb, "Individuen_Naechte", tabColour="#d4a843")
    out_ind <- individ_nacht_tab[order(igel, tage_seit), .(
      Igel=igel, Datum=as.character(datum), Tage_seit=tage_seit,
      N_Fixes=n_fixes,
      Waldweg_pct=round(pct_waldweg,1), Waldrand_pct=round(pct_waldrand,1),
      Edge_pct=round(pct_edge,1), Wald_pct=round(pct_wald,1),
      Acker_pct=round(pct_acker,1)
    )]
    safe_write(wb, "Individuen_Naechte", out_ind,
               info="Naechtliche Habitatnutzung pro Tier (% der Fixes in dieser Nacht)")
  }

  # ── Sheet 7: Rohdaten Sample ─────────────────────────────────
  if (!is.null(habitat_fixes) && nrow(habitat_fixes) > 0) {
    addWorksheet(wb, "Rohdaten_Fixes", tabColour="#636363")
    out_fixes <- habitat_fixes[, .(
      Igel=igel, Datum=as.character(datum),
      Tage_seit=if ("tage_seit" %in% names(habitat_fixes)) tage_seit else NA_integer_,
      X_UTM32N=round(x,1), Y_UTM32N=round(y,1), Habitat=habitat_grob,
      Phase=if ("phase2" %in% names(habitat_fixes)) as.character(phase2) else NA_character_
    )]
    if (nrow(out_fixes) > 30000) out_fixes <- out_fixes[sample(.N, 30000)]
    safe_write(wb, "Rohdaten_Fixes", out_fixes,
               info=sprintf("Sample von %d Nacht-Fixes mit Habitatklassifikation (max. 30.000)", nrow(out_fixes)))
  }

  excel_pfad <- file.path(out_ordner, "Block5_Uebersicht.xlsx")
  saveWorkbook(wb, excel_pfad, overwrite=TRUE)
  cat(sprintf("Excel gespeichert: %s\n\n", basename(excel_pfad)))
}

saveRDS(list(
  wetter_dt        = wetter_dt,
  nacht_wetter     = nacht_wetter,
  fixes_wetter     = fixes_wetter,
  habitat_fixes    = habitat_fixes,
  habitat_polygone = habitat_polygone
), file.path(out_ordner, "Block5_Ergebnisse.rds"))
cat("RDS gespeichert\n\n")

# ══════════════════════════════════════════════════════════════
# 9. PUBLIKATIONSFIGUR + TABELLE — Habitat Selection
# ══════════════════════════════════════════════════════════════
cat("====================================================\n")
cat("9. Publikationsfigur — Habitat Selection\n")
cat("====================================================\n\n")

if (exists("sel_dt") && !is.null(sel_dt) &&
    requireNamespace("ggplot2",   quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE)) {
  tryCatch({

    # ── Farbpalette (blau / grau) ─────────────────────────────
    col_avail    <- "#B4B2A9"   # hellgrau — verfügbar
    col_selected <- "#000000"   # blau     — bevorzugt
    col_avoided  <- "#B4B2A9"   # dunkelgrau — gemieden

    # ── X-Achsenlabels (Vollnamen, 45° rotiert) ──────────────
    hab_wrap <- c(
      "Edge habitat"      = "Edge habitat",
      "Forest interior"   = "Forest interior",
      "Agricultural land" = "Agricultural land",
      "Grassland"         = "Grassland"
    )

    # ── Daten vorbereiten ─────────────────────────────────────
    main_hab <- c("Edge habitat", "Forest interior",
                  "Agricultural land", "Grassland")
    sel_pub <- sel_dt[habitat %in% main_hab, .(
      hab      = factor(habitat, levels = main_hab),
      avail100 = round(avail_pct * 100, 1),
      use100   = round(use_pct   * 100, 1),
      ratio    = round(selektionsratio, 2)
    )]
    # Character-Spalte für fill (robuster als logical)
    sel_pub[, sel_dir := ifelse(ratio > 1, "Selected", "Avoided")]
    sel_pub[, sel_dir := factor(sel_dir, levels = c("Selected","Avoided"))]
    # Mindest-Balkenhöhe für Grassland (ratio ≈ 0)
    sel_pub[, ratio_display := pmax(ratio, 0.04)]

    # Long format für gepaarte Balken
    sel_long <- melt(
      sel_pub[, .(hab, avail100, use100)],
      id.vars       = "hab",
      measure.vars  = c("avail100", "use100"),
      variable.name = "typ", value.name = "pct"
    )
    sel_long[, typ_label := factor(
      ifelse(typ == "avail100", "Available", "Used"),
      levels = c("Available", "Used")
    )]
    # Labels: kleine Werte gesondert behandeln
    sel_long[, bar_label := ifelse(pct < 0.05, "< 0.1%",
                                   paste0(pct, "%"))]
    # Grassland-Used: Label nach rechts verschieben (hjust=0) um Überlappung
    # zu vermeiden → wird als separater geom_text Layer ergänzt
    sel_long_main  <- sel_long[!(hab == "Grassland" & typ == "use100")]
    sel_long_grass <- sel_long[ (hab == "Grassland" & typ == "use100")]

    # ── Panel (a): Gepaarte Balken ────────────────────────────
    p_A <- ggplot(sel_long,
                  aes(x = hab, y = pct, fill = typ_label)) +
      geom_col(position = position_dodge(width = 0.7),
               width = 0.65, color = "white", linewidth = 0.25) +
      # Alle Labels ausser Grassland-Used
      geom_text(data = sel_long_main,
                aes(label = bar_label),
                position = position_dodge(width = 0.7),
                vjust = -0.4, size = 2.9, color = "grey25") +
      # Grassland-Used: nach rechts versetzt (hjust = 0)
      geom_text(data = sel_long_grass,
                aes(label = bar_label),
                position = position_dodge(width = 0.7),
                vjust = -0.4, hjust = -0.1, size = 2.9, color = "grey25") +
      scale_x_discrete(labels = hab_wrap) +
      scale_fill_manual(
        values = c("Available" = "#B4B2A9", "Used" = "#000000"),
        name   = NULL
      ) +
      scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        expand = expansion(mult = c(0, 0.2))
      ) +
      labs(x = NULL, y = "Proportion (%)", tag = "(a)") +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "top",
        legend.key.size = unit(0.4, "cm"),
        legend.text     = element_text(size = 10),
        axis.text.x     = element_text(size = 10, color = "grey20",
                                        angle = 45, hjust = 1, vjust = 1),
        axis.text.y     = element_text(size = 10),
        axis.title.y    = element_text(size = 11),
        plot.tag        = element_text(size = 13, face = "bold"),
        plot.margin     = margin(5, 5, 10, 5),
        plot.background = element_rect(colour = "black",
                                       fill = NA, linewidth = 0.6)
      )

    # ── Panel (b): Selektionsratio ────────────────────────────
    # "no preference" rechts oben über Grassland platzieren
    # (dort ist der Balken fast 0 → viel freier Raum)
    p_B <- ggplot(sel_pub,
                  aes(x = hab, y = ratio_display, fill = sel_dir)) +
      geom_hline(yintercept = 1, linetype = "dashed",
                 color = "grey55", linewidth = 0.9) +
      geom_col(width = 0.65, color = "white", linewidth = 0.25) +
      geom_text(aes(
        y     = ratio_display,
        label = ifelse(ratio < 0.01, "< 0.01", sprintf("%.2f", ratio))
      ), vjust = -0.4, size = 3.2, color = "grey55") +
      # "no preference" über dem Grassland-Balken (x=4), knapp über Linie
      annotate("text", x = 3.5, y = 1.07,
               label = "— no preference", hjust = 1, vjust = 0,
               size = 2.8, color = "grey55", fontface = "italic") +
      scale_x_discrete(labels = hab_wrap) +
      scale_fill_manual(
        values = c("Selected" = "#000000", "Avoided" = "#888780"),
        name   = NULL
      ) +
      scale_y_continuous(
        limits = c(0, 2.2),
        breaks = seq(0, 2, 0.5),
        expand = expansion(mult = c(0, 0.08))
      ) +
      labs(x = NULL, y = "Selection ratio (used / available)",
           tag = "(b)") +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "top",
        legend.key.size = unit(0.4, "cm"),
        legend.text     = element_text(size = 10),
        axis.text.x     = element_text(size = 10, color = "grey20",
                                        angle = 45, hjust = 1, vjust = 1),
        axis.text.y     = element_text(size = 10),
        axis.title.y    = element_text(size = 11),
        plot.tag        = element_text(size = 13, face = "bold"),
        plot.margin     = margin(5, 15, 10, 5),
        plot.background = element_rect(colour = "black",
                                       fill = NA, linewidth = 0.6)
      )

    # ── Patchwork kombinieren ─────────────────────────────────
    p_hab_pub <- p_A + p_B +
      plot_layout(ncol = 2, widths = c(1.1, 1))

    # ── Speichern (höher wegen 2-zeiliger x-Achse) ───────────
    fig_png <- file.path(out_ordner, "Fig_Habitat_Selection.png")
    fig_pdf <- file.path(out_ordner, "Fig_Habitat_Selection.pdf")
    ggsave(fig_png, p_hab_pub, width = 20, height = 12,
           dpi = 300, units = "cm")
    ggsave(fig_pdf, p_hab_pub, width = 20, height = 12,
           units = "cm")
    cat(sprintf("  ✓ Fig_Habitat_Selection.png gespeichert\n"))
    cat(sprintf("  ✓ Fig_Habitat_Selection.pdf gespeichert\n\n"))

    # ── Tabelle drucken ───────────────────────────────────────
    tab_sel <- sel_pub[order(match(hab, main_hab)), .(
      Habitat           = as.character(hab),
      `Available (%)`   = avail100,
      `Used (%)`        = ifelse(use100 < 0.1, "< 0.1",
                                 as.character(use100)),
      `Selection ratio` = sprintf("%.2f", ratio),
      Preference        = ifelse(selected, "Selected", "Avoided")
    )]
    cat("  Habitat selection table (for manuscript):\n")
    print(as.data.frame(tab_sel), row.names = FALSE)
    cat(paste0(
      "\n  Note: inference via per-animal Wilcoxon (n = 25 hedgehogs); ",
      "pooled fix-level chi-square removed (pseudoreplicated)\n",
      "  Edge temporal: W = 68, p = 0.705, n = 17 hedgehogs\n\n"
    ))

  }, error = function(e) {
    cat(sprintf("  [WARN] Publikationsfigur Habitat: %s\n\n",
                conditionMessage(e)))
  })
} else {
  cat("  sel_dt nicht verfuegbar — Abschnitt uebersprungen\n\n")
}

if (requireNamespace("officer", quietly=TRUE) &&
    requireNamespace("flextable", quietly=TRUE)) {

  doc <- read_docx()
  doc <- doc |>
    body_add_par("Block 5: Umweltdaten — Wetter & Habitat", style="heading 1") |>
    body_add_par(paste0("Erstellt: ", format(Sys.Date(), "%d.%m.%Y"))) |>
    body_add_par("") |>
    body_add_par("1. Wetterdaten (DWD)", style="heading 2") |>
    body_add_par(
      if (!is.null(wetter_dt))
        sprintf(paste0(
          "Taegliche Klimadaten wurden von der DWD-Station %s ",
          "(Entfernung: %.1f km) fuer den Beobachtungszeitraum %s bis %s geladen. ",
          "Die Tagesmitteltemperatur lag zwischen %.1f und %.1f Grad C."),
          station_name, station_dist,
          datum_range[1], datum_range[2],
          min(wetter_dt$temp_mittel, na.rm=TRUE),
          max(wetter_dt$temp_mittel, na.rm=TRUE))
      else "DWD-Wetterdaten konnten nicht geladen werden."
    ) |>
    body_add_par("") |>
    body_add_par("2. Habitatnutzung & Bewegung", style="heading 2") |>
    body_add_par(paste0(
      "Analysiert wurden Habitatnutzung (Edge habitat = Waldrand+Waldweg kombiniert, Forest interior, Agricultural land). ",
      "naechtliche Aktivitaetsradien und Distanzen zur Auswilderungsvoliere."
    )) |>
    body_add_par("")

  for (plot_name in c("temperaturverlauf.png", "wetter_homerange.png",
                       "habitat_karte.png", "habitat_nutzung.png",
                       "habitat_nutzung_vs_verfuegbar.png")) {
    pfad <- file.path(out_ordner, plot_name)
    if (file.exists(pfad)) {
      h_cm <- if (grepl("wetter_homerange", plot_name)) 13 else 10
      doc <- body_add_img(doc, src=pfad, width=17/2.54, height=h_cm/2.54)
      doc <- body_add_par(doc, "")
    }
  }

  word_pfad <- file.path(out_ordner, "Block5_Bericht.docx")
  print(doc, target=word_pfad)
  cat(sprintf("Word gespeichert: %s\n\n", basename(word_pfad)))
}

cat("====================================================\n")
cat("Block 5 abgeschlossen!\n")
cat("====================================================\n")
