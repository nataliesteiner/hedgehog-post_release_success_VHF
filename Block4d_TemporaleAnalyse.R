# ==============================================================
# Block 4d — temporal movement analysis
# ==============================================================
# Project:  Hedgehog VHF telemetry, Sachsenhagen, Lower Saxony
# Author:   Natalie Steiner
#
# --------------------------------------------------------------
# METHODS
# --------------------------------------------------------------
#
# Data:
#   VHF multilateration data (tRackIT system), sub-minute
#   resolution, several weeks post-release. Only night fixes
#   (Night == TRUE) are analysed (nocturnal species).
#
# 1. Nightly movement metrics
#
#    Path length (m):
#      Sum of all Euclidean distances between consecutive
#      fixes of one night = actual distance travelled

# ── Pakete ────────────────────────────────────────────────────
pakete <- c("adehabitatHR", "sp", "sf", "data.table", "ggplot2",
            "patchwork", "openxlsx", "readxl", "lubridate", "viridis",
            "officer", "flextable")
neu <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) install.packages(neu)

suppressPackageStartupMessages({
  library(adehabitatHR); library(sp);        library(sf)
  library(data.table);   library(ggplot2);   library(patchwork)
  library(openxlsx);     library(readxl);    library(lubridate)
  library(viridis);      library(officer);   library(flextable)
})
cat("Alle Pakete geladen\n\n")

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

# ── 1. Einstellungen ─────────────────────────────────────────
projekt_root  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
gpkg_ordner   <- file.path(projekt_root, "data", "kernel_files")
meta_datei    <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
output_ordner <- file.path(projekt_root, "output", "Block4d_Temporal")
dir.create(output_ordner, showWarnings = FALSE, recursive = TRUE)

FENSTER_N   <- 3      # Nächte pro gleitendem Hotspot-Fenster
HOTSPOT_PCT <- 10     # Obere X% der KDE = Hotspot
KDE_GRID    <- 100    # KDE-Rastergröße
KDE_H       <- "href" # KDE-Bandbreite

release_x <- 514743.7
release_y <- 5805363.8

cat("Output:", output_ordner, "\n\n")

# ── 2. Metadaten ─────────────────────────────────────────────
# Direkte Spaltennamen aus data_igel.xlsx (verifiziert):
# individual, sex, date_release, time_reha, diagnosis_main, tagging_date
meta_raw <- as.data.table(read_excel(meta_datei))

meta <- data.table(
  igel         = trimws(meta_raw$individual),       # trailing spaces entfernen
  sex          = meta_raw$sex,                        # "Male" / "Female"
  release_date = as.Date(meta_raw$date_release),     # readxl konvertiert Datum automatisch
  time_reha    = as.numeric(meta_raw$time_reha),     # Tage in Rehabilitation
  diagnose     = meta_raw$diagnosis_main,             # "orphan", "parasites", etc.
  tagging_date = as.Date(meta_raw$tagging_date)      # Besenderungsdatum
)
meta <- meta[!is.na(igel) & igel != "" & !is.na(release_date)]

# ── 3. GPKG laden ────────────────────────────────────────────
gpkg_dateien <- list.files(gpkg_ordner, pattern = "\\.gpkg$",
                            full.names = TRUE, recursive = TRUE)
gpkg_dateien <- gpkg_dateien[!grepl("\\(\\d+\\)\\.gpkg$", gpkg_dateien)]

# GPKG-Spalten (verifiziert durch direkte Dateianalyse):
#   Individual Name = "Igel1"  (sauberer Name, kein Filename-Parsing noetig)
#   Date = "2024-09-20"        (Datum als Text)
#   Time = "09:30:05"          (Uhrzeit als Text)
#   Night = numeric (> 0 = Nacht, NA/0 = Tag)
#   _time existiert in SQLite, wird aber von sf nicht korrekt importiert
#          → Date + Time Kombination verwenden

lade_fixes <- function(pfad) {
  # Originale Logik (funktioniert): st_read → as.data.table → nacht filtern →
  # st_as_sf (mit gespeichertem CRS) → st_coordinates.
  # Verbesserungen: einmaliges st_read (nicht mehr doppelt), sichtbare Fehlermeldungen.
  tryCatch({
    layer  <- st_layers(pfad)$name[1]
    sf_all <- st_read(pfad, layer = layer, quiet = TRUE)  # einmaliger Lesevorgang
    crs_pfad <- st_crs(sf_all)                             # CRS merken (kein zweites st_read)
    dat    <- as.data.table(sf_all)                        # sf → data.table (Spaltennamen bleiben)

    if (!"Night" %in% names(dat)) {
      message(sprintf("  [SKIP] %s: Keine 'Night'-Spalte", basename(pfad)))
      return(NULL)
    }

    # Night ist numerisch: > 0 = Nacht, NA oder 0 = Tag
    nacht <- dat[!is.na(Night) & Night > 0]
    if (nrow(nacht) < 15) {
      message(sprintf("  [SKIP] %s: nur %d Nacht-Fixes (<15)", basename(pfad), nrow(nacht)))
      return(NULL)
    }

    # Koordinaten extrahieren — Spaltenname dynamisch erkennen
    # (tRackIT GPKG hat "geometry", nicht "geom"; R-Konvention variiert)
    geom_col <- attr(sf_all, "sf_column")  # z.B. "geometry"
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
    # Fallback: Koordinaten direkt aus Geometrie-Objekten (wie Block5)
    if (!xy_ok && geom_col %in% names(nacht)) {
      geom_list <- nacht[[geom_col]]
      nacht[, x := vapply(geom_list, function(g) {
        v <- tryCatch(as.numeric(g), error=function(e) c(NA_real_,NA_real_))
        if (length(v)>=1) v[1] else NA_real_
      }, numeric(1))]
      nacht[, y := vapply(geom_list, function(g) {
        v <- tryCatch(as.numeric(g), error=function(e) c(NA_real_,NA_real_))
        if (length(v)>=2) v[2] else NA_real_
      }, numeric(1))]
    }

    # Robuster Spaltenname-Lookup: sf/data.table kann Leerzeichen in
    # Spaltennamen zu Punkten konvertieren ("Individual Name" → "Individual.Name").
    # grep() findet die Spalte unabhaengig davon welche Variante vorliegt.
    igel_col  <- grep("^Individual[ .]Name$",
                       names(nacht), value = TRUE)[1]
    date_col  <- grep("^Date$",  names(nacht), value = TRUE)[1]
    time_col  <- grep("^Time$",  names(nacht), value = TRUE)[1]

    if (is.na(igel_col)) {
      # Fallback: Tiernamen direkt aus GPKG-Dateinamen extrahieren
      igel_col_val <- sub("^(Igel\\d+)_.*$", "\\1", basename(pfad))
      message(sprintf("  [WARN] %s: 'Individual Name'-Spalte nicht gefunden, nutze '%s' aus Dateinamen",
                      basename(pfad), igel_col_val))
      nacht[, igel := igel_col_val]
    } else {
      nacht[, igel := trimws(nacht[[igel_col]])]
    }

    nacht[, datum := as.Date(if (!is.na(date_col)) nacht[[date_col]] else NA_character_)]
    nacht[, datetime := with_tz(
      as.POSIXct(paste(if (!is.na(date_col)) nacht[[date_col]] else "1970-01-01",
                       if (!is.na(time_col)) nacht[[time_col]] else "00:00:00"),
                 format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      "Europe/Berlin"
    )]

    nacht <- nacht[!is.na(x) & !is.na(y) & x != 0 & y != 0]
    if (nrow(nacht) < 15) {
      message(sprintf("  [SKIP] %s: zu wenige valide Koordinaten", basename(pfad)))
      return(NULL)
    }
    nacht[, .(igel, datum, datetime, x, y, Night)]

  }, error = function(e) {
    message(sprintf("  [FEHLER] %s: %s", basename(pfad), conditionMessage(e)))
    NULL
  })
}

cat(sprintf("Lade GPKG-Fixes aus %d Dateien...\n", length(gpkg_dateien)))
alle_fixes_liste <- lapply(gpkg_dateien, lade_fixes)
n_geladen <- sum(!sapply(alle_fixes_liste, is.null))
cat(sprintf("  %d/%d GPKG-Dateien erfolgreich geladen\n", n_geladen, length(gpkg_dateien)))

alle_fixes <- rbindlist(Filter(Negate(is.null), alle_fixes_liste), fill = TRUE)

# BUG-FIX: Hard stop wenn keine Fixes geladen — verhindert leere Outputs
# ohne jegliche Fehlermeldung (passierte bei stillem st_read-Fehler)
if (nrow(alle_fixes) == 0) {
  stop(paste0(
    "\nKeine GPKG-Fixes geladen! Moegliche Ursachen:\n",
    "  (a) sf-Paket-Problem: install.packages('sf') neu ausfuehren\n",
    "  (b) GPKG-Pfad falsch: ", gpkg_ordner, "\n",
    "  (c) Alle Tiere haben < 15 Nacht-Fixes (Night > 0)\n",
    "  (d) Fehlermeldungen oben pruefen!\n"
  ))
}

alle_fixes <- merge(alle_fixes, meta[, .(igel, release_date, sex, time_reha, diagnose, tagging_date)],
                    by = "igel", all.x = TRUE)
alle_fixes[, tage_seit := as.numeric(datum - release_date)]

# Diagnostik: Merge-Ergebnis pruefen
n_kein_release <- alle_fixes[is.na(release_date), uniqueN(igel)]
if (n_kein_release > 0) {
  warning(sprintf(
    "%d Igel ohne release_date nach Merge (igel-Namen-Mismatch?):\n  %s",
    n_kein_release,
    paste(alle_fixes[is.na(release_date), unique(igel)], collapse=", ")
  ))
}

alle_fixes <- alle_fixes[!is.na(tage_seit) & tage_seit >= 0]
setorder(alle_fixes, igel, datetime)

if (nrow(alle_fixes) == 0) {
  stop(paste0(
    "\nalle_fixes nach tage_seit-Filter leer!\n",
    "  Wahrscheinlich: release_date NA fuer alle Tiere nach Merge.\n",
    "  Igel-Namen in GPKG vs. Excel pruefen (Gross-/Kleinschreibung, Leerzeichen).\n"
  ))
}

n_mit_zeit <- alle_fixes[!is.na(datetime), .N]
cat(sprintf("  %d Fixes, %d Igel | Timestamps: %d/%d gueltig (%.0f%%)\n\n",
            nrow(alle_fixes), uniqueN(alle_fixes$igel),
            n_mit_zeit, nrow(alle_fixes),
            100 * n_mit_zeit / nrow(alle_fixes)))

# ══════════════════════════════════════════════════════════════
# HILFSFUNKTIONEN
# ══════════════════════════════════════════════════════════════

# Tortuosität: Straightness Index (Netto / Weglänge)
tortuositaet <- function(x, y) {
  if (length(x) < 2) return(NA_real_)
  steps <- sqrt(diff(x)^2 + diff(y)^2)
  total <- sum(steps, na.rm = TRUE)
  if (total == 0) return(NA_real_)
  netto <- sqrt((tail(x,1) - x[1])^2 + (tail(y,1) - y[1])^2)
  round(netto / total, 4)
}

# Mittlerer absoluter Drehwinkel (°)
mean_turn <- function(x, y) {
  if (length(x) < 3) return(NA_real_)
  angles <- atan2(diff(y), diff(x))
  turns  <- abs(diff(angles))
  turns  <- ifelse(turns > pi, 2*pi - turns, turns)
  round(mean(turns, na.rm = TRUE) * 180 / pi, 2)
}

# KDE-Raster für Hotspot-Analyse
# BUG-FIX: raster::raster() und raster::as.data.frame() entfernt —
# das raster-Paket war nicht in der pakete-Liste und ist veraltet (seit terra).
# Ersetzt durch direkte sp-Konvertierung via SpatialPixelsDataFrame@data + coordinates(),
# kein zusätzliches Paket nötig.
kde_raster_dt <- function(x, y, grid = KDE_GRID, h = KDE_H) {
  if (length(x) < 10) return(NULL)
  tryCatch({
    sp_obj <- SpatialPoints(cbind(x, y),
                             proj4string = CRS("+proj=utm +zone=32 +datum=WGS84"))
    ud     <- kernelUD(sp_obj, h = h, grid = grid)
    spdf   <- as(ud, "SpatialPixelsDataFrame")
    coords <- sp::coordinates(spdf)
    dt     <- data.table(x = coords[, 1], y = coords[, 2],
                          dichte = spdf@data[[1]])
    dt[!is.na(dichte)]
  }, error = function(e) NULL)
}

hotspot_zentroid <- function(kde_dt, pct = HOTSPOT_PCT) {
  if (is.null(kde_dt) || nrow(kde_dt) == 0) return(c(NA_real_, NA_real_))
  hot <- kde_dt[dichte >= quantile(dichte, 1 - pct/100, na.rm = TRUE)]
  if (nrow(hot) == 0) return(c(NA_real_, NA_real_))
  c(weighted.mean(hot$x, hot$dichte), weighted.mean(hot$y, hot$dichte))
}

raster_overlap <- function(dt1, dt2) {
  # Hilfsfunktion: konvertiert beliebiges Objekt sicher zu data.table
  # Behandelt: data.table, data.frame, list(data.table) (list-column-Extraktion
  # in manchen data.table-Versionen gibt list-verpackte Objekte zurueck)
  to_dt <- function(x) {
    if (is.null(x)) return(NULL)
    # list-verpackt? (Listenspalten-Extraktion gibt manchmal list(dt) statt dt)
    if (!is.data.table(x) && is.list(x) && !is.data.frame(x)) {
      if (length(x) == 1L) x <- x[[1L]] else return(NULL)
    }
    if (is.data.table(x))  return(x)
    if (is.data.frame(x))  return(as.data.table(x))
    return(NULL)
  }

  dt1 <- to_dt(dt1); if (is.null(dt1)) return(NA_real_)
  dt2 <- to_dt(dt2); if (is.null(dt2)) return(NA_real_)

  needed <- c("x","y","dichte")
  if (!all(needed %in% names(dt1)) || !all(needed %in% names(dt2))) return(NA_real_)

  r    <- 50L
  dt1c <- copy(dt1); dt1c[, `:=`(xr = round(x/r)*r, yr = round(y/r)*r)]
  dt2c <- copy(dt2); dt2c[, `:=`(xr = round(x/r)*r, yr = round(y/r)*r)]
  m    <- merge(dt1c[, .(xr, yr, d1 = dichte)],
                dt2c[, .(xr, yr, d2 = dichte)], by = c("xr","yr"))
  if (nrow(m) == 0L) return(0)
  s1 <- sum(m$d1, na.rm=TRUE); s2 <- sum(m$d2, na.rm=TRUE)
  if (s1 == 0 || s2 == 0) return(0)
  m[, d1n := d1/s1][, d2n := d2/s2]
  round(sum(sqrt(m$d1n * m$d2n), na.rm = TRUE), 4)
}

# ══════════════════════════════════════════════════════════════
# HAUPTSCHLEIFE
# ══════════════════════════════════════════════════════════════
cat("====================================\n")
cat("Starte Analyse pro Igel...\n\n")

nacht_metriken_alle  <- list()
hotspot_verlauf_alle <- list()
igel_liste <- sort(unique(alle_fixes$igel))

for (igel_i in igel_liste) {
  cat(sprintf("-- %s --\n", igel_i))
  fixes_i   <- alle_fixes[igel == igel_i]
  naechte   <- sort(unique(fixes_i$datum))
  n_naechte <- length(naechte)

  if (n_naechte < 3) {
    cat(sprintf("  ! Nur %d Naechte -- ueberspringe\n\n", n_naechte))
    next
  }

  # ── A. Nächtliche Bewegungsmetriken ───────────────────────
  metriken <- lapply(naechte, function(nd) {
    f <- fixes_i[datum == nd]
    setorder(f, datetime)
    n <- nrow(f)
    if (n < 3) return(NULL)

    steps <- sqrt(diff(f$x)^2 + diff(f$y)^2)

    # BUG-FIX: Wenn datetime = NA (GPKG-Timestamps nicht parsbar),
    # wird der Zeitfilter übersprungen und ein Distanzfilter verwendet.
    # Ohne diesen Fix war weg_gesamt_m = 0 für alle Tiere.
    dt_sec <- tryCatch(
      as.numeric(diff(f$datetime), units = "secs"),
      error = function(e) rep(NA_real_, length(steps))
    )
    hat_zeit <- !all(is.na(dt_sec)) && any(dt_sec > 0, na.rm = TRUE)

    if (hat_zeit) {
      # Schritte mit Zeitlücke > 30 min ausschließen
      valid <- !is.na(steps) & !is.na(dt_sec) & dt_sec > 0 & dt_sec < 1800
    } else {
      # Kein Zeitfilter verfügbar: Distanzfilter (> 500 m = nicht plausibel)
      valid <- !is.na(steps) & steps < 500
    }
    steps_v <- steps[valid]

    weg_m      <- sum(steps_v, na.rm = TRUE)
    netto_m    <- sqrt((tail(f$x,1) - f$x[1])^2 + (tail(f$y,1) - f$y[1])^2)
    max_dist   <- max(sqrt((f$x - f$x[1])^2 + (f$y - f$y[1])^2), na.rm = TRUE)
    dist_rel   <- sqrt((median(f$x) - release_x)^2 + (median(f$y) - release_y)^2)
    tort       <- tortuositaet(f$x, f$y)
    mt         <- mean_turn(f$x, f$y)
    m_step     <- if (length(steps_v) > 0) mean(steps_v, na.rm=TRUE) else NA_real_
    sd_step    <- if (length(steps_v) > 1) sd(steps_v,   na.rm=TRUE) else NA_real_
    zx <- mean(f$x); zy <- mean(f$y)
    r95 <- quantile(sqrt((f$x - zx)^2 + (f$y - zy)^2), 0.95, na.rm=TRUE)

    data.table(
      igel             = igel_i,
      datum            = nd,
      tage_seit        = unique(f$tage_seit)[1],
      n_fixes          = n,
      weg_gesamt_m     = round(weg_m, 1),
      netto_m          = round(netto_m, 1),
      max_dist_start_m = round(max_dist, 1),
      dist_release_m   = round(dist_rel, 1),
      tortuositaet     = tort,
      mean_drehwinkel  = mt,
      mean_step_m      = round(m_step, 1),
      sd_step_m        = round(sd_step, 1),
      radius_95_m      = round(as.numeric(r95), 1)
    )
  })

  metriken_dt <- rbindlist(Filter(Negate(is.null), metriken), fill = TRUE)
  nacht_metriken_alle[[igel_i]] <- metriken_dt
  cat(sprintf("  %d Naechte | Weg: %.0f m (Median) | Distanz Ausw.: %.0f m (Median)\n",
              nrow(metriken_dt),
              median(metriken_dt$weg_gesamt_m, na.rm=TRUE),
              median(metriken_dt$dist_release_m, na.rm=TRUE)))

  # ── B. Temporale Hotspot-Analyse ──────────────────────────
  if (n_naechte >= FENSTER_N) {
    hotspot_verlauf <- lapply(seq_len(n_naechte - FENSTER_N + 1), function(i) {
      fn <- naechte[i:(i + FENSTER_N - 1)]
      f_fen <- fixes_i[datum %in% fn]
      if (nrow(f_fen) < 15) return(NULL)
      kdt  <- kde_raster_dt(f_fen$x, f_fen$y)
      if (is.null(kdt)) return(NULL)
      zent <- hotspot_zentroid(kdt)
      d_rel <- sqrt((zent[1]-release_x)^2 + (zent[2]-release_y)^2)
      n_hot <- kdt[dichte >= quantile(dichte, 1-HOTSPOT_PCT/100, na.rm=TRUE), .N]
      zm    <- if (nrow(kdt)>1) min(diff(sort(unique(kdt$x))), na.rm=TRUE) else 50

      data.table(
        igel          = igel_i,
        fenster_start = fn[1],
        fenster_mitte = fn[ceiling(FENSTER_N/2)],
        fenster_ende  = tail(fn,1),
        tage_seit     = fixes_i[datum == fn[ceiling(FENSTER_N/2)], first(tage_seit)],
        n_fixes_fen   = nrow(f_fen),
        hotspot_x     = round(zent[1], 1),
        hotspot_y     = round(zent[2], 1),
        hotspot_dist_release_m = round(d_rel, 1),
        hotspot_ha    = round(n_hot * zm^2 / 10000, 3),
        kde_dt_ref    = list(kdt)
      )
    })

    hs_dt <- rbindlist(Filter(Negate(is.null), hotspot_verlauf), fill=TRUE)
    if (nrow(hs_dt) >= 2) {
      hs_dt[, stabilitaet := NA_real_]
      n_ok <- 0L
      for (k in 2:nrow(hs_dt)) {
        ov <- tryCatch(
          raster_overlap(hs_dt$kde_dt_ref[[k-1L]], hs_dt$kde_dt_ref[[k]]),
          error = function(e) {
            # Erste Fehlermeldung pro Tier anzeigen (Diagnose)
            if (k == 2L)
              message(sprintf("  [WARN] Stabilitaet Fehler %s (k=2): %s | class=%s",
                              igel_i, conditionMessage(e),
                              paste(class(hs_dt$kde_dt_ref[[1L]]), collapse=",")))
            NA_real_
          })
        hs_dt[k, stabilitaet := ov]
        if (!is.na(ov)) n_ok <- n_ok + 1L
      }
      if (n_ok > 0L)
        cat(sprintf("  Stabilitaet: %d/%d Werte berechnet\n", n_ok, nrow(hs_dt)-1L))
      else if (nrow(hs_dt) >= 2)
        cat(sprintf("  [WARN] Stabilitaet: 0/%d Werte (Fehler oben pruefen)\n", nrow(hs_dt)-1L))
      hs_dt[, hotspot_shift_m := c(NA_real_,
               sqrt(diff(hotspot_x)^2 + diff(hotspot_y)^2))]
    }
    hs_dt[, kde_dt_ref := NULL]
    hotspot_verlauf_alle[[igel_i]] <- hs_dt
    cat(sprintf("  %d Hotspot-Fenster\n", nrow(hs_dt)))
  }
  cat("\n")
}

# ── Zusammenführen ────────────────────────────────────────────
nacht_tab   <- rbindlist(nacht_metriken_alle,  fill=TRUE)
hotspot_tab <- rbindlist(hotspot_verlauf_alle, fill=TRUE)
cat("====================================\n")
cat(sprintf("%d Nacht-Eintraege | %d Hotspot-Fenster\n\n",
            nrow(nacht_tab), nrow(hotspot_tab)))

# ══════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════
cat("Erstelle Plots...\n\n")

theme_bew <- theme_minimal(base_size=10) +
  theme(plot.title=element_text(face="bold",size=10),
        plot.subtitle=element_text(size=8,color="grey50"),
        panel.grid.minor=element_blank())

for (igel_i in igel_liste) {
  nm_i <- nacht_tab[igel == igel_i]
  hs_i <- if (igel_i %in% names(hotspot_verlauf_alle)) hotspot_verlauf_alle[[igel_i]] else NULL
  if (nrow(nm_i) == 0) next
  hl <- nrow(nm_i) >= 5

  # BUG-FIX: geom_smooth in ggplot2 >= 3.4.0 erfordert explizites formula-Argument.
  # Ohne formula = y ~ x wirft ggplot2 intern einen object$timestamp-Fehler.
  # Zusaetzlich: ggsave in tryCatch → ein fehlerhafter Plot stoppt nicht den ganzen Loop.
  sm <- function(...) geom_smooth(formula = y ~ x, ...)

  p1 <- ggplot(nm_i, aes(x=tage_seit, y=weg_gesamt_m)) +
    geom_col(fill="#3A86FF", alpha=0.7) +
    { if(hl) sm(method="loess",se=TRUE,color="#E8A020",linewidth=1.2,fill="#E8A020",alpha=0.2) } +
    labs(title="Taegliche Weglaenge", subtitle="Gesamtstrecke pro Nacht",
         x="Tage nach Auswilderung", y="Weglaenge (m)") + theme_bew

  p2 <- ggplot(nm_i, aes(x=tage_seit, y=dist_release_m)) +
    geom_line(color="#E63946",linewidth=1) + geom_point(color="#E63946",size=2,alpha=0.8) +
    { if(hl) sm(method="loess",se=TRUE,color="#457B9D",linewidth=1,fill="#457B9D",alpha=0.2) } +
    labs(title="Distanz zum Auswilderungspunkt", subtitle="Median-Position der Nacht",
         x="Tage nach Auswilderung", y="Distanz (m)") + theme_bew

  p3 <- ggplot(nm_i[!is.na(tortuositaet)], aes(x=tage_seit, y=tortuositaet)) +
    geom_line(color="#2D6A4F",linewidth=1) + geom_point(color="#2D6A4F",size=2,alpha=0.8) +
    { if(hl && nm_i[!is.na(tortuositaet), .N] >= 5)
        sm(method="loess",se=TRUE,color="#52B788",linewidth=1,fill="#52B788",alpha=0.2) } +
    geom_hline(yintercept=0.5,linetype="dashed",color="grey60") +
    scale_y_continuous(limits=c(0,1)) +
    labs(title="Tortuositaet (Straightness Index)",
         subtitle="0 = kurvenreich | 1 = geradlinig | Igel typisch < 0.10",
         x="Tage nach Auswilderung", y="Tortuositaet") + theme_bew

  p4 <- ggplot(nm_i[!is.na(mean_drehwinkel)], aes(x=tage_seit, y=mean_drehwinkel)) +
    geom_line(color="#9B2226",linewidth=1) + geom_point(color="#9B2226",size=2,alpha=0.8) +
    { if(hl && nm_i[!is.na(mean_drehwinkel), .N] >= 5)
        sm(method="loess",se=TRUE,color="#AE2012",linewidth=1,fill="#AE2012",alpha=0.2) } +
    labs(title="Mittlerer Drehwinkel",
         subtitle="180 deg = geradlinig | ~90 deg = Suchbewegung",
         x="Tage nach Auswilderung", y="Drehwinkel (deg)") + theme_bew

  p5 <- ggplot(nm_i[!is.na(mean_step_m)], aes(x=tage_seit, y=mean_step_m)) +
    geom_ribbon(aes(ymin=pmax(0,mean_step_m-sd_step_m), ymax=mean_step_m+sd_step_m),
                fill="#4361EE", alpha=0.15) +
    geom_line(color="#4361EE",linewidth=1) +
    labs(title="Mittlere Schrittlaenge", subtitle="+- SD (blaues Band)",
         x="Tage nach Auswilderung", y="Schrittlaenge (m)") + theme_bew

  p6 <- ggplot(nm_i[!is.na(radius_95_m)], aes(x=tage_seit, y=radius_95_m)) +
    geom_col(fill="#7B2D8B", alpha=0.65) +
    { if(hl) sm(method="loess",se=FALSE,color="#C77DFF",linewidth=1.2) } +
    labs(title="Aktivitaetsradius (95. Pz.)",
         subtitle="Radius um naechtlichen Schwerpunkt",
         x="Tage nach Auswilderung", y="Radius (m)") + theme_bew

  if (!is.null(hs_i) && nrow(hs_i) >= 2) {
    p7 <- ggplot(hs_i[!is.na(stabilitaet)], aes(x=tage_seit, y=stabilitaet)) +
      geom_line(color="#F4A261",linewidth=1.2) + geom_point(color="#E76F51",size=2.5) +
      { if(hs_i[!is.na(stabilitaet), .N] >= 5)
          sm(method="loess",se=TRUE,color="#264653",linewidth=1,fill="#264653",alpha=0.15) } +
      scale_y_continuous(limits=c(0,1)) +
      labs(title="Hotspot-Stabilitaet",
           subtitle=paste0("KDE-Ueberlappung ",FENSTER_N,"-Naechte-Fenster"),
           x="Tage nach Auswilderung", y="Stabilitaetsindex") + theme_bew

    p8 <- ggplot(hs_i[!is.na(hotspot_shift_m)], aes(x=tage_seit, y=hotspot_shift_m)) +
      geom_col(fill="#F4A261", alpha=0.75) +
      labs(title="Hotspot-Verschiebung (m)",
           x="Tage nach Auswilderung", y="Verschiebung (m)") + theme_bew

    p9 <- ggplot(hs_i[!is.na(hotspot_x)]) +
      geom_path(aes(x=hotspot_x,y=hotspot_y,color=tage_seit), linewidth=1,
                arrow=arrow(length=unit(0.15,"cm"),type="closed")) +
      geom_point(aes(x=hotspot_x,y=hotspot_y,color=tage_seit), size=2.5) +
      # annotate("point") durch geom_point() ersetzt — stabiler in allen ggplot2-Versionen
      geom_point(data=data.frame(x=release_x, y=release_y),
                 aes(x=x, y=y), shape=23, size=4, fill="red", color="white",
                 inherit.aes=FALSE) +
      scale_color_viridis_c(option="plasma",name="Tage") +
      labs(title="Hotspot-Pfad", subtitle="Bewegung des Nutzungsschwerpunkts",
           x="Easting", y="Northing") + theme_bew + coord_equal()

    p_g <- (p1|p2|p3)/(p4|p5|p6)/(p7|p8|p9) +
      plot_annotation(title=paste("Block 4d:",igel_i),
                      subtitle=sprintf("%d Naechte | %d Hotspot-Fenster",nrow(nm_i),nrow(hs_i)),
                      theme=theme(plot.title=element_text(face="bold",size=13),
                                  plot.subtitle=element_text(size=9,color="grey40")))
    tryCatch(
      ggsave(file.path(output_ordner,paste0(igel_i,"_temporal.png")),
             p_g, width=18, height=14, dpi=150),
      error = function(e)
        cat(sprintf("  ! Plot-Fehler %s (9-Panel): %s\n", igel_i, conditionMessage(e)))
    )
  } else {
    p_g <- (p1|p2)/(p3|p4)/(p5|p6) +
      plot_annotation(title=paste("Block 4d:",igel_i),
                      subtitle=sprintf("%d Naechte",nrow(nm_i)),
                      theme=theme(plot.title=element_text(face="bold",size=13),
                                  plot.subtitle=element_text(size=9,color="grey40")))
    tryCatch(
      ggsave(file.path(output_ordner,paste0(igel_i,"_temporal.png")),
             p_g, width=14, height=14, dpi=150),
      error = function(e)
        cat(sprintf("  ! Plot-Fehler %s (6-Panel): %s\n", igel_i, conditionMessage(e)))
    )
  }
  cat(sprintf("  Plot: %s_temporal.png\n", igel_i))
}

# ── Übersichtsplot ────────────────────────────────────────────
if (nrow(nacht_tab) > 0) {
  igel_genug <- nacht_tab[, .N, by=igel][N>=5, igel]
  ns <- nacht_tab[igel %in% igel_genug]
  ns[, phase := factor(
    ifelse(tage_seit <= 7, "Woche 1\n(direkt nach\nAuswilderung)",
                           "Ab Woche 2\n(etabliert)"),
    levels=c("Woche 1\n(direkt nach\nAuswilderung)","Ab Woche 2\n(etabliert)"))]

  pu1 <- ggplot(ns,aes(x=tage_seit,y=weg_gesamt_m,color=igel,group=igel)) +
    geom_line(alpha=0.5,linewidth=0.7) +
    geom_smooth(aes(group=1),method="loess",formula=y~x,se=TRUE,color="black",linewidth=1.5,fill="grey80") +
    scale_color_viridis_d(option="turbo",guide="none") +
    labs(title="Weglaenge pro Nacht -- alle Igel",x="Tage",y="Weglaenge (m)") + theme_bew

  pu2 <- ggplot(ns,aes(x=tage_seit,y=dist_release_m,color=igel,group=igel)) +
    geom_line(alpha=0.5,linewidth=0.7) +
    geom_smooth(aes(group=1),method="loess",formula=y~x,se=TRUE,color="black",linewidth=1.5,fill="grey80") +
    scale_color_viridis_d(option="turbo",guide="none") +
    labs(title="Distanz Auswilderung -- alle Igel",x="Tage",y="Distanz (m)") + theme_bew

  pu3 <- ggplot(ns[!is.na(tortuositaet)],aes(x=tage_seit,y=tortuositaet,color=igel,group=igel)) +
    geom_line(alpha=0.4,linewidth=0.6) +
    geom_smooth(aes(group=1),method="loess",formula=y~x,se=TRUE,color="#2D6A4F",linewidth=1.5,
                fill="#95D5B2",alpha=0.3) +
    scale_color_viridis_d(option="turbo",guide="none") +
    scale_y_continuous(limits=c(0,1)) +
    labs(title="Tortuositaet -- alle Igel",x="Tage",y="Tortuositaet") + theme_bew

  pu4 <- ggplot(ns,aes(x=phase,y=weg_gesamt_m,fill=phase)) +
    geom_boxplot(alpha=0.7,outlier.shape=21,outlier.size=2) +
    geom_jitter(width=0.15,alpha=0.4,size=1) +
    scale_fill_manual(values=c("#ADB5BD","#4361EE"),guide="none") +
    labs(title="Weglaenge: frueh vs. etabliert",x=NULL,y="Weglaenge (m)") + theme_bew

  pu5 <- ggplot(ns[!is.na(tortuositaet)],aes(x=phase,y=tortuositaet,fill=phase)) +
    geom_boxplot(alpha=0.7,outlier.shape=21) +
    geom_jitter(width=0.15,alpha=0.4,size=1) +
    scale_fill_manual(values=c("#ADB5BD","#2D6A4F"),guide="none") +
    scale_y_continuous(limits=c(0,1)) +
    labs(title="Tortuositaet: frueh vs. etabliert",x=NULL,y="Tortuositaet") + theme_bew

  pu6 <- ggplot(ns,aes(x=phase,y=dist_release_m,fill=phase)) +
    geom_boxplot(alpha=0.7,outlier.shape=21) +
    geom_jitter(width=0.15,alpha=0.4,size=1) +
    scale_fill_manual(values=c("#ADB5BD","#E63946"),guide="none") +
    labs(title="Distanz Auswilderung: frueh vs. etabliert",x=NULL,y="Distanz (m)") + theme_bew

  p_ueb <- (pu1|pu2|pu3)/(pu4|pu5|pu6) +
    plot_annotation(title="Block 4d -- Alle Igel: Bewegungsveraenderung ueber Zeit",
                    subtitle=sprintf("N = %d Igel | Schwarze Linie = Gesamt-Trend (loess)",
                                     length(igel_genug)),
                    theme=theme(plot.title=element_text(face="bold",size=14),
                                plot.subtitle=element_text(size=10,color="grey40")))

  ggsave(file.path(output_ordner,"alle_igel_temporal_uebersicht.png"),
         p_ueb, width=18, height=12, dpi=150)
  cat("  Uebersichtsplot gespeichert\n")
}

# ══════════════════════════════════════════════════════════════
# EXCEL-EXPORT
# ══════════════════════════════════════════════════════════════
cat("\nErstelle Excel-Bericht...\n")

wb <- createWorkbook()
mk_s <- function(fg="FFFFFF",fc="000000",bold=FALSE,halign="center",size=10) {
  createStyle(fontName="Arial",fontSize=size,fontColour=paste0("#",fc),
               fgFill=paste0("#",fg),textDecoration=if(bold)"bold" else NULL,
               halign=halign,valign="center",
               border="TopBottomLeftRight",borderColour="#AABFD4")
}
hdr  <- mk_s("2C5F8A","FFFFFF",bold=TRUE)
alt1 <- mk_s("F5F9FF"); alt2 <- mk_s("FFFFFF")

# Sheet 1: Naechtliche Metriken
addWorksheet(wb,"Naechtliche Metriken")
lbl1 <- c("Igel","Datum","Tage seit Ausw.","N Fixes","Weglaenge (m)",
           "Netto-Disp. (m)","Max Dist Start (m)","Dist Ausw. (m)",
           "Tortuositaet","Drehwinkel (deg)","Schrittlaenge (m)",
           "SD Schritt (m)","Aktivitaetsradius 95% (m)")
writeData(wb,"Naechtliche Metriken",as.data.table(t(lbl1)),startRow=1,colNames=FALSE)
addStyle(wb,"Naechtliche Metriken",hdr,rows=1,cols=1:13,gridExpand=TRUE)
writeData(wb,"Naechtliche Metriken",nacht_tab,startRow=2,colNames=FALSE)
for (z in seq_len(nrow(nacht_tab)))
  addStyle(wb,"Naechtliche Metriken",if(z%%2==0)alt1 else alt2,rows=z+1,cols=1:13,gridExpand=TRUE)
setColWidths(wb,"Naechtliche Metriken",cols=1:13,
              widths=c(10,12,14,9,13,14,16,14,12,14,15,12,20))
freezePane(wb,"Naechtliche Metriken",firstRow=TRUE)

# Sheet 2: Hotspot-Verlauf
if (nrow(hotspot_tab) > 0) {
  addWorksheet(wb,"Hotspot-Verlauf")
  lbl2 <- c("Igel","Fenster Start","Fenster Mitte","Fenster Ende",
             "Tage seit Ausw.","N Fixes","Hotspot X","Hotspot Y",
             "Dist Ausw. (m)","Hotspot Flaeche (ha)","Stabilitaetsindex","Verschiebung (m)")
  hs_exp <- hotspot_tab[,.(igel,fenster_start,fenster_mitte,fenster_ende,
                             tage_seit,n_fixes_fen,hotspot_x,hotspot_y,
                             hotspot_dist_release_m,hotspot_ha,stabilitaet,hotspot_shift_m)]
  writeData(wb,"Hotspot-Verlauf",as.data.table(t(lbl2)),startRow=1,colNames=FALSE)
  addStyle(wb,"Hotspot-Verlauf",hdr,rows=1,cols=1:12,gridExpand=TRUE)
  writeData(wb,"Hotspot-Verlauf",hs_exp,startRow=2,colNames=FALSE)
  for (z in seq_len(nrow(hs_exp)))
    addStyle(wb,"Hotspot-Verlauf",if(z%%2==0)alt1 else alt2,rows=z+1,cols=1:12,gridExpand=TRUE)
  setColWidths(wb,"Hotspot-Verlauf",cols=1:12,widths=c(10,14,14,14,14,10,12,12,16,16,14,14))
  freezePane(wb,"Hotspot-Verlauf",firstRow=TRUE)
}

# Sheet 3: Zusammenfassung
addWorksheet(wb,"Zusammenfassung")
if (nrow(nacht_tab) > 0) {
  zus <- nacht_tab[,.(
    n_naechte      = .N,
    mean_weg_m     = round(mean(weg_gesamt_m,  na.rm=TRUE),1),
    sd_weg_m       = round(sd(weg_gesamt_m,    na.rm=TRUE),1),
    max_weg_m      = round(max(weg_gesamt_m,   na.rm=TRUE),1),
    mean_dist_m    = round(mean(dist_release_m,na.rm=TRUE),1),
    max_dist_m     = round(max(dist_release_m, na.rm=TRUE),1),
    mean_tort      = round(mean(tortuositaet,  na.rm=TRUE),3),
    mean_radius_m  = round(mean(radius_95_m,   na.rm=TRUE),1),
    trend_weg      = { x<-.SD$tage_seit;y<-.SD$weg_gesamt_m;v<-complete.cases(x,y)
                       if(sum(v)>=3) round(coef(lm(y[v]~x[v]))[2],2) else NA_real_ },
    trend_tort     = { x<-.SD$tage_seit;y<-.SD$tortuositaet;v<-complete.cases(x,y)
                       if(sum(v)>=3) round(coef(lm(y[v]~x[v]))[2],4) else NA_real_ }
  ),by=igel]
  zus[, interpret_weg  := fcase(!is.na(trend_weg) & trend_weg >  5, "Weglaenge nimmt zu",
                                  !is.na(trend_weg) & trend_weg < -5, "Weglaenge nimmt ab",
                                  default = "stabil")]
  zus[, interpret_tort := fcase(!is.na(trend_tort) & trend_tort >  0.005, "geradliniger",
                                  !is.na(trend_tort) & trend_tort < -0.005, "mehr Kurven",
                                  default = "stabil")]
  lbl3 <- c("Igel","N Naechte","Weg Mittel (m)","SD Weg","Max Weg",
             "Dist Ausw. Mittel (m)","Dist Ausw. Max (m)",
             "Tortuositaet Mittel","Aktivitaetsradius Mittel (m)",
             "Trend Weg","Trend Tort.","Interpret Weg","Interpret Tort.")
  writeData(wb,"Zusammenfassung",as.data.table(t(lbl3)),startRow=1,colNames=FALSE)
  addStyle(wb,"Zusammenfassung",hdr,rows=1,cols=1:13,gridExpand=TRUE)
  writeData(wb,"Zusammenfassung",zus,startRow=2,colNames=FALSE)
  for (z in seq_len(nrow(zus)))
    addStyle(wb,"Zusammenfassung",if(z%%2==0)alt1 else alt2,rows=z+1,cols=1:13,gridExpand=TRUE)
  setColWidths(wb,"Zusammenfassung",cols=1:13,
                widths=c(10,10,14,10,10,18,16,18,20,12,12,18,16))
  freezePane(wb,"Zusammenfassung",firstRow=TRUE)
}

saveWorkbook(wb, file.path(output_ordner,"temporal_ergebnisse.xlsx"), overwrite=TRUE)
cat("Excel gespeichert\n")

saveRDS(list(nacht_metriken=nacht_tab, hotspot_verlauf=hotspot_tab),
        file.path(output_ordner,"temporal_ergebnisse.rds"))
cat("RDS gespeichert\n")

# ══════════════════════════════════════════════════════════════
# ABSCHNITT 8: GRUPPENVERGLEICHE — Geschlecht, Reha-Dauer, Saison
# ══════════════════════════════════════════════════════════════
# Fragestellungen:
#   - Haben Maennchen groessere Aktionsraeume als Weibchen?
#   - Zeigen laenger rehabilitierte Tiere andere Bewegungsmuster?
#   - Unterscheidet sich das Verhalten zwischen Sommer (Foragieren)
#     und Herbst (Vorbereitung auf Winterschlaf)?
# ══════════════════════════════════════════════════════════════
cat("\n====================================\n")
cat("Abschnitt 8: Gruppenvergleiche\n")
cat("====================================\n\n")

# Sicherheitscheck: nacht_tab muss datum-Spalte haben
if (nrow(nacht_tab) == 0 || !"datum" %in% names(nacht_tab)) {
  cat("  Hinweis: nacht_tab leer oder datum fehlt.\n")
  cat("  Spaltennamen:", paste(names(nacht_tab), collapse=", "), "\n\n")
  nacht_tab <- data.table()
} else {
  nacht_tab[, datum := as.Date(datum)]  # sicherstellen dass datum ein Date ist
}

if (nrow(nacht_tab) > 0) {
# Saison aus Monat ableiten (biologisch sinnvoller als Kalender)
nacht_tab[, monat := month(datum)]
nacht_tab[, saison := fcase(
  monat %in% 4:5,   "Fruehjahr (Apr-Mai)",
  monat %in% 6:9,   "Sommer (Jun-Sep)",
  monat %in% 10:11, "Herbst/Vorwintersch. (Okt-Nov)",
  default           = "Winter (Dez-Marz)"
)]
nacht_tab[, saison := factor(saison, levels = c(
  "Fruehjahr (Apr-Mai)", "Sommer (Jun-Sep)",
  "Herbst/Vorwintersch. (Okt-Nov)", "Winter (Dez-Marz)"))]

# Reha-Dauer-Gruppe: kurz (< Median) vs. lang (>= Median)
if ("time_reha" %in% names(nacht_tab) && nacht_tab[!is.na(time_reha), .N] > 0) {
  med_reha <- median(nacht_tab$time_reha, na.rm = TRUE)
  nacht_tab[, reha_gruppe := fifelse(
    is.na(time_reha), NA_character_,
    fifelse(time_reha < med_reha,
            sprintf("Kurz (< %d Tage)", round(med_reha)),
            sprintf("Lang (>= %d Tage)", round(med_reha)))
  )]
  cat(sprintf("  Reha-Dauer Median: %.0f Tage\n", med_reha))
}

hat_sex   <- "sex"       %in% names(nacht_tab) && nacht_tab[!is.na(sex), .N] > 0
hat_reha  <- "reha_gruppe" %in% names(nacht_tab) && nacht_tab[!is.na(reha_gruppe), .N] > 0
hat_alter <- "alter"     %in% names(nacht_tab) && nacht_tab[!is.na(alter), .N] > 0
hat_saison <- nacht_tab[!is.na(saison) & saison != "Winter (Dez-Marz)", .N] > 0

# Farbpaletten
farbe_sex   <- c("m" = "#4575b4", "M" = "#4575b4", "maennlich" = "#4575b4",
                  "f" = "#d6604d", "F" = "#d6604d", "weiblich" = "#d6604d",
                  "w" = "#d6604d", "W" = "#d6604d")
farbe_saison <- c("Fruehjahr (Apr-Mai)"          = "#74add1",
                   "Sommer (Jun-Sep)"              = "#f46d43",
                   "Herbst/Vorwintersch. (Okt-Nov)"= "#2d6a4f",
                   "Winter (Dez-Marz)"             = "#4575b4")
farbe_reha  <- c("#fee090", "#d73027")

plot_liste_gruppe <- list()

# ── 8a: Geschlechtsvergleich ──────────────────────────────────
if (hat_sex) {
  sex_tab <- nacht_tab[!is.na(sex) & !is.na(radius_95_m)]
  n_sex   <- sex_tab[, .N, by = sex]
  cat("  Geschlechtsverteilung:
"); print(n_sex)

  # Wilcoxon-Test Geschlecht x Radius
  sex_vals <- unique(sex_tab$sex)
  if (length(sex_vals) == 2) {
    wx <- tryCatch(
      wilcox.test(radius_95_m ~ sex, data = sex_tab),
      error = function(e) NULL)
    if (!is.null(wx))
      cat(sprintf("  Wilcoxon Geschlecht x Radius: W=%.0f, p=%.4f\n",
                  wx$statistic, wx$p.value))
  }

  p_sex_rad <- ggplot(sex_tab, aes(x = sex, y = radius_95_m, fill = sex)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
    scale_fill_manual(values = farbe_sex, guide = "none") +
    labs(title = "Aktivitaetsradius nach Geschlecht",
         subtitle = "Maennchen haben typischerweise groessere Aktionsraeume",
         x = NULL, y = "95%-Aktivitaetsradius (m)") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))

  p_sex_dist <- ggplot(sex_tab[!is.na(dist_release_m)],
                        aes(x = sex, y = dist_release_m, fill = sex)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.5) +
    scale_fill_manual(values = farbe_sex, guide = "none") +
    labs(title = "Distanz Auswilderungsgehege nach Geschlecht",
         x = NULL, y = "Distanz (m)") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))

  p_sex_weg <- ggplot(sex_tab[!is.na(weg_gesamt_m) & weg_gesamt_m > 0],
                       aes(x = sex, y = weg_gesamt_m, fill = sex)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.5) +
    scale_fill_manual(values = farbe_sex, guide = "none") +
    labs(title = "Weglaenge nach Geschlecht",
         x = NULL, y = "Weglaenge (m)") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))

  p_sex <- p_sex_rad + p_sex_dist + p_sex_weg +
    plot_annotation(title = "Gruppenvergleich: Geschlecht",
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave(file.path(output_ordner, "gruppe_geschlecht.png"),
         p_sex, width = 15, height = 6, dpi = 150)
  plot_liste_gruppe[["sex"]] <- p_sex
  cat("  gruppe_geschlecht.png gespeichert\n")
}

# ── 8b: Saison-Vergleich ──────────────────────────────────────
if (hat_saison) {
  saison_tab <- nacht_tab[!is.na(saison) & !is.na(radius_95_m)]
  saison_n   <- saison_tab[, .N, by = saison][order(saison)]
  cat("\n  Saisonal verfuegbare Naechte:\n"); print(saison_n)

  p_sais_rad <- ggplot(saison_tab, aes(x = saison, y = radius_95_m, fill = saison)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.6) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.45) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
    scale_fill_manual(values = farbe_saison, guide = "none") +
    labs(title = "Aktivitaetsradius je Saison",
         subtitle = "Mehr Nahrungssuche im Sommer? Mehr Bewegung vor Winterschlaf?",
         x = NULL, y = "95%-Aktivitaetsradius (m)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 20, hjust = 1, size = 9))

  p_sais_weg <- ggplot(saison_tab[!is.na(weg_gesamt_m) & weg_gesamt_m > 0],
                        aes(x = saison, y = weg_gesamt_m, fill = saison)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.6) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.45) +
    scale_fill_manual(values = farbe_saison, guide = "none") +
    labs(title = "Weglaenge je Saison",
         x = NULL, y = "Weglaenge (m)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 20, hjust = 1, size = 9))

  # Zeitverlauf: Median-Radius pro Monat (alle Tiere)
  monat_tab <- nacht_tab[!is.na(radius_95_m), .(
    median_rad  = median(radius_95_m, na.rm = TRUE),
    median_weg  = median(weg_gesamt_m, na.rm = TRUE),
    n_naechte   = .N
  ), by = .(monat, saison)]
  setorder(monat_tab, monat)

  p_sais_trend <- ggplot(monat_tab[n_naechte >= 2],
                          aes(x = monat, y = median_rad, color = saison, size = n_naechte)) +
    geom_line(color = "grey40", linewidth = 0.8) +
    geom_point() +
    scale_color_manual(values = farbe_saison, name = "Saison") +
    scale_size_continuous(name = "N Naechte", range = c(2, 7)) +
    scale_x_continuous(breaks = 1:12,
                        labels = c("Jan","Feb","Mar","Apr","Mai","Jun",
                                   "Jul","Aug","Sep","Okt","Nov","Dez")) +
    labs(title = "Monatlicher Median-Aktivitaetsradius",
         subtitle = "Groesse = Anzahl Beobachtungsnaechte",
         x = NULL, y = "Median 95%-Radius (m)") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))

  p_saison <- (p_sais_rad + p_sais_weg) / p_sais_trend +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(title = "Gruppenvergleich: Saison",
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave(file.path(output_ordner, "gruppe_saison.png"),
         p_saison, width = 15, height = 10, dpi = 150)
  plot_liste_gruppe[["saison"]] <- p_saison
  cat("  gruppe_saison.png gespeichert\n")
}

# ── 8c: Reha-Dauer-Vergleich ──────────────────────────────────
if (hat_reha) {
  reha_tab <- nacht_tab[!is.na(reha_gruppe) & !is.na(radius_95_m)]

  wx_reha <- tryCatch(
    wilcox.test(radius_95_m ~ reha_gruppe, data = reha_tab),
    error = function(e) NULL)
  if (!is.null(wx_reha))
    cat(sprintf("  Wilcoxon Reha-Dauer x Radius: W=%.0f, p=%.4f\n",
                wx_reha$statistic, wx_reha$p.value))

  # Korrelation Reha-Dauer (numerisch) x Radius
  if ("time_reha" %in% names(reha_tab)) {
    cor_reha <- tryCatch(
      cor.test(reha_tab$time_reha, reha_tab$radius_95_m, method = "spearman"),
      error = function(e) NULL)
    if (!is.null(cor_reha))
      cat(sprintf("  Spearman Reha-Dauer x Radius: rho=%.3f, p=%.4f\n",
                  cor_reha$estimate, cor_reha$p.value))
  }

  p_reha_rad <- ggplot(reha_tab, aes(x = reha_gruppe, y = radius_95_m, fill = reha_gruppe)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
    scale_fill_manual(values = setNames(farbe_reha, unique(reha_tab$reha_gruppe)),
                      guide = "none") +
    labs(title = "Aktivitaetsradius nach Reha-Dauer",
         subtitle = "Laengere Reha = staerkere Habituation = andere Bewegungsmuster?",
         x = NULL, y = "95%-Aktivitaetsradius (m)") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))

  # Streudiagramm: Reha-Dauer (numerisch) x Median-Radius pro Tier
  if ("time_reha" %in% names(reha_tab)) {
    per_tier <- reha_tab[, .(
      med_rad  = median(radius_95_m, na.rm = TRUE),
      time_reha = first(time_reha),
      sex       = first(sex)
    ), by = igel]

    p_reha_scatter <- ggplot(per_tier[!is.na(time_reha)],
                              aes(x = time_reha, y = med_rad)) +
      geom_point(aes(color = if (hat_sex) sex else "alle"), size = 3, alpha = 0.8) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30",
                  fill = "grey80", alpha = 0.2, linewidth = 0.9) +
      { if (hat_sex) scale_color_manual(values = farbe_sex, name = "Geschlecht")
        else scale_color_manual(values = c("alle" = "#4575b4"), guide = "none") } +
      labs(title = "Reha-Dauer vs. median. Aktivitaetsradius (pro Tier)",
           x = "Reha-Dauer (Tage)", y = "Median 95%-Radius (m)") +
      theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))
  } else {
    p_reha_scatter <- ggplot() + theme_void()
  }

  p_reha <- p_reha_rad + p_reha_scatter +
    plot_annotation(title = "Gruppenvergleich: Rehabilitationsdauer",
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave(file.path(output_ordner, "gruppe_reha_dauer.png"),
         p_reha, width = 13, height = 6, dpi = 150)
  plot_liste_gruppe[["reha"]] <- p_reha
  cat("  gruppe_reha_dauer.png gespeichert\n")
}

  cat(sprintf("\n  Gruppenvergleiche abgeschlossen: %d Plots erstellt\n\n",
              length(plot_liste_gruppe)))
} # end if nrow(nacht_tab) > 0



# ══════════════════════════════════════════════════════════════
# ABSCHNITT 9: SIGNALVERLUST & DISPERSAL-ERKENNUNG
# ══════════════════════════════════════════════════════════════
# Methode: Direkte Nutzung der Activity-Klassifikationsdaten
#   (data/activity/classification_active_passive_IgelX_*.csv)
#   Jede Zeile = 1 Detektionsereignis mit Klassifikation:
#     "a" = aktiv (Tier bewegt sich)
#     "p" = passiv (Tier sitzt still / moeglicherweise Sender)
#
# Unterscheidung Sender abgefallen vs. Winterschlaf:
#   Sender abgefallen: nur "p", HOHE Detektionsrate (Sender auf Boden)
#   Winterschlaf:      nur "p", NIEDRIGE Rate (Antenne unter Erde)
#                      + saisonales Muster (Okt-Nov -> 0% aktiv -> Marz)
#   Abgewandert:       Signal verschwindet ausserhalb Winterschlafzeit
# ══════════════════════════════════════════════════════════════
cat("====================================\n")
cat("Abschnitt 9: Signalverlust & Dispersal\n")
cat("====================================\n\n")

# path.expand() noetig: list.files() expandiert ~ nicht automatisch
activity_ordner <- path.expand(file.path(projekt_root, "data", "activity"))

if (!dir.exists(activity_ordner)) {
  cat(sprintf("  ! Activity-Ordner nicht gefunden: %s\n", activity_ordner))
  cat("  Bitte Pfad in projekt_root pruefen\n\n")
} else {
  cat(sprintf("  Activity-Ordner: %s\n", activity_ordner))
  cat(sprintf("  Gefundene Dateien: %d\n",
              length(list.files(activity_ordner, pattern = "\\.csv$"))))
}

# Activity-CSV laden (eine pro Tier)
# Beruecksichtigt:
#   - Igel10/11: 10 Spalten (nutzt pred_nested_loio_wmv statt smoothed_wmv)
#   - Igel3: hat Frequenz im Dateinamen (Igel3_150.080_...)
#
# BUG-FIX 1: act[[col]] statt .SD[[col]] — expliziter und stabiler in allen
#   data.table-Versionen. .SD[[col]] in := kann in aelteren Versionen leer sein.
# BUG-FIX 2: Diagnostik-Output wenn Datei nicht gefunden oder kein Ergebnis.
# BUG-FIX 3: Guard danach fuer leeres activity_raw (kein Spaltennamen-Crash).
lade_activity <- function(igel_name) {
  f <- list.files(activity_ordner,
                  pattern = paste0("classification_active_passive_", igel_name, "[_.]"),
                  full.names = TRUE)
  if (length(f) == 0) {
    message(sprintf("    Keine Activity-Datei fuer %s in: %s", igel_name, activity_ordner))
    return(NULL)
  }

  tryCatch({
    act <- fread(f[1], showProgress = FALSE, na.strings = c("", "NA"))

    # Zeitstempel-Spalte finden
    t_col <- intersect(c("X_time", "_time", "time", "timestamp"), names(act))
    if (length(t_col) == 0) {
      message(sprintf("    %s: Keine Zeitstempel-Spalte gefunden (erwartet: _time)", igel_name))
      return(NULL)
    }

    # Klassifikations-Spalte: bevorzuge smoothed_wmv (gefiltert), sonst wmv
    k_candidates <- names(act)[grepl("smoothed_wmv|_wmv$|smoothed$", names(act))]
    if (length(k_candidates) == 0)
      k_candidates <- names(act)[grepl("pred_nested", names(act))]
    if (length(k_candidates) == 0) {
      message(sprintf("    %s: Keine Klassifikationsspalte gefunden", igel_name))
      return(NULL)
    }
    k_col <- tail(k_candidates, 1)  # letzter Treffer = am staerksten verarbeitet

    # BUG-FIX: act[[col]] statt .SD[[col]] — direkter Spaltenzugriff, keine
    # Abhaengigkeit vom data.table-internen .SD-Verhalten in := Ausdruecken
    time_vals  <- as.character(act[[t_col[1]]])
    klasse_vals <- act[[k_col]]

    result <- data.table(
      igel  = igel_name,
      datum = as.Date(substr(time_vals, 1L, 10L)),
      aktiv = (klasse_vals == "a")
    )

    result <- result[!is.na(datum) & !is.na(aktiv)]

    if (nrow(result) == 0L) {
      message(sprintf("    %s: 0 valide Zeilen nach Filter (Datum/Aktiv NA?)", igel_name))
      return(NULL)
    }
    result

  }, error = function(e) {
    message(sprintf("    ! lade_activity(%s): %s", igel_name, conditionMessage(e)))
    NULL
  })
}

cat("  Lade Activity-Daten fuer alle Tiere...\n")

# Diagnostik: welche igel-Namen sind in alle_fixes?
igel_fuer_activity <- sort(unique(alle_fixes$igel))
cat(sprintf("  Suche Activity-Daten fuer %d Tiere: %s\n",
            length(igel_fuer_activity),
            paste(igel_fuer_activity, collapse = ", ")))

activity_liste <- Filter(Negate(is.null),
                          lapply(igel_fuer_activity, lade_activity))

activity_raw <- if (length(activity_liste) > 0) {
  rbindlist(activity_liste, fill = TRUE)
} else {
  data.table(igel = character(), datum = as.Date(character()),
             aktiv = logical())
}

# BUG-FIX: Guard fuer leeres activity_raw — verhindert "object 'igel' not found"
# wenn lade_activity() fuer alle Tiere NULL zurueckgibt
if (nrow(activity_raw) == 0L || !"igel" %in% names(activity_raw)) {
  cat(sprintf("  [WARN] Keine Activity-Daten geladen (%d/%d Tiere). Abschnitt 9 wird uebersprungen.\n",
              length(activity_liste), length(igel_fuer_activity)))
  cat("  Tipp: Pruefe ob die Meldungen oben (Keine Activity-Datei / Fehler) Hinweise geben.\n\n")
  activity_daily  <- data.table(igel = character(), datum = as.Date(character()),
                                 n_gesamt = integer(), n_aktiv = integer(), pct_aktiv = numeric())
  activity_monthly <- data.table(igel = character(), monat = character(),
                                  pct_aktiv_monat = numeric(), n_tage = integer())
  signal_analyse  <- data.table(igel = character(), first_date = as.Date(character()),
                                 last_date = as.Date(character()), klassifikation = character())
} else {
  cat(sprintf("  %d Tiere mit Activity-Daten | %s bis %s\n\n",
              uniqueN(activity_raw$igel),
              min(activity_raw$datum, na.rm = TRUE),
              max(activity_raw$datum, na.rm = TRUE)))

# Tagesweise Zusammenfassung: % aktiv pro Tag
activity_daily <- activity_raw[, .(
  n_gesamt  = .N,
  n_aktiv   = sum(aktiv, na.rm = TRUE),
  pct_aktiv = round(100 * sum(aktiv, na.rm = TRUE) / .N, 1)
), by = .(igel, datum)]

# Monatliche Zusammenfassung (fuer Winterschlaf-Erkennung)
activity_monthly <- activity_daily[, .(
  pct_aktiv_monat = round(mean(pct_aktiv, na.rm = TRUE), 1),
  n_tage          = .N
), by = .(igel, monat = format(datum, "%Y-%m"))]

# Pro Tier: Signalverlust-Klassifikation
signal_analyse <- activity_daily[, {
  setorder(.SD, datum)
  first_date <- min(datum)
  last_date  <- max(datum)

  # Winterschlaf-Signal: mind. 2 aufeinanderfolgende Monate mit 0% aktiv
  monat_stats <- .SD[, .(pct = mean(pct_aktiv)), by = .(m = format(datum, "%Y-%m"))]
  setorder(monat_stats, m)
  passive_monate <- monat_stats[pct < 2, .N]  # Monate mit < 2% Aktivitaet

  # Konsekutive passive Monate pruefen
  hat_winterschlaf <- FALSE
  if (nrow(monat_stats) >= 2) {
    runs <- rle(monat_stats$pct < 2)
    hat_winterschlaf <- any(runs$values & runs$lengths >= 2)
  }

  # Wiederaufwachen nach Winterschlaf?
  hat_reaktivierung <- hat_winterschlaf && {
    letzter_monat_pct <- tail(monat_stats$pct, 1)
    letzter_monat_pct > 20  # letzter Monat wieder aktiv
  }

  # Endphase: letzte 14 Tage
  letzte14 <- .SD[datum >= last_date - 14]
  pct_ende <- round(mean(letzte14$pct_aktiv, na.rm = TRUE), 1)
  n_ende   <- sum(letzte14$n_gesamt)

  # Detektionsrate: frueher vs. spaeter (Trend)
  n_monate <- nrow(monat_stats)
  rate_ratio <- if (n_monate >= 2) {
    frueh <- mean(head(monat_stats$pct, ceiling(n_monate/2)), na.rm = TRUE)
    spaet <- mean(tail(monat_stats$pct, ceiling(n_monate/2)), na.rm = TRUE)
    if (frueh > 0) round(spaet / frueh, 2) else NA_real_
  } else NA_real_

  analyse_ende <- max(activity_raw$datum, na.rm = TRUE)
  tage_seit_letztem <- as.numeric(analyse_ende - last_date)

  # Klassifikation mit biologischer Begruendung
  klassifikation <- fcase(
    # Noch aktiv: letzter Datenpunkt juenger als 14 Tage
    tage_seit_letztem <= 14,
      "Noch aktiv (aktuell getrackt)",
    # Winterschlaf mit Reaktivierung = Tier hat ueberwintert und lebt noch
    hat_reaktivierung,
      "Ueberwinterung + Reaktivierung (Tier lebt)",
    # Winterschlaf ohne Reaktivierung = Tier nicht mehr detektiert nach Winterschlaf
    hat_winterschlaf,
      "Winterschlaf (kein Signal mehr danach)",
    # Nur passiv am Ende, hohe Detektionsrate = Sender abgefallen
    pct_ende < 2 & n_ende > 5000,
      "Sender abgefallen (nur passiv, hohe Rate)",
    # Signal weg, kein saisonales Muster = abgewandert
    tage_seit_letztem > 30 & !hat_winterschlaf,
      "Abgewandert / Signal verloren",
    default = "Unbekannt"
  )

  .(
    first_date           = first_date,
    last_date            = last_date,
    tracking_tage        = as.numeric(last_date - first_date),
    tage_seit_letztem    = round(tage_seit_letztem),
    passive_monate       = passive_monate,
    hat_winterschlaf     = hat_winterschlaf,
    hat_reaktivierung    = hat_reaktivierung,
    pct_aktiv_gesamt     = round(mean(pct_aktiv, na.rm = TRUE), 1),
    pct_aktiv_ende14     = pct_ende,
    rate_ratio           = rate_ratio,
    klassifikation       = klassifikation
  )
}, by = igel]

setorder(signal_analyse, last_date)
cat("  Signalverlust-Klassifikation:\n\n")
print(signal_analyse[, .(igel, first_date, last_date, hat_winterschlaf,
                          hat_reaktivierung, pct_aktiv_ende14, klassifikation)])
cat("\n")
cat("  Zusammenfassung:\n")
print(signal_analyse[, .N, by = klassifikation][order(-N)])
cat("\n")

# ── Plot S1: Monatliche Aktivitaetsrate (Heatmap) ─────────────
# Zeigt Winterschlaf und Reaktivierung sehr klar
activity_monthly_plot <- activity_monthly[igel %in% signal_analyse$igel]
activity_monthly_plot[, monat_date := as.Date(paste0(monat, "-01"))]
activity_monthly_plot <- merge(activity_monthly_plot,
                                signal_analyse[, .(igel, klassifikation)],
                                by = "igel", all.x = TRUE)

p_heatmap <- ggplot(activity_monthly_plot,
                     aes(x = monat_date, y = reorder(igel, monat_date),
                         fill = pct_aktiv_monat)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#4575b4", mid = "#fee090", high = "#d73027",
                        midpoint = 30, name = "% aktiv",
                        limits = c(0, 100), na.value = "grey90") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  labs(title = "Monatliche Aktivitaetsrate aller Igel",
       subtitle = "Blau = passiv/Winterschlaf | Gelb-Rot = aktiv | Grau = keine Daten",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title      = element_text(face = "bold", size = 11),
        axis.text.y     = element_text(size = 8),
        legend.position = "right",
        panel.grid      = element_blank())

# ── Plot S2: Zeitstrahl mit Klassifikation ────────────────────
sig_plot <- signal_analyse[order(first_date)]
sig_plot[, igel_f := factor(igel, levels = igel)]

farben_sig <- c(
  "Noch aktiv (aktuell getrackt)"           = "#2d6a4f",
  "Ueberwinterung + Reaktivierung (Tier lebt)" = "#4575b4",
  "Winterschlaf (kein Signal mehr danach)"   = "#74add1",
  "Sender abgefallen (nur passiv, hohe Rate)"= "#e76f51",
  "Abgewandert / Signal verloren"            = "#d73027",
  "Unbekannt"                                = "grey60"
)

p_zeitstrahl <- ggplot(sig_plot, aes(y = igel_f)) +
  geom_segment(aes(x = first_date, xend = last_date,
                   y = igel_f, yend = igel_f,
                   color = klassifikation), linewidth = 3.5, alpha = 0.85) +
  geom_point(aes(x = last_date, color = klassifikation),
             size = 3.5, shape = 21, fill = "white", stroke = 1.8) +
  scale_color_manual(values = farben_sig, name = NULL) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Tracking-Zeitraeume und Signal-Status",
       subtitle = "Balken = Trackingperiode aus Activity-Daten | Kreis = letzter Detektionszeitpunkt",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        axis.text.x   = element_text(angle = 30, hjust = 1))

p_signal_gesamt <- p_zeitstrahl / p_heatmap +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title    = "Signalverlust-Analyse: Wann und wie endete das Tracking?",
    subtitle = "Basierend auf Activity-Klassifikationsdaten (aktiv/passiv pro Tier und Tag)",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(output_ordner, "signal_klassifikation.png"),
       p_signal_gesamt, width = 16, height = 12, dpi = 150)
cat("  signal_klassifikation.png gespeichert\n")

# Excel: Signalanalyse-Tabelle
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb_sig <- tryCatch(
    loadWorkbook(file.path(output_ordner, "temporal_ergebnisse.xlsx")),
    error = function(e) {
      wb <- createWorkbook()
      wb
    }
  )
  if (!"Signalverlust" %in% names(wb_sig)) addWorksheet(wb_sig, "Signalverlust")
  hdr_s <- createStyle(fontName="Arial", fontSize=10, fontColour="white",
                        fgFill="#2C5F8A", halign="center", textDecoration="bold",
                        border="TopBottomLeftRight", borderColour="#AABFD4")
  lbl_s <- c("Igel","Erster Tag","Letzter Tag","Tracking (Tage)",
             "Tage seit letzt. Detektion","Passive Monate",
             "Winterschlaf","Reaktivierung",
             "% aktiv gesamt","% aktiv letzte 14 Tage","Rate-Verhaeltnis","Klassifikation")
  writeData(wb_sig, "Signalverlust", as.data.table(t(lbl_s)),
            startRow=1, colNames=FALSE)
  addStyle(wb_sig, "Signalverlust", hdr_s, rows=1, cols=1:12, gridExpand=TRUE)
  writeData(wb_sig, "Signalverlust", signal_analyse, startRow=2, colNames=FALSE)
  setColWidths(wb_sig, "Signalverlust", cols=1:12,
               widths=c(10,12,12,14,20,14,12,12,16,20,14,32))
  saveWorkbook(wb_sig, file.path(output_ordner,"temporal_ergebnisse.xlsx"),
               overwrite=TRUE)
  cat("  Signalverlust-Sheet gespeichert\n")
}

cat("\n====================================\n")
cat("Abschnitt 9 abgeschlossen\n")
cat("====================================\n\n")

} # end else (activity_raw nicht leer)

# ══════════════════════════════════════════════════════════════
# WORD-BERICHT (Methodik + Ergebnisse)
# ══════════════════════════════════════════════════════════════
cat("\nErstelle Word-Bericht...\n")

tryCatch({
  # Kennzahlen berechnen
  if (!exists("zus")) {
    zus <- nacht_tab[, .(
      n_naechte     = .N,
      mean_weg_m    = round(mean(weg_gesamt_m,   na.rm=TRUE), 1),
      mean_dist_m   = round(mean(dist_release_m, na.rm=TRUE), 1),
      max_dist_m    = round(max(dist_release_m,  na.rm=TRUE), 1),
      mean_tort     = round(mean(tortuositaet,   na.rm=TRUE), 3),
      mean_radius_m = round(mean(radius_95_m,    na.rm=TRUE), 1)
    ), by = igel]
  }

  n_igel_r <- uniqueN(nacht_tab$igel)
  med_weg  <- round(median(nacht_tab$weg_gesamt_m,   na.rm=TRUE))
  med_dist <- round(median(nacht_tab$dist_release_m, na.rm=TRUE))
  max_dist <- round(max(nacht_tab$dist_release_m,    na.rm=TRUE))
  med_tort <- round(median(nacht_tab$tortuositaet,   na.rm=TRUE), 3)
  med_rad  <- round(median(nacht_tab$radius_95_m,    na.rm=TRUE))
  dist_w1  <- round(median(nacht_tab[tage_seit <= 7, dist_release_m], na.rm=TRUE))
  dist_w2  <- round(median(nacht_tab[tage_seit  > 7, dist_release_m], na.rm=TRUE))

  # Hilfsfunktionen
  h1 <- function(txt) officer::fpar(
    officer::ftext(txt, officer::fp_text(bold=TRUE, font.size=14, color="#2C5F8A",
                                          font.family="Arial")),
    fp_p = officer::fp_par(space_before=240, space_after=80,
                             border.bottom=officer::fp_border(color="#2C5F8A", width=1.5)))
  h2 <- function(txt) officer::fpar(
    officer::ftext(txt, officer::fp_text(bold=TRUE, font.size=12, color="#1a1a1a",
                                          font.family="Arial")),
    fp_p = officer::fp_par(space_before=160, space_after=60))
  txt <- function(t) officer::fpar(
    officer::ftext(t, officer::fp_text(font.size=11, font.family="Arial")),
    fp_p = officer::fp_par(space_before=0, space_after=80, line_spacing=1.2))
  sp  <- function() officer::fpar(officer::ftext(" "),
    fp_p = officer::fp_par(space_before=0, space_after=40))

  doc <- officer::read_docx()

  # Titel
  doc <- officer::body_add_fpar(doc, officer::fpar(
    officer::ftext("Temporale Bewegungsanalyse rehabilitierter Igel",
                   officer::fp_text(bold=TRUE, font.size=18, color="#2C5F8A",
                                     font.family="Arial")),
    fp_p = officer::fp_par(text.align="center", space_before=360, space_after=100)))
  doc <- officer::body_add_fpar(doc, officer::fpar(
    officer::ftext(paste0("Igelbesenderung Sachsenhagen | ",
                           format(Sys.Date(), "%B %Y")),
                   officer::fp_text(font.size=11, color="#666666",
                                     font.family="Arial", italic=TRUE)),
    fp_p = officer::fp_par(text.align="center", space_after=60)))
  doc <- officer::body_add_fpar(doc, officer::fpar(
    officer::ftext("Natalie Steiner, ITAW, Stiftung Tieraerztliche Hochschule Hannover",
                   officer::fp_text(font.size=11, color="#666666",
                                     font.family="Arial", italic=TRUE)),
    fp_p = officer::fp_par(text.align="center", space_after=400)))
  doc <- officer::body_add_break(doc)

  # ── 1. METHODIK ──────────────────────────────────────────────
  doc <- officer::body_add_fpar(doc, h1("1. Methodik"))

  doc <- officer::body_add_fpar(doc, h2("Daten und Studiendesign"))
  doc <- officer::body_add_fpar(doc, txt(paste0(
    "Die Daten stammen aus dem automatisierten VHF-Telemetriesystem tRackIT, ",
    "das am Auswilderungsgelaende Sachsenhagen installiert wurde. Das System ",
    "erfasst Positionen rehabilitierter Igel (Erinaceus europaeus) mit sub-minutlicher ",
    "Zeitaufloesung ueber mehrere Wochen nach der Auswilderung. ",
    "In dieser Analyse wurden ausschliesslich Nacht-Fixes verwendet ",
    "(Spalte Night == TRUE), da Igel streng nachtaktiv sind und tagsueber ",
    "in der Regel keine Bewegungsaktivitaet zeigen. ",
    sprintf("Insgesamt wurden %d Individuen mit %d Nacht-Datensaetzen analysiert.",
            n_igel_r, nrow(nacht_tab)))))

  doc <- officer::body_add_fpar(doc, sp())
  doc <- officer::body_add_fpar(doc, h2("Berechnete Metriken pro Nacht"))

  doc <- officer::body_add_fpar(doc, txt(paste0(
    "Weglaenge (m): Die pro Nacht zurueckgelegte Gesamtstrecke als Summe aller ",
    "Euklidischen Distanzen zwischen aufeinanderfolgenden Fixes. Schritte mit ",
    "Zeitluecken > 30 Minuten (Empfangsluecken) werden ausgeschlossen.")))

  doc <- officer::body_add_fpar(doc, txt(paste0(
    "Tortuositaet (Straightness Index, 0-1): Verhaeltnis zwischen der Luftlinie ",
    "(erster bis letzter Fix der Nacht) und der tatsaechlichen Weglaenge. ",
    "Ein Wert nahe 1 bedeutet geradlinige Bewegung, Werte nahe 0 beschreiben ",
    "eine stark gewundene, kurvenreiche Bewegung, wie sie beim Nahrungssuchen ",
    "typisch ist (Benhamou 2004). Igel zeigen typischerweise Werte unter 0.10.")))

  doc <- officer::body_add_fpar(doc, txt(paste0(
    "Aktivitaetsradius 95% (m): Das 95. Perzentil aller Distanzen vom naechtlichen ",
    "Schwerpunkt (Zentroid aller Nacht-Fixes). Gibt an, in welchem Radius um ",
    "den Schwerpunkt 95% der Aktivitaet stattfand.")))

  doc <- officer::body_add_fpar(doc, txt(paste0(
    "Distanz zum Auswilderungsgehege (m): Euklidische Distanz des Nacht-Zentroids ",
    "(Medianposition aller Fixes der Nacht) zum Auswilderungsgehege. ",
    "Zeigt, wie weit sich ein Tier im Laufe der Zeit vom Ausgangspunkt entfernt.")))

  doc <- officer::body_add_fpar(doc, sp())
  doc <- officer::body_add_fpar(doc, h2("Hotspot-Analyse"))
  doc <- officer::body_add_fpar(doc, txt(paste0(
    "Mit einem gleitenden ", FENSTER_N, "-Naechte-Fenster wurde fuer jedes ",
    "Zeitfenster ein KDE-Raster berechnet (R-Paket adehabitatHR, Bandbreite href). ",
    "Als Hotspot gilt der gewichtete Schwerpunkt der oberen ", HOTSPOT_PCT,
    "% der Dichtezellen. Der Stabilitaetsindex beschreibt die raeumliche ",
    "Ueberlappung aufeinanderfolgender Fenster (0 = vollstaendig verlagert, ",
    "1 = unveraendert). Koordinatensystem: UTM Zone 32N (EPSG:25832).")))

  doc <- officer::body_add_break(doc)

  # ── 2. ERGEBNISSE ──────────────────────────────────────────
  doc <- officer::body_add_fpar(doc, h1("2. Ergebnisse"))

  doc <- officer::body_add_fpar(doc, h2("Distanz zum Auswilderungsgehege"))
  doc <- officer::body_add_fpar(doc, txt(sprintf(paste0(
    "Die Tiere blieben insgesamt in unmittelbarer Naehe des Auswilderungsgeheges. ",
    "Die mediane Nacht-Zentroid-Distanz betrug %d m (Bereich: %d - %d m). ",
    "Kein Tier entfernte sich weiter als %d m vom Ausgangspunkt. ",
    "In der ersten Woche nach Auswilderung lag die mediane Distanz bei %d m, ",
    "ab Woche zwei bei %d m. Dieses Muster deutet auf eine graduelle Ausdehnung ",
    "des Aktionsraums hin, ohne jedoch zu einer weitraeumigen Dispersion zu fuehren."),
    med_dist,
    round(min(nacht_tab$dist_release_m, na.rm=TRUE)),
    max_dist, max_dist, dist_w1, dist_w2)))

  doc <- officer::body_add_fpar(doc, sp())
  doc <- officer::body_add_fpar(doc, h2("Tortuositaet"))
  doc <- officer::body_add_fpar(doc, txt(sprintf(paste0(
    "Der mediane Straightness Index betrug %.3f (Bereich: %.3f - %.3f). ",
    "Alle Tiere zeigten Werte weit unter 0.10, was ein hochgradig kurvenreiches ",
    "Bewegungsmuster belegt. Dies ist typisch fuer die Nahrungssuchstrategie ",
    "insektivorer Sauger und deutet auf normales Foragierverhalten ",
    "bereits kurz nach der Auswilderung hin."),
    med_tort,
    round(min(nacht_tab$tortuositaet, na.rm=TRUE), 3),
    round(max(nacht_tab$tortuositaet, na.rm=TRUE), 3))))

  doc <- officer::body_add_fpar(doc, sp())
  doc <- officer::body_add_fpar(doc, h2("Aktivitaetsradius"))
  doc <- officer::body_add_fpar(doc, txt(sprintf(paste0(
    "Der mediane naechtliche Aktivitaetsradius (95. Pz.) betrug %d m ",
    "(Bereich: %d - %d m). Die Variabilitaet zwischen Individuen war hoch, ",
    "was individuelle Unterschiede in der Raumnutzung widerspiegelt."),
    med_rad,
    round(min(nacht_tab$radius_95_m, na.rm=TRUE)),
    round(max(nacht_tab$radius_95_m, na.rm=TRUE)))))

  doc <- officer::body_add_fpar(doc, sp())
  doc <- officer::body_add_fpar(doc, h2("Eingewoehnungsphase vs. etablierte Phase"))
  doc <- officer::body_add_fpar(doc, txt(paste0(
    "Der Vergleich zwischen Woche 1 (Tage 0-7) und den Folgewochen (ab Tag 8) ",
    "zeigt keine systematische Veraenderung der raeumlichen Nutzungsmuster. ",
    "Distanz und Aktivitaetsradius blieben zwischen den Phasen vergleichbar. ",
    "Dies deutet darauf hin, dass die Tiere ihr Aktionsareal bereits ",
    "ab der ersten Nacht etablierten und kein laengerer Eingewoehnungszeitraum ",
    "beobachtet werden konnte.")))

  doc <- officer::body_add_fpar(doc, sp())
  doc <- officer::body_add_fpar(doc, h2("Individuelle Ergebnisse"))
  doc <- officer::body_add_fpar(doc, txt(paste0(
    "Tabelle 1 zeigt die Zusammenfassung der Bewegungsmetriken pro Individuum. ",
    "Dist. Ausw. = mittlere und maximale Distanz zum Auswilderungsgehege; ",
    "Tortuositaet = Straightness Index (0-1); Radius = 95. Pz. der Distanzen ",
    "vom naechtlichen Schwerpunkt.")))
  doc <- officer::body_add_fpar(doc, sp())

  # Ergebnistabelle
  tab_df <- as.data.frame(zus[order(igel), .(
    Igel             = igel,
    `N Naechte`      = n_naechte,
    `Weglaenge (m)`  = mean_weg_m,
    `Dist. Ausw. (m)`= mean_dist_m,
    `Max Dist. (m)`  = max_dist_m,
    `Tortuositaet`   = mean_tort,
    `Radius 95% (m)` = mean_radius_m
  )])

  ft <- flextable::flextable(tab_df) |>
    flextable::theme_vanilla() |>
    flextable::bg(bg = "#2C5F8A", part = "header") |>
    flextable::color(color = "white", part = "header") |>
    flextable::bold(part = "header") |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(j = 1, align = "left", part = "all") |>
    flextable::autofit()
  doc <- flextable::body_add_flextable(doc, ft)

  docx_pfad <- file.path(output_ordner, "Block4d_Methodik_Ergebnisse.docx")
  print(doc, target = docx_pfad)
  cat(sprintf("  Word-Bericht gespeichert: %s\n", basename(docx_pfad)))

}, error = function(e) {
  cat(sprintf("  ! Word-Bericht fehlgeschlagen: %s\n", conditionMessage(e)))
  cat("  Tipp: officer und flextable muessen installiert sein.\n")
})

# ══════════════════════════════════════════════════════════════
# ERGEBNIS-ZUSAMMENFASSUNG (wird nach jedem Lauf automatisch ausgegeben)
# ══════════════════════════════════════════════════════════════
cat("\n====================================\n")
cat("ERGEBNIS-ZUSAMMENFASSUNG Block 4d\n")
cat("====================================\n\n")

if (nrow(nacht_tab) > 0) {
  cat(sprintf("Analysierte Individuen:  %d Igel, %d Nacht-Eintraege\n\n",
              uniqueN(nacht_tab$igel), nrow(nacht_tab)))

  cat("── Weglaenge pro Nacht ──────────────────\n")
  cat(sprintf("  Median (alle):   %.0f m\n",  median(nacht_tab$weg_gesamt_m, na.rm=TRUE)))
  cat(sprintf("  Bereich:         %.0f – %.0f m\n",
              min(nacht_tab$weg_gesamt_m, na.rm=TRUE),
              max(nacht_tab$weg_gesamt_m, na.rm=TRUE)))
  cat(sprintf("  Woche 1:         %.0f m (Median)\n",
              nacht_tab[tage_seit<=7, median(weg_gesamt_m, na.rm=TRUE)]))
  cat(sprintf("  Ab Woche 2:      %.0f m (Median)\n\n",
              nacht_tab[tage_seit>7,  median(weg_gesamt_m, na.rm=TRUE)]))

  cat("── Distanz zum Auswilderungsgehege ──────\n")
  cat(sprintf("  Median (alle):   %.0f m\n",  median(nacht_tab$dist_release_m, na.rm=TRUE)))
  cat(sprintf("  Bereich:         %.0f – %.0f m\n",
              min(nacht_tab$dist_release_m, na.rm=TRUE),
              max(nacht_tab$dist_release_m, na.rm=TRUE)))
  cat(sprintf("  Woche 1:         %.0f m | Ab Woche 2: %.0f m\n\n",
              nacht_tab[tage_seit<=7, median(dist_release_m, na.rm=TRUE)],
              nacht_tab[tage_seit>7,  median(dist_release_m, na.rm=TRUE)]))

  cat("── Tortuositaet (Straightness Index) ────\n")
  cat(sprintf("  Median:          %.3f\n",  median(nacht_tab$tortuositaet, na.rm=TRUE)))
  cat(sprintf("  Bereich:         %.3f – %.3f\n",
              min(nacht_tab$tortuositaet, na.rm=TRUE),
              max(nacht_tab$tortuositaet, na.rm=TRUE)))
  cat("  (Werte < 0.10 = typisch kurvenreiches Foragierverhalten)\n\n")

  cat("── Aktivitaetsradius (95. Pz.) ──────────\n")
  cat(sprintf("  Median:          %.0f m  (%.0f – %.0f m)\n\n",
              median(nacht_tab$radius_95_m, na.rm=TRUE),
              min(nacht_tab$radius_95_m, na.rm=TRUE),
              max(nacht_tab$radius_95_m, na.rm=TRUE)))

  # Individuelle Trends
  tr <- nacht_tab[,.(
    tw = { x<-tage_seit;y<-weg_gesamt_m;v<-complete.cases(x,y)
           if(sum(v)>=3) coef(lm(y[v]~x[v]))[2] else NA_real_ },
    td = { x<-tage_seit;y<-dist_release_m;v<-complete.cases(x,y)
           if(sum(v)>=3) coef(lm(y[v]~x[v]))[2] else NA_real_ }
  ),by=igel]
  cat("── Individuelle Weglaenge-Trends ────────\n")
  cat(sprintf("  Zunehmend (> +5 m/Tag):  %d Igel\n", tr[!is.na(tw)&tw> 5,.N]))
  cat(sprintf("  Abnehmend (< -5 m/Tag):  %d Igel\n", tr[!is.na(tw) & tw < -5, .N]))
  cat(sprintf("  Stabil:                  %d Igel\n\n", tr[is.na(tw) | (tw >= -5 & tw <= 5), .N]))
  cat("── Individuelle Distanz-Trends ──────────\n")
  cat(sprintf("  Zunehmend (> +2 m/Tag):  %d Igel\n", tr[!is.na(td) & td >  2, .N]))
  cat(sprintf("  Abnehmend (< -2 m/Tag):  %d Igel\n", tr[!is.na(td) & td < -2, .N]))
  cat(sprintf("  Stabil:                  %d Igel\n",  tr[is.na(td)|(td>=-2&td<=2),.N]))
}

cat("\n====================================\n")
cat("Block 4d abgeschlossen!\n")
cat(sprintf("  Output: %s\n", normalizePath(output_ordner)))
cat("====================================\n")
