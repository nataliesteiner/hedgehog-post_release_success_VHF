# ==============================================================
# Block 4b — individual analysis & aKDE comparison
# ==============================================================
# Project:  Hedgehog VHF telemetry, Lower Saxony
# Author:   Natalie Steiner
#
# Goals:
#  1. KDE vs. aKDE (autocorrelated KDE, ctmm package) — comparison
#  2. Time resolution: raw vs. 5-min vs. 20-min bins (optional)
#  3. Individual analysis: cumulative HR, night radius, distance
#
# Requirement: Block4_Kernel_HomeRange.R must run first
# Input:  output/Block4_Kernel/kernel_ergebnisse.rds
# Output: output/Block4b_Einzeltier/
# ==============================================================

# ── EINSTELLUNGEN (hier anpassen) ──────────────────────────────
projekt_root        <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
release_punkt       <- data.frame(x = 514743.7, y = 5805363.8, label = "Auswilderungsgehege")
ziel_epsg           <- 25832

kde_grid            <- 100    # KDE-Gitterauflösung (kleiner = schneller; Standard: 100–150)
min_fixes_kde       <- 10     # Mindest-Fixes für KDE
min_fixes_akde      <- 30     # Mindest-Fixes für aKDE
zeitaufloesung_an   <- FALSE  # TRUE = Raw/5-min/20-min Vergleich aktivieren (langsam!)
kum_max_schritte    <- 15     # Max. Datenpunkte für kumulative HR-Kurve

# OSM-Hintergrundkarte (TRUE = schöner, braucht Internet & maptiles-Paket)
# Bei Problemen auf FALSE setzen → läuft stabiler und schneller
osm_tiles_an        <- TRUE
osm_zoom            <- 15     # Zoom-Stufe (14 = schneller, 16 = mehr Details)

# ── 0. Pakete ──────────────────────────────────────────────────
pakete <- c("sf", "adehabitatHR", "sp", "data.table", "ggplot2", "ggspatial",
            "patchwork", "scales", "viridis", "lubridate",
            "officer", "flextable", "openxlsx")
fehlend <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(fehlend) > 0) install.packages(fehlend)
invisible(lapply(c("sf","adehabitatHR","sp","data.table","ggplot2","ggspatial",
                   "patchwork","scales","viridis","lubridate"),
                 library, character.only = TRUE))
for (p in c("ctmm","officer","flextable","openxlsx"))
  if (requireNamespace(p, quietly=TRUE)) library(p, character.only=TRUE)

ctmm_verfuegbar <- requireNamespace("ctmm", quietly = TRUE)

# maptiles installieren falls nötig (löst terra/rosm-Konflikt)
if (osm_tiles_an && !requireNamespace("maptiles", quietly = TRUE)) {
  message("Installiere maptiles...")
  install.packages("maptiles", quiet = TRUE)
}
if (osm_tiles_an && requireNamespace("maptiles", quietly = TRUE))
  library(maptiles)

osm_tiles_verfuegbar <- osm_tiles_an && requireNamespace("maptiles", quietly = TRUE)
if (osm_tiles_an && !osm_tiles_verfuegbar)
  message("maptiles nicht verfuegbar → Karten ohne Hintergrund.")

# suncalc entfernen falls aus vorherigem Block geladen
try(detach("package:suncalc", unload = TRUE, character.only = TRUE), silent = TRUE)
try(unloadNamespace("suncalc"), silent = TRUE)

# ctmm überschreibt ggplot2::annotate — das hier stellt die richtige Version wieder her
annotate <- ggplot2::annotate
cat(sprintf("Pakete geladen | ctmm: %s | OSM-Tiles: %s\n\n",
            ctmm_verfuegbar, osm_tiles_verfuegbar))

# ── 1. Pfade ───────────────────────────────────────────────────
rds_pfad   <- file.path(projekt_root, "output", "Block4_Kernel", "kernel_ergebnisse.rds")
out_ordner <- file.path(projekt_root, "output", "Block4b_Einzeltier")
dir.create(out_ordner, showWarnings = FALSE, recursive = TRUE)

# ── 2. Daten laden ─────────────────────────────────────────────
if (!file.exists(rds_pfad))
  stop("Block4-RDS nicht gefunden:\n  ", rds_pfad,
       "\nBitte zuerst Block4_Kernel_HomeRange.R ausführen!")

block4     <- readRDS(rds_pfad)
alle_fixes <- block4$alle_fixes
kde_dt     <- block4$kde_gesamt
setorder(alle_fixes, igel, datetime)

cat(sprintf("Geladen: %d Fixes, %d Igel\n\n", nrow(alle_fixes), alle_fixes[, uniqueN(igel)]))

# ── 3. Hilfsfunktionen ─────────────────────────────────────────

# Standard-KDE — kernel.area() gibt Hektar zurück (unout="ha" ist Standard in adehabitatHR)
berechne_kde <- function(fixes_dt, h = "href", grid = kde_grid, min_fixes = min_fixes_kde) {
  if (is.null(fixes_dt) || nrow(fixes_dt) < min_fixes) return(NULL)
  tryCatch({
    sp_obj <- SpatialPoints(
      coords      = as.matrix(fixes_dt[, .(x, y)]),
      proj4string = CRS(SRS_string = paste0("EPSG:", ziel_epsg))
    )
    kde <- kernelUD(sp_obj, h = h, grid = grid)
    list(
      a95ha   = kernel.area(kde, percent = 95)[1],
      a50ha   = kernel.area(kde, percent = 50)[1],
      poly_95 = tryCatch(st_as_sf(getverticeshr(kde, 95)), error = function(e) NULL),
      poly_50 = tryCatch(st_as_sf(getverticeshr(kde, 50)), error = function(e) NULL)
    )
  }, error = function(e) NULL)
}

# UTM 32N → WGS84 für ctmm
utm_zu_wgs84 <- function(x_vec, y_vec, epsg_in = 25832) {
  pts <- st_as_sf(data.frame(x = x_vec, y = y_vec), coords = c("x","y"), crs = epsg_in)
  co  <- st_coordinates(st_transform(pts, 4326))
  data.frame(lon = co[,1], lat = co[,2])
}

# aKDE-Fläche aus AKDE-Objekt (m² → ha)
get_akde_area <- function(UD, level_ud = 0.95) {
  tryCatch({
    # units=FALSE erzwingt SI-Einheiten (m²) — ctmm 1.x skaliert sonst auto auf km²
    # was nach Division durch 10000 zu 0.00 führt
    s        <- summary(UD, level.UD = level_ud, units = FALSE)
    ci       <- s$CI
    est_col  <- grep("^est$|^50%$", colnames(ci), value = TRUE)[1]
    if (is.na(est_col)) est_col <- colnames(ci)[ceiling(ncol(ci)/2)]
    area_row <- grep("area|Area", rownames(ci), ignore.case = TRUE, value = TRUE)
    if (length(area_row) == 0) area_row <- rownames(ci)[1]
    val <- as.numeric(ci[area_row[1], est_col])
    if (!is.na(val) && val > 0) val / 10000 else NA_real_  # m² → ha
  }, error = function(e) NA_real_)
}

# aKDE-Polygon extrahieren
# BUG-FIX: ctmm:::SpatialPolygonsDataFrame.UD / .AKDE sind interne (nicht-exportierte)
# Funktionen und brechen mit neuen ctmm-Versionen. Neue Strategie:
#   Methode 1: exportierter S3-Generic as() → SpatialPolygonsDataFrame (ctmm >= 1.2)
#   Methode 2: interne .UD-Funktion (Fallback für ältere ctmm-Versionen)
#   Methode 3: interne .AKDE-Funktion (letzter Fallback)
extract_akde_poly <- function(UD, level_ud = 0.95, crs_out = 25832) {

  # Methode 1: exportierter Generic as(..., "SpatialPolygonsDataFrame")
  result <- tryCatch({
    sp_obj <- as(UD, "SpatialPolygonsDataFrame")
    sp_est <- sp_obj[grepl("est|50%", sp_obj$name, ignore.case = TRUE), ]
    if (nrow(sp_est) == 0) sp_est <- sp_obj[ceiling(nrow(sp_obj) / 2), ]
    st_transform(st_as_sf(sp_est), crs = crs_out)
  }, error = function(e) NULL)
  if (!is.null(result) && nrow(result) > 0) return(result)

  # Methode 2: interne UD-Funktion (ältere ctmm-Versionen)
  result <- tryCatch({
    sp_obj <- ctmm:::SpatialPolygonsDataFrame.UD(
      UD, level.UD = level_ud, level = 0.95
    )
    sp_est <- sp_obj[grepl("est|50%", sp_obj$name, ignore.case = TRUE), ]
    if (nrow(sp_est) == 0) sp_est <- sp_obj[ceiling(nrow(sp_obj) / 2), ]
    st_transform(st_as_sf(sp_est), crs = crs_out)
  }, error = function(e) NULL)
  if (!is.null(result) && nrow(result) > 0) return(result)

  # Methode 3: AKDE-spezifische interne Funktion (letzter Fallback)
  tryCatch({
    sp_obj <- ctmm:::SpatialPolygonsDataFrame.AKDE(UD, level.UD = level_ud, level = 0.95)
    sp_est <- sp_obj[grepl("est|50%", sp_obj$name, ignore.case = TRUE), ]
    if (nrow(sp_est) == 0) sp_est <- sp_obj[2, ]
    st_transform(st_as_sf(sp_est), crs = crs_out)
  }, error = function(e) NULL)
}

# Zeitliches Binning
bin_fixes_temporal <- function(fixes_dt, bin_min = 5) {
  bin_sek <- as.integer(bin_min * 60L)
  dt_c    <- copy(fixes_dt)
  dt_c[, time_bin := as.POSIXct(floor(as.numeric(datetime) / bin_sek) * bin_sek,
                                  origin = "1970-01-01", tz = "Europe/Berlin")]
  meta_cols <- intersect(names(fixes_dt),
    c("alter","sex","saison_auswild","diag_gruppe","release_date","time_reha","epsg"))
  binned <- dt_c[, .(x=mean(x,na.rm=TRUE), y=mean(y,na.rm=TRUE),
                      datetime=min(datetime), datum=first(datum),
                      tage_seit=first(tage_seit), Night=first(Night)),
                  by = .(igel, time_bin)]
  binned[, time_bin := NULL]
  if (length(meta_cols) > 0) {
    meta_extra <- unique(fixes_dt[, c("igel", meta_cols), with = FALSE])
    binned <- merge(binned, meta_extra, by = "igel", all.x = TRUE)
  }
  binned
}

# ── 4. ZEITAUFLÖSUNGS-VERGLEICH (optional) ─────────────────────
if (zeitaufloesung_an) {
  cat("== 4. Zeitauflösungs-Vergleich ==\n")
  fixes_5min  <- bin_fixes_temporal(alle_fixes, 5)
  fixes_20min <- bin_fixes_temporal(alle_fixes, 20)

  kde_fuer_aufloesung <- function(fixes_dt, label) {
    rbindlist(lapply(fixes_dt[, unique(igel)], function(ig) {
      fi  <- fixes_dt[igel == ig]
      res <- berechne_kde(fi)
      if (is.null(res)) return(NULL)
      data.table(igel=ig, aufloesung=label, n_fixes=nrow(fi),
                 a95ha=res$a95ha, a50ha=res$a50ha)
    }), fill = TRUE)
  }

  aufloesung_dt <- rbindlist(list(
    kde_fuer_aufloesung(alle_fixes,  "Raw"),
    kde_fuer_aufloesung(fixes_5min,  "5-min"),
    kde_fuer_aufloesung(fixes_20min, "20-min")
  ))
  aufloesung_dt[, aufloesung := factor(aufloesung, c("Raw","5-min","20-min"))]

  za_wide <- dcast(aufloesung_dt[aufloesung != "Raw"], igel ~ aufloesung, value.var = "a95ha")
  setnames(za_wide, c("5-min","20-min"), c("a95_5min","a95_20min"))
  za_wide <- za_wide[!is.na(a95_5min) & !is.na(a95_20min)]
  za_wide[, pct_diff := round(100 * (a95_20min - a95_5min) / a95_5min, 1)]

  p_za <- ggplot(za_wide, aes(x=a95_5min, y=a95_20min, label=igel)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
    geom_point(aes(color=igel), size=3, show.legend=FALSE) +
    geom_text(size=2.8, vjust=-0.6) +
    scale_color_viridis_d(option="turbo") +
    labs(title="Zeitauflösung: 5-min vs. 20-min Bins",
         x="95%-KDE bei 5-min (ha)", y="95%-KDE bei 20-min (ha)") +
    theme_bw(base_size=12)

  ggsave(file.path(out_ordner, "vergleich_zeitaufloesung.png"),
         p_za, width=10, height=7, dpi=150)
  cat("  vergleich_zeitaufloesung.png gespeichert\n\n")
} else {
  cat("== 4. Zeitauflösungs-Vergleich übersprungen (zeitaufloesung_an = FALSE) ==\n\n")
  aufloesung_dt <- data.table()
}

# ── 5. aKDE MIT ctmm ─────────────────────────────────────────
cat("== 5. aKDE (ctmm) ==\n")
akde_ergebnisse <- list()

if (!ctmm_verfuegbar) {
  cat("  ctmm nicht installiert — aKDE übersprungen\n\n")
} else {
  for (igel_i in alle_fixes[, unique(igel)]) {
    fixes_i <- alle_fixes[igel == igel_i][order(datetime)]
    n_i     <- nrow(fixes_i)
    cat(sprintf("  %s (%d Fixes) ... ", igel_i, n_i))

    if (n_i < min_fixes_akde) { cat("zu wenige Fixes\n"); next }

    tryCatch({
      wgs84  <- utm_zu_wgs84(fixes_i$x, fixes_i$y)
      tel_df <- data.frame(
        `individual.local.identifier` = igel_i,
        timestamp       = as.POSIXct(fixes_i$datetime, tz = "Europe/Berlin"),
        `location.long` = wgs84$lon,
        `location.lat`  = wgs84$lat,
        check.names = FALSE
      )
      tel   <- as.telemetry(tel_df, timezone = "Europe/Berlin")
      GUESS <- ctmm.guess(tel, interactive = FALSE)

      # ctmm.fit statt ctmm.select — deutlich schneller, kein hängen
      FIT   <- tryCatch(ctmm.fit(tel, CTMM = GUESS, trace = 0),
                        error = function(e) GUESS)

      modell_name <- tryCatch(summary(FIT)$name, error = function(e) "Unbekannt")
      UD  <- akde(tel, FIT)
      a95 <- get_akde_area(UD, 0.95)
      a50 <- get_akde_area(UD, 0.50)
      ess <- tryCatch(summary(UD)$DOF["area"], error = function(e) NA_real_)

      cat(sprintf("95%%=%.2f ha | ESS=%.1f [%s]\n", a95, ess, modell_name))

      akde_ergebnisse[[igel_i]] <- list(
        igel=igel_i, modell=modell_name, n_fixes=n_i,
        a95ha=a95, a50ha=a50, ess_area=ess,
        poly_95=extract_akde_poly(UD, 0.95),
        poly_50=extract_akde_poly(UD, 0.50),
        UD_objekt=UD
      )
    }, error = function(e) cat(sprintf("Fehler: %s\n", conditionMessage(e))))
  }
  cat(sprintf("\naKDE berechnet: %d von %d Igeln\n\n",
              length(akde_ergebnisse), alle_fixes[, uniqueN(igel)]))
}

# ── 6. aKDE-ZUSAMMENFASSUNG ───────────────────────────────────
# Standard-KDE wird nicht mehr verwendet (Primaermetrik: aKDE).
# Dieser Abschnitt fasst die aKDE-Ergebnisse zusammen.
cat("== 6. aKDE-Zusammenfassung ==\n")

if (length(akde_ergebnisse) > 0) {
  akde_summary <- rbindlist(lapply(akde_ergebnisse, function(x)
    data.table(igel      = x$igel,
               modell    = x$modell,
               akde_95ha = x$a95ha,
               akde_50ha = x$a50ha,
               ess_area  = x$ess_area,
               n_fixes   = x$n_fixes)))

  # Mit Metadaten aus kde_dt anreichern
  vergleich_dt <- merge(
    kde_dt[, .(igel, n_naechte)],
    akde_summary,
    by = "igel", all.x = TRUE
  )

  cat(sprintf("  aKDE erfolgreich: %d Tiere\n", sum(!is.na(akde_summary$akde_95ha))))
  cat(sprintf("  Median aKDE 95%%: %.2f ha (Range: %.2f - %.2f ha)\n",
              median(akde_summary$akde_95ha, na.rm=TRUE),
              min(akde_summary$akde_95ha, na.rm=TRUE),
              max(akde_summary$akde_95ha, na.rm=TRUE)))
  cat(sprintf("  Median ESS: %.0f\n\n", median(akde_summary$ess_area, na.rm=TRUE)))
} else {
  akde_summary <- data.table()
  vergleich_dt <- kde_dt[, .(igel, n_naechte)]
  vergleich_dt[, c("modell","akde_95ha","akde_50ha","ess_area","n_fixes") := NA]
  cat("  Kein aKDE berechnet.\n\n")
}

# ── 7. EINZELTIER-ANALYSE ─────────────────────────────────────
cat("== 7. Einzeltier-Analyse ==\n")

einzeltier_results <- list()

# maptiles arbeitet direkt mit terra — kein Vorab-Test nötig

for (igel_i in alle_fixes[, unique(igel)]) {
  fixes_i     <- alle_fixes[igel == igel_i][order(datetime)]
  n_i         <- nrow(fixes_i)
  naechte_srt <- sort(fixes_i[, unique(datum)])
  n_naechte_i <- length(naechte_srt)
  cat(sprintf("  %s (%d Fixes, %d Nächte)\n", igel_i, n_i, n_naechte_i))

  kde_i  <- berechne_kde(fixes_i)
  akde_i <- akde_ergebnisse[[igel_i]]
  fixes_i[, nacht_nr := as.integer(factor(datum, levels=naechte_srt))]

  # 7a: Kumulative HR — beschleunigt: max kum_max_schritte gleichmäßig verteilte Punkte
  idx_kum      <- unique(round(seq(1, n_naechte_i, length.out=min(n_naechte_i, kum_max_schritte))))
  kumulativ_dt <- rbindlist(lapply(idx_kum, function(k) {
    fi_k <- fixes_i[datum %in% naechte_srt[1:k]]
    if (nrow(fi_k) < 5) return(NULL)
    res_k <- berechne_kde(fi_k)
    if (is.null(res_k)) return(NULL)
    data.table(nacht_nr=k, datum=naechte_srt[k],
               tage_seit=fi_k[datum==naechte_srt[k], mean(tage_seit)],
               n_fixes=nrow(fi_k), a95ha=res_k$a95ha, a50ha=res_k$a50ha)
  }), fill=TRUE)

  # 7b: Nacht-Metriken
  nacht_metriken_dt <- rbindlist(lapply(naechte_srt, function(d) {
    fi_n <- fixes_i[datum == d]
    if (nrow(fi_n) < 3) return(NULL)
    cx_n <- mean(fi_n$x); cy_n <- mean(fi_n$y)
    dz   <- sqrt((fi_n$x - cx_n)^2 + (fi_n$y - cy_n)^2)
    data.table(datum=d, tage_seit=fi_n[1,tage_seit], nacht_nr=fi_n[1,nacht_nr],
               cx=cx_n, cy=cy_n, n_fixes_nacht=nrow(fi_n),
               radius_95pct_m=quantile(dz,0.95), radius_mean_m=mean(dz),
               dist_release_m=sqrt((cx_n-release_punkt$x)^2+(cy_n-release_punkt$y)^2))
  }), fill=TRUE)

  # Bounding box
  x_pad <- max((max(fixes_i$x)-min(fixes_i$x))*0.2, 150)
  y_pad <- max((max(fixes_i$y)-min(fixes_i$y))*0.2, 150)
  xlims <- c(min(fixes_i$x,release_punkt$x)-x_pad, max(fixes_i$x,release_punkt$x)+x_pad)
  ylims <- c(min(fixes_i$y,release_punkt$y)-y_pad, max(fixes_i$y,release_punkt$y)+y_pad)

  # Zoom automatisch an Ausdehnung anpassen (vermeidet zu viele Kacheln bei kleinen Gebieten)
  # BUG-FIX: dplyr::case_when() entfernt — dplyr ist nicht in der pakete-Liste.
  # Ersetzt durch nested ifelse() in Base R.
  ausdehnung_m <- max(diff(xlims), diff(ylims))
  zoom_auto <- ifelse(ausdehnung_m > 5000, 14L,
               ifelse(ausdehnung_m > 2000, 15L,
               ifelse(ausdehnung_m > 500,  16L, 17L)))

  # sf-Objekte für ggplot
  fixes_sf   <- st_as_sf(fixes_i, coords = c("x", "y"), crs = 25832)
  release_sf <- st_as_sf(release_punkt, coords = c("x", "y"), crs = 25832)

  # OSM-Hintergrund via maptiles (kein rosm/terra-Konflikt)
  osm_layer <- if (osm_tiles_verfuegbar) {
    tryCatch({
      aoi <- st_bbox(c(xmin=xlims[1], xmax=xlims[2], ymin=ylims[1], ymax=ylims[2]),
                     crs = st_crs(25832)) |> st_as_sfc() |> st_transform(4326)
      tiles <- maptiles::get_tiles(aoi, provider = "OpenStreetMap",
                                   zoom = zoom_auto, crop = TRUE)
      list(layer_spatial(tiles, alpha = 0.65))
    }, error = function(e) {
      message(sprintf("  [%s] Tiles nicht geladen: %s", igel_i, conditionMessage(e)))
      list()
    })
  } else {
    list()
  }
  # Plot 1: Karte
  p_karte <- ggplot() +
    osm_layer +
    geom_sf(data=fixes_sf, aes(color=nacht_nr),
            size=1.8, alpha=0.7, shape=16) +
    scale_color_viridis_c(name="Nacht\n(Nr.)", option="plasma", end=0.9) +
    {if (!is.null(kde_i) && !is.null(kde_i$poly_95))
      geom_sf(data=kde_i$poly_95, inherit.aes=FALSE,
              fill=NA, color="#2166ac", linewidth=0.9, linetype="dashed")} +
    {if (!is.null(kde_i) && !is.null(kde_i$poly_50))
      geom_sf(data=kde_i$poly_50, inherit.aes=FALSE,
              fill="#2166ac", alpha=0.18, color="#2166ac", linewidth=1.1)} +
    {if (!is.null(akde_i) && !is.null(akde_i$poly_95))
      geom_sf(data=akde_i$poly_95, inherit.aes=FALSE,
              fill=NA, color="#d6604d", linewidth=0.9, linetype="dotted")} +
    geom_sf(data=release_sf, shape=23, size=5,
            fill="yellow", color="black", stroke=1.5, inherit.aes=FALSE) +
    ggplot2::annotate("text", x=release_punkt$x, y=release_punkt$y,
             label="Auswild.", size=2.8, color="grey20",
             vjust=2.5) +
    coord_sf(xlim=xlims, ylim=ylims, crs=st_crs(25832), expand=FALSE) +
    labs(title=paste0(igel_i, " — Nacht-Aktionsareal"),
         subtitle=sprintf("KDE: 95%%=%.2f ha | 50%%=%.2f ha | %d Fixes | %d Nächte",
           if(!is.null(kde_i)) kde_i$a95ha else NA,
           if(!is.null(kde_i)) kde_i$a50ha else NA, n_i, n_naechte_i),
         x="UTM32N Easting (m)", y="UTM32N Northing (m)") +
    theme_bw(base_size=11) +
    theme(plot.title=element_text(face="bold"),
          axis.title=element_blank())

  # Plot 2: Kumulative HR
  p_kumulativ <- if (nrow(kumulativ_dt) >= 2) {
    ggplot(kumulativ_dt, aes(x=nacht_nr)) +
      geom_ribbon(aes(ymin=a50ha, ymax=a95ha), fill="#2166ac", alpha=0.12) +
      geom_line(aes(y=a95ha), color="#2166ac", linewidth=1.0) +
      geom_point(aes(y=a95ha), color="#2166ac", size=2.5) +
      geom_line(aes(y=a50ha), color="#d6604d", linewidth=0.8, linetype="dashed") +
      geom_point(aes(y=a50ha), color="#d6604d", size=2.0) +
      labs(title="Kumulative Home Range",
           subtitle="Blau=95% | Rot=50%",
           x="Kumulierte Nächte", y="Areal (ha)") +
      theme_bw(base_size=10) + theme(plot.title=element_text(face="bold", size=10.5))
  } else {
    ggplot() + labs(title="Zu wenige Daten") + theme_bw(base_size=10)
  }

  # Plot 3: Nacht-Radius
  p_radius <- if (nrow(nacht_metriken_dt) >= 2) {
    ggplot(nacht_metriken_dt, aes(x=tage_seit)) +
      geom_col(aes(y=radius_95pct_m), fill="#74c476", alpha=0.75) +
      geom_line(aes(y=radius_mean_m), color="#238b45", linewidth=0.9) +
      geom_point(aes(y=radius_mean_m), color="#238b45", size=2) +
      scale_y_continuous(labels=scales::label_number(suffix=" m")) +
      labs(title="Nacht-Aktivitätsradius", subtitle="Balken=95% | Linie=Mittelwert",
           x="Tage seit Auswilderung", y="Distanz vom Nachtzentroid") +
      theme_bw(base_size=10) + theme(plot.title=element_text(face="bold", size=10.5))
  } else {
    ggplot() + labs(title="Zu wenige Daten") + theme_bw(base_size=10)
  }

  # Plot 4: Distanz zum Auswilderungsgehege
  p_distanz <- if (nrow(nacht_metriken_dt) >= 2) {
    p_d <- ggplot(nacht_metriken_dt, aes(x=tage_seit, y=dist_release_m)) +
      geom_line(color="#9e9ac8", linewidth=0.8, alpha=0.7) +
      geom_point(color="#756bb1", size=2.5)
    if (nrow(nacht_metriken_dt) >= 5)
      p_d <- p_d + geom_smooth(method="loess", span=0.8, se=TRUE,
                                color="#3f007d", fill="#bcbddc", alpha=0.25)
    p_d +
      scale_y_continuous(labels=scales::label_number(suffix=" m")) +
      labs(title="Distanz zum Auswilderungsgehege",
           x="Tage seit Auswilderung", y="Distanz (m)") +
      theme_bw(base_size=10) + theme(plot.title=element_text(face="bold", size=10.5))
  } else {
    ggplot() + labs(title="Zu wenige Daten") + theme_bw(base_size=10)
  }

  p_kombi <- (p_karte | (p_kumulativ / p_radius / p_distanz)) +
    plot_layout(widths=c(1.5,1)) +
    plot_annotation(
      title    = paste0("Einzeltier-Analyse: ", igel_i),
      subtitle = sprintf("KDE: 95%%=%.2f ha, 50%%=%.2f ha | %d Nächte | %d Fixes",
        if(!is.null(kde_i)) kde_i$a95ha else NA,
        if(!is.null(kde_i)) kde_i$a50ha else NA, n_naechte_i, n_i),
      theme = theme(plot.title=element_text(face="bold", size=13),
                    plot.subtitle=element_text(size=9, color="grey30"))
    )

  png_name <- paste0("einzeltier_", gsub("[^a-zA-Z0-9]","_",igel_i), ".png")
  ggsave(file.path(out_ordner, png_name), p_kombi, width=18, height=9, dpi=150)
  cat(sprintf("    %s gespeichert\n", png_name))

  einzeltier_results[[igel_i]] <- list(
    n_fixes=n_i, n_naechte=n_naechte_i,
    kde_95ha    = if(!is.null(kde_i)) kde_i$a95ha  else NA_real_,
    kde_50ha    = if(!is.null(kde_i)) kde_i$a50ha  else NA_real_,
    akde_95ha   = if(!is.null(akde_i)) akde_i$a95ha  else NA_real_,
    akde_modell = if(!is.null(akde_i)) akde_i$modell else NA_character_,
    kumulativ_dt=kumulativ_dt, nacht_metriken=nacht_metriken_dt
  )
}

cat(sprintf("\nEinzeltier-Analyse: %d Igel\n\n", length(einzeltier_results)))

# ── 8. ÜBERSICHTS-PLOTS ───────────────────────────────────────
et_summary <- rbindlist(lapply(names(einzeltier_results), function(ig) {
  x  <- einzeltier_results[[ig]]
  nm <- x$nacht_metriken
  data.table(
    igel=ig, n_fixes=x$n_fixes, n_naechte=x$n_naechte,
    kde_95ha=x$kde_95ha, kde_50ha=x$kde_50ha,
    akde_95ha=x$akde_95ha, akde_modell=x$akde_modell,
    median_dist_m   = if(!is.null(nm)&&nrow(nm)>0) round(median(nm$dist_release_m, na.rm=TRUE)) else NA_real_,
    final_dist_m    = if(!is.null(nm)&&nrow(nm)>0) round(nm[.N, dist_release_m]) else NA_real_,
    median_radius_m = if(!is.null(nm)&&nrow(nm)>0) round(median(nm$radius_95pct_m, na.rm=TRUE)) else NA_real_
  )
}), fill=TRUE)

if ("alter" %in% names(kde_dt))
  et_summary <- merge(et_summary,
    kde_dt[, .(igel, alter, sex, saison_auswild, time_reha)], by="igel", all.x=TRUE)

# Distanz-Übersicht
if (et_summary[!is.na(median_dist_m), .N] >= 3) {
  et_ord <- et_summary[!is.na(median_dist_m)][order(median_dist_m)]
  et_ord[, igel_f := factor(igel, levels=igel)]

  p_et_dist <- ggplot(et_ord, aes(x=igel_f, y=median_dist_m/1000)) +
    geom_col(aes(fill=if("alter"%in%names(et_ord)) alter else "all"),
             width=0.7, alpha=0.85, show.legend="alter"%in%names(et_ord)) +
    geom_point(aes(y=final_dist_m/1000), shape=18, size=3, color="black") +
    geom_hline(yintercept=median(et_ord$median_dist_m,na.rm=TRUE)/1000,
               linetype="dashed", color="grey40") +
    scale_fill_manual(values=c("Jungtier"="#fdae61","Adult"="#4dac26","all"="#6baed6"),
                      name="Altersklasse") +
    coord_flip() +
    labs(title="Distanz zum Auswilderungsgehege",
         subtitle="Balken=Median | Raute=letzte Nacht",
         x=NULL, y="Distanz (km)") +
    theme_bw(base_size=12) + theme(legend.position="bottom")

  ggsave(file.path(out_ordner, "einzeltier_distanz_uebersicht.png"),
         p_et_dist, width=9, height=7, dpi=150)
  cat("  einzeltier_distanz_uebersicht.png gespeichert\n")
}

# Kumulative HR — alle Tiere
kum_alle <- rbindlist(lapply(names(einzeltier_results), function(ig) {
  x <- einzeltier_results[[ig]]$kumulativ_dt
  if (!is.null(x) && nrow(x)>0) { x[, igel:=ig]; x }
}), fill=TRUE)

if (nrow(kum_alle) >= 5) {
  p_kum_all <- ggplot(kum_alle, aes(x=nacht_nr, y=a95ha)) +
    geom_line(color="#2166ac", linewidth=0.8) +
    geom_point(color="#2166ac", size=1.5) +
    geom_line(aes(y=a50ha), color="#d6604d", linewidth=0.6, linetype="dashed") +
    facet_wrap(~igel, scales="free_y", ncol=5) +
    labs(title="Kumulative Home Range (alle Igel)",
         subtitle="Blau=95% | Rot gestrichelt=50%",
         x="Nächte (kumuliert)", y="Areal (ha)") +
    theme_bw(base_size=9) +
    theme(strip.text=element_text(size=8,face="bold"), panel.spacing=unit(0.3,"lines"))

  ggsave(file.path(out_ordner, "einzeltier_kumulativ_alle.png"), p_kum_all,
         width=16, height=max(4, ceiling(length(alle_fixes[,unique(igel)])/5)*3.5), dpi=150)
  cat("  einzeltier_kumulativ_alle.png gespeichert\n\n")
}

# ── 9. aKDE PAARWEISE OVERLAP-ANALYSE ────────────────────────
# Für Publikation: Überlappung der aKDE-Streifgebiete gleichzeitig
# freigelassener Paare. Methode: Jaccard-Index + Overlap-Koeffizient
# über sf-Polygonschnitt der aKDE-95%-Konturen.
# HINWEIS: Block4_Kernel verwendet Standard-KDE-Polygone für diese Analyse.
# Hier wird korrekt aKDE verwendet — Ergebnisse für Paper bevorzugen.
cat("== 9. aKDE Paarweise Overlap ==\n")

akde_overlap_dt <- NULL

if (length(akde_ergebnisse) >= 2 && requireNamespace("sf", quietly = TRUE)) {

  # Metadaten laden um gleichzeitig freigelassene Paare zu finden
  meta_4b <- tryCatch({
    m <- as.data.table(readxl::read_excel(
      file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")))
    data.table(
      igel         = trimws(m$individual),
      release_date = as.Date(m$date_release)
    )
  }, error = function(e) NULL)

  sf_ueberlap <- function(p1, p2) {
    p1u <- tryCatch(st_union(p1), error = function(e) p1)
    p2u <- tryCatch(st_union(p2), error = function(e) p2)
    inter <- tryCatch(st_intersection(p1u, p2u), error = function(e) NULL)
    a_i  <- if (is.null(inter) || length(inter) == 0) 0 else
              as.numeric(st_area(inter))
    a1 <- as.numeric(st_area(p1u))
    a2 <- as.numeric(st_area(p2u))
    list(
      jaccard = round(if ((a1+a2-a_i) > 0) a_i/(a1+a2-a_i) else 0, 3),
      oc      = round(if (min(a1,a2)  > 0) a_i/min(a1,a2)  else 0, 3),
      a_inter_ha = round(a_i / 10000, 3)
    )
  }

  if (!is.null(meta_4b)) {
    rd_counts <- meta_4b[!is.na(release_date), .N, by = release_date][N >= 2]

    if (nrow(rd_counts) > 0) {
      overlap_liste <- list()
      for (rd in rd_counts$release_date) {
        igel_rd <- meta_4b[release_date == rd, igel]
        combs   <- combn(igel_rd, 2, simplify = FALSE)
        for (cp in combs) {
          ia <- cp[1]; ib <- cp[2]
          pl <- paste0(ia, "_", ib)
          # Nur Paare für die aKDE vorhanden ist
          hat_akde <- ia %in% names(akde_ergebnisse) &&
                      ib %in% names(akde_ergebnisse) &&
                      !is.null(akde_ergebnisse[[ia]]$poly_95) &&
                      !is.null(akde_ergebnisse[[ib]]$poly_95)

          cat(sprintf("  Paar %s: aKDE verfuegbar: %s\n", pl, hat_akde))

          if (!hat_akde) {
            cat(sprintf("  [WARN] %s: kein aKDE-Polygon — Paar uebersprungen\n", pl))
            next
          }
          poly_a_95 <- akde_ergebnisse[[ia]]$poly_95
          poly_b_95 <- akde_ergebnisse[[ib]]$poly_95

          tryCatch({
            ov95 <- sf_ueberlap(akde_ergebnisse[[ia]]$poly_95,
                                 akde_ergebnisse[[ib]]$poly_95)
            ov50 <- if (!is.null(akde_ergebnisse[[ia]]$poly_50) &&
                         !is.null(akde_ergebnisse[[ib]]$poly_50))
              sf_ueberlap(akde_ergebnisse[[ia]]$poly_50,
                           akde_ergebnisse[[ib]]$poly_50)
            else list(jaccard=NA, oc=NA, a_inter_ha=NA)

            overlap_liste[[pl]] <- data.table(
              Paar              = pl,
              Igel_A            = ia,
              Igel_B            = ib,
              Release_Datum     = format(as.Date(rd, origin="1970-01-01"), "%d.%m.%Y"),
              aKDE_95_Jaccard   = ov95$jaccard,
              aKDE_95_OC        = ov95$oc,
              aKDE_95_Schnitt_ha= ov95$a_inter_ha,
              aKDE_50_Jaccard   = ov50$jaccard,
              aKDE_50_OC        = ov50$oc,
              Interpretation    = fcase(
                ov95$jaccard >= 0.5, "Stark ueberlappend",
                ov95$jaccard >= 0.2, "Maessig ueberlappend",
                default = "Kaum ueberlappend")
            )
            cat(sprintf("  %s: Jaccard_95=%.3f | OC_95=%.3f | Schnitt=%.3f ha → %s\n",
                        pl, ov95$jaccard, ov95$oc, ov95$a_inter_ha,
                        overlap_liste[[pl]]$Interpretation))

            # Karte Paar
            poly_95_both <- do.call(rbind, list(
              cbind(akde_ergebnisse[[ia]]$poly_95, igel=ia),
              cbind(akde_ergebnisse[[ib]]$poly_95, igel=ib)
            ))
            p_paar <- ggplot() +
              geom_sf(data=poly_95_both, aes(fill=igel, color=igel),
                      alpha=0.25, linewidth=1.1) +
              {if (!is.null(akde_ergebnisse[[ia]]$poly_50) &&
                    !is.null(akde_ergebnisse[[ib]]$poly_50))
                geom_sf(data=do.call(rbind, list(
                  cbind(akde_ergebnisse[[ia]]$poly_50, igel=ia),
                  cbind(akde_ergebnisse[[ib]]$poly_50, igel=ib)
                )), aes(fill=igel, color=igel), alpha=0.4, linewidth=0.7)} +
              scale_fill_manual(values=c("#2166ac","#d6604d"), name="Igel") +
              scale_color_manual(values=c("#2166ac","#d6604d"), guide="none") +
              labs(
                title    = sprintf("aKDE: %s vs. %s", ia, ib),
                subtitle = sprintf("95%%-Jaccard=%.3f | OC=%.3f | Schnitt=%.3f ha",
                                   ov95$jaccard, ov95$oc, ov95$a_inter_ha),
                caption  = "Aussen: 95%-aKDE | Innen: 50%-aKDE | Methode: ctmm"
              ) +
              theme_bw(base_size=11) +
              theme(plot.title    = element_text(face="bold"),
                    plot.subtitle = element_text(size=9, color="grey40"),
                    plot.caption  = element_text(size=8, color="grey50"))
            ggsave(file.path(out_ordner, paste0("akde_paar_", pl, ".png")),
                   p_paar, width=10, height=7, dpi=150)
            cat(sprintf("    akde_paar_%s.png gespeichert\n", pl))

          }, error = function(e)
            cat(sprintf("  ! Fehler bei %s: %s\n", pl, conditionMessage(e))))
        }
      }

      if (length(overlap_liste) > 0) {
        akde_overlap_dt <- rbindlist(overlap_liste, fill=TRUE)
        cat("\n  aKDE Paarweise Overlap Ergebnisse:\n")
        print(akde_overlap_dt[, .(Paar, aKDE_95_Jaccard, aKDE_95_OC,
                                   aKDE_95_Schnitt_ha, Interpretation)])
      }
    } else {
      cat("  Keine gleichzeitigen Auswilderungen (≥2 Igel selber Tag)\n")
    }
  }
} else {
  cat("  aKDE-Polygone nicht verfuegbar — Block4b mit ctmm neu ausfuehren\n")
}
cat("\n")

# ── 9b. POPULATIONSKARTEN auf ATKIS-Landcover (3 Varianten) ───────────────────
# A: ATKIS-Hintergrund + gepoolte Populations-UD (95% Linie / 50% Flaeche)
# B: ATKIS-Hintergrund + alle 50%-Cores in EINER neutralen Farbe (Dichte)
# C: Small multiples — ein Panel je Tier (95% Linie + 50% Flaeche)
# Hintergrund: ATKIS Basis-DLM (habitat_sachsenhagen_atkis.gpkg) in EPSG:25832,
# kein Internet/Satellit noetig. Plot in WGS84 (Grad-Achsen).
cat("== 9b. Populationskarten (ATKIS, 3 Varianten) ==\n")

if (length(akde_ergebnisse) >= 2) {

  suppressWarnings(suppressMessages({ library(sf); library(ggplot2) }))

  ## (1) aKDE-Polygone aller Tiere sammeln (UTM 25832) -------------------------
  igel_sorted <- names(akde_ergebnisse)
  igel_sorted <- igel_sorted[order(as.numeric(gsub("\\D", "", igel_sorted)))]

  pl95 <- list(); pl50 <- list()
  for (ig in igel_sorted) {
    p95 <- akde_ergebnisse[[ig]]$poly_95
    p50 <- akde_ergebnisse[[ig]]$poly_50
    if (!is.null(p95) && inherits(p95, "sf")) { p95$igel <- ig; pl95[[ig]] <- p95[, c("igel", "geometry")] }
    if (!is.null(p50) && inherits(p50, "sf")) { p50$igel <- ig; pl50[[ig]] <- p50[, c("igel", "geometry")] }
  }

  if (length(pl95) >= 2) {

    alle_95 <- do.call(rbind, pl95); if (is.na(st_crs(alle_95))) st_crs(alle_95) <- 25832
    alle_50 <- do.call(rbind, pl50); if (is.na(st_crs(alle_50))) st_crs(alle_50) <- 25832
    alle_95$igel <- factor(alle_95$igel, levels = igel_sorted)
    alle_50$igel <- factor(alle_50$igel, levels = igel_sorted)

    # Plot-Ausschnitt (UTM) mit 60 m Puffer um 95%-Polygone
    bb_utm <- st_bbox(alle_95)
    bb_utm["xmin"] <- bb_utm["xmin"] - 60; bb_utm["ymin"] <- bb_utm["ymin"] - 60
    bb_utm["xmax"] <- bb_utm["xmax"] + 60; bb_utm["ymax"] <- bb_utm["ymax"] + 60

    ## (2) ATKIS-Landcover laden, auf Englisch mappen, Farben -----------------
    atkis_pfad <- file.path(projekt_root, "data", "excel_files",
                            "habitat_sachsenhagen_atkis.gpkg")
    hl_labels <- c(Wald = "Forest", Gehoelz_Strauch = "Shrubland / hedgerow", Acker = "Arable land",
                   Wiese = "Grassland", Gewaesser = "Water body", Siedlung = "Built-up area")
    hl_cols <- c("Forest" = "#6E9E73", "Shrubland / hedgerow" = "#AECB9A", "Arable land" = "#E6DAB0",
                 "Grassland" = "#D2E8BC", "Water body" = "#A6CEE3", "Built-up area" = "#D8D2C6",
                 "Other" = "#E3E0D8")
    hab_w <- NULL
    if (file.exists(atkis_pfad)) {
      hab <- tryCatch(st_read(atkis_pfad, quiet = TRUE), error = function(e) NULL)
      if (!is.null(hab)) {
        if (is.na(st_crs(hab))) st_crs(hab) <- 25832
        cls <- unname(hl_labels[as.character(hab$habitat)]); cls[is.na(cls)] <- "Other"
        hab$cls <- factor(cls, levels = names(hl_cols))
        hab <- suppressWarnings(st_crop(hab, bb_utm))
        hab_w <- st_transform(hab, 4326)
      }
    } else {
      cat("  ATKIS-GPKG nicht gefunden — Hintergrund bleibt leer:\n  ", atkis_pfad, "\n")
    }

    ## (3) ATKIS-Wege (optional, duenne Linien) ------------------------------
    wege_w <- NULL
    wege_pfad <- file.path(projekt_root, "data", "excel_files",
                           "atkis_wege_sachsenhagen.gpkg")
    if (file.exists(wege_pfad)) {
      wege <- tryCatch(st_read(wege_pfad, quiet = TRUE), error = function(e) NULL)
      if (!is.null(wege)) {
        if (is.na(st_crs(wege))) st_crs(wege) <- 25832
        wege <- suppressWarnings(st_crop(wege, bb_utm))
        wege_w <- st_transform(wege, 4326)
      }
    }

    ## (4) aKDE-Polygone + Gehege nach WGS84, Plot-Limits --------------------
    alle_95_w <- st_transform(alle_95, 4326)
    alle_50_w <- st_transform(alle_50, 4326)
    enc_w     <- st_transform(st_as_sf(release_punkt, coords = c("x", "y"), crs = 25832), 4326)
    bb_w <- st_bbox(alle_95_w)
    dx <- (bb_w["xmax"] - bb_w["xmin"]) * 0.06; dy <- (bb_w["ymax"] - bb_w["ymin"]) * 0.06
    xlim_w <- c(bb_w["xmin"] - dx, bb_w["xmax"] + dx)
    ylim_w <- c(bb_w["ymin"] - dy, bb_w["ymax"] + dy)

    ## (5) Gepoolte Populations-UD (KDE ueber ALLE Fixes) -------------------
    pooled_95_w <- NULL; pooled_50_w <- NULL
    if (requireNamespace("adehabitatHR", quietly = TRUE) &&
        requireNamespace("sp", quietly = TRUE)) {
      pooled <- tryCatch({
        xy  <- as.data.frame(alle_fixes[, .(x, y)])
        spp <- sp::SpatialPoints(xy, proj4string = sp::CRS(SRS_string = "EPSG:25832"))
        ud  <- adehabitatHR::kernelUD(spp, h = "href", grid = 200, extent = 1.2)
        v95 <- sf::st_as_sf(adehabitatHR::getverticeshr(ud, 95))
        v50 <- sf::st_as_sf(adehabitatHR::getverticeshr(ud, 50))
        if (is.na(st_crs(v95))) st_crs(v95) <- 25832
        if (is.na(st_crs(v50))) st_crs(v50) <- 25832
        list(p95 = st_transform(v95, 4326), p50 = st_transform(v50, 4326))
      }, error = function(e) { message("  Pooled-KDE fehlgeschlagen: ", conditionMessage(e)); NULL })
      if (!is.null(pooled)) { pooled_95_w <- pooled$p95; pooled_50_w <- pooled$p50 }
    }
    if (is.null(pooled_95_w)) {                 # Fallback: Vereinigung der Einzel-UDs
      pooled_95_w <- st_union(alle_95_w); pooled_50_w <- st_union(alle_50_w)
      cat("  (Pooled-KDE nicht verfuegbar — Vereinigung der Einzel-UDs als Huelle)\n")
    }

    ## Gemeinsame Hintergrund-Layer + Theme ---------------------------------
    base_layers <- function() {
      ll <- list()
      if (!is.null(hab_w))  ll <- c(ll, list(geom_sf(data = hab_w,  aes(fill = cls), colour = NA, alpha = 0.90)))
      if (!is.null(wege_w)) ll <- c(ll, list(geom_sf(data = wege_w, colour = "grey35", linewidth = 0.25)))
      ll
    }
    common_theme <- theme_bw(base_size = 10) +
      theme(legend.position = "right",
            plot.title    = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 8.5, colour = "grey30"),
            plot.caption  = element_text(size = 6.5, colour = "grey50", hjust = 1),
            panel.grid    = element_line(colour = "grey90", linewidth = 0.2),
            axis.text     = element_text(size = 8))
    grad_x <- scale_x_continuous(labels = function(x) sprintf("%.3f°E", x))
    grad_y <- scale_y_continuous(labels = function(y) sprintf("%.3f°N", y))
    scale_atkis <- scale_fill_manual(name = "Land cover (ATKIS)", values = hl_cols, drop = TRUE)
    n_pop <- length(igel_sorted)

    ## ── MAP A: ATKIS + gepoolte Populations-UD ───────────────────────────
    p_A <- ggplot() + base_layers() +
      geom_sf(data = pooled_95_w, fill = NA, colour = "#7F0000",
              linewidth = 0.8, linetype = "longdash") +
      geom_sf(data = pooled_50_w, fill = "#B2182B", colour = "#7F0000",
              alpha = 0.45, linewidth = 0.5) +
      geom_sf(data = enc_w, shape = 23, size = 3.6, fill = "#FFD700",
              colour = "black", stroke = 1.1) +
      scale_atkis +
      coord_sf(xlim = xlim_w, ylim = ylim_w, expand = FALSE) + grad_x + grad_y +
      ggspatial::annotation_scale(location = "bl", width_hint = 0.18,
              bar_cols = c("black", "white"), text_cex = 0.8) +
      ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
              height = unit(0.9, "cm"), width = unit(0.6, "cm"),
              style = ggspatial::north_arrow_fancy_orienteering()) +
      labs(title = "Population space use on ATKIS land cover — Sachsenhagen",
           subtitle = sprintf("Pooled utilization distribution | dashed: 95%% UD | filled: 50%% core | N = %d", n_pop),
           x = "Longitude", y = "Latitude",
           caption = "Land cover: ATKIS Basis-DLM, © LGLN (2026), CC BY 4.0") +
      common_theme
    ggsave(file.path(out_ordner, "akde_map_A_atkis_pooled.png"), p_A,
           width = 20, height = 18, units = "cm", dpi = 300, bg = "white")
    ggsave(file.path(out_ordner, "akde_map_A_atkis_pooled.pdf"), p_A,
           width = 20, height = 18, units = "cm", bg = "white")
    cat("  Map A (ATKIS + pooled UD) gespeichert\n")

    ## ── MAP B: ATKIS + alle 50%-Cores in EINER Farbe (Dichte) ────────────
    p_B <- ggplot() + base_layers() +
      geom_sf(data = alle_50_w, fill = "#3B0F70", colour = "#3B0F70",
              alpha = 0.16, linewidth = 0.2) +
      geom_sf(data = enc_w, shape = 23, size = 3.6, fill = "#FFD700",
              colour = "black", stroke = 1.1) +
      scale_atkis +
      coord_sf(xlim = xlim_w, ylim = ylim_w, expand = FALSE) + grad_x + grad_y +
      ggspatial::annotation_scale(location = "bl", width_hint = 0.18,
              bar_cols = c("black", "white"), text_cex = 0.8) +
      ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
              height = unit(0.9, "cm"), width = unit(0.6, "cm"),
              style = ggspatial::north_arrow_fancy_orienteering()) +
      labs(title = "Overlap of individual core areas — Sachsenhagen",
           subtitle = sprintf("All 50%% aKDE cores in one colour (transparent) | darker = more overlap | N = %d", n_pop),
           x = "Longitude", y = "Latitude",
           caption = "Land cover: ATKIS Basis-DLM, © LGLN (2026), CC BY 4.0") +
      common_theme
    ggsave(file.path(out_ordner, "akde_map_B_atkis_cores.png"), p_B,
           width = 20, height = 18, units = "cm", dpi = 300, bg = "white")
    ggsave(file.path(out_ordner, "akde_map_B_atkis_cores.pdf"), p_B,
           width = 20, height = 18, units = "cm", bg = "white")
    cat("  Map B (ATKIS + 50%-Cores neutral) gespeichert\n")

    ## ── MAP C: Small multiples (ein Panel je Tier) ───────────────────────
    p_C <- ggplot() +
      { if (!is.null(hab_w)) geom_sf(data = hab_w, aes(fill = cls), colour = NA, alpha = 0.85) } +
      geom_sf(data = alle_95_w, fill = NA, colour = "#08519C", linewidth = 0.4) +
      geom_sf(data = alle_50_w, fill = "#08519C", colour = "#08519C",
              alpha = 0.40, linewidth = 0.2) +
      geom_sf(data = enc_w, shape = 23, size = 1.5, fill = "#FFD700",
              colour = "black", stroke = 0.5) +
      scale_atkis +
      facet_wrap(~ igel, ncol = 5,
                 labeller = as_labeller(function(x) sub("^Igel", "H", x))) +
      coord_sf(xlim = xlim_w, ylim = ylim_w, expand = FALSE) +
      labs(title = "Individual nocturnal home ranges — Sachsenhagen",
           subtitle = "95% UD (outline) and 50% core (filled) per individual on ATKIS land cover",
           x = NULL, y = NULL,
           caption = "Land cover: ATKIS Basis-DLM, © LGLN (2026), CC BY 4.0") +
      theme_bw(base_size = 8) +
      theme(legend.position = "bottom",
            axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_rect(fill = "grey92", colour = NA),
            strip.text = element_text(face = "bold", size = 8),
            plot.title = element_text(face = "bold", size = 11))
    ggsave(file.path(out_ordner, "akde_map_C_smallmultiples.png"), p_C,
           width = 24, height = 26, units = "cm", dpi = 300, bg = "white")
    ggsave(file.path(out_ordner, "akde_map_C_smallmultiples.pdf"), p_C,
           width = 24, height = 26, units = "cm", bg = "white")
    cat("  Map C (small multiples) gespeichert\n\n")

  } else {
    cat("  Zu wenig aKDE-Polygone fuer Populationskarten\n\n")
  }
} else {
  cat("  Keine aKDE-Ergebnisse — Populationskarten uebersprungen\n\n")
}

# ── 10. EXPORT ─────────────────────────────────────────────────
cat("== 10. Export ==\n")

saveRDS(list(
  einzeltier_results=einzeltier_results, akde_ergebnisse=akde_ergebnisse,
  vergleich_dt=vergleich_dt, aufloesung_dt=aufloesung_dt, et_summary=et_summary,
  akde_overlap=akde_overlap_dt   # aKDE Paaranalyse — für Paper bevorzugen
), file.path(out_ordner, "Block4b_Ergebnisse.rds"))
cat("  RDS gespeichert\n")

if (requireNamespace("openxlsx", quietly=TRUE)) {
  wb   <- createWorkbook()
  st_h <- createStyle(fontName="Arial", fontSize=11, fontColour="white",
                      fgFill="#2C5F8A", halign="center", textDecoration="bold")
  st_n <- createStyle(fontName="Arial", fontSize=10)
  st_2 <- createStyle(fontName="Arial", fontSize=10, numFmt="0.00", halign="center")
  st_0 <- createStyle(fontName="Arial", fontSize=10, numFmt="0",    halign="center")

  hdr <- function(wb, sh, cols, row=1) {
    writeData(wb, sh, as.data.frame(t(cols)), startRow=row, colNames=FALSE)
    addStyle(wb, sh, st_h, rows=row, cols=seq_along(cols), gridExpand=TRUE)
  }

  # Sheet 1: Einzeltier
  addWorksheet(wb, "Einzeltier", tabColour="#2C5F8A")
  writeData(wb, "Einzeltier", "Block 4b: Einzeltier-Übersicht", startRow=1)
  hdr(wb, "Einzeltier",
      c("Igel","N Fixes","N Nächte","KDE 95% (ha)","KDE 50% (ha)",
        "aKDE 95% (ha)","aKDE-Modell","Med. Dist. (m)","Letzte Dist. (m)","Med. Radius (m)"),
      row=3)
  xl1 <- et_summary[order(igel), .(igel, n_fixes, n_naechte,
    round(kde_95ha,2), round(kde_50ha,2), round(akde_95ha,2),
    akde_modell, median_dist_m, final_dist_m, median_radius_m)]
  writeData(wb, "Einzeltier", xl1, startRow=4, colNames=FALSE)
  addStyle(wb, "Einzeltier", st_n, rows=4:(3+nrow(xl1)), cols=1,    gridExpand=TRUE)
  addStyle(wb, "Einzeltier", st_0, rows=4:(3+nrow(xl1)), cols=2:3,  gridExpand=TRUE)
  addStyle(wb, "Einzeltier", st_2, rows=4:(3+nrow(xl1)), cols=4:6,  gridExpand=TRUE)
  addStyle(wb, "Einzeltier", st_0, rows=4:(3+nrow(xl1)), cols=8:10, gridExpand=TRUE)

  # Sheet 2: aKDE-Uebersicht (Primaermetrik)
  addWorksheet(wb, "aKDE_Uebersicht", tabColour="#4dac26")
  writeData(wb, "aKDE_Uebersicht", "aKDE — Homerange-Groessen (ctmm, Primaermetrik)", startRow=1)
  hdr(wb, "aKDE_Uebersicht",
      c("Igel","N Fixes","N Naechte","aKDE 95% (ha)","aKDE 50% (ha)","Modell","ESS"),
      row=3)
  xl2 <- vergleich_dt[order(igel), .(
    igel, n_fixes, n_naechte,
    round(akde_95ha, 2), round(akde_50ha, 2),
    modell, round(ess_area, 1)
  )]
  writeData(wb, "aKDE_Uebersicht", xl2, startRow=4, colNames=FALSE)
  addStyle(wb, "aKDE_Uebersicht", st_n, rows=4:(3+nrow(xl2)), cols=c(1,6), gridExpand=TRUE)
  addStyle(wb, "aKDE_Uebersicht", st_0, rows=4:(3+nrow(xl2)), cols=2:3,    gridExpand=TRUE)
  addStyle(wb, "aKDE_Uebersicht", st_2, rows=4:(3+nrow(xl2)), cols=c(4,5,7), gridExpand=TRUE)

  # Sheet 3: aKDE Paaranalyse (falls vorhanden)
  if (!is.null(akde_overlap_dt) && nrow(akde_overlap_dt) > 0) {
    addWorksheet(wb, "aKDE_Paaranalyse", tabColour="#e76f51")
    hdr(wb, "aKDE_Paaranalyse",
        c("Paar","Igel A","Igel B","Datum","Jaccard 95%","OC 95%",
          "Schnitt 95% (ha)","Jaccard 50%","OC 50%","Interpretation"),
        row=1)
    writeData(wb, "aKDE_Paaranalyse", akde_overlap_dt, startRow=2, colNames=FALSE)
    addStyle(wb, "aKDE_Paaranalyse", st_n,
             rows=2:(1+nrow(akde_overlap_dt)), cols=1:10, gridExpand=TRUE)
    addStyle(wb, "aKDE_Paaranalyse", st_2,
             rows=2:(1+nrow(akde_overlap_dt)), cols=5:9, gridExpand=TRUE)
  }

  saveWorkbook(wb, file.path(out_ordner,"Block4b_Uebersicht.xlsx"), overwrite=TRUE)
  cat("  Excel gespeichert\n")
}

# ── 10. Abschluss ─────────────────────────────────────────────
cat("\n== Block 4b abgeschlossen! ==\n")
cat("Output:", out_ordner, "\n")
