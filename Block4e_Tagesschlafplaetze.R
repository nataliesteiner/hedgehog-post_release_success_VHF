# ==============================================================
# Block 4e — daytime resting-site (nest) analysis
# ==============================================================
# Project:  Hedgehog VHF telemetry, Lower Saxony
# Author:   Natalie Steiner | TiHo Hannover
#
# Questions:
#   1. Where are the daytime resting sites?
#      (centroid of daytime fixes per calendar day)
#   2. How far are resting sites from the release enclosure?
#   3. How often do animals change resting site?
#      (shift between consecutive days)
#   4. Does the switching frequency change over time?
#      (LME: nest-switch distance ~ days since release)
#
# Methods:
#   Daytime resting site = median centroid of all daytime fixes
#   (Night = 0) on a calendar day. Nest switch = distance between
#   consecutive daytime resting sites. Threshold 50 m
#   (based on tRackIT positional accuracy ~20-50 m).

# ── 0. Pakete ──────────────────────────────────────────────────
pakete <- c("sf","data.table","ggplot2","patchwork","scales",
            "lubridate","readxl","openxlsx","nlme","viridis")
fehlend <- pakete[!sapply(pakete, requireNamespace, quietly=TRUE)]
if (length(fehlend) > 0) install.packages(fehlend)

suppressPackageStartupMessages({
  library(sf); library(data.table); library(ggplot2)
  library(patchwork); library(scales); library(lubridate)
  library(readxl); library(lme4); library(lmerTest); library(viridis)
})
if (requireNamespace("openxlsx")) library(openxlsx)
cat("Pakete geladen\n\n")

# ── 1. Pfade ───────────────────────────────────────────────────
projekt_root <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
gpkg_ordner  <- file.path(projekt_root, "data", "kernel_files")
meta_datei   <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
out_ordner   <- file.path(projekt_root, "output", "Block4e_Schlafplaetze")
dir.create(out_ordner, showWarnings=FALSE, recursive=TRUE)

ziel_epsg  <- 32632
release_x  <- 514743.7
release_y  <- 5805363.8
WECHSEL_SCHWELLE_M <- 50  # Mindestdistanz fuer einen echten Nestwechsel

cat("Output:", out_ordner, "\n\n")

# ── 2. Metadaten laden ─────────────────────────────────────────
meta <- tryCatch({
  m <- as.data.table(read_excel(path.expand(meta_datei), sheet="Raw_data"))
  data.table(
    igel         = trimws(as.character(m$individual)),
    release_date = as.Date(m$date_release),
    sex          = as.character(m$sex)
  )
}, error=function(e) {
  cat("[WARN] Metadaten:", conditionMessage(e), "\n")
  data.table(igel=character(), release_date=as.Date(NA), sex=character())
})
cat(sprintf("Metadaten: %d Igel\n\n", nrow(meta)))

# ── 3. GPKG-Dateien laden — NUR Tagesfixes ────────────────────
# Tagesschlafplatz = Zentroid der Tagesfixes (Night = 0)
# Tagesfixes repraesentieren den Aufenthaltsort im Ruhezustand
cat("Lade Tagesfixes aus GPKG...\n")

gpkg_files <- list.files(path.expand(gpkg_ordner), pattern="\\.gpkg$",
                          full.names=TRUE)
gpkg_files <- gpkg_files[!grepl("\\(1\\)", gpkg_files)]

tag_fixes_liste <- lapply(gpkg_files, function(f) {
  igel_i <- sub("_Sachsenhagen_.*\\.gpkg$", "", basename(f))
  tryCatch({
    gdf <- st_read(f, quiet=TRUE)

    # Nur Tagesfixes (Night = 0 oder NA)
    if ("Night" %in% names(gdf)) {
      gdf <- gdf[is.na(gdf$Night) | gdf$Night == 0, ]
    }
    gdf <- gdf[!st_is_empty(gdf) & !is.na(st_geometry(gdf)), ]
    if (nrow(gdf) < 3) return(NULL)

    # Koordinaten extrahieren
    geom_col <- attr(gdf, "sf_column")
    if (is.null(geom_col)) geom_col <- "geometry"
    coords <- tryCatch(
      st_coordinates(gdf),
      error=function(e) NULL
    )
    if (is.null(coords) || anyNA(coords[,1])) {
      geom_list <- gdf[[geom_col]]
      x_v <- vapply(geom_list, function(g) { v <- tryCatch(as.numeric(g),error=function(e) c(NA_real_,NA_real_)); if(length(v)>=1) v[1] else NA_real_ }, numeric(1))
      y_v <- vapply(geom_list, function(g) { v <- tryCatch(as.numeric(g),error=function(e) c(NA_real_,NA_real_)); if(length(v)>=2) v[2] else NA_real_ }, numeric(1))
    } else {
      x_v <- coords[,1]; y_v <- coords[,2]
    }

    # Zeitstempel
    tc <- intersect(c("X_time","_time","Date"), names(gdf))[1]
    if (!is.na(tc) && tc == "X_time") {
      ts <- with_tz(as.POSIXct(gdf[[tc]], tz="UTC"), "Europe/Berlin")
    } else if (!is.na(tc) && tc == "Date") {
      ts <- as.POSIXct(as.character(gdf[[tc]]), tz="Europe/Berlin")
    } else {
      ts <- rep(as.POSIXct(NA), nrow(gdf))
    }

    # Release-Datum
    rel <- meta[igel == igel_i, release_date]
    rel <- if (length(rel) && !is.na(rel[1])) rel[1] else as.Date(NA)

    dt <- data.table(
      igel         = igel_i,
      datetime     = ts,
      datum        = as.Date(ts, tz="Europe/Berlin"),
      x            = x_v,
      y            = y_v,
      release_date = rel,
      tage_seit    = as.integer(as.Date(ts, tz="Europe/Berlin") - rel)
    )
    dt <- dt[!is.na(x) & !is.na(y) & x != 0 & y != 0]
    dt
  }, error=function(e) {
    cat(sprintf("  [WARN] %s: %s\n", igel_i, conditionMessage(e)))
    NULL
  })
})

tag_fixes <- rbindlist(Filter(Negate(is.null), tag_fixes_liste), fill=TRUE)
tag_fixes <- merge(tag_fixes, meta[, .(igel, sex)], by="igel", all.x=TRUE)
tag_fixes <- tag_fixes[!is.na(x) & !is.na(y) & !is.na(datum)]

cat(sprintf("  %d Tagesfixes, %d Igel\n\n", nrow(tag_fixes), tag_fixes[, uniqueN(igel)]))

if (nrow(tag_fixes) == 0) stop("Keine Tagesfixes gefunden — GPKG-Laden pruefen.")

# ── 4. Tagesschlafplatz pro Tag berechnen ─────────────────────
cat("Berechne Tagesschlafplaetze (Median-Zentroid pro Tag)...\n")

schlaf_dt <- tag_fixes[, .(
  nest_x     = median(x, na.rm=TRUE),
  nest_y     = median(y, na.rm=TRUE),
  n_fixes    = .N,
  tage_seit  = first(na.omit(tage_seit))
), by = .(igel, datum)]

# Nur Tage mit genuegend Fixes
schlaf_dt <- schlaf_dt[n_fixes >= 3]

# Distanz zur Auswilderungsvoliere
schlaf_dt[, dist_voliere_m := sqrt((nest_x - release_x)^2 + (nest_y - release_y)^2)]

# Release-Datum ergaenzen
schlaf_dt <- merge(schlaf_dt, meta[, .(igel, release_date, sex)], by="igel", all.x=TRUE)

# ── N-at-risk-Cutoff (Survivorship Bias Korrektur) ─────────────
MIN_N_AT_RISK <- 5

n_risk_per_day <- schlaf_dt[!is.na(tage_seit) & tage_seit >= 0,
                              .(n_igel = uniqueN(igel)), by = tage_seit]
cutoff_tag <- n_risk_per_day[n_igel >= MIN_N_AT_RISK, max(tage_seit)]
cat(sprintf("  N-at-risk-Cutoff: Tag %d (letzter Tag mit N >= %d Tieren)\n",
            cutoff_tag, MIN_N_AT_RISK))

cat("  N-at-risk pro Tag:\n")
print(n_risk_per_day[tage_seit <= cutoff_tag + 3][order(tage_seit)])

schlaf_ok  <- schlaf_dt[!is.na(tage_seit) & tage_seit >= 0 & tage_seit <= cutoff_tag]
cat(sprintf("\n  Nach Cutoff: %d Tagesschlafplaetze (%d -> %d Tage)\n\n",
            nrow(schlaf_ok), schlaf_ok[, min(tage_seit)], schlaf_ok[, max(tage_seit)]))

cat(sprintf("  %d Tagesschlafplaetze bei %d Igeln (Median %d Tage/Tier)\n\n",
            nrow(schlaf_dt),
            schlaf_dt[, uniqueN(igel)],
            as.integer(median(schlaf_dt[, .N, by=igel]$N))))

# ── 5. Nestwechsel berechnen ───────────────────────────────────
cat("Berechne Nestwechsel zwischen aufeinanderfolgenden Tagen...\n")

wechsel_dt <- rbindlist(lapply(schlaf_dt[, unique(igel)], function(ig) {
  d <- schlaf_dt[igel == ig][order(datum)]
  if (nrow(d) < 2) return(NULL)

  d[, wechsel_m := c(NA_real_,
    sqrt(diff(nest_x)^2 + diff(nest_y)^2))]
  d[, tage_intervall := c(NA_integer_, as.integer(diff(datum)))]

  d[!is.na(wechsel_m) & tage_intervall <= 2]
}), fill=TRUE)

wechsel_dt[, echter_wechsel := wechsel_m > WECHSEL_SCHWELLE_M]

cat(sprintf("  %d Tages-zu-Tages-Vergleiche bei %d Igeln\n",
            nrow(wechsel_dt), wechsel_dt[, uniqueN(igel)]))
cat(sprintf("  Echter Nestwechsel (>%d m): %d von %d (%.0f%%)\n\n",
            WECHSEL_SCHWELLE_M,
            sum(wechsel_dt$echter_wechsel, na.rm=TRUE),
            nrow(wechsel_dt),
            100*mean(wechsel_dt$echter_wechsel, na.rm=TRUE)))

# ── 6. Deskriptive Statistik ───────────────────────────────────
cat("=== Deskriptive Statistik ===\n")

tier_stats <- schlaf_dt[, .(
  N_Tage         = .N,
  Median_Dist_m  = round(median(dist_voliere_m, na.rm=TRUE)),
  Max_Dist_m     = round(max(dist_voliere_m, na.rm=TRUE))
), by = igel]

tier_wechsel <- wechsel_dt[, .(
  N_Vergleiche       = .N,
  Median_Wechsel_m   = round(median(wechsel_m, na.rm=TRUE)),
  Pct_Echter_Wechsel = round(100*mean(echter_wechsel, na.rm=TRUE), 1),
  N_Unique_Nester    = .N
), by = igel]

tier_summary <- merge(tier_stats, tier_wechsel, by="igel", all.x=TRUE)
tier_summary <- merge(tier_summary, meta[, .(igel, sex)], by="igel", all.x=TRUE)
setorder(tier_summary, igel)

cat(sprintf("Population (n=%d Igel):\n", nrow(tier_summary)))
cat(sprintf("  Median Nest-Distanz zur Voliere: %d m (Range: %d-%d m)\n",
            as.integer(median(tier_summary$Median_Dist_m, na.rm=TRUE)),
            as.integer(min(tier_summary$Median_Dist_m, na.rm=TRUE)),
            as.integer(max(tier_summary$Median_Dist_m, na.rm=TRUE))))
cat(sprintf("  Median taegliche Nestwechsel-Distanz: %d m\n",
            as.integer(median(wechsel_dt$wechsel_m, na.rm=TRUE))))
cat(sprintf("  Anteil echter Nestwechsel (>%d m): %.0f%%\n",
            WECHSEL_SCHWELLE_M,
            100*mean(wechsel_dt$echter_wechsel, na.rm=TRUE)))
print(tier_summary)
cat("\n")

# ── 6b. Sexvergleich ──────────────────────────────────────────
cat("=== Sexvergleich: Nest-Distanz zur Voliere ===\n")
sex_dist   <- tier_summary[!is.na(sex) & !is.na(Median_Dist_m)]
sex_male   <- sex_dist[sex == "Male",   Median_Dist_m]
sex_female <- sex_dist[sex == "Female", Median_Dist_m]

cat(sprintf("  Male   n=%d: median=%d m (range %d-%d m)\n",
            length(sex_male), as.integer(median(sex_male)),
            min(sex_male), max(sex_male)))
cat(sprintf("  Female n=%d: median=%d m (range %d-%d m)\n",
            length(sex_female), as.integer(median(sex_female)),
            min(sex_female), max(sex_female)))

if (length(sex_male) >= 3 && length(sex_female) >= 3) {
  wt <- wilcox.test(sex_female, sex_male, alternative = "two.sided", exact = FALSE)
  cat(sprintf("  Wilcoxon rank-sum: W=%.1f, p=%.4f (%s)\n",
              wt$statistic, wt$p.value,
              ifelse(wt$p.value < 0.05, "signifikant", "nicht signifikant")))
}

cat("\n=== Sexvergleich: Nestwechsel-Distanz ===\n")
wechsel_sex <- wechsel_dt[!is.na(wechsel_m) & !is.na(sex)]
for (s in c("Male","Female")) {
  vals <- wechsel_sex[sex == s, wechsel_m]
  cat(sprintf("  %s n=%d shifts: median=%.0f m\n", s, length(vals), median(vals)))
}
wm <- wechsel_sex[sex == "Male",   wechsel_m]
wf <- wechsel_sex[sex == "Female", wechsel_m]
if (length(wm) >= 5 && length(wf) >= 5) {
  wt2 <- wilcox.test(wf, wm, alternative = "two.sided", exact = FALSE)
  cat(sprintf("  Wilcoxon: W=%.1f, p=%.4f\n", wt2$statistic, wt2$p.value))
}
cat("\n")

# ── 7. LME: Veraendert sich Wechseldistanz ueber die Zeit? ────
cat("LME: Nestwechsel-Distanz ~ Tage seit Auswilderung...\n")

lme_data <- wechsel_dt[!is.na(wechsel_m) & !is.na(tage_seit) & !is.na(igel) &
                        tage_seit >= 0 & tage_seit <= cutoff_tag]
cat(sprintf("  LME-Daten: %d Vergleiche bei %d Igeln (Tag 0-%d, N-at-risk >= %d)\n",
            nrow(lme_data), lme_data[, uniqueN(igel)], cutoff_tag, MIN_N_AT_RISK))

if (nrow(lme_data) >= 20 && lme_data[, uniqueN(igel)] >= 5) {
  lme_res <- tryCatch(
    lmer(wechsel_m ~ tage_seit + (1|igel), data=lme_data,
         REML=TRUE, na.action=na.omit),
    error=function(e) { cat("[WARN] LME:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(lme_res)) {
    s <- summary(lme_res)
    cat(sprintf("  Intercept: %.1f m (SE: %.1f)\n", coef(s)[1,1], coef(s)[1,2]))
    cat(sprintf("  Slope (tage_seit): %.2f m/Tag (SE: %.2f, p=%.4f)\n",
                coef(s)[2,1], coef(s)[2,2], coef(s)[2,5]))
    cat(sprintf("  Interpretation: Nestwechsel-Distanz %s signifikant ueber die Zeit\n\n",
                ifelse(coef(s)[2,5] < 0.05,
                       ifelse(coef(s)[2,1] < 0, "NIMMT AB", "NIMMT ZU"),
                       "veraendert sich NICHT")))
  }
} else {
  cat("  Zu wenig Daten fuer LME\n\n")
  lme_res <- NULL
}

# ── 8. Plots ──────────────────────────────────────────────────
cat("Erstelle Plots...\n")

farben_sex <- c("Male"="#4575b4","Female"="#d73027")

p_dist_zeit <- ggplot(schlaf_ok, aes(x=tage_seit, y=dist_voliere_m)) +
  geom_point(aes(color=sex), size=1.5, alpha=0.4) +
  geom_smooth(method="loess", span=0.6, se=TRUE,
              color="grey30", fill="grey80", linewidth=1.0) +
  geom_hline(yintercept=100, linetype="dashed", color="#2166ac", linewidth=0.7) +
  scale_color_manual(values=farben_sex, name="Sex",
                     guide=guide_legend(override.aes=list(size=3,alpha=1))) +
  scale_y_continuous(labels=label_number(suffix=" m")) +
  labs(title="Day nest distance from release enclosure over time",
       subtitle="Dashed line = 100 m | LOESS smooth with 95% CI",
       x="Days since release", y="Distance to enclosure (m)") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold"), legend.position="right")

p_wechsel_zeit <- ggplot(
  wechsel_dt[!is.na(tage_seit) & tage_seit >= 0 & tage_seit <= cutoff_tag],
  aes(x=tage_seit, y=wechsel_m)) +
  geom_hline(yintercept=WECHSEL_SCHWELLE_M, linetype="dashed",
             color="grey50", linewidth=0.7) +
  geom_point(aes(color=echter_wechsel), size=1.8, alpha=0.5) +
  geom_smooth(method="loess", span=0.7, se=TRUE,
              color="grey30", fill="grey80", linewidth=1.0) +
  { if (!is.null(lme_res)) {
      b <- coef(summary(lme_res))[,1]
      geom_abline(slope=b[2], intercept=b[1],
                  color="#c0392b", linewidth=1.0, linetype="dotted")
  }} +
  scale_color_manual(values=c("TRUE"="#c0392b","FALSE"="#4575b4"),
                     labels=c("TRUE"=paste0(">",WECHSEL_SCHWELLE_M,"m (nest change)"),
                              "FALSE"=paste0("<=",WECHSEL_SCHWELLE_M,"m (same site)")),
                     name=NULL) +
  scale_y_continuous(labels=label_number(suffix=" m"), limits=c(0, NA)) +
  labs(title="Day-to-day nest shift distance over time",
       subtitle=paste0("Dashed = ", WECHSEL_SCHWELLE_M,
                       "m threshold | Red dotted = LME trend"),
       x="Days since release", y="Nest shift distance (m)") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")

igel_ord <- tier_summary[order(Median_Dist_m), igel]
schlaf_plot <- copy(schlaf_dt)
schlaf_plot[, igel_f := factor(igel, levels=igel_ord)]

p_dist_tier <- ggplot(schlaf_plot, aes(x=igel_f, y=dist_voliere_m, fill=sex)) +
  geom_boxplot(alpha=0.7, outlier.shape=21, outlier.size=1.5, width=0.7) +
  geom_hline(yintercept=c(100,250), linetype=c("dashed","dotted"),
             color=c("#2166ac","grey50"), linewidth=0.7) +
  scale_fill_manual(values=farben_sex, name="Sex") +
  scale_y_continuous(labels=label_number(suffix=" m")) +
  coord_flip() +
  annotate("text", x=0.5, y=105, label="100 m", size=3, color="#2166ac", hjust=0) +
  annotate("text", x=0.5, y=255, label="250 m detection limit", size=3, color="grey50", hjust=0) +
  labs(title="Day nest distance per individual",
       subtitle="Sorted by median distance | Dashed = 100m | Dotted = detection limit",
       x=NULL, y="Distance to release enclosure (m)") +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="right")

p_karte_nester <- ggplot(schlaf_dt, aes(x=nest_x, y=nest_y, color=tage_seit)) +
  geom_point(size=1.5, alpha=0.6) +
  geom_point(aes(x=release_x, y=release_y), color="red", shape=8,
             size=6, stroke=2, inherit.aes=FALSE) +
  scale_color_viridis_c(option="plasma", name="Days since\nrelease") +
  facet_wrap(~igel, ncol=5, scales="free") +
  labs(title="Individual day nest locations over time",
       subtitle="Red star = release enclosure | Colour = days since release",
       x="UTM32N Easting (m)", y="UTM32N Northing (m)") +
  theme_bw(base_size=8) +
  theme(plot.title=element_text(face="bold", size=11),
        strip.text=element_text(face="bold", size=7),
        axis.text=element_text(size=5),
        legend.position="right")

p_kombi <- (p_dist_zeit + p_wechsel_zeit) /
            p_dist_tier +
  plot_annotation(
    title    = "Day nest analysis — post-release hedgehogs",
    subtitle = sprintf("N = %d hedgehogs | Day 0-%d (N-at-risk >= %d) | Threshold: >%d m",
                       schlaf_dt[, uniqueN(igel)], cutoff_tag, MIN_N_AT_RISK, WECHSEL_SCHWELLE_M),
    theme    = theme(plot.title=element_text(face="bold", size=13))
  )

ggsave(file.path(out_ordner, "tagesschlafplaetze_uebersicht.png"),
       p_kombi, width=14, height=12, dpi=150)
ggsave(file.path(out_ordner, "tagesschlafplaetze_karte.png"),
       p_karte_nester, width=16, height=12, dpi=150)
cat("  tagesschlafplaetze_uebersicht.png / _karte.png\n\n")

# ── 9. Excel Export ───────────────────────────────────────────
cat("Speichere Ergebnisse...\n")
wb <- createWorkbook()

addWorksheet(wb, "Zusammenfassung", tabColour="#2C5F8A")
writeData(wb, "Zusammenfassung",
          data.frame(Info=paste0(
            "Tagesschlafplatz-Analyse | Schwellenwert Nestwechsel: >",
            WECHSEL_SCHWELLE_M, " m | Positionen: Median-Zentroid der Tagesfixes (Night=0) pro Kalendertag")),
          startRow=1, colNames=FALSE)
writeData(wb, "Zusammenfassung", as.data.frame(tier_summary), startRow=2, colNames=TRUE)

addWorksheet(wb, "Schlafplaetze_taeglich", tabColour="#1a6b1a")
out_nester <- schlaf_dt[order(igel, datum), .(
  Igel=igel, Datum=as.character(datum), Tage_seit=tage_seit,
  Nest_X=round(nest_x,1), Nest_Y=round(nest_y,1),
  Dist_Voliere_m=round(dist_voliere_m,1), N_Fixes=n_fixes
)]
writeData(wb, "Schlafplaetze_taeglich", as.data.frame(out_nester), startRow=1, colNames=TRUE)

addWorksheet(wb, "Nestwechsel", tabColour="#c97a1e")
out_wechsel <- wechsel_dt[order(igel, datum), .(
  Igel=igel, Datum=as.character(datum), Tage_seit=tage_seit,
  Wechsel_m=round(wechsel_m,1), Tage_Intervall=tage_intervall,
  Echter_Wechsel=echter_wechsel
)]
writeData(wb, "Nestwechsel", as.data.frame(out_wechsel), startRow=1, colNames=TRUE)

if (!is.null(lme_res)) {
  addWorksheet(wb, "LME_Ergebnis", tabColour="#d73027")
  s <- summary(lme_res)
  lme_tab <- data.frame(
    Term    = rownames(coef(s)),
    Estimate= round(coef(s)[,1],3),
    SE      = round(coef(s)[,2],3),
    t_value = round(coef(s)[,4],3),
    p_value = round(coef(s)[,5],4)
  )
  writeData(wb, "LME_Ergebnis",
            data.frame(Info="LME: Nestwechsel-Distanz (m) ~ Tage seit Auswilderung + (1|Igel)"),
            startRow=1, colNames=FALSE)
  writeData(wb, "LME_Ergebnis", lme_tab, startRow=2, colNames=TRUE)
}

# ── 10. Habitatklassifikation (ATKIS Basis-DLM) ───────────────
cat("Klassifiziere Habitattyp der Tagesschlafplaetze...\n")
cat("  Datenquelle: ATKIS Basis-DLM (LGLN Niedersachsen, CC BY 4.0)\n")

# ATKIS Basis-DLM — amtliche Geobasisdaten
# Einmalig Habitat_Download_ATKIS.R ausfuehren, um das GPKG zu erzeugen.
hab_pfad <- file.path(projekt_root, "data", "excel_files", "habitat_sachsenhagen_atkis.gpkg")

if (!file.exists(path.expand(hab_pfad))) {
  cat("  ATKIS-GPKG nicht gefunden — Habitat-Analyse uebersprungen\n")
  cat("  Fehlend:", hab_pfad, "\n")
  cat("  Loesung: Habitat_Download_ATKIS.R einmalig ausfuehren.\n")
} else {

  hab_sf_e <- tryCatch(
    st_transform(st_read(path.expand(hab_pfad), quiet=TRUE), crs=ziel_epsg),
    error=function(e) { cat(sprintf("  [FEHLER] GPKG laden: %s\n", conditionMessage(e))); NULL })

  if (!is.null(hab_sf_e)) {
    hab_sf_e <- hab_sf_e[st_is_valid(hab_sf_e), ]
    nester_sf_e <- st_as_sf(schlaf_ok, coords=c("nest_x","nest_y"), crs=ziel_epsg)

    # Basis-Klassifikation (ATKIS-Klassen: Wald, Acker, Wiese, Gehoelz_Strauch)
    schlaf_ok[, nest_hab_roh := "Offenland"]
    for (hab in c("Gehoelz_Strauch","Gewaesser","Wald","Wiese","Acker")) {
      hab_poly <- hab_sf_e[hab_sf_e$habitat == hab, ]
      if (nrow(hab_poly) == 0) next
      in_hab <- tryCatch(
        suppressMessages(lengths(st_intersects(nester_sf_e, st_union(hab_poly))) > 0),
        error=function(e) rep(FALSE, nrow(schlaf_ok)))
      schlaf_ok[in_hab & nest_hab_roh == "Offenland", nest_hab_roh := hab]
    }

    # Waldrand-Zone (±20 m) — konsistent mit Block5_Umweltdaten.R
    wald_u <- tryCatch(st_union(hab_sf_e[hab_sf_e$habitat == "Wald", ]), error=function(e) NULL)
    if (!is.null(wald_u)) {
      wald_rand_e <- tryCatch(
        st_make_valid(st_difference(st_buffer(wald_u, 20), st_buffer(wald_u, -20))),
        error=function(e) NULL)
      if (!is.null(wald_rand_e)) {
        in_rand <- tryCatch(
          suppressMessages(lengths(st_intersects(nester_sf_e, wald_rand_e)) > 0),
          error=function(e) rep(FALSE, nrow(schlaf_ok)))
        schlaf_ok[in_rand, nest_hab_roh := "Waldrand"]
      }
    }

    # Zusammenfassen (konsistent mit Block5: Edge = Waldrand + Gehoelz_Strauch)
    # Waldweg wird bei Tagesschlafen NICHT beruecksichtigt:
    # Schlafende Tiere liegen nicht auf Waldwegen (Punkt-Koordinate, kein Naehepuffer)
    schlaf_ok[, nest_habitat := fcase(
      nest_hab_roh %in% c("Waldrand","Gehoelz_Strauch"), "Edge habitat",
      nest_hab_roh == "Wald",  "Forest interior",
      nest_hab_roh == "Acker", "Agricultural land",
      nest_hab_roh == "Wiese", "Grassland",
      default = "Other"
    )]

    schlaf_ok[, phase := fifelse(tage_seit <= 5, "Early (day 1-5)", "Late (day 6+)")]

    cat("\n  === Habitatverteilung Tagesschlafplaetze (ATKIS) ===\n")
    hab_vert_e <- schlaf_ok[, .N, by=nest_habitat]
    hab_vert_e[, use_pct := round(100*N/sum(N),1)]
    print(hab_vert_e[order(-use_pct)])

    # ── Verfuegbarkeit (600 m Radius, ATKIS) ─────────────────────
    release_sf_e <- st_as_sf(data.frame(x=release_x, y=release_y),
                              coords=c("x","y"), crs=ziel_epsg)
    study_buf_e  <- st_buffer(release_sf_e, 600)
    hab_av_e     <- suppressMessages(st_intersection(hab_sf_e, study_buf_e))
    hab_av_e$area_ha <- as.numeric(st_area(hab_av_e)) / 10000
    avail_r <- as.data.table(hab_av_e)[, .(avail_ha=sum(area_ha)), by=habitat]

    hab_map2 <- c(
      "Wald"           = "Forest interior",
      "Acker"          = "Agricultural land",
      "Wiese"          = "Grassland",
      "Gehoelz_Strauch"= "Edge habitat"
    )
    avail_r[, hab_grob := hab_map2[habitat]]
    avail_g <- avail_r[!is.na(hab_grob), .(avail_ha=sum(avail_ha)), by=hab_grob]

    # Waldrand-Flaeche zu Edge habitat hinzufuegen, von Forest interior abziehen
    edge_ha_e2 <- 0
    if (!is.null(wald_rand_e)) {
      rand_s <- tryCatch(suppressMessages(st_intersection(wald_rand_e, study_buf_e)),
                          error=function(e) NULL)
      if (!is.null(rand_s)) edge_ha_e2 <- as.numeric(st_area(rand_s)) / 10000
    }
    if ("Forest interior" %in% avail_g$hab_grob)
      avail_g[hab_grob=="Forest interior", avail_ha := pmax(0, avail_ha - edge_ha_e2)]
    if ("Edge habitat" %in% avail_g$hab_grob)
      avail_g[hab_grob=="Edge habitat", avail_ha := avail_ha + edge_ha_e2]
    else
      avail_g <- rbind(avail_g, data.table(hab_grob="Edge habitat", avail_ha=edge_ha_e2))

    avail_g[, avail_pct := avail_ha / sum(avail_ha)]
    setnames(avail_g, "hab_grob", "nest_habitat")

    uva_n <- merge(hab_vert_e, avail_g[, .(nest_habitat, avail_pct)],
                   by="nest_habitat", all=TRUE)
    uva_n[is.na(use_pct), use_pct := 0]
    uva_n[is.na(N), N := 0]
    uva_n[is.na(avail_pct), avail_pct := 0.001]
    uva_n[, avail_pct_pct     := round(avail_pct*100, 1)]
    uva_n[, selektionsratio   := round(use_pct/100 / avail_pct, 2)]
    uva_n[, praeferenz        := fcase(
      selektionsratio>1.2, "Preferred",
      selektionsratio<0.8, "Avoided",
      default = "Neutral")]
    setorder(uva_n, -selektionsratio)

    cat("\n  Use vs. Availability (ATKIS):\n")
    print(uva_n[, .(nest_habitat, N_Nester=N, Genutzt_pct=use_pct,
                    Verfuegbar_pct=avail_pct_pct, Ratio=selektionsratio)])

    # ── Per-Tier-Wilcoxon-Test (Hauptanalyse) ────────────────────
    wilcox_nest_dt <- NULL
    tryCatch({
      n_tiere_n <- length(unique(schlaf_ok$igel))
      cat(sprintf("\n  === Per-Tier-Wilcoxon (Tagnester, n=%d Tiere) ===\n", n_tiere_n))

      pt_nest     <- schlaf_ok[, .N, by=.(igel, nest_habitat)]
      pt_nest_tot <- schlaf_ok[, .(N_total=.N), by=igel]
      pt_nest     <- merge(pt_nest, pt_nest_tot, by="igel")
      pt_nest[, prop := N / N_total]

      hab_test_n2 <- c("Edge habitat","Forest interior","Agricultural land","Grassland")
      wilcox_nest_dt <- rbindlist(lapply(hab_test_n2, function(h) {
        avail_p <- uva_n[nest_habitat==h, avail_pct]
        if (length(avail_p)==0 || is.na(avail_p)) return(NULL)
        props_with <- pt_nest[nest_habitat==h, prop]
        all_tiere_n <- unique(schlaf_ok$igel)
        n_without   <- length(setdiff(all_tiere_n, pt_nest[nest_habitat==h, igel]))
        props_all   <- c(props_with, rep(0, n_without))
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
      cat("\n  Per-Tier-Wilcoxon-Ergebnisse (Tagesnester, ATKIS):\n")
      print(wilcox_nest_dt)
    }, error=function(e) {
      cat(sprintf("  [WARN] Wilcoxon: %s\n", conditionMessage(e)))
      wilcox_nest_dt <<- NULL
    })

    # ── Excel-Sheets ──────────────────────────────────────────────
    tryCatch({
      addWorksheet(wb, "Nest_Habitat_Rohdaten", tabColour="#c97a1e")
      out_r <- as.data.frame(schlaf_ok[order(igel,datum),
        .(Igel=igel, Datum=as.character(datum), Tage_seit=tage_seit, Phase=phase,
          Nest_X=round(nest_x,1), Nest_Y=round(nest_y,1),
          Dist_Voliere_m=round(dist_voliere_m,1),
          Habitat_roh=nest_hab_roh, Habitat=nest_habitat, N_Fixes=n_fixes)])
      writeData(wb, "Nest_Habitat_Rohdaten",
                data.frame(Info=paste0(
                  "Habitatklassifikation ATKIS | Edge=Waldrand+-20m+Gehoelz_Strauch | Tag 0-",cutoff_tag)),
                startRow=1, colNames=FALSE)
      writeData(wb, "Nest_Habitat_Rohdaten", out_r, startRow=2, colNames=TRUE)
      cat(sprintf("    Nest_Habitat_Rohdaten: %d Zeilen\n", nrow(out_r)))
    }, error=function(e) cat(sprintf("    [FEHLER] Rohdaten: %s\n", conditionMessage(e))))

    tryCatch({
      addWorksheet(wb, "Nest_Habitat_UvA", tabColour="#d4a843")
      writeData(wb, "Nest_Habitat_UvA",
                data.frame(Info=paste0(
                  "Use vs. Availability (ATKIS) | Selektionsratio>1=bevorzugt | Tag 0-",cutoff_tag)),
                startRow=1, colNames=FALSE)
      writeData(wb, "Nest_Habitat_UvA",
                as.data.frame(uva_n[order(-selektionsratio),
                  .(Habitat=nest_habitat, N_Nester=N,
                    Genutzt_pct=use_pct, Verfuegbar_pct=avail_pct_pct,
                    Selektionsratio=selektionsratio, Praeferenz=praeferenz)]),
                startRow=2, colNames=TRUE)
    }, error=function(e) cat(sprintf("    [FEHLER] UvA: %s\n", conditionMessage(e))))

    if (!is.null(wilcox_nest_dt) && nrow(wilcox_nest_dt) > 0) {
      tryCatch({
        addWorksheet(wb, "Nest_Hab_Wilcoxon", tabColour="#2166ac")
        writeData(wb, "Nest_Hab_Wilcoxon",
                  data.frame(Info=paste0(
                    "HAUPTANALYSE: Per-Tier-Wilcoxon | ATKIS Basis-DLM Habitatdaten | ",
                    "n = Tiere mit Nestdaten | *** p<0.001 | ** p<0.01 | * p<0.05 | n.s. n.s.")),
                  startRow=1, colNames=FALSE)
        writeData(wb, "Nest_Hab_Wilcoxon", as.data.frame(wilcox_nest_dt),
                  startRow=2, colNames=TRUE)
        cat(sprintf("    Nest_Hab_Wilcoxon: %d Habitate\n", nrow(wilcox_nest_dt)))
      }, error=function(e) cat(sprintf("    [FEHLER] Wilcoxon-Sheet: %s\n", conditionMessage(e))))
    }

    tryCatch({
      addWorksheet(wb, "Nest_Habitat_Phase", tabColour="#f7c948")
      ph_wide <- dcast(
        schlaf_ok[!is.na(phase), .(N=.N), by=.(igel,phase,nest_habitat)][
          , pct:=round(100*N/sum(N),1), by=.(igel,phase)],
        igel + nest_habitat ~ phase, value.var="pct", fill=0)
      writeData(wb, "Nest_Habitat_Phase",
                data.frame(Info="% Nester pro Habitat und Phase | Frueh=Tag 1-5 | Spaet=Tag 6+"),
                startRow=1, colNames=FALSE)
      writeData(wb, "Nest_Habitat_Phase", as.data.frame(ph_wide), startRow=2, colNames=TRUE)
    }, error=function(e) cat(sprintf("    [FEHLER] Phase: %s\n", conditionMessage(e))))

  } # Ende if (!is.null(hab_sf_e))
} # Ende if (file.exists)

saveWorkbook(wb, file.path(out_ordner, "Block4e_Schlafplaetze_Ergebnisse.xlsx"), overwrite=TRUE)
cat("  Block4e_Schlafplaetze_Ergebnisse.xlsx\n\n")

cat("====================================================\n")
cat("Block 4e abgeschlossen!\n")
cat("====================================================\n")
cat("Outputs:\n")
cat("  tagesschlafplaetze_uebersicht.png\n")
cat("  tagesschlafplaetze_karte.png\n")
cat("  Block4e_Schlafplaetze_Ergebnisse.xlsx\n\n")
cat(sprintf("Kernbefunde:\n"))
cat(sprintf("  N Tiere analysiert: %d\n", schlaf_dt[, uniqueN(igel)]))
cat(sprintf("  Median Nest-Distanz: %d m zur Voliere\n",
            as.integer(median(schlaf_dt$dist_voliere_m, na.rm=TRUE))))
cat(sprintf("  Anteil echter Nestwechsel (>%d m): %.0f%%\n",
            WECHSEL_SCHWELLE_M, 100*mean(wechsel_dt$echter_wechsel, na.rm=TRUE)))
