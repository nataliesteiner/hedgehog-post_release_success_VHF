# ==============================================================================
#  Block6_Dispersal.R
#  Hedgehog VHF telemetry, Wildtierstation Sachsenhagen
#
#  Question: How far and how fast do the hedgehogs move away from the
#            release enclosure after release?
#
#  Hypothesis: In the first days after release (soft-release phase) the
#              animals stay near the enclosure and disperse only later.
#
#  Data: VHF multilaterations (GeoPackages), all 26 hedgehogs
#
#  Method adapted from:
#    Edible-dormouse project (Steiner N.) — 05_softrelease.R
#    Same site, same antenna infrastructure
#
#  Input:
#    data/kernel_files/multilateration_April2026/*.multilaterations.gpkg
#    data/excel_files/data_igel.xlsx

# ── 0. Pakete ──────────────────────────────────────────────────────────────────

pakete <- c("sf", "lme4", "lmerTest", "ggplot2", "dplyr", "tidyr",
            "purrr", "lubridate", "readxl", "openxlsx", "scales",
            "mgcv",     # GAM-Smoother (Plot E)
            "ggdist",   # Raincloud-Plot (Plot F + G)
            "ggrepel"   # Text-Labels anti-overlap (Plot G, I, J, K)
            )
neu <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) {
  cat("Installiere fehlende Pakete:", paste(neu, collapse = ", "), "\n")
  install.packages(neu)
}

suppressPackageStartupMessages({
  library(sf)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lubridate)
  library(readxl)
  library(openxlsx)
  library(scales)
})

cat("✓ Pakete geladen\n\n")

# ── 1. Pfade ───────────────────────────────────────────────────────────────────
# Projektwurzel — einzige Zeile, die du ggf. anpassen musst
projekt_root <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"

gpkg_ordner  <- file.path(projekt_root, "data", "kernel_files")
meta_datei   <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
out_ordner   <- file.path(projekt_root, "output", "Block6_Dispersal")
dir.create(out_ordner, showWarnings = FALSE, recursive = TRUE)

cat("📁 GPKG-Ordner:", gpkg_ordner, "\n")
cat("📁 Output:     ", out_ordner,  "\n\n")

# ── 2. Auswilderungsort (Voliere) ──────────────────────────────────────────────
VOLIERE_LAT      <- 52.398608
VOLIERE_LON      <-  9.216834
VOLIERE_RADIUS_M <- 50    # Visueller Richtwert: ~50 m Voliere-Durchmesser

voliere_wgs <- st_sfc(st_point(c(VOLIERE_LON, VOLIERE_LAT)), crs = 4326)
voliere_utm <- st_transform(voliere_wgs, crs = 25832)
vol_xy      <- st_coordinates(voliere_utm)

DAYS_MAX <- 30L   # Analysefenster: bis 30 Tage nach Auswilderung

cat(sprintf("📍 Voliere UTM: X = %.1f, Y = %.1f\n", vol_xy[1,1], vol_xy[1,2]))

# ── 3. Metadaten laden ─────────────────────────────────────────────────────────
cat("Lade Metadaten aus", basename(meta_datei), "...\n")

meta_raw <- read_excel(meta_datei)
names(meta_raw) <- tolower(gsub("[^a-zA-Z0-9]", "_", names(meta_raw)))

# Spaltenname suchen — gibt character(0) wenn nicht gefunden
get_col <- function(df, patterns) {
  found <- unlist(lapply(patterns, function(p)
    grep(p, names(df), ignore.case = TRUE, value = TRUE)))
  if (length(found) > 0) found[1] else NA_character_
}

col_rel    <- get_col(meta_raw, c("date_release", "release_date"))
col_sr     <- get_col(meta_raw, c("soft_release_start", "soft_release"))
col_last   <- get_col(meta_raw, c("date_last_signal", "last_signal"))
col_indiv  <- get_col(meta_raw, c("^individual$", "individual"))
col_sex    <- get_col(meta_raw, c("^sex$", "sex"))
col_diag   <- get_col(meta_raw, c("diagnosis_main", "diagnosis"))
col_weight <- get_col(meta_raw, c("tagging_weight", "weight"))
col_id     <- get_col(meta_raw, c("^id$"))

cat("  Spalten gefunden:\n")
cat(sprintf("    individual:         %s\n", col_indiv))
cat(sprintf("    date_release:       %s\n", col_rel))
cat(sprintf("    soft_release_start: %s\n", col_sr))
cat(sprintf("    date_last_signal:   %s\n", col_last))
cat(sprintf("    sex:                %s\n", col_sex))
cat(sprintf("    diagnosis:          %s\n", col_diag))

meta <- meta_raw %>%
  transmute(
    tier_id        = if (!is.na(col_id))     .data[[col_id]]    else row_number(),
    tier_label     = trimws(.data[[col_indiv]]),
    date_release   = as.Date(.data[[col_rel]]),
    soft_release_start = if (!is.na(col_sr))   as.Date(.data[[col_sr]])   else as.Date(NA),
    date_last_signal   = if (!is.na(col_last)) as.Date(.data[[col_last]]) else as.Date(NA),
    sex            = if (!is.na(col_sex))    .data[[col_sex]]   else NA_character_,
    diagnosis      = if (!is.na(col_diag))   .data[[col_diag]]  else NA_character_,
    tagging_weight = if (!is.na(col_weight)) as.numeric(.data[[col_weight]]) else NA_real_,
    sr_dauer_tage  = as.numeric(date_release - soft_release_start),
    # tag_found: 1 = Sender physisch auf dem Gelände wiedergefunden, 0 = nicht gefunden
    # Quelle: manuell eingetragen in data_igel.xlsx (Spalte tag_found)
    tag_found = if ("tag_found" %in% names(meta_raw)) as.integer(.data[["tag_found"]]) else NA_integer_,
    # longterm_outcome: "alive" / "dead" / "unknown" — manuell eingetragen
    longterm_outcome = if ("longterm_outcome" %in% names(meta_raw)) as.character(.data[["longterm_outcome"]]) else "unknown"
  ) %>%
  filter(!is.na(tier_label))
# Tiere ohne date_release werden nicht ausgeschlossen — ihr erster Fix-Tag
# wird später als Ersatz-Tag-0 verwendet (Fallback, mit Kennzeichnung)

cat(sprintf("\n  %d Tiere in Metadaten:\n", nrow(meta)))
print(meta %>% select(tier_label, date_release, sr_dauer_tage, sex) %>% as.data.frame())

# ── 4. GeoPackages laden ───────────────────────────────────────────────────────
cat("\n── GeoPackages laden ────────────────────────────────────────\n")

# recursive = TRUE: findet GPKGs auch in Unterordnern (z.B. multilateration_April2026/)
gpkg_dateien <- list.files(gpkg_ordner,
                            pattern    = "\\.gpkg$",
                            full.names = TRUE,
                            recursive  = TRUE)
cat(sprintf("  %d GeoPackage-Dateien gefunden.\n", length(gpkg_dateien)))

extract_label_gpkg <- function(fp) {
  bn <- basename(fp)
  m  <- regmatches(bn, regexpr("Igel\\d+", bn, perl = TRUE))
  if (length(m) > 0) return(m[1])
  sub("_.*", "", bn)
}

find_time_col <- function(sf_obj) {
  candidates <- c("X_time", "_time", "time", "timestamp", "Timestamp", "datetime")
  found <- intersect(candidates, names(sf_obj))
  if (length(found) > 0) return(found[1])
  if (all(c("Date", "Time") %in% names(sf_obj))) return("Date_Time_combined")
  NA_character_
}

alle_fixes <- map_dfr(gpkg_dateien, function(fp) {
  tryCatch({
    sf_obj <- st_read(fp, quiet = TRUE)
    lbl    <- extract_label_gpkg(fp)

    sc_col  <- intersect(c("Station Count", "station_count", "StationCount"), names(sf_obj))[1]
    sc_vals <- if (!is.na(sc_col)) as.integer(sf_obj[[sc_col]]) else rep(NA_integer_, nrow(sf_obj))
    if (!is.na(sc_col)) {
      sf_obj  <- sf_obj[!is.na(sc_vals) & sc_vals >= 2, ]
      sc_vals <- sc_vals[!is.na(sc_vals) & sc_vals >= 2]
    }

    sf_obj <- sf_obj[!st_is_empty(sf_obj), ]
    if (nrow(sf_obj) == 0) return(NULL)

    if (is.na(st_crs(sf_obj))) sf_obj <- st_set_crs(sf_obj, 4326)
    if (st_crs(sf_obj)$epsg != 25832)
      sf_obj <- st_transform(sf_obj, crs = 25832)
    coords <- st_coordinates(sf_obj)

    time_col <- find_time_col(sf_obj)
    if (is.na(time_col)) {
      cat(sprintf("  [WARN] Keine Zeitspalte in %s — übersprungen\n", basename(fp)))
      return(NULL)
    }
    if (time_col == "Date_Time_combined") {
      ts_raw <- as.POSIXct(paste(sf_obj$Date, sf_obj$Time), tz = "UTC")
    } else {
      ts_raw <- sf_obj[[time_col]]
      if (!inherits(ts_raw, "POSIXct"))
        ts_raw <- as.POSIXct(as.character(ts_raw), tz = "UTC")
    }

    sc_export <- if (!is.na(sc_col)) as.integer(sf_obj[[sc_col]]) else NA_integer_

    tibble::tibble(
      tier_label    = lbl,
      timestamp     = with_tz(ts_raw, "Europe/Berlin"),
      x_utm         = coords[, 1],
      y_utm         = coords[, 2],
      station_count = if (length(sc_export) == nrow(sf_obj)) sc_export else NA_integer_
    ) %>%
      filter(!is.na(x_utm), !is.na(y_utm), !is.na(timestamp)) %>%
      mutate(Date = as.Date(timestamp), Hour = hour(timestamp))

  }, error = function(e) {
    cat(sprintf("  [WARN] Fehler in %s: %s\n", basename(fp), conditionMessage(e)))
    NULL
  })
})

cat(sprintf("\n  Roh gesamt: %d Fixes von %d Tieren\n",
            nrow(alle_fixes), n_distinct(alle_fixes$tier_label)))

# ── 5. Metadaten joinen + Tage-seit-Auswilderung ──────────────────────────────
alle_fixes <- alle_fixes %>%
  left_join(meta %>% select(tier_label, date_release, date_last_signal,
                             sex, diagnosis, sr_dauer_tage),
            by = "tier_label")

# Fallback: Tiere ohne date_release in Metadaten → ersten Fix-Tag als Tag 0 verwenden
tiere_ohne_datum <- alle_fixes %>%
  filter(is.na(date_release)) %>%
  group_by(tier_label) %>%
  summarise(date_release_fallback = min(Date, na.rm = TRUE), .groups = "drop")

if (nrow(tiere_ohne_datum) > 0) {
  cat(sprintf("  [INFO] %d Tier(e) ohne date_release in Metadaten — erster Fix-Tag als Tag 0:\n",
              nrow(tiere_ohne_datum)))
  print(tiere_ohne_datum)
  alle_fixes <- alle_fixes %>%
    left_join(tiere_ohne_datum, by = "tier_label") %>%
    mutate(date_release = if_else(is.na(date_release), date_release_fallback, date_release),
           release_fallback = !is.na(date_release_fallback)) %>%
    select(-date_release_fallback)
} else {
  alle_fixes$release_fallback <- FALSE
}

alle_fixes <- alle_fixes %>%
  filter(!is.na(date_release)) %>%
  mutate(
    days_since_rel = as.numeric(Date - date_release),
    nach_release   = days_since_rel >= -2 & days_since_rel <= DAYS_MAX,
    vor_last_sig   = is.na(date_last_signal) | Date <= date_last_signal
  ) %>%
  filter(nach_release, vor_last_sig) %>%
  select(-nach_release, -vor_last_sig)

cat(sprintf("  Nach Zeitfilter: %d Fixes von %d Tieren\n",
            nrow(alle_fixes), n_distinct(alle_fixes$tier_label)))

# ── 6. Distanz zur Voliere ─────────────────────────────────────────────────────
alle_fixes <- alle_fixes %>%
  mutate(dist_voliere_m = sqrt((x_utm - vol_xy[1,1])^2 + (y_utm - vol_xy[1,2])^2))

cat(sprintf("\n  Distanz zur Voliere:\n"))
cat(sprintf("    Median: %.0f m\n",  median(alle_fixes$dist_voliere_m, na.rm = TRUE)))
cat(sprintf("    Mittel: %.0f m\n",  mean(alle_fixes$dist_voliere_m, na.rm = TRUE)))
cat(sprintf("    Max:    %.0f m\n",  max(alle_fixes$dist_voliere_m, na.rm = TRUE)))

# ── 7. Tägliche Mediandistanz pro Tier ────────────────────────────────────────
daily_dist <- alle_fixes %>%
  group_by(tier_label, sex, diagnosis, sr_dauer_tage, days_since_rel) %>%
  summarise(
    dist_median = median(dist_voliere_m, na.rm = TRUE),
    dist_q25    = quantile(dist_voliere_m, 0.25, na.rm = TRUE),
    dist_q75    = quantile(dist_voliere_m, 0.75, na.rm = TRUE),
    n_fixes     = n(),
    .groups     = "drop"
  )

mean_dist <- daily_dist %>%
  group_by(days_since_rel) %>%
  summarise(
    dist_mean = mean(dist_median, na.rm = TRUE),
    dist_se   = sd(dist_median, na.rm = TRUE) / sqrt(n()),
    n_tiere   = n(),
    .groups   = "drop"
  )

# ── 8. Theme ───────────────────────────────────────────────────────────────────
COL_IGEL   <- "#4C9A52"
COL_M      <- "#2166AC"
COL_F      <- "#D6604D"
MIN_N_PLOT <- 3L

theme_igel <- theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle    = element_text(size = 11, color = "grey45", margin = margin(b = 8)),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 11),
    axis.line        = element_line(color = "grey40"),
    panel.grid.major = element_line(color = "grey93", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 15, 8, 10)
  )

# ── 9. Plot A — Zeitverlauf Distanz ───────────────────────────────────────────
cat("\n── Erstelle Plots ───────────────────────────────────────────\n")

p6a <- ggplot() +
  annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = Inf,
           fill = COL_IGEL, alpha = 0.07) +
  annotate("text", x = 3.5, y = Inf, vjust = 1.8, size = 3.3,
           label = "Erste Woche", color = COL_IGEL, fontface = "italic") +
  geom_hline(yintercept = VOLIERE_RADIUS_M,
             linetype = "dashed", color = "grey60", linewidth = 0.7) +
  annotate("label", x = DAYS_MAX, y = VOLIERE_RADIUS_M,
           label = paste0("Volieren-Radius (~", VOLIERE_RADIUS_M, " m)"),
           hjust = 1, size = 3.2, color = "grey50",
           fill = "white", label.padding = unit(2, "pt")) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.8) +
  annotate("text", x = 0.3, y = Inf, vjust = 1.8, hjust = 0,
           label = "Auswilderung", size = 3.3, color = "grey35", fontface = "italic") +
  geom_line(data = daily_dist,
            aes(x = days_since_rel, y = dist_median, group = tier_label, color = tier_label),
            alpha = 0.35, linewidth = 0.7) +
  geom_ribbon(data = mean_dist %>% filter(n_tiere >= MIN_N_PLOT),
              aes(x = days_since_rel,
                  ymin = pmax(0, dist_mean - dist_se),
                  ymax = dist_mean + dist_se),
              fill = COL_IGEL, alpha = 0.30, color = NA) +
  geom_line(data = mean_dist %>% filter(n_tiere >= MIN_N_PLOT),
            aes(x = days_since_rel, y = dist_mean),
            color = COL_IGEL, linewidth = 2.2, lineend = "round") +
  geom_point(data = mean_dist %>% filter(n_tiere >= MIN_N_PLOT),
             aes(x = days_since_rel, y = dist_mean),
             color = COL_IGEL, size = 2, shape = 16) +
  { if (any(mean_dist$n_tiere > 0 & mean_dist$n_tiere < MIN_N_PLOT))
      geom_line(data = mean_dist %>% filter(n_tiere > 0 & n_tiere < MIN_N_PLOT),
                aes(x = days_since_rel, y = dist_mean),
                color = "grey70", linewidth = 1.2, linetype = "dotted") } +
  scale_color_viridis_d(option = "D", end = 0.88, guide = "none") +
  scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                     expand = expansion(mult = c(0.04, 0.06))) +
  scale_x_continuous(breaks = seq(0, DAYS_MAX, by = 5),
                     expand = expansion(add = c(0.5, 0.5))) +
  labs(title    = "Distanz zur Auswilderungsvoliere nach Auswilderung",
       subtitle = sprintf("n = %d Igel  |  Dicke Linie = Gruppenmedian ± SE  |  Dünne Linien = Einzeltiere  |  Gestrichelt = n < %d",
                          n_distinct(daily_dist$tier_label), MIN_N_PLOT),
       x = "Tage nach Auswilderung", y = "Distanz zur Voliere (m)") +
  theme_igel

ggsave(file.path(out_ordner, "06a_distanz_voliere_zeitverlauf.png"), p6a,
       width = 12, height = 6, dpi = 200, bg = "white")
cat("  ✓ 06a_distanz_voliere_zeitverlauf.png\n")

# ── 10. Plot B — Presentation-Version (Englisch) ──────────────────────────────
n_annot <- mean_dist %>%
  filter(n_tiere > 0, days_since_rel >= 0,
         days_since_rel %% 5 == 0 | days_since_rel == 1) %>%
  mutate(label = paste0("n=", n_tiere))

x_hi <- max(mean_dist %>% filter(n_tiere >= MIN_N_PLOT) %>% pull(days_since_rel), na.rm = TRUE) + 1

p6b <- ggplot() +
  annotate("rect", xmin = 0, xmax = min(7, x_hi), ymin = -Inf, ymax = Inf,
           fill = COL_IGEL, alpha = 0.07) +
  annotate("text", x = 3.5, y = Inf, vjust = 1.8, hjust = 0.5,
           label = "Early\nphase", size = 3.3, color = "grey50", fontface = "italic") +
  geom_hline(yintercept = VOLIERE_RADIUS_M,
             linetype = "dashed", color = "grey65", linewidth = 0.7) +
  annotate("label", x = x_hi, y = VOLIERE_RADIUS_M,
           label = paste0("Aviary boundary (~", VOLIERE_RADIUS_M, " m)"),
           hjust = 1, size = 3.2, color = "grey50",
           fill = "white", label.padding = unit(2, "pt")) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.8) +
  annotate("text", x = 0.3, y = Inf, vjust = 1.8, hjust = 0,
           label = "Release", size = 3.5, color = "grey35", fontface = "italic") +
  geom_line(data = daily_dist %>% filter(days_since_rel >= -1, days_since_rel <= x_hi),
            aes(x = days_since_rel, y = dist_median, group = tier_label),
            color = COL_IGEL, alpha = 0.15, linewidth = 0.6) +
  geom_ribbon(data = mean_dist %>% filter(n_tiere >= MIN_N_PLOT, days_since_rel <= x_hi),
              aes(x = days_since_rel,
                  ymin = pmax(0, dist_mean - dist_se),
                  ymax = dist_mean + dist_se),
              fill = COL_IGEL, alpha = 0.28, color = NA) +
  geom_line(data = mean_dist %>% filter(n_tiere >= MIN_N_PLOT, days_since_rel <= x_hi),
            aes(x = days_since_rel, y = dist_mean),
            color = COL_IGEL, linewidth = 2.3, lineend = "round") +
  geom_point(data = mean_dist %>% filter(n_tiere >= MIN_N_PLOT, days_since_rel <= x_hi),
             aes(x = days_since_rel, y = dist_mean),
             color = COL_IGEL, size = 2.5, shape = 16) +
  { if (any(mean_dist$n_tiere > 0 & mean_dist$n_tiere < MIN_N_PLOT))
      geom_line(data = mean_dist %>% filter(n_tiere > 0 & n_tiere < MIN_N_PLOT, days_since_rel <= x_hi),
                aes(x = days_since_rel, y = dist_mean),
                color = "grey75", linewidth = 1.2, linetype = "dotted") } +
  geom_text(data = n_annot %>% filter(days_since_rel <= x_hi),
            aes(x = days_since_rel, y = 0, label = label),
            vjust = 1.8, size = 3.0, color = "grey55") +
  scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                     expand = expansion(mult = c(0.08, 0.05))) +
  scale_x_continuous(breaks = seq(0, x_hi, by = 5),
                     expand = expansion(add = c(0.5, 0.8))) +
  labs(title    = sprintf("Post-release dispersal from the release aviary  (n = %d hedgehogs)",
                           n_distinct(daily_dist$tier_label)),
       subtitle = sprintf("Bold line = group mean ± SE  |  Faint lines = individual animals  |  Dotted grey = n < %d animals",
                           MIN_N_PLOT),
       x = "Days since release", y = "Distance to release aviary (m)") +
  theme_classic(base_size = 14) +
  theme(plot.title       = element_text(face = "bold", size = 16, hjust = 0),
        plot.subtitle    = element_text(size = 12, color = "grey45", margin = margin(b = 10)),
        axis.title       = element_text(size = 13),
        axis.text        = element_text(size = 12),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(12, 18, 8, 12))

ggsave(file.path(out_ordner, "06b_distanz_voliere_pres.png"), p6b,
       width = 11, height = 5.5, dpi = 200, bg = "white")
cat("  ✓ 06b_distanz_voliere_pres.png\n")

# ── 10b. Plot — PUBLICATION VERSION (LME-modelled trend + spread) ─────────────
# Population trend from the linear mixed model (dist ~ day + (1+day|individual)),
# 95% CI from fixed-effects covariance; light band = interquartile spread across
# individuals; the modelled line is solid while >=5 animals contributed and
# dashed/greyed below that (survivorship). N-at-risk shown along the lower axis.
cat("== 10b. Publication dispersal figure ==\n")

lme_data_pub <- daily_dist %>% filter(days_since_rel >= 0, !is.na(dist_median))

if (n_distinct(lme_data_pub$tier_label) >= 3 && nrow(lme_data_pub) >= 10) {

  spread_pub <- lme_data_pub %>%
    group_by(days_since_rel) %>%
    summarise(q25 = quantile(dist_median, 0.25, na.rm = TRUE),
              q75 = quantile(dist_median, 0.75, na.rm = TRUE),
              n   = dplyr::n(), .groups = "drop")

  N_SOLID  <- 5L
  x_hi_pub <- max(spread_pub$days_since_rel[spread_pub$n >= MIN_N_PLOT], na.rm = TRUE)

  # LME population-level prediction + 95% CI (fixed effects only)
  m_pub <- tryCatch(
    lme4::lmer(dist_median ~ days_since_rel + (1 + days_since_rel | tier_label),
               data = lme_data_pub, REML = TRUE,
               control = lme4::lmerControl(optimizer = "bobyqa")),
    error = function(e)
      lme4::lmer(dist_median ~ days_since_rel + (1 | tier_label),
                 data = lme_data_pub, REML = TRUE))

  nd      <- data.frame(days_since_rel = seq(0, x_hi_pub, by = 1))
  mm      <- model.matrix(~ days_since_rel, nd)
  nd$fit  <- as.numeric(mm %*% lme4::fixef(m_pub))
  vc      <- as.matrix(vcov(m_pub))
  nd$se   <- sqrt(rowSums((mm %*% vc) * mm))
  nd$lo   <- nd$fit - 1.96 * nd$se
  nd$hi   <- nd$fit + 1.96 * nd$se
  nd      <- merge(nd, spread_pub[, c("days_since_rel", "n")],
                   by = "days_since_rel", all.x = TRUE)
  nd_solid <- nd[nd$n >= N_SOLID, ]
  nd_faint <- nd[nd$n <  N_SOLID, ]
  beta_day <- as.numeric(lme4::fixef(m_pub)["days_since_rel"])

  iqr_pub <- spread_pub %>% filter(n >= MIN_N_PLOT, days_since_rel <= x_hi_pub)
  nrisk   <- spread_pub %>% filter(days_since_rel %in% seq(0, x_hi_pub, by = 5))

  col_line <- "#2E7D32"; col_ci <- "#66BB6A"; col_iqr <- "#BCDFBC"

  # y-Achse ab 0; Oberkante am Datenbereich. N-at-risk-Zeile in der unteren Bandzone.
  q_hi    <- max(iqr_pub$q75, na.rm = TRUE)
  y_top   <- ceiling((q_hi + 15) / 10) * 10
  y_nrisk <- y_top * 0.06

  p6_pub <- ggplot() +
    annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = Inf, fill = "#E8F3E8") +
    annotate("text", x = 0.4, y = y_top * 0.985,
             label = "Supplementary feeding\n(days 0–7)",
             hjust = 0, vjust = 1, size = 3, colour = "#5C8A5C",
             fontface = "italic", lineheight = 0.95) +
    geom_ribbon(data = iqr_pub, aes(x = days_since_rel, ymin = q25, ymax = q75),
                fill = col_iqr, alpha = 0.7) +
    geom_ribbon(data = nd, aes(x = days_since_rel, ymin = lo, ymax = hi),
                fill = col_ci, alpha = 0.55) +
    geom_line(data = nd_solid, aes(x = days_since_rel, y = fit),
              colour = col_line, linewidth = 1.7) +
    { if (nrow(nd_faint) > 0)
        geom_line(data = rbind(tail(nd_solid, 1), nd_faint),
                  aes(x = days_since_rel, y = fit),
                  colour = col_line, linewidth = 1.2, linetype = "22", alpha = 0.55) } +
    geom_hline(yintercept = VOLIERE_RADIUS_M, colour = "grey60",
               linetype = "21", linewidth = 0.5) +
    annotate("text", x = x_hi_pub, y = VOLIERE_RADIUS_M + 7,
             label = sprintf("Enclosure vicinity (~%d m)", VOLIERE_RADIUS_M),
             hjust = 1, vjust = 0, size = 2.9, colour = "grey55") +
    geom_text(data = nrisk, aes(x = days_since_rel, y = y_nrisk, label = n),
              size = 3, colour = "grey55") +
    annotate("text", x = -1.3, y = y_nrisk, label = "N =", hjust = 1, size = 3,
             colour = "grey55", fontface = "bold") +
    scale_x_continuous(breaks = seq(0, x_hi_pub, 5)) +
    scale_y_continuous(breaks = seq(0, y_top, 50), labels = label_number(suffix = " m")) +
    coord_cartesian(xlim = c(-2, x_hi_pub + 0.5),
                    ylim = c(0, y_top), clip = "off") +
    labs(x = "Days since release", y = "Distance to release site") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title    = element_text(size = 12),
          axis.text     = element_text(size = 11),
          plot.margin   = margin(10, 16, 8, 12))

  ggsave(file.path(out_ordner, "06b_distanz_voliere_publication.png"), p6_pub,
         width = 11, height = 5.4, dpi = 300, bg = "white")
  ggsave(file.path(out_ordner, "06b_distanz_voliere_publication.pdf"), p6_pub,
         width = 11, height = 5.4, bg = "white")
  cat("  ✓ 06b_distanz_voliere_publication.png/.pdf\n")
} else {
  cat("  Zu wenig Daten fuer Publikationsabbildung\n")
}


# ── 11. Plot C — Vergleich Männchen vs. Weibchen ──────────────────────────────
if (!all(is.na(daily_dist$sex))) {

  daily_sex <- daily_dist %>%
    filter(!is.na(sex), tolower(sex) %in% c("male","female","männlich","weiblich","m","f","w")) %>%
    mutate(sex_de = case_when(
      tolower(sex) %in% c("male", "m", "männlich")   ~ "Männchen",
      tolower(sex) %in% c("female", "f", "w", "weiblich") ~ "Weibchen",
      TRUE ~ sex
    ))

  if (n_distinct(daily_sex$sex_de) >= 2) {
    mean_sex <- daily_sex %>%
      group_by(sex_de, days_since_rel) %>%
      summarise(dist_mean = mean(dist_median, na.rm = TRUE),
                dist_se   = sd(dist_median, na.rm = TRUE) / sqrt(n()),
                n_tiere   = n(), .groups = "drop")

    sex_farben <- c("Männchen" = COL_M, "Weibchen" = COL_F)
    n_m <- n_distinct(daily_sex %>% filter(sex_de == "Männchen") %>% pull(tier_label))
    n_w <- n_distinct(daily_sex %>% filter(sex_de == "Weibchen") %>% pull(tier_label))

    p6c <- ggplot(mean_sex %>% filter(n_tiere >= MIN_N_PLOT),
                  aes(x = days_since_rel, y = dist_mean, color = sex_de, fill = sex_de)) +
      annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = Inf,
               fill = "grey90", alpha = 0.5) +
      geom_hline(yintercept = VOLIERE_RADIUS_M,
                 linetype = "dashed", color = "grey60", linewidth = 0.7) +
      geom_vline(xintercept = 0, color = "grey35", linewidth = 0.8) +
      annotate("text", x = 0.3, y = Inf, vjust = 1.8, hjust = 0,
               label = "Auswilderung", size = 3.2, color = "grey35", fontface = "italic") +
      geom_ribbon(aes(ymin = pmax(0, dist_mean - dist_se), ymax = dist_mean + dist_se),
                  alpha = 0.20, color = NA) +
      geom_line(linewidth = 2, lineend = "round") +
      geom_point(size = 2.5, shape = 16) +
      scale_color_manual(values = sex_farben,
                         labels = c("Männchen" = sprintf("Männchen (n=%d)", n_m),
                                    "Weibchen" = sprintf("Weibchen (n=%d)", n_w))) +
      scale_fill_manual(values = sex_farben, guide = "none") +
      scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                         expand = expansion(mult = c(0.04, 0.06))) +
      scale_x_continuous(breaks = seq(0, DAYS_MAX, by = 5)) +
      labs(title    = "Distanz zur Voliere — Männchen vs. Weibchen",
           subtitle = "Mittelwert ± SE | Nur Tage mit ≥ 3 Tieren",
           x = "Tage nach Auswilderung", y = "Distanz zur Voliere (m)", color = NULL) +
      theme_igel + theme(legend.position = "bottom")

    ggsave(file.path(out_ordner, "06c_distanz_sex_vergleich.png"), p6c,
           width = 11, height = 6, dpi = 200, bg = "white")
    cat("  ✓ 06c_distanz_sex_vergleich.png\n")
  } else {
    cat("  [INFO] Nicht genug Geschlechtskategorien für Vergleichsplot.\n")
  }
}

# ── 12. Plot D — Boxplot: Medianabstand pro Tier ──────────────────────────────
animal_summary <- daily_dist %>%
  filter(days_since_rel >= 0) %>%
  group_by(tier_label, sex) %>%
  summarise(median_dist = median(dist_median, na.rm = TRUE),
            max_dist    = max(dist_median, na.rm = TRUE),
            n_tage      = n(), .groups = "drop") %>%
  mutate(sex_de = case_when(
    tolower(sex) %in% c("male", "m", "männlich")        ~ "Männchen",
    tolower(sex) %in% c("female", "f", "w", "weiblich") ~ "Weibchen",
    TRUE ~ "Unbekannt"
  ))

p6d <- ggplot(animal_summary, aes(x = 1, y = median_dist)) +
  geom_hline(yintercept = VOLIERE_RADIUS_M,
             linetype = "dashed", color = "grey60", linewidth = 0.8) +
  annotate("label", x = 1, y = VOLIERE_RADIUS_M,
           label = paste0("Volieren-Radius (~", VOLIERE_RADIUS_M, " m)"),
           size = 3.2, color = "grey50", fill = "white", label.padding = unit(2, "pt")) +
  geom_boxplot(fill = COL_IGEL, alpha = 0.65, outlier.shape = NA,
               width = 0.45, color = "grey30", linewidth = 0.6) +
  geom_jitter(width = 0.12, size = 3, alpha = 0.75, color = "grey25", shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5.5, color = "black") +
  scale_x_continuous(limits = c(0.5, 1.5), breaks = NULL, labels = NULL) +
  scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                     expand = expansion(mult = c(0.03, 0.10))) +
  labs(title    = "Medianabstand zur Auswilderungsvoliere",
       subtitle = sprintf("Pro Tier: Median der täglichen Mediandistanzen  |  n = %d Igel  |  ◆ = Gesamtmittel",
                          nrow(animal_summary)),
       x = NULL, y = "Median-Distanz zur Voliere (m)") +
  theme_igel

ggsave(file.path(out_ordner, "06d_distanz_boxplot.png"), p6d,
       width = 6, height = 6, dpi = 200, bg = "white")
cat("  ✓ 06d_distanz_boxplot.png\n")

# ── 12b. Plot E — GAMM-Smoother ───────────────────────────────────────────────
cat("\n── Plot E: GAMM-Smoother ────────────────────────────────────\n")

if (requireNamespace("mgcv", quietly = TRUE)) {
  gam_data <- daily_dist %>%
    filter(days_since_rel >= 0, days_since_rel <= DAYS_MAX, !is.na(dist_median))

  p6e <- ggplot(gam_data, aes(x = days_since_rel, y = dist_median)) +
    annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = Inf,
             fill = COL_IGEL, alpha = 0.06) +
    geom_vline(xintercept = 0, color = "grey35", linewidth = 0.8) +
    annotate("text", x = 0.3, y = Inf, vjust = 1.6, hjust = 0,
             label = "Auswilderung", size = 3.2, color = "grey35", fontface = "italic") +
    annotate("text", x = 3.5, y = Inf, vjust = 3.5, size = 3.2,
             label = "Erste Woche", color = COL_IGEL, fontface = "italic") +
    geom_hline(yintercept = VOLIERE_RADIUS_M,
               linetype = "dashed", color = "grey65", linewidth = 0.7) +
    annotate("label", x = DAYS_MAX, y = VOLIERE_RADIUS_M,
             label = paste0("Volieren-Radius (~", VOLIERE_RADIUS_M, " m)"),
             hjust = 1, size = 3.1, color = "grey50",
             fill = "white", label.padding = unit(2, "pt")) +
    geom_line(aes(group = tier_label), color = COL_IGEL, alpha = 0.12, linewidth = 0.55) +
    stat_summary(fun = mean, geom = "point",
                 color = "grey60", size = 1.8, alpha = 0.6, shape = 16) +
    stat_smooth(method = "gam", formula = y ~ s(x, k = 6),
                color = COL_IGEL, fill = COL_IGEL, alpha = 0.22,
                linewidth = 2.2, se = TRUE) +
    scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                       expand = expansion(mult = c(0.04, 0.06))) +
    scale_x_continuous(breaks = seq(0, DAYS_MAX, by = 5),
                       expand = expansion(add = c(0.3, 0.5))) +
    labs(title    = "Distanz zur Auswilderungsvoliere — GAM-Smoother",
         subtitle = sprintf("n = %d Igel  |  Linie = GAM-Smooth (k=6) mit 95%%–CI  |  Punkte = Tages-Gruppenmittel  |  Dünne Linien = Einzeltiere",
                             n_distinct(gam_data$tier_label)),
         x = "Tage nach Auswilderung", y = "Distanz zur Voliere (m)") +
    theme_igel

  ggsave(file.path(out_ordner, "06e_distanz_gam_smoother.png"), p6e,
         width = 12, height = 6, dpi = 200, bg = "white")
  cat("  ✓ 06e_distanz_gam_smoother.png\n")

  p6e_en <- p6e +
    labs(title    = sprintf("Post-release dispersal — GAM smoother  (n = %d hedgehogs)",
                             n_distinct(gam_data$tier_label)),
         subtitle = "Line = GAM smooth (k=6) with 95% CI  |  Dots = daily group mean  |  Faint lines = individual animals",
         x = "Days since release", y = "Distance to release aviary (m)") +
    annotate("text", x = 0.3, y = Inf, vjust = 1.6, hjust = 0,
             label = "Release", size = 3.2, color = "grey35", fontface = "italic") +
    annotate("text", x = 3.5, y = Inf, vjust = 3.5, size = 3.2,
             label = "First week", color = COL_IGEL, fontface = "italic") +
    theme_classic(base_size = 14) +
    theme(plot.title       = element_text(face = "bold", size = 16, hjust = 0),
          plot.subtitle    = element_text(size = 12, color = "grey45", margin = margin(b = 10)),
          axis.title       = element_text(size = 13),
          axis.text        = element_text(size = 12),
          panel.grid.major = element_line(color = "grey93", linewidth = 0.4),
          panel.grid.minor = element_blank(),
          plot.margin      = margin(12, 18, 8, 12))

  ggsave(file.path(out_ordner, "06e_distanz_gam_pres.png"), p6e_en,
         width = 11, height = 5.5, dpi = 200, bg = "white")
  cat("  ✓ 06e_distanz_gam_pres.png\n")
} else {
  cat("  [INFO] mgcv nicht verfügbar — Plot E übersprungen.\n")
}

# ── 12c. Plot F — Raincloud-Plot ──────────────────────────────────────────────
cat("\n── Plot F: Raincloud-Plot ───────────────────────────────────\n")

if (requireNamespace("ggdist", quietly = TRUE)) {
  library(ggdist)

  p6f <- ggplot(animal_summary, aes(x = 1, y = median_dist)) +
    geom_hline(yintercept = VOLIERE_RADIUS_M,
               linetype = "dashed", color = "grey65", linewidth = 0.8) +
    annotate("label", x = 1.45, y = VOLIERE_RADIUS_M,
             label = paste0("~", VOLIERE_RADIUS_M, " m"), hjust = 1,
             size = 3.0, color = "grey55", fill = "white", label.padding = unit(1.5, "pt")) +
    stat_halfeye(adjust = 0.8, width = 0.45, .width = 0, point_colour = NA,
                 fill = COL_IGEL, alpha = 0.55, justification = -0.25) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                 color = "grey30", linewidth = 0.7, alpha = 0.8) +
    geom_jitter(width = 0.06, size = 3.2, alpha = 0.80, color = COL_IGEL, shape = 16) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 6, color = "black") +
    scale_x_continuous(limits = c(0.5, 1.7), breaks = NULL, labels = NULL) +
    scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                       expand = expansion(mult = c(0.03, 0.08))) +
    labs(title    = "Medianabstand zur Auswilderungsvoliere",
         subtitle = sprintf("Pro Tier: Median aller Tages-Mediandistanzen  |  n = %d Igel  |  ◆ = Gesamtmittel  |  Kurve = Dichte",
                             nrow(animal_summary)),
         x = NULL, y = "Median-Distanz zur Voliere (m)") +
    theme_igel + theme(axis.ticks.x = element_blank())

  ggsave(file.path(out_ordner, "06f_distanz_raincloud.png"), p6f,
         width = 6, height = 7, dpi = 200, bg = "white")
  cat("  ✓ 06f_distanz_raincloud.png\n")
} else {
  cat("  [INFO] ggdist nicht verfügbar — Plot F übersprungen.\n")
}

# ── 12d. Plot G — Raincloud nach Diagnose-Kategorie ───────────────────────────
cat("\n── Plot G: Diagnose-Vergleich ───────────────────────────────\n")

animal_diag <- animal_summary %>%
  left_join(meta %>% select(tier_label, diagnosis), by = "tier_label") %>%
  mutate(
    diag_kat = case_when(
      grepl("blind",       diagnosis, ignore.case = TRUE)                            ~ "Blindheit",
      grepl("orph",        diagnosis, ignore.case = TRUE)                            ~ "Waisling",
      grepl("endoparasit", diagnosis, ignore.case = TRUE) &
        !grepl("blind|orphan|orph", diagnosis, ignore.case = TRUE)                  ~ "Endoparasiten",
      grepl("parasite",    diagnosis, ignore.case = TRUE) &
        !grepl("blind|orphan|orph|endo", diagnosis, ignore.case = TRUE)             ~ "Parasiten",
      grepl("trauma",      diagnosis, ignore.case = TRUE) &
        !grepl("blind",    diagnosis, ignore.case = TRUE)                            ~ "Trauma",
      grepl("fungal",      diagnosis, ignore.case = TRUE)                            ~ "Pilzerkrankung",
      TRUE                                                                            ~ "Sonstige"
    ),
    diag_kat = factor(diag_kat, levels = c("Waisling","Parasiten","Endoparasiten",
                                            "Trauma","Pilzerkrankung","Sonstige","Blindheit"))
  )

diag_farben <- c("Waisling"       = "#4393C3",
                  "Parasiten"      = "#74C476",
                  "Endoparasiten"  = "#31A354",
                  "Trauma"         = "#F4A582",
                  "Pilzerkrankung" = "#DFC27D",
                  "Sonstige"       = "#B2ABD2",
                  "Blindheit"      = "#D6604D")

n_diag <- animal_diag %>% count(diag_kat) %>%
  mutate(label = paste0(diag_kat, "\n(n=", n, ")"))
diag_labels <- setNames(n_diag$label, n_diag$diag_kat)

if (requireNamespace("ggdist", quietly = TRUE)) {
  p6g <- ggplot(animal_diag,
                aes(x = diag_kat, y = median_dist, fill = diag_kat, color = diag_kat)) +
    geom_hline(yintercept = VOLIERE_RADIUS_M,
               linetype = "dashed", color = "grey65", linewidth = 0.7) +
    annotate("label", x = 0.5, y = VOLIERE_RADIUS_M,
             label = paste0("~", VOLIERE_RADIUS_M, " m"),
             hjust = 0, size = 3.0, color = "grey55",
             fill = "white", label.padding = unit(1.5, "pt")) +
    stat_halfeye(adjust = 0.9, width = 0.45, .width = 0, point_colour = NA,
                 alpha = 0.40, justification = -0.2) +
    geom_boxplot(width = 0.18, outlier.shape = NA, color = "grey25",
                 linewidth = 0.6, alpha = 0.0, fill = "white") +
    geom_jitter(width = 0.08, size = 3, alpha = 0.85, shape = 16) +
    # ggrepel: Guard gegen fehlende Installation
    { if (requireNamespace("ggrepel", quietly = TRUE))
        ggrepel::geom_text_repel(
          aes(label = tier_label), size = 2.8, color = "grey30",
          box.padding = 0.3, point.padding = 0.2,
          max.overlaps = 8, segment.color = "grey70", segment.size = 0.3,
          show.legend = FALSE)
      else
        geom_text(aes(label = tier_label), vjust = -1.0, size = 2.8,
                  color = "grey30", show.legend = FALSE) } +
    stat_summary(fun = mean, geom = "point",
                 shape = 18, size = 5, color = "black", show.legend = FALSE) +
    scale_fill_manual(values  = diag_farben, guide = "none") +
    scale_color_manual(values = diag_farben, guide = "none") +
    scale_x_discrete(labels = diag_labels) +
    scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                       expand = expansion(mult = c(0.03, 0.10))) +
    labs(title    = "Medianabstand zur Voliere nach Diagnose-Kategorie",
         subtitle = sprintf("n = %d Igel  |  ◆ = Kategoriemittel  |  Kurve = Dichte  |  Rot-Orange = Blindheit",
                             nrow(animal_diag)),
         x = NULL, y = "Median-Distanz zur Voliere (m)") +
    theme_igel + theme(axis.text.x = element_text(size = 10, lineheight = 1.2))

  ggsave(file.path(out_ordner, "06g_distanz_diagnose.png"), p6g,
         width = 11, height = 7, dpi = 200, bg = "white")
  cat("  ✓ 06g_distanz_diagnose.png\n")
} else {
  p6g_simple <- ggplot(animal_diag,
                        aes(x = diag_kat, y = median_dist, fill = diag_kat, color = diag_kat)) +
    geom_hline(yintercept = VOLIERE_RADIUS_M,
               linetype = "dashed", color = "grey65", linewidth = 0.7) +
    geom_boxplot(alpha = 0.55, outlier.shape = NA,
                 color = "grey30", linewidth = 0.6, width = 0.5) +
    geom_jitter(width = 0.12, size = 3, alpha = 0.85, shape = 16) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "black") +
    scale_fill_manual(values  = diag_farben, guide = "none") +
    scale_color_manual(values = diag_farben, guide = "none") +
    scale_x_discrete(labels = diag_labels) +
    scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                       expand = expansion(mult = c(0.03, 0.10))) +
    labs(title    = "Medianabstand zur Voliere nach Diagnose-Kategorie",
         subtitle = sprintf("n = %d Igel  |  ◆ = Kategoriemittel", nrow(animal_diag)),
         x = NULL, y = "Median-Distanz zur Voliere (m)") +
    theme_igel + theme(axis.text.x = element_text(size = 10, lineheight = 1.2))

  ggsave(file.path(out_ordner, "06g_distanz_diagnose.png"), p6g_simple,
         width = 11, height = 6, dpi = 200, bg = "white")
  cat("  ✓ 06g_distanz_diagnose.png (Fallback ohne ggdist)\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# ── Persistenz-Analyse ────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Persistenz-Analyse: Wann verlassen Igel das Gelände? ───\n")

ABGANG_SCHWELLE_M <- 200L
ABGANG_MIN_TAGE   <- 5L
RANGE_SCHWELLE_M  <- 250L

endpunkte <- daily_dist %>%
  filter(days_since_rel >= 0) %>%
  group_by(tier_label, sex) %>%
  summarise(letzter_tag     = max(days_since_rel),
            enddistanz      = dist_median[which.max(days_since_rel)],
            max_dist_gesamt = max(dist_median, na.rm = TRUE),
            .groups = "drop")

slope_letzte5_df <- daily_dist %>%
  filter(days_since_rel >= 0) %>%
  left_join(endpunkte %>% select(tier_label, letzter_tag), by = "tier_label") %>%
  filter(days_since_rel >= letzter_tag - 4) %>%
  group_by(tier_label) %>%
  summarise(slope_letzte5 = if (n() >= 2)
              coef(lm(dist_median ~ days_since_rel))[["days_since_rel"]]
            else NA_real_, .groups = "drop")

# ── tag_found aus Metadaten einmergen ─────────────────────────────────────────
endpunkte <- endpunkte %>%
  left_join(slope_letzte5_df, by = "tier_label") %>%
  left_join(meta %>% select(tier_label, tag_found, longterm_outcome), by = "tier_label")

# ══════════════════════════════════════════════════════════════════════════════
# NEUE KLASSIFIKATIONSLOGIK (basierend auf sub-minutlicher Auflösung)
# ══════════════════════════════════════════════════════════════════════════════
# Begründung:
#   Das tRackIT-System liefert sub-minutliche Positionsdaten. Eine Abwanderung
#   über den 250m-Detektionsbereich würde als graduelle Distanzzunahme sichtbar
#   sein (Tier läuft ~0.5 m/s → 300m Strecke = ~10 min = Hunderte Fixes).
#   Der Signalverlust bei Enddistanz <200m ist daher KEIN Abwanderungsindikator,
#   sondern weist auf plötzlichen Senderverlust (Detachment) oder Tod am Ort hin.
#
# Kategorien:
#   A: "Sender auf Gelände gefunden"
#      → tag_found = 1: Sender physisch auf dem Gelände wiedergefunden.
#        Beweist dass der Sender vor Ort abgefallen ist. Tier-Verbleib danach unklar.
#
#   B: "Aktive Abwanderung"
#      → tag_found = 0 UND Enddistanz ≥200m UND ≥5 Tage Tracking
#        Graduelle Bewegung zur Detektionsgrenze erkennbar. Sender nicht auf Gelände
#        → Tier hat den 250m-Bereich aktiv verlassen.
#
#   C: "Plötzlicher Signalverlust nahe Zentrum"
#      → tag_found = 0 UND Enddistanz <200m UND ≥5 Tage Tracking
#        Kein gradueller Distanzanstieg vor Signalverlust. Sender nicht auf Gelände
#        gefunden (ggf. in Vegetation, von Tier mitgenommen, oder nicht gesucht).
#
#   D: "Kurz-Tracking, nicht klassifizierbar"
#      → tag_found = 0 UND <5 Tage Tracking
#        Zu kurze Trackingdauer für gesicherte Aussage.

endpunkte <- endpunkte %>%
  mutate(
    abgang_typ = case_when(
      # Kategorie A: Sender physisch auf Gelände gefunden
      !is.na(tag_found) & tag_found == 1
        ~ "A: Sender auf Gelände gefunden",
      # Kategorie B: Aktive Abwanderung (graduelle Bewegung zur Grenze)
      (is.na(tag_found) | tag_found == 0) &
        letzter_tag >= ABGANG_MIN_TAGE &
        enddistanz  >= ABGANG_SCHWELLE_M
        ~ "B: Aktive Abwanderung\n(Enddist. ≥200 m, kein Sender)",
      # Kategorie C: Multilateralisation verloren (letzte Position nahe Zentrum)
      # Begruendung: tRackIT unterscheidet SIGNAL (Activity-CSV, eine Antenne reicht)
      # und POSITION (GPKG, braucht >=3 Antennen gleichzeitig). Signalverlust in
      # GPKG bedeutet NICHT zwingend dass das Tier weg ist — Antennengeometrie oder
      # Vegetation koennen Multilateralisation verhindern waehrend Signal weitergeht.
      # Activity-CSV zeigt bei manchen Tieren noch Aktivitaet nach letztem GPKG-Fix.
      # Sender nicht auf Gelände gefunden. Tier-Verbleib: UNKLAR.
      (is.na(tag_found) | tag_found == 0) &
        letzter_tag >= ABGANG_MIN_TAGE &
        enddistanz  <  ABGANG_SCHWELLE_M
        ~ "C: Multilateral. verloren\n(Pos. <200 m, Verbleib unklar)",
      # Kategorie D: Kurz-Tracking
      TRUE
        ~ "D: Kurz-Tracking\n(<5 Tage, kein Sender)"
    ),
    abgang_typ = factor(abgang_typ, levels = c(
      "A: Sender auf Gelände gefunden",
      "B: Aktive Abwanderung\n(Enddist. ≥200 m, kein Sender)",
      "C: Multilateral. verloren\n(Pos. <200 m, Verbleib unklar)",
      "D: Kurz-Tracking\n(<5 Tage, kein Sender)"))
  )

cat("\n  ═══════════════════════════════════════════════════════════\n")
cat("  Abgangstyp-Klassifikation — revidierte Logik\n")
cat("  ═══════════════════════════════════════════════════════════\n")
cat("  Grundlage: Sub-minutliche Positionsdaten + Geländesuche\n")
cat("  A: Sender physisch auf Gelände gefunden\n")
cat("     → Detachment bestätigt; Tier-Verbleib nach Detachment unklar\n")
cat("  B: Gelände verlassen, letzte Pos ≥200m von Voliere\n")
cat("     → Tier war bereits an Geländegrenze, dann Signal weg\n")
cat("  C: Multilateral. verloren, letzte Pos <200m — Verbleib unklar\n")
cat("     → Gelände abgesucht, Sender nicht gefunden = Sender weg\n")
cat("     → Activity-CSV zeigt teils weiteres Signal; Sender nicht auffindbar\n")
cat("  D: <5 Tage Tracking, kein Sender = nicht klassifizierbar\n")
cat("  ═══════════════════════════════════════════════════════════\n\n")
print(endpunkte %>%
        select(tier_label, letzter_tag, enddistanz, slope_letzte5,
               tag_found, longterm_outcome, abgang_typ) %>%
        arrange(abgang_typ, letzter_tag) %>% as.data.frame())

cat("\n  Zusammenfassung:\n")
print(table(endpunkte$abgang_typ))
cat("\n  Bestätigte Outcomes:\n")
print(endpunkte %>% filter(!is.na(longterm_outcome) & longterm_outcome != "unknown") %>%
        select(tier_label, abgang_typ, longterm_outcome) %>% as.data.frame())

abgang_farben <- c(
  "A: Sender auf Gelände gefunden"                         = "#4393C3",
  "B: Aktive Abwanderung\n(Enddist. ≥200 m, kein Sender)" = "#D6604D",
  "C: Multilateral. verloren\n(Pos. <200 m, Verbleib unklar)"= "#F4A582",
  "D: Kurz-Tracking\n(<5 Tage, kein Sender)"               = "#BDBDBD"
)

# ── Plot H — Retention-Kurve ──────────────────────────────────────────────────
retention <- mean_dist %>%
  filter(days_since_rel >= 0) %>%
  select(days_since_rel, n_tiere)
n_tag0 <- retention %>% filter(days_since_rel == 0) %>% pull(n_tiere)
retention <- retention %>%
  mutate(pct_aktiv = n_tiere / n_tag0 * 100, verluste = n_tag0 - n_tiere)

p6h <- ggplot(retention, aes(x = days_since_rel, y = n_tiere)) +
  annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = Inf,
           fill = "#F1A340", alpha = 0.08) +
  annotate("text", x = 3.5, y = Inf, vjust = 1.8, size = 3.2,
           label = "Erste Woche", color = "#F1A340", fontface = "italic") +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.8) +
  geom_step(color = COL_IGEL, linewidth = 1.8, direction = "hv") +
  geom_point(color = COL_IGEL, size = 3, shape = 16) +
  geom_text(aes(label = n_tiere), vjust = -0.8, size = 3.2, color = "grey40") +
  scale_y_continuous(name     = "Anzahl noch aktiver Tiere",
                     limits   = c(0, n_tag0 * 1.12),
                     breaks   = seq(0, n_tag0, by = 5),
                     sec.axis = sec_axis(transform = ~ . / n_tag0 * 100,
                                         name   = "Anteil noch aktiver Tiere (%)",
                                         breaks = seq(0, 100, by = 25),
                                         labels = function(x) paste0(x, "%"))) +
  scale_x_continuous(breaks = seq(0, DAYS_MAX, by = 5),
                     expand = expansion(add = c(0.3, 0.5))) +
  labs(title    = "Retention-Kurve: Wann verschwinden Tiere aus dem Datensatz?",
       subtitle = sprintf("n = %d Tiere bei Auswilderung  |  Abnahme = Senderverlust ODER Abwanderung außer Reichweite",
                           n_tag0),
       x = "Tage nach Auswilderung") +
  theme_igel

ggsave(file.path(out_ordner, "06h_retention_kurve.png"), p6h,
       width = 12, height = 5.5, dpi = 200, bg = "white")
cat("  ✓ 06h_retention_kurve.png\n")

# ── Plot I — Enddistanz ───────────────────────────────────────────────────────
p6i <- ggplot(endpunkte,
              aes(x = letzter_tag, y = enddistanz, fill = abgang_typ, color = abgang_typ)) +
  geom_hline(yintercept = RANGE_SCHWELLE_M,
             linetype = "dashed", color = "firebrick", linewidth = 0.8, alpha = 0.7) +
  annotate("label", x = DAYS_MAX, y = RANGE_SCHWELLE_M,
           label = paste0("Detektionsgrenze (~", RANGE_SCHWELLE_M, " m)"),
           hjust = 1, size = 3.1, color = "firebrick",
           fill = "white", label.padding = unit(2, "pt")) +
  geom_hline(yintercept = ABGANG_SCHWELLE_M,
             linetype = "dotted", color = "grey50", linewidth = 0.7) +
  annotate("text", x = 0.5, y = ABGANG_SCHWELLE_M + 6,
           label = paste0("Schwelle Abwanderung (", ABGANG_SCHWELLE_M, " m)"),
           hjust = 0, size = 3.0, color = "grey50") +
  geom_point(size = 5, shape = 21, alpha = 0.85, color = "grey20") +
  { if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(aes(label = tier_label), size = 3.0, color = "grey20",
                                box.padding = 0.35, point.padding = 0.3,
                                max.overlaps = 15, segment.color = "grey70",
                                segment.size = 0.3, show.legend = FALSE)
    else
      geom_text(aes(label = tier_label), vjust = -1.0, size = 2.8, color = "grey20") } +
  scale_fill_manual(values  = abgang_farben, name = "Abgangstyp") +
  scale_color_manual(values = abgang_farben, guide = "none") +
  scale_x_continuous(breaks = seq(0, max(endpunkte$letzter_tag) + 2, by = 5),
                     expand = expansion(add = c(0.5, 1.5))) +
  scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                     expand = expansion(mult = c(0.03, 0.10))) +
  labs(title    = "Enddistanz vs. Trackingdauer: Abwanderung oder Senderverlust?",
       subtitle = "Jeder Punkt = 1 Tier  |  Rot = Abwanderung vermutet (>200 m)  |  Blau = Verbleib / Senderverlust",
       x = "Letzter Tracking-Tag (Tage nach Auswilderung)",
       y = "Distanz am letzten Tracking-Tag (m)") +
  theme_igel + theme(legend.position = "bottom",
                     legend.text     = element_text(size = 9, lineheight = 1.1),
                     legend.key.size = unit(0.8, "lines"))

ggsave(file.path(out_ordner, "06i_enddistanz_analyse.png"), p6i,
       width = 11, height = 7, dpi = 200, bg = "white")
cat("  ✓ 06i_enddistanz_analyse.png\n")

# ── Plot J — Trajektorien nach Abgangstyp ─────────────────────────────────────
traj_data <- daily_dist %>%
  filter(days_since_rel >= 0) %>%
  left_join(endpunkte %>% select(tier_label, abgang_typ, letzter_tag), by = "tier_label")

endpunkt_seg <- traj_data %>%
  group_by(tier_label, abgang_typ) %>%
  filter(days_since_rel >= letzter_tag - 2) %>%
  ungroup()

p6j <- ggplot() +
  annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = Inf,
           fill = "#F1A340", alpha = 0.06) +
  geom_hline(yintercept = RANGE_SCHWELLE_M,
             linetype = "dashed", color = "firebrick", linewidth = 0.7, alpha = 0.5) +
  annotate("text", x = DAYS_MAX, y = RANGE_SCHWELLE_M + 5,
           label = paste0("Detektionsgrenze ~", RANGE_SCHWELLE_M, " m"),
           hjust = 1, size = 2.9, color = "firebrick") +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.8) +
  geom_line(data = traj_data,
            aes(x = days_since_rel, y = dist_median, group = tier_label, color = abgang_typ),
            linewidth = 0.65, alpha = 0.50) +
  geom_line(data = endpunkt_seg,
            aes(x = days_since_rel, y = dist_median, group = tier_label, color = abgang_typ),
            linewidth = 2.0, alpha = 0.90) +
  geom_point(data = endpunkte,
             aes(x = letzter_tag, y = enddistanz, fill = abgang_typ),
             size = 4, shape = 21, color = "grey20", alpha = 0.90) +
  { if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(
        data = endpunkte %>% filter(grepl("Abwanderung", as.character(abgang_typ))),
        aes(x = letzter_tag, y = enddistanz, label = tier_label),
        size = 3.0, color = "grey20", nudge_x = 0.5, nudge_y = 8,
        segment.color = "grey60", segment.size = 0.3,
        max.overlaps = 10, show.legend = FALSE) } +
  scale_color_manual(values  = abgang_farben, name = "Abgangstyp") +
  scale_fill_manual(values   = abgang_farben, guide = "none") +
  scale_y_continuous(labels = label_number(suffix = " m"), limits = c(0, NA),
                     expand = expansion(mult = c(0.03, 0.08))) +
  scale_x_continuous(breaks = seq(0, DAYS_MAX, by = 5),
                     expand = expansion(add = c(0.3, 1.5))) +
  labs(title    = "Individuelle Trajektorien — Abgangstyp hervorgehoben",
       subtitle = "Dicker Strich = letzte 3 Tracking-Tage  |  Punkt = Endpunkt  |  Rot = Abwanderung  |  Blau = Verbleib/Senderverlust",
       x = "Tage nach Auswilderung", y = "Distanz zur Voliere (m)") +
  theme_igel + theme(legend.position = "bottom",
                     legend.text     = element_text(size = 9, lineheight = 1.1),
                     legend.key.size = unit(0.8, "lines"))

ggsave(file.path(out_ordner, "06j_trajektorien_abgang.png"), p6j,
       width = 13, height = 7, dpi = 200, bg = "white")
cat("  ✓ 06j_trajektorien_abgang.png\n")

# ── Plot K — Slope letzte 5 Tage ──────────────────────────────────────────────
slope_data <- endpunkte %>%
  filter(!is.na(slope_letzte5), letzter_tag >= ABGANG_MIN_TAGE) %>%
  mutate(richtung = case_when(
    slope_letzte5 >  5  ~ "Zunehmend\n(Tier entfernt sich)",
    slope_letzte5 < -5  ~ "Abnehmend\n(Tier nähert sich)",
    TRUE                ~ "Stabil\n(±5 m/Tag)"
  ),
  richtung = factor(richtung, levels = c("Abnehmend\n(Tier nähert sich)",
                                          "Stabil\n(±5 m/Tag)",
                                          "Zunehmend\n(Tier entfernt sich)")))

richtung_farben <- c("Abnehmend\n(Tier nähert sich)"   = "#4393C3",
                      "Stabil\n(±5 m/Tag)"              = COL_IGEL,
                      "Zunehmend\n(Tier entfernt sich)" = "#D6604D")

p6k <- ggplot(slope_data,
              aes(x = reorder(tier_label, slope_letzte5), y = slope_letzte5, fill = richtung)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.8) +
  geom_hline(yintercept = c(-5, 5), color = "grey70", linewidth = 0.5, linetype = "dotted") +
  geom_col(alpha = 0.85, color = "grey30", linewidth = 0.4, width = 0.7) +
  { if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(
        aes(label = sprintf("%+.1f", slope_letzte5),
            vjust = ifelse(slope_letzte5 >= 0, -0.4, 1.4)),
        size = 3.0, color = "grey25", direction = "y",
        max.overlaps = 20, show.legend = FALSE)
    else
      geom_text(aes(label = sprintf("%+.1f", slope_letzte5),
                    vjust = ifelse(slope_letzte5 >= 0, -0.4, 1.4)),
                size = 2.8, color = "grey25") } +
  scale_fill_manual(values = richtung_farben, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(sprintf("%+.0f", x), " m/Tag")) +
  labs(title    = "Distanzveränderung in den letzten 5 Tracking-Tagen",
       subtitle = sprintf("Nur Tiere mit ≥ %d Tracking-Tagen (n = %d)  |  Positive Werte = Tier entfernte sich zuletzt",
                           ABGANG_MIN_TAGE, nrow(slope_data)),
       x = NULL, y = "Distanz-Slope letzte 5 Tage (m/Tag)") +
  coord_flip() +
  theme_igel + theme(legend.position = "bottom",
                     legend.text     = element_text(size = 9, lineheight = 1.1),
                     axis.text.y     = element_text(size = 10))

ggsave(file.path(out_ordner, "06k_slope_letzte_woche.png"), p6k,
       width = 10, height = max(5, nrow(slope_data) * 0.42 + 2),
       dpi = 200, bg = "white")
cat("  ✓ 06k_slope_letzte_woche.png\n")

# Zusammenfassung
cat("\n  ── Persistenz-Zusammenfassung ───────────────────────────\n")
tab <- endpunkte %>% count(abgang_typ)
for (i in seq_len(nrow(tab)))
  cat(sprintf("    %-45s %d Tiere\n",
              gsub("\n", " ", as.character(tab$abgang_typ[i])), tab$n[i]))
n_abwanderer <- sum(grepl("Abwanderung", as.character(endpunkte$abgang_typ)))
cat(sprintf("\n    ► %d von %d Tieren verlassen das Gebiet (Enddistanz > %d m)\n",
            n_abwanderer, nrow(endpunkte), ABGANG_SCHWELLE_M))
slope_pos <- slope_data %>% filter(slope_letzte5 > 5)
cat(sprintf("    ► %d Tiere zeigten zunehmende Distanz in der letzten Woche\n", nrow(slope_pos)))

# ── Excel-Export Persistenz ────────────────────────────────────────────────────
cat("\n── Excel-Export: Persistenz-Analyse ────────────────────────\n")

wb_p <- createWorkbook()

add_p_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeDataTable(wb, name, as.data.frame(data),
                 tableStyle = "TableStyleMedium9", withFilter = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}

abgang_summary <- endpunkte %>%
  mutate(abgang_typ_clean = gsub("\n", " ", as.character(abgang_typ))) %>%
  select(Tier = tier_label, Geschlecht = sex, Letzter_Tag = letzter_tag,
         Enddistanz_m = enddistanz, MaxDistanz_m = max_dist_gesamt,
         Slope_letzte5 = slope_letzte5, Abgangstyp = abgang_typ_clean) %>%
  mutate(across(where(is.numeric), ~ round(., 1))) %>%
  arrange(Abgangstyp, Letzter_Tag)

add_p_sheet(wb_p, "01_Abgangstypen", abgang_summary)
add_p_sheet(wb_p, "02_Retention_Kurve",
            retention %>% transmute(Tag_nach_Auswilderung = days_since_rel,
                                     N_aktiv = n_tiere,
                                     Anteil_aktiv_pct = round(pct_aktiv, 1),
                                     Kumulativer_Verlust = verluste))
if (nrow(slope_data) > 0)
  add_p_sheet(wb_p, "03_Slope_letzte5Tage",
              slope_data %>%
                mutate(richtung_clean = gsub("\n", " ", as.character(richtung))) %>%
                select(Tier = tier_label, Geschlecht = sex,
                       Letzter_Tag = letzter_tag, Enddistanz_m = enddistanz,
                       Slope_m_pro_Tag = slope_letzte5, Richtung = richtung_clean) %>%
                mutate(across(where(is.numeric), ~ round(., 2))) %>%
                arrange(Slope_m_pro_Tag))

add_p_sheet(wb_p, "04_Trajektorien_Rohdaten",
            traj_data %>%
              mutate(abgang_typ_clean = gsub("\n", " ", as.character(abgang_typ))) %>%
              select(Tier = tier_label, Geschlecht = sex, Tag = days_since_rel,
                     Distanz_Median_m = dist_median, Distanz_Q25_m = dist_q25,
                     Distanz_Q75_m = dist_q75, N_Fixes = n_fixes,
                     Letzter_Tag = letzter_tag, Abgangstyp = abgang_typ_clean) %>%
              mutate(across(where(is.numeric), ~ round(., 1))) %>%
              arrange(Tier, Tag))

zusammenfassung_df <- data.frame(
  Kennzahl = c("Tiere gesamt (Tag 0)", "Tiere noch aktiv Tag 5",
                "Tiere noch aktiv Tag 10", "Tiere noch aktiv Tag 20",
                paste0("Abwanderung vermutet (Enddistanz >", ABGANG_SCHWELLE_M, " m, ≥5 Tage)"),
                paste0("Verbleib/Senderverlust (Enddistanz ≤", ABGANG_SCHWELLE_M, " m, ≥5 Tage)"),
                paste0("Kurz-Tracking (<", ABGANG_MIN_TAGE, " Tage)"),
                "Tiere mit zunehm. Distanztrend letzte 5 Tage (>+5 m/Tag)",
                "Tiere mit stabiler Distanz letzte 5 Tage (±5 m/Tag)",
                "Tiere mit abnehm. Distanztrend letzte 5 Tage (<-5 m/Tag)",
                "Median Enddistanz alle Tiere [m]", "Median Letzter-Tracking-Tag [Tage]"),
  Wert = c(
    n_tag0,
    { r <- retention %>% filter(days_since_rel == 5) %>% pull(n_tiere); if(length(r)==0) NA else r },
    { r <- retention %>% filter(days_since_rel == 10) %>% pull(n_tiere); if(length(r)==0) NA else r },
    { r <- retention %>% filter(days_since_rel == 20) %>% pull(n_tiere); if(length(r)==0) NA else r },
    sum(grepl("Abwanderung",   as.character(endpunkte$abgang_typ))),
    sum(grepl("Verbleib",      as.character(endpunkte$abgang_typ))),
    sum(grepl("Kurz-Tracking", as.character(endpunkte$abgang_typ))),
    if (nrow(slope_data) > 0) sum(slope_data$slope_letzte5 >  5, na.rm = TRUE) else NA,
    if (nrow(slope_data) > 0) sum(abs(slope_data$slope_letzte5) <= 5, na.rm = TRUE) else NA,
    if (nrow(slope_data) > 0) sum(slope_data$slope_letzte5 < -5, na.rm = TRUE) else NA,
    round(median(endpunkte$enddistanz, na.rm = TRUE), 1),
    round(median(endpunkte$letzter_tag, na.rm = TRUE), 1)
  ), stringsAsFactors = FALSE)
add_p_sheet(wb_p, "00_Zusammenfassung", zusammenfassung_df)

legende_p <- data.frame(
  Begriff = c("Letzter_Tag","Enddistanz_m","MaxDistanz_m","Slope_letzte5",
               "Abgangstyp: Abwanderung","Abgangstyp: Verbleib / Senderverlust",
               "Abgangstyp: Kurz-Tracking","Detektionsgrenze","Schwelle Abwanderung"),
  Erklaerung = c(
    "Letzter Tag mit mindestens 1 Fix (Tage nach Auswilderung)",
    "Mediandistanz zur Voliere am letzten Tracking-Tag [m]",
    "Maximale Mediandistanz über die gesamte Trackingperiode [m]",
    "Steigung der linearen Regression Distanz ~ Tag über die letzten 5 Tracking-Tage [m/Tag]",
    paste0("Enddistanz > ", ABGANG_SCHWELLE_M, " m UND Tracking ≥ ", ABGANG_MIN_TAGE, " Tage"),
    paste0("Enddistanz ≤ ", ABGANG_SCHWELLE_M, " m UND Tracking ≥ ", ABGANG_MIN_TAGE, " Tage"),
    paste0("Tracking < ", ABGANG_MIN_TAGE, " Tage — zu kurz für Klassifikation"),
    paste0(RANGE_SCHWELLE_M, " m — ab diesem Abstand werden Tiere kaum noch detektiert"),
    paste0(ABGANG_SCHWELLE_M, " m — heuristische Grenze Abwanderung vs. Verbleib")),
  stringsAsFactors = FALSE)
addWorksheet(wb_p, "Legende")
writeData(wb_p, "Legende", legende_p)
s_hdr_p <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white",
                        fgFill = "#7B4F2E", textDecoration = "bold",
                        halign = "center", border = "Bottom")
addStyle(wb_p, "Legende", s_hdr_p, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb_p, "Legende", cols = 1, widths = 30)
setColWidths(wb_p, "Legende", cols = 2, widths = 75)

saveWorkbook(wb_p, file.path(out_ordner, "06_persistenz_ergebnisse.xlsx"), overwrite = TRUE)
cat("  ✓ 06_persistenz_ergebnisse.xlsx\n")

# ── 13. LME-Modell ────────────────────────────────────────────────────────────
cat("\n── LME-Modell: Zeittrend ────────────────────────────────────\n")

lme_data <- daily_dist %>% filter(days_since_rel >= 0, !is.na(dist_median))

mod_zeit <- tryCatch(
  lmer(dist_median ~ days_since_rel + (1 + days_since_rel | tier_label),
       data = lme_data, REML = TRUE, control = lmerControl(optimizer = "bobyqa")),
  error = function(e) {
    cat("  [WARN] LME mit random slope fehlgeschlagen:", conditionMessage(e), "\n")
    tryCatch(
      lmer(dist_median ~ days_since_rel + (1 | tier_label), data = lme_data, REML = TRUE),
      error = function(e2) { cat("  [WARN] LME fehlgeschlagen:", conditionMessage(e2), "\n"); NULL })
  })

if (!is.null(mod_zeit)) {
  coef_tab <- as.data.frame(summary(mod_zeit)$coefficients)
  print(round(coef_tab, 4))
  slope <- coef_tab["days_since_rel", "Estimate"]
  se    <- coef_tab["days_since_rel", "Std. Error"]
  pval  <- coef_tab["days_since_rel", "Pr(>|t|)"]
  cat(sprintf("\n  ► Slope: %.2f m/Tag (SE = %.2f, p = %.4f)\n", slope, se, pval))
  if (pval < 0.05) {
    if (slope > 0) cat("  → Signifikante Zunahme der Distanz (Dispersion).\n")
    else           cat("  → Signifikante Abnahme der Distanz (Kontraktion).\n")
  } else {
    cat("  → Kein signifikanter Zeittrend (p ≥ 0.05).\n")
  }
}

mod_sex <- NULL
if (n_distinct(na.omit(lme_data$sex)) >= 2) {
  lme_sex_data <- lme_data %>% filter(!is.na(sex)) %>% mutate(sex = factor(sex))
  mod_sex <- tryCatch(
    lmer(dist_median ~ days_since_rel * sex + (1 | tier_label),
         data = lme_sex_data, REML = TRUE, control = lmerControl(optimizer = "bobyqa")),
    error = function(e) { cat("  [WARN] Sex-Modell:", conditionMessage(e), "\n"); NULL })
  if (!is.null(mod_sex)) {
    cat("\n  Geschlechtseffekt-Modell:\n")
    print(round(as.data.frame(summary(mod_sex)$coefficients), 4))
  }
}

# ── 14. Excel-Export LME ──────────────────────────────────────────────────────
cat("\n── Excel-Export LME ─────────────────────────────────────────\n")

wb <- createWorkbook()
s_hdr  <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white",
                       fgFill = "#2E6B34", textDecoration = "bold",
                       halign = "center", border = "Bottom")
s_bold <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold")

add_sheet_styled <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  writeDataTable(wb, sheet_name, as.data.frame(data),
                 tableStyle = "TableStyleMedium7", withFilter = TRUE)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(data)), widths = "auto")
}

write_lme_sheet <- function(wb, sheet_name, mod, titel) {
  addWorksheet(wb, sheet_name)
  if (is.null(mod)) {
    writeData(wb, sheet_name, data.frame(Hinweis = paste(titel, "— Modell nicht berechnet.")))
    return(invisible())
  }
  coefs <- as.data.frame(summary(mod)$coefficients) %>%
    tibble::rownames_to_column("Term") %>%
    mutate(across(where(is.numeric), ~ round(., 4)))
  writeData(wb, sheet_name, data.frame(Abschnitt = paste0("Modell: ", titel)), startRow = 1)
  writeData(wb, sheet_name, coefs, startRow = 2)
  addStyle(wb, sheet_name, s_hdr,  rows = 2, cols = seq_len(ncol(coefs)), gridExpand = TRUE)
  addStyle(wb, sheet_name, s_bold, rows = 1, cols = 1)
  re_row <- nrow(coefs) + 5
  writeData(wb, sheet_name, data.frame(Abschnitt = "Zufallseffekte (Varianz)"), startRow = re_row)
  re_dt <- as.data.frame(VarCorr(mod)) %>% mutate(across(where(is.numeric), ~ round(., 4)))
  writeData(wb, sheet_name, re_dt, startRow = re_row + 1)
  setColWidths(wb, sheet_name, cols = 1:max(ncol(coefs), ncol(re_dt)), widths = "auto")
}

add_sheet_styled(wb, "Tagesdistanzen_Tier",
                 daily_dist %>%
                   select(tier_label, sex, days_since_rel, dist_median, dist_q25, dist_q75, n_fixes) %>%
                   mutate(across(where(is.numeric), ~ round(., 1))))
add_sheet_styled(wb, "Gruppenstatistik",
                 mean_dist %>% mutate(across(where(is.numeric), ~ round(., 1))))
add_sheet_styled(wb, "ProTier_Zusammenfassung",
                 animal_summary %>%
                   left_join(meta %>% select(tier_label, date_release, sr_dauer_tage,
                                              diagnosis, tagging_weight), by = "tier_label") %>%
                   mutate(across(where(is.numeric), ~ round(., 1))))
write_lme_sheet(wb, "LME_Zeittrend",  mod_zeit, "Distanz ~ Tage seit Auswilderung")
write_lme_sheet(wb, "LME_Geschlecht", mod_sex,  "Distanz ~ Tage × Geschlecht")

legend_df <- data.frame(
  Begriff = c("days_since_rel","dist_median","dist_q25 / dist_q75",
               "n_fixes","sr_dauer_tage","Intercept (LME)","days_since_rel (LME)","Volieren-Radius"),
  Erklaerung = c(
    "Tage seit Auswilderung (0 = Auswilderungstag)",
    "Mediandistanz aller Multilateralisierungen eines Tieres an einem Tag zur Voliere [m]",
    "25. / 75. Perzentile der Tages-Fixes [m]",
    "Anzahl VHF-Multilateralisierungen an diesem Tag",
    "Dauer der Soft-Release-Phase in der Voliere (date_release - soft_release_start)",
    "Geschätzte mittlere Distanz am Auswilderungstag (Tag 0)",
    "Veränderung der Distanz pro Tag (positiv = Zunahme = Dispersion)",
    paste0(VOLIERE_RADIUS_M, " m — ungefährer Volieren-Durchmesser")),
  stringsAsFactors = FALSE)
addWorksheet(wb, "Legende")
writeData(wb, "Legende", legend_df)
addStyle(wb, "Legende", s_hdr, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "Legende", cols = 1, widths = 25)
setColWidths(wb, "Legende", cols = 2, widths = 70)

saveWorkbook(wb, file.path(out_ordner, "06_lme_ergebnisse.xlsx"), overwrite = TRUE)
cat("  ✓ 06_lme_ergebnisse.xlsx\n")

# ── 15. RDS-Export ────────────────────────────────────────────────────────────
saveRDS(
  list(alle_fixes    = alle_fixes,
       daily_dist    = daily_dist,
       mean_dist     = mean_dist,
       animal_summary= animal_summary,
       endpunkte     = endpunkte,
       meta          = meta,
       mod_zeit      = mod_zeit,
       mod_sex       = mod_sex,
       voliere_xy    = vol_xy),
  file = file.path(out_ordner, "06_dispersal_ergebnisse.rds"))
cat("  ✓ 06_dispersal_ergebnisse.rds\n")

# ══════════════════════════════════════════════════════════════════════════════
# 16. Word-Methodenbericht
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Word-Methodenbericht ─────────────────────────────────────\n")

if (!requireNamespace("officer",   quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")

if (requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)) {
  library(officer)
  library(flextable)

  make_ft <- function(df, hl_rows = NULL) {
    brd <- fp_border(color = "#BFBFBF", width = 0.5)
    ft  <- flextable(df) |>
      font(fontname = "Arial", part = "all") |>
      fontsize(size = 10, part = "all") |>
      bold(part = "header") |>
      color(part = "header", color = "white") |>
      bg(part = "header", bg = "#2E6B34") |>
      hline(border = brd, part = "all") |>
      vline(border = brd, part = "all") |>
      hline_top(border = brd, part = "header") |>
      bg(part = "body", bg = "white") |>
      set_table_properties(layout = "autofit")
    if (!is.null(hl_rows) && length(hl_rows) > 0)
      ft <- bg(ft, i = hl_rows, bg = "#E2EFDA", part = "body")
    ft
  }

  doc <- read_docx()

  # Titel
  doc <- doc |>
    body_add_par("Block 6: Dispersal-Analyse — Raumnutzung nach Auswilderung",
                 style = "heading 1") |>
    body_add_par("Igelbesenderung Niedersachsen — TiHo Hannover | Natalie Steiner",
                 style = "Normal") |>
    body_add_par(paste0("Erstellt: ", format(Sys.Date(), "%d. %B %Y")), style = "Normal") |>
    body_add_par("", style = "Normal")

  # 1. Fragestellung
  doc <- doc |>
    body_add_par("1. Fragestellung", style = "heading 2") |>
    body_add_par(paste0(
      "Wie weit und wie schnell entfernen sich ausgewilderte Igel von der Auswilderungsvoliere? ",
      "Gibt es eine initiale Soft-Release-Phase mit reduzierter Mobilität? ",
      "Wann verlassen Tiere das Detektionsgebiet — durch Abwanderung oder Senderverlust? ",
      "Unterscheidet sich die Dispersal-Dynamik zwischen Männchen und Weibchen?"
    ), style = "Normal") |>
    body_add_par("", style = "Normal")

  # 2. Abgrenzung zu Block 4
  doc <- doc |>
    body_add_par("2. Abgrenzung zu Block 4 (Kernel Home Range)", style = "heading 2") |>
    body_add_par(paste0(
      "Block 6 ist komplementär zu Block 4 und nicht redundant: ",
      "Block 4 berechnet statische Raumnutzungsmetriken über den gesamten Beobachtungszeitraum ",
      "(KDE-Aktionsareal, max. Nachtdistanz, Rückkehrrate). ",
      "Block 6 analysiert die zeitliche Dynamik der Distanz zur Voliere in den ersten 30 Tagen — ",
      "speziell die Soft-Release-Hypothese, die Retention-Kurve und die Endpunkt-Klassifikation."
    ), style = "Normal") |>
    body_add_par("", style = "Normal")

  # 3. Methodik
  doc <- doc |>
    body_add_par("3. Methodik", style = "heading 2") |>
    body_add_par("3.1 Datenbasis", style = "heading 3") |>
    body_add_par(paste0(
      "Datenquelle: VHF-Multilateralisierungen (GeoPackages) aus dem Ordner data/kernel_files/. ",
      "Filter: Fixes mit Station Count ≥ 2 (mindestens 2 Antennen für Multilateralisation). ",
      "Koordinatensystem: UTM Zone 32N (EPSG 25832). ",
      "Analysefenster: Tag 0 bis Tag ", DAYS_MAX, " nach Auswilderung. ",
      "Tägliche Aggregation: Mediandistanz aller Fixes eines Tieres pro Tag zur Voliere."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("3.2 Volieren-Koordinaten", style = "heading 3") |>
    body_add_par(paste0(
      "Auswilderungsvoliere Wildtierstation Sachsenhagen: ",
      sprintf("Lat = %.6f, Lon = %.6f (WGS84). ", VOLIERE_LAT, VOLIERE_LON),
      sprintf("UTM: X = %.1f, Y = %.1f (EPSG 25832). ", vol_xy[1,1], vol_xy[1,2]),
      "Volieren-Radius: ~", VOLIERE_RADIUS_M, " m (visueller Richtwert, kein Ausschluss-Kriterium)."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("3.3 Abgangstyp-Klassifikation", style = "heading 3") |>
    body_add_par(paste0(
      "Heuristische Klassifikation in drei Kategorien anhand Enddistanz und Trackingdauer: ",
      "(1) Abwanderung: Enddistanz > ", ABGANG_SCHWELLE_M, " m UND Tracking ≥ ", ABGANG_MIN_TAGE, " Tage. ",
      "(2) Verbleib/Senderverlust: Enddistanz ≤ ", ABGANG_SCHWELLE_M, " m UND Tracking ≥ ", ABGANG_MIN_TAGE, " Tage. ",
      "(3) Kurz-Tracking: < ", ABGANG_MIN_TAGE, " Tage — zu kurz für Klassifikation. ",
      "Detektionsgrenze: ~", RANGE_SCHWELLE_M, " m (ab diesem Abstand zur nächsten Antenne kaum noch Signal)."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("3.4 Statistisches Modell", style = "heading 3") |>
    body_add_par(paste0(
      "Lineares gemischtes Modell (LME, lme4): Distanz ~ Tage seit Auswilderung + (1 + Tage | Tier). ",
      "Random slopes und random intercepts pro Tier. ",
      "Fallback auf random-intercept-only Modell bei Konvergenzproblemen. ",
      "Optimizer: bobyqa (robuster als Standard-REML für unbalancierte Daten). ",
      "GAM-Smoother (mgcv, k=6): nicht-parametrische Alternative für Visualisierung ohne Modellannahmen."
    ), style = "Normal") |>
    body_add_par("", style = "Normal")

  # 4. Ergebnisse - Übersichtstabelle
  doc <- doc |>
    body_add_par("4. Ergebnisse — Übersicht", style = "heading 2")

  erg_df <- data.frame(
    Kennzahl = c("Analysierte Tiere", "Analysefenster", "Gesamtdistanz Median",
                 "Pro-Tier-Median (Median)", "Tiere Abwanderung vermutet",
                 "Tiere Verbleib/Senderverlust", "Tiere Kurz-Tracking",
                 "LME Slope (m/Tag)", "LME p-Wert"),
    Wert = c(
      n_distinct(daily_dist$tier_label),
      paste0("Tag 0 – ", DAYS_MAX),
      paste0(round(median(alle_fixes$dist_voliere_m, na.rm=TRUE)), " m"),
      paste0(round(median(animal_summary$median_dist, na.rm=TRUE)), " m"),
      n_abwanderer,
      sum(grepl("Verbleib", as.character(endpunkte$abgang_typ))),
      sum(grepl("Kurz-Tracking", as.character(endpunkte$abgang_typ))),
      if (!is.null(mod_zeit)) {
        ct <- as.data.frame(summary(mod_zeit)$coefficients)
        round(ct["days_since_rel","Estimate"], 2)
      } else "n.v.",
      if (!is.null(mod_zeit)) {
        ct <- as.data.frame(summary(mod_zeit)$coefficients)
        round(ct["days_since_rel","Pr(>|t|)"], 4)
      } else "n.v."
    ), stringsAsFactors = FALSE)

  doc <- doc |>
    body_add_flextable(make_ft(erg_df)) |>
    body_add_par("", style = "Normal")

  # Abgangstyp-Tabelle
  doc <- doc |>
    body_add_par("4.1 Abgangstyp-Klassifikation pro Tier", style = "heading 3")

  abgang_bericht <- endpunkte %>%
    mutate(Abgangstyp = gsub("\n", " ", as.character(abgang_typ))) %>%
    select(Tier = tier_label, Geschlecht = sex, `Letzter Tag` = letzter_tag,
           `Enddistanz (m)` = enddistanz, `Max. Distanz (m)` = max_dist_gesamt,
           `Slope (m/Tag)` = slope_letzte5, Abgangstyp) %>%
    mutate(across(where(is.numeric), ~ round(., 1))) %>%
    arrange(Abgangstyp, `Letzter Tag`) %>%
    as.data.frame()

  abgang_hl <- which(grepl("Abwanderung", abgang_bericht$Abgangstyp))

  doc <- doc |>
    body_add_flextable(make_ft(abgang_bericht, hl_rows = abgang_hl)) |>
    body_add_par("Grün = Abwanderung vermutet (Enddistanz > 200 m).",
                 style = "Normal") |>
    body_add_par("", style = "Normal")

  # 5. Interpretation
  doc <- doc |>
    body_add_par("5. Interpretation", style = "heading 2") |>
    body_add_par(paste0(
      "Der LME-Zeittrend zeigt, ob es eine systematische Veränderung der Distanz zur Voliere ",
      "über die 30 Tage gibt. Ein nicht-signifikanter Slope (p ≥ 0.05) bedeutet, dass die ",
      "Tiere im Mittel keine gerichtete Dispersion zeigen — was mit der Central-Place-Foraging-",
      "Hypothese aus Block 4 übereinstimmt (Igel kehren zur Voliere zurück)."
    ), style = "Normal") |>
    body_add_par(paste0(
      "Die Retention-Kurve zeigt wann Tiere aus dem Datensatz verschwinden. ",
      "Da Senderverlust und echte Abwanderung nicht immer unterscheidbar sind, ",
      "wird die Enddistanz als Proxy verwendet: Tiere mit Enddistanz > 200 m ",
      "haben das Detektionsgebiet wahrscheinlich verlassen, während Tiere mit ",
      "geringer Enddistanz eher einen Senderverlust in der Nähe der Voliere erlitten."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("6. Referenzen", style = "heading 2") |>
    body_add_par(paste0(
      "Bates, D., Mächler, M., Bolker, B., Walker, S. (2015). Fitting Linear Mixed-Effects ",
      "Models Using lme4. Journal of Statistical Software, 67(1), 1–48."
    ), style = "Normal") |>
    body_add_par(paste0(
      "Wood, S.N. (2017). Generalized Additive Models: An Introduction with R. ",
      "2nd edition. Chapman and Hall/CRC."
    ), style = "Normal")

  bericht_pfad <- file.path(out_ordner, "Block6_Dispersal_Methodenbericht.docx")
  print(doc, target = bericht_pfad)
  cat("✓ Block6_Dispersal_Methodenbericht.docx gespeichert\n")
} else {
  cat("  officer/flextable nicht verfügbar — Word-Bericht übersprungen\n")
}

# ── 17. Abschluss-Zusammenfassung ─────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════\n")
cat("  ERGEBNISSE — Dispersal-Analyse\n")
cat("══════════════════════════════════════════════════════════════\n\n")
cat(sprintf("  Tiere analysiert:          %d\n", n_distinct(daily_dist$tier_label)))
cat(sprintf("  Multilateralisierungen:    %d\n", nrow(alle_fixes)))
cat(sprintf("  Analysefenster:            Tag 0 bis Tag %d\n", DAYS_MAX))
cat(sprintf("\n  Distanz zur Voliere (alle Tiere, alle Tage):\n"))
cat(sprintf("    Median:  %.0f m\n", median(alle_fixes$dist_voliere_m, na.rm = TRUE)))
cat(sprintf("    Mittel:  %.0f m\n", mean(alle_fixes$dist_voliere_m, na.rm = TRUE)))
cat(sprintf("    Max:     %.0f m\n", max(alle_fixes$dist_voliere_m, na.rm = TRUE)))
cat(sprintf("\n  Pro-Tier-Median:\n"))
cat(sprintf("    Median:  %.0f m\n", median(animal_summary$median_dist, na.rm = TRUE)))
cat(sprintf("    Min:     %.0f m\n", min(animal_summary$median_dist, na.rm = TRUE)))
cat(sprintf("    Max:     %.0f m\n", max(animal_summary$median_dist, na.rm = TRUE)))
if (!is.null(mod_zeit)) {
  coef_tab <- as.data.frame(summary(mod_zeit)$coefficients)
  cat(sprintf("\n  LME Zeittrend:\n"))
  cat(sprintf("    Intercept: %.1f m  (Tag 0)\n",     coef_tab["(Intercept)",    "Estimate"]))
  cat(sprintf("    Slope:     %.2f m/Tag  (p = %.4f)\n",
              coef_tab["days_since_rel", "Estimate"],
              coef_tab["days_since_rel", "Pr(>|t|)"]))
}
cat(sprintf("\n  Output gespeichert in:\n  %s\n", out_ordner))
cat("══════════════════════════════════════════════════════════════\n")
cat("\n✓ Block6_Dispersal.R fertig.\n")
