# ==============================================================
# Chronobiological analysis — hedgehog activity after release
# ==============================================================
# Project:  Hedgehog VHF telemetry, Lower Saxony
# Author:   Natalie Steiner
#
# Methods (classical chronobiology, per individual):
#   1. Actograms  — double-plotted, sunrise/sunset
#   2. Cosinor    — mesor, amplitude, acrophase per individual
#                   and per time window (rolling weekly window)
#   3. Rayleigh   — circular test for rhythmicity per individual
#
# Rationale: GAMMs pooling all hedgehogs assume a shared rhythm
# structure. It is chronobiologically more correct to first estimate
# rhythm parameters per individual (cosinor, Rayleigh) and only then
# summarise across individuals.
#
# Requirement: Block0_Datenpipeline.R has been run and has saved
#   output/Block0_Pipeline/chrono_minuten24h.rds.
# ==============================================================
# IMPORTANT: always run the script from line 1 (Ctrl+Alt+R).
# ==============================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(circular)

cat("✓ Alle Pakete geladen\n\n")

# Hilfsfunktion für Word-Bericht (muss früh definiert sein)
if (!requireNamespace("officer",   quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
library(officer)
library(flextable)

add_img_safe <- function(doc, pfad, b = 15, h = 10) {
  if (file.exists(pfad))
    body_add_img(doc, pfad, width = b/2.54, height = h/2.54)
  else
    body_add_par(doc, paste0("[Bild fehlt: ", basename(pfad), "]"), style = "Normal")
}

# ──────────────────────────────────────────────────────────────
# EINSTELLUNGEN
# ──────────────────────────────────────────────────────────────

# Projektwurzel — einzige Zeile die du ggf. anpassen musst
projekt_root   <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"

output_ordner  <- file.path(projekt_root, "output", "Block3_Chronobiologie")
chrono_ordner  <- output_ordner   # Plots direkt in Block3_Chronobiologie/
dir.create(output_ordner, showWarnings = FALSE, recursive = TRUE)

# Sonnenaufgang/untergang für Hannover (52.40°N, 9.72°E)
# Wird für Actogramm-Überlagerung gebraucht
lat_h <- 52.397
lon_h <-  9.217

# Rollierendes Zeitfenster für Cosinor (in Tagen)
fenster_tage <- 7    # Fenstergröße
schritt_tage <- 3    # Schrittweite (Überlappung erlaubt)
min_stunden  <- 12   # Mindestanzahl Stunden mit Daten pro Fenster für Cosinor-Fit

# ── Einschluss-Schwellenwert: Mindestanzahl aktiver Nachtminuten ──────────────
# Begruendung: Konsistent mit Block 2 (>= 60 Nachtminuten als Mindestdatenbasis).
# Rayleigh-Test mit < 60 Datenpunkten produziert trivial signifikante Ergebnisse
# wenn alle Minuten zufaellig auf denselben Kreissektor fallen (Small-Sample-Artefakt).
# Fisher (1993), Zar (2010): mind. 1 Datenpunkt/h empfohlen fuer zirkulaere Statistik.
MIN_AKTIV_MIN <- 60L   # Mindestanzahl aktiver Nachtminuten fuer Einschluss in Analyse

# ──────────────────────────────────────────────────────────────
# DATEN LADEN
# ──────────────────────────────────────────────────────────────

cat("Lade 24h-Minutendaten...\n")

chrono_rds <- file.path(projekt_root, "output", "Block0_Pipeline", "chrono_minuten24h.rds")
if (!file.exists(chrono_rds)) {
  stop(paste0(
    "Datei '", chrono_rds, "' nicht gefunden.\n",
    "Bitte zuerst Block0_Datenpipeline.R ausführen."))
}

dt <- readRDS(chrono_rds)
setDT(dt)

cat("Geladene Minuten:", format(nrow(dt), big.mark = "'"), "\n")
cat("Igel:            ", nlevels(dt$igel), "\n")
cat("Datumsbereich:   ", format(min(dt$datum)), "–", format(max(dt$datum)), "\n\n")

# ── Einschluss-Filter: nur Tiere mit >= MIN_AKTIV_MIN aktiven Nachtminuten ───
# Gilt fuer Rayleigh-Test, Cosinor und alle weiteren Analysen in Block 3.
aktiv_pro_igel <- dt[aktiv == 1, .(n_aktiv_min = .N), by = igel]
igel_ok        <- aktiv_pro_igel[n_aktiv_min >= MIN_AKTIV_MIN, igel]
igel_liste_alle <- levels(dt$igel)   # alle Tiere im Datensatz
igel_liste      <- as.character(igel_ok)  # nur Tiere mit genuegend Daten

# ── Zusaetzlicher Ausschluss: unzureichende Tracking-Dauer ───────────────────
# H26 (Igel26) wurde nur 1 Tag getrackt. Ein Rayleigh-/Cosinor-Ergebnis aus einem
# einzigen Tag ist nicht belastbar. Ausschluss haelt Block 3 konsistent mit allen
# anderen Analysen (NAF, IS), aus denen H26 bereits wegen 1-Tages-Tracking
# ausgeschlossen ist (siehe Manuskript 3.5.1).
IGEL_EXCL_KURZ <- c("Igel26")
igel_liste     <- setdiff(igel_liste, IGEL_EXCL_KURZ)

# Ausgeschlossene Tiere dokumentieren
igel_excl <- setdiff(igel_liste_alle, igel_liste)
cat(sprintf("\n  Einschluss-Filter (>= %d aktive Nachtminuten):\n", MIN_AKTIV_MIN))
cat(sprintf("  Eingeschlossen: %d Tiere\n", length(igel_liste)))
if (length(igel_excl) > 0) {
  for (ig in igel_excl) {
    n_ig <- aktiv_pro_igel[igel == ig, n_aktiv_min]
    if (length(n_ig) == 0) n_ig <- 0L
    cat(sprintf("  [AUSGESCHLOSSEN] %s: %d aktive Min. (< %d)\n",
                ig, n_ig, MIN_AKTIV_MIN))
  }
}
cat("\n")

# Sonnenauf/-untergang laden (suncalc)
if (!requireNamespace("suncalc", quietly = TRUE)) install.packages("suncalc")
library(suncalc)

alle_daten   <- unique(dt$datum)
sonnen_dt    <- getSunlightTimes(
  date = alle_daten,
  lat  = lat_h, lon = lon_h,
  tz   = "Europe/Berlin",
  keep = c("sunrise", "sunset")
)
setDT(sonnen_dt)
sonnen_dt[, datum     := as.Date(date)]
sonnen_dt[, aufgang_h := as.numeric(format(sunrise, "%H")) +
                         as.numeric(format(sunrise, "%M")) / 60]
sonnen_dt[, untergang_h := as.numeric(format(sunset, "%H")) +
                           as.numeric(format(sunset, "%M")) / 60]

cat("✓ Sonnenauf/-untergang geladen\n\n")

# ══════════════════════════════════════════════════════════════
# 1. ACTOGRAMME
# ══════════════════════════════════════════════════════════════
# Standard in der Chronobiologie: doppelt aufgetragen (Double-Plot)
# Jede Zeile = 48h, verschoben um 24h → Phasendrift sofort sichtbar
# Dunkel = aktiv, Hell = passiv
# Gelbe Linie = Sonnenuntergang, Hellblaue Linie = Sonnenaufgang
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("1. Actogramme\n")
cat("════════════════════════════════════\n\n")

# Aggregation: 10-Minuten-Bins (für übersichtliche Actogramme)
dt[, bin10 := floor(stunde * 6) / 6]   # auf 10-Min-Raster runden

actogramm_data <- dt[, .(aktiv_anteil = mean(aktiv)),
                     by = .(igel, datum, tage_seit, bin10)]
actogramm_data <- merge(actogramm_data, sonnen_dt[, .(datum, aufgang_h, untergang_h)],
                        by = "datum", all.x = TRUE)

# Doppelt-Plot-Funktion
plot_actogramm <- function(igel_name) {

  ad <- actogramm_data[igel == igel_name]
  if (nrow(ad) == 0) return(NULL)

  t_range <- ad[, range(tage_seit)]
  t_seq   <- t_range[1]:t_range[2]

  # Double-Plot: Tag t erscheint in Zeile t (x: 0–24h) UND in Zeile t-1 (x: 24–48h)
  # → Zeile i zeigt Stunden 0–24 von Tag i UND 0–24 von Tag i+1 als 24–48
  dp_links  <- ad[, .(tage_seit, bin10, aktiv_anteil, aufgang_h, untergang_h)]
  dp_rechts <- copy(dp_links)
  dp_rechts[, bin10     := bin10 + 24]
  dp_rechts[, aufgang_h  := aufgang_h  + 24]
  dp_rechts[, untergang_h := untergang_h + 24]
  # Rechte Hälfte: Tag t erscheint eine Zeile HÖHER (tage_seit - 1)
  dp_rechts[, tage_seit := tage_seit - 1]

  dp_all <- rbind(dp_links, dp_rechts)
  dp_all <- dp_all[tage_seit %in% t_seq]

  # Sonnenlinien pro Tag (Mittelwert, da saisonal leicht variiert)
  sonnen_igel <- merge(
    data.table(tage_seit = t_seq),
    ad[, .(aufgang_h = mean(aufgang_h, na.rm = TRUE),
           untergang_h = mean(untergang_h, na.rm = TRUE)),
       by = tage_seit],
    by = "tage_seit", all.x = TRUE
  )
  # Rechte Hälfte der Sonnenlinien
  sonnen_rechts <- copy(sonnen_igel)
  sonnen_rechts[, aufgang_h   := aufgang_h   + 24]
  sonnen_rechts[, untergang_h := untergang_h + 24]
  sonnen_rechts[, tage_seit   := tage_seit - 1]
  sonnen_dp <- rbind(sonnen_igel, sonnen_rechts)
  sonnen_dp <- sonnen_dp[tage_seit %in% t_seq]

  ggplot(dp_all, aes(x = bin10, y = tage_seit, fill = aktiv_anteil)) +
    geom_tile(width = 1/6, height = 0.9) +
    # Sonnenuntergang (Nachtbeginn)
    geom_segment(data = sonnen_dp,
                 aes(x = untergang_h, xend = untergang_h,
                     y = tage_seit - 0.5, yend = tage_seit + 0.5),
                 color = "#F4A261", linewidth = 0.7, inherit.aes = FALSE) +
    # Sonnenaufgang (Nachtende)
    geom_segment(data = sonnen_dp,
                 aes(x = aufgang_h, xend = aufgang_h,
                     y = tage_seit - 0.5, yend = tage_seit + 0.5),
                 color = "#90E0EF", linewidth = 0.7, inherit.aes = FALSE) +
    scale_fill_gradient(low = "white", high = "#1a2744",
                        name = "Aktivität", limits = c(0, 1),
                        breaks = c(0, 0.5, 1),
                        labels = c("0%", "50%", "100%")) +
    scale_x_continuous(
      limits = c(0, 48),
      breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48),
      labels = c("00:00","06:00","12:00","18:00",
                 "00:00","06:00","12:00","18:00","00:00")
    ) +
    scale_y_reverse(breaks = pretty_breaks(n = 8)) +
    labs(
      title    = paste0("Actogramm — ", igel_name, " (Double-Plot)"),
      subtitle = "Dunkel = aktiv | Orange = Sonnenuntergang | Blau = Sonnenaufgang\nJede Zeile = 48h (überlappend) → Phasendrift sichtbar als Diagonalen",
      x = "Uhrzeit",
      y = "Tage nach Auswilderung"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold"),
      panel.grid    = element_blank(),
      axis.text.x   = element_text(size = 7),
      legend.position = "right"
    )
}

cat("Erstelle Actogramme für", length(igel_liste), "Igel...\n")

for (ig in igel_liste) {
  p <- plot_actogramm(ig)
  if (!is.null(p)) {
    pfad <- file.path(chrono_ordner, paste0("actogramm_", ig, ".png"))
    png(pfad, width = 1400, height = max(600, 100 * uniqueN(dt[igel == ig]$tage_seit)),
        res = 130)
    print(p)
    dev.off()
  }
}
cat("✓ Actogramme gespeichert in:", chrono_ordner, "\n\n")

# Übersichts-Actogramm: alle Igel nebeneinander (kleine Version)
cat("Erstelle Übersichts-Actogramm (alle Igel)...\n")

p_alle_acto <- ggplot(dt[tage_seit <= 21,
                          .(aktiv_anteil = mean(aktiv)),
                          by = .(igel, tage_seit, bin30 = floor(stunde * 2) / 2)],
                      aes(x = bin30, y = tage_seit, fill = aktiv_anteil)) +
  geom_tile(width = 0.5, height = 0.9) +
  facet_wrap(~ igel, ncol = 5) +
  scale_fill_gradient(low = "white", high = "#1a2744",
                      name = "Aktiv", limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                     labels = c("0h","6h","12h","18h","0h")) +
  scale_y_reverse(breaks = c(1, 7, 14, 21)) +
  labs(
    title    = "Actogramme — alle Igel, Tage 1–21 (Übersicht)",
    subtitle = "Dunkel = aktiv | Erste 21 Tage nach Auswilderung",
    x = "Uhrzeit", y = "Tag"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title    = element_text(face = "bold"),
    panel.grid    = element_blank(),
    strip.text    = element_text(face = "bold", size = 7),
    legend.position = "right"
  )

png(file.path(chrono_ordner, "actogramme_alle_uebersicht.png"),
    width = 2400, height = 1800, res = 130)
print(p_alle_acto)
dev.off()
cat("✓ Übersichts-Actogramm gespeichert\n\n")

# ══════════════════════════════════════════════════════════════
# 2. COSINOR-ANALYSE
# ══════════════════════════════════════════════════════════════
# Für jeden Igel × Zeitfenster:
#   Fit: y_h = Mesor + β·cos(2π·h/24) + γ·sin(2π·h/24) + ε
#   → Mesor    = mittleres Aktivitätsniveau
#   → Amplitude = sqrt(β² + γ²) = Stärke des 24h-Rhythmus
#   → Akrophase = atan2(-γ, β) × 24/(2π)  [in Stunden: peak-Zeit]
#
# Aggregation auf Stundenmittelwerte (robuster als Minuten)
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("2. Cosinor-Analyse\n")
cat("════════════════════════════════════\n\n")

# Stündliche Aggregation
dt_std <- dt[, .(
  aktiv_anteil = mean(aktiv),
  n_min        = .N
), by = .(igel, tage_seit, stunde_int = floor(stunde))]

# Cosinor für ein Zeitfenster fitten (Hilfsfunktion)
fit_cosinor <- function(dat) {
  # dat: data.table mit Spalten stunde_int, aktiv_anteil
  if (nrow(dat) < min_stunden || uniqueN(dat$stunde_int) < min_stunden) return(NULL)

  dat <- dat[!is.na(aktiv_anteil)]
  if (nrow(dat) < min_stunden) return(NULL)

  h <- dat$stunde_int
  y <- dat$aktiv_anteil

  # Lineares Cosinor-Modell
  cos_h <- cos(2 * pi * h / 24)
  sin_h <- sin(2 * pi * h / 24)

  fit <- tryCatch(lm(y ~ cos_h + sin_h), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  cf <- coef(fit)
  if (length(cf) < 3 || any(is.na(cf))) return(NULL)

  mesor     <- cf[1]                           # Intercept = Mesor
  beta      <- cf[2]                           # Koeffizient cos
  gamma     <- cf[3]                           # Koeffizient sin
  amplitude <- sqrt(beta^2 + gamma^2)          # Amplitude

  # Akrophase in Stunden (0–24): Peak der Aktivität
  acrophase_rad <- atan2(-gamma, beta)          # Winkel in Rad
  acrophase_h   <- acrophase_rad * 24 / (2*pi) # Umrechnung in Stunden
  if (acrophase_h < 0) acrophase_h <- acrophase_h + 24

  # R² des Fits
  r2 <- summary(fit)$r.squared

  # Signifikanz des Gesamtmodells (F-Test)
  fstat <- summary(fit)$fstatistic
  p_val <- if (!is.null(fstat)) {
    pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  } else NA_real_

  list(
    mesor      = as.numeric(mesor),
    amplitude  = as.numeric(amplitude),
    acrophase  = as.numeric(acrophase_h),
    r2         = as.numeric(r2),
    p_cosinor  = as.numeric(p_val),
    n_stunden  = nrow(dat)
  )
}

# Cosinor pro Igel × rollendes Zeitfenster
cat("Fitte Cosinor-Modelle (rollendes Fenster:", fenster_tage,
    "Tage, Schritt:", schritt_tage, "Tage)...\n\n")

cosinor_results <- list()

for (ig in igel_liste) {
  dt_ig   <- dt_std[igel == ig]
  t_range <- dt_ig[, range(tage_seit)]

  # Igel mit weniger Tagen als Fenstergröße überspringen
  if (t_range[2] - t_range[1] + 1 < fenster_tage) next

  starts  <- seq(t_range[1], t_range[2] - fenster_tage + 1, by = schritt_tage)

  if (length(starts) == 0) next

  for (t0 in starts) {
    t1  <- t0 + fenster_tage - 1
    dat <- dt_ig[tage_seit >= t0 & tage_seit <= t1]
    res <- fit_cosinor(dat)
    if (!is.null(res)) {
      cosinor_results[[paste0(ig, "_t", t0)]] <- c(
        list(igel = ig, t_start = t0, t_mitte = (t0 + t1) / 2),
        res
      )
    }
  }
}

# Cosinor auch für jeden Igel GESAMT (ein einzelner Wert pro Tier)
cosinor_gesamt <- list()
for (ig in igel_liste) {
  dat <- dt_std[igel == ig]
  res <- fit_cosinor(dat)
  if (!is.null(res)) {
    cosinor_gesamt[[ig]] <- c(list(igel = ig), res)
  }
}

cosinor_dt        <- rbindlist(cosinor_results, fill = TRUE)
cosinor_gesamt_dt <- rbindlist(cosinor_gesamt,  fill = TRUE)

# ── Dynamischer N-at-risk-Cutoff für rollende Zeitreihe ──────────
# Analog zur GAMM-Analyse: nur Zeitfenster mit >= n_min_cosinor
# Igeln mit gültigem Fit werden im Trend dargestellt.
n_min_cosinor <- 3   # Mindestanzahl Igel pro Zeitfenster

if (nrow(cosinor_dt) > 0) {
  # Anzahl Igel mit gültigem Fit pro Fenstermitte
  n_igel_pro_fenster <- cosinor_dt[, .(n_igel = uniqueN(igel)), by = t_mitte]
  cosinor_dt <- merge(cosinor_dt, n_igel_pro_fenster, by = "t_mitte", all.x = TRUE)

  # Cutoff = letzte Fenstermitte mit noch >= n_min_cosinor Igeln
  cutoff_cosinor <- n_igel_pro_fenster[n_igel >= n_min_cosinor, max(t_mitte)]
  cosinor_dt_plot <- cosinor_dt[t_mitte <= cutoff_cosinor]

  cat(sprintf(
    "N-at-risk-Cutoff rollender Cosinor: Tag %.1f (letzte Fenstermitte mit N ≥ %d Igeln)\n",
    cutoff_cosinor, n_min_cosinor))
  cat(sprintf("  Fenster im Plot: %d (von %d gesamt)\n\n",
              nrow(cosinor_dt_plot), nrow(cosinor_dt)))
} else {
  cosinor_dt_plot <- cosinor_dt
  cutoff_cosinor  <- NA_real_
}

cat("Cosinor-Fits berechnet:\n")
cat("  Rollend:    ", nrow(cosinor_dt), "Fenster über", length(igel_liste), "Igel\n")
cat("  Gesamt:     ", nrow(cosinor_gesamt_dt), "Igel\n\n")

# Übersicht Gesamt-Cosinor
cat("── Cosinor-Parameter (gesamt, pro Igel) ──\n")
print(cosinor_gesamt_dt[, .(
  igel,
  Mesor      = round(mesor      * 100, 1),
  Amplitude  = round(amplitude  * 100, 1),
  Akrophase  = round(acrophase, 1),
  R2         = round(r2, 3),
  p          = round(p_cosinor, 3)
)][order(Akrophase)])
cat("\nMesor und Amplitude in % | Akrophase = Stunde des Aktivitätsgipfels\n\n")

# ── Plots Cosinor ──────────────────────────────────────────────

# Plot C1: Akrophase aller Igel (Polarplot)
# Zeigt auf dem "Zifferblatt" wann jeder Igel seinen Aktivitätsgipfel hat
cosinor_gesamt_dt[, acro_rad := acrophase * 2 * pi / 24]

p_polar <- ggplot(cosinor_gesamt_dt,
                  aes(x = acro_rad, y = amplitude * 100,
                      color = igel, label = igel)) +
  geom_segment(aes(xend = acro_rad, yend = 0), linewidth = 1.2, alpha = 0.7) +
  geom_point(size = 3) +
  geom_text(hjust = -0.2, size = 2.8, alpha = 0.9) +
  # Markierungen für Nacht (21:00–05:00 = typisch Igel)
  annotate("rect",
           xmin = 21 * 2 * pi / 24, xmax = 24 * 2 * pi / 24,
           ymin = 0, ymax = Inf, alpha = 0.08, fill = "#1a2744") +
  annotate("rect",
           xmin = 0, xmax = 5 * 2 * pi / 24,
           ymin = 0, ymax = Inf, alpha = 0.08, fill = "#1a2744") +
  coord_polar(start = 0) +
  scale_x_continuous(
    limits  = c(0, 2 * pi),
    breaks  = (0:23) * 2 * pi / 24,
    labels  = sprintf("%02d:00", 0:23),
    expand  = c(0, 0)
  ) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  labs(
    title    = "Akrophase aller Igel — wann ist der Aktivitätsgipfel?",
    subtitle = "Jeder Strich = ein Igel | Länge = Amplitude (Rhythmusstärke)\nBlauer Sektor = typische Nachtphase",
    x = NULL, y = "Amplitude (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold"),
    axis.text.x   = element_text(size = 7),
    panel.grid.major = element_line(color = "grey85")
  )

# Plot C2: Akrophase über Zeit (rollend) pro Igel
if (nrow(cosinor_dt_plot) > 0) {

  # N-at-risk pro Fenstermitte für Beschriftung
  n_label_acro <- cosinor_dt_plot[, .(n_igel = uniqueN(igel)), by = t_mitte]

  p_acro_zeit <- ggplot(cosinor_dt_plot[p_cosinor < 0.05],  # nur signifikante Fits
                        aes(x = t_mitte, y = acrophase,
                            color = igel, group = igel)) +
    # Nachtbereiche
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = 21, ymax = 24, alpha = 0.07, fill = "#1a2744") +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin =  0, ymax =  5, alpha = 0.07, fill = "#1a2744") +
    geom_hline(yintercept = c(21, 5), linetype = "dashed",
               color = "#1a2744", alpha = 0.4) +
    # Cutoff-Linie
    {if (!is.na(cutoff_cosinor))
      geom_vline(xintercept = cutoff_cosinor, linetype = "dotted",
                 color = "grey40", linewidth = 0.8)} +
    geom_line(alpha = 0.5, linewidth = 0.8) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(aes(group = 1), method = "loess", span = 0.5,
                color = "black", linewidth = 1.3, se = TRUE,
                fill = "grey30", alpha = 0.15) +
    # N-at-risk am unteren Rand
    geom_text(data = n_label_acro,
              aes(x = t_mitte, y = 0.3, label = paste0("N=", n_igel)),
              inherit.aes = FALSE, size = 2.5, color = "grey50", vjust = 0) +
    scale_y_continuous(limits = c(0, 24),
                       breaks = seq(0, 24, 3),
                       labels = sprintf("%02d:00", seq(0, 24, 3))) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title    = "Akrophase über Zeit — verschiebt sich der Aktivitätsgipfel?",
      subtitle = paste0(
        "Nur signifikante Cosinor-Fits (p < 0.05) | Fenster: ", fenster_tage, " Tage\n",
        "Blauer Bereich = Nacht | Schwarz = Pop.-Trend | N = Anzahl Igel je Fenster",
        if (!is.na(cutoff_cosinor)) paste0(" | gestrichelt = N<", n_min_cosinor, "-Cutoff") else ""),
      x = "Tage nach Auswilderung",
      y = "Akrophase (Uhrzeit des Aktivitätsgipfels)",
      color = "Igel"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  # Plot C3: Amplitude über Zeit
  p_ampl_zeit <- ggplot(cosinor_dt_plot[p_cosinor < 0.05],
                        aes(x = t_mitte, y = amplitude * 100,
                            color = igel, group = igel)) +
    # Cutoff-Linie
    {if (!is.na(cutoff_cosinor))
      geom_vline(xintercept = cutoff_cosinor, linetype = "dotted",
                 color = "grey40", linewidth = 0.8)} +
    geom_line(alpha = 0.5, linewidth = 0.8) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(aes(group = 1), method = "loess", span = 0.5,
                color = "black", linewidth = 1.3, se = TRUE,
                fill = "grey30", alpha = 0.15) +
    # N-at-risk am unteren Rand
    geom_text(data = n_label_acro,
              aes(x = t_mitte, y = 0.2, label = paste0("N=", n_igel)),
              inherit.aes = FALSE, size = 2.5, color = "grey50", vjust = 0) +
    scale_color_viridis_d(option = "turbo") +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(
      title    = "Amplitude (Rhythmusstärke) über Zeit",
      subtitle = paste0(
        "Hohe Amplitude = stark ausgeprägter 24h-Rhythmus | N = Anzahl Igel je Fenster\n",
        "Schwarz = Pop.-Trend (LOESS)",
        if (!is.na(cutoff_cosinor)) paste0(" | gestrichelt = N<", n_min_cosinor, "-Cutoff") else ""),
      x = "Tage nach Auswilderung",
      y = "Amplitude (Prozentpunkte)",
      color = "Igel"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  png(file.path(chrono_ordner, "cosinor_akrophase_zeitverlauf.png"),
      width = 1600, height = 900, res = 130)
  print(p_acro_zeit)
  dev.off()

  png(file.path(chrono_ordner, "cosinor_amplitude_zeitverlauf.png"),
      width = 1600, height = 900, res = 130)
  print(p_ampl_zeit)
  dev.off()
}

# Plot C4: Cosinor-Kurven pro Igel (gefittete 24h-Profile)
cat("Erstelle gefittete Cosinor-Profile...\n")

cosinor_kurven <- lapply(igel_liste, function(ig) {
  res <- cosinor_gesamt[[ig]]
  if (is.null(res)) return(NULL)
  h_seq <- seq(0, 23.9, by = 0.1)
  y_fit <- res$mesor +
    res$amplitude * cos(2 * pi * h_seq / 24 - res$acrophase * 2 * pi / 24)
  data.table(igel = ig, stunde = h_seq, fit = y_fit * 100,
             p_sig = res$p_cosinor < 0.05)
})
cosinor_kurven_dt <- rbindlist(cosinor_kurven[!sapply(cosinor_kurven, is.null)])

# Beobachtete Stundenmittelwerte (Population)
dt_std_pop <- dt_std[, .(
  aktiv_pct = mean(aktiv_anteil, na.rm = TRUE) * 100
), by = .(igel, stunde_int)]

p_kurven <- ggplot() +
  geom_line(data = dt_std_pop,
            aes(x = stunde_int, y = aktiv_pct, group = igel),
            color = "grey70", alpha = 0.5, linewidth = 0.6) +
  geom_line(data = cosinor_kurven_dt[p_sig == TRUE],
            aes(x = stunde, y = fit, color = igel, group = igel),
            linewidth = 1.0, alpha = 0.8) +
  # Tagesdurchschnitt (thick black)
  geom_line(data = dt_std_pop[, .(aktiv_pct = mean(aktiv_pct)), by = stunde_int],
            aes(x = stunde_int, y = aktiv_pct),
            color = "black", linewidth = 1.4) +
  annotate("rect", xmin = -0.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf, fill = "#1a2744", alpha = 0.06) +
  annotate("rect", xmin = 20.5, xmax = 24,
           ymin = -Inf, ymax = Inf, fill = "#1a2744", alpha = 0.06) +
  scale_x_continuous(breaks = seq(0, 23, 3),
                     labels = sprintf("%02d:00", seq(0, 23, 3))) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  labs(
    title    = "Cosinor-Profile aller Igel — gefittete 24h-Kurven",
    subtitle = "Farbige Linien = Cosinor-Fit je Igel (nur p < 0.05)\nGraue Linien = beobachtete Stundenmittel | Schwarz = Populationsmittel",
    x = "Uhrzeit", y = "Aktivität (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

png(file.path(chrono_ordner, "cosinor_profile_alle.png"),
    width = 1400, height = 800, res = 130)
print(p_kurven)
dev.off()

png(file.path(chrono_ordner, "cosinor_akrophase_polar.png"),
    width = 1000, height = 1000, res = 130)
print(p_polar)
dev.off()

cat("✓ Cosinor-Plots gespeichert\n\n")

# ══════════════════════════════════════════════════════════════
# 3. RAYLEIGH-TEST (zirkuläre Statistik)
# ══════════════════════════════════════════════════════════════
# Fragestellung: Ist die Aktivität eines Igels signifikant auf
# eine Tageszeit konzentriert (= rhythmisch)?
#
# H0: gleichmäßige Verteilung über 24h (= arrhythmisch)
# H1: Aktivität konzentriert sich auf bestimmte Tageszeit
#
# Maße:
#   R̄  (mean resultant length, 0–1): Rhythmusstärke
#       0 = völlig gleichmäßig, 1 = alle Aktivität zu einem Zeitpunkt
#   μ   (mean direction): mittlere Aktivitätszeit
#   p   (Rayleigh-p): Signifikanz der Rhythmizität
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("3. Rayleigh-Test (Zirkuläre Statistik)\n")
cat("════════════════════════════════════\n\n")

# Alle aktiven Minuten → Zeitpunkt als Winkel
# θ = 2π × stunde / 24  (0 Uhr = 0, 24 Uhr = 2π)

rayleigh_res <- lapply(igel_liste, function(ig) {
  aktiv_min <- dt[igel == ig & aktiv == 1, stunde]

  # Einschluss-Filter bereits oben angewendet (igel_liste gefiltert).
  # Alle Tiere in igel_liste haben >= MIN_AKTIV_MIN aktive Minuten.

  theta  <- circular(aktiv_min * 2 * pi / 24,
                     type  = "angles",
                     units = "radians",
                     modulo = "2pi")

  # Rayleigh-Test
  rt     <- rayleigh.test(theta)

  # Mittlere Richtung (= mittlere Aktivitätszeit)
  mu_rad <- mean(theta)
  mu_h   <- as.numeric(mu_rad) * 24 / (2 * pi)
  if (mu_h < 0) mu_h <- mu_h + 24

  # Mean resultant length (Konzentration)
  rho    <- rho.circular(theta)

  data.table(
    igel          = ig,
    n_aktiv_min   = length(aktiv_min),
    rho           = round(as.numeric(rho), 3),
    mu_h          = round(mu_h, 2),
    rayleigh_p    = round(rt$p.value, 4),
    rhythmisch    = rt$p.value < 0.05
  )
})

rayleigh_dt <- rbindlist(rayleigh_res[!sapply(rayleigh_res, is.null)])
rayleigh_dt <- rayleigh_dt[order(rayleigh_p)]

cat("── Rayleigh-Test Ergebnisse ──\n")
print(rayleigh_dt[, .(
  igel,
  `N aktive Min.` = n_aktiv_min,
  `ρ̄ (Rhythmusstärke)` = rho,
  `μ (mittl. Aktivitätszeit)` = sprintf("%05.2f Uhr", mu_h),
  `p (Rayleigh)` = rayleigh_p,
  `rhythmisch?`  = ifelse(rhythmisch, "✓ ja (p<0.05)", "✗ nein")
)])

cat("\n")
n_rhythmisch <- sum(rayleigh_dt$rhythmisch, na.rm = TRUE)
cat("Rhythmische Igel (p < 0.05):", n_rhythmisch, "von", nrow(rayleigh_dt), "\n")
cat("Median ρ̄:", round(median(rayleigh_dt$rho, na.rm = TRUE), 3),
    "(0 = arrhythmisch, 1 = perfekter Rhythmus)\n\n")

# ══════════════════════════════════════════════════════════════
# 3b. SONNENUNTERGANGS-RELATIVE TIMING-METRIKEN
# ══════════════════════════════════════════════════════════════
# Begruendung: Die Acrophase/μ in Uhrzeit sind saisonabhaengig
# (Sonnenuntergang variiert zwischen den Auswilderungs-Monaten).
# Hier werden Cosinor-Acrophase und Rayleigh-μ RELATIV ZUM
# SONNENUNTERGANG berechnet ("Stunden nach Sonnenuntergang"),
# konsistent mit dem Aktivitaetsbeginn und mit Vazquez et al. (2019).
# Additiv — bestehende (Uhrzeit-)Auswertungen bleiben unveraendert.
# ──────────────────────────────────────────────────────────────
cat("════════════════════════════════════\n")
cat("3b. Timing relativ zum Sonnenuntergang\n")
cat("════════════════════════════════════\n\n")

# Sonnenuntergang pro Tag an die Minutendaten anfuegen
dt_rel <- merge(dt, sonnen_dt[, .(datum, untergang_h)], by = "datum", all.x = TRUE)
dt_rel[, stunde_rel := (stunde - untergang_h) %% 24]   # 0–24 h nach Sonnenuntergang

# Rayleigh-μ relativ zum Sonnenuntergang (pro Tier)
mu_rel_dt <- rbindlist(lapply(igel_liste, function(ig) {
  hr <- dt_rel[igel == ig & aktiv == 1 & !is.na(stunde_rel), stunde_rel]
  if (length(hr) < MIN_AKTIV_MIN) return(NULL)
  th <- circular(hr * 2 * pi / 24, type = "angles", units = "radians", modulo = "2pi")
  m  <- as.numeric(mean(th)) * 24 / (2 * pi); if (m < 0) m <- m + 24
  data.table(igel = ig,
             mu_rel_h  = round(m, 2),
             rho_rel   = round(as.numeric(rho.circular(th)), 3))
}))

# Cosinor-Acrophase relativ zum Sonnenuntergang (pro Tier, GESAMT)
# fit_cosinor erwartet Spalten 'stunde_int' + 'aktiv_anteil' → rel-Bin so benennen
dt_std_rel <- dt_rel[!is.na(stunde_rel),
                     .(aktiv_anteil = mean(aktiv), n_min = .N),
                     by = .(igel, tage_seit, stunde_int = floor(stunde_rel))]
acro_rel_dt <- rbindlist(lapply(igel_liste, function(ig) {
  res <- fit_cosinor(dt_std_rel[igel == ig])
  if (is.null(res)) return(NULL)
  data.table(igel = ig, acrophase_rel_h = round(res$acrophase, 2))
}))

timing_rel <- merge(mu_rel_dt, acro_rel_dt, by = "igel", all = TRUE)
setorder(timing_rel, igel)
cat("── Timing relativ zum Sonnenuntergang (Stunden nach SU) ──\n")
print(timing_rel)
cat(sprintf("\nMedian Acrophase (Cosinor) rel. Sonnenuntergang: +%.2f h (n = %d)\n",
            median(timing_rel$acrophase_rel_h, na.rm = TRUE),
            sum(!is.na(timing_rel$acrophase_rel_h))))
cat(sprintf("Median μ (mittl. Aktivitätszeit) rel. Sonnenuntergang: +%.2f h\n",
            median(timing_rel$mu_rel_h, na.rm = TRUE)))

# speichern (Base-R, keine zusaetzliche Abhaengigkeit)
write.csv(timing_rel,
          file.path(output_ordner, "timing_relative_to_sunset.csv"),
          row.names = FALSE)
cat("✓ timing_relative_to_sunset.csv gespeichert\n\n")

# Plot R1: Rayleigh-Ergebnis — Polarplot (publikationsreif, englisch)
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

# Polarplot auf Sonnenuntergangs-relative μ umstellen (siehe Abschnitt 3b).
# Winkel = Stunden nach Sonnenuntergang (konsistent mit Onset & berichteten Werten).
rayleigh_dt <- merge(rayleigh_dt, mu_rel_dt[, .(igel, mu_rel_h)], by = "igel", all.x = TRUE)
rayleigh_dt[, mu_rad := mu_rel_h * 2 * pi / 24]

# Blinde Tiere identifizieren.
# Hinweis: chrono_meta wird erst weiter unten (Abschnitt Metadaten-Merge) gebaut.
# Daher exists()-Guard, damit der Lauf von oben nicht abbricht; als Fallback die
# bekannten Augen-Pathologie-Tiere (H22, H23, H26), sodass der Rayleigh-Plot sie
# auch dann korrekt einfärbt, wenn chrono_meta hier noch nicht existiert.
blind_igel_vec <- if (exists("chrono_meta") && !is.null(chrono_meta) &&
                      "Diagnosis_group" %in% names(chrono_meta))
  chrono_meta[grepl("(?i)blind", Diagnosis_group), igel] else c("Igel22", "Igel23", "Igel26")

rayleigh_pub <- copy(rayleigh_dt)
rayleigh_pub[, group    := ifelse(igel %in% blind_igel_vec, "Visually impaired", "Sighted")]
# Kürze Label: Igel1 → H1  (platzsparend im Plot)
rayleigh_pub[, label_id := paste0("H", gsub("Igel", "", igel))]
# Blinde Tiere vollständig beschriften, sichtige Tiere nur als Nummer
rayleigh_pub[, plot_label := ifelse(group == "Blind", label_id, label_id)]

p_rayleigh_polar <- ggplot(rayleigh_pub,
                            aes(x      = mu_rad,
                                y      = rho,
                                colour = group,
                                shape  = group,
                                label  = plot_label)) +
  # Nachtlängen-Bandbreite im Messzeitraum (relativ zum Sonnenuntergang):
  #   hellgrau = längste Nacht (~12.2 h nach SU, Ende Sept.)
  #   dunkelgrau = kürzeste Nacht (~7.2 h nach SU, 21. Juni; immer dunkel)
  annotate("rect",
           xmin = 0, xmax = 12.2 * 2*pi/24,
           ymin = 0, ymax = Inf, fill = "grey70", alpha = 0.22) +
  annotate("rect",
           xmin = 0, xmax = 7.2 * 2*pi/24,
           ymin = 0, ymax = Inf, fill = "grey40", alpha = 0.28) +
  # Linien vom Ursprung
  geom_segment(aes(xend = mu_rad, yend = 0),
               linewidth = 0.9, alpha = 0.80) +
  # Punkte an der Spitze
  geom_point(size = 2.8) +
  # Labels — ggrepel verhindert Überlappung
  ggrepel::geom_text_repel(
    size          = 2.3,
    max.overlaps  = 40,
    show.legend   = FALSE,
    segment.size  = 0.3,
    segment.alpha = 0.45,
    box.padding   = 0.20,
    point.padding = 0.05,
    force         = 2.0,
    min.segment.length = 0.2
  ) +
  coord_polar(start = 0) +
  # Achse = Stunden nach Sonnenuntergang (alle 3 h)
  scale_x_continuous(
    limits = c(0, 2 * pi),
    breaks = (0:7) * 2 * pi / 8,
    labels = paste0("+", seq(0, 21, 3), " h"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(0.25, 0.50, 0.75, 1.00),
    labels = c("0.25", "0.50", "0.75", "1.00")
  ) +
  scale_colour_manual(
    values = c("Sighted" = "#2C3E6B", "Visually impaired" = "#C0392B"),
    name   = "Diagnosis"
  ) +
  scale_shape_manual(
    values = c("Sighted" = 19, "Visually impaired" = 17),
    name   = "Diagnosis"
  ) +
  labs(
    x = NULL,
    y = expression(bar(rho) ~ "(mean resultant length)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(size = 9,  colour = "grey30"),
    axis.text.y      = element_text(size = 7,  colour = "grey45"),
    axis.title.y     = element_text(size = 10),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

# Plot R2: ρ̄ und μ als Balkenplot (Supplement / Explorer-Plot)
rayleigh_bar_dt <- merge(rayleigh_dt, rayleigh_pub[, .(igel, group, label_id)],
                         by = "igel", all.x = TRUE)
rayleigh_bar_dt <- rayleigh_bar_dt[order(mu_rel_h)]
rayleigh_bar_dt[, label_f := factor(label_id, levels = label_id)]

p_rayleigh_bar <- ggplot(rayleigh_bar_dt,
                         aes(y = label_f, x = rho, fill = group)) +
  geom_col(alpha = 0.82, width = 0.75) +
  geom_text(aes(x     = rho + 0.01,
                label = sprintf("μ = +%.1f h", mu_rel_h)),
            hjust = 0, size = 2.9, colour = "grey30") +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(
    values = c("Sighted" = "#2C3E6B", "Visually impaired" = "#C0392B"),
    name   = "Diagnosis"
  ) +
  scale_x_continuous(
    limits = c(0, max(rayleigh_bar_dt$rho, na.rm = TRUE) * 1.45),
    labels = function(x) round(x, 2)
  ) +
  labs(
    x = expression(bar(rho) ~ "(mean resultant length)"),
    y = NULL,
    caption = "Labels show mean activity time (μ) in hours after sunset. Sorted by μ."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold"),
    plot.caption    = element_text(colour = "grey50", size = 8)
  )

png(file.path(chrono_ordner, "rayleigh_polar.png"),
    width = 1000, height = 1000, res = 130)
print(p_rayleigh_polar)
dev.off()

png(file.path(chrono_ordner, "rayleigh_balken.png"),
    width = 1400, height = 900, res = 130)
print(p_rayleigh_bar)
dev.off()

cat("✓ Rayleigh-Plots gespeichert\n\n")

# ══════════════════════════════════════════════════════════════
# 4. AKTIVITÄTSBEGINN RELATIV ZU SONNENUNTERGANG
# ══════════════════════════════════════════════════════════════
# Chronobiologisches Maß: Wann beginnt jeder Igel abends
# seine Aktivität relativ zum Sonnenuntergang?
# Positiv = nach Sonnenuntergang aktiv, negativ = vorher.
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("4. Aktivitätsbeginn relativ zu Sonnenuntergang\n")
cat("════════════════════════════════════\n\n")

if (!requireNamespace("suncalc", quietly = TRUE)) install.packages("suncalc")

# Merge Minutendaten × Sonnenzeiten
onset_dt2 <- merge(
  dt[aktiv == 1 & !is.na(stunde)],
  sonnen_dt[, .(datum, untergang_h)],
  by = "datum", all.x = TRUE
)

# Abendliches Suchfenster: Sonnenuntergang −1h bis Mitternacht
onset_dt2 <- onset_dt2[
  !is.na(untergang_h) &
  stunde >= (untergang_h - 1) &
  stunde <  24
]

# Erster aktiver Zeitpunkt pro Igel × Nacht
onset_nacht <- onset_dt2[, .(
  onset_h     = min(stunde, na.rm = TRUE),
  untergang_h = first(untergang_h)
), by = .(igel, datum, tage_seit)]

onset_nacht[, onset_min := (onset_h - untergang_h) * 60]
onset_nacht <- onset_nacht[onset_min >= -60]  # Ausreißer: > 1h vor SU entfernen

# Zusammenfassung pro Igel
onset_summary <- onset_nacht[, .(
  med_onset_min = median(onset_min, na.rm = TRUE),
  q25_onset     = quantile(onset_min, 0.25, na.rm = TRUE),
  q75_onset     = quantile(onset_min, 0.75, na.rm = TRUE),
  n_naechte     = .N
), by = igel]
onset_summary <- onset_summary[order(med_onset_min)]
onset_summary[, igel_f := factor(igel, levels = igel)]

cat("Aktivitätsbeginn (Median Minuten nach Sonnenuntergang):\n")
cat("  Populationsmedian:", round(median(onset_nacht$onset_min, na.rm=TRUE), 0), "min\n")
cat("  Bereich:", round(min(onset_summary$med_onset_min),0),
    "bis", round(max(onset_summary$med_onset_min),0), "min\n\n")

# ── Plot O1: Boxplot Onset pro Igel ──────────────────────────
# Merge für Boxplot: onset_nacht braucht igel_f
onset_nacht_plot <- merge(onset_nacht,
                          onset_summary[, .(igel, igel_f)],
                          by = "igel", all.x = TRUE)

p_onset_box <- ggplot(onset_nacht_plot, aes(x = igel_f, y = onset_min)) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 1) +
  annotate("text", x = Inf, y = 4,
           label = "↑ nach SU", size = 2.8, color = "grey40",
           hjust = 1.1, fontface = "italic") +
  geom_boxplot(fill = "#A8D8E8", color = "#2C3E6B",
               outlier.size = 1, outlier.alpha = 0.5, width = 0.6) +
  geom_point(data = onset_summary,
             aes(x = igel_f, y = med_onset_min),
             shape = 18, size = 3.5, color = "#C0392B",
             inherit.aes = FALSE) +
  scale_y_continuous(labels = function(x) paste0(ifelse(x >= 0, "+", ""), x, " min")) +
  labs(
    title    = "Aktivitätsbeginn: Wann werden die Igel abends aktiv?",
    subtitle = paste0(
      "Minuten nach Sonnenuntergang (0 = exakt SU) | Raute = Median\n",
      "Sortiert nach medianer Onset-Zeit | Boxplot = IQR ± 1.5×IQR | n = Beobachtungsnächte"),
    x = NULL, y = "Aktivitätsbeginn (min nach Sonnenuntergang)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

png(file.path(chrono_ordner, "onset_boxplot.png"), width = 1600, height = 800, res = 130)
print(p_onset_box)
dev.off()

# ── Plot O2: Onset über Zeit (erste 14 Tage, Populationstrend) ──
onset_pop_tag <- onset_nacht[tage_seit <= 14, .(
  med_min = median(onset_min, na.rm = TRUE),
  q25     = quantile(onset_min, 0.25, na.rm = TRUE),
  q75     = quantile(onset_min, 0.75, na.rm = TRUE),
  n       = .N
), by = tage_seit]

p_onset_zeit <- ggplot() +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.8) +
  geom_jitter(data = onset_nacht[tage_seit <= 14],
              aes(x = tage_seit, y = onset_min, color = igel),
              width = 0.2, alpha = 0.35, size = 1.8) +
  geom_ribbon(data = onset_pop_tag,
              aes(x = tage_seit, ymin = q25, ymax = q75),
              fill = "black", alpha = 0.12, inherit.aes = FALSE) +
  geom_line(data = onset_pop_tag,
            aes(x = tage_seit, y = med_min),
            color = "black", linewidth = 1.4, inherit.aes = FALSE) +
  geom_text(data = onset_pop_tag,
            aes(x = tage_seit,
                y = min(onset_nacht[tage_seit <= 14, onset_min], na.rm=TRUE) - 15,
                label = n),
            size = 2.3, color = "grey55") +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  scale_x_continuous(breaks = c(1, 3, 7, 14)) +
  scale_y_continuous(labels = function(x) paste0(ifelse(x >= 0, "+", ""), x, " min")) +
  labs(
    title    = "Aktivitätsbeginn über die ersten 14 Tage nach Auswilderung",
    subtitle = "Minuten nach Sonnenuntergang | Schwarz = Median ± IQR | N = Igel-Nächte",
    x = "Tage nach Auswilderung",
    y = "Aktivitätsbeginn (min nach Sonnenuntergang)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

png(file.path(chrono_ordner, "onset_zeitverlauf.png"), width = 1500, height = 800, res = 130)
print(p_onset_zeit)
dev.off()

cat("✓ Onset-Plots gespeichert\n\n")

# ══════════════════════════════════════════════════════════════
# 5. METADATEN-VERKNÜPFUNG
# ══════════════════════════════════════════════════════════════
# Frage: Unterscheiden sich chronobiologische Parameter zwischen
# Diagnosegruppen, Jahreszeiten, Altersklassen, Geschlechtern?
# Zentrales Ergebnis: Diagnose "Blindness" → gestörte Rhythmizität
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("5. Metadaten-Verknüpfung\n")
cat("════════════════════════════════════\n\n")

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

meta_pfad <- file.path(projekt_root, "output", "Block1_Metadaten", "Block1_Uebersicht.xlsx")
if (!file.exists(meta_pfad)) {
  cat("⚠ Block1_Uebersicht.xlsx nicht gefunden — Metadaten-Analyse übersprungen\n\n")
  meta_dt     <- NULL
  chrono_meta <- NULL
} else {
  meta_raw <- read.xlsx(meta_pfad, sheet = "All animals")
  setDT(meta_raw)
  setnames(meta_raw, names(meta_raw),
           gsub("[^A-Za-z0-9_]", "_", names(meta_raw)))
  # Igel-Spalte
  igel_col <- names(meta_raw)[grepl("(?i)animal.*id|igel", names(meta_raw))][1]
  setnames(meta_raw, igel_col, "igel")
  meta_raw[, igel := as.character(igel)]

  # Relevante Spalten auswählen
  wanted <- c("igel", "Sex", "Age_class", "Diagnosis_group",
              "Season", "Rel__weight__g_", "Rehab__d_", "Tracking__d_")
  wanted <- intersect(wanted, names(meta_raw))
  meta_dt <- meta_raw[, ..wanted]

  # Chrono-Parameter zusammenführen
  chrono_merge <- merge(cosinor_gesamt_dt, rayleigh_dt, by = "igel", all = TRUE)
  chrono_merge <- merge(chrono_merge,
                        onset_summary[, .(igel, med_onset_min, n_naechte)],
                        by = "igel", all.x = TRUE)
  chrono_meta  <- merge(chrono_merge, meta_dt, by = "igel", all.x = TRUE)

  # Fallback: 'ist_blind' IMMER als Spalte definieren, damit die Metadaten-Plots
  # (Jahreszeit, Alter/Geschlecht) nicht abstuerzen, falls keine passende
  # Diagnosespalte im Merge vorhanden ist (Spaltenname-Mismatch).
  .diag_sp <- intersect(c("Diagnosis_group", "Diagnosis", "diagnosis_main", "diagnosis"),
                        names(chrono_meta))
  if (length(.diag_sp) > 0) {
    chrono_meta[, ist_blind := grepl("(?i)blind", get(.diag_sp[1]))]
  } else {
    chrono_meta[, ist_blind := igel %in% c("Igel22", "Igel23", "Igel26")]
  }

  cat("Metadaten geladen:", nrow(meta_dt), "Igel\n")
  cat("Merged chrono × meta:", nrow(chrono_meta), "Zeilen\n\n")

  # Diagnosegruppen-Ausgabe
  if ("Diagnosis_group" %in% names(chrono_meta)) {
    cat("── Rayleigh ρ̄ nach Diagnosegruppe ──\n")
    print(chrono_meta[, .(
      N      = .N,
      med_rho     = round(median(rho,      na.rm=TRUE), 3),
      med_acro    = round(median(acrophase, na.rm=TRUE), 1),
      med_onset   = round(median(med_onset_min, na.rm=TRUE), 0)
    ), by = Diagnosis_group][order(-med_rho)])
    cat("\n")
  }

  # ── Plot M1: ρ̄ nach Diagnosegruppe (Punkte + Boxplot) ───────
  diag_col <- "Diagnosis_group"
  if (diag_col %in% names(chrono_meta) && !all(is.na(chrono_meta[[diag_col]]))) {

    # Blinde Igel hervorheben
    chrono_meta[, ist_blind := grepl("(?i)blind", get(diag_col))]

    p_meta_diag <- ggplot(chrono_meta[!is.na(get(diag_col))],
                          aes_string(x = diag_col, y = "rho",
                                     color = "ist_blind", label = "igel")) +
      geom_boxplot(aes_string(group = diag_col), fill = "grey92",
                   color = "grey60", outlier.shape = NA, width = 0.5) +
      geom_jitter(size = 3.5, width = 0.15, alpha = 0.85) +
      {if (requireNamespace("ggrepel", quietly=TRUE))
         ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20, show.legend = FALSE)
       else
         geom_text(vjust = -0.9, size = 2.8, show.legend = FALSE)} +
      geom_hline(yintercept = median(chrono_meta$rho, na.rm=TRUE),
                 linetype = "dashed", color = "grey40") +
      scale_color_manual(
        values  = c("FALSE" = "#2C3E6B", "TRUE" = "#C0392B"),
        labels  = c("FALSE" = "Other diagnosis", "TRUE" = "Visually impaired"),
        name    = ""
      ) +
      labs(
        title    = "Rhythmusstärke (ρ̄) nach Diagnosegruppe",
        subtitle = paste0("Rayleigh ρ̄: 0 = arrhythmisch, 1 = perfekter Rhythmus\n",
                          "Gestrichelt = Populationsmedian | Rot = Blindheitsdiagnose"),
        x = "Diagnosegruppe", y = expression(bar(rho) ~ "(Rhythmusstärke)")
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title      = element_text(face = "bold"),
            legend.position = "bottom",
            axis.text.x     = element_text(angle = 25, hjust = 1))

    png(file.path(chrono_ordner, "meta_rho_diagnose.png"),
        width = 1400, height = 900, res = 130)
    print(p_meta_diag)
    dev.off()
    cat("✓ Plot: ρ̄ nach Diagnose\n")
  }

  # ── Plot M2: Akrophase nach Jahreszeit ───────────────────────
  if ("Season" %in% names(chrono_meta)) {

    p_meta_season <- ggplot(chrono_meta[!is.na(Season) & !is.na(acrophase)],
                            aes(x = Season, y = acrophase,
                                color = ist_blind, label = igel)) +
      geom_boxplot(aes(group = Season), fill = "grey92",
                   color = "grey60", outlier.shape = NA, width = 0.5) +
      geom_jitter(size = 3.5, width = 0.15, alpha = 0.85) +
      {if (requireNamespace("ggrepel", quietly=TRUE))
         ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20, show.legend = FALSE)
       else
         geom_text(vjust = -0.9, size = 2.8, show.legend = FALSE)} +
      # Typische Nachtphase
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = 21, ymax = 24, fill = "#1a2744", alpha = 0.08) +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin =  0, ymax =  5, fill = "#1a2744", alpha = 0.08) +
      scale_y_continuous(breaks = seq(0, 24, 3),
                         labels = sprintf("%02d:00", seq(0, 24, 3)),
                         limits = c(0, 24)) +
      scale_color_manual(
        values = c("FALSE" = "#2C3E6B", "TRUE" = "#C0392B"),
        labels = c("FALSE" = "Andere Diagnose", "TRUE" = "Blind"), name = ""
      ) +
      labs(
        title    = "Cosinor-Akrophase (Aktivitätsgipfel) nach Jahreszeit",
        subtitle = "Blauer Bereich = typische Nachtphase (21:00–05:00) | Rot = Blindheitsdiagnose",
        x = "Jahreszeit der Auswilderung", y = "Akrophase (Uhrzeit)"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom")

    png(file.path(chrono_ordner, "meta_akrophase_saison.png"),
        width = 1300, height = 900, res = 130)
    print(p_meta_season)
    dev.off()
    cat("✓ Plot: Akrophase nach Jahreszeit\n")
  }

  # ── Plot M3: ρ̄ nach Alter und Geschlecht ────────────────────
  if (all(c("Age_class", "Sex") %in% names(chrono_meta))) {

    p_meta_age <- ggplot(chrono_meta[!is.na(Age_class) & !is.na(Sex)],
                         aes(x = Age_class, y = rho,
                             color = ist_blind, shape = Sex, label = igel)) +
      geom_boxplot(aes(group = Age_class), fill = "grey92",
                   color = "grey60", outlier.shape = NA, width = 0.5) +
      geom_jitter(size = 3.5, width = 0.15, alpha = 0.85) +
      {if (requireNamespace("ggrepel", quietly=TRUE))
         ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20, show.legend = FALSE)
       else
         geom_text(vjust = -0.9, size = 2.8, show.legend = FALSE)} +
      scale_color_manual(
        values = c("FALSE" = "#2C3E6B", "TRUE" = "#C0392B"),
        labels = c("FALSE" = "Other diagnosis", "TRUE" = "Visually impaired"), name = "Diagnosis"
      ) +
      scale_shape_manual(values = c("Male" = 16, "Female" = 17),
                         name = "Geschlecht") +
      labs(
        title    = "Rhythmusstärke (ρ̄) nach Altersklasse und Geschlecht",
        subtitle = "Kreis = männlich, Dreieck = weiblich | Rot = Blindheitsdiagnose",
        x = "Altersklasse", y = expression(bar(rho) ~ "(Rhythmusstärke)")
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom")

    png(file.path(chrono_ordner, "meta_rho_alter_sex.png"),
        width = 1200, height = 900, res = 130)
    print(p_meta_age)
    dev.off()
    cat("✓ Plot: ρ̄ nach Alter/Geschlecht\n")
  }

  cat("\n")
}

# ══════════════════════════════════════════════════════════════
# 6. ZUSAMMENFASSUNG & WORD-BERICHT
# ══════════════════════════════════════════════════════════════

cat("════════════════════════════════════\n")
cat("Zusammenfassung\n")
cat("════════════════════════════════════\n\n")

cat("── Cosinor Gesamt-Übersicht ──\n")
cosinor_sig <- cosinor_gesamt_dt[p_cosinor < 0.05]
cat("Signifikante Cosinor-Fits (p < 0.05):",
    nrow(cosinor_sig), "von", nrow(cosinor_gesamt_dt), "Igeln\n")
if (nrow(cosinor_sig) > 0) {
  cat("Mittlere Akrophase (signifikante Fits):",
      round(mean(cosinor_sig$acrophase), 1), "Uhr\n")
  cat("Mittlere Amplitude:", round(mean(cosinor_sig$amplitude * 100), 1), "%\n")
}
cat("\n")

# Ergebnisse speichern
saveRDS(list(
  cosinor_rollend  = cosinor_dt,
  cosinor_gesamt   = cosinor_gesamt_dt,
  rayleigh         = rayleigh_dt,
  onset_summary    = onset_summary,
  chrono_meta      = chrono_meta
), file.path(chrono_ordner, "chrono_ergebnisse.rds"))

cat("✓ Ergebnisse gespeichert: Chronobiologie/chrono_ergebnisse.rds\n\n")

# ── Word-Bericht ───────────────────────────────────────────────
cat("Erstelle Word-Bericht...\n")

doc <- read_docx() |>
  body_set_default_section(prop_section(
    page_margins = page_mar(top = 2, bottom = 2, left = 3, right = 2.5),
    page_size    = page_size(width = 21/2.54, height = 29.7/2.54)
  ))

# Titel
doc <- doc |>
  body_add_par("Chronobiologische Analyse — Igelbesenderung Niedersachsen",
               style = "heading 1") |>
  body_add_par(paste("Erstellt:", format(Sys.time(), "%d.%m.%Y %H:%M")),
               style = "Normal") |>
  body_add_par(paste0(
    "Methoden: Actogramme (Double-Plot), Cosinor-Analyse (Mesor/Amplitude/Akrophase), ",
    "Rayleigh-Test (zirkuläre Statistik). Analysen pro Individuum — keine Poolung ",
    "über Individuen vor der Parameterextraktion."),
    style = "Normal") |>
  body_add_break()

# 1. Actogramme
doc <- body_add_par(doc, "1  Actogramme", style = "heading 1")
doc <- body_add_par(doc, paste0(
  "Doppelt aufgetragene Actogramme (Double-Plot) zeigen die 24h-Aktivität jedes ",
  "Igels über alle Beobachtungstage. Jede Zeile entspricht 48h (doppelt dargestellt), ",
  "sodass eine eventuelle Phasendrift des Rhythmus als Diagonale erkennbar wird. ",
  "Orange Linie = Sonnenuntergang, blau = Sonnenaufgang."),
  style = "Normal")
doc <- body_add_par(doc, "Übersicht — alle Igel, erste 21 Tage", style = "heading 2")
doc <- add_img_safe(doc, file.path(chrono_ordner, "actogramme_alle_uebersicht.png"),
                    b = 18, h = 14)
doc <- body_add_par(doc,
  "Einzelne Actogramme für jeden Igel sind im Ordner Chronobiologie/ gespeichert.",
  style = "Normal")
doc <- body_add_break(doc)

# 2. Cosinor
doc <- body_add_par(doc, "2  Cosinor-Analyse", style = "heading 1")
doc <- body_add_par(doc, paste0(
  "Die Cosinor-Analyse (Halberg et al. 1967) fittet pro Individuum eine Cosinus-Kurve ",
  "an die stündlichen Aktivitätsmittelwerte: y = Mesor + Amplitude × cos(2πt/24 − Akrophase). ",
  "Mesor = mittleres Aktivitätsniveau; Amplitude = Stärke des 24h-Rhythmus; ",
  "Akrophase = Uhrzeit des Aktivitätsgipfels. Nur Fits mit p < 0.05 (F-Test) ",
  "werden als signifikant gewertet."),
  style = "Normal")

# Tabelle Gesamt-Cosinor
tbl_cos <- cosinor_gesamt_dt[, .(
  Igel       = igel,
  Mesor      = paste0(round(mesor*100, 1), "%"),
  Amplitude  = paste0(round(amplitude*100, 1), "%"),
  Akrophase  = sprintf("%05.2f Uhr", acrophase),
  R2         = round(r2, 3),
  p          = ifelse(p_cosinor < 0.001, "< 0.001",
               ifelse(p_cosinor < 0.01,  "< 0.01",
               ifelse(p_cosinor < 0.05,  "< 0.05",
                      as.character(round(p_cosinor, 3))))),
  signifikant = ifelse(p_cosinor < 0.05, "✓", "")
)]
ft_cos <- flextable(as.data.frame(tbl_cos)) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#2C3E6B") |>
  color(part = "header", color = "white") |>
  bg(i = which(tbl_cos$signifikant == "✓"), bg = "#EAF4FB") |>
  bg(i = seq(1, nrow(tbl_cos), 2), bg = "#F8F8F8") |>
  fontsize(size = 9, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  autofit()
doc <- body_add_par(doc, "Cosinor-Parameter pro Igel (Gesamtbeobachtungszeitraum)",
                    style = "heading 2")
doc <- body_add_flextable(doc, ft_cos)

doc <- body_add_par(doc, "Akrophase — Polarplot", style = "heading 2")
doc <- add_img_safe(doc, file.path(chrono_ordner, "cosinor_akrophase_polar.png"), b=12, h=12)

doc <- body_add_par(doc, "Gefittete 24h-Profile aller Igel", style = "heading 2")
doc <- add_img_safe(doc, file.path(chrono_ordner, "cosinor_profile_alle.png"), b=16, h=9)

if (file.exists(file.path(chrono_ordner, "cosinor_akrophase_zeitverlauf.png"))) {
  doc <- body_add_par(doc, "Akrophase über Zeit (rollendes Fenster)", style = "heading 2")
  doc <- body_add_par(doc, paste0(
    "Zeigt ob sich der Zeitpunkt des Aktivitätsgipfels über die Wochen nach der ",
    "Auswilderung verschiebt. Jede Linie = ein Igel. Nur signifikante Fenster (p < 0.05)."),
    style = "Normal")
  doc <- add_img_safe(doc, file.path(chrono_ordner, "cosinor_akrophase_zeitverlauf.png"),
                      b=16, h=9)
  doc <- body_add_par(doc, "Amplitude über Zeit", style = "heading 2")
  doc <- add_img_safe(doc, file.path(chrono_ordner, "cosinor_amplitude_zeitverlauf.png"),
                      b=16, h=9)
}
doc <- body_add_break(doc)

# 3. Rayleigh
doc <- body_add_par(doc, "3  Rayleigh-Test", style = "heading 1")
doc <- body_add_par(doc, paste0(
  "Der Rayleigh-Test (Rayleigh 1880, Fisher 1993) prüft mit zirkulärer Statistik, ",
  "ob die Aktivitätszeitpunkte eines Individuums signifikant auf eine Tageszeit ",
  "konzentriert sind (H0: gleichmäßige Verteilung = arrhythmisch). ",
  "ρ̄ (mean resultant length) misst die Konzentration (0 = gleichmäßig, 1 = alle ",
  "Aktivität zu einem Zeitpunkt). μ ist die mittlere Aktivitätszeit."),
  style = "Normal")

tbl_ray <- rayleigh_dt[, .(
  Igel                   = igel,
  `N aktive Min.`        = format(n_aktiv_min, big.mark = "'"),
  `ρ̄`                   = rho,
  `μ (mittl. Zeit)`      = sprintf("%05.2f Uhr", mu_h),
  `p (Rayleigh)`         = ifelse(rayleigh_p < 0.001, "< 0.001",
                           ifelse(rayleigh_p < 0.01,  "< 0.01",
                           ifelse(rayleigh_p < 0.05,  "< 0.05",
                                  as.character(rayleigh_p)))),
  `Rhythmisch?`          = ifelse(rhythmisch, "✓ ja", "✗ nein")
)]
ft_ray <- flextable(as.data.frame(tbl_ray)) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#2C3E6B") |>
  color(part = "header", color = "white") |>
  bg(i = which(tbl_ray$`Rhythmisch?` == "✓ ja"), bg = "#EAF4FB") |>
  bg(i = seq(1, nrow(tbl_ray), 2), bg = "#F8F8F8") |>
  fontsize(size = 9, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  autofit()
doc <- body_add_flextable(doc, ft_ray)

doc <- body_add_par(doc, "Rayleigh-Test — Polarplot", style = "heading 2")
doc <- add_img_safe(doc, file.path(chrono_ordner, "rayleigh_polar.png"), b=12, h=12)
doc <- body_add_par(doc, "Rayleigh-Test — Balkendiagramm", style = "heading 2")
doc <- add_img_safe(doc, file.path(chrono_ordner, "rayleigh_balken.png"), b=16, h=9)
doc <- body_add_break(doc)

# 4. Onset
doc <- body_add_par(doc, "4  Aktivitätsbeginn relativ zum Sonnenuntergang", style = "heading 1")
doc <- body_add_par(doc, paste0(
  "Für jede Nacht und jeden Igel wurde der erste Zeitpunkt ermittelt, zu dem das Tier ",
  "im abendlichen Suchfenster (Sonnenuntergang −1h bis Mitternacht) aktiv war. ",
  "Das Ergebnis wird in Minuten nach Sonnenuntergang ausgedrückt (0 = exakt Sonnenuntergang, ",
  "positiv = später, negativ = früher). Werte unter −60 min wurden als Ausreißer entfernt. ",
  "Gesunder Igelrhythmus: Aktivitätsbeginn kurz nach Sonnenuntergang erwartet."),
  style = "Normal")

if (nrow(onset_summary) > 0) {
  # Onset-Tabelle
  tbl_onset <- onset_summary[, .(
    Igel            = igel,
    `Median (min)`  = round(med_onset_min, 0),
    `Q25 (min)`     = round(q25_onset, 0),
    `Q75 (min)`     = round(q75_onset, 0),
    `N Nächte`      = n_naechte
  )]
  ft_onset <- flextable(as.data.frame(tbl_onset)) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#2C3E6B") |>
    color(part = "header", color = "white") |>
    bg(i = seq(1, nrow(tbl_onset), 2), bg = "#F8F8F8") |>
    fontsize(size = 9, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    autofit()
  doc <- body_add_par(doc, "Onset-Zeiten je Igel (sortiert nach Median)", style = "heading 2")
  doc <- body_add_flextable(doc, ft_onset)
}
doc <- body_add_par(doc, "Onset-Boxplot je Igel", style = "heading 2")
doc <- add_img_safe(doc, file.path(chrono_ordner, "onset_boxplot.png"), b=17, h=8)
doc <- body_add_par(doc, "Onset-Verlauf (erste 14 Tage)", style = "heading 2")
doc <- add_img_safe(doc, file.path(chrono_ordner, "onset_zeitverlauf.png"), b=16, h=8)
doc <- body_add_break(doc)

# 5. Metadaten
doc <- body_add_par(doc, "5  Metadaten-Verknüpfung", style = "heading 1")
doc <- body_add_par(doc, paste0(
  "Die chronobiologischen Parameter (ρ̄, Akrophase, Onset) wurden mit den ",
  "Metadaten aus Block 1 (Diagnosegruppe, Jahreszeit, Altersklasse, Geschlecht) verknüpft. ",
  "Zentrales Ergebnis: Drei Igel mit der Diagnose 'Blindness' (Igel 22, 23, 26) zeigen ",
  "deutlich abweichende Rhythmizität. Blinde Säugetiere können ihren circadianen Rhythmus ",
  "nicht über den Hell-Dunkel-Zyklus synchronisieren (fehlender Lichtzeitgeber), was zu ",
  "schwächerer oder zeitlich verschobener Rhythmizität führt — ein in der chronobiologischen ",
  "Literatur bekanntes Phänomen (sogenannte 'free-running rhythms')."),
  style = "Normal")

if (file.exists(file.path(chrono_ordner, "meta_rho_diagnose.png"))) {
  doc <- body_add_par(doc, "Rhythmusstärke (ρ̄) nach Diagnosegruppe", style = "heading 2")
  doc <- body_add_par(doc, paste0(
    "Igel 22 und 23 (Blindness) zeigen die niedrigsten ρ̄-Werte (0.309 bzw. 0.417) ",
    "im gesamten Datensatz. Igel 26 (ebenfalls blind) hat weniger Beobachtungsnächte ",
    "und zeigt noch eine relativ normale Rhythmizität. Die Diagnosegruppe 'Blindness' ",
    "weist systematisch niedrigere ρ̄-Werte auf als alle anderen Gruppen."),
    style = "Normal")
  doc <- add_img_safe(doc, file.path(chrono_ordner, "meta_rho_diagnose.png"), b=15, h=9)
}

if (file.exists(file.path(chrono_ordner, "meta_akrophase_saison.png"))) {
  doc <- body_add_par(doc, "Akrophase nach Jahreszeit der Auswilderung", style = "heading 2")
  doc <- body_add_par(doc, paste0(
    "Die Hauptgruppe der Igel zeigt unabhängig von der Jahreszeit der Auswilderung ",
    "Aktivitätsgipfel in der typischen Nachtphase (21:00–05:00 Uhr). Saisonale ",
    "Unterschiede im Sonnenuntergang (z.B. Sommer: ~21:30 Uhr, Herbst: ~18:30 Uhr) ",
    "könnten die absolute Akrophase leicht verschieben. Die beiden blinden Ausreißer ",
    "(Igel 22, 23) stammen beide aus der Sommerkohorte."),
    style = "Normal")
  doc <- add_img_safe(doc, file.path(chrono_ordner, "meta_akrophase_saison.png"), b=14, h=9)
}

if (file.exists(file.path(chrono_ordner, "meta_rho_alter_sex.png"))) {
  doc <- body_add_par(doc, "Rhythmusstärke nach Altersklasse und Geschlecht", style = "heading 2")
  doc <- body_add_par(doc, paste0(
    "Erwachsene und juvenile Igel zeigen vergleichbare Rhythmusstärken. ",
    "Systematische Unterschiede zwischen den Geschlechtern sind nicht erkennbar. ",
    "Die Blindheitsdiagnose ist der stärkste Prädiktor für reduzierte Rhythmizität."),
    style = "Normal")
  doc <- add_img_safe(doc, file.path(chrono_ordner, "meta_rho_alter_sex.png"), b=13, h=9)
}
doc <- body_add_break(doc)

# Speichern
bericht_pfad <- file.path(output_ordner, "Block3_Chronobiologie_Bericht.docx")
print(doc, target = bericht_pfad)

cat("✓ Word-Bericht gespeichert:", basename(bericht_pfad), "\n\n")

# ══════════════════════════════════════════════════════════════
# 7. EXCEL-ÜBERSICHT (automatisch, wird bei jedem Lauf neu erzeugt)
# ══════════════════════════════════════════════════════════════
cat("Erstelle Excel-Übersicht...\n")

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

wb <- createWorkbook()

# ── Hilfsfunktion: Header-Stil ──────────────────────────────
stil_header <- createStyle(
  fgFill = "#2C3E6B", fontColour = "white",
  textDecoration = "bold", halign = "center",
  border = "Bottom", borderColour = "#AAAAAA"
)
stil_grau <- createStyle(fgFill = "#F2F2F2")
stil_blau <- createStyle(fgFill = "#EAF4FB")
stil_zahl <- createStyle(numFmt = "0.000", halign = "right")
stil_pct  <- createStyle(numFmt = "0.0%",  halign = "right")

# ────────────────────────────────────────────────────────────
# Sheet 1: Kennzahlen (Übersicht)
# ────────────────────────────────────────────────────────────
addWorksheet(wb, "Kennzahlen")

n_chrono_gesamt <- length(igel_liste_alle)  # alle Tiere im chrono-Datensatz
n_ausgeschlossen <- length(igel_excl)
n_eingeschlossen <- length(igel_liste)      # nach >= MIN_AKTIV_MIN-Filter
n_cosinor_sig    <- sum(cosinor_gesamt_dt$p_cosinor < 0.05, na.rm = TRUE)
n_rayleigh_sig   <- sum(rayleigh_dt$rayleigh_p < 0.001, na.rm = TRUE)
rho_median       <- round(median(rayleigh_dt$rho, na.rm = TRUE), 3)
rho_min          <- round(min(rayleigh_dt$rho, na.rm = TRUE), 3)
rho_max          <- round(max(rayleigh_dt$rho, na.rm = TRUE), 3)

kennzahlen <- data.frame(
  Kennzahl = c(
    paste0("Block 3 — Chronobiologie  |  VHF-Igel-Monitoring  |  TiHo Hannover"),
    paste0("Erstellt: ", format(Sys.time(), "%d.%m.%Y %H:%M"),
           " | Einschluss-Schwellenwert: >=", MIN_AKTIV_MIN, " aktive Nachtminuten"),
    "Kennzahl", # Header-Zeile
    "Tiere im Chrono-Datensatz (gesamt)",
    paste0("  davon ausgeschlossen (<", MIN_AKTIV_MIN, " aktive Min.): ",
           paste(igel_excl, collapse = ", ")),
    "Eingeschlossene Tiere (Analyse)",
    "Igel mit signifikantem Cosinor-Fit (p<0.05)",
    "Igel rhythmisch laut Rayleigh (p<0.001)",
    "Median ρ̄ (Rhythmusstärke)",
    "Median Akrophase (Aktivitätsgipfel)",
    "Median Aktivitätsbeginn nach SU",
    if ("Diagnosis_group" %in% names(chrono_meta) && !is.null(chrono_meta))
      paste0("Igel mit Diagnose 'Blindness'") else character(0),
    "⚠  Wichtiger Befund"
  ),
  Wert = c(
    "", "", "Wert",
    as.character(n_chrono_gesamt),
    as.character(n_ausgeschlossen),
    as.character(n_eingeschlossen),
    paste0(n_cosinor_sig, " / ", n_eingeschlossen),
    paste0(n_rayleigh_sig, " / ", n_eingeschlossen),
    as.character(rho_median),
    {mu_pop <- rayleigh_dt[, median(mu_h, na.rm = TRUE)]; sprintf("%02d:%02d Uhr",
      floor(mu_pop), round((mu_pop - floor(mu_pop)) * 60))},
    {if (nrow(onset_summary) > 0)
      paste0("+", round(median(onset_summary$med_onset_min, na.rm=TRUE), 0), " min")
      else "—"},
    if ("Diagnosis_group" %in% names(chrono_meta) && !is.null(chrono_meta)) {
      blind_igel <- chrono_meta[grepl("(?i)blind", Diagnosis_group), igel]
      paste0(length(blind_igel), " (", paste(blind_igel, collapse = ", "), ")")
    } else character(0),
    "Blinde Igel (Igel 22, 23) zeigen schwache/verschobene Rhythmizität (free-running rhythms)"
  ),
  `Einheit / Anmerkung` = c(
    "", "", "Einheit / Anmerkung",
    paste0("chrono_minuten24h.rds"),
    paste0("<", MIN_AKTIV_MIN, " aktive Nachtminuten"),
    "mit Rayleigh-Test und Cosinor",
    paste0(round(n_cosinor_sig / n_eingeschlossen * 100, 1), "%"),
    "100% — alle eingeschlossenen Igel",
    paste0("Range: ", rho_min, " – ", rho_max),
    "Hauptgruppe 22:00–02:00",
    "Population (SU = Sonnenuntergang)",
    if ("Diagnosis_group" %in% names(chrono_meta) && !is.null(chrono_meta))
      "Alle zeigen abweichende Muster" else character(0),
    "Blindheit = stärkster Prädiktor für gestörte Rhythmizität"
  ),
  stringsAsFactors = FALSE
)

writeData(wb, "Kennzahlen", kennzahlen, startRow = 1, colNames = FALSE)
addStyle(wb, "Kennzahlen", stil_header, rows = 3, cols = 1:3, gridExpand = TRUE)
setColWidths(wb, "Kennzahlen", cols = 1:3, widths = c(55, 30, 45))
mergeCells(wb, "Kennzahlen", cols = 1:3, rows = 1)
mergeCells(wb, "Kennzahlen", cols = 1:3, rows = 2)
addStyle(wb, "Kennzahlen",
  createStyle(textDecoration = "bold", fontSize = 13), rows = 1, cols = 1)

# ────────────────────────────────────────────────────────────
# Sheet 2: Cosinor-Ergebnisse
# ────────────────────────────────────────────────────────────
addWorksheet(wb, "Cosinor")

cos_export <- cosinor_gesamt_dt[, .(
  Igel        = igel,
  `Mesor (%)`      = round(mesor * 100, 2),
  `Amplitude (%)`  = round(amplitude * 100, 2),
  `Akrophase (h)`  = round(acrophase, 2),
  `R²`        = round(r2, 4),
  `p (F-Test)`     = round(p_cosinor, 4),
  `Signifikant p<0.05` = ifelse(p_cosinor < 0.05, "Ja", "Nein")
)]

title_cos <- data.frame(
  V1 = "Cosinor-Analyse — Zirkulare Rhythmuskennzahlen pro Igel",
  stringsAsFactors = FALSE)
writeData(wb, "Cosinor", title_cos, startRow = 1, colNames = FALSE)
addStyle(wb, "Cosinor",
  createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
writeData(wb, "Cosinor", as.data.frame(cos_export), startRow = 3, colNames = TRUE)
addStyle(wb, "Cosinor", stil_header, rows = 3, cols = 1:7, gridExpand = TRUE)
# Abwechselnde Zeilenfärbung + blau für signifikante
for (i in seq_len(nrow(cos_export))) {
  r <- i + 3
  if (cos_export$`Signifikant p<0.05`[i] == "Ja") {
    addStyle(wb, "Cosinor", stil_blau, rows = r, cols = 1:7, gridExpand = TRUE, stack = TRUE)
  } else if (i %% 2 == 0) {
    addStyle(wb, "Cosinor", stil_grau, rows = r, cols = 1:7, gridExpand = TRUE, stack = TRUE)
  }
}
setColWidths(wb, "Cosinor", cols = 1:7, widths = c(10, 13, 15, 15, 8, 12, 18))
freezePane(wb, "Cosinor", firstRow = TRUE, firstActiveRow = 4)

# ────────────────────────────────────────────────────────────
# Sheet 3: Rayleigh-Test
# ────────────────────────────────────────────────────────────
addWorksheet(wb, "Rayleigh")

ray_export <- rayleigh_dt[, .(
  Igel                   = igel,
  `N aktive Min.`        = n_aktiv_min,
  `ρ̄ (Rhythmusstärke)` = rho,
  `μ (mittl. Aktivitätszeit)` = sprintf("%05.2f Uhr", mu_h),
  `p (Rayleigh)`         = ifelse(rayleigh_p == 0, "< 0.0001",
                            ifelse(rayleigh_p < 0.001, "< 0.001",
                            ifelse(rayleigh_p < 0.01, "< 0.01",
                                   as.character(round(rayleigh_p, 4))))),
  `Rhythmisch? (p<0.001)` = ifelse(rayleigh_p < 0.001, "Ja", "Nein")
)]

title_ray <- data.frame(
  V1 = "Rayleigh-Test — Zirkuläre Statistik (Rhythmizität pro Igel)",
  stringsAsFactors = FALSE)
writeData(wb, "Rayleigh", title_ray, startRow = 1, colNames = FALSE)
addStyle(wb, "Rayleigh",
  createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
writeData(wb, "Rayleigh", as.data.frame(ray_export), startRow = 3, colNames = TRUE)
addStyle(wb, "Rayleigh", stil_header, rows = 3, cols = 1:6, gridExpand = TRUE)
for (i in seq_len(nrow(ray_export))) {
  r <- i + 3
  if (ray_export$`Rhythmisch? (p<0.001)`[i] == "Ja") {
    addStyle(wb, "Rayleigh", stil_blau, rows = r, cols = 1:6, gridExpand = TRUE, stack = TRUE)
  } else if (i %% 2 == 0) {
    addStyle(wb, "Rayleigh", stil_grau, rows = r, cols = 1:6, gridExpand = TRUE, stack = TRUE)
  }
}
# Legende
leg_row <- nrow(ray_export) + 5
writeData(wb, "Rayleigh",
  data.frame(Legende = c(
    "ρ̄ ≥ 0.80: stark rhythmisch",
    "ρ̄ 0.60–0.79: rhythmisch",
    "ρ̄ 0.40–0.59: schwach rhythmisch",
    "ρ̄ < 0.40: sehr schwach rhythmisch"
  )), startRow = leg_row, colNames = FALSE)
setColWidths(wb, "Rayleigh", cols = 1:6,
  widths = c(10, 14, 22, 24, 15, 22))
freezePane(wb, "Rayleigh", firstRow = TRUE, firstActiveRow = 4)

# ────────────────────────────────────────────────────────────
# Sheet 4: Aktivitätsbeginn
# ────────────────────────────────────────────────────────────
addWorksheet(wb, "Aktivitaetsbeginn")

if (nrow(onset_summary) > 0) {
  onset_export <- onset_summary[, .(
    Igel                = igel,
    `Median (min nach SU)` = round(med_onset_min, 0),
    `Q25 (min)`         = round(q25_onset, 0),
    `Q75 (min)`         = round(q75_onset, 0),
    `N Nächte`          = n_naechte
  )]
  title_onset <- data.frame(V1 =
    "Aktivitätsbeginn relativ zu Sonnenuntergang (SU) — Median pro Igel",
    stringsAsFactors = FALSE)
  writeData(wb, "Aktivitaetsbeginn", title_onset, startRow = 1, colNames = FALSE)
  addStyle(wb, "Aktivitaetsbeginn",
    createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
  writeData(wb, "Aktivitaetsbeginn", as.data.frame(onset_export),
            startRow = 3, colNames = TRUE)
  addStyle(wb, "Aktivitaetsbeginn", stil_header, rows = 3, cols = 1:5, gridExpand = TRUE)
  for (i in seq_len(nrow(onset_export))) {
    if (i %% 2 == 0)
      addStyle(wb, "Aktivitaetsbeginn", stil_grau, rows = i + 3, cols = 1:5,
               gridExpand = TRUE, stack = TRUE)
  }
  setColWidths(wb, "Aktivitaetsbeginn", cols = 1:5, widths = c(10, 22, 12, 12, 12))
}

# ────────────────────────────────────────────────────────────
# Sheet 5: Metadaten-Verknüpfung
# ────────────────────────────────────────────────────────────
addWorksheet(wb, "Metadaten_Verknuepfung")

if (!is.null(chrono_meta) && nrow(chrono_meta) > 0) {
  meta_export_cols <- intersect(
    c("igel", "Sex", "Age_class", "Diagnosis_group", "Season",
      "Rel__weight__g_", "Rehab__d_", "Tracking__d_",
      "rho", "mu_h", "acrophase", "amplitude", "p_cosinor", "med_onset_min"),
    names(chrono_meta))
  meta_export <- as.data.frame(chrono_meta[, ..meta_export_cols])
  title_meta <- data.frame(V1 =
    "Metadaten × Chronobiologie — eine Zeile pro Igel",
    stringsAsFactors = FALSE)
  writeData(wb, "Metadaten_Verknuepfung", title_meta, startRow = 1, colNames = FALSE)
  addStyle(wb, "Metadaten_Verknuepfung",
    createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
  writeData(wb, "Metadaten_Verknuepfung", meta_export, startRow = 3, colNames = TRUE)
  addStyle(wb, "Metadaten_Verknuepfung", stil_header,
           rows = 3, cols = seq_along(meta_export_cols), gridExpand = TRUE)
  for (i in seq_len(nrow(meta_export))) {
    if (i %% 2 == 0)
      addStyle(wb, "Metadaten_Verknuepfung", stil_grau,
               rows = i + 3, cols = seq_along(meta_export_cols),
               gridExpand = TRUE, stack = TRUE)
  }
  setColWidths(wb, "Metadaten_Verknuepfung",
               cols = seq_along(meta_export_cols),
               widths = rep(14, length(meta_export_cols)))
}

xlsx_pfad <- file.path(output_ordner, "Block3_Chronobiologie_Uebersicht.xlsx")
saveWorkbook(wb, xlsx_pfad, overwrite = TRUE)
cat("✓ Excel-Übersicht gespeichert:", basename(xlsx_pfad), "\n\n")

cat("── Erzeugte Dateien ──\n")
cat("  Chronobiologie/actogramm_<Igel>.png         — je ein Actogramm\n")
cat("  Chronobiologie/actogramme_alle_uebersicht.png\n")
cat("  Chronobiologie/cosinor_profile_alle.png\n")
cat("  Chronobiologie/cosinor_akrophase_polar.png\n")
cat("  Chronobiologie/cosinor_akrophase_zeitverlauf.png\n")
cat("  Chronobiologie/cosinor_amplitude_zeitverlauf.png\n")
cat("  Chronobiologie/rayleigh_polar.png\n")
cat("  Chronobiologie/rayleigh_balken.png\n")
cat("  Chronobiologie/onset_boxplot.png\n")
cat("  Chronobiologie/onset_zeitverlauf.png\n")
cat("  Chronobiologie/meta_rho_diagnose.png\n")
cat("  Chronobiologie/meta_akrophase_saison.png\n")
cat("  Chronobiologie/meta_rho_alter_sex.png\n")
cat("  Chronobiologie/chrono_ergebnisse.rds\n")
cat("  Block3_Chronobiologie_Bericht.docx\n")
cat("  Block3_Chronobiologie_Uebersicht.xlsx   ← wird bei jedem Lauf neu erzeugt\n")
cat("\n✓ Chronobiologische Analyse abgeschlossen!\n")

sink(file.path(output_ordner, "chrono_kennzahlen.txt"))
print(cosinor_gesamt_dt)
print(rayleigh_dt)
if (!is.null(onset_summary)) print(onset_summary[, .(igel, med_onset_min, n_naechte)])
sink()

# ══════════════════════════════════════════════════════════════
# 8. PUBLIKATIONS-FIGUREN (englisch, druckfertig)
# ══════════════════════════════════════════════════════════════
cat("════════════════════════════════════\n")
cat("8. Publikations-Figuren\n")
cat("════════════════════════════════════\n\n")

pub_ordner <- file.path(output_ordner, "pub_figures")
dir.create(pub_ordner, showWarnings = FALSE)

# ── Robuste PDF-Ausgabe ───────────────────────────────────────────────────────
# Manche Windows-R-Installationen haben keine funktionierende eingebaute Cairo-DLL
# ("failed to load cairo DLL"). Das Paket 'Cairo' bringt sein eigenes Cairo mit
# und rendert Unicode-Sonderzeichen (z. B. ρ̄) korrekt. Fallback: Basis-pdf().
open_pdf <- function(path, width, height) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoPDF(file = path, width = width, height = height)
  } else {
    message("Hinweis: Paket 'Cairo' nicht installiert – nutze Basis-pdf(). ",
            "Sonderzeichen wie ρ̄ werden dort evtl. nicht korrekt dargestellt. ",
            "Fuer beste PDF-Qualitaet: install.packages(\"Cairo\").")
    pdf(file = path, width = width, height = height)
  }
}

# ── Fig. 2A: Übersichts-Actogramm (Supplement-Qualität, englisch) ─────────────
# Verbesserte Version des Übersichtsplots: englisch, größere Schrift,
# lesbare Achsen. Bleibt Supplement, da 25+ Tiere immer gedrängt wirken.
# ─────────────────────────────────────────────────────────────────────────────

# Tierkürzel: Igel1 → H1
dt_pub <- copy(dt)
dt_pub[, hedgehog := paste0("H", gsub("Igel", "", igel))]
dt_pub[, hedgehog := factor(hedgehog,
  levels = paste0("H", gsub("Igel", "", levels(dt$igel))))]

# Nur Tiere in igel_liste (nach >=60-min-Filter)
igel_liste_pub <- paste0("H", gsub("Igel", "", igel_liste))
dt_pub_filt <- dt_pub[hedgehog %in% igel_liste_pub]

actogramm_pub2 <- dt_pub_filt[, .(aktiv_anteil = mean(aktiv)),
                                by = .(hedgehog, datum, tage_seit,
                                       bin30 = floor(stunde * 2) / 2)]
actogramm_pub2 <- merge(actogramm_pub2,
                         sonnen_dt[, .(datum, aufgang_h, untergang_h)],
                         by = "datum", all.x = TRUE)

p_acto_supp <- ggplot(actogramm_pub2[tage_seit <= 21],
                       aes(x = bin30, y = tage_seit, fill = aktiv_anteil)) +
  geom_tile(width = 0.5, height = 0.9) +
  # Sonnenuntergang
  geom_segment(
    data = actogramm_pub2[tage_seit <= 21,
                           .(untergang_h = mean(untergang_h, na.rm = TRUE)),
                           by = .(hedgehog, tage_seit)],
    aes(x = untergang_h, xend = untergang_h,
        y = tage_seit - 0.5, yend = tage_seit + 0.5),
    colour = "#F4A261", linewidth = 0.5, inherit.aes = FALSE
  ) +
  # Sonnenaufgang
  geom_segment(
    data = actogramm_pub2[tage_seit <= 21,
                           .(aufgang_h = mean(aufgang_h, na.rm = TRUE)),
                           by = .(hedgehog, tage_seit)],
    aes(x = aufgang_h, xend = aufgang_h,
        y = tage_seit - 0.5, yend = tage_seit + 0.5),
    colour = "#90E0EF", linewidth = 0.5, inherit.aes = FALSE
  ) +
  facet_wrap(~ hedgehog, ncol = 5) +
  scale_fill_gradient(low = "white", high = "#1a2744",
                      name = "Activity", limits = c(0, 1),
                      breaks = c(0, 0.5, 1),
                      labels = c("0%", "50%", "100%")) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                     labels = c("00:00", "06:00", "12:00", "18:00", "00:00")) +
  scale_y_reverse(breaks = c(1, 7, 14, 21)) +
  labs(
    x = "Time of day",
    y = "Days post-release"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid    = element_blank(),
    strip.text    = element_text(face = "bold", size = 8),
    axis.text.x   = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y   = element_text(size = 7),
    legend.position = "right",
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 7)
  )

open_pdf(file.path(pub_ordner, "FigS1_actogram_overview.pdf"),
    width = 16.25, height = 12.5)
print(p_acto_supp)
dev.off()
cat("✓ Supplementary actogram overview gespeichert\n")

# ── Fig. 2 Panel: 4 repräsentative Actogramme (Haupttext) ────────────────────
# Auswahl:
#   H14 — stark rhythmisch, sichtig      (ρ̄ = 0.809)
#   H9  — moderat rhythmisch, sichtig    (ρ̄ = 0.712)
#   H22 — blind, tagesaktiv              (ρ̄ = 0.309)
#   H23 — blind, verschobene Phase       (ρ̄ = 0.417)
# ─────────────────────────────────────────────────────────────────────────────

# Hilfsfunktion: Englisches Double-Plot-Actogramm
plot_actogramm_pub <- function(igel_name,
                               label_title,
                               label_rho,
                               highlight_blind = FALSE,
                               max_tage        = NULL,
                               legend_pos      = "none") {

  ad <- actogramm_data[igel == igel_name]
  if (nrow(ad) == 0) return(NULL)

  # Optionale Begrenzung auf ersten max_tage Tage
  if (!is.null(max_tage)) ad <- ad[tage_seit <= max_tage]
  if (nrow(ad) == 0) return(NULL)

  t_range <- ad[, range(tage_seit)]
  t_seq   <- t_range[1]:t_range[2]

  dp_links  <- ad[, .(tage_seit, bin10, aktiv_anteil, aufgang_h, untergang_h)]
  dp_rechts <- copy(dp_links)
  dp_rechts[, bin10      := bin10      + 24]
  dp_rechts[, aufgang_h   := aufgang_h  + 24]
  dp_rechts[, untergang_h := untergang_h + 24]
  dp_rechts[, tage_seit  := tage_seit - 1]

  dp_all <- rbind(dp_links, dp_rechts)
  dp_all <- dp_all[tage_seit %in% t_seq]

  sonnen_igel <- merge(
    data.table(tage_seit = t_seq),
    ad[, .(aufgang_h    = mean(aufgang_h,   na.rm = TRUE),
           untergang_h  = mean(untergang_h, na.rm = TRUE)),
       by = tage_seit],
    by = "tage_seit", all.x = TRUE
  )
  sonnen_rechts <- copy(sonnen_igel)
  sonnen_rechts[, aufgang_h   := aufgang_h   + 24]
  sonnen_rechts[, untergang_h := untergang_h + 24]
  sonnen_rechts[, tage_seit   := tage_seit - 1]
  sonnen_dp <- rbind(sonnen_igel, sonnen_rechts)
  sonnen_dp <- sonnen_dp[tage_seit %in% t_seq]

  # Titelfarbe: rot für blinde Tiere
  title_col <- if (highlight_blind) "#C0392B" else "black"

  ggplot(dp_all, aes(x = bin10, y = tage_seit, fill = aktiv_anteil)) +
    # Mitternachts-Trennlinie zwischen linker (Tag d) und rechter (Tag d+1) Hälfte
    geom_vline(xintercept = 24, colour = "grey80",
               linewidth = 0.4, linetype = "22") +
    geom_tile(width = 1/6, height = 0.95) +
    geom_segment(data = sonnen_dp,
                 aes(x = untergang_h, xend = untergang_h,
                     y = tage_seit - 0.5, yend = tage_seit + 0.5),
                 colour = "#F4A261", linewidth = 0.9, inherit.aes = FALSE) +
    geom_segment(data = sonnen_dp,
                 aes(x = aufgang_h, xend = aufgang_h,
                     y = tage_seit - 0.5, yend = tage_seit + 0.5),
                 colour = "#90E0EF", linewidth = 0.9, inherit.aes = FALSE) +
    scale_fill_gradient(low = "white", high = "#1a2744",
                        name   = "Activity",
                        limits = c(0, 1),
                        breaks = c(0, 0.5, 1),
                        labels = c("0 %", "50 %", "100 %")) +
    scale_x_continuous(
      limits = c(0, 48),
      breaks = c(0, 12, 24, 36, 48),
      labels = c("00:00", "12:00", "00:00", "12:00", "00:00"),
      expand = c(0, 0)
    ) +
    # Ganzzahlige Tag-Achse alle 2 Tage (statt pretty_breaks → keine Halb-Tage)
    scale_y_reverse(breaks = seq(0, 60, 2), expand = c(0, 0.5)) +
    labs(
      title    = label_title,
      subtitle = bquote(bar(rho) == .(label_rho)),
      x        = "Time of day",
      y        = "Days post-release"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", size = 15,
                                      colour = title_col),
      plot.subtitle    = element_text(size = 12, colour = "grey30"),
      panel.grid       = element_blank(),
      axis.text.x      = element_text(size = 12),
      axis.text.y      = element_text(size = 12),
      axis.title       = element_text(size = 13),
      legend.position   = legend_pos,
      legend.title      = element_text(size = 11, face = "bold"),
      legend.text       = element_text(size = 10),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width  = unit(0.5, "cm"),
      legend.margin     = margin(4, 8, 4, 8),
      plot.background   = element_rect(colour = "black", fill = NA,
                                       linewidth = 0.6),
      plot.margin       = margin(30, 24, 8, 24)
    )
}

# Gemeinsames Zeitfenster für alle 4 Panels (erste N Tage) — eine Zahl ändern reicht
TAGE_FENSTER_FIG2 <- 10

# Die 4 ausgewählten Tiere (alle auf dasselbe Fenster begrenzt → gleiche Höhe, besser lesbar)
p_h14 <- plot_actogramm_pub("Igel14", "H14 — sighted, strongly rhythmic",
                             label_rho = 0.809,
                             max_tage  = TAGE_FENSTER_FIG2)
p_h9  <- plot_actogramm_pub("Igel9",  "H9 — sighted, moderately rhythmic",
                             label_rho = 0.712,
                             max_tage  = TAGE_FENSTER_FIG2)
p_h22 <- plot_actogramm_pub("Igel22", "H22 — visually impaired, near-diurnal",
                             label_rho = 0.309, highlight_blind = TRUE,
                             max_tage  = TAGE_FENSTER_FIG2)
# H23 bekommt die Legende — wird von patchwork gesammelt und unten platziert
p_h23 <- plot_actogramm_pub("Igel23", "H23 — visually impaired, phase-shifted",
                             label_rho  = 0.417, highlight_blind = TRUE,
                             max_tage   = TAGE_FENSTER_FIG2,
                             legend_pos = "bottom")

# Kombinieren mit patchwork (2×2 Grid, gemeinsame Legende)
p_acto_4panel <- (p_h14 + p_h9 + p_h22 + p_h23) +
  patchwork::plot_layout(ncol = 2, guides = "collect") +
  patchwork::plot_annotation(
    tag_levels = "a",
    tag_suffix = ")",
    caption = paste0(
      "Dark fill = active | Orange line = sunset | Blue line = sunrise\n",
      "Double-plot: each row spans 48 h (day d and d+1 side by side)"),
    theme = theme(
      plot.caption = element_text(size = 9, colour = "grey40",
                                  hjust = 0, lineheight = 1.3)
    )
  ) &
  theme(
    plot.tag          = element_text(face = "bold", size = 15),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = 11, face = "bold",
                                     vjust = 0.8),
    legend.text       = element_text(size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(3.0, "cm")
  )

open_pdf(file.path(pub_ordner, "Fig2A_actogram_4panel.pdf"),
    width  = 16.36,
    height = 9.77)   # 10-Tage-Fenster → flachere Panels, ausgewogeneres Seitenverhältnis
print(p_acto_4panel)
dev.off()
cat("✓ Fig 2A — 4-Panel-Actogramm gespeichert\n")

# ── Fig. 2B: Rayleigh-Polarplot (Publikationsversion) ────────────────────────
open_pdf(file.path(pub_ordner, "Fig2B_rayleigh_polar.pdf"),
    width = 8, height = 8)
print(p_rayleigh_polar)
dev.off()
cat("✓ Fig 2B — Rayleigh-Polarplot gespeichert\n")

cat("\n✓ Alle Publikations-Figuren gespeichert in:", pub_ordner, "\n\n")