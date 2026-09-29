# ==============================================================
# Block 2d — release-effect analysis
# ==============================================================
# Project:  Hedgehog VHF telemetry, Wildtierstation Sachsenhagen
# Author:   Natalie Steiner
#
# Question: Does behaviour in the first 2-3 nights after release
#           differ significantly from later nights?
#
# Data (from Block0_Datenpipeline.R):
#   gamm_tagesdaten.rds   — per hedgehog x night: % active minutes (pct_aktiv 0-1)
#   chrono_minuten24h.rds — per hedgehog x minute (24h): active 0/1, hour 0-24
#   Classification: pred_nested_loio_smoothed_wmv (all hedgehogs incl. 1-6)
#
# Analyses:
#   1. Night activity phase 1 (day 1-3) vs. phase 2 (day 4+)
#      -> Wilcoxon signed-rank test (paired), spaghetti plot, pair plot
#   2. Activity onset relative to sunset (first 2 weeks)
#      -> first active time per night after sunset
#   3. Between-individual variability over time
#      -> SD of night activity across hedgehogs per day
#   4. Trajectory-pattern classification

# ──────────────────────────────────────────────────────────────
# 0. PAKETE
# ──────────────────────────────────────────────────────────────

pakete <- c("data.table", "ggplot2", "patchwork", "scales",
            "suncalc", "officer", "flextable", "openxlsx")
neu    <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) install.packages(neu)

invisible(lapply(pakete, library, character.only = TRUE))

# ggrepel: optional (nur für Plot 4 — Labels ohne Überlappung)
use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (use_ggrepel) library(ggrepel)

cat("✓ Pakete geladen\n\n")

# ──────────────────────────────────────────────────────────────
# 1. EINSTELLUNGEN
# ──────────────────────────────────────────────────────────────

projekt_root  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
output_ordner <- file.path(projekt_root, "output", "Block2d_ReleaseEffect")
dir.create(output_ordner, showWarnings = FALSE, recursive = TRUE)

# Koordinaten Wildtierstation Sachsenhagen
lat_h <- 52.397
lon_h <-  9.217

# Phasendefinition
phase1_bis <- 3L      # Tage 1–3  = Release-Phase
n_min_igel <- 5L      # Populations-Cutoff: min. N Igel pro Tag

# Ausschlusskriterien (identisch mit Block2_GAMM.R)
min_naechte   <- 3
min_nacht_min <- 60
min_pct_aktiv <- 1    # %

# ──────────────────────────────────────────────────────────────
# 2. DATEN LADEN
# ──────────────────────────────────────────────────────────────

cat("Lade Datensätze...\n")

tage_rds <- file.path(projekt_root, "output", "Block0_Pipeline", "gamm_tagesdaten.rds")
min_rds  <- file.path(projekt_root, "output", "Block0_Pipeline", "chrono_minuten24h.rds")

if (!file.exists(tage_rds)) stop("gamm_tagesdaten.rds fehlt — Block0 ausführen.")
if (!file.exists(min_rds))  stop("chrono_minuten24h.rds fehlt — Block0 ausführen.")

dt_tage <- readRDS(tage_rds); setDT(dt_tage)
dt_chrono <- readRDS(min_rds); setDT(dt_chrono)

# igel immer als character (robust gegen factor / character)
dt_tage[,   igel := as.character(igel)]
dt_chrono[, igel := as.character(igel)]

cat("Tagesdaten:    ", nrow(dt_tage),   "Zeilen\n")
cat("Minutendaten:  ", format(nrow(dt_chrono), big.mark = "'"), "Zeilen\n\n")

# ──────────────────────────────────────────────────────────────
# 3. QC — AUSSCHLUSSKRITERIEN (wie Block2_GAMM.R)
# ──────────────────────────────────────────────────────────────

igel_qc <- dt_tage[tageszeit == "Nacht", .(
  n_naechte    = uniqueN(datum),
  n_nacht_min  = sum(n_min, na.rm = TRUE),
  mean_aktiv_n = mean(pct_aktiv * 100, na.rm = TRUE)
), by = igel]

igel_qc[, ausschluss := fcase(
  n_naechte    < min_naechte,   paste0("< ", min_naechte, " Nächte"),
  n_nacht_min  < min_nacht_min, paste0("< ", min_nacht_min, " Nachtminuten"),
  mean_aktiv_n < min_pct_aktiv, "< 1% Nachtaktivität",
  default = NA_character_
)]

igel_ok <- igel_qc[is.na(ausschluss), as.character(igel)]

if (length(igel_ok) == 0) stop("Kein Igel erfüllt die QC-Kriterien.")
cat("QC: ", length(igel_ok), "Igel eingeschlossen\n")

dt_tage   <- dt_tage[igel %in% igel_ok]
dt_chrono <- dt_chrono[igel %in% igel_ok]

# ──────────────────────────────────────────────────────────────
# 4. POPULATIONS-CUTOFF (N ≥ 5 Igel pro Tag, wie Block2_GAMM.R)
# ──────────────────────────────────────────────────────────────

n_pro_tag  <- dt_tage[tageszeit == "Nacht" & !is.na(pct_aktiv),
                       .(n_igel = uniqueN(igel)), by = tage_seit]
if (n_pro_tag[, max(n_igel, na.rm = TRUE)] < n_min_igel)
  stop("N-at-risk-Schwelle nie erreicht — n_min_igel verringern.")

cutoff_tag <- n_pro_tag[n_igel >= n_min_igel, max(tage_seit)]
cat("Populations-Cutoff: Tag", cutoff_tag, "(N ≥", n_min_igel, "Igel)\n\n")

# ──────────────────────────────────────────────────────────────
# 5. NACHT-TAGESDATEN (Basis für Analysen 1, 3, 4)
# ──────────────────────────────────────────────────────────────

dt_nacht <- dt_tage[
  tageszeit == "Nacht" &
  !is.na(pct_aktiv) &
  tage_seit >= 1L &
  tage_seit <= cutoff_tag
]
dt_nacht[, pct_pct := pct_aktiv * 100]   # Anteil → Prozent

# Phasenzuordnung
dt_nacht[, phase := ifelse(tage_seit <= phase1_bis, "Release (Tag 1-3)", "Etablierung (Tag 4+)")]
dt_nacht[, phase := factor(phase, levels = c("Release (Tag 1-3)", "Etablierung (Tag 4+)"))]

# Igel mit Daten in BEIDEN Phasen
phase_check <- dt_nacht[, .(hat_p1 = any(tage_seit <= phase1_bis),
                             hat_p2 = any(tage_seit >  phase1_bis)), by = igel]
igel_beide  <- phase_check[hat_p1 == TRUE & hat_p2 == TRUE, igel]
cat("Igel mit Daten in beiden Phasen:", length(igel_beide), "\n\n")

# ──────────────────────────────────────────────────────────────
# 6. SONNENUNTERGANG (für Analyse 2)
# ──────────────────────────────────────────────────────────────

alle_daten <- sort(unique(dt_chrono$datum))
sonnen_dt  <- as.data.table(getSunlightTimes(
  date = alle_daten, lat = lat_h, lon = lon_h,
  tz = "Europe/Berlin", keep = c("sunrise", "sunset")))
sonnen_dt[, datum       := as.Date(date)]
sonnen_dt[, aufgang_h   := as.numeric(format(sunrise, "%H")) +
                           as.numeric(format(sunrise, "%M")) / 60]
sonnen_dt[, untergang_h := as.numeric(format(sunset,  "%H")) +
                           as.numeric(format(sunset,  "%M")) / 60]
cat("✓ Sonnenzeiten geladen\n\n")

# ══════════════════════════════════════════════════════════════
# ANALYSE 1: NACHTAKTIVITÄT — RELEASE-PHASE vs. ETABLIERUNG
# ══════════════════════════════════════════════════════════════

cat("════════ Analyse 1: Release vs. Etablierung ════════\n\n")

# Mittlere Nachtaktivität je Igel × Phase
a1_phase <- dt_nacht[igel %in% igel_beide, .(
  mittel = mean(pct_pct, na.rm = TRUE),
  sd     = sd(pct_pct,   na.rm = TRUE),
  n      = .N
), by = .(igel, phase)]

# Auf breites Format für paarweisen Test
a1_wide <- dcast(a1_phase, igel ~ phase, value.var = "mittel")
# Spaltennamen sicher umbenennen (unabhängig von der Schreibweise)
sp_names <- names(a1_wide)
sp_p1    <- sp_names[grepl("Release",    sp_names)]
sp_p2    <- sp_names[grepl("Etablierung", sp_names)]

if (length(sp_p1) == 1 && length(sp_p2) == 1) {
  setnames(a1_wide, c(sp_p1, sp_p2), c("p1", "p2"))
} else {
  stop("dcast-Spalten nicht erkannt — Phase-Labels prüfen.")
}
a1_wide <- a1_wide[!is.na(p1) & !is.na(p2)]

# Wilcoxon Signed-Rank Test
wt       <- wilcox.test(a1_wide$p1, a1_wide$p2, paired = TRUE, exact = FALSE)
med_diff <- median(a1_wide$p1 - a1_wide$p2, na.rm = TRUE)
cat("Wilcoxon (N =", nrow(a1_wide), "Igel):\n")
cat("  Median Differenz (Release – Etablierung):", round(med_diff, 1), "PP\n")
cat("  p =", format.pval(wt$p.value, digits = 3), "\n\n")

# Richtungs-Flag für Paarplot
a1_wide[, richtung := ifelse(p2 > p1, "gestiegen", "gesunken")]

# ── Plot 1: Spaghetti-Plot ────────────────────────────────────
n_ar <- dt_nacht[, .(n = uniqueN(igel)), by = tage_seit]
pop  <- dt_nacht[, .(m = mean(pct_pct, na.rm=TRUE),
                      se = sd(pct_pct, na.rm=TRUE) / sqrt(.N)), by = tage_seit]

p1_spaghetti <- ggplot() +
  annotate("rect", xmin = 0.5, xmax = phase1_bis + 0.5,
           ymin = -Inf, ymax = Inf, fill = "#F4A261", alpha = 0.13) +
  geom_vline(xintercept = phase1_bis + 0.5,
             linetype = "dashed", color = "#C0392B", linewidth = 0.7) +
  annotate("text", x = (1 + phase1_bis) / 2, y = 98,
           label = "Release\n(Tag 1–3)", size = 2.9, color = "#C0392B",
           fontface = "italic", hjust = 0.5, vjust = 1) +
  geom_line(data = dt_nacht,
            aes(x = tage_seit, y = pct_pct, group = igel, color = igel),
            alpha = 0.35, linewidth = 0.65) +
  geom_ribbon(data = pop, aes(x = tage_seit, ymin = m - se, ymax = m + se),
              fill = "black", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = pop, aes(x = tage_seit, y = m),
            color = "black", linewidth = 1.4, inherit.aes = FALSE) +
  geom_text(data = n_ar, aes(x = tage_seit, y = -6, label = n),
            size = 2.3, color = "grey55") +
  annotate("text", x = 0.3, y = -6, label = "N=",
           size = 2.3, color = "grey55", hjust = 1) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  scale_x_continuous(breaks = sort(unique(c(1, 3, 7, 14, 21, cutoff_tag)))) +
  scale_y_continuous(limits = c(-9, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(title    = "Nachtaktivität aller Igel — Tage 1 bis 27 nach Auswilderung",
       subtitle = paste0("Orange = Release-Phase (Tag 1–", phase1_bis, ") | ",
                         "Jede Linie = ein Igel | Schwarz = Pop.-Mittel ± SE\n",
                         "N = Anzahl Igel mit Daten je Tag"),
       x = "Tage nach Auswilderung", y = "Nachtaktivität (%)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ── Plot 2: Paarvergleich ─────────────────────────────────────
med_p1 <- median(a1_wide$p1, na.rm = TRUE)
med_p2 <- median(a1_wide$p2, na.rm = TRUE)
y_max  <- max(c(a1_wide$p1, a1_wide$p2), na.rm = TRUE)

p2_paar <- ggplot() +
  geom_segment(data = a1_wide,
               aes(x = 1, xend = 2, y = p1, yend = p2, color = richtung),
               linewidth = 0.9, alpha = 0.65,
               arrow = arrow(length = unit(0.07, "inches"), ends = "last")) +
  geom_point(data = a1_wide, aes(x = 1, y = p1), size = 2.8,
             color = "#F4A261", alpha = 0.85) +
  geom_point(data = a1_wide, aes(x = 2, y = p2), size = 2.8,
             color = "#2C3E6B", alpha = 0.85) +
  # Mediane
  annotate("point", x = 1, y = med_p1, size = 5, shape = 18, color = "black") +
  annotate("point", x = 2, y = med_p2, size = 5, shape = 18, color = "black") +
  annotate("segment", x = 1, xend = 2, y = med_p1, yend = med_p2,
           linewidth = 1.6, color = "black", linetype = "dashed") +
  annotate("text", x = 1.5, y = y_max * 1.06,
           label = paste0("Wilcoxon p ",
                          ifelse(wt$p.value < 0.001, "< 0.001",
                                 paste0("= ", round(wt$p.value, 3)))),
           size = 3.4, fontface = "italic") +
  scale_color_manual(values = c("gestiegen" = "#2C7BB6", "gesunken" = "#D7191C"),
                     name = "Aktivität ab Tag 4:") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c(paste0("Tag 1–", phase1_bis, "\n(Release)"),
                                paste0("Tag ", phase1_bis + 1, "+\n(Etablierung)")),
                     limits = c(0.65, 2.35)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title    = paste0("Release-Effekt: Nachtaktivität — Phase 1 vs. Phase 2\n",
                          "N = ", nrow(a1_wide), " Igel mit Daten in beiden Phasen"),
       subtitle = "Jede Linie = ein Igel | Raute = Median",
       x = NULL, y = "Mittlere Nachtaktivität (%)") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "bottom")

# Speichern
png(file.path(output_ordner, "p1_spaghetti.png"),     width=1600, height=900, res=130)
print(p1_spaghetti); dev.off()
png(file.path(output_ordner, "p2_paarvergleich.png"), width=1000, height=950, res=130)
print(p2_paar); dev.off()
cat("✓ Plots 1 & 2 gespeichert\n\n")

# ══════════════════════════════════════════════════════════════
# ANALYSE 2: AKTIVITÄTSBEGINN RELATIV ZU SONNENUNTERGANG
# ══════════════════════════════════════════════════════════════

cat("════════ Analyse 2: Aktivitätsbeginn ════════\n\n")

# Merge: Minuten × Sonnenzeiten
dt_cs <- merge(dt_chrono[tage_seit >= 1 & tage_seit <= min(14L, cutoff_tag)],
               sonnen_dt[, .(datum, aufgang_h, untergang_h)],
               by = "datum", all.x = TRUE)

# Erster aktiver Zeitpunkt pro Igel × Nacht (abendliches Fenster: SU-1h bis Mitternacht)
onset_dt <- dt_cs[
  aktiv == 1 &
  !is.na(untergang_h) &
  stunde >= (untergang_h - 1) &
  stunde <  24,
  .(onset_h = min(stunde, na.rm = TRUE),
    untergang_h = first(untergang_h)),
  by = .(igel, datum, tage_seit)
]
onset_dt[, onset_min := (onset_h - untergang_h) * 60]  # Minuten nach Sonnenuntergang
onset_dt <- onset_dt[onset_min >= -60]  # entferne extreme Ausreißer

if (nrow(onset_dt) > 0) {
  cat("Aktivitätsbeginn (Median Minuten nach Sonnenuntergang):\n")
  cat("  Release-Phase  :", round(median(onset_dt[tage_seit <= phase1_bis, onset_min], na.rm=TRUE), 0), "min\n")
  cat("  Etablierungsph.:", round(median(onset_dt[tage_seit >  phase1_bis, onset_min], na.rm=TRUE), 0), "min\n\n")

  # Populations-Zusammenfassung pro Tag
  pop_onset <- onset_dt[, .(med_min = median(onset_min, na.rm=TRUE),
                              q25 = quantile(onset_min, 0.25, na.rm=TRUE),
                              q75 = quantile(onset_min, 0.75, na.rm=TRUE),
                              n   = .N), by = tage_seit]

  p3_onset <- ggplot() +
    annotate("rect", xmin = 0.5, xmax = phase1_bis + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#F4A261", alpha = 0.13) +
    geom_vline(xintercept = phase1_bis + 0.5,
               linetype = "dashed", color = "#C0392B", linewidth = 0.7) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.8) +
    annotate("text", x = min(14, cutoff_tag) + 0.2, y = 2,
             label = "↑ nach SU", size = 2.6, color = "grey40",
             hjust = 0, fontface = "italic") +
    geom_jitter(data = onset_dt, aes(x = tage_seit, y = onset_min, color = igel),
                width = 0.18, alpha = 0.45, size = 1.9) +
    geom_ribbon(data = pop_onset, aes(x = tage_seit, ymin = q25, ymax = q75),
                fill = "black", alpha = 0.12, inherit.aes = FALSE) +
    geom_line(data = pop_onset, aes(x = tage_seit, y = med_min),
              color = "black", linewidth = 1.4, inherit.aes = FALSE) +
    geom_text(data = pop_onset,
              aes(x = tage_seit,
                  y = min(onset_dt$onset_min, na.rm = TRUE) - 12,
                  label = n),
              size = 2.3, color = "grey55") +
    scale_color_viridis_d(option = "turbo", guide = "none") +
    scale_x_continuous(breaks = sort(unique(c(1, 3, 7, min(14, cutoff_tag))))) +
    scale_y_continuous(labels = function(x) paste0(ifelse(x >= 0, "+", ""), x, " min")) +
    labs(title    = "Aktivitätsbeginn: Wann werden die Igel abends aktiv?",
         subtitle = paste0("Minuten nach Sonnenuntergang (0 = exakt Sonnenuntergang)\n",
                           "Orange = Release-Phase | Schwarz = Median ± IQR | N = Igel-Nächte"),
         x = "Tage nach Auswilderung",
         y = "Aktivitätsbeginn (Minuten nach Sonnenuntergang)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  png(file.path(output_ordner, "p3_aktivitaetsbeginn.png"), width=1500, height=900, res=130)
  print(p3_onset); dev.off()
  cat("✓ Plot 3 gespeichert\n\n")
} else {
  cat("⚠ Keine Onset-Daten verfügbar — Plot 3 wird übersprungen\n\n")
  pop_onset <- data.table()
}

# ══════════════════════════════════════════════════════════════
# ANALYSE 3: VARIABILITÄT ZWISCHEN INDIVIDUEN ÜBER ZEIT
# ══════════════════════════════════════════════════════════════

cat("════════ Analyse 3: Variabilität über Zeit ════════\n\n")

variab <- dt_nacht[, .(sd  = sd(pct_pct,  na.rm = TRUE),
                        m   = mean(pct_pct, na.rm = TRUE),
                        n   = .N),
                   by = tage_seit]

cat("SD Nachtaktivität:\n")
cat("  Release-Phase  :", round(variab[tage_seit <= phase1_bis, mean(sd, na.rm=TRUE)], 1), "PP\n")
cat("  Etablierungsph.:", round(variab[tage_seit >  phase1_bis, mean(sd, na.rm=TRUE)], 1), "PP\n\n")

p4_variab <- ggplot(variab, aes(x = tage_seit)) +
  annotate("rect", xmin = 0.5, xmax = phase1_bis + 0.5,
           ymin = -Inf, ymax = Inf, fill = "#F4A261", alpha = 0.13) +
  geom_vline(xintercept = phase1_bis + 0.5,
             linetype = "dashed", color = "#C0392B", linewidth = 0.7) +
  geom_col(aes(y = sd, fill = tage_seit <= phase1_bis), alpha = 0.75, width = 0.85) +
  geom_smooth(aes(y = sd), method = "loess", span = 0.55,
              color = "black", linewidth = 1.2, se = TRUE,
              fill = "grey30", alpha = 0.15, inherit.aes = TRUE) +
  geom_text(aes(y = sd + 0.9, label = n), size = 2.3, color = "grey50") +
  scale_fill_manual(values = c("TRUE"  = "#F4A261", "FALSE" = "#2C3E6B"),
                    labels = c("TRUE"  = paste0("Tag 1–", phase1_bis),
                               "FALSE" = paste0("Tag ", phase1_bis + 1, "+")),
                    name = "Phase") +
  scale_x_continuous(breaks = sort(unique(c(1, 3, 7, 14, 21, cutoff_tag)))) +
  labs(title    = "Variabilität der Nachtaktivität zwischen Individuen",
       subtitle = paste0("SD aller Igel pro Tag | Hohe SD = große Unterschiede zwischen Tieren\n",
                         "N über Balken = Anzahl Igel"),
       x = "Tage nach Auswilderung",
       y = "SD Nachtaktivität (Prozentpunkte)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

png(file.path(output_ordner, "p4_variabilitaet.png"), width=1500, height=800, res=130)
print(p4_variab); dev.off()
cat("✓ Plot 4 gespeichert\n\n")

# ══════════════════════════════════════════════════════════════
# ANALYSE 4: VERLAUFSMUSTER-KLASSIFIKATION
# ══════════════════════════════════════════════════════════════

cat("════════ Analyse 4: Verlaufsmuster ════════\n\n")

norm_dt <- dt_nacht[igel %in% igel_beide, .(
  p1 = mean(pct_pct[tage_seit <= phase1_bis], na.rm = TRUE),
  p2 = mean(pct_pct[tage_seit >  phase1_bis], na.rm = TRUE),
  n1 = sum(tage_seit <= phase1_bis),
  n2 = sum(tage_seit >  phase1_bis)
), by = igel]
norm_dt <- norm_dt[!is.na(p1) & !is.na(p2) & n1 >= 1 & n2 >= 2]
norm_dt[, delta := p1 - p2]

# Klassifikation (positiv = Release höher, negativ = Etablierung höher)
norm_dt[, typ := fcase(
  delta < -15,       "Zunehmend",     # in Etablierung aktiver
  delta >  15,       "Abnehmend",     # in Release aktiver
  default = "Stabil"
)]
norm_dt[, typ := factor(typ, levels = c("Zunehmend", "Stabil", "Abnehmend"))]

cat("Klassifikation (N =", nrow(norm_dt), "Igel):\n")
print(norm_dt[order(delta), .(Igel = igel,
                               `Phase1 (%)` = round(p1, 1),
                               `Phase2 (%)` = round(p2, 1),
                               `Δ (PP)` = round(delta, 1),
                               Typ = as.character(typ))])
cat("\nTyp-Verteilung:\n")
print(norm_dt[, .N, by = typ][order(typ)])
cat("\n")

p5_norm <- ggplot(norm_dt, aes(x = p1, y = p2, color = typ, label = igel)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 1) +
  geom_abline(slope = 1, intercept = c(-15, 15),
              linetype = "dotted", color = "grey75", linewidth = 0.7) +
  geom_point(size = 4, alpha = 0.88) +
  {if (use_ggrepel) geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE)
   else             geom_text(vjust = -0.9, size = 2.9, show.legend = FALSE)} +
  scale_color_manual(
    values = c("Zunehmend" = "#2C7BB6", "Stabil" = "#1A7A4A", "Abnehmend" = "#D7191C"),
    labels = c("Zunehmend" = paste0("Zunehmend (Phase 2 > Phase 1, > 15 PP)"),
               "Stabil"    = "Stabil (Differenz ≤ 15 PP)",
               "Abnehmend" = paste0("Abnehmend (Phase 1 > Phase 2, > 15 PP)")),
    name = "Verlaufsmuster") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  labs(title    = "Verlaufsmuster je Individuum: Release vs. Etablierungsphase",
       subtitle = paste0("x = Phase 1 (Tag 1–", phase1_bis, ") | y = Phase 2 (Tag ",
                         phase1_bis + 1, "+)\n",
                         "Gestrichelt = kein Unterschied | Gepunktet = ±15 PP"),
       x = paste0("Nachtaktivität Tag 1–", phase1_bis, " (%)"),
       y = paste0("Nachtaktivität Tag ",   phase1_bis + 1, "+ (%)")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
        legend.text = element_text(size = 9))

png(file.path(output_ordner, "p5_verlaufsmuster.png"), width=1200, height=1050, res=130)
print(p5_norm); dev.off()
cat("✓ Plot 5 gespeichert\n\n")

# ══════════════════════════════════════════════════════════════
# EXCEL-ÜBERSICHT
# ══════════════════════════════════════════════════════════════

cat("Erstelle Excel-Übersicht...\n")

wb <- createWorkbook()

# Stil-Definitionen
h_style  <- createStyle(fontColour = "#FFFFFF", fgFill = "#2C3E6B", textDecoration = "bold",
                         halign = "center", border = "Bottom", borderColour = "#FFFFFF")
z_style  <- createStyle(fgFill = "#F0F4FA")
ok_style <- createStyle(fgFill = "#D5E8D4", halign = "center")
warn_style <- createStyle(fgFill = "#FFE6CC", halign = "center")

# ── Sheet 1: Übersicht pro Igel ──────────────────────────────
addWorksheet(wb, "Verlaufsmuster")
df1 <- as.data.frame(norm_dt[order(delta), .(
  Igel            = igel,
  `Aktivität Phase 1 (%)` = round(p1, 1),
  `Aktivität Phase 2 (%)` = round(p2, 1),
  `Differenz (PP)`        = round(delta, 1),
  `Nächte Phase 1`        = n1,
  `Nächte Phase 2`        = n2,
  `Verlaufsmuster`        = as.character(typ)
)])
writeData(wb, "Verlaufsmuster", df1, headerStyle = h_style)
# Zebrastreifen
for (r in seq(2, nrow(df1) + 1, 2)) addStyle(wb, "Verlaufsmuster", z_style, rows = r, cols = 1:ncol(df1), gridExpand = TRUE)
# Farbige Typen
for (r in seq_len(nrow(df1))) {
  col_s <- switch(df1$Verlaufsmuster[r],
    "Zunehmend" = createStyle(fgFill = "#DAE8FC", halign = "center"),
    "Abnehmend" = createStyle(fgFill = "#F8CECC", halign = "center"),
    "Stabil"    = createStyle(fgFill = "#D5E8D4", halign = "center"),
    NULL)
  if (!is.null(col_s)) addStyle(wb, "Verlaufsmuster", col_s, rows = r + 1, cols = 7)
}
setColWidths(wb, "Verlaufsmuster", cols = 1:ncol(df1), widths = "auto")

# ── Sheet 2: Nachtaktivität pro Igel × Nacht ────────────────
addWorksheet(wb, "Nacht_Aktivitaet")
df2 <- as.data.frame(dt_nacht[order(igel, tage_seit), .(
  Igel            = igel,
  Datum           = as.character(datum),
  `Tage seit AW`  = tage_seit,
  Phase           = as.character(phase),
  `% aktiv`       = round(pct_pct, 1),
  `N Minuten`     = n_min
)])
writeData(wb, "Nacht_Aktivitaet", df2, headerStyle = h_style)
setColWidths(wb, "Nacht_Aktivitaet", cols = 1:ncol(df2), widths = "auto")

# ── Sheet 3: Aktivitätsbeginn ────────────────────────────────
if (nrow(onset_dt) > 0) {
  addWorksheet(wb, "Aktivitaetsbeginn")
  df3 <- as.data.frame(onset_dt[order(igel, tage_seit), .(
    Igel            = igel,
    Datum           = as.character(datum),
    `Tage seit AW`  = tage_seit,
    Phase           = ifelse(tage_seit <= phase1_bis, "Release", "Etablierung"),
    `Beginn (Uhr)`  = round(onset_h, 2),
    `Min nach SU`   = round(onset_min, 0),
    `Sonnenuntergang (Uhr)` = round(untergang_h, 2)
  )])
  writeData(wb, "Aktivitaetsbeginn", df3, headerStyle = h_style)
  setColWidths(wb, "Aktivitaetsbeginn", cols = 1:ncol(df3), widths = "auto")
}

# ── Sheet 4: Variabilität pro Tag ────────────────────────────
addWorksheet(wb, "Variabilitaet_pro_Tag")
df4 <- as.data.frame(variab[order(tage_seit), .(
  `Tage seit AW`   = tage_seit,
  Phase            = ifelse(tage_seit <= phase1_bis, "Release", "Etablierung"),
  `N Igel`         = n,
  `Mittel % aktiv` = round(m, 1),
  `SD (PP)`        = round(sd, 1)
)])
writeData(wb, "Variabilitaet_pro_Tag", df4, headerStyle = h_style)
setColWidths(wb, "Variabilitaet_pro_Tag", cols = 1:ncol(df4), widths = "auto")

# ── Sheet 5: Kennzahlen ───────────────────────────────────────
addWorksheet(wb, "Kennzahlen")
onset_med_p1 <- if (nrow(onset_dt) > 0)
  paste0(round(median(onset_dt[tage_seit <= phase1_bis, onset_min], na.rm=TRUE), 0), " min") else "–"
onset_med_p2 <- if (nrow(onset_dt) > 0)
  paste0(round(median(onset_dt[tage_seit >  phase1_bis, onset_min], na.rm=TRUE), 0), " min") else "–"

df5 <- data.frame(
  Kennzahl = c(
    "Analysefenster",
    "Klassifikationsbasis",
    "Eingeschlossene Igel",
    "Igel mit Daten in beiden Phasen",
    "Release-Phase",
    "Etablierungsphase",
    "Mittl. Nachtaktivität Phase 1 (%)",
    "Mittl. Nachtaktivität Phase 2 (%)",
    "Mediandifferenz Phase1 – Phase2 (PP)",
    "Wilcoxon p-Wert",
    "Interpretation Wilcoxon",
    "SD Phase 1 (PP)",
    "SD Phase 2 (PP)",
    "Aktivitätsbeginn nach SU — Phase 1",
    "Aktivitätsbeginn nach SU — Phase 2",
    "Verlaufsmuster: Zunehmend (N)",
    "Verlaufsmuster: Stabil (N)",
    "Verlaufsmuster: Abnehmend (N)"
  ),
  Wert = c(
    paste0("Tag 1–", cutoff_tag, " (N ≥ ", n_min_igel, " Igel)"),
    "pred_nested_loio_smoothed_wmv (alle Igel inkl. 1–6)",
    as.character(length(igel_ok)),
    as.character(length(igel_beide)),
    paste0("Tag 1–", phase1_bis),
    paste0("Tag ", phase1_bis + 1, " bis ", cutoff_tag),
    paste0(round(mean(dt_nacht[tage_seit <= phase1_bis, pct_pct], na.rm=TRUE), 1)),
    paste0(round(mean(dt_nacht[tage_seit >  phase1_bis, pct_pct], na.rm=TRUE), 1)),
    as.character(round(med_diff, 1)),
    format.pval(wt$p.value, digits = 3),
    ifelse(wt$p.value < 0.05,
           "Signifikanter Unterschied zwischen Phasen",
           "Kein signifikanter Unterschied"),
    as.character(round(variab[tage_seit <= phase1_bis, mean(sd, na.rm=TRUE)], 1)),
    as.character(round(variab[tage_seit >  phase1_bis, mean(sd, na.rm=TRUE)], 1)),
    onset_med_p1,
    onset_med_p2,
    as.character(norm_dt[typ == "Zunehmend", .N]),
    as.character(norm_dt[typ == "Stabil",    .N]),
    as.character(norm_dt[typ == "Abnehmend", .N])
  ),
  stringsAsFactors = FALSE
)
writeData(wb, "Kennzahlen", df5, headerStyle = h_style)
for (r in seq(2, nrow(df5) + 1, 2))
  addStyle(wb, "Kennzahlen", z_style, rows = r, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "Kennzahlen", cols = 1:2, widths = c(45, 30))

excel_pfad <- file.path(output_ordner, "Block2d_ReleaseEffect_Uebersicht.xlsx")
saveWorkbook(wb, excel_pfad, overwrite = TRUE)
cat("✓ Excel gespeichert\n\n")

# ══════════════════════════════════════════════════════════════
# WORD-METHODENBERICHT
# ══════════════════════════════════════════════════════════════

cat("Erstelle Word-Methodenbericht...\n")

library(officer)
library(flextable)

add_img <- function(doc, pfad, b = 15, h = 9) {
  if (file.exists(pfad))
    body_add_img(doc, pfad, width = b / 2.54, height = h / 2.54)
  else
    body_add_par(doc, paste0("[Bild fehlt: ", basename(pfad), "]"), style = "Normal")
}

# Flextable-Hilfsfunktion
make_ft <- function(df) {
  flextable(as.data.frame(df)) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#2C3E6B") |>
    color(part = "header", color = "white") |>
    bg(i = seq(1, nrow(df), 2), bg = "#F0F4FA") |>
    fontsize(size = 9, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    autofit()
}

doc <- read_docx() |>
  body_set_default_section(prop_section(
    page_size    = page_size(width = 21 / 2.54, height = 29.7 / 2.54),
    page_margins = page_mar(top = 2.5, bottom = 2.5, left = 3, right = 2.5)
  ))

# ── Titelseite ────────────────────────────────────────────────
doc <- doc |>
  body_add_par("Block 2d — Release-Effekt-Analyse", style = "heading 1") |>
  body_add_par("Methodenbericht", style = "heading 2") |>
  body_add_par(paste("Erstellt:", format(Sys.time(), "%d.%m.%Y %H:%M")), style = "Normal") |>
  body_add_par(paste0("Projekt: Igelbesenderung Niedersachsen — Wildtierstation Sachsenhagen"),
               style = "Normal") |>
  body_add_break()

# ── 1. Hintergrund & Fragestellung ───────────────────────────
doc <- doc |>
  body_add_par("1  Hintergrund und Fragestellung", style = "heading 1") |>
  body_add_par(paste0(
    "In der Wildtiertelemetrie ist der sogenannte 'Release-Effekt' ein bekanntes Phänomen: ",
    "Tiere zeigen direkt nach der Freilassung häufig atypisches Verhalten, das durch Stress, ",
    "Desorientierung in einem unbekannten Gebiet oder die veränderte Situation nach der ",
    "Gefangenschaft bedingt sein kann. Für rehabilitierte Wildtiere ist dieser Effekt ",
    "besonders relevant, da die Station wissen muss, ob und wann sich die Tiere auf ein ",
    "stabiles, artgerechtes Verhaltensmuster einpendeln."),
    style = "Normal") |>
  body_add_par(paste0(
    "Die vorliegende Analyse (Block 2d) untersucht, ob sich die nächtliche Aktivität der ",
    "besenderten Igel in den ersten ", phase1_bis, " Tagen nach Auswilderung (Release-Phase) ",
    "signifikant von den folgenden Tagen (Etablierungsphase: Tag ", phase1_bis + 1,
    "+) unterscheidet. Zusätzlich wird der Zeitpunkt des abendlichen Aktivitätsbeginns ",
    "relativ zum Sonnenuntergang, die Variabilität zwischen Individuen und die ",
    "Klassifikation der Verlaufsmuster je Tier analysiert."),
    style = "Normal")

# ── 2. Datenbasis ─────────────────────────────────────────────
doc <- doc |>
  body_add_break() |>
  body_add_par("2  Datenbasis und Klassifikationsmethode", style = "heading 1") |>
  body_add_par("2.1  Verwendete Dateien", style = "heading 2") |>
  body_add_par(paste0(
    "Block 2d verwendet zwei Datensätze, die von Block0_Datenpipeline.R erzeugt werden:"),
    style = "Normal")

df_data <- data.frame(
  Datei = c("gamm_tagesdaten.rds", "chrono_minuten24h.rds"),
  Inhalt = c(
    "1 Zeile pro Igel × Datum × Tageszeit; enthält Anzahl Minuten (n_min) und Anteil aktiver Minuten (pct_aktiv, 0–1)",
    "1 Zeile pro Igel × Minute (24h); enthält stunde (kontinuierlich 0–24) und aktiv (0/1)"
  ),
  `Verwendung in Block 2d` = c(
    "Analysen 1, 3, 4 (Nachtaktivität, Variabilität, Verlaufsmuster)",
    "Analyse 2 (Aktivitätsbeginn relativ zu Sonnenuntergang)"
  ),
  stringsAsFactors = FALSE
)
doc <- body_add_flextable(doc, make_ft(df_data))

doc <- doc |>
  body_add_par("2.2  Klassifikationsmethode", style = "heading 2") |>
  body_add_par(paste0(
    "Die Aktivitätsklassifikation aller Igel — einschließlich Igel 1–6 — basiert auf der Spalte ",
    "'pred_nested_loio_smoothed_wmv' aus den Rohdaten-CSVs. Diese Klassifikation verwendet ",
    "die gewichtete Bewegungsvarianz (Weighted Movement Variance, WMV) des Sendesignals als ",
    "Indikator für körperliche Aktivität. Über ein gleitendes 5-Minuten-Fenster wird das ",
    "Signal geglättet (smoothed WMV), bevor es mittels eines Leave-One-Individual-Out-Ansatzes ",
    "(LOIO) in 'aktiv' (a) und 'passiv' (p) klassifiziert wird. ",
    "Diese Methode reduziert Fehlklassifikationen durch kurzzeitige Signalschwankungen und ist ",
    "für alle Tiere einheitlich angewendet worden."),
    style = "Normal") |>
  body_add_par(paste0(
    "Hinweis: In Block 2c (HMM-Analyse) wird für die Validierungsplots ein separater WMV-",
    "Referenzdatensatz geladen. Da für Igel 1–6 das Format dieser Referenzdatei abweicht, ",
    "wird dort auf 'pred_nested_loio_mv' zurückgegriffen. Diese abweichende Behandlung ",
    "betrifft ausschließlich die Vergleichsplots in Block 2c — nicht die Grundklassifikation, ",
    "die in allen anderen Blöcken (2, 2b, 2d, 3) verwendet wird."),
    style = "Normal") |>
  body_add_par("2.3  Qualitätskontrolle und Zeitfenster", style = "heading 2") |>
  body_add_par(paste0(
    "Die QC-Kriterien sind identisch mit Block 2 (GAMM): Ausschluss bei weniger als ",
    min_naechte, " Nächten mit Nachtdaten, weniger als ", min_nacht_min,
    " absoluten Nachtminuten, oder unter ", min_pct_aktiv,
    "% mittlerer Nachtaktivität (statisches Signal / toter Sender). ",
    "Das Analysefenster endet dynamisch am letzten Tag, an dem noch mindestens ",
    n_min_igel, " Igel Daten liefern (N-at-risk-Cutoff = Tag ", cutoff_tag, ")."),
    style = "Normal")

# ── 3. Methoden ───────────────────────────────────────────────
doc <- doc |>
  body_add_break() |>
  body_add_par("3  Statistische Methoden", style = "heading 1") |>
  body_add_par("3.1  Phasenvergleich: Wilcoxon Signed-Rank Test", style = "heading 2") |>
  body_add_par(paste0(
    "Für jeden Igel wird die mittlere nächtliche Aktivität (% aktiver Minuten) in der ",
    "Release-Phase (Tag 1–", phase1_bis, ") und der Etablierungsphase (Tag ",
    phase1_bis + 1, "+) berechnet. Da dieselben Tiere in beiden Phasen beobachtet werden, ",
    "wird ein paarweiser Wilcoxon Signed-Rank Test (Wilcoxon 1945) angewendet. ",
    "Dieser Test ist nicht-parametrisch und eignet sich für kleine Stichprobengrößen und ",
    "nicht normalverteilte Daten. Er prüft, ob die Mediane beider Phasen systematisch ",
    "verschieden sind (H0: kein Unterschied)."),
    style = "Normal") |>
  body_add_par("3.2  Aktivitätsbeginn relativ zu Sonnenuntergang", style = "heading 2") |>
  body_add_par(paste0(
    "Für jede Nacht und jeden Igel wird der erste Zeitpunkt (stunde) ermittelt, zu dem der ",
    "Igel nach Sonnenuntergang aktiv ist (aktiv = 1 und stunde ≥ Sonnenuntergang − 1h). ",
    "Der Sonnenuntergang wird tagesgenau über das R-Paket 'suncalc' für den Standort ",
    "Sachsenhagen (Breitengrad ", lat_h, "°N, Längengrad ", lon_h, "°E) berechnet. ",
    "Das Ergebnis wird in Minuten nach Sonnenuntergang ausgedrückt: positiv = nach Sonnenuntergang aktiv, ",
    "negativ = bereits vor Sonnenuntergang aktiv. Werte unter −60 Minuten werden als Ausreißer entfernt."),
    style = "Normal") |>
  body_add_par("3.3  Variabilität zwischen Individuen", style = "heading 2") |>
  body_add_par(paste0(
    "Die Standardabweichung (SD) der Nachtaktivität zwischen allen Igeln wird pro Tag berechnet. ",
    "Eine hohe SD zeigt große Unterschiede zwischen Individuen, eine niedrige SD zeigt ",
    "ähnliches Verhalten aller Tiere. Ein abnehmender Trend würde auf eine Konvergenz ",
    "der Verhaltensmuster im Verlauf der Beobachtungszeit hinweisen."),
    style = "Normal") |>
  body_add_par("3.4  Verlaufsmuster-Klassifikation", style = "heading 2") |>
  body_add_par(paste0(
    "Jedes Tier wird nach der Differenz (Δ) zwischen seiner mittleren Nachtaktivität in ",
    "Phase 1 und Phase 2 klassifiziert: 'Zunehmend' wenn Phase 2 um mehr als 15 Prozentpunkte ",
    "höher liegt (Δ < −15), 'Abnehmend' wenn Phase 1 um mehr als 15 Prozentpunkte höher liegt ",
    "(Δ > 15), und 'Stabil' bei einer Differenz ≤ 15 Prozentpunkte. Der Schwellenwert von ",
    "15 Prozentpunkten entspricht einem biologisch relevanten Unterschied bei VHF-Telemetrie-Daten."),
    style = "Normal")

# ── 4. Ergebnisse ─────────────────────────────────────────────
doc <- doc |>
  body_add_break() |>
  body_add_par("4  Ergebnisse", style = "heading 1") |>
  body_add_par("4.1  Kennzahlen-Übersicht", style = "heading 2") |>
  body_add_flextable(make_ft(df5)) |>
  body_add_par("4.2  Nachtaktivität über die Zeit (Spaghetti-Plot)", style = "heading 2") |>
  body_add_par(paste0(
    "Der Spaghetti-Plot zeigt die individuelle Nachtaktivität aller ", length(igel_ok),
    " Igel über den gesamten Beobachtungszeitraum. Die orange markierte Zone umfasst die ",
    "Release-Phase (Tag 1–", phase1_bis, "). Die schwarze Linie zeigt den Populationsmittelwert ± SE. ",
    "N (unten) gibt die Anzahl Igel mit Daten pro Tag an."),
    style = "Normal")
doc <- add_img(doc, file.path(output_ordner, "p1_spaghetti.png"), b = 17, h = 9)

doc <- doc |>
  body_add_par("4.3  Paarvergleich Release-Phase vs. Etablierungsphase", style = "heading 2") |>
  body_add_par(paste0(
    "Jede Linie verbindet denselben Igel in beiden Phasen. Blaue Linien zeigen Tiere, ",
    "deren Nachtaktivität in der Etablierungsphase höher war; rote Linien zeigen das Gegenteil. ",
    "Die Rauten markieren die Mediane. ",
    "Wilcoxon Signed-Rank Test (N = ", nrow(a1_wide), " Igel): ",
    "Mediandifferenz = ", round(med_diff, 1), " Prozentpunkte, ",
    "p = ", format.pval(wt$p.value, digits = 3), ". ",
    ifelse(wt$p.value < 0.05,
           "Die Release-Phase unterscheidet sich signifikant von der Etablierungsphase.",
           "Es besteht kein signifikanter Unterschied zwischen den Phasen.")),
    style = "Normal")
doc <- add_img(doc, file.path(output_ordner, "p2_paarvergleich.png"), b = 12, h = 11)
doc <- body_add_break(doc)

if (file.exists(file.path(output_ordner, "p3_aktivitaetsbeginn.png"))) {
  doc <- doc |>
    body_add_par("4.4  Aktivitätsbeginn relativ zu Sonnenuntergang", style = "heading 2") |>
    body_add_par(paste0(
      "Dargestellt ist der erste abendliche Aktivitätszeitpunkt je Igel und Nacht in Minuten ",
      "nach Sonnenuntergang. Der Median der Release-Phase lag bei ", onset_med_p1, " (nach SU), ",
      "der Etablierungsphase bei ", onset_med_p2, ". ",
      "Der Nullpunkt entspricht dem exakten Sonnenuntergang (tagesgenau berechnet via suncalc). ",
      "Negative Werte bedeuten: Igel war bereits vor Sonnenuntergang aktiv."),
      style = "Normal")
  doc <- add_img(doc, file.path(output_ordner, "p3_aktivitaetsbeginn.png"), b = 16, h = 9)
}

doc <- doc |>
  body_add_par("4.5  Variabilität zwischen Individuen", style = "heading 2") |>
  body_add_par(paste0(
    "Die SD der Nachtaktivität zwischen allen Igeln pro Tag zeigt, ob die Streuung im Verlauf ",
    "der Beobachtungszeit ab- oder zunimmt. In der Release-Phase betrug die mittlere SD ",
    round(variab[tage_seit <= phase1_bis, mean(sd, na.rm=TRUE)], 1),
    " Prozentpunkte, in der Etablierungsphase ",
    round(variab[tage_seit >  phase1_bis, mean(sd, na.rm=TRUE)], 1),
    " Prozentpunkte."),
    style = "Normal")
doc <- add_img(doc, file.path(output_ordner, "p4_variabilitaet.png"), b = 16, h = 8)
doc <- body_add_break(doc)

doc <- doc |>
  body_add_par("4.6  Verlaufsmuster je Individuum", style = "heading 2") |>
  body_add_par(paste0(
    "Von ", nrow(norm_dt), " Igeln mit ausreichend Daten in beiden Phasen wurden ",
    norm_dt[typ == "Zunehmend", .N], " als 'Zunehmend' (Aktivität steigt ab Tag 4), ",
    norm_dt[typ == "Stabil",    .N], " als 'Stabil' und ",
    norm_dt[typ == "Abnehmend", .N], " als 'Abnehmend' klassifiziert. ",
    "Jeder Datenpunkt entspricht einem Igel. Die gestrichelte Diagonale zeigt den Bereich ",
    "ohne Unterschied zwischen den Phasen; die gepunkteten Linien markieren die ±15 PP-Schwellen."),
    style = "Normal")
doc <- add_img(doc, file.path(output_ordner, "p5_verlaufsmuster.png"), b = 14, h = 12)

# ── 5. Interpretation ─────────────────────────────────────────
doc <- doc |>
  body_add_break() |>
  body_add_par("5  Biologische Interpretation und Empfehlungen", style = "heading 1") |>
  body_add_par(paste0(
    "Ein signifikanter Release-Effekt (p < 0.05) würde bedeuten, dass die Tiere direkt nach ",
    "der Auswilderung systematisch mehr oder weniger aktiv sind als in den Folgewochen — ",
    "ein Hinweis auf Stress, Orientierungslosigkeit oder veränderte Aktivitätsmuster ",
    "unmittelbar post-release. Ein nicht signifikantes Ergebnis zeigt, dass die Tiere ",
    "bereits in der ersten Nacht nach der Auswilderung ein ähnliches Aktivitätsniveau zeigen ",
    "wie in den Folgewochen, was auf eine schnelle Adaptation hindeutet."),
    style = "Normal") |>
  body_add_par(paste0(
    "Die Variabilität zwischen Individuen (Analyse 3) gibt Aufschluss darüber, ob einzelne ",
    "Tiere mit der Auswilderung 'Ausreißer' sind oder ob das gesamte Kollektiv ähnlich reagiert. ",
    "Die Verlaufsmuster-Klassifikation (Analyse 4) erlaubt es der Station, Tiere zu identifizieren, ",
    "die ihren Aktivitätslevel nach der Release-Phase noch deutlich verändern — ",
    "potenziell ein Indikator für Rehabilitationserfolg oder -bedarf."),
    style = "Normal") |>
  body_add_par(paste0(
    "Der Aktivitätsbeginn relativ zu Sonnenuntergang (Analyse 2) ist chronobiologisch besonders ",
    "relevant: gesunde Igel sollten kurz nach Sonnenuntergang aktiv werden. Wenn Tiere in der ",
    "Release-Phase deutlich später oder früher aktiv sind als in der Etablierungsphase, ",
    "kann das auf circadiane Desynchronisation durch die Stresssituation der Freilassung ",
    "oder durch veränderte Lichtverhältnisse im Freiland vs. Station hinweisen."),
    style = "Normal")

# Speichern
doc_pfad <- file.path(output_ordner, "Block2d_ReleaseEffect_Methodenbericht.docx")
print(doc, target = doc_pfad)
cat("✓ Word-Bericht gespeichert\n\n")

cat("══════════════════════════════════════════\n")
cat("✓ Block 2d abgeschlossen!\n")
cat("Dateien in:", output_ordner, "\n")
cat("  p1_spaghetti.png\n")
cat("  p2_paarvergleich.png\n")
cat("  p3_aktivitaetsbeginn.png\n")
cat("  p4_variabilitaet.png\n")
cat("  p5_verlaufsmuster.png\n")
cat("  Block2d_ReleaseEffect_Uebersicht.xlsx\n")
cat("  Block2d_ReleaseEffect_Methodenbericht.docx\n")
cat("══════════════════════════════════════════\n")
