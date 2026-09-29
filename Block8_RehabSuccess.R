# ==============================================================================
# Block 8 — rehabilitation-success synthesis
#
# Question: Which factors influence behavioural normalisation after release?
#           Do hedgehogs reach nocturnal activity patterns like wild animals?
#
# PRIMARY METRIC: night-activity proportion (analysis E)
#   Proportion of activity during the astronomical night (sunset to sunrise,
#   Sachsenhagen 52.397 N 9.217 E). Direct, model-free measure. Robust to the
#   ~20-minute smoothing in pred_nested_loio_smoothed_wmv.
#
# SECONDARY METRIC: IS — Interdaily Stability (analyses A-D)
#   Measures synchronisation of the 24h rhythm. Significant correlation with
#   weight gain (rho=0.57, p=0.026) justifies use as a secondary measure.
#   Only animals with >= MIN_TAGE_IS_IV days included.
#
# IV — Intradaily Variability: NOT USED
#   IV was excluded after the empirical comparison (Block0_wmv_Vergleich.R)

# ── 0. Pakete ──────────────────────────────────────────────────────────────────

pakete <- c("data.table", "ggplot2", "lubridate", "scales", "patchwork",
            "readxl", "openxlsx",   # openxlsx statt writexl (einheitliches Styling)
            "mgcv",                  # GAM-Smoother (geom_smooth method="gam")
            "officer", "flextable")  # Word-Berichte

neu <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) {
  cat("Installiere fehlende Pakete:", paste(neu, collapse = ", "), "\n")
  install.packages(neu)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(lubridate)
  library(scales)
  library(patchwork)
  library(readxl)
  library(openxlsx)
  library(mgcv)
})
cat("✓ Pakete geladen\n\n")

# ── Einstellungen ──────────────────────────────────────────────────────────────

projekt_root  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
output_ordner <- file.path(projekt_root, "output", "Block8_RehabSuccess")
dir.create(output_ordner, showWarnings = FALSE, recursive = TRUE)

is_rds     <- file.path(projekt_root, "output", "Block5_IS_IV", "is_iv_ergebnisse.rds")
chrono_rds <- file.path(projekt_root, "output", "Block0_Pipeline", "chrono_minuten24h.rds")
meta_pfad  <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")

epoch_min       <- 60L
is_fenster_tage <- 7L
min_tage_window <- 3L

# Mindest-Tracking-Dauer für IS/IV-Analysen
MIN_TAGE_IS_IV  <- 7L

# Referenz Nachtaktivität
# Koordinaten Studiengebiet Sachsenhagen (identisch mit Block3_Chronobiologie.R)
STUDIE_LAT <- 52.397
STUDIE_LON <-  9.217

# Nachtfenster: astronomisch (per suncalc, berechnet nach Laden der Minutendaten)
# Kein fixer NACHT_START/NACHT_ENDE mehr — wird als dt$ist_nacht gesetzt.
REFERENZ_NACHT <- 0.85   # interner Schwellenwert (kein externer Literaturwert — nur für interne Bewertungslogik)

# Mindest-Gesamtaktivität für Winterschlaf-Filter
# Fenster mit Nacht+Tag-Aktivität < 20% gelten als Torpor/Winterschlaf und werden
# aus der Primärmetrik-Analyse ausgeschlossen.
# Konsequenz: Igel1 (3 Tracking-Tage, <10% Gesamtaktivität) vollständig ausgeschlossen;
#             Igel5: Winterschlaf-Phasen (Tage ~26–165 und ~211+) werden gefiltert.
AKTIV_SCHWELLE <- 0.20

# ── Farben ─────────────────────────────────────────────────────────────────────
COL_IGEL    <- "#5F9E6E"
COL_BLOCK8  <- "#4C6A9C"   # Blau als Block-8-Akzentfarbe (Berichte)

diag_farben <- c(
  "Waisling"       = "#4393C3",
  "Parasiten"      = "#74C476",
  "Trauma"         = "#F4A582",
  "Pilzerkrankung" = "#DFC27D",
  "Blindheit"      = "#D6604D",
  "Sonstige"       = "#B2ABD2"
)

sex_farben <- c("Male" = "#4393C3", "Female" = "#D6604D")

cat("══════════════════════════════════════════════════════════════\n")
cat("  Block 8 — Rehabilitations-Erfolgs-Synthese\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 1. DATEN LADEN
# ══════════════════════════════════════════════════════════════════════════════

if (!file.exists(is_rds))
  stop("is_iv_ergebnisse.rds nicht gefunden. Bitte Block5_IS_IV.R ausführen.")
is_ergebnisse <- readRDS(is_rds)
is_gesamt     <- as.data.table(is_ergebnisse$gesamt)
is_zeit       <- as.data.table(is_ergebnisse$zeitverlauf)
cat("IS/IV Gesamt:", nrow(is_gesamt), "Tiere |",
    nrow(is_gesamt[!is.na(IS)]), "mit gültigem IS\n")
cat("IS/IV Zeitverlauf:", nrow(is_zeit[!is.na(IS)]), "Datenpunkte\n")

if (!file.exists(chrono_rds))
  stop("chrono_minuten24h.rds nicht gefunden. Bitte Block0_Pipeline.R ausführen.")
dt <- readRDS(chrono_rds)
setDT(dt)
cat("Aktivitätsdaten:", format(nrow(dt), big.mark = "'"), "Minuten\n")

# ── Astronomisches Nachtfenster per Datum (suncalc, Sachsenhagen) ─────────────
if (!requireNamespace("suncalc", quietly = TRUE)) install.packages("suncalc")
library(suncalc)

datums_vec <- sort(unique(dt$datum))
sonnen_dt  <- as.data.table(getSunlightTimes(
  date = datums_vec,
  lat  = STUDIE_LAT,
  lon  = STUDIE_LON,
  keep = c("sunrise", "sunset"),
  tz   = "Europe/Berlin"
))
sonnen_dt[, datum     := as.Date(date)]
sonnen_dt[, sunrise_h := as.numeric(format(sunrise, "%H", tz = "Europe/Berlin")) +
                         as.numeric(format(sunrise, "%M", tz = "Europe/Berlin")) / 60]
sonnen_dt[, sunset_h  := as.numeric(format(sunset,  "%H", tz = "Europe/Berlin")) +
                         as.numeric(format(sunset,  "%M", tz = "Europe/Berlin")) / 60]
sonnen_dt <- sonnen_dt[, .(datum, sunrise_h, sunset_h)]

dt <- merge(dt, sonnen_dt, by = "datum", all.x = TRUE)
# ist_nacht = TRUE wenn stunde >= Sonnenuntergang ODER stunde < Sonnenaufgang
dt[, ist_nacht := stunde >= sunset_h | stunde < sunrise_h]

cat(sprintf(
  "Astronomisches Nachtfenster: Ø SU %.1f Uhr, Ø SA %.1f Uhr (Sachsenhagen, n=%d Daten)\n",
  mean(sonnen_dt$sunset_h, na.rm = TRUE),
  mean(sonnen_dt$sunrise_h, na.rm = TRUE),
  length(datums_vec)
))

# ── Metadaten aus Excel laden ─────────────────────────────────────────────────
meta_raw <- as.data.table(read_excel(meta_pfad))
setnames(meta_raw, names(meta_raw), tolower(trimws(gsub("\\s+", "_", names(meta_raw)))))

meta_voll <- data.table(
  igel           = trimws(as.character(meta_raw[["individual"]])),
  sex            = trimws(as.character(meta_raw[["sex"]])),
  diagnosis      = trimws(as.character(meta_raw[["diagnosis_main"]])),
  time_reha      = suppressWarnings(as.numeric(meta_raw[["time_reha"]])),
  tagging_weight = suppressWarnings(as.numeric(meta_raw[["tagging_weight"]])),
  weight_entry   = suppressWarnings(as.numeric(meta_raw[["weigth_entry"]])),
  weight_gain    = suppressWarnings(as.numeric(meta_raw[["weight_gain"]])),
  tagging_period = suppressWarnings(as.numeric(meta_raw[["tagging_period"]])),
  date_release   = suppressWarnings(as.Date(meta_raw[["date_release"]]))
)

meta_voll[!grepl("^Igel", igel, ignore.case = TRUE), igel := paste0("Igel", igel)]
meta_voll[, alter := ifelse(!is.na(weight_entry) & weight_entry < 300, "Jungtier", "Adult")]
meta_voll[, monat_auswild := month(date_release)]
meta_voll[, saison_auswild := fcase(
  monat_auswild %in% c(3, 4, 5),   "Frühling",
  monat_auswild %in% c(6, 7, 8),   "Sommer",
  monat_auswild %in% c(9, 10, 11), "Herbst",
  default = NA_character_
)]
meta_voll[, c("date_release", "monat_auswild") := NULL]

# ── Diagnose-Kategorisierung ──────────────────────────────────────────────────
meta_voll[, diag_kat := fcase(
  grepl("blind",     diagnosis, ignore.case = TRUE),                          "Blindheit",
  grepl("orph",      diagnosis, ignore.case = TRUE),                          "Waisling",
  grepl("fungal",    diagnosis, ignore.case = TRUE),                          "Pilzerkrankung",
  grepl("endoparasit", diagnosis, ignore.case = TRUE),                        "Endoparasiten",
  grepl("parasite",  diagnosis, ignore.case = TRUE) &
    !grepl("blind|orph|endo", diagnosis, ignore.case = TRUE),                 "Parasiten",
  grepl("trauma",    diagnosis, ignore.case = TRUE) &
    !grepl("blind",  diagnosis, ignore.case = TRUE),                          "Trauma",
  default = "Sonstige"
)]
meta_voll[, diag_kat := factor(diag_kat,
  levels = c("Waisling", "Parasiten", "Trauma", "Pilzerkrankung", "Blindheit", "Sonstige"))]

cat("\nDiagnose-Kategorien:\n")
print(meta_voll[, .N, by = diag_kat][order(diag_kat)])
cat("\n")

# ── IS/IV mit Metadaten verknüpfen ────────────────────────────────────────────
meta_join <- meta_voll[, .(igel, sex, alter, diag_kat, time_reha,
                            tagging_weight, weight_gain, tagging_period, saison_auswild)]

is_voll <- merge(
  is_gesamt[, .(igel, n_tage, IS, IV, RA)],
  meta_join, by = "igel", all.x = TRUE
)
is_zeit_voll <- merge(
  is_zeit[, .(igel, tage_seit, IS, IV, n_tage)],
  meta_join, by = "igel", all.x = TRUE
)

cat("IS-Datensatz (gesamt):", nrow(is_voll), "Tiere |",
    nrow(is_voll[!is.na(IS)]), "mit gültigem IS\n\n")

# ── Mindest-Tracking-Filter ───────────────────────────────────────────────────
is_voll_analyse <- is_voll[!is.na(n_tage) & n_tage >= MIN_TAGE_IS_IV]
n_ausgeschlossen <- nrow(is_voll) - nrow(is_voll_analyse)

# ── IS-Artefakt-Flag: Tiere mit unzuverlässigem globalen IS ──────────────────
# Igel5: Zwei Winterschlaf-Perioden (Tage ~26–165 und ~211+) über 221 Tracking-Tage.
#        Der globale IS (0.094) spiegelt hauptsächlich die Winterschlaf-Dynamik wider,
#        nicht die circadiane Stabilität im aktiven Zustand.
# Igel1: Nur 3 Tracking-Tage mit < 5% Gesamtaktivität → IS (0.072) statistisch instabil.
# Diese Tiere bleiben im IS-Zeitverlauf (Analyse C), werden aber aus den
# Gruppenvergleichen (A, B, D) ausgeschlossen.
IS_ARTEFAKT_TIERE <- c("Igel5", "Igel1")
is_voll_analyse_clean <- is_voll_analyse[!(igel %in% IS_ARTEFAKT_TIERE)]
cat(sprintf("IS-Artefakt-Ausschluss (Analyse A/B/D): %s\n",
            paste(IS_ARTEFAKT_TIERE, collapse = ", ")))

cat(sprintf("IS/IV-Filter (>= %d Tage): %d von %d Tieren eingeschlossen",
            MIN_TAGE_IS_IV, nrow(is_voll_analyse), nrow(is_voll)))
if (n_ausgeschlossen > 0) {
  ausgeschl_tiere <- is_voll[is.na(n_tage) | n_tage < MIN_TAGE_IS_IV, igel]
  cat(sprintf(" | %d ausgeschlossen: %s",
              n_ausgeschlossen, paste(ausgeschl_tiere, collapse = ", ")))
}
cat("\n")
cat(sprintf("IS-Artefakt-Ausschluss (Analysen A/B/D): %d Tiere (%s)\n  Grund: Igel5 = Winterschlaf über gesamte Tracking-Dauer; Igel1 = 3 Tage/< 5%% Aktivität\n",
            length(IS_ARTEFAKT_TIERE), paste(IS_ARTEFAKT_TIERE, collapse = ", ")))
cat(sprintf("IS-Analyse-Stichprobe (bereinigt): %d Tiere\n\n", nrow(is_voll_analyse_clean)))

# ══════════════════════════════════════════════════════════════════════════════
# ROLLING NACHTAKTIVITÄTS-ANTEIL (für Analyse E)
# ══════════════════════════════════════════════════════════════════════════════
cat("Berechne Rolling Nachtaktivitäts-Anteil …\n")

igel_liste <- dt[, levels(droplevels(factor(igel)))]

nacht_zeit <- rbindlist(lapply(igel_liste, function(ig) {
  sub  <- dt[igel == ig & !is.na(aktiv)]
  if (sub[, uniqueN(datum)] < min_tage_window) return(NULL)
  tage <- sort(unique(sub$tage_seit))
  halb <- floor(is_fenster_tage / 2L)

  rbindlist(lapply(tage, function(t) {
    sub_f  <- sub[tage_seit >= (t - halb) & tage_seit <= (t + halb)]
    n_t    <- sub_f[, uniqueN(datum)]
    if (n_t < min_tage_window) return(NULL)
    gesamt <- sub_f[, .N]
    nacht   <- sub_f[ist_nacht == TRUE,  sum(aktiv, na.rm = TRUE)]
    tag_akt <- sub_f[ist_nacht == FALSE, sum(aktiv, na.rm = TRUE)]
    if (gesamt == 0L) return(NULL)
    data.table(
      igel         = ig,
      tage_seit    = t,
      n_tage       = n_t,
      nacht_anteil = nacht / gesamt,
      tag_anteil   = tag_akt / gesamt
    )
  }), fill = TRUE)
}), fill = TRUE)

nacht_voll <- merge(nacht_zeit, meta_join, by = "igel", all.x = TRUE)

# ── Abgeleitete Metriken ──────────────────────────────────────────────────────
# gesamt_aktiv = Anteil aktiver Minuten (Tag+Nacht) an allen detektierten Minuten
# nacht_von_aktiv = Anteil der AKTIVEN Zeit der nachts stattfindet
#   → im Gegensatz zu nacht_anteil (= aktive Nachtmin / ALLE detektierten Min)
#     enthält nacht_von_aktiv keine passiven Ruhephasen im Nenner
nacht_voll[, gesamt_aktiv    := nacht_anteil + tag_anteil]
nacht_voll[, nacht_von_aktiv := fifelse(gesamt_aktiv > 0.01,
                                         nacht_anteil / gesamt_aktiv,
                                         NA_real_)]

# Winterschlaf-Flag: Gesamtaktivität < AKTIV_SCHWELLE → Torpor/Hibernation-Phase
nacht_voll[, winterschlaf := gesamt_aktiv < AKTIV_SCHWELLE]

# Gefilterte Stichprobe für Primärmetrik-Analysen
nacht_voll_aktiv <- nacht_voll[winterschlaf == FALSE]

cat("Nachtaktivitäts-Datenpunkte (gesamt):", nrow(nacht_voll[!is.na(nacht_anteil)]), "\n")
n_ws_excl <- sum(nacht_voll$winterschlaf, na.rm = TRUE)
cat(sprintf("Winterschlaf-Filter (%d%% Schwelle): %d Fenster ausgeschlossen\n",
            as.integer(AKTIV_SCHWELLE * 100), n_ws_excl))
igel_ws_ausgeschl <- setdiff(unique(nacht_voll$igel), unique(nacht_voll_aktiv$igel))
if (length(igel_ws_ausgeschl) > 0)
  cat("  Vollständig ausgeschlossen:", paste(igel_ws_ausgeschl, collapse = ", "),
      "(< 20% Gesamtaktivität in allen Fenstern)\n")
cat("Datenpunkte nach Filter:     ", nrow(nacht_voll_aktiv), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSE A — IS nach Diagnose-Kategorie
# ══════════════════════════════════════════════════════════════════════════════
cat("════════════════════\nAnalyse A: IS nach Diagnose\n")

# IS-Artefakt-Tiere (Igel5 Winterschlaf, Igel1 Kurztracking) ausgeschlossen
dat_a  <- is_voll_analyse_clean[!is.na(IS) & !is.na(diag_kat)]
n_diag <- dat_a[, .N, by = diag_kat]
kw_IS  <- tryCatch(kruskal.test(IS ~ diag_kat, data = dat_a), error = function(e) NULL)

dat_a <- merge(dat_a, n_diag, by = "diag_kat")
dat_a[, diag_label := paste0(as.character(diag_kat), "\n(n=", N, ")")]
dat_a[, diag_label := factor(diag_label,
  levels = dat_a[, .(diag_kat, diag_label)][order(diag_kat)][!duplicated(diag_kat)]$diag_label)]

p_ann_a <- if (!is.null(kw_IS))
  sprintf("Kruskal-Wallis: χ²=%.2f, df=%d, p=%.3f",
          kw_IS$statistic, kw_IS$parameter, kw_IS$p.value) else ""

p8a <- ggplot(dat_a, aes(x = diag_label, y = IS, fill = diag_kat)) +
  geom_violin(alpha = 0.35, trim = FALSE, linewidth = 0.4,
              data = ~ subset(., diag_kat %in% dat_a[, .N, by = diag_kat][N > 1]$diag_kat)) +
  geom_boxplot(width = 0.25, alpha = 0.7, outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(color = diag_kat), width = 0.08, size = 2.8, alpha = 0.9) +
  scale_fill_manual(values = diag_farben, guide = "none") +
  scale_color_manual(values = diag_farben, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title   = "A  IS nach Diagnose-Kategorie",
       x = NULL, y = "IS — Interdaily Stability (0–1)",
       caption = p_ann_a) +
  theme_minimal(base_size = 11) +
  theme(plot.title         = element_text(face = "bold", size = 12),
        axis.text.x        = element_text(size = 9, lineheight = 1.1),
        plot.caption       = element_text(size = 8, color = "grey40"),
        panel.grid.major.x = element_blank())

ggsave(file.path(output_ordner, "08a_IS_diagnose.png"), p8a,
       width = 8, height = 5.3, dpi = 170, bg = "white")
cat("✓ 08a_IS_diagnose.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSE B — IS nach Geschlecht
# ══════════════════════════════════════════════════════════════════════════════
cat("Analyse B: IS nach Geschlecht\n")

dat_b <- is_voll_analyse_clean[!is.na(IS) & !is.na(sex)]
n_sex <- dat_b[, .N, by = sex]
wt_IS <- tryCatch(wilcox.test(IS ~ sex, data = dat_b), error = function(e) NULL)

sex_n_labels <- setNames(paste0(n_sex$sex, "\n(n=", n_sex$N, ")"), n_sex$sex)
p_ann_b <- if (!is.null(wt_IS))
  sprintf("Wilcoxon: W=%.0f, p=%.3f", wt_IS$statistic, wt_IS$p.value) else ""

p8b <- ggplot(dat_b, aes(x = sex, y = IS, fill = sex)) +
  geom_violin(alpha = 0.35, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(color = sex), width = 0.07, size = 3, alpha = 0.9) +
  scale_fill_manual(values  = sex_farben, guide = "none") +
  scale_color_manual(values = sex_farben, guide = "none") +
  scale_x_discrete(labels = sex_n_labels) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title   = "B  IS nach Geschlecht",
       x = NULL, y = "IS — Interdaily Stability (0–1)",
       caption = p_ann_b) +
  theme_minimal(base_size = 11) +
  theme(plot.title         = element_text(face = "bold", size = 12),
        axis.text.x        = element_text(size = 10),
        plot.caption       = element_text(size = 8, color = "grey40"),
        panel.grid.major.x = element_blank())

ggsave(file.path(output_ordner, "08b_IS_sex.png"), p8b,
       width = 5.3, height = 4.7, dpi = 170, bg = "white")
cat("✓ 08b_IS_sex.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSE C — IS Zeitverlauf nach Diagnose (Spaghetti + GAM)
# ══════════════════════════════════════════════════════════════════════════════
cat("Analyse C: IS Zeitverlauf nach Diagnose\n")

dat_c    <- is_zeit_voll[!is.na(IS) & !is.na(diag_kat)]
igel_n4  <- dat_c[, .N, by = igel][N >= 4]$igel
dat_c    <- dat_c[igel %in% igel_n4]

if (nrow(dat_c) > 0) {
  p8c <- ggplot(dat_c, aes(x = tage_seit, y = IS, color = diag_kat, group = igel)) +
    geom_line(alpha = 0.35, linewidth = 0.6) +
    geom_smooth(aes(group = diag_kat), method = "gam",
                formula = y ~ s(x, bs = "cs", k = 4),
                se = TRUE, linewidth = 1.2, alpha = 0.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    annotate("text", x = max(dat_c$tage_seit) * 0.95, y = 0.52,
             label = "IS = 0.5 (Schwellenwert)", size = 2.8,
             color = "grey40", hjust = 1, fontface = "italic") +
    scale_color_manual(values = diag_farben, name = "Diagnose") +
    scale_x_continuous(name = "Tage seit Auswilderung") +
    scale_y_continuous(name = "IS — Interdaily Stability (0–1)",
                       limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title    = "C  IS-Zeitverlauf nach Diagnose-Kategorie",
         subtitle = "Dünne Linien: Einzeltiere  |  Breite Linien: GAM-Trend pro Gruppe (±SE)",
         caption  = paste0("Gestrichelte Linie: IS = 0.5  |  n = ",
                            length(unique(dat_c$igel)), " Igel")) +
    theme_minimal(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", size = 12),
          plot.subtitle   = element_text(size = 9, color = "grey40"),
          plot.caption    = element_text(size = 8, color = "grey40"),
          legend.position = "right")

  ggsave(file.path(output_ordner, "08c_IS_zeitverlauf_diagnose.png"), p8c,
         width = 10.6, height = 5.3, dpi = 170, bg = "white")
  cat("✓ 08c_IS_zeitverlauf_diagnose.png\n")
} else {
  cat("⚠ Zu wenige Datenpunkte für Zeitverlauf-Plot.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSE D — Reha-Korrelationen
# ══════════════════════════════════════════════════════════════════════════════
cat("Analyse D: Reha-Korrelationen\n")

dat_d <- is_voll_analyse_clean[!is.na(IS)]

korr_vars <- list(
  list(var = "time_reha",      xlab = "Rehab.-Dauer (Tage)",         title = "IS ~ Reha-Dauer"),
  list(var = "tagging_weight", xlab = "Auswilderungsgewicht (g)",    title = "IS ~ Auswilderungsgewicht"),
  list(var = "weight_gain",    xlab = "Gewichtszunahme in Reha (g)", title = "IS ~ Gewichtszunahme")
)

# Bug-Fix: aes_string() ist in neuem ggplot2 deprecated → aes(.data[[xvar]])
plot_korr <- function(dat, xvar, xlab, title) {
  sub <- dat[!is.na(get(xvar))]
  if (nrow(sub) < 4) return(ggplot() + labs(title = paste0(title, "\n(zu wenige Daten)")))
  ct <- tryCatch(cor.test(sub[[xvar]], sub$IS, method = "spearman", exact = FALSE),
                 error = function(e) NULL)
  r_lab <- if (!is.null(ct))
    sprintf("ρ = %.2f\np = %.3f\nn = %d", ct$estimate, ct$p.value, nrow(sub)) else ""

  ggplot(sub, aes(x = .data[[xvar]], y = IS)) +
    geom_point(aes(color = diag_kat), size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                color = "grey30", fill = "grey70", linewidth = 0.9, alpha = 0.25) +
    scale_color_manual(values = diag_farben, name = "Diagnose") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.6,
             label = r_lab, size = 3.2, color = "grey20") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = title, x = xlab, y = "IS (0–1)") +
    theme_minimal(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", size = 10),
          legend.position = "none")
}

plots_d <- lapply(korr_vars, function(v) plot_korr(dat_d, v$var, v$xlab, v$title))

# Bug-Fix: ggpubr guard — plot_layout(guides="collect") reicht, kein ggpubr nötig
p8d <- (plots_d[[1]] | plots_d[[2]] | plots_d[[3]]) + plot_layout(guides = "collect")

ggsave(file.path(output_ordner, "08d_reha_korrelationen.png"), p8d,
       width = 11.8, height = 4.7, dpi = 170, bg = "white")
cat("✓ 08d_reha_korrelationen.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSE E — Nachtaktivitäts-Anteil über Zeit (PRIMÄRMETRIK)
# ══════════════════════════════════════════════════════════════════════════════
cat("Analyse E: Nachtaktivitäts-Normalisierung (Primärmetrik: nacht_von_aktiv)\n")

# Primärmetrik: nacht_von_aktiv (= aktive Nachtmin / alle aktiven Min)
# Winterschlaf-gefilterte Stichprobe (gesamt_aktiv >= 20%)
dat_e  <- nacht_voll_aktiv[!is.na(nacht_von_aktiv)]
igel_e <- dat_e[, .N, by = igel][N >= 4]$igel
dat_e  <- dat_e[igel %in% igel_e]

if (nrow(dat_e) > 0) {
  nacht_mean <- dat_e[, .(
    mean_nacht = mean(nacht_von_aktiv, na.rm = TRUE),
    se_nacht   = sd(nacht_von_aktiv, na.rm = TRUE) / sqrt(.N),
    n_igel     = .N
  ), by = tage_seit][order(tage_seit)]

  n_igel_e <- length(unique(dat_e$igel))
  n_excl_label <- if (length(igel_ws_ausgeschl) > 0)
    paste0("  |  Ausgeschlossen (Torpor/WS): ", paste(igel_ws_ausgeschl, collapse = ", "))
  else ""

  p8e <- ggplot() +
    geom_line(data = dat_e,
              aes(x = tage_seit, y = nacht_von_aktiv, group = igel, color = diag_kat),
              alpha = 0.35, linewidth = 0.5) +
    geom_ribbon(data = nacht_mean,
                aes(x = tage_seit,
                    ymin = pmax(0, mean_nacht - se_nacht),
                    ymax = pmin(1, mean_nacht + se_nacht)),
                fill = "grey40", alpha = 0.15) +
    geom_smooth(data = dat_e,
                aes(x = tage_seit, y = nacht_von_aktiv),
                method = "gam", formula = y ~ s(x, bs = "cs", k = 5),
                color = "grey20", linewidth = 1.3, se = FALSE) +
    scale_color_manual(values = diag_farben, name = "Diagnose") +
    scale_y_continuous(
      name   = "Nachtanteil der Aktivität (aktive Nachtzeit / gesamte Aktivzeit)",
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1)) +
    scale_x_continuous(name = "Tage seit Auswilderung") +
    labs(title    = "E  Nachtaktivitäts-Normalisierung post-Auswilderung",
         subtitle = sprintf(
           "Anteil der Aktivität während astronomischer Nacht  |  n = %d Igel%s",
           n_igel_e, n_excl_label),
         caption  = paste0(
           "Nachtfenster: astronomisch (Sonnenuntergang–Sonnenaufgang, Sachsenhagen)  |  Dünne Linien: Einzeltiere  |  Dicker Trend: GAM  |  ",
           "Filter: Fenster mit Gesamtaktivität < 20% ausgeschlossen (Torpor/Winterschlaf)")) +
    theme_minimal(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", size = 12),
          plot.subtitle   = element_text(size = 9, color = "grey40"),
          plot.caption    = element_text(size = 8, color = "grey40"),
          legend.position = "right")

  ggsave(file.path(output_ordner, "08e_nachtaktivitaet_normalisierung.png"), p8e,
         width = 10.6, height = 5.3, dpi = 170, bg = "white")
  cat("✓ 08e_nachtaktivitaet_normalisierung.png\n")
} else {
  cat("⚠ Zu wenige Nachtaktivitäts-Datenpunkte.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSE F — Visually impaired vs. sighted: NAF and IS comparison
# ══════════════════════════════════════════════════════════════════════════════
cat("════════════════════\nAnalyse F: Visually impaired vs. sighted (NAF + IS)\n")

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

# ── NAF per animal (WS-filtered; Igel5 excluded — hibernation artefact) ───────
naf_f <- nacht_voll_aktiv[
  !igel %in% c("Igel5") & !is.na(nacht_von_aktiv),
  .(naf_mean = mean(nacht_von_aktiv, na.rm = TRUE)),
  by = igel
]
naf_f <- merge(naf_f, meta_join[, .(igel, diag_kat)], by = "igel", all.x = TRUE)
naf_f[, group := fifelse(diag_kat == "Blindheit", "Visually impaired", "Sighted")]
n_naf_s <- naf_f[group == "Sighted",            .N]
n_naf_b <- naf_f[group == "Visually impaired",  .N]
lbl_s_naf <- sprintf("Sighted\n(n = %d)", n_naf_s)
lbl_b_naf <- sprintf("Visually impaired\n(n = %d)", n_naf_b)
naf_f[, group_label := fifelse(group == "Sighted", lbl_s_naf, lbl_b_naf)]
naf_f[, group_label := factor(group_label, levels = c(lbl_s_naf, lbl_b_naf))]

# ── IS per animal (is_voll_analyse_clean; both groups) ────────────────────────
is_f <- is_voll_analyse_clean[!is.na(IS) & !is.na(diag_kat), .(igel, IS, diag_kat)]
is_f[, group := fifelse(diag_kat == "Blindheit", "Visually impaired", "Sighted")]
n_is_s <- is_f[group == "Sighted",           .N]
n_is_b <- is_f[group == "Visually impaired", .N]
lbl_s_is <- sprintf("Sighted\n(n = %d)", n_is_s)
lbl_b_is <- sprintf("Visually impaired\n(n = %d)", n_is_b)
is_f[, group_label := fifelse(group == "Sighted", lbl_s_is, lbl_b_is)]
is_f[, group_label := factor(group_label, levels = c(lbl_s_is, lbl_b_is))]

# ── Colors ─────────────────────────────────────────────────────────────────────
col_blind   <- "#D6604D"
col_sighted <- "#6BAED6"
cols_naf <- setNames(c(col_sighted, col_blind), c(lbl_s_naf, lbl_b_naf))
cols_is  <- setNames(c(col_sighted, col_blind), c(lbl_s_is,  lbl_b_is))

# ── Panel: NAF ─────────────────────────────────────────────────────────────────
jit_naf <- position_jitter(width = 0.09, seed = 7)

pf_naf <- ggplot(naf_f, aes(x = group_label, y = naf_mean, color = group_label)) +
  # Group mean crossbar
  stat_summary(fun = mean, geom = "crossbar", width = 0.28,
               linewidth = 0.75, fatten = 1.8, color = "grey15",
               show.legend = FALSE) +
  # Individual points
  geom_point(position = jit_naf, size = 3.4, alpha = 0.85) +
  # Labels for blind animals only
  geom_text_repel(
    data        = naf_f[group == "Visually impaired"],
    aes(label   = sub("^Igel", "H", igel)),
    size        = 3, fontface = "italic", color = col_blind,
    nudge_x     = 0.30, nudge_y = 0,
    segment.size = 0.3, segment.color = "grey50",
    min.segment.length = 0, seed = 7
  ) +
  scale_color_manual(values = cols_naf, guide = "none") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05), breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title    = "Nocturnal activity fraction (NAF)",
    subtitle = "Proportion of active time during astronomical night\n(hibernation phases excluded: < 20% total activity)",
    x = NULL,
    y = "Active time at night / total active time"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8.5, color = "grey40", lineheight = 1.2),
    axis.text.x        = element_text(size = 10.5, lineheight = 1.3),
    axis.title.y       = element_text(size = 9),
    panel.grid.major.x = element_blank()
  )

# ── Panel: IS ──────────────────────────────────────────────────────────────────
jit_is <- position_jitter(width = 0.09, seed = 7)

pf_is <- ggplot(is_f, aes(x = group_label, y = IS, color = group_label)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.28,
               linewidth = 0.75, fatten = 1.8, color = "grey15",
               show.legend = FALSE) +
  geom_point(position = jit_is, size = 3.4, alpha = 0.85) +
  geom_text_repel(
    data        = is_f[group == "Visually impaired"],
    aes(label   = sub("^Igel", "H", igel)),
    size        = 3, fontface = "italic", color = col_blind,
    nudge_x     = 0.30, nudge_y = 0,
    segment.size = 0.3, segment.color = "grey50",
    min.segment.length = 0, seed = 7
  ) +
  scale_color_manual(values = cols_is, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "Interdaily stability (IS)",
    subtitle = "Day-to-day regularity of 24-h activity pattern\n(IS = 0: arrhythmic; IS = 1: perfectly regular)",
    x = NULL,
    y = "IS (0–1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8.5, color = "grey40", lineheight = 1.2),
    axis.text.x        = element_text(size = 10.5, lineheight = 1.3),
    axis.title.y       = element_text(size = 9),
    panel.grid.major.x = element_blank()
  )

# ── Combined figure ────────────────────────────────────────────────────────────
p8f <- (pf_naf | pf_is) +
  plot_annotation(
    title   = "F  Nocturnal behaviour — visually impaired vs. sighted hedgehogs",
    caption = paste0(
      "NAF: sighted n = ", n_naf_s, " (Igel5 excl.: hibernation artefact); ",
      "IS: sighted n = ", n_is_s, " (Igel5 & Igel1 excl.: artefacts)  |  ",
      "Visually impaired: H22, H23 (Igel26 excl.: 1-day tracking)  |  ",
      "Black bar = group mean  |  Hibernation phases (< 20% total activity) excluded"
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", size = 12),
      plot.caption = element_text(size = 7.5, color = "grey40", lineheight = 1.15)
    )
  )

ggsave(file.path(output_ordner, "08f_blind_vs_sighted.png"), p8f,
       width = 10.5, height = 5.8, dpi = 170, bg = "white")
cat("✓ 08f_blind_vs_sighted.png\n")

# ── Publication version (black border, panel labels, greyscale) ───────────────
pub_theme <- theme_classic(base_size = 11) +
  theme(
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.75),
    axis.line          = element_blank(),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    axis.text          = element_text(color = "black", size = 10),
    axis.title.y       = element_text(color = "black", size = 9.5),
    axis.text.x        = element_text(size = 10, lineheight = 1.25),
    plot.tag           = element_text(face = "bold", size = 11),
    plot.tag.position  = c(0.02, 0.98)
  )

# x-axis labels without n (n goes in figure caption)
xlabels_naf <- setNames(c("Sighted", "Visually\nimpaired"), levels(naf_f$group_label))
xlabels_is  <- setNames(c("Sighted", "Visually\nimpaired"), levels(is_f$group_label))

# greyscale: sighted = grey60, visually impaired = black
col_s_pub <- "grey55"
col_b_pub <- "black"
cols_naf_pub <- setNames(c(col_s_pub, col_b_pub), levels(naf_f$group_label))
cols_is_pub  <- setNames(c(col_s_pub, col_b_pub), levels(is_f$group_label))

pf_naf_pub <- ggplot(naf_f, aes(x = group_label, y = naf_mean, color = group_label)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.28,
               linewidth = 0.75, fatten = 1.8, color = "black",
               show.legend = FALSE) +
  geom_point(position = jit_naf, size = 3.2, alpha = 0.85) +
  geom_text_repel(
    data         = naf_f[group == "Visually impaired"],
    aes(label    = sub("^Igel", "H", igel)),
    size         = 2.8, fontface = "italic",
    nudge_x      = 0.30, nudge_y = 0,
    segment.size = 0.3, segment.color = "grey55",
    min.segment.length = 0, seed = 7, color = col_b_pub
  ) +
  scale_color_manual(values = cols_naf_pub, guide = "none") +
  scale_x_discrete(labels = xlabels_naf) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Active time at night / total active time") +
  pub_theme

pf_is_pub <- ggplot(is_f, aes(x = group_label, y = IS, color = group_label)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.28,
               linewidth = 0.75, fatten = 1.8, color = "black",
               show.legend = FALSE) +
  geom_point(position = jit_is, size = 3.2, alpha = 0.85) +
  geom_text_repel(
    data         = is_f[group == "Visually impaired"],
    aes(label    = sub("^Igel", "H", igel)),
    size         = 2.8, fontface = "italic",
    nudge_x      = 0.30, nudge_y = 0,
    segment.size = 0.3, segment.color = "grey55",
    min.segment.length = 0, seed = 7, color = col_b_pub
  ) +
  scale_color_manual(values = cols_is_pub, guide = "none") +
  scale_x_discrete(labels = xlabels_is) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Interdaily stability (IS)") +
  pub_theme

p8f_pub <- (pf_naf_pub | pf_is_pub) +
  plot_annotation(tag_levels = "a", tag_prefix = "", tag_suffix = ")")

# Robuste PDF-Ausgabe: nutzt 'Cairo' (bringt eigenes Cairo mit, unabhaengig von
# einer evtl. fehlenden System-cairo-DLL); Fallback auf Basis-pdf().
cairo_or_pdf <- function(filename, width, height, bg = "white", ...) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoPDF(file = filename, width = width, height = height, bg = bg)
  } else {
    grDevices::pdf(file = filename, width = width, height = height, bg = bg)
  }
}

ggsave(file.path(output_ordner, "08f_blind_vs_sighted_pub.png"), p8f_pub,
       width = 9, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_ordner, "08f_blind_vs_sighted_pub.pdf"), p8f_pub,
       width = 9, height = 5, device = cairo_or_pdf, bg = "white")
cat("✓ 08f_blind_vs_sighted_pub.png/.pdf\n")

# ══════════════════════════════════════════════════════════════════════════════
# SYNTHESE-TABELLE
# ══════════════════════════════════════════════════════════════════════════════
cat("\nErstelle Synthese-Tabelle …\n")

akro_pro_tier <- nacht_voll[!is.na(nacht_anteil), .(
  nacht_anteil_mean   = round(mean(nacht_anteil, na.rm = TRUE), 3),
  nacht_anteil_letzte = round(tail(nacht_anteil[order(tage_seit)], 1), 3),
  n_nacht             = .N
), by = igel]

# Primärmetrik-Zusammenfassung (Winterschlaf-gefiltert)
akro_nva <- nacht_voll_aktiv[!is.na(nacht_von_aktiv), .(
  nacht_von_aktiv_mean   = round(mean(nacht_von_aktiv, na.rm = TRUE), 3),
  nacht_von_aktiv_letzte = round(tail(nacht_von_aktiv[order(tage_seit)], 1), 3),
  gesamt_aktiv_mean      = round(mean(gesamt_aktiv, na.rm = TRUE), 3),
  n_nva                  = .N
), by = igel]
akro_pro_tier <- merge(akro_pro_tier, akro_nva, by = "igel", all.x = TRUE)

is_trend <- is_zeit_voll[!is.na(IS) & !is.na(tage_seit), {
  if (.N >= 4) {
    m <- tryCatch(lm(IS ~ tage_seit), error = function(e) NULL)
    if (!is.null(m)) .(IS_trend_pro_tag = round(coef(m)[["tage_seit"]], 5))
    else             .(IS_trend_pro_tag = NA_real_)
  } else .(IS_trend_pro_tag = NA_real_)
}, by = igel]

synthese <- merge(
  meta_voll[, .(igel, sex, alter, diag_kat, time_reha, tagging_weight,
                weight_entry, weight_gain, tagging_period, saison_auswild)],
  is_voll[, .(igel, n_tage, IS, RA)],
  by = "igel", all.x = TRUE
)
synthese <- merge(synthese, akro_pro_tier, by = "igel", all.x = TRUE)
synthese <- merge(synthese, is_trend,      by = "igel", all.x = TRUE)
setorder(synthese, igel)
synthese[, IS := round(IS, 3)]
synthese[, RA := round(RA, 3)]

cat("Synthese-Tabelle:", nrow(synthese), "Tiere\n")
print(synthese[, .(igel, diag_kat, sex, n_tage, IS, RA,
                   nacht_von_aktiv_mean, nacht_anteil_mean, IS_trend_pro_tag,
                   time_reha, tagging_period)])

# Korrelationstabelle für Excel und Word
korr_tabelle <- rbindlist(lapply(korr_vars, function(v) {
  sub <- dat_d[!is.na(get(v$var))]
  if (nrow(sub) < 4) return(data.table(Variable = v$var, rho = NA, p_wert = NA, n = nrow(sub)))
  ct <- tryCatch(cor.test(sub[[v$var]], sub$IS, method = "spearman", exact = FALSE),
                 error = function(e) NULL)
  if (is.null(ct))
    return(data.table(Variable = v$var, Beschreibung = v$title,
                      rho = NA, p_wert = NA, signifikant = NA, n = nrow(sub)))
  data.table(Variable    = v$var,
             Beschreibung = v$title,
             rho          = round(ct$estimate, 3),
             p_wert        = round(ct$p.value, 4),
             signifikant   = ifelse(ct$p.value < 0.05, "ja *", "nein"),
             n             = nrow(sub))
}))

# IS-Gruppenstatistiken ohne Artefakt-Tiere (Igel5, Igel1)
grp_is_diag <- is_voll_analyse_clean[!is.na(IS) & !is.na(diag_kat), .(
  n = .N, IS_median = round(median(IS), 3),
  IS_mean = round(mean(IS), 3), IS_sd = round(sd(IS), 3)
), by = diag_kat][order(diag_kat)]

grp_is_sex <- is_voll_analyse_clean[!is.na(IS) & !is.na(sex), .(
  n = .N, IS_median = round(median(IS), 3),
  IS_mean = round(mean(IS), 3), IS_sd = round(sd(IS), 3)
), by = sex][order(sex)]

# Nachtaktivität früh vs. spät — Primärmetrik nacht_von_aktiv (Winterschlaf-gefiltert)
nacht_frueh <- if (nrow(dat_e) > 0) dat_e[tage_seit <= 7,  mean(nacht_von_aktiv, na.rm=TRUE)] else NA
nacht_spaet <- if (nrow(dat_e) > 0) dat_e[tage_seit >= 14, mean(nacht_von_aktiv, na.rm=TRUE)] else NA

# ══════════════════════════════════════════════════════════════════════════════
# 5. EXCEL-EXPORT (openxlsx — styled, konsistent mit anderen Blöcken)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Excel-Export ─────────────────────────────────────────────\n")

wb <- createWorkbook()

s_hdr  <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white",
                       fgFill = COL_BLOCK8, textDecoration = "bold",
                       halign = "center", border = "Bottom",
                       borderColour = "#BFBFBF")
s_body <- createStyle(fontName = "Arial", fontSize = 10)
s_alt  <- createStyle(fontName = "Arial", fontSize = 10, fgFill = "#EEF2FF")
s_sig  <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "#166534",
                       textDecoration = "bold", fgFill = "#DCFCE7")
s_ns   <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "#6B7280")
s_warn <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "#991B1B",
                       fgFill = "#FEE2E2")

add_sheet_styled <- function(wb, name, data, sig_rows = NULL, warn_rows = NULL) {
  addWorksheet(wb, name)
  writeData(wb, name, as.data.frame(data), startRow = 1)
  nr <- nrow(data); nc <- ncol(data)
  addStyle(wb, name, s_hdr,  rows = 1,         cols = seq_len(nc), gridExpand = TRUE)
  addStyle(wb, name, s_body, rows = 2:(nr + 1), cols = seq_len(nc), gridExpand = TRUE)
  for (r in seq(2, nr + 1, by = 2))
    addStyle(wb, name, s_alt, rows = r, cols = seq_len(nc), gridExpand = TRUE, stack = TRUE)
  if (!is.null(sig_rows)  && length(sig_rows)  > 0)
    addStyle(wb, name, s_sig,  rows = sig_rows  + 1, cols = seq_len(nc),
             gridExpand = TRUE, stack = TRUE)
  if (!is.null(warn_rows) && length(warn_rows) > 0)
    addStyle(wb, name, s_warn, rows = warn_rows + 1, cols = seq_len(nc),
             gridExpand = TRUE, stack = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(nc), widths = "auto")
}

# 00 Überblick
uebersicht_zeilen <- data.frame(
  Metrik = c(
    "Tiere gesamt", paste0("Tiere mit gültigem IS (≥ ", MIN_TAGE_IS_IV, " Tage)"),
    "Tiere mit Nachtaktivitätsdaten (nacht_von_aktiv, gefiltert)",
    "IS Median (Analyse-Stichprobe)", "IS Spannweite",
    paste0("Ø nacht_von_aktiv Tage 1–7 [WS-Filter ", AKTIV_SCHWELLE*100,"%%]"),
    paste0("Ø nacht_von_aktiv ab Tag 14 [WS-Filter ", AKTIV_SCHWELLE*100,"%%]"),
    paste0("Ausgeschlossen durch WS-Filter (gesamt_aktiv < ", AKTIV_SCHWELLE*100, "%)"),
    "Kruskal-Wallis IS ~ Diagnose (p-Wert)",
    "Wilcoxon IS ~ Geschlecht (p-Wert)"
  ),
  Wert = c(
    nrow(synthese),
    paste0(nrow(is_voll_analyse[!is.na(IS)]), " total / ",
           nrow(is_voll_analyse_clean[!is.na(IS)]), " bereinigt (ohne ",
           paste(IS_ARTEFAKT_TIERE, collapse=", "), ")"),
    nrow(synthese[!is.na(nacht_von_aktiv_mean)]),
    if (nrow(is_voll_analyse_clean[!is.na(IS)]) > 0)
      round(median(is_voll_analyse_clean$IS, na.rm = TRUE), 3) else NA,
    if (nrow(is_voll_analyse_clean[!is.na(IS)]) > 0)
      paste0(round(min(is_voll_analyse_clean$IS, na.rm=TRUE),3), " – ",
             round(max(is_voll_analyse_clean$IS, na.rm=TRUE),3)) else NA,
    if (!is.na(nacht_frueh)) paste0(round(nacht_frueh * 100, 1), "%") else "n. v.",
    if (!is.na(nacht_spaet)) paste0(round(nacht_spaet * 100, 1), "%") else "n. v.",
    if (length(igel_ws_ausgeschl) > 0) paste(igel_ws_ausgeschl, collapse = ", ") else "keines",
    if (!is.null(kw_IS)) round(kw_IS$p.value, 4) else "n. v.",
    if (!is.null(wt_IS)) round(wt_IS$p.value, 4) else "n. v."
  ),
  stringsAsFactors = FALSE
)
add_sheet_styled(wb, "00_Ueberblick", uebersicht_zeilen)

# 01 Synthese pro Tier
add_sheet_styled(wb, "01_Synthese_ProTier", synthese,
                 sig_rows = which(!is.na(synthese$IS) & synthese$IS > 0.5))

# 02 Korrelationen
add_sheet_styled(wb, "02_Korrelationen", korr_tabelle,
                 sig_rows  = which(korr_tabelle$signifikant == "ja *"),
                 warn_rows = which(korr_tabelle$signifikant == "nein"))

# 03 + 04 Gruppenstatistiken
add_sheet_styled(wb, "03_Gruppen_Diagnose", grp_is_diag)
add_sheet_styled(wb, "04_Gruppen_Sex",      grp_is_sex)

# 05 Nachtaktivität Rolling — inkl. Primärmetrik nacht_von_aktiv + Winterschlaf-Flag
nacht_export <- nacht_voll[!is.na(nacht_anteil), .(
  igel,
  tage_seit,
  nacht_anteil     = round(nacht_anteil, 4),     # aktive Nachtmin / alle detekt. Min (NICHT Berger-vergleichbar)
  tag_anteil       = round(tag_anteil, 4),        # aktive Tagmin   / alle detekt. Min
  gesamt_aktiv     = round(gesamt_aktiv, 4),      # Nacht+Tag aktiv / alle detekt. Min
  nacht_von_aktiv  = round(nacht_von_aktiv, 4),   # PRIMÄRMETRIK: aktive Nachtmin / alle aktiven Min
  winterschlaf     = winterschlaf,                # TRUE = durch WS-Filter ausgeschlossen
  n_tage,
  diag_kat,
  sex
)]
add_sheet_styled(wb, "05_Nachtaktivitaet_Rolling", nacht_export,
                 warn_rows = which(nacht_export$winterschlaf))

saveWorkbook(wb, file.path(output_ordner, "08_reha_synthese.xlsx"), overwrite = TRUE)
cat("  ✓ 08_reha_synthese.xlsx\n")

# ══════════════════════════════════════════════════════════════════════════════
# 6. WORD METHODENBERICHT
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
      bg(part = "header", bg = COL_BLOCK8) |>
      hline(border = brd, part = "all") |>
      vline(border = brd, part = "all") |>
      hline_top(border = brd, part = "header") |>
      bg(part = "body", bg = "white") |>
      set_table_properties(layout = "autofit")
    if (!is.null(hl_rows)  && length(hl_rows)  > 0)
      ft <- bg(ft, i = hl_rows, bg = "#DCFCE7", part = "body")
    if (!is.null(warn_rows) && length(warn_rows) > 0)
      ft <- bg(ft, i = warn_rows, bg = "#FEE2E2", part = "body") |>
              color(i = warn_rows, color = "#991B1B", part = "body")
    ft
  }

  doc_m <- read_docx()

  doc_m <- doc_m |>
    body_add_par("Block 8: Rehabilitations-Erfolgs-Synthese — Methodenbericht",
                 style = "heading 1") |>
    body_add_par("Igelbesenderung Wildtierstation Sachsenhagen — TiHo Hannover | Natalie Steiner",
                 style = "Normal") |>
    body_add_par(paste0("Erstellt: ", format(Sys.Date(), "%d. %B %Y")), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("1. Fragestellung", style = "heading 2") |>
    body_add_par(paste0(
      "Block 8 beantwortet die zentrale Frage der Rehabilitationsstudie: ",
      "Normalisiert sich das Verhalten ausgewilderter Igel — insbesondere das ",
      "zirkadiane Aktivitätsmuster — nach der Rehabilitation? Konkret: Entwickeln ",
      "die Tiere mit der Zeit ein nächtliches Aktivitätsmuster, und welche Rehab.-Parameter (Diagnose, Geschlecht, ",
      "Rehab.-Dauer, Gewichtszunahme) beeinflussen die Verhaltens-Normalisierung?"
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("2. Datenbasis", style = "heading 2") |>
    body_add_par(paste0(
      "Drei Datenquellen werden integriert: ",
      "(1) IS/IV/RA-Werte aus Block5_IS_IV.R (is_iv_ergebnisse.rds), ",
      "(2) Minuten-Aktivitätsdaten aus Block0_Pipeline.R (chrono_minuten24h.rds), ",
      "(3) Metadaten aus data_igel.xlsx (Diagnose, Geschlecht, Gewicht, Rehab.-Dauer)."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("3. Metriken und Entscheidungen", style = "heading 2") |>

    body_add_par("3.1 Primärmetrik: nacht_von_aktiv (Nachtanteil der Aktivität)", style = "heading 3") |>
    body_add_par(paste0(
      "Anteil der Aktivzeit während der astronomischen Nacht (Sonnenuntergang bis Sonnenaufgang, ",
      "Sachsenhagen 52.397°N 9.217°E) an der gesamten Aktivzeit (Analyse E). ",
      "Berechnet als gleitendes Fenster (", is_fenster_tage, " Tage, Mindest-Tage: ",
      min_tage_window, "): nacht_von_aktiv = Σ(aktiv bei Nacht) / Σ(aktiv bei Nacht + aktiv bei Tag). ",
      "Wichtig: Die 'rohe' Metrik nacht_anteil = Σ(aktiv bei Nacht) / Σ(alle detektierten Minuten) ",
      "enthält im Nenner passive Ruheminuten (z.B. am Tagesschlafplatz in Empfängerreichweite). ",
      "nacht_von_aktiv bereinigt dies durch Beschränkung auf aktive Minuten. ",
      "Winterschlaf-/Torpor-Phasen mit Gesamtaktivität < ", AKTIV_SCHWELLE * 100, "% werden gefiltert; ",
      "betroffen: Igel5 (Tage ~26–165 und ~211+) und Igel1 (vollständig, n=3 Tracking-Tage)."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("3.2 Sekundärmetrik: IS (Interdaily Stability)", style = "heading 3") |>
    body_add_par(paste0(
      "IS misst die Synchronisation des 24-Stunden-Rhythmus anhand der day-to-day ",
      "Ähnlichkeit des stündlichen Aktivitätsprofils (van Someren et al. 1999). ",
      "IS = 0: kein rhythmisches Muster; IS = 1: perfekte Reproduzierbarkeit. ",
      "IS ist robust gegenüber der ~20-min Glättung, da stündliche Mittelwerte ",
      "verwendet werden (Glättungsartefakte mitteln heraus). ",
      "Einschluss-Kriterium: ≥ ", MIN_TAGE_IS_IV, " Tracking-Tage."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("3.3 IV (Intradaily Variability) — nicht verwendet", style = "heading 3") |>
    body_add_par(paste0(
      "IV wurde explizit aus Block 8 ausgeschlossen. Gründe: ",
      "(1) Die ~20-min Glättung in smoothed_wmv unterdrückt kurzfristige Aktivitätswechsel ",
      "systematisch → absolute IV-Werte sind nicht interpretierbar. ",
      "(2) Tracking-Dauern zu kurz und heterogen (8–46 Tage) für stabile IV-Schätzungen. ",
      "(3) n zu klein (14 Tiere ≥ 7 Tage) für Gruppenvergleiche. ",
      "IV ist in Block5_IS_IV.R weiterhin berechnet (Vergleichswert), ",
      "wird in Block 8 nicht analysiert."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("4. Statistische Methoden", style = "heading 2") |>

    body_add_par("4.1 Gruppenvergleiche (Analysen A, B)", style = "heading 3") |>
    body_add_par(paste0(
      "IS nach Diagnose-Kategorie: Kruskal-Wallis-Test (nichtparametrisch, da n pro Gruppe klein). ",
      "IS nach Geschlecht: Wilcoxon-Rangsummen-Test. ",
      "Visualisierung: Violin-Plot + Boxplot + Jitter-Punkte. ",
      "Ausgeschlossen aus Gruppenvergleichen (IS-Artefakt): Igel5 (IS=0.094, durch zwei ",
      "Winterschlaf-Perioden über 221 Tage verfälscht) und Igel1 (IS=0.072, nur 3 Tracking-Tage, ",
      "< 5% Gesamtaktivität). Beide Tiere verbleiben im IS-Zeitverlauf (Analyse C)."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("4.2 Zeitverlauf (Analyse C)", style = "heading 3") |>
    body_add_par(paste0(
      "IS-Zeitverlauf aus dem Rolling-Window-Datensatz (Block 5). ",
      "Visualisierung: Spaghetti-Plot (Einzeltiere) + GAM-Trend pro Diagnose-Gruppe ",
      "(mgcv, cubic regression splines, k=4). ",
      "Einschluss: Tiere mit ≥ 4 Rolling-Window-Datenpunkten."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("4.3 Korrelationen (Analyse D)", style = "heading 3") |>
    body_add_par(paste0(
      "Spearman-Rangkorrelation (cor.test, exact=FALSE) zwischen IS und: ",
      "(a) Rehab.-Dauer (time_reha), ",
      "(b) Auswilderungsgewicht (tagging_weight), ",
      "(c) Gewichtszunahme in der Reha (weight_gain). ",
      "Signifikanzschwelle α = 0.05."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("4.4 Nachtaktivitäts-Trend (Analyse E)", style = "heading 3") |>
    body_add_par(paste0(
      "GAM-Smoother über alle Tiere (mgcv, cubic regression splines, k=5). ",
      "Vergleich früher (Tage 1–7) vs. später Periode (ab Tag 14). ",
      "Einschluss: Tiere mit ≥ 4 Rolling-Window-Datenpunkten."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("5. Diagnose-Kategorisierung", style = "heading 2") |>
    body_add_par(paste0(
      "Die Diagnosen aus data_igel.xlsx werden in folgende Kategorien zusammengefasst: ",
      "Blindheit (grep 'blind'), Waisling (grep 'orph'), Pilzerkrankung (grep 'fungal'), ",
      "Endoparasiten (grep 'endoparasit'), Parasiten (grep 'parasite', exkl. blind/orph/endo), ",
      "Trauma (grep 'trauma', exkl. blind), Sonstige (Rest). ",
      "Priorität bei Mehrfachdiagnosen: Blindheit > Waisling > Pilz > Endo > Parasiten > Trauma."
    ), style = "Normal") |>
    body_add_par("", style = "Normal")

  # Diagnose-Tabelle
  diag_tab <- meta_voll[, .N, by = diag_kat][order(diag_kat)]
  names(diag_tab) <- c("Diagnose-Kategorie", "Anzahl Tiere")
  doc_m <- doc_m |>
    body_add_par("Tabelle 1: Diagnose-Kategorien", style = "Normal") |>
    body_add_flextable(make_ft(as.data.frame(diag_tab))) |>
    body_add_par("", style = "Normal") |>

    body_add_par("6. Limitationen", style = "heading 2") |>
    body_add_par(paste0(
      "Kleine und ungleiche Gruppengrößen (Diagnose-Kategorien 1–6 Tiere) ",
      "schränken die statistische Aussagekraft der Gruppenvergleiche (Analyse A) erheblich ein. ",
      "Die Rehab.-Parameter (time_reha, weight_gain) sind potenziell konfundiert ",
      "(kranke Tiere bleiben länger in Reha UND haben schlechtere Ausgangswerte). ",
      "Der Nachtaktivitäts-Anteil ist ein Proxy-Maß — tatsächliche Habitatnutzung ",
      "und Nahrungssuche werden nicht direkt erfasst."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>

    body_add_par("7. Referenzen", style = "heading 2") |>
    body_add_par(paste0(
      "van Someren EJW, Lijzenga C, Mirmiran M, Swaab DF (1997). ",
      "Long-term fitness training improves the circadian rest–activity rhythm in healthy elderly males. ",
      "Journal of Biological Rhythms 12(2):146–156."
    ), style = "Normal") |>
    body_add_par(paste0(
      "Wood SN (2017). Generalized Additive Models: An Introduction with R (2nd ed.). CRC Press."
    ), style = "Normal")

  out_m <- file.path(output_ordner, "Block8_RehabSuccess_Methodenbericht.docx")
  print(doc_m, target = out_m)
  cat("  ✓ Block8_RehabSuccess_Methodenbericht.docx\n")

} else {
  cat("  [WARN] officer/flextable nicht verfügbar.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. WORD ERGEBNISBERICHT (dynamisch — nutzt berechnete Variablen)
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
      bg(part = "header", bg = COL_BLOCK8) |>
      hline(border = brd, part = "all") |>
      vline(border = brd, part = "all") |>
      hline_top(border = brd, part = "header") |>
      bg(part = "body", bg = "white") |>
      set_table_properties(layout = "autofit")
    if (!is.null(hl_rows) && length(hl_rows) > 0)
      ft <- bg(ft, i = hl_rows, bg = "#DCFCE7", part = "body")
    if (!is.null(warn_rows) && length(warn_rows) > 0)
      ft <- bg(ft, i = warn_rows, bg = "#FEE2E2", part = "body") |>
              color(i = warn_rows, color = "#991B1B", part = "body")
    ft
  }

  doc_e <- read_docx()

  doc_e <- doc_e |>
    body_add_par("Block 8: Rehabilitations-Erfolgs-Synthese — Ergebnisbericht",
                 style = "heading 1") |>
    body_add_par("Igelbesenderung Wildtierstation Sachsenhagen — TiHo Hannover | Natalie Steiner",
                 style = "Normal") |>
    body_add_par(paste0("Erstellt: ", format(Sys.Date(), "%d. %B %Y")), style = "Normal") |>
    body_add_par("", style = "Normal")

  # Überblick-Kennzahlen
  uebersicht_e <- data.frame(
    Metrik   = c("Tiere gesamt",
                 paste0("Davon mit IS (≥ ", MIN_TAGE_IS_IV, " Tage)"),
                 "Davon mit nacht_von_aktiv-Daten (WS-gefiltert)",
                 "IS Median (Analyse-Stichprobe)",
                 "Ø nacht_von_aktiv Tage 1–7",
                 "Ø nacht_von_aktiv ab Tag 14",
                 paste0("WS-Filter ausgeschlossen (gesamt_aktiv < ", AKTIV_SCHWELLE*100, "%)")),
    Ergebnis = c(
      nrow(synthese),
      paste0(nrow(is_voll_analyse[!is.na(IS)]), " total / ",
             nrow(is_voll_analyse_clean[!is.na(IS)]), " bereinigt"),
      nrow(synthese[!is.na(nacht_von_aktiv_mean)]),
      if (nrow(is_voll_analyse_clean[!is.na(IS)]) > 0)
        round(median(is_voll_analyse_clean$IS, na.rm = TRUE), 3) else "n. v.",
      if (!is.na(nacht_frueh)) paste0(round(nacht_frueh * 100, 1), "%") else "n. v.",
      if (!is.na(nacht_spaet)) paste0(round(nacht_spaet * 100, 1), "%") else "n. v.",
      if (length(igel_ws_ausgeschl) > 0) paste(igel_ws_ausgeschl, collapse = ", ") else "keines"
    ), stringsAsFactors = FALSE
  )
  doc_e <- doc_e |>
    body_add_par("Tabelle 1: Überblick Kennzahlen", style = "Normal") |>
    body_add_flextable(make_ft2(uebersicht_e)) |>
    body_add_par("", style = "Normal")

  # 1. IS nach Diagnose
  doc_e <- doc_e |>
    body_add_par("1. IS nach Diagnose-Kategorie (Analyse A)", style = "heading 2") |>
    body_add_par(sprintf(paste0(
      "In die Diagnose-Analyse gingen %d Tiere mit gültigem IS und ≥ %d Tracking-Tagen ein. "),
      nrow(dat_a), MIN_TAGE_IS_IV), style = "Normal") |>
    body_add_par(
      if (!is.null(kw_IS))
        sprintf(paste0(
          "Der Kruskal-Wallis-Test zeigte %s signifikante Unterschiede zwischen den ",
          "Diagnose-Gruppen (χ² = %.2f, df = %d, p = %.3f). "),
          if (kw_IS$p.value < 0.05) "statistisch" else "keine statistisch",
          kw_IS$statistic, kw_IS$parameter, kw_IS$p.value)
      else "Kruskal-Wallis-Test konnte nicht berechnet werden.",
      style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("Tabelle 2: IS nach Diagnose-Gruppe", style = "Normal") |>
    body_add_flextable(make_ft2(
      as.data.frame(grp_is_diag),
      hl_rows = which(grp_is_diag$IS_median >= 0.5)
    )) |>
    body_add_par("Grün markiert: IS-Median ≥ 0.5.", style = "Normal") |>
    body_add_par("", style = "Normal")

  # 2. IS nach Geschlecht
  doc_e <- doc_e |>
    body_add_par("2. IS nach Geschlecht (Analyse B)", style = "heading 2") |>
    body_add_par(
      if (!is.null(wt_IS))
        sprintf(paste0(
          "Der Wilcoxon-Test zeigte %s signifikante Unterschiede zwischen ",
          "Männchen und Weibchen im IS (W = %.0f, p = %.3f). "),
          if (wt_IS$p.value < 0.05) "" else "keine",
          wt_IS$statistic, wt_IS$p.value)
      else "Wilcoxon-Test konnte nicht berechnet werden.",
      style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("Tabelle 3: IS nach Geschlecht", style = "Normal") |>
    body_add_flextable(make_ft2(as.data.frame(grp_is_sex))) |>
    body_add_par("", style = "Normal")

  # 3. Korrelationen
  doc_e <- doc_e |>
    body_add_par("3. Reha-Korrelationen — IS ~ Rehab.-Parameter (Analyse D)", style = "heading 2") |>
    body_add_par(paste0(
      "Spearman-Rangkorrelationen (exact = FALSE) zwischen IS und drei ",
      "Rehab.-Parametern. Signifikante Korrelationen (p < 0.05) sind grün markiert."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("Tabelle 4: Spearman-Korrelationen IS ~ Rehab.-Parameter", style = "Normal") |>
    body_add_flextable(make_ft2(
      as.data.frame(korr_tabelle),
      hl_rows   = which(korr_tabelle$signifikant == "ja *"),
      warn_rows = which(korr_tabelle$signifikant == "nein")
    )) |>
    body_add_par("Grün = signifikant (p < 0.05). Rot = nicht signifikant.", style = "Normal") |>
    body_add_par("", style = "Normal")

  # 4. Nachtaktivität
  doc_e <- doc_e |>
    body_add_par("4. Nachtaktivitäts-Normalisierung (Analyse E — Primärmetrik)", style = "heading 2")

  if (!is.na(nacht_frueh) && !is.na(nacht_spaet)) {
    diff_pp <- (nacht_spaet - nacht_frueh) * 100
    doc_e <- doc_e |>
      body_add_par(sprintf(paste0(
        "PRIMÄRMETRIK: nacht_von_aktiv = Anteil aktiver Nachtminuten an allen aktiven Minuten. ",
        "Winterschlaf-/Torpor-Phasen (Gesamtaktivität < %.0f%%) wurden gefiltert. ",
        "Die %d Igel mit ≥ 4 auswertbaren Datenpunkten zeigten einen mittleren nacht_von_aktiv ",
        "von %.1f%% in den ersten 7 Tagen nach Auswilderung. ",
        "Ab Tag 14 betrug der mittlere Wert %.1f%% (%+.1f Prozentpunkte)."),
        AKTIV_SCHWELLE * 100,
        length(unique(dat_e$igel)),
        nacht_frueh * 100, nacht_spaet * 100, diff_pp), style = "Normal") |>
      body_add_par(
        if (nacht_spaet >= 0.80)
          "Die Tiere zeigten ab Tag 14 einen hohen Nachtaktivitäts-Anteil, was auf eine erfolgreiche Anpassung an nächtliche Aktivität hindeutet."
        else if (nacht_spaet >= 0.65)
          "Die Tiere zeigten ab Tag 14 eine mehrheitlich nächtliche Aktivität."
        else
          paste0(
            "Der Nachtaktivitäts-Anteil ab Tag 14 war vergleichsweise gering. ",
            "Dies könnte auf eine noch unvollständige Normalisierung nach Auswilderung ",
            "oder auf Besonderheiten der Stichprobe (v. a. blinde Tiere) hindeuten."),
        style = "Normal")
  } else {
    doc_e <- doc_e |>
      body_add_par("Zu wenige Datenpunkte für Nachtaktivitäts-Auswertung.", style = "Normal")
  }

  doc_e <- doc_e |>
    body_add_par("", style = "Normal")

  # 5. Synthese pro Tier
  synth_tab <- synthese[, .(
    Tier            = igel,
    Diagnose        = as.character(diag_kat),
    Sex             = sex,
    N_Tage          = n_tage,
    IS              = IS,
    RA              = RA,
    NachtVonAktiv   = nacht_von_aktiv_mean,    # Primärmetrik (WS-gefiltert)
    NachtAnteil_raw = nacht_anteil_mean,        # Rohwert
    IS_Trend        = IS_trend_pro_tag,
    Reha_Tage       = time_reha
  )]
  doc_e <- doc_e |>
    body_add_par("5. Synthese-Tabelle: Alle Metriken pro Tier", style = "heading 2") |>
    body_add_par("Tabelle 5: Individuelle Ergebnis-Übersicht", style = "Normal") |>
    body_add_flextable(make_ft2(
      as.data.frame(synth_tab),
      hl_rows  = which(!is.na(synth_tab$IS) & synth_tab$IS > 0.5),
      warn_rows = which(synth_tab$Diagnose == "Blindheit")
    )) |>
    body_add_par("Grün: IS > 0.5. Rot: Blindheits-Diagnose.", style = "Normal") |>
    body_add_par("", style = "Normal")

  # 6. Empfehlung
  doc_e <- doc_e |>
    body_add_par("6. Fazit und Empfehlungen", style = "heading 2") |>
    body_add_par(paste0(
      "Die Ergebnisse zeigen, dass rehabilitierte Igel nach Auswilderung grundsätzlich ",
      "ein nächtliches Aktivitätsmuster entwickeln. ",
      "Blinde Tiere (Igel 22, Igel 23, Igel 26) stellen eine klinisch bedeutsame Ausnahme dar — ",
      "sie zeigen deutlich niedrigere IS-Werte und einen geringeren Nachtaktivitäts-Anteil, ",
      "konsistent mit circadianer Desorganisation ohne phototischen Zeitgeber. ",
      "Für künftige Studien empfiehlt sich: (a) Mindesttracking-Dauer ≥ 14 Tage pro Tier, ",
      "(b) standardisierte Gewichtsdokumentation (Eingangsgewicht, Entlassungsgewicht), ",
      "(c) Trennung Diagnosen-Kategorien bei Mehrfachdiagnosen."
    ), style = "Normal")

  out_e <- file.path(output_ordner, "Block8_RehabSuccess_Ergebnisbericht.docx")
  print(doc_e, target = out_e)
  cat("  ✓ Block8_RehabSuccess_Ergebnisbericht.docx\n")

} else {
  cat("  [WARN] officer/flextable nicht verfügbar.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# ZUSAMMENFASSUNG (Textdatei wie vorher)
# ══════════════════════════════════════════════════════════════════════════════
sink(file.path(output_ordner, "08_kennzahlen.txt"))
cat("══════════════════════════════════════════════════════════════\n")
cat("Block 8 — Rehabilitations-Erfolgs-Synthese\n")
cat(format(Sys.time(), "Erstellt: %Y-%m-%d %H:%M"), "\n")
cat("══════════════════════════════════════════════════════════════\n\n")
cat("── Stichprobe ──\n")
cat("Tiere gesamt:                    ", nrow(synthese), "\n")
cat("Tiere mit gültigem IS (≥", MIN_TAGE_IS_IV, "Tage): ", nrow(synthese[!is.na(IS)]), "\n")
cat("Tiere mit nacht_von_aktiv-Daten: ", nrow(synthese[!is.na(nacht_von_aktiv_mean)]),
    "(WS-Filter:", AKTIV_SCHWELLE * 100, "%)\n")
if (length(igel_ws_ausgeschl) > 0)
  cat("WS-Filter ausgeschlossen:        ", paste(igel_ws_ausgeschl, collapse=", "), "\n")
cat("HINWEIS: IV nicht analysiert — siehe Skript-Header.\n\n")
cat("── IS Diagnose-Gruppen ──\n")
print(grp_is_diag)
cat("\n── IS Geschlecht ──\n")
print(grp_is_sex)
cat("\n── Reha-Korrelationen (Spearman) ──\n")
print(korr_tabelle)
cat("\n── Nachtaktivität (PRIMÄRMETRIK: nacht_von_aktiv, WS-gefiltert) ──\n")
cat(sprintf("  Metrik: aktive Nachtmin / alle aktiven Min\n"))
cat(sprintf("  Winterschlaf-Filter: gesamt_aktiv < %.0f%%\n", AKTIV_SCHWELLE * 100))
if (length(igel_ws_ausgeschl) > 0)
  cat("  Ausgeschlossen:", paste(igel_ws_ausgeschl, collapse=", "), "\n")
if (!is.na(nacht_frueh) && !is.na(nacht_spaet)) {
  cat(sprintf("  Ø nacht_von_aktiv Tage 1–7:  %.1f%%\n", nacht_frueh * 100))
  cat(sprintf("  Ø nacht_von_aktiv ab Tag 14: %.1f%%\n", nacht_spaet * 100))
  cat(sprintf("  Veränderung:                 %+.1f Prozentpunkte\n",
              (nacht_spaet - nacht_frueh) * 100))
}
cat("\n══════════════════════════════════════════════════════════════\n")
sink()
cat("  ✓ 08_kennzahlen.txt\n")

cat("\n══════════════════════════════════════════════════════════════\n")
cat("Block 8 abgeschlossen!\n")
cat("  Primärmetrik:      nacht_von_aktiv (Analyse E)\n")
cat("                     = aktive Nachtmin / alle aktiven Min\n")
cat("                     WS-Filter:", AKTIV_SCHWELLE * 100, "% Gesamtaktivität\n")
if (length(igel_ws_ausgeschl) > 0)
  cat("                     WS-ausgeschlossen:", paste(igel_ws_ausgeschl, collapse=", "), "\n")
cat("  Sekundärmetrik:    IS (Analysen A–D)\n")
cat("  IS-Artefakt-Tiere: Igel5 (WS) + Igel1 (Kurztracking) aus A/B/D ausgeschlossen\n")
cat("  IV entfernt:       ~20-min Glättung → absolute Werte uninterpretierbar\n")
cat(sprintf("  Tage-Filter:       %d Tiere ausgeschlossen (< %d Tage)\n",
            n_ausgeschlossen, MIN_TAGE_IS_IV))
cat(sprintf("  Output:            %s\n", output_ordner))
cat("══════════════════════════════════════════════════════════════\n")
