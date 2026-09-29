# ==============================================================
# wmv vs. smoothed_wmv — method comparison
# ==============================================================
# Purpose: Empirical comparison of the two tRackIT classifications
#          pred_nested_loio_wmv          (300s rolling window)
#          pred_nested_loio_smoothed_wmv (300s + ~20 min smoothing)
#
# Questions:
#   1. How high is the agreement (concordance)?
#   2. How do IS and IV differ?
#   3. Which variable correlates more strongly with rehab success?
#   4. Does wmv look biologically more meaningful or noisier?
#
# Author: Natalie Steiner
# ==============================================================

pakete <- c("data.table", "lubridate", "ggplot2", "suncalc",
            "readxl", "scales", "patchwork", "writexl")
neu <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) install.packages(neu)
invisible(lapply(pakete, library, character.only = TRUE))

# ──────────────────────────────────────────────────────────────
# EINSTELLUNGEN
# ──────────────────────────────────────────────────────────────
projekt_root  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
daten_ordner  <- file.path(projekt_root, "data", "activity")
meta_datei    <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
output_ordner <- file.path(projekt_root, "output", "Block0_Pipeline", "wmv_Vergleich")
dir.create(output_ordner, recursive = TRUE, showWarnings = FALSE)

standort_lat <- 52.39729710523643
standort_lon <-  9.216876248766871
zeitzone     <- "Europe/Berlin"

# 3 repräsentative Tiere für Zeitreihenplot (werden automatisch gewählt falls leer)
beispiel_tiere <- c()  # z.B. c("Igel5", "Igel11", "Igel14") oder leer lassen

# ──────────────────────────────────────────────────────────────
# HILFSFUNKTIONEN
# ──────────────────────────────────────────────────────────────
modal_val <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_integer_)
  as.integer(names(sort(table(x), decreasing = TRUE))[1L])
}

# IS  (Interdaily Stability)
calc_IS <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 24) return(NA_real_)
  p   <- 24L
  xbar <- mean(x)
  if (xbar == 0 || var(x) == 0) return(NA_real_)
  h   <- n %/% p
  xh  <- sapply(0:(p - 1), function(i) mean(x[seq(i + 1, h * p, by = p)]))
  (n * sum((xh - xbar)^2)) / (p * sum((x - xbar)^2))
}

# IV  (Intradaily Variability)
calc_IV <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) return(NA_real_)
  xbar <- mean(x)
  if (xbar == 0 || var(x) == 0) return(NA_real_)
  (n * sum(diff(x)^2)) / ((n - 1) * sum((x - xbar)^2))
}

# ──────────────────────────────────────────────────────────────
# DATEN EINLESEN
# ──────────────────────────────────────────────────────────────
cat("── Einlesen der Rohdaten ──\n")

csv_dateien <- list.files(daten_ordner, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_dateien) == 0) stop("Keine CSV-Dateien in: ", daten_ordner)

meta_pfad <- path.expand(meta_datei)
meta <- as.data.table(read_xlsx(meta_pfad))

# Igel-Name aus Dateiname extrahieren
igel_name_aus_pfad <- function(pfad) {
  basis <- tools::file_path_sans_ext(basename(pfad))
  # "classification_active_passive_Igel5_300s_window_plus_filter" → "Igel5"
  m <- regmatches(basis, regexpr("Igel\\s*\\d+", basis, ignore.case = TRUE))
  if (length(m) == 0) return(basis)
  gsub("\\s+", "", m)  # "Igel 5" → "Igel5"
}

# Speicherschonend: jede CSV einzeln einlesen und SOFORT auf Minuten reduzieren.
# So werden nie alle Rohzeilen gleichzeitig gehalten — wichtig für sehr große
# Tiere wie Igel5 (~5,18 Mio. Rohzeilen inkl. Winterschlaf).
dt_list <- lapply(csv_dateien, function(f) {
  ig <- igel_name_aus_pfad(f)

  # Nur den Header prüfen, ob beide Klassifikations-Spalten vorhanden sind
  kopf <- names(fread(file = path.expand(f), nrows = 0L))
  has_wmv      <- "pred_nested_loio_wmv"          %in% kopf
  has_smoothed <- "pred_nested_loio_smoothed_wmv" %in% kopf
  if (!has_wmv || !has_smoothed) {
    cat("  WARNUNG:", ig, "— fehlende Spalte(n), wird übersprungen\n")
    return(NULL)
  }

  # Nur die drei benötigten Spalten laden (spart Speicher & Zeit)
  d <- fread(file = path.expand(f), encoding = "UTF-8",
             select = c("_time",
                        "pred_nested_loio_wmv",
                        "pred_nested_loio_smoothed_wmv"))

  # _time liegt vor als "2024-09-30 19:44:33.120311+00:00" (Leerzeichen + Offset,
  # inkl. Sekundenbruchteile) — robust mit lubridate parsen (wie Block0/Block2).
  d[, datetime := ymd_hms(`_time`, tz = "UTC", quiet = TRUE)]
  d <- d[!is.na(datetime)]
  d[, datetime := with_tz(datetime, zeitzone)]
  d <- d[!is.na(pred_nested_loio_wmv) & !is.na(pred_nested_loio_smoothed_wmv)]
  if (nrow(d) == 0L) {
    cat("  WARNUNG:", ig, "— keine gültigen Zeilen, wird übersprungen\n")
    return(NULL)
  }

  # Binäre Aktivität (a=1, p=0)
  d[, aktiv_wmv      := fifelse(pred_nested_loio_wmv          == "a", 1L, 0L)]
  d[, aktiv_smoothed := fifelse(pred_nested_loio_smoothed_wmv == "a", 1L, 0L)]

  # Minutenreduktion (Modalwert) bereits pro Tier
  d[, datetime_min := floor_date(datetime, "1 minute")]
  dmin <- d[, .(
    aktiv_wmv      = modal_val(aktiv_wmv),
    aktiv_smoothed = modal_val(aktiv_smoothed),
    n              = .N
  ), by = datetime_min]
  dmin[, igel := ig]

  cat(sprintf("  %-8s  Rohzeilen: %9d  →  Minuten: %7d\n",
              ig, nrow(d), nrow(dmin)))
  rm(d); invisible(gc(verbose = FALSE))

  dmin[, .(igel, datetime_min, aktiv_wmv, aktiv_smoothed, n)]
})

dt_min <- rbindlist(Filter(Negate(is.null), dt_list))
if (nrow(dt_min) == 0) stop("Keine verwertbaren Aktivitätsdaten eingelesen.")

dt_min[, date   := as.Date(datetime_min, tz = zeitzone)]
dt_min[, stunde := hour(datetime_min)]

cat("\n  Tiere eingelesen:      ", uniqueN(dt_min$igel), "\n")
cat("  Minuten nach Reduktion:", nrow(dt_min), "\n\n")

# ──────────────────────────────────────────────────────────────
# 1. KONKORDANZ
# ──────────────────────────────────────────────────────────────
cat("── 1. Konkordanz ──\n")

konkordanz <- dt_min[!is.na(aktiv_wmv) & !is.na(aktiv_smoothed), .(
  n_gesamt     = .N,
  n_gleich     = sum(aktiv_wmv == aktiv_smoothed),
  n_wmv_aktiv_smoothed_passiv = sum(aktiv_wmv == 1L & aktiv_smoothed == 0L),
  n_wmv_passiv_smoothed_aktiv = sum(aktiv_wmv == 0L & aktiv_smoothed == 1L)
), by = igel]

konkordanz[, pct_konkordanz          := round(n_gleich / n_gesamt * 100, 1)]
konkordanz[, pct_wmv_aktiver         := round(n_wmv_aktiv_smoothed_passiv / n_gesamt * 100, 2)]
konkordanz[, pct_smoothed_aktiver    := round(n_wmv_passiv_smoothed_aktiv / n_gesamt * 100, 2)]

gesamt_konk <- dt_min[!is.na(aktiv_wmv) & !is.na(aktiv_smoothed),
  .(pct = round(sum(aktiv_wmv == aktiv_smoothed) / .N * 100, 1))]

cat("  Gesamtkondordanz:", gesamt_konk$pct, "%\n")
cat("  Spannweite je Tier:", min(konkordanz$pct_konkordanz), "–",
    max(konkordanz$pct_konkordanz), "%\n\n")

print(konkordanz[order(pct_konkordanz)])

# ──────────────────────────────────────────────────────────────
# 2. IS und IV für beide Variablen
# ──────────────────────────────────────────────────────────────
cat("\n── 2. IS und IV je Tier ──\n")

metriken <- dt_min[, {
  # Stündliche Zeitreihe für IS/IV
  stunden_ts_wmv      <- tapply(aktiv_wmv,      paste(date, stunde), mean, na.rm = TRUE)
  stunden_ts_smoothed <- tapply(aktiv_smoothed, paste(date, stunde), mean, na.rm = TRUE)

  .(
    IS_wmv      = calc_IS(as.numeric(stunden_ts_wmv)),
    IS_smoothed = calc_IS(as.numeric(stunden_ts_smoothed)),
    IV_wmv      = calc_IV(as.numeric(stunden_ts_wmv)),
    IV_smoothed = calc_IV(as.numeric(stunden_ts_smoothed)),
    n_tage      = uniqueN(date)
  )
}, by = igel]

metriken[, IS_diff := IS_wmv      - IS_smoothed]
metriken[, IV_diff := IV_wmv      - IV_smoothed]

cat("\n  IS — Mittelwert wmv:", round(mean(metriken$IS_wmv,      na.rm=TRUE), 3),
    "| smoothed:", round(mean(metriken$IS_smoothed, na.rm=TRUE), 3))
cat("\n  IV — Mittelwert wmv:", round(mean(metriken$IV_wmv,      na.rm=TRUE), 3),
    "| smoothed:", round(mean(metriken$IV_smoothed, na.rm=TRUE), 3), "\n")

print(metriken[order(igel)])

# ──────────────────────────────────────────────────────────────
# 3. KORRELATION MIT GEWICHTSZUNAHME
# ──────────────────────────────────────────────────────────────
cat("\n── 3. Korrelation mit Gewichtszunahme ──\n")

# Flexibler Spaltenzugriff
find_col <- function(dt, pattern) {
  hits <- grep(pattern, names(dt), ignore.case = TRUE, value = TRUE)
  if (length(hits) == 0) return(NULL)
  hits[1]
}

# Tier-ID und Gewichtszunahme aus data_igel.xlsx:
#   individual (Tiername), weight_gain (Zunahme in g).
# Fallback: aus weigth_entry / tagging_weight berechnen, falls weight_gain fehlt.
col_ig   <- find_col(meta, "^individual$|^igel$|^id$|^tier|^name")
col_gain <- find_col(meta, "^weight_gain$|weight_gain|gewichtszunahme|zunahme")
col_ein  <- find_col(meta, "weigth_entry|weight_entry|eingang|intake")
col_aus  <- find_col(meta, "tagging_weight|weight_release|ausgang")

if (!is.null(col_ig) && (!is.null(col_gain) || (!is.null(col_ein) && !is.null(col_aus)))) {
  meta_gew <- meta[, .(
    igel            = trimws(as.character(get(col_ig))),
    gewicht_zunahme = if (!is.null(col_gain))
                        suppressWarnings(as.numeric(get(col_gain)))
                      else
                        suppressWarnings(as.numeric(get(col_aus)) - as.numeric(get(col_ein)))
  )]
  meta_gew[!grepl("^Igel", igel, ignore.case = TRUE), igel := paste0("Igel", igel)]

  dat_korr <- merge(metriken, meta_gew, by = "igel")
  dat_korr <- dat_korr[!is.na(gewicht_zunahme) & !is.na(IS_wmv)]

  if (nrow(dat_korr) >= 5) {
    r_IS_wmv      <- cor.test(dat_korr$IS_wmv,      dat_korr$gewicht_zunahme, method = "spearman")
    r_IS_smoothed <- cor.test(dat_korr$IS_smoothed, dat_korr$gewicht_zunahme, method = "spearman")
    r_IV_wmv      <- cor.test(dat_korr$IV_wmv,      dat_korr$gewicht_zunahme, method = "spearman")
    r_IV_smoothed <- cor.test(dat_korr$IV_smoothed, dat_korr$gewicht_zunahme, method = "spearman")

    cat(sprintf("  IS ~ Gewichtszunahme:  wmv ρ=%.3f p=%.3f | smoothed ρ=%.3f p=%.3f\n",
                r_IS_wmv$estimate, r_IS_wmv$p.value,
                r_IS_smoothed$estimate, r_IS_smoothed$p.value))
    cat(sprintf("  IV ~ Gewichtszunahme:  wmv ρ=%.3f p=%.3f | smoothed ρ=%.3f p=%.3f\n",
                r_IV_wmv$estimate, r_IV_wmv$p.value,
                r_IV_smoothed$estimate, r_IV_smoothed$p.value))
  } else {
    cat("  Zu wenige Tiere mit vollständigen Gewichtsdaten (n =", nrow(dat_korr), ")\n")
  }
} else {
  cat("  Metadaten-Spalten nicht erkannt — Korrelation übersprungen\n")
  cat("  Gefundene Spalten:", paste(names(meta), collapse = ", "), "\n")
}

# ──────────────────────────────────────────────────────────────
# 4. GRAFIKEN
# ──────────────────────────────────────────────────────────────
cat("\n── 4. Grafiken erstellen ──\n")

# ── Plot A: Konkordanz je Tier ──
p_konk <- ggplot(konkordanz, aes(x = reorder(igel, pct_konkordanz), y = pct_konkordanz)) +
  geom_col(fill = "#2E75B6", width = 0.7) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_text(aes(label = paste0(pct_konkordanz, "%")),
            hjust = -0.1, size = 3, color = "grey30") +
  coord_flip(ylim = c(80, 102)) +
  labs(title = "A — Konkordanz wmv vs. smoothed_wmv je Tier",
       subtitle = "Rote Linie = 95%-Schwelle",
       x = NULL, y = "Übereinstimmung (%)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# ── Plot B: IV-Vergleich Boxplot ──
iv_long <- melt(
  metriken[, .(igel, IV_wmv, IV_smoothed)],
  id.vars = "igel",
  variable.name = "Methode",
  value.name = "IV"
)
iv_long[, Methode := fifelse(Methode == "IV_wmv", "wmv (300s)", "smoothed_wmv (300s+20min)")]

p_iv <- ggplot(iv_long[!is.na(IV)], aes(x = Methode, y = IV, fill = Methode)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  geom_line(aes(group = igel), color = "grey60", alpha = 0.5, linewidth = 0.4,
            data = iv_long[!is.na(IV)]) +
  scale_fill_manual(values = c("wmv (300s)" = "#70AD47",
                                "smoothed_wmv (300s+20min)" = "#2E75B6")) +
  labs(title = "B — Intradaily Variability (IV)",
       subtitle = "Linien verbinden dasselbe Tier",
       x = NULL, y = "IV") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# ── Plot C: IS-Vergleich Boxplot ──
is_long <- melt(
  metriken[, .(igel, IS_wmv, IS_smoothed)],
  id.vars = "igel",
  variable.name = "Methode",
  value.name = "IS"
)
is_long[, Methode := fifelse(Methode == "IS_wmv", "wmv (300s)", "smoothed_wmv (300s+20min)")]

p_is <- ggplot(is_long[!is.na(IS)], aes(x = Methode, y = IS, fill = Methode)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  geom_line(aes(group = igel), color = "grey60", alpha = 0.5, linewidth = 0.4,
            data = is_long[!is.na(IS)]) +
  scale_fill_manual(values = c("wmv (300s)" = "#70AD47",
                                "smoothed_wmv (300s+20min)" = "#2E75B6")) +
  labs(title = "C — Interdaily Stability (IS)",
       subtitle = "Linien verbinden dasselbe Tier",
       x = NULL, y = "IS") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# ── Plot D: Zeitreihe für 3 Beispieltiere ──
if (length(beispiel_tiere) == 0) {
  alle_igel <- metriken[!is.na(IS_wmv) & n_tage >= 7, igel]
  if (length(alle_igel) >= 3) {
    set.seed(42)
    beispiel_tiere <- sample(alle_igel, 3)
  } else {
    beispiel_tiere <- alle_igel
  }
}
cat("  Zeitreihenplot für:", paste(beispiel_tiere, collapse = ", "), "\n")

ts_data <- dt_min[igel %in% beispiel_tiere & !is.na(aktiv_wmv) & !is.na(aktiv_smoothed)]

# Stündliche Mittelwerte für bessere Lesbarkeit
ts_stunde <- ts_data[, .(
  aktiv_wmv      = mean(aktiv_wmv,      na.rm = TRUE),
  aktiv_smoothed = mean(aktiv_smoothed, na.rm = TRUE)
), by = .(igel, date, stunde)]
ts_stunde[, datetime_h := as.POSIXct(paste(date, sprintf("%02d:00:00", stunde)),
                                       tz = zeitzone)]

ts_long <- melt(ts_stunde, id.vars = c("igel", "datetime_h"),
                measure.vars = c("aktiv_wmv", "aktiv_smoothed"),
                variable.name = "Methode", value.name = "aktiv")
ts_long[, Methode := fifelse(Methode == "aktiv_wmv", "wmv (300s)", "smoothed_wmv (300s+20min)")]

p_ts <- ggplot(ts_long, aes(x = datetime_h, y = aktiv, color = Methode)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  facet_wrap(~ igel, ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("wmv (300s)" = "#70AD47",
                                 "smoothed_wmv (300s+20min)" = "#2E75B6")) +
  scale_x_datetime(date_labels = "%d.%m", date_breaks = "3 days") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title = "D — Zeitreihenvergleich (stündliche Mittelwerte)",
       x = NULL, y = "Aktivitätsanteil", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#D5E8F0"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ── Plots speichern ──
ggsave(file.path(output_ordner, "Vergleich_Konkordanz.png"),  p_konk, width=8,  height=5,  dpi=200)
ggsave(file.path(output_ordner, "Vergleich_IV.png"),          p_iv,   width=6,  height=5,  dpi=200)
ggsave(file.path(output_ordner, "Vergleich_IS.png"),          p_is,   width=6,  height=5,  dpi=200)
ggsave(file.path(output_ordner, "Vergleich_Zeitreihe.png"),   p_ts,   width=12, height=8,  dpi=200)

# Kombinierter Überblicksplot (A+B+C)
p_kombi <- (p_konk | p_iv | p_is) +
  plot_annotation(
    title    = "wmv vs. smoothed_wmv — Methodenvergleich",
    subtitle = "Grün = wmv (300s)  |  Blau = smoothed_wmv (300s + ~20min Glättung)",
    theme    = theme(plot.title    = element_text(size=14, face="bold", color="#1F4E79"),
                     plot.subtitle = element_text(size=10, color="grey40"))
  )
ggsave(file.path(output_ordner, "Vergleich_Uebersicht.png"), p_kombi, width=16, height=6, dpi=200)

cat("  Plots gespeichert in:", output_ordner, "\n")

# ──────────────────────────────────────────────────────────────
# EXCEL EXPORT
# ──────────────────────────────────────────────────────────────
cat("\n── Excel Export ──\n")

sheets <- list(
  Konkordanz   = as.data.frame(konkordanz),
  IS_IV        = as.data.frame(metriken)
)

write_xlsx(sheets, file.path(output_ordner, "wmv_Vergleich_Ergebnisse.xlsx"))
cat("  Excel gespeichert\n")

# ──────────────────────────────────────────────────────────────
# ZUSAMMENFASSUNG
# ──────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════\n")
cat("ZUSAMMENFASSUNG\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("Gesamtkondordanz:     %.1f%%\n", gesamt_konk$pct))
cat(sprintf("IV wmv      Median:   %.3f\n", median(metriken$IV_wmv,      na.rm=TRUE)))
cat(sprintf("IV smoothed Median:   %.3f\n", median(metriken$IV_smoothed, na.rm=TRUE)))
cat(sprintf("IV-Differenz Median:  %.3f (wmv - smoothed)\n", median(metriken$IV_diff, na.rm=TRUE)))
cat(sprintf("IS wmv      Median:   %.3f\n", median(metriken$IS_wmv,      na.rm=TRUE)))
cat(sprintf("IS smoothed Median:   %.3f\n", median(metriken$IS_smoothed, na.rm=TRUE)))
cat(sprintf("IS-Differenz Median:  %.3f (wmv - smoothed)\n\n", median(metriken$IS_diff, na.rm=TRUE)))

cat("Ausgabedateien:\n")
cat(" ", file.path(output_ordner, "Vergleich_Uebersicht.png"), "\n")
cat(" ", file.path(output_ordner, "Vergleich_Zeitreihe.png"), "\n")
cat(" ", file.path(output_ordner, "wmv_Vergleich_Ergebnisse.xlsx"), "\n")
cat("══════════════════════════════════════════════════════\n")
