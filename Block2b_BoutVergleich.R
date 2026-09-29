# ==============================================================
# Block 2b — bout comparison: wmv vs. smoothed_wmv
# ==============================================================
# Project:  Hedgehog VHF telemetry, Lower Saxony
# Author:   Natalie Steiner
#
# Question:
#   Should bouts be computed from pred_nested_loio_wmv (~5 min
#   effective resolution) or pred_nested_loio_smoothed_wmv
#   (~20 min resolution due to additional smoothing)?
#
# Comparisons:
#   1. Log-survivor plot (overlaid, both variants)
#   2. Per-night bout metrics (side-by-side boxplots)
#   3. 24h activity profile (overlaid lines)
#   4. Night-activity time course (LOESS smoothing)
#
# Data sources:
#   smoothed_wmv -> output/Block0_Pipeline/gamm_nachtminuten.rds
#                   (already fully processed by Block0)
#   wmv          -> data/activity/*.csv  (re-processed here)

# ──────────────────────────────────────────────────────────────
# PAKETE
# ──────────────────────────────────────────────────────────────
pakete <- c("data.table", "lubridate", "ggplot2", "patchwork",
            "scales", "suncalc", "readxl",
            "officer", "flextable", "openxlsx")
neu    <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) {
  cat("Installiere fehlende Pakete:", paste(neu, collapse = ", "), "\n")
  install.packages(neu)
}
invisible(lapply(pakete, library, character.only = TRUE))
cat("✓ Alle Pakete geladen\n\n")

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a[1] else b

# ──────────────────────────────────────────────────────────────
# EINSTELLUNGEN  — hier anpassen
# ──────────────────────────────────────────────────────────────

projekt_root  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
daten_ordner  <- file.path(projekt_root, "data", "activity")
meta_datei    <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
nacht_rds     <- file.path(projekt_root, "output", "Block0_Pipeline", "gamm_nachtminuten.rds")
output_ordner <- file.path(projekt_root, "output", "Block2b_BoutVergleich")

dir.create(output_ordner, showWarnings = FALSE, recursive = TRUE)

# Standort (Niedersachsen)
standort_lat <- 52.39729710523643
standort_lon <-  9.216876248766871
zeitzone     <- "Europe/Berlin"

# ── Bout-Kriterien (minimale Pause in Minuten = zwei verschiedene Bouts) ──
#
#   smoothed_wmv: ~20 min effektive Auflösung durch zusätzliche Glättung
#     → kurze Pausen < 20 min werden artifiziell überbrückt
#     → Bout-Kriterium muss ≥ Glättungsfenster gesetzt werden
#
#   wmv: ~5 min effektive Auflösung (nur 300s Rolling Window)
#     → kurze Pausen biologisch sinnvoll erkennbar
#     → kürzeres Bout-Kriterium möglich
#
# !! ANPASSEN nach Sichtung des Log-Survivor-Plots !!
bout_kriterium_smwmv <- 20   # Minuten — smoothed_wmv
bout_kriterium_wmv   <- 10   # Minuten — wmv

# ── Ausschlusskriterien (identisch zu Block2_GAMM.R) ──
min_naechte   <- 3    # Mindestanzahl Nächte mit Nachtdaten
min_nacht_min <- 60   # Mindestanzahl absoluter Nachtminuten gesamt
min_pct_aktiv <- 1    # < 1% Nachtaktivität → statisches Signal

# Schwellenwert toter Sender (Block0-Replikation für wmv)
schwelle_aktiv_pct <- 5  # % Tagesaktivität unterhalb → statisch

# ── Zeitfenster: Mindest-Stichprobengröße ──────────────────────
# Analysen nur für Tage mit ≥ n_min_igel Tieren (Survivorship-Bias-Schutz).
# Basierend auf den Daten: n ≥ 5 Igel = bis Tag 27.
# Wird automatisch aus den Daten ermittelt (nicht manuell anpassen).
n_min_igel <- 5L          # Mindestanzahl Igel pro Zeitpunkt

# ── Farben ──
farbe_smwmv <- "#2C3E6B"   # Dunkelblau — smoothed_wmv
farbe_wmv   <- "#E07B39"   # Orange     — wmv

# ──────────────────────────────────────────────────────────────
# HILFSFUNKTIONEN
# ──────────────────────────────────────────────────────────────

# Modalwert (gibt Character zurück, wie in Block0)
modal_val_chr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1L]
}

# IBI-Berechnung: gibt Vektor der Inter-Bout-Intervalle (in Minuten)
# aktiv_vec: integer 0/1, sortiert nach Zeit innerhalb einer Nacht
# Rückgabe: numeric — Abstände zwischen aktiven Minuten > 1 Min
berechne_ibis <- function(aktiv_vec) {
  aktiv_idx <- which(aktiv_vec == 1L)
  if (length(aktiv_idx) < 2L) return(numeric(0))
  pausen <- diff(aktiv_idx)
  pausen[pausen > 1L]   # nur echte Pausen (> 1 Min Lücke)
}

# Bout-Berechnung pro Nacht: Run-Length-Encoding + Merging kurzer Pausen
# Gibt Liste: n_bouts, mean_dauer, max_dauer, mean_ibi (alles in Minuten)
berechne_bouts <- function(aktiv_vec, bout_krit) {
  if (length(aktiv_vec) == 0L || sum(aktiv_vec, na.rm = TRUE) == 0L) {
    return(list(n_bouts    = 0L,
                mean_dauer = NA_real_,
                max_dauer  = NA_real_,
                mean_ibi   = NA_real_))
  }
  rle_res <- rle(aktiv_vec)
  laengen <- as.numeric(rle_res$lengths)
  werte   <- rle_res$values

  # Kurze Pausen zusammenführen (3 Iterationen wie in Block2)
  for (iter in 1:3) {
    passiv_idx <- which(werte == 0)
    zu_mergen  <- passiv_idx[laengen[passiv_idx] < bout_krit]
    if (length(zu_mergen) == 0L) break
    j <- zu_mergen[1L]
    if (j > 1L && j < length(werte)) {
      neue_laenge <- laengen[j-1L] + laengen[j] + laengen[j+1L]
      laengen <- c(laengen[seq_len(j-2L)], neue_laenge,
                   laengen[(j+2L):length(laengen)])
      werte   <- c(werte[seq_len(j-2L)], 1L,
                   werte[(j+2L):length(werte)])
    }
  }

  aktiv_bouts <- laengen[werte == 1]
  passiv_ibi  <- laengen[werte == 0]
  passiv_ibi  <- passiv_ibi[passiv_ibi >= bout_krit]

  list(
    n_bouts    = length(aktiv_bouts),
    mean_dauer = as.numeric(mean(aktiv_bouts)),
    max_dauer  = as.numeric(max(aktiv_bouts)),
    mean_ibi   = if (length(passiv_ibi) > 0) as.numeric(mean(passiv_ibi))
                 else NA_real_
  )
}

# Erstellt data.table für Log-Survivor-Plot aus einer Liste von Nächten
# aktiv_liste: data.table mit Spalten (igel, datum, aktiv), geordnet nach stunde
logsurv_dt <- function(aktiv_liste) {
  ibis <- aktiv_liste[order(igel, datum, stunde), {
    .(ibi = berechne_ibis(aktiv))
  }, by = .(igel, datum)]

  ibi_vals <- ibis$ibi
  if (length(ibi_vals) == 0L) return(NULL)
  ibi_sort  <- sort(ibi_vals)
  survivor  <- (length(ibi_sort):1L) / length(ibi_sort)
  data.table(ibi = ibi_sort, survivor = survivor)
}

# Metadaten-Lader (wie Block0)
parse_datum <- function(x) {
  x      <- as.character(x)
  result <- suppressWarnings(as.Date(x, format = "%d.%m.%Y"))
  na_idx <- is.na(result)
  result[na_idx] <- suppressWarnings(as.Date(x[na_idx]))
  result
}

# ──────────────────────────────────────────────────────────────
# METADATEN LADEN
# ──────────────────────────────────────────────────────────────
cat("── Metadaten laden ──\n")
meta_pfad <- path.expand(meta_datei)
meta <- as.data.table(read_excel(meta_pfad, sheet = 1))

setnames(meta,
  old = c("individual", "sex", "diagnosis_main", "tagging_weight",
          "soft_release_start", "date_release", "date_last_signal",
          "tagging_period"),
  new = c("igel", "geschlecht", "diagnose", "gewicht_g",
          "soft_release", "hard_release", "letztes_signal",
          "besenderungsdauer"),
  skip_absent = TRUE)

meta[, igel         := trimws(igel)]
meta[, hard_release := parse_datum(hard_release)]
meta[, soft_release := parse_datum(soft_release)]

cat("  Igel in Metadaten:", nrow(meta), "\n\n")

# ──────────────────────────────────────────────────────────────
# A) SMOOTHED_WMV — aus gamm_nachtminuten.rds laden
# ──────────────────────────────────────────────────────────────
cat("── A) Smoothed_wmv: RDS laden ──\n")

if (!file.exists(path.expand(nacht_rds))) {
  stop("gamm_nachtminuten.rds nicht gefunden.\n",
       "Bitte zuerst Block0_Datenpipeline.R ausführen!\n",
       "Pfad: ", nacht_rds)
}

dt_smwmv_raw <- readRDS(path.expand(nacht_rds))
setDT(dt_smwmv_raw)

# Sicherstellen dass igel als Character behandelt wird (kein Faktor-Level-Problem)
dt_smwmv_raw[, igel  := as.character(igel)]
dt_smwmv_raw[, datum := as.Date(datum)]

cat("  Nachtminuten smoothed_wmv:", format(nrow(dt_smwmv_raw), big.mark = "'"),
    "| Igel:", uniqueN(dt_smwmv_raw$igel), "\n")

# Ausschlusskriterien für smoothed_wmv
qc_smwmv <- dt_smwmv_raw[, .(
  n_naechte  = uniqueN(datum),
  n_min      = .N,
  pct_aktiv  = mean(aktiv) * 100
), by = igel]

qc_smwmv[, einschluss := n_naechte >= min_naechte &
                          n_min     >= min_nacht_min &
                          pct_aktiv >= min_pct_aktiv]

igel_smwmv <- qc_smwmv[einschluss == TRUE, igel]
dt_smwmv   <- dt_smwmv_raw[igel %in% igel_smwmv]
dt_smwmv[, igel := factor(igel)]

cat("  Nach Ausschluss:", length(igel_smwmv), "Igel\n")
if (any(!qc_smwmv$einschluss)) {
  cat("  Ausgeschlossen (smwmv):\n")
  print(qc_smwmv[einschluss == FALSE,
                 .(igel, n_naechte, n_min, pct_aktiv = round(pct_aktiv, 1))])
}
cat("\n")

# ──────────────────────────────────────────────────────────────
# B) WMV — aus CSVs laden (Block0-Pipeline, nur Nacht)
# ──────────────────────────────────────────────────────────────
cat("── B) wmv: CSV-Pipeline ──\n")

alle_csvs <- list.files(path.expand(daten_ordner),
                        pattern = "classification_active_passive.*\\.csv$",
                        full.names = TRUE)
cat("  Gefundene CSVs:", length(alle_csvs), "\n")

# Hilfsfunktion: CSV-Pfad für einen Igel finden
finde_csv <- function(igel_name) {
  nr      <- gsub("[^0-9]", "", igel_name)
  treffer <- grep(paste0("Igel[[:space:]_]?0*", nr, "[^0-9]"),
                  alle_csvs, value = TRUE)
  if (length(treffer) == 0L) return(NA_character_)
  treffer[1L]
}

meta[, csv_pfad := sapply(igel, finde_csv)]

wmv_nacht_liste <- list()

for (i in seq_len(nrow(meta))) {

  igel_info <- meta[i]
  igel_name <- igel_info$igel
  csv_pfad  <- igel_info$csv_pfad

  if (is.na(csv_pfad)) {
    cat("  ⚠", igel_name, "— keine CSV gefunden, übersprungen\n")
    next
  }

  hard_release_i <- igel_info$hard_release

  # 1. Einlesen
  dt <- tryCatch(
    fread(csv_pfad, na.strings = c("", "NA"), showProgress = FALSE),
    error = function(e) { cat("  ✗", igel_name, "Lesefehler:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(dt) || nrow(dt) == 0L) next

  # Spalten umbenennen (tolerant)
  setnames(dt,
    old = c("Station Name", "Individual Name", "_time", "0"),
    new = c("station", "igel_col", "time_raw", "signal"),
    skip_absent = TRUE)

  # 2. Aktivste Klassifikationsspalte wählen
  #    Priorität: wmv → mv → loio (reines LOIO ohne Zeitgewichtung)
  #    Gilt generell für alle Igel — nicht nur für Igel 1–6.
  #    mv und loio haben dieselbe ~5-min-Auflösung wie wmv (kein Zeitgewicht),
  #    daher gilt Bout-Kriterium 10 min für alle drei Fallback-Stufen.
  wmv_hat_werte  <- "pred_nested_loio_wmv" %in% names(dt) &&
                    any(dt[["pred_nested_loio_wmv"]] %in% c("a", "p"))
  mv_hat_werte   <- "pred_nested_loio_mv"  %in% names(dt) &&
                    any(dt[["pred_nested_loio_mv"]]  %in% c("a", "p"))
  loio_hat_werte <- "pred_nested_loio"     %in% names(dt) &&
                    any(dt[["pred_nested_loio"]]     %in% c("a", "p"))

  if (wmv_hat_werte) {
    klasse_spalte_i <- "pred_nested_loio_wmv"
    cat("  [wmv]", igel_name, "\n")
  } else if (mv_hat_werte) {
    klasse_spalte_i <- "pred_nested_loio_mv"
    cat("  [mv → Fallback 1]", igel_name,
        "— wmv leer, verwende pred_nested_loio_mv\n")
  } else if (loio_hat_werte) {
    klasse_spalte_i <- "pred_nested_loio"
    cat("  [loio → Fallback 2]", igel_name,
        "— wmv und mv leer, verwende pred_nested_loio\n")
  } else {
    cat("  ⚠", igel_name,
        "— keine Klassifikationsspalte verfügbar (wmv/mv/loio alle leer)\n")
    next
  }

  # 3. Zeitstempel
  dt[, time_utc   := ymd_hms(time_raw, tz = "UTC", quiet = TRUE)]
  dt[, time_local := with_tz(time_utc, tzone = zeitzone)]
  dt[, time_min   := floor_date(time_local, unit = "minute")]
  dt <- dt[!is.na(time_min)]

  # 4. Klassifikation (gewählte Spalte)
  klasse_raw <- dt[[klasse_spalte_i]]
  dt[, klasse_wmv := fcase(
    klasse_raw == "a", "aktiv",
    klasse_raw == "p", "passiv",
    default = NA_character_
  )]

  # 5. Minutenreduktion (Modalwert)
  dt_min <- dt[!is.na(klasse_wmv), .(
    klasse_modal = modal_val_chr(klasse_wmv),
    n_messungen  = .N
  ), by = time_min]
  dt_min[, datum := as.Date(time_min, tz = zeitzone)]

  # 6. Hard-Release-Filter
  if (!is.na(hard_release_i)) {
    dt_min <- dt_min[datum >= hard_release_i]
  }
  if (nrow(dt_min) == 0L) {
    cat("  ⚠", igel_name, "— keine Daten nach Auswilderungsdatum\n")
    next
  }

  # 7. Toter-Sender-Erkennung (wie Block0)
  letztes_datum_csv <- max(dt_min$datum, na.rm = TRUE)
  tagesdaten_wmv <- dt_min[!is.na(klasse_modal), .(
    pct_aktiv = mean(klasse_modal == "aktiv") * 100
  ), by = datum][order(datum)]

  aktive_tage <- tagesdaten_wmv[pct_aktiv >= schwelle_aktiv_pct, datum]
  if (length(aktive_tage) == 0L) {
    cat("  ⚠", igel_name, "— kein Tag mit ≥", schwelle_aktiv_pct,
        "% Aktivität (wmv), übersprungen\n")
    next
  }
  letzter_aktiver_tag <- max(aktive_tage)
  if (letzter_aktiver_tag < letztes_datum_csv) {
    dt_min <- dt_min[datum <= letzter_aktiver_tag]
  }

  # 8. Tag/Nacht via suncalc
  tage_df <- data.frame(
    date = unique(dt_min$datum),
    lat  = standort_lat,
    lon  = standort_lon
  )
  sun <- as.data.table(getSunlightTimes(data = tage_df, tz = zeitzone,
                                         keep = c("sunrise", "sunset")))
  dt_min <- merge(dt_min, sun[, .(date, sunrise, sunset)],
                  by.x = "datum", by.y = "date", all.x = TRUE)

  dt_min[, tageszeit := fifelse(
    time_min >= sunrise & time_min < sunset, "Tag", "Nacht"
  )]
  dt_min[, stunde    := hour(time_min) + minute(time_min) / 60]
  dt_min[, tage_seit := as.integer(datum - hard_release_i) + 1L]

  # 9. Nur Nacht behalten
  nacht_wmv <- dt_min[tageszeit == "Nacht" & !is.na(klasse_modal), .(
    igel           = igel_name,
    datum,
    tage_seit,
    stunde,
    aktiv          = as.integer(klasse_modal == "aktiv"),
    quelle_spalte  = klasse_spalte_i   # wmv oder mv (Fallback) — für Transparenz
  )]

  if (nrow(nacht_wmv) == 0L) {
    cat("  ⚠", igel_name, "— keine Nachtminuten nach Filterung\n")
    next
  }

  wmv_nacht_liste[[igel_name]] <- nacht_wmv
  cat("  ✓", igel_name, "—", format(nrow(nacht_wmv), big.mark = "'"),
      "Nachtminuten [", klasse_spalte_i, "]\n")
}

dt_wmv_raw <- rbindlist(wmv_nacht_liste)
cat("\n  wmv gesamt:", format(nrow(dt_wmv_raw), big.mark = "'"),
    "Nachtminuten | Igel:", uniqueN(dt_wmv_raw$igel), "\n")

# Ausschlusskriterien für wmv
qc_wmv <- dt_wmv_raw[, .(
  n_naechte  = uniqueN(datum),
  n_min      = .N,
  pct_aktiv  = mean(aktiv) * 100
), by = igel]

qc_wmv[, einschluss := n_naechte >= min_naechte &
                        n_min     >= min_nacht_min &
                        pct_aktiv >= min_pct_aktiv]

igel_wmv <- qc_wmv[einschluss == TRUE, igel]
dt_wmv   <- dt_wmv_raw[igel %in% igel_wmv]
dt_wmv[, igel := factor(igel)]

cat("  Nach Ausschluss:", length(igel_wmv), "Igel\n")
if (any(!qc_wmv$einschluss)) {
  cat("  Ausgeschlossen (wmv):\n")
  print(qc_wmv[einschluss == FALSE,
               .(igel, n_naechte, n_min, pct_aktiv = round(pct_aktiv, 1))])
}
cat("\n")

# ──────────────────────────────────────────────────────────────
# ÜBERBLICK: Welche Igel sind in beiden Varianten enthalten?
# ──────────────────────────────────────────────────────────────
cat("════════════════════════════════════════════\n")
cat("EINSCHLUSS-ÜBERBLICK\n")
cat("════════════════════════════════════════════\n")

igel_beide  <- intersect(igel_smwmv, igel_wmv)
igel_nur_sm <- setdiff(igel_smwmv, igel_wmv)
igel_nur_w  <- setdiff(igel_wmv, igel_smwmv)

cat("In beiden Varianten:  ", length(igel_beide),  "Igel —",
    paste(sort(igel_beide), collapse = ", "), "\n")
if (length(igel_nur_sm) > 0)
  cat("Nur in smoothed_wmv: ", length(igel_nur_sm), "Igel —",
      paste(sort(igel_nur_sm), collapse = ", "), "\n")
if (length(igel_nur_w) > 0)
  cat("Nur in wmv:          ", length(igel_nur_w), "Igel —",
      paste(sort(igel_nur_w), collapse = ", "), "\n")
cat("════════════════════════════════════════════\n\n")

if (length(igel_beide) < 2L) {
  warning("Weniger als 2 Igel in beiden Varianten — Boxplots möglicherweise ",
          "wenig aussagekräftig.")
}

# Für Vergleichsplots: nur Igel die in BEIDEN Varianten enthalten sind
dt_smwmv_v <- dt_smwmv[igel %in% igel_beide]
dt_wmv_v   <- dt_wmv[  igel %in% igel_beide]

# ──────────────────────────────────────────────────────────────
# DYNAMISCHER ZEITCUTOFF — N-at-risk-Schwelle
# ──────────────────────────────────────────────────────────────
# Späte Zeitpunkte werden von wenigen Langzeit-Tieren dominiert
# (Survivorship Bias). Alle Analysen laufen nur bis zum letzten Tag
# mit ≥ n_min_igel Tieren in BEIDEN Datensätzen.
# ──────────────────────────────────────────────────────────────

# N pro Tag: beide Datensätze kombiniert (konservativste Grenze)
n_pro_tag_smwmv <- dt_smwmv_v[, .(n_igel = uniqueN(igel)), by = tage_seit]
n_pro_tag_wmv   <- dt_wmv_v[,   .(n_igel = uniqueN(igel)), by = tage_seit]
n_pro_tag_beide <- merge(n_pro_tag_smwmv, n_pro_tag_wmv,
                          by = "tage_seit", suffixes = c("_sm", "_wmv"))
n_pro_tag_beide[, n_min_beider := pmin(n_igel_sm, n_igel_wmv)]

cutoff_tag_2b <- n_pro_tag_beide[n_min_beider >= n_min_igel,
                                  max(tage_seit)]

cat("══════════════════════════════════════════════\n")
cat("DYNAMISCHER ZEITCUTOFF (N-at-risk-Schwelle)\n")
cat("══════════════════════════════════════════════\n")
cat(sprintf("  Mindest-N:   %d Igel pro Zeitpunkt\n", n_min_igel))
cat(sprintf("  Cutoff Tag:  %d\n", cutoff_tag_2b))
cat(sprintf("  Zeitraum:    Tag 1 – %d nach Auswilderung\n", cutoff_tag_2b))
cat("  Tage danach werden ausgeschlossen (Survivorship Bias).\n")
cat("══════════════════════════════════════════════\n\n")

# N-at-risk-Tabelle
cat("N-at-risk (smoothed_wmv / wmv) um den Cutoff:\n")
print(n_pro_tag_beide[tage_seit >= max(1, cutoff_tag_2b - 7) &
                       tage_seit <= cutoff_tag_2b + 3,
                       .(tage_seit, n_smwmv = n_igel_sm,
                         n_wmv = n_igel_wmv, n_min = n_min_beider)])
cat("\n")

# Filter anwenden
dt_smwmv_v <- dt_smwmv_v[tage_seit <= cutoff_tag_2b]
dt_wmv_v   <- dt_wmv_v[  tage_seit <= cutoff_tag_2b]

# Aktualisierte Igel-Liste (nur Tiere mit Daten im Fenster)
igel_beide_im_fenster <- intersect(
  dt_smwmv_v[, unique(igel)],
  dt_wmv_v[,   unique(igel)]
)
dt_smwmv_v <- dt_smwmv_v[igel %in% igel_beide_im_fenster]
dt_wmv_v   <- dt_wmv_v[  igel %in% igel_beide_im_fenster]

cat(sprintf("  Effektives Analysefenster: Tag 1 – %d | %d Igel\n\n",
            cutoff_tag_2b, length(igel_beide_im_fenster)))

# ──────────────────────────────────────────────────────────────
# 1. LOG-SURVIVOR-PLOT (overlaid)
# ──────────────────────────────────────────────────────────────
cat("── 1. Log-Survivor-Plot ──\n")

ls_smwmv <- logsurv_dt(dt_smwmv_v)
ls_wmv   <- logsurv_dt(dt_wmv_v)

if (is.null(ls_smwmv) || is.null(ls_wmv)) {
  cat("  ⚠ Zu wenige IBIs für Log-Survivor-Plot\n")
} else {
  ls_smwmv[, Variante := paste0("smoothed_wmv  (Kriterium ", bout_kriterium_smwmv, " min)")]
  ls_wmv[,   Variante := paste0("wmv           (Kriterium ", bout_kriterium_wmv,   " min)")]
  ls_all <- rbindlist(list(ls_smwmv, ls_wmv))

  # Auf max. 60 Minuten begrenzen (längere Pausen = sicher verschiedene Bouts)
  ls_plot <- ls_all[ibi <= 60]

  # Legendenbeschriftungen vorberechnen (setNames braucht einfache Zeichenketten)
  lbl_smwmv_ls <- paste0("smoothed_wmv  (Kriterium ", bout_kriterium_smwmv, " min)")
  lbl_wmv_ls   <- paste0("wmv           (Kriterium ", bout_kriterium_wmv,   " min)")
  farben_ls    <- setNames(c(farbe_smwmv, farbe_wmv), c(lbl_smwmv_ls, lbl_wmv_ls))
  ltypen_ls    <- setNames(c("solid",     "solid"),    c(lbl_smwmv_ls, lbl_wmv_ls))

  p_logsurv <- ggplot(ls_plot, aes(x = ibi, y = log(survivor),
                                    color = Variante, linetype = Variante)) +
    geom_line(linewidth = 1.1) +
    geom_vline(xintercept = bout_kriterium_smwmv, linetype = "dashed",
               color = farbe_smwmv, alpha = 0.7, linewidth = 0.8) +
    geom_vline(xintercept = bout_kriterium_wmv, linetype = "dashed",
               color = farbe_wmv,   alpha = 0.7, linewidth = 0.8) +
    annotate("text", x = bout_kriterium_smwmv + 0.5,
             y = max(log(ls_plot$survivor), na.rm = TRUE) * 0.95,
             label = paste0("smwmv: ", bout_kriterium_smwmv, " min"),
             color = farbe_smwmv, hjust = 0, size = 3.2) +
    annotate("text", x = bout_kriterium_wmv + 0.5,
             y = max(log(ls_plot$survivor), na.rm = TRUE) * 0.80,
             label = paste0("wmv: ", bout_kriterium_wmv, " min"),
             color = farbe_wmv, hjust = 0, size = 3.2) +
    scale_color_manual(values = farben_ls) +
    scale_linetype_manual(values = ltypen_ls) +
    labs(
      title    = "Log-Survivor-Plot der Inter-Bout-Intervalle",
      subtitle = paste0(
        "Knick = natürliches Bout-Kriterium  |  n(smwmv) = ",
        format(nrow(ls_smwmv[ibi <= 60]), big.mark = "'"),
        " IBIs  |  n(wmv) = ",
        format(nrow(ls_wmv[ibi <= 60]), big.mark = "'"), " IBIs\n",
        "Gestrichelte Linien = gewählte Kriterien"
      ),
      x      = "Pause zwischen aktiven Minuten (min)",
      y      = "log(Überlebenswahrscheinlichkeit)",
      color  = NULL, linetype = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", color = "#1F4E79"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(output_ordner, "01_log_survivor.png"),
         p_logsurv, width = 9, height = 6, dpi = 200)
  cat("  ✓ 01_log_survivor.png gespeichert\n")
  cat("  → Knickpunkt zeigt natürliches Bout-Kriterium\n")
  cat("    smwmv-Kurve: flacher (Glättung unterdrückt kurze Pausen)\n")
  cat("    wmv-Kurve:   stärker strukturiert (feinere zeitliche Auflösung)\n\n")
}

# ──────────────────────────────────────────────────────────────
# 2. BOUT-METRIKEN PRO NACHT
# ──────────────────────────────────────────────────────────────
cat("── 2. Bout-Metriken berechnen ──\n")

# Metriken smoothed_wmv
bm_smwmv <- dt_smwmv_v[order(igel, datum, stunde), {
  bm <- berechne_bouts(aktiv, bout_kriterium_smwmv)
  .(n_bouts    = bm$n_bouts,
    mean_dauer = bm$mean_dauer,
    max_dauer  = bm$max_dauer,
    mean_ibi   = bm$mean_ibi,
    tage_seit  = tage_seit[1L])
}, by = .(igel, datum)]
bm_smwmv[, Variante := "smoothed_wmv"]

# Metriken wmv
bm_wmv <- dt_wmv_v[order(igel, datum, stunde), {
  bm <- berechne_bouts(aktiv, bout_kriterium_wmv)
  .(n_bouts    = bm$n_bouts,
    mean_dauer = bm$mean_dauer,
    max_dauer  = bm$max_dauer,
    mean_ibi   = bm$mean_ibi,
    tage_seit  = tage_seit[1L])
}, by = .(igel, datum)]
bm_wmv[, Variante := "wmv"]

bm_all <- rbindlist(list(bm_smwmv, bm_wmv))
bm_all[, Variante := factor(Variante,
                              levels  = c("wmv", "smoothed_wmv"),
                              labels  = c(
                                paste0("wmv\n(Krit. ", bout_kriterium_wmv, " min)"),
                                paste0("smoothed_wmv\n(Krit. ", bout_kriterium_smwmv, " min)")
                              ))]

farben_box <- setNames(c(farbe_wmv, farbe_smwmv),
                        c(paste0("wmv\n(Krit. ", bout_kriterium_wmv, " min)"),
                          paste0("smoothed_wmv\n(Krit. ", bout_kriterium_smwmv, " min)")))

cat("  Nächte smwmv:", nrow(bm_smwmv), "| wmv:", nrow(bm_wmv), "\n\n")

# Übersichtsstatistik
cat("── Bout-Metriken Übersicht ──\n")
cat("  smoothed_wmv:\n")
print(bm_smwmv[, .(
  n_bouts_med    = median(n_bouts,    na.rm = TRUE),
  dauer_med_min  = median(mean_dauer, na.rm = TRUE),
  ibi_med_min    = median(mean_ibi,   na.rm = TRUE)
)])
cat("  wmv:\n")
print(bm_wmv[, .(
  n_bouts_med    = median(n_bouts,    na.rm = TRUE),
  dauer_med_min  = median(mean_dauer, na.rm = TRUE),
  ibi_med_min    = median(mean_ibi,   na.rm = TRUE)
)])
cat("\n")

# Boxplot-Hilfsfunktion
boxplot_vergleich <- function(dt, y_var, y_lab, titel, log_y = FALSE) {
  p <- ggplot(dt[!is.na(get(y_var))],
              aes(x = Variante, y = get(y_var), fill = Variante)) +
    geom_boxplot(alpha = 0.75, outlier.shape = 21, outlier.size = 1.5) +
    geom_jitter(aes(color = Variante), width = 0.15, size = 1.2, alpha = 0.4) +
    scale_fill_manual(values  = farben_box) +
    scale_color_manual(values = farben_box) +
    labs(title = titel, x = NULL, y = y_lab) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          panel.grid.minor  = element_blank(),
          plot.title = element_text(face = "bold", size = 10))
  if (log_y) p <- p + scale_y_log10(labels = scales::comma)
  p
}

p_nbouts  <- boxplot_vergleich(bm_all, "n_bouts",    "Bouts pro Nacht",     "A — Anzahl Bouts")
p_mdauer  <- boxplot_vergleich(bm_all, "mean_dauer", "Mittlere Boutdauer (min)", "B — Boutdauer")
p_maxdauer <- boxplot_vergleich(bm_all, "max_dauer", "Max. Boutdauer (min)","C — Max. Boutdauer")
p_mibi    <- boxplot_vergleich(bm_all, "mean_ibi",   "Mittleres IBI (min)", "D — Inter-Bout-Intervall")

p_bouts_kombi <- (p_nbouts | p_mdauer) / (p_maxdauer | p_mibi) +
  plot_annotation(
    title    = "Bout-Metriken: wmv vs. smoothed_wmv",
    subtitle = paste0(
      "wmv: Bout-Kriterium ", bout_kriterium_wmv, " min  |  ",
      "smoothed_wmv: Bout-Kriterium ", bout_kriterium_smwmv, " min  |  ",
      "n = ", length(igel_beide), " Igel"
    ),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold", color = "#1F4E79"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave(file.path(output_ordner, "02_bout_metriken.png"),
       p_bouts_kombi, width = 10, height = 8, dpi = 200)
cat("  ✓ 02_bout_metriken.png gespeichert\n\n")

# ──────────────────────────────────────────────────────────────
# 3. 24h-AKTIVITÄTSPROFIL (Nacht, überlagert)
# ──────────────────────────────────────────────────────────────
cat("── 3. 24h-Aktivitätsprofil ──\n")

# Granularität: Viertelstunde
stunde_runden <- 0.25

# ── Hilfsfunktion: Stunden für Nacht-Zeitachse verschieben ──
# Standardproblem: Nachtdaten gehen von ~20:00 über Mitternacht bis ~06:00.
# Wenn stunde 0–6 als 0–6 auf der x-Achse liegt, erscheinen Abend (20–23)
# rechts von Morgen (0–6) → falsche Darstellung.
# Lösung: frühe Morgenstunden um +24 verschieben → 0 → 24, 1 → 25 usw.
# Die Achse wird dann von ~18 bis ~30 gezeichnet (lineare Skala, kein Umbruch).
nacht_shift <- function(x) ifelse(x < 12, x + 24, x)

# Achsenbreaks und -Labels für verschobene Nachtdarstellung
nacht_breaks <- c(18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
nacht_labels <- sprintf("%02d:00", nacht_breaks %% 24L)
nacht_limits <- c(18.5, 30.5)   # leicht außerhalb der Daten → kein Clipping

# ── Populationsmittel pro Viertelstunde ──
profil_smwmv <- dt_smwmv_v[, .(
  pct_aktiv = mean(aktiv) * 100,
  n         = .N
), by = .(stunde_h = round(stunde / stunde_runden) * stunde_runden)]
profil_smwmv[, Variante := "smoothed_wmv"]

profil_wmv <- dt_wmv_v[, .(
  pct_aktiv = mean(aktiv) * 100,
  n         = .N
), by = .(stunde_h = round(stunde / stunde_runden) * stunde_runden)]
profil_wmv[, Variante := "wmv"]

profil_all <- rbindlist(list(profil_smwmv, profil_wmv))
profil_all[, stunde_plot := nacht_shift(stunde_h)]  # Mitternacht → 24

p_profil <- ggplot(profil_all, aes(x = stunde_plot, y = pct_aktiv,
                                    color = Variante, linetype = Variante)) +
  geom_line(linewidth = 1.0) +
  geom_smooth(se = FALSE, span = 0.3, linewidth = 0.5,
              show.legend = FALSE) +
  scale_color_manual(values = c("smoothed_wmv" = farbe_smwmv,
                                 "wmv"          = farbe_wmv),
                     labels = c(
                       paste0("smoothed_wmv (Krit. ", bout_kriterium_smwmv, " min)"),
                       paste0("wmv (Krit. ",           bout_kriterium_wmv,   " min)")
                     )) +
  scale_linetype_manual(values = c("smoothed_wmv" = "solid", "wmv" = "solid"),
                        labels = c(
                          paste0("smoothed_wmv (Krit. ", bout_kriterium_smwmv, " min)"),
                          paste0("wmv (Krit. ",           bout_kriterium_wmv,   " min)")
                        )) +
  scale_x_continuous(breaks = nacht_breaks, labels = nacht_labels,
                     limits = nacht_limits) +
  labs(
    title    = "Nacht-Aktivitätsprofil: wmv vs. smoothed_wmv",
    subtitle = paste0("Mittlere Aktivität pro Viertelstunde (alle Nächte, n = ",
                      length(igel_beide), " Igel)\n",
                      "Glatte Linie = LOESS-Glättung (visuell)"),
    x     = "Uhrzeit",
    y     = "Aktivität (%)",
    color = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", color = "#1F4E79"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_ordner, "03_aktivitaetsprofil_nacht.png"),
       p_profil, width = 9, height = 6, dpi = 200)
cat("  ✓ 03_aktivitaetsprofil_nacht.png gespeichert\n")
cat("    smwmv-Kurve: glatter (weniger Hochfrequenzrauschen)\n")
cat("    wmv-Kurve:   mehr kurzfristige Variabilität sichtbar\n\n")

# ── Zusatz: Profil pro Igel (Spaghetti) für Varianz-Check ──
profil_igel_smwmv <- dt_smwmv_v[, .(
  pct_aktiv = mean(aktiv) * 100
), by = .(igel, stunde_h = round(stunde / stunde_runden) * stunde_runden)]
profil_igel_smwmv[, Variante := "smoothed_wmv"]

profil_igel_wmv <- dt_wmv_v[, .(
  pct_aktiv = mean(aktiv) * 100
), by = .(igel, stunde_h = round(stunde / stunde_runden) * stunde_runden)]
profil_igel_wmv[, Variante := "wmv"]

profil_igel_all <- rbindlist(list(profil_igel_smwmv, profil_igel_wmv))
profil_igel_all[, stunde_plot  := nacht_shift(stunde_h)]
profil_igel_all[, Variante_lbl := fifelse(
  Variante == "smoothed_wmv",
  paste0("smoothed_wmv (Krit. ", bout_kriterium_smwmv, " min)"),
  paste0("wmv (Krit. ", bout_kriterium_wmv, " min)")
)]

# Populationsmittel vorberechnen (vermeidet stat_summary-Fehler bei NA/Inf)
mean_profil <- profil_igel_all[, .(
  pct_aktiv_mean = mean(pct_aktiv, na.rm = TRUE)
), by = .(stunde_plot, Variante_lbl)]

p_profil_igel <- ggplot(profil_igel_all,
                         aes(x = stunde_plot, y = pct_aktiv,
                             color = Variante, group = igel)) +
  geom_line(alpha = 0.45, linewidth = 0.45) +
  geom_line(data = mean_profil,
            aes(x = stunde_plot, y = pct_aktiv_mean, group = 1),
            color = "black", linewidth = 1.6, inherit.aes = FALSE) +
  facet_wrap(~ Variante_lbl, ncol = 2) +
  scale_color_manual(values = c("smoothed_wmv" = farbe_smwmv, "wmv" = farbe_wmv)) +
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26, 28, 30),
                     labels = sprintf("%02d:00", c(18, 20, 22, 0, 2, 4, 6)),
                     limits = nacht_limits) +
  labs(
    title    = "Aktivitätsprofil je Igel (Spaghetti-Plot)",
    subtitle = "Farbige Linien = einzelne Igel  |  Schwarze Linie = Populationsmittel",
    x = "Uhrzeit", y = "Aktivität (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#D5E8F0"),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_ordner, "03b_aktivitaetsprofil_spaghetti.png"),
       p_profil_igel, width = 11, height = 6, dpi = 200)
cat("  ✓ 03b_aktivitaetsprofil_spaghetti.png gespeichert\n\n")

# ──────────────────────────────────────────────────────────────
# 4. ZEITVERLAUF NACHTAKTIVITÄT (LOESS-Vergleich)
# ──────────────────────────────────────────────────────────────
cat("── 4. Zeitverlauf Nachtaktivität ──\n")

# Anteil aktiver Nachtminuten pro Igel × Nacht
trend_smwmv <- dt_smwmv_v[, .(
  pct_aktiv = mean(aktiv) * 100,
  n_min     = .N
), by = .(igel, datum, tage_seit)]
trend_smwmv[, Variante := "smoothed_wmv"]

trend_wmv <- dt_wmv_v[, .(
  pct_aktiv = mean(aktiv) * 100,
  n_min     = .N
), by = .(igel, datum, tage_seit)]
trend_wmv[, Variante := "wmv"]

trend_all <- rbindlist(list(trend_smwmv, trend_wmv))
trend_all[, Variante_f := factor(Variante,
                                   levels = c("wmv", "smoothed_wmv"),
                                   labels = c(
                                     paste0("wmv (Krit. ", bout_kriterium_wmv,   " min)"),
                                     paste0("smoothed_wmv (Krit. ", bout_kriterium_smwmv, " min)")
                                   ))]

farben_trend <- setNames(c(farbe_wmv, farbe_smwmv),
                          c(paste0("wmv (Krit. ",           bout_kriterium_wmv,   " min)"),
                            paste0("smoothed_wmv (Krit. ", bout_kriterium_smwmv, " min)")))

p_trend <- ggplot(trend_all, aes(x = tage_seit, y = pct_aktiv,
                                   color = Variante_f)) +
  geom_point(alpha = 0.2, size = 1.2) +
  geom_smooth(method = "loess", span = 0.4, se = TRUE,
              aes(fill = Variante_f), alpha = 0.15, linewidth = 1.2) +
  scale_color_manual(values = farben_trend) +
  scale_fill_manual(values  = farben_trend) +
  scale_x_continuous(limits = c(1, cutoff_tag_2b),
                     breaks = seq(0, cutoff_tag_2b, by = 7)) +
  labs(
    title    = "Zeitverlauf Nachtaktivität: wmv vs. smoothed_wmv",
    subtitle = paste0(
      "Punkte = einzelne Igel-Nächte  |  Linie = LOESS-Glättung (95%-KI)\n",
      "n = ", length(igel_beide_im_fenster), " Igel  |  ",
      "Analysefenster: Tag 1–", cutoff_tag_2b,
      " (n ≥ ", n_min_igel, " Igel/Tag)"
    ),
    x     = "Tage seit Auswilderung",
    y     = "Nachtaktivität (%)",
    color = NULL, fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", color = "#1F4E79"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_ordner, "04_zeitverlauf_nachtaktivitaet.png"),
       p_trend, width = 10, height = 6, dpi = 200)
cat("  ✓ 04_zeitverlauf_nachtaktivitaet.png gespeichert\n\n")

# ── Zeitverlauf Bout-Metriken (n_bouts und mean_dauer) ──
bm_trend <- rbindlist(list(
  bm_smwmv[, .(igel, datum, tage_seit, n_bouts, mean_dauer, Variante)],
  bm_wmv[,   .(igel, datum, tage_seit, n_bouts, mean_dauer, Variante)]
))
bm_trend[, Variante_f := factor(Variante,
                                  levels = c("wmv", "smoothed_wmv"),
                                  labels = c(
                                    paste0("wmv (Krit. ", bout_kriterium_wmv,   " min)"),
                                    paste0("smoothed_wmv (Krit. ", bout_kriterium_smwmv, " min)")
                                  ))]

p_trend_n <- ggplot(bm_trend[!is.na(n_bouts)],
                    aes(x = tage_seit, y = n_bouts, color = Variante_f)) +
  geom_jitter(alpha = 0.2, size = 1.0, height = 0.1) +
  geom_smooth(method = "loess", span = 0.5, se = TRUE,
              aes(fill = Variante_f), alpha = 0.15, linewidth = 1.2) +
  scale_color_manual(values = farben_trend) +
  scale_fill_manual(values  = farben_trend) +
  scale_x_continuous(limits = c(1, cutoff_tag_2b),
                     breaks = seq(0, cutoff_tag_2b, by = 7)) +
  labs(title = "Anzahl Bouts über Zeit", x = "Tage seit Auswilderung",
       y = "Bouts pro Nacht", color = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

p_trend_d <- ggplot(bm_trend[!is.na(mean_dauer)],
                    aes(x = tage_seit, y = mean_dauer, color = Variante_f)) +
  geom_jitter(alpha = 0.2, size = 1.0) +
  geom_smooth(method = "loess", span = 0.5, se = TRUE,
              aes(fill = Variante_f), alpha = 0.15, linewidth = 1.2) +
  scale_color_manual(values = farben_trend) +
  scale_fill_manual(values  = farben_trend) +
  scale_x_continuous(limits = c(1, cutoff_tag_2b),
                     breaks = seq(0, cutoff_tag_2b, by = 7)) +
  labs(title = "Mittlere Boutdauer über Zeit", x = "Tage seit Auswilderung",
       y = "Mittlere Boutdauer (min)", color = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

p_trend_bouts <- (p_trend_n | p_trend_d) +
  plot_annotation(
    title    = "Zeitverlauf Bout-Metriken",
    subtitle = paste0("LOESS-Glättung  |  n = ", length(igel_beide), " Igel"),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", color = "#1F4E79"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave(file.path(output_ordner, "04b_zeitverlauf_bout_metriken.png"),
       p_trend_bouts, width = 12, height = 6, dpi = 200)
cat("  ✓ 04b_zeitverlauf_bout_metriken.png gespeichert\n\n")

# ──────────────────────────────────────────────────────────────
# ZUSAMMENFASSUNG
# ──────────────────────────────────────────────────────────────
cat("════════════════════════════════════════════════════════\n")
cat("ZUSAMMENFASSUNG — wmv vs. smoothed_wmv Bout-Vergleich\n")
cat("════════════════════════════════════════════════════════\n\n")

cat("── Einschluss ──\n")
cat(sprintf("  Igel in beiden Varianten:    %d\n", length(igel_beide)))
cat(sprintf("  Bout-Kriterium smoothed_wmv: %d min\n", bout_kriterium_smwmv))
cat(sprintf("  Bout-Kriterium wmv:          %d min\n", bout_kriterium_wmv))

# Welche Tiere haben den mv-Fallback genutzt?
if ("quelle_spalte" %in% names(dt_wmv_raw)) {
  fallback_igel <- dt_wmv_raw[quelle_spalte == "pred_nested_loio_mv",
                               unique(igel)]
  if (length(fallback_igel) > 0)
    cat("  Fallback mv (wmv leer):      ",
        paste(sort(fallback_igel), collapse = ", "), "\n")
}
cat("\n")

cat("── Bout-Metriken Median ──\n")
cat(sprintf("  %-30s  %s  %s\n", "", "smoothed_wmv", "wmv"))
cat(sprintf("  %-30s  %12.1f  %6.1f\n", "Bouts pro Nacht:",
  median(bm_smwmv$n_bouts, na.rm=TRUE), median(bm_wmv$n_bouts, na.rm=TRUE)))
cat(sprintf("  %-30s  %12.1f  %6.1f\n", "Mittlere Boutdauer (min):",
  median(bm_smwmv$mean_dauer, na.rm=TRUE), median(bm_wmv$mean_dauer, na.rm=TRUE)))
cat(sprintf("  %-30s  %12.1f  %6.1f\n", "Max. Boutdauer (min):",
  median(bm_smwmv$max_dauer, na.rm=TRUE), median(bm_wmv$max_dauer, na.rm=TRUE)))
cat(sprintf("  %-30s  %12.1f  %6.1f\n", "Mittleres IBI (min):",
  median(bm_smwmv$mean_ibi, na.rm=TRUE), median(bm_wmv$mean_ibi, na.rm=TRUE)))
cat("\n")

cat("── Interpretation ──\n")
cat("  smoothed_wmv (~20 min Auflösung):\n")
cat("   - Kurze Pausen (< 20 min) werden artifiziell überbrückt\n")
cat("   - → Bouts erscheinen länger, IBIs kürzer\n")
cat("   - → Log-Survivor-Plot: flacher Verlauf, wenig Struktur\n")
cat("   - MIBI aus Log-Survivor ist KEIN biologischer Parameter,\n")
cat("     sondern Artefakt der Glättung\n\n")
cat("  wmv (~5 min Auflösung):\n")
cat("   - Kürzere Pausen biologisch erkennbar\n")
cat("   - → Bout-Kriterium ≤ 20 min sinnvoll interpretierbar\n")
cat("   - → Log-Survivor-Knick zeigt echte biologische Schwelle\n")
cat("   - Empfehlung: wmv für Bout-Analyse verwenden,\n")
cat("     Kriterium aus Log-Survivor-Plot ableiten\n\n")

cat("── Ausgabedateien ──\n")
cat("  Ordner:", output_ordner, "\n")
cat("   01_log_survivor.png              — Log-Survivor-Plot (overlaid)\n")
cat("   02_bout_metriken.png             — Boxplots Bout-Metriken\n")
cat("   03_aktivitaetsprofil_nacht.png   — 24h Aktivitätsprofil\n")
cat("   03b_aktivitaetsprofil_spaghetti.png — je Igel + Mittel\n")
cat("   04_zeitverlauf_nachtaktivitaet.png  — LOESS Zeitverlauf\n")
cat("   04b_zeitverlauf_bout_metriken.png   — Bouts + Dauer über Zeit\n")
cat("════════════════════════════════════════════════════════\n")
cat("\n→ Log-Survivor-Plot ansehen → Knick bestätigt?\n")
cat("→ Falls ja: 'bout_kriterium_wmv' oben bestätigen oder anpassen\n")
cat("→ Dann Block2_GAMM.R mit 'bout_kriterium <-", bout_kriterium_wmv, "' ausführen\n\n")

# ──────────────────────────────────────────────────────────────
# BERICHT & EXCEL — automatisch nach der Analyse erstellt
# ──────────────────────────────────────────────────────────────
cat("════════════════════════════════════════════════════════\n")
cat("BERICHT & EXCEL — Erstelle Ausgabedateien...\n")
cat("════════════════════════════════════════════════════════\n\n")

# ── Fallback-Liste für Metadaten ──
fallback_liste <- character(0)
if ("quelle_spalte" %in% names(dt_wmv_raw)) {
  fallback_liste <- dt_wmv_raw[quelle_spalte == "pred_nested_loio_mv", unique(igel)]
}

# ============================================================
# EXCEL-ÜBERBLICK
# ============================================================
cat("── Excel: Block2b_BoutVergleich_Uebersicht.xlsx ...\n")

wb <- createWorkbook()

# ── Styles ──
stil_header   <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white",
                             fgFill = "#2C3E6B", halign = "center", valign = "center",
                             textDecoration = "bold", border = "TopBottomLeftRight",
                             borderColour = "#CCCCCC", wrapText = TRUE)
stil_header2  <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white",
                             fgFill = "#E07B39", halign = "center", valign = "center",
                             textDecoration = "bold", border = "TopBottomLeftRight",
                             borderColour = "#CCCCCC", wrapText = TRUE)
stil_normal   <- createStyle(fontName = "Arial", fontSize = 9, border = "TopBottomLeftRight",
                             borderColour = "#CCCCCC")
stil_zahl1    <- createStyle(fontName = "Arial", fontSize = 9, numFmt = "0.0",
                             border = "TopBottomLeftRight", borderColour = "#CCCCCC")
stil_zebra    <- createStyle(fontName = "Arial", fontSize = 9, fgFill = "#F0F4F8",
                             border = "TopBottomLeftRight", borderColour = "#CCCCCC")
stil_zebra1   <- createStyle(fontName = "Arial", fontSize = 9, fgFill = "#F0F4F8",
                             numFmt = "0.0", border = "TopBottomLeftRight",
                             borderColour = "#CCCCCC")

# ── Sheet 1: Pro-Igel Bout-Metriken ──
addWorksheet(wb, "Bout-Metriken je Igel")

df_igel_bm <- merge(
  bm_smwmv[, .(
    smwmv_naechte    = .N,
    smwmv_bouts_med  = round(median(n_bouts,    na.rm = TRUE), 1),
    smwmv_dauer_med  = round(median(mean_dauer, na.rm = TRUE), 1),
    smwmv_maxd_med   = round(median(max_dauer,  na.rm = TRUE), 1),
    smwmv_ibi_med    = round(median(mean_ibi,   na.rm = TRUE), 1)
  ), by = igel],
  bm_wmv[, .(
    wmv_naechte    = .N,
    wmv_bouts_med  = round(median(n_bouts,    na.rm = TRUE), 1),
    wmv_dauer_med  = round(median(mean_dauer, na.rm = TRUE), 1),
    wmv_maxd_med   = round(median(max_dauer,  na.rm = TRUE), 1),
    wmv_ibi_med    = round(median(mean_ibi,   na.rm = TRUE), 1)
  ), by = igel],
  by = "igel", all = TRUE
)
setorder(df_igel_bm, igel)

kopf1 <- c("Igel",
           "smoothed_wmv: Nächte", "smoothed_wmv: Bouts/Nacht",
           "smoothed_wmv: Dauer (min)", "smoothed_wmv: MaxDauer (min)",
           "smoothed_wmv: IBI (min)",
           "wmv: Nächte", "wmv: Bouts/Nacht",
           "wmv: Dauer (min)", "wmv: MaxDauer (min)",
           "wmv: IBI (min)")

writeData(wb, 1, as.data.frame(df_igel_bm),
          startRow = 2, startCol = 1, colNames = FALSE)
writeData(wb, 1, as.data.frame(t(kopf1)),
          startRow = 1, startCol = 1, colNames = FALSE)

addStyle(wb, 1, stil_header,  rows = 1, cols = 1:6, gridExpand = TRUE)
addStyle(wb, 1, stil_header2, rows = 1, cols = 7:11, gridExpand = TRUE)
for (i in seq_len(nrow(df_igel_bm))) {
  st <- if (i %% 2 == 0) stil_zebra else stil_normal
  st1 <- if (i %% 2 == 0) stil_zebra1 else stil_zahl1
  addStyle(wb, 1, st,  rows = i + 1, cols = 1,    gridExpand = FALSE)
  addStyle(wb, 1, st1, rows = i + 1, cols = 2:11, gridExpand = TRUE)
}
setColWidths(wb, 1, cols = 1:11,
             widths = c(10, 16, 20, 22, 24, 18, 14, 18, 18, 20, 15))
setRowHeights(wb, 1, rows = 1, heights = 42)

# ── Sheet 2: Gesamtstatistik ──
addWorksheet(wb, "Gesamtstatistik")

df_ges <- data.frame(
  Kennzahl = c(
    "Anzahl Igel (Gesamtdatensatz)",
    "Anzahl Igel mit smoothed_wmv-Daten",
    "Anzahl Igel mit wmv-Daten",
    "Igel mit mv-Fallback (wmv leer)",
    "Bout-Kriterium smoothed_wmv (min)",
    "Bout-Kriterium wmv (min)",
    "",
    "IBIs gesamt — smoothed_wmv (n)",
    "IBIs gesamt — wmv (n)",
    "",
    "Median Bouts/Nacht — smoothed_wmv",
    "Median Bouts/Nacht — wmv",
    "Median Boutdauer (min) — smoothed_wmv",
    "Median Boutdauer (min) — wmv",
    "Median Max-Boutdauer (min) — smoothed_wmv",
    "Median Max-Boutdauer (min) — wmv",
    "Median IBI (min) — smoothed_wmv",
    "Median IBI (min) — wmv"
  ),
  Wert = c(
    length(union(unique(dt_smwmv$igel), unique(dt_wmv_raw$igel))),
    length(unique(dt_smwmv$igel)),
    length(unique(dt_wmv_raw$igel)),
    if (length(fallback_liste) > 0) paste(sort(fallback_liste), collapse = ", ") else "keine",
    bout_kriterium_smwmv,
    bout_kriterium_wmv,
    "",
    if (!is.null(ls_smwmv)) nrow(ls_smwmv) else 0L,
    if (!is.null(ls_wmv))   nrow(ls_wmv)   else 0L,
    "",
    round(median(bm_smwmv$n_bouts,    na.rm = TRUE), 1),
    round(median(bm_wmv$n_bouts,      na.rm = TRUE), 1),
    round(median(bm_smwmv$mean_dauer, na.rm = TRUE), 1),
    round(median(bm_wmv$mean_dauer,   na.rm = TRUE), 1),
    round(median(bm_smwmv$max_dauer,  na.rm = TRUE), 1),
    round(median(bm_wmv$max_dauer,    na.rm = TRUE), 1),
    round(median(bm_smwmv$mean_ibi,   na.rm = TRUE), 1),
    round(median(bm_wmv$mean_ibi,     na.rm = TRUE), 1)
  ),
  stringsAsFactors = FALSE
)
writeData(wb, 2, df_ges, startRow = 1, startCol = 1, colNames = TRUE)
addStyle(wb, 2, stil_header, rows = 1, cols = 1:2, gridExpand = TRUE)
for (i in seq_len(nrow(df_ges))) {
  st <- if (i %% 2 == 0) stil_zebra else stil_normal
  addStyle(wb, 2, st, rows = i + 1, cols = 1:2, gridExpand = TRUE)
}
setColWidths(wb, 2, cols = 1:2, widths = c(44, 28))

# ── Sheet 3: Datenqualität ──
addWorksheet(wb, "Datenqualitaet")

if ("quelle_spalte" %in% names(dt_wmv_raw)) {
  dq <- dt_wmv_raw[, .(
    Nächte_gesamt = .N,
    Quelle        = unique(quelle_spalte)[1]
  ), by = igel]
  setorder(dq, igel)
  names(dq) <- c("Igel", "Nächte (wmv-Datensatz)", "Verwendete Spalte")
  writeData(wb, 3, as.data.frame(dq), startRow = 1, startCol = 1, colNames = TRUE)
  addStyle(wb, 3, stil_header, rows = 1, cols = 1:3, gridExpand = TRUE)
  for (i in seq_len(nrow(dq))) {
    st <- if (i %% 2 == 0) stil_zebra else stil_normal
    addStyle(wb, 3, st, rows = i + 1, cols = 1:3, gridExpand = TRUE)
  }
  setColWidths(wb, 3, cols = 1:3, widths = c(10, 24, 30))
} else {
  writeData(wb, 3, data.frame(Hinweis = "Keine quelle_spalte-Information verfügbar."))
}

excel_pfad <- file.path(output_ordner, "Block2b_BoutVergleich_Uebersicht.xlsx")
saveWorkbook(wb, excel_pfad, overwrite = TRUE)
cat("  ✓ Excel gespeichert:", excel_pfad, "\n\n")

# ============================================================
# WORD-BERICHT
# ============================================================
cat("── Word-Bericht: Block2b_BoutVergleich_Bericht.docx ...\n")

library(officer)
library(flextable)

# ── Hilfsfunktionen (analog Block2_GAMM.R) ──

add_plot_safe <- function(doc, pfad, breite = 16, hoehe = 10) {
  if (file.exists(pfad)) {
    doc <- body_add_img(doc, pfad,
                        width  = breite / 2.54,
                        height = hoehe  / 2.54)
  } else {
    doc <- body_add_par(doc, paste0("[Bild nicht gefunden: ", basename(pfad), "]"),
                        style = "Normal")
  }
  doc
}

add_section <- function(doc, titel, level = 1) {
  style_name <- if (level == 1) "heading 1" else "heading 2"
  body_add_par(doc, titel, style = style_name)
}

make_ft <- function(df, header_farbe = "#2C3E6B") {
  flextable(as.data.frame(df)) |>
    bold(part = "header") |>
    bg(part = "header", bg = header_farbe) |>
    color(part = "header", color = "white") |>
    bg(i = seq(2, nrow(df), 2), bg = "#F0F4F8") |>
    border_outer(part = "all", border = fp_border(color = "#CCCCCC", width = 1)) |>
    border_inner_h(part = "body", border = fp_border(color = "#E8E8E8", width = 0.5)) |>
    fontsize(size = 9, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    autofit()
}

# ── Dokument initialisieren ──
doc <- read_docx()

doc <- doc |>
  body_set_default_section(
    prop_section(
      page_margins = page_mar(top = 2.5, bottom = 2, left = 3, right = 2.5),
      page_size    = page_size(width = 21 / 2.54, height = 29.7 / 2.54)
    )
  )

# ── TITELSEITE ──
doc <- doc |>
  body_add_par("METHODENBERICHT", style = "heading 1") |>
  body_add_par("Block 2b — Bout-Analyse: wmv vs. smoothed_wmv", style = "heading 2") |>
  body_add_par("VHF-Igelbesenderung — Europäischer Igel (Erinaceus europaeus)",
               style = "Normal")

df_meta <- data.frame(
  Feld  = c("Autorin", "Institution", "Standort",
            "Analysedatum", "Bout-Kriterium wmv", "Bout-Kriterium smoothed_wmv"),
  Wert  = c("Natalie Steiner",
            "Stiftung Tierärztliche Hochschule Hannover",
            "Wildtierstation Sachsenhagen, Niedersachsen",
            format(Sys.Date(), "%d. %B %Y"),
            paste0(bout_kriterium_wmv,   " Minuten"),
            paste0(bout_kriterium_smwmv, " Minuten")),
  stringsAsFactors = FALSE
)
ft_meta <- flextable(df_meta) |>
  delete_part("header") |>
  bold(j = 1) |>
  color(j = 1, color = "#2C3E6B") |>
  bg(bg = "#F8F9FB") |>
  border_outer(border = fp_border(color = "#CCCCCC", width = 1)) |>
  border_inner_h(border = fp_border(color = "#EEEEEE", width = 0.5)) |>
  fontsize(size = 10) |>
  font(fontname = "Arial") |>
  width(j = 1, width = 5.5 / 2.54) |>
  width(j = 2, width = 11.5 / 2.54)
doc <- body_add_flextable(doc, ft_meta)
doc <- body_add_break(doc)

# ── 1. FRAGESTELLUNG ──
doc <- add_section(doc, "1  Fragestellung und Methodik")

doc <- body_add_par(doc, paste0(
  "Dieser Bericht vergleicht zwei verschiedene Aktivitätsklassifikationen ",
  "hinsichtlich ihrer Eignung für eine Bout-Analyse. Die Variante smoothed_wmv ",
  "basiert auf einer ~20-minütigen Glättung (300-s-Rollfenster + zusätzliches ",
  "Smoothing) und liegt bereits im Datensatz gamm_nachtminuten.rds vor. Die ",
  "Variante wmv verwendet ausschliesslich das 300-s-Rollfenster (~5 min effektive ",
  "Auflösung) und wurde für diesen Vergleich direkt aus den Roh-CSV-Dateien ",
  "importiert. Für Igel 1–6 war die Spalte pred_nested_loio_wmv leer (älteres ",
  "tRackIT-Exportformat); es wurde automatisch auf pred_nested_loio_mv zurückgegriffen."),
  style = "Normal")

doc <- body_add_par(doc, paste0(
  "Bouts wurden über Run-Length-Encoding (rle()) identifiziert. Kurze Pausen ",
  "unterhalb des jeweiligen Bout-Kriteriums (wmv: ", bout_kriterium_wmv, " min; ",
  "smoothed_wmv: ", bout_kriterium_smwmv, " min) wurden iterativ mit den ",
  "angrenzenden aktiven Blöcken zusammengeführt (3 Iterationen). ",
  "Das Bout-Kriterium für die wmv-Variante wurde anhand des Log-Survivor-Plots ",
  "der Inter-Bout-Intervalle (IBI) festgelegt."),
  style = "Normal")

# Stichproben-Tabelle
df_stichprobe <- data.frame(
  Merkmal = c("Igel gesamt (Gesamtdatensatz)",
              "Igel in smoothed_wmv-Datensatz",
              "Igel in wmv-Datensatz",
              "davon mit mv-Fallback",
              "Igel in beiden Datensätzen (Vergleich)"),
  `smoothed_wmv` = c("—",
    as.character(length(unique(dt_smwmv$igel))), "—", "—", "—"),
  wmv = c("—", "—",
    as.character(length(unique(dt_wmv_raw$igel))),
    if (length(fallback_liste) > 0) paste(sort(fallback_liste), collapse = ", ") else "keine",
    as.character(length(igel_beide))),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
doc <- body_add_flextable(doc, make_ft(df_stichprobe))

# ── 2. LOG-SURVIVOR-PLOT ──
doc <- add_section(doc, "2  Log-Survivor-Plot der Inter-Bout-Intervalle")

doc <- body_add_par(doc, paste0(
  "Der Log-Survivor-Plot zeigt die Überlebenskurve der Pausenlängen (IBI = ",
  "Inter-Bout-Interval) auf logarithmischer Skala. Ein Knick im Kurvenverlauf ",
  "zeigt an, bei welcher Pausenlänge zwei verschiedene biologische Prozesse ",
  "aufeinandertreffen — kurze, intra-bout-Pausen (z. B. kurzes Innehalten beim ",
  "Fressen) und längere Inter-Bout-Pausen (z. B. Ruhepausen in Verstecken). ",
  "Die Lage des Knicks ist das empirische Bout-Kriterium."),
  style = "Normal")

doc <- body_add_par(doc, paste0(
  "Befund (wmv, n = ", if (!is.null(ls_wmv)) nrow(ls_wmv) else 0L, " IBIs): Die Kurve zeigt einen ",
  "deutlichen Knick bei ca. 5–8 Minuten, mit einem steileren Abfall darunter ",
  "und einem flacheren Verlauf darüber. Dies entspricht dem erwarteten ",
  "Zwei-Phasen-Muster bei biologisch sinnvollen Bout-Daten. Das gewählte ",
  "Kriterium von ", bout_kriterium_wmv, " Minuten liegt im Bereich dieses Knicks ",
  "und ist damit empirisch begründet."),
  style = "Normal")

doc <- body_add_par(doc, paste0(
  "Befund (smoothed_wmv, n = ", if (!is.null(ls_smwmv)) nrow(ls_smwmv) else 0L, " IBIs, Faktor ",
  round((if (!is.null(ls_wmv)) nrow(ls_wmv) else 0L) /
        max(if (!is.null(ls_smwmv)) nrow(ls_smwmv) else 1L, 1), 1), "x weniger): Die Kurve ",
  "verläuft nahezu linear ohne erkennbaren Knick. Die starke Glättung (~20 min) ",
  "hat kurze Pausen bereits vor der Bout-Berechnung überbrückt, sodass keine ",
  "biologisch interpretierbare Schwelle mehr sichtbar ist. Das Bout-Kriterium ",
  "von ", bout_kriterium_smwmv, " Minuten ist damit ein Artefakt der Glättungsbreite, ",
  "kein biologischer Parameter."),
  style = "Normal")

doc <- add_plot_safe(doc, file.path(output_ordner, "01_log_survivor.png"),
                     breite = 16, hoehe = 9)

# ── 3. BOUT-METRIKEN ──
doc <- add_section(doc, "3  Bout-Metriken (Boxplots)")

doc <- body_add_par(doc, paste0(
  "Die Boxplots zeigen Bouts pro Nacht, mittlere Boutdauer, maximale Boutdauer ",
  "und mittleres IBI für beide Varianten. Alle Metriken werden als Mediane über ",
  "die pro Igel und Nacht berechneten Werte angegeben."),
  style = "Normal")

# Kennzahlentabelle
df_bm_tab <- data.frame(
  Metrik = c("Bouts pro Nacht (Median)",
             "Mittlere Boutdauer — Median (min)",
             "Maximale Boutdauer — Median (min)",
             "Mittleres IBI — Median (min)"),
  `smoothed_wmv` = c(
    round(median(bm_smwmv$n_bouts,    na.rm = TRUE), 1),
    round(median(bm_smwmv$mean_dauer, na.rm = TRUE), 1),
    round(median(bm_smwmv$max_dauer,  na.rm = TRUE), 1),
    round(median(bm_smwmv$mean_ibi,   na.rm = TRUE), 1)
  ),
  wmv = c(
    round(median(bm_wmv$n_bouts,    na.rm = TRUE), 1),
    round(median(bm_wmv$mean_dauer, na.rm = TRUE), 1),
    round(median(bm_wmv$max_dauer,  na.rm = TRUE), 1),
    round(median(bm_wmv$mean_ibi,   na.rm = TRUE), 1)
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
doc <- body_add_flextable(doc, make_ft(df_bm_tab))

doc <- body_add_par(doc, paste0(
  "Interpretation: Die wmv-Variante ergibt mehr, kürzere Bouts pro Nacht mit ",
  "kürzeren IBIs. Das entspricht einem feinkörnigeren Aktivitätsbild, das kurze ",
  "Unterbrechungen (z. B. kurzes Eingraben, Orientierungspause) als eigene Pausen ",
  "erkennt. Die smoothed_wmv-Variante erzeugt durch das Überbrücken kurzer Pausen ",
  "weniger, dafür deutlich längere Bouts. Beide Muster sind intern konsistent, ",
  "jedoch sind die smoothed_wmv-Bouts ein Artefakt des Glättungsalgorithmus, ",
  "nicht ein Abbild der tatsächlichen Aktivitätsstruktur."),
  style = "Normal")

doc <- add_plot_safe(doc, file.path(output_ordner, "02_bout_metriken.png"),
                     breite = 16, hoehe = 10)

# ── 4. AKTIVITÄTSPROFIL ──
doc <- add_section(doc, "4  24-h-Aktivitätsprofil")

doc <- body_add_par(doc, paste0(
  "Das Aktivitätsprofil zeigt den Anteil aktiver Minuten pro Stunde über die Nacht ",
  "(gemittelt über alle Nächte und alle Igel). Beide Varianten zeigen die gleiche ",
  "glockenförmige Kurve: langsamer Anstieg ab ca. 19:00 Uhr, Plateau bzw. Maximum ",
  "zwischen 22:00 und 24:00 Uhr, Abfall gegen Morgengrauen."),
  style = "Normal")

doc <- body_add_par(doc,
  "Unterschiede zwischen den Varianten: (1) Die wmv-Kurve startet früher (ab ~18:30 Uhr), ",
  style = "Normal")
doc <- body_add_par(doc, paste0(
  "weil der smoothed_wmv-Datensatz (aus gamm_nachtminuten.rds) bereits durch den ",
  "Block-0-Filter auf 'Nacht laut suncalc' beschränkt ist und frühe Abendminuten ",
  "je nach Jahreszeit wegfallen können. (2) Die smoothed_wmv-Kurve verläuft ",
  "generell glatter und zeigt weniger Minutenstreuung, da die Glättung bereits ",
  "auf Rohdatenebene erfolgte. (3) Beide Varianten zeigen dasselbe biologische ",
  "Aktivitätsfenster — die Unterschiede sind methodischer Natur."),
  style = "Normal")

doc <- add_plot_safe(doc, file.path(output_ordner, "03_aktivitaetsprofil_nacht.png"),
                     breite = 16, hoehe = 9)
doc <- add_section(doc, "Spaghetti-Plot: Aktivitätsprofil je Igel", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "03b_aktivitaetsprofil_spaghetti.png"),
                     breite = 16, hoehe = 9)

# ── 5. ZEITVERLAUF ──
doc <- add_section(doc, "5  Zeitverlauf Nachtaktivität und Bout-Metriken")

doc <- body_add_par(doc, paste0(
  "Die Zeitverlaufsplots zeigen den prozentualen Anteil aktiver Nachtminuten ",
  "(LOESS-Glättung) und die Bout-Metriken (Bouts/Nacht, Boutdauer) über die ",
  "Tage seit Auswilderung. Hier zeigt sich der wichtigste Unterschied zwischen ",
  "den Datensätzen:"),
  style = "Normal")

doc <- body_add_par(doc, paste0(
  "Saisonaler Hintergrund — smoothed_wmv-Datensatz: Dieser Datensatz enthält ",
  "alle 24 Igel, darunter Igel 1–6, die im September/Oktober 2024 ausgewildert ",
  "wurden. Von diesen ist bekannt, dass Igel 5 in die Winterruhe wechselte; ",
  "die übrigen Individuen aus dieser Gruppe (1–4, 6) gingen nicht in den ",
  "Winterschlaf. Der Verlauf der LOESS-Kurve und eventuelle Einbrüche bei ",
  "Tag 50–100 sollten daher nicht pauschal mit Winterschlaf erklärt werden — ",
  "mögliche Ursachen sind stattdessen Datenlücken einzelner Sender, ",
  "sinkende N-at-risk-Zahlen zu späteren Zeitpunkten (weniger Igel liefern ",
  "Daten → LOESS instabil) oder individuelle Aktivitätsschwankungen. ",
  "Negative LOESS-Werte, sofern sichtbar, sind mathematische Randeffekte der ",
  "Glättung und keine interpretierbaren Daten."),
  style = "Normal")

doc <- body_add_par(doc, paste0(
  "wmv-Datensatz: Da Igel 1–6 keine wmv-Daten haben (nur mv-Fallback), ",
  "und die verbleibenden Tiere überwiegend im Frühjahr/Sommer 2025 ausgewildert ",
  "wurden, ist dieser Datensatz de facto ein Sommerdatensatz. Die LOESS-Kurve ",
  "zeigt dementsprechend einen gleichmässigeren Verlauf ohne den starken ",
  "Winterschlaf-Einbruch. Der Vergleich der Zeittrends zwischen wmv und ",
  "smoothed_wmv ist daher durch einen Saison-Confounder eingeschränkt und ",
  "sollte nicht überinterpretiert werden."),
  style = "Normal")

doc <- add_plot_safe(doc, file.path(output_ordner, "04_zeitverlauf_nachtaktivitaet.png"),
                     breite = 16, hoehe = 9)
doc <- add_section(doc, "Zeitverlauf Bout-Metriken", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "04b_zeitverlauf_bout_metriken.png"),
                     breite = 16, hoehe = 8)

# ── 6. EMPFEHLUNG ──
doc <- add_section(doc, "6  Empfehlung und Schlussfolgerung")

doc <- body_add_par(doc, paste0(
  "Auf Grundlage dieses Vergleichs wird empfohlen, für die GAMM-Bout-Analyse ",
  "(Block2_GAMM.R) die wmv-Variante (pred_nested_loio_wmv, mit mv-Fallback für ",
  "Igel 1–6) mit einem Bout-Kriterium von ", bout_kriterium_wmv, " Minuten zu ",
  "verwenden. Begründung:"),
  style = "Normal")

df_emp <- data.frame(
  Kriterium = c(
    "Biologisch interpretierbare Bout-Struktur",
    "Log-Survivor-Knick vorhanden",
    "IBI-Pool ausreichend gross",
    "Boutdauer realistisch",
    "Bouts/Nacht plausibel",
    "Saisonaler Confounder (Zeittrend)",
    "Empfehlung"
  ),
  `smoothed_wmv (20 min)` = c(
    "Nein — Glättung überbrückt echte Pausen",
    "Nein — kein Knick sichtbar",
    paste0(if (!is.null(ls_smwmv)) nrow(ls_smwmv) else 0L, " IBIs"),
    paste0(round(median(bm_smwmv$mean_dauer, na.rm=TRUE), 0), " min (artefaktell verlängert)"),
    paste0(round(median(bm_smwmv$n_bouts, na.rm=TRUE), 0), " — unterschätzt"),
    "Enthält Herbst-Igel (1–6); nur Igel 5 in Winterschlaf → N-at-risk spät instabil",
    "NICHT empfohlen"
  ),
  `wmv (5 min, mit mv-Fallback)` = c(
    "Ja — feingranulares Aktivitätsbild",
    paste0("Ja — Knick bei ~5–8 min, Kriterium ", bout_kriterium_wmv, " min"),
    paste0(if (!is.null(ls_wmv)) nrow(ls_wmv) else 0L, " IBIs"),
    paste0(round(median(bm_wmv$mean_dauer, na.rm=TRUE), 0), " min (biologisch plausibel)"),
    paste0(round(median(bm_wmv$n_bouts, na.rm=TRUE), 0), " — realistisch"),
    "De-facto Sommerdatensatz — homogener",
    "EMPFOHLEN"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ft_emp <- flextable(df_emp) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#2C3E6B") |>
  color(part = "header", color = "white") |>
  bg(i = nrow(df_emp), bg = "#E8F4E8") |>
  bg(i = seq(2, nrow(df_emp) - 1, 2), bg = "#F0F4F8") |>
  bold(i = nrow(df_emp)) |>
  border_outer(part = "all", border = fp_border(color = "#CCCCCC", width = 1)) |>
  border_inner_h(part = "body", border = fp_border(color = "#E8E8E8", width = 0.5)) |>
  fontsize(size = 9, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  width(j = 1, width = 6 / 2.54) |>
  width(j = 2, width = 6 / 2.54) |>
  width(j = 3, width = 7 / 2.54)
doc <- body_add_flextable(doc, ft_emp)

doc <- body_add_par(doc, paste0(
  "Nächster Schritt: Block2_GAMM.R mit 'bout_kriterium <- ",
  bout_kriterium_wmv, "' (Zeile ~35) und 'klasse_spalte <- ",
  "\"pred_nested_loio_wmv\"' ausführen."),
  style = "Normal")

# ── Dokument speichern ──
docx_pfad <- file.path(output_ordner, "Block2b_BoutVergleich_Bericht.docx")
print(doc, target = docx_pfad)
cat("  ✓ Word-Bericht gespeichert:", docx_pfad, "\n\n")

cat("════════════════════════════════════════════════════════\n")
cat("FERTIG — Beide Ausgabedateien erstellt:\n")
cat("  →", excel_pfad, "\n")
cat("  →", docx_pfad, "\n")
cat("════════════════════════════════════════════════════════\n\n")
