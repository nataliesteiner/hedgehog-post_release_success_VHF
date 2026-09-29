# ==============================================================
# Hierarchical additive mixed model (GAMM) — hedgehog activity
# ==============================================================
# Project:  Hedgehog VHF telemetry, Lower Saxony
# Author:   Natalie Steiner
# Question: Does night activity change over time after release?
#           Are there individual differences?
#
# Method:   GAMM via mgcv::bam()
#   Option A: daily level  — proportion of active night minutes ~ s(days_since)
#   Option B: minute level — active/passive ~ s(days_since) x s(hour)
#   Bout:     per-night bout metrics ~ s(days_since)
#
# Note: Block0_Datenpipeline.R does NOT need to be re-run.
#   Data are loaded directly from data/activity/*.csv (wmv pipeline).
#   Rationale: Block2b showed that smoothed_wmv (~20 min resolution)
#   produces artefactual bouts; wmv (~5 min) + 10 min criterion is preferred.
# ==============================================================
# IMPORTANT: always run the script from line 1 (Ctrl+Alt+R).
# ==============================================================

pakete <- c("mgcv", "data.table", "ggplot2", "patchwork", "scales",
            "officer", "flextable", "openxlsx",
            "lubridate", "suncalc", "readxl")
neu    <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) {
  cat("Installiere fehlende Pakete:", paste(neu, collapse = ", "), "\n")
  install.packages(neu)
}

library(mgcv)
library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(openxlsx)
library(lubridate)
library(suncalc)
library(readxl)

# NULL-coalescing helper (not in base R)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# Modalwert-Funktion (wie Block0/Block2b)
modal_val_chr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1L]
}

# Datumsparser (wie Block0/Block2b — toleriert TT.MM.JJJJ und JJJJ-MM-TT)
parse_datum <- function(x) {
  x      <- as.character(x)
  result <- suppressWarnings(as.Date(x, format = "%d.%m.%Y"))
  na_idx <- is.na(result)
  result[na_idx] <- suppressWarnings(as.Date(x[na_idx]))
  result
}

# gratia für schöne GAMM-Plots (optional aber empfohlen)
gratia_verfuegbar <- requireNamespace("gratia", quietly = TRUE)
if (gratia_verfuegbar) library(gratia)

cat("✓ Alle Pakete geladen\n\n")

# ──────────────────────────────────────────────────────────────
# EINSTELLUNGEN
# ──────────────────────────────────────────────────────────────

# Projektwurzel — einzige Zeile die du ggf. anpassen musst
projekt_root  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"

output_ordner <- file.path(projekt_root, "output", "Block2_GAMM")
dir.create(output_ordner, showWarnings = FALSE, recursive = TRUE)

# ── wmv-Datenquelle (CSV-Pipeline, wie in Block2b_BoutVergleich.R) ──────────
# HINWEIS: Block2b hat gezeigt, dass smoothed_wmv (~20 Min Auflösung) keine
# biologisch sinnvolle Bout-Struktur ergibt (kein Log-Survivor-Knick).
# Empfehlung: wmv (~5 Min) mit Bout-Kriterium 10 Min verwenden.
# Diese Pipeline ersetzt das direkte Laden aus gamm_nachtminuten.rds
# (das smoothed_wmv enthält) und lädt wmv/mv direkt aus den CSVs.
# Block0_Datenpipeline.R bleibt UNVERÄNDERT.
daten_ordner       <- file.path(projekt_root, "data", "activity")
meta_datei         <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")

# Standort (Niedersachsen) — für suncalc Tag/Nacht-Bestimmung
standort_lat       <- 52.39729710523643
standort_lon       <-  9.216876248766871
zeitzone           <- "Europe/Berlin"

# Schwellenwert Toter-Sender-Erkennung (wie Block0)
schwelle_aktiv_pct <- 5   # % Tagesaktivität unterhalb → statisches Signal

# ── Ausschlusskriterien ──────────────────────────────────────
#
# WICHTIG: Ein prozentualer Nachtdaten-Schwellenwert (z.B. < 30%)
# würde Sommer-Auswilderungen systematisch benachteiligen, weil die
# Nacht im Juli (~52°N) biologisch nur ~8h lang ist → max. ~33% der
# Gesamtminuten fallen in die Nacht. Stattdessen: absoluter Schwellenwert
# für Nachtminuten (unabhängig von Tageslänge).
#
min_naechte      <- 3     # Mindestanzahl Nächte mit Nachtdaten
min_nacht_min    <- 60    # Mindestanzahl absoluter Nachtminuten gesamt
                          # (60 Min = grobe Mindestbasis; bei Bedarf anpassen)
min_pct_aktiv    <- 1     # Igel mit < 1% Nachtaktivität → statisches Signal

# ── Visualisierungs-Cutoff (nur für Plots, NICHT fürs Modell) ────────────────
# Ab Tag 15 sinkt N auf ≤10 und Igel22+23 dominieren die Schätzung.
# Modelle laufen weiterhin auf allen Daten; nur die Zeitachse in den Grafiken
# wird auf diesen Wert begrenzt bzw. der Bereich danach grau markiert.
heatmap_vis_cutoff <- 14

farben <- c("Nacht" = "#2C3E6B", "Tag" = "#E8A87C")

# ──────────────────────────────────────────────────────────────
# DATEN LADEN — wmv-Pipeline (ersetzt readRDS gamm_nachtminuten.rds)
# ──────────────────────────────────────────────────────────────
# Empfehlung aus Block2b_BoutVergleich: wmv statt smoothed_wmv
# verwenden, da smoothed_wmv artifizielle Bout-Überbrückung erzeugt
# (~20 Min Auflösung → kein Log-Survivor-Knick sichtbar).
# Klassifikationsspalte: pred_nested_loio_wmv (Fallback: pred_nested_loio_mv)
# Bout-Kriterium: 10 Minuten (unverändert, war schon korrekt).
# ──────────────────────────────────────────────────────────────

cat("Lade wmv-Daten aus CSVs (statt smoothed_wmv aus gamm_nachtminuten.rds)...\n")

# ── Metadaten laden (Hard-Release-Datum je Igel) ──────────────
meta_pfad <- path.expand(meta_datei)
if (!file.exists(meta_pfad)) {
  stop("Metadatei nicht gefunden: ", meta_pfad, "\n",
       "Bitte meta_datei-Pfad in den EINSTELLUNGEN prüfen.")
}
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

# ── CSV-Dateien lokalisieren ──────────────────────────────────
alle_csvs <- list.files(path.expand(daten_ordner),
                        pattern = "classification_active_passive.*\\.csv$",
                        full.names = TRUE)
cat("  Gefundene CSVs:", length(alle_csvs), "\n")

if (length(alle_csvs) == 0L) {
  stop("Keine classification_active_passive*.csv gefunden in:\n  ", daten_ordner,
       "\nBitte daten_ordner in den EINSTELLUNGEN prüfen.")
}

# Hilfsfunktion: CSV-Pfad für einen Igel finden
finde_csv <- function(igel_name) {
  nr      <- gsub("[^0-9]", "", igel_name)
  treffer <- grep(paste0("Igel[[:space:]_]?0*", nr, "[^0-9]"),
                  alle_csvs, value = TRUE)
  if (length(treffer) == 0L) return(NA_character_)
  treffer[1L]
}
meta[, csv_pfad := sapply(igel, finde_csv)]

# ── CSV-Pipeline: alle Minuten laden (Tag + Nacht) für dt_tage ──
# Nachtminuten werden separat gefiltert für dt_nacht (Modell B)
wmv_alle_liste  <- list()   # alle Minuten (für dt_tage-Aggregation)
wmv_nacht_liste <- list()   # nur Nachtminuten (dt_nacht für Modell B)

for (i in seq_len(nrow(meta))) {

  igel_info <- meta[i]
  igel_name <- igel_info$igel
  csv_pfad  <- igel_info$csv_pfad

  if (is.na(csv_pfad)) {
    cat("  ⚠", igel_name, "— keine CSV gefunden, übersprungen\n")
    next
  }

  hard_release_i <- igel_info$hard_release

  # 1. CSV einlesen
  dt <- tryCatch(
    fread(csv_pfad, na.strings = c("", "NA"), showProgress = FALSE),
    error = function(e) {
      cat("  ✗", igel_name, "Lesefehler:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0L) next

  # Spalten tolerant umbenennen
  setnames(dt,
    old = c("Station Name", "Individual Name", "_time", "0"),
    new = c("station", "igel_col", "time_raw", "signal"),
    skip_absent = TRUE)

  # 2. Klassifikationsspalte wählen
  #    Priorität: wmv → mv → loio (reines LOIO ohne Zeitgewichtung)
  #    Gilt für alle Igel — nicht nur für Igel 1–6.
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

  # 4. Klassifikation binär codieren
  klasse_raw <- dt[[klasse_spalte_i]]
  dt[, klasse_wmv := fcase(
    klasse_raw == "a", "aktiv",
    klasse_raw == "p", "passiv",
    default = NA_character_
  )]

  # 5. Minutenreduktion: Modalwert pro Minute
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
        "% Aktivität, übersprungen\n")
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
  dt_min[, igel      := igel_name]
  dt_min[, aktiv     := as.integer(klasse_modal == "aktiv")]

  # 9a. Alle Minuten sammeln (für dt_tage-Aggregation)
  wmv_alle_liste[[igel_name]] <- dt_min[!is.na(klasse_modal),
    .(igel, datum, tage_seit, tageszeit, stunde, aktiv)]

  # 9b. Nur Nachtminuten (für dt_nacht / Modell B)
  nacht_rows <- dt_min[tageszeit == "Nacht" & !is.na(klasse_modal),
    .(igel, datum, tage_seit, stunde, aktiv)]
  if (nrow(nacht_rows) == 0L) {
    cat("  ⚠", igel_name, "— keine Nachtminuten nach Filterung\n")
    next
  }
  wmv_nacht_liste[[igel_name]] <- nacht_rows
  cat("  ✓", igel_name, "—", format(nrow(nacht_rows), big.mark = "'"),
      "Nachtminuten [", klasse_spalte_i, "]\n")
}

# ── Zusammenführen ────────────────────────────────────────────
wmv_alle  <- rbindlist(wmv_alle_liste)
dt_nacht  <- rbindlist(wmv_nacht_liste)
dt_nacht[, igel := factor(igel)]

cat("\nwmv gesamt:", format(nrow(dt_nacht), big.mark = "'"),
    "Nachtminuten | Igel:", uniqueN(dt_nacht$igel), "\n")

# ── dt_tage aufbauen (Tagesaggregation, wie Block0-Struktur) ──
# Benötigt für Ausschlussprüfung (Zeilen tageszeit=="Nacht") und Modell A.
# pct_aktiv ist ein Anteil 0–1 (wie in gamm_tagesdaten.rds).
dt_tage <- wmv_alle[, .(
  n_min     = .N,
  pct_aktiv = mean(aktiv, na.rm = TRUE)
), by = .(igel, datum, tageszeit, tage_seit)]
dt_tage[, igel := factor(igel)]

cat("Tagesdaten:    ", nrow(dt_tage), "Zeilen (Tag+Nacht aggregiert)\n")
cat("Nachtminuten:  ", format(nrow(dt_nacht), big.mark = "'"), "Zeilen\n")
cat("Igel gesamt:   ", nlevels(dt_tage$igel), "\n\n")

# ──────────────────────────────────────────────────────────────
# AUSSCHLUSSKRITERIEN
# ──────────────────────────────────────────────────────────────

cat("── Ausschlussprüfung ──\n")

# Pro Igel: Anzahl Nächte, absolute Nachtminuten, mittlere Aktivität nachts
igel_qc <- dt_tage[, .(
  n_naechte      = uniqueN(datum[tageszeit == "Nacht"]),
  n_nacht_min    = sum(n_min[tageszeit == "Nacht"]),          # absolut, nicht %
  pct_nacht_dat  = sum(n_min[tageszeit == "Nacht"]) /
                   sum(n_min) * 100,                           # nur zur Info
  mean_aktiv_n   = mean(pct_aktiv[tageszeit == "Nacht"] * 100, na.rm = TRUE)
), by = igel]

# Ausschlussentscheidung — KEIN prozentualer Nachtdaten-Schwellenwert,
# da dieser Sommer-Auswilderungen systematisch benachteiligt
igel_qc[, ausschluss_grund := fcase(
  n_naechte  < min_naechte,  paste0("< ", min_naechte, " Nächte (n=", n_naechte, ")"),
  n_nacht_min < min_nacht_min, paste0("< ", min_nacht_min, " Nachtminuten gesamt (n=",
                                       n_nacht_min, ")"),
  mean_aktiv_n < min_pct_aktiv, "< 1% Nachtaktivität (statisches Signal / toter Sender)",
  default = NA_character_
)]

igel_ausschluss <- igel_qc[!is.na(ausschluss_grund)]
igel_einschluss <- igel_qc[ is.na(ausschluss_grund)]

cat("\nAusgeschlossen (", nrow(igel_ausschluss), "Igel):\n")
if (nrow(igel_ausschluss) > 0)
  print(igel_ausschluss[, .(igel, n_naechte, pct_nacht_dat = round(pct_nacht_dat,1),
                              mean_aktiv_n = round(mean_aktiv_n,1), ausschluss_grund)])

cat("\nEingeschlossen (", nrow(igel_einschluss), "Igel):\n")
print(igel_einschluss[, .(igel, n_naechte, pct_nacht_dat = round(pct_nacht_dat,1),
                            mean_aktiv_n = round(mean_aktiv_n,1))])

# Datensätze filtern
igel_ok <- igel_einschluss$igel
dt_tage_ok  <- dt_tage[ igel %in% igel_ok]
dt_nacht_ok <- dt_nacht[igel %in% igel_ok]

# Faktor neu setzen (nur eingeschlossene Igel)
dt_tage_ok[,  igel := droplevels(igel)]
dt_nacht_ok[, igel := droplevels(igel)]

cat("\n→", nlevels(dt_tage_ok$igel), "Igel für GAMM-Analyse\n\n")

# ──────────────────────────────────────────────────────────────
# DYNAMISCHER ZEITCUTOFF — N-AT-RISK-SCHWELLE
# ──────────────────────────────────────────────────────────────
# Problem: Späte Zeitpunkte werden von nur wenigen Langzeit-
# beobachteten dominiert (Survivorship Bias). Der "Populations-
# trend" ab Tag X basiert evtl. nur auf 2–3 Individuen.
#
# Lösung: Alle Modelle laufen nur bis zu dem Tag, an dem noch
# mindestens n_min_igel Igel Daten liefern. Danach ist kein
# sinnvoller Populationsschluss mehr möglich.
# ──────────────────────────────────────────────────────────────

n_min_igel <- 5   # Mindestanzahl Igel pro Zeitpunkt
                  # → anpassen falls nötig (bei kleiner Stichprobe evtl. 4)
                  # Hinweis: Ab Tag ~15 sinkt N auf ≤10 (Survivorship Bias).
                  # Die Modelle nutzen aber alle verfügbaren Daten für
                  # stabilere Schätzungen; in der Heatmap-Visualisierung
                  # wird der Bereich >14 Tage separat kommentiert.

# Pro Tag: wie viele eingeschlossene Igel haben Nachtdaten?
n_pro_tag <- dt_tage_ok[tageszeit == "Nacht" & !is.na(pct_aktiv),
                         .(n_igel = uniqueN(igel)), by = tage_seit]

# Letzter Tag mit N ≥ Schwelle
tage_ausreichend <- n_pro_tag[n_igel >= n_min_igel, tage_seit]

if (length(tage_ausreichend) == 0) {
  stop(paste0("Kein einziger Tag hat ≥ ", n_min_igel,
              " Igel mit Daten. Bitte n_min_igel verringern."))
}

cutoff_tag <- max(tage_ausreichend)

cat("══════════════════════════════════════════════\n")
cat("DYNAMISCHER ZEITCUTOFF (N-at-risk-Schwelle)\n")
cat("══════════════════════════════════════════════\n")
cat("Mindest-N:          ", n_min_igel, "Igel pro Zeitpunkt\n")
cat("Gesamter Zeitraum:  ", dt_tage_ok[, min(tage_seit)], "–",
                           dt_tage_ok[, max(tage_seit)], "Tage\n")
cat("Cutoff bei Tag:     ", cutoff_tag, "\n")
cat("Ausgeschlossene Tage:", dt_tage_ok[tage_seit > cutoff_tag, uniqueN(tage_seit)],
    "(zu wenig Igel)\n")

# Welche Igel fallen nach dem Cutoff raus (nur Daten nach cutoff_tag)?
igel_nur_spaet <- dt_tage_ok[tage_seit > cutoff_tag,
                               .(max_tag = max(tage_seit)), by = igel][
                 !igel %in% dt_tage_ok[tage_seit <= cutoff_tag]$igel]

if (nrow(igel_nur_spaet) > 0) {
  cat("Hinweis: Igel ohne Daten ≤ cutoff_tag (tragen nicht zu Modellen bei):\n")
  print(igel_nur_spaet)
}
cat("══════════════════════════════════════════════\n\n")

# N-at-risk-Tabelle ausgeben
cat("N-at-risk um den Cutoff:\n")
print(n_pro_tag[tage_seit >= (cutoff_tag - 5) & tage_seit <= (cutoff_tag + 5)])
cat("\n")

# Alle Modelldaten auf Zeitfenster ≤ cutoff_tag beschränken
dt_tage_ok  <- dt_tage_ok[ tage_seit <= cutoff_tag]
dt_nacht_ok <- dt_nacht_ok[tage_seit <= cutoff_tag]

# Faktor-Levels bereinigen (Igel ohne Daten im Fenster entfernen)
igel_im_fenster <- dt_tage_ok[!is.na(pct_aktiv), unique(igel)]
dt_tage_ok  <- dt_tage_ok[ igel %in% igel_im_fenster]
dt_nacht_ok <- dt_nacht_ok[igel %in% igel_im_fenster]
dt_tage_ok[,  igel := droplevels(factor(igel))]
dt_nacht_ok[, igel := droplevels(factor(igel))]

cat("Effektives Analysefenster: Tag", dt_tage_ok[, min(tage_seit)],
    "–", dt_tage_ok[, max(tage_seit)], "\n")
cat("Igel im Analysefenster:   ", nlevels(dt_tage_ok$igel), "\n\n")

# ──────────────────────────────────────────────────────────────
# BEOBACHTUNGSABDECKUNG — WER WAR WANN IM DATENSATZ?
# ──────────────────────────────────────────────────────────────
# WICHTIG: Nicht alle Igel wurden gleich lang beobachtet.
# Spät-Zeitpunkte werden von wenigen Langzeitbeobachteten dominiert
# (Survivorship Bias). Das muss transparent gemacht werden.
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("Beobachtungsabdeckung (Survivorship-Check)\n")
cat("════════════════════════════════════\n\n")

# Pro Igel: erster und letzter Beobachtungstag
nacht_roh_all <- dt_tage_ok[tageszeit == "Nacht" & !is.na(pct_aktiv)]
abdeckung <- nacht_roh_all[, .(
  t_start  = min(tage_seit),
  t_ende   = max(tage_seit),
  n_naechte = .N
), by = igel][order(t_ende - t_start, decreasing = TRUE)]

abdeckung[, beob_dauer := t_ende - t_start + 1]

cat("Beobachtungsdauer je Igel:\n")
print(abdeckung[, .(igel, t_start, t_ende, beob_dauer, n_naechte)])
cat("\nMittelwert Beobachtungsdauer:", round(mean(abdeckung$beob_dauer), 1), "Tage\n")
cat("Median:                      ", round(median(abdeckung$beob_dauer), 1), "Tage\n")
cat("Max:                         ", max(abdeckung$beob_dauer), "Tage\n")
cat("Min:                         ", min(abdeckung$beob_dauer), "Tage\n\n")

# N-at-risk: Wie viele Igel pro Zeitfenster?
t_max <- max(abdeckung$t_ende)
zeitfenster <- seq(0, t_max, by = 7)
n_at_risk <- sapply(zeitfenster, function(t) {
  sum(abdeckung$t_start <= t & abdeckung$t_ende >= t)
})
n_risk_dt <- data.table(
  tage_seit = zeitfenster,
  n_igel    = n_at_risk
)
cat("N-at-risk pro Woche:\n")
print(n_risk_dt)
cat("\n")

# ── Plot COV-1: Gantt-Chart — wer war wann im Datensatz ──
# Igel sortiert nach Beobachtungsdauer (längste oben)
abdeckung[, igel_sort := factor(igel, levels = abdeckung[order(beob_dauer)]$igel)]

p_gantt <- ggplot(abdeckung) +
  geom_segment(aes(x = t_start, xend = t_ende,
                   y = igel_sort, yend = igel_sort,
                   color = beob_dauer),
               linewidth = 5, lineend = "round") +
  geom_point(aes(x = t_start, y = igel_sort), color = "white", size = 2) +
  geom_point(aes(x = t_ende,  y = igel_sort), shape = 4,
             color = "grey20", size = 2.5, stroke = 1.2) +
  scale_color_gradientn(
    colours = c("#f7f7f7", "#92c5de", "#0571b0"),
    name    = "Obs. duration\n(days)"
  ) +
  scale_x_continuous(breaks = seq(0, t_max, 7),
                     labels = paste0("T", seq(0, t_max, 7))) +
  labs(
    title    = "Observation window per hedgehog (Gantt chart)",
    subtitle = "White dot = start | × = end | Colour = duration\nFew hedgehogs observed long-term → Survivorship bias!",
    x = "Days since release",
    y = "Hedgehog"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "right"
  )

# ── Plot COV-2: N-at-risk über Zeit ──
p_n_risk <- ggplot(n_risk_dt, aes(x = tage_seit, y = n_igel)) +
  geom_col(fill = "#2C3E6B", alpha = 0.75, width = 5) +
  geom_text(aes(label = n_igel), vjust = -0.4, size = 3.2, color = "grey20") +
  scale_x_continuous(breaks = zeitfenster,
                     labels = paste0("T", zeitfenster)) +
  scale_y_continuous(breaks = 0:max(n_at_risk),
                     limits = c(0, max(n_at_risk) + 1)) +
  labs(
    title    = "N-at-risk: How many hedgehogs per time window?",
    subtitle = "Bar = number of hedgehogs with data in that week\nLate time points: population-level conclusions are limited!",
    x = "Days since release",
    y = "Number of hedgehogs"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

png(file.path(output_ordner, "abdeckung_gantt.png"),
    width = 1600, height = 1000, res = 130)
print(p_gantt)
dev.off()

png(file.path(output_ordner, "abdeckung_n_risk.png"),
    width = 1400, height = 700, res = 130)
print(p_n_risk)
dev.off()

cat("✓ Abdeckungs-Plots gespeichert:\n")
cat("  abdeckung_gantt.png   — Beobachtungsfenster je Igel\n")
cat("  abdeckung_n_risk.png  — N-at-risk pro Zeitfenster\n\n")

# ── Beobachtungsgewichte für Modell A & Frühphasen-GAMM ──
# Igel mit kurzer Beobachtungszeit liefern nur für frühe Zeitpunkte Daten.
# Damit ihr Beitrag nicht durch lange-beobachtete Igel überschattet wird,
# berechnen wir Gewichte: kürzere Beobachtungsdauer → höheres Gewicht
# pro Datenpunkt (da sie seltener im Datensatz vertreten sind).
#
# Gewicht = 1 / (Beobachtungsdauer^0.5)  → moderate Anpassung
# (^0.5 statt ^1 verhindert extreme Gewichte)
abdeckung[, gewicht_raw := 1 / sqrt(beob_dauer)]
abdeckung[, gewicht     := gewicht_raw / mean(gewicht_raw)]  # Mittelwert = 1

cat("Beobachtungsgewichte (kürzere Beobachtung → höheres Gewicht):\n")
print(abdeckung[, .(igel, beob_dauer, gewicht = round(gewicht, 3))][order(beob_dauer)])
cat("\n")

# ──────────────────────────────────────────────────────────────
# EXPLORATIVE VISUALISIERUNG — ROHDATEN VOR DEM MODELL
# ──────────────────────────────────────────────────────────────
# Diese Plots zeigen die Rohdaten OHNE Modellannahmen.
# So sieht man direkt: Gibt es überhaupt einen sichtbaren Trend?
# ──────────────────────────────────────────────────────────────

cat("── Explorative Rohdaten-Plots ──\n")

nacht_roh <- dt_tage_ok[tageszeit == "Nacht" & !is.na(pct_aktiv)]
nacht_roh[, pct_aktiv_pct := pct_aktiv * 100]

# ── Plot E1: Spaghetti-Plot — jede Linie = ein Igel ──
# Zeigt sofort: Sind die Tiere von Anfang an nachtaktiv?
# Gibt es einen steigenden/fallenden Trend?
p_spaghetti <- ggplot(nacht_roh,
                      aes(x = tage_seit, y = pct_aktiv_pct,
                          color = igel, group = igel)) +
  # Grauer Hintergrund fuer Bereich nach Tag 14 (Survivorship Bias Zone)
  annotate("rect",
           xmin = heatmap_vis_cutoff + 0.5,
           xmax = max(nacht_roh$tage_seit, na.rm = TRUE) + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.35) +
  annotate("text",
           x = heatmap_vis_cutoff + 1, y = 95,
           label = "N \u2264 10\n(Survivorship Bias)",
           hjust = 0, vjust = 1, size = 2.8, color = "grey40") +
  geom_vline(xintercept = heatmap_vis_cutoff + 0.5,
             linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_line(alpha = 0.45, linewidth = 0.8) +
  geom_smooth(aes(group = 1),
              method  = "loess", span = 0.4,
              color   = "black", linewidth = 1.5,
              se      = TRUE, fill = "grey30", alpha = 0.2) +
  geom_hline(yintercept = 50, linetype = "dashed",
             color = "firebrick", linewidth = 0.8) +
  annotate("text", x = min(nacht_roh$tage_seit, na.rm = TRUE),
           y = 52, label = "50%-Schwelle", hjust = 0,
           color = "firebrick", size = 3) +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title    = "Nocturnal activity of all hedgehogs — raw data (spaghetti plot)",
    subtitle = paste0("Each line = one hedgehog | Black line = LOESS trend | ",
                      "Grey = N \u2264 10 (from day ", heatmap_vis_cutoff + 1, ")"),
    x        = "Days since release",
    y        = "Proportion of active night minutes (%)",
    color    = "Hedgehog"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

# ── Plot E2: Facetten-Plot — jeder Igel einzeln mit Trendlinie ──
p_facetten <- ggplot(nacht_roh,
                     aes(x = tage_seit, y = pct_aktiv_pct)) +
  geom_point(color = "#2C3E6B", alpha = 0.5, size = 1.2) +
  geom_smooth(method = "loess", span = 0.7,
              color = "#E8A87C", fill = "#E8A87C", alpha = 0.3,
              linewidth = 1) +
  geom_hline(yintercept = 50, linetype = "dashed",
             color = "firebrick", linewidth = 0.5) +
  facet_wrap(~ igel, scales = "free_x") +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Nocturnal activity per hedgehog — raw data with individual trend",
    subtitle = "Orange line = individual LOESS trend | Red = 50% threshold",
    x        = "Days since release",
    y        = "% active night minutes"
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.title     = element_text(face = "bold"),
        strip.text     = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEF2F7", color = NA))

png(file.path(output_ordner, "explorativ_spaghetti.png"),
    width = 1600, height = 900, res = 130)
print(p_spaghetti)
dev.off()

png(file.path(output_ordner, "explorativ_facetten.png"),
    width = 2000, height = 1600, res = 130)
print(p_facetten)
dev.off()
cat("✓ Explorative Plots gespeichert:\n")
cat("  explorativ_spaghetti.png  — alle Igel, Populationstrend\n")
cat("  explorativ_facetten.png   — jeder Igel einzeln\n\n")

# ── Deskriptive Statistik pro Igel ──
cat("── Deskriptiver Überblick (Nachtaktivität) ──\n")
desc_igel <- nacht_roh[, .(
  n_Nächte     = .N,
  Tage_Bereich = paste0(min(tage_seit), "–", max(tage_seit)),
  Mittelwert   = round(mean(pct_aktiv_pct), 1),
  Median       = round(median(pct_aktiv_pct), 1),
  SD           = round(sd(pct_aktiv_pct), 1),
  Min          = round(min(pct_aktiv_pct), 1),
  Max          = round(max(pct_aktiv_pct), 1)
), by = igel][order(Mittelwert, decreasing = TRUE)]
print(desc_igel)
cat("\n")

# Früh- vs. Spätphase — falls Daten reichen
grenze <- 14   # erste 2 Wochen vs. danach
frueh_spaet <- nacht_roh[, phase := ifelse(tage_seit <= grenze, "Früh (≤14d)", "Spät (>14d)")]
frueh_spaet_tab <- frueh_spaet[, .(
  n            = .N,
  Mittelwert   = round(mean(pct_aktiv_pct), 1),
  Median       = round(median(pct_aktiv_pct), 1),
  SD           = round(sd(pct_aktiv_pct), 1)
), by = phase]

cat("── Früh- vs. Spätphase (Populationsmittel) ──\n")
print(frueh_spaet_tab)

# Wilcoxon-Test Früh vs. Spät
if (all(c("Früh (≤14d)", "Spät (>14d)") %in% frueh_spaet$phase)) {
  wt <- wilcox.test(
    nacht_roh[tage_seit <= grenze,  pct_aktiv_pct],
    nacht_roh[tage_seit  > grenze,  pct_aktiv_pct],
    exact = FALSE
  )
  cat("\nWilcoxon-Test Früh vs. Spät: W =", round(wt$statistic, 1),
      ", p =", round(wt$p.value, 3), "\n")
  cat(if (wt$p.value < 0.05)
        "→ SIGNIFIKANTER Unterschied zwischen früher und später Phase!\n"
      else
        "→ Kein signifikanter Unterschied zwischen früher und später Phase.\n")
}
cat("\n")

# ──────────────────────────────────────────────────────────────
# ERSTE-TAGE-ANALYSE — TAGE 1 BIS 5 ISOLIERT
# ──────────────────────────────────────────────────────────────
# Warum isoliert?
#   - In den ersten Tagen sind fast alle Igel noch im Datensatz
#     → kein Survivorship Bias, maximale Vergleichbarkeit
#   - Genau hier würde eine Anpassungszeit sichtbar sein
#   - Danach (Wochen) dominieren nur noch Langzeitbeobachtete
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("Erste-Tage-Analyse (Tage 1–5)\n")
cat("════════════════════════════════════\n\n")

erste_tage_max <- 5

# Nocturnal Index für alle Tage (wird später auch in NI-Abschnitt gebraucht)
ni_alle <- dcast(dt_tage_ok[!is.na(pct_aktiv)],
                 igel + datum + tage_seit ~ tageszeit,
                 value.var = "pct_aktiv",
                 fun.aggregate = mean)

if ("Nacht" %in% names(ni_alle) && "Tag" %in% names(ni_alle)) {
  ni_alle[, nocturnal_index := (Nacht - Tag) * 100]
} else {
  ni_alle[, nocturnal_index := NA_real_]
}

# Nur Tage 1–5 und nur ganzzahlige Tage (Tag 0 = Auswilderungstag, evtl. unvollständig)
et <- ni_alle[tage_seit >= 1 & tage_seit <= erste_tage_max & !is.na(nocturnal_index)]
et[, tag_f := factor(paste0("Tag ", tage_seit), levels = paste0("Tag ", 1:erste_tage_max))]

# N pro Tag (Frühphase — eigene Variable, überschreibt nicht das GAMM n_pro_tag)
n_pro_tag_et <- et[, .(n_igel = uniqueN(igel)), by = tag_f][order(tag_f)]
cat("N Igel pro Tag:\n")
print(n_pro_tag_et)
cat("\n")

# Deskriptive Statistik pro Tag
desc_et <- et[, .(
  n         = uniqueN(igel),
  median_NI = round(median(nocturnal_index, na.rm = TRUE), 1),
  mean_NI   = round(mean(nocturnal_index,   na.rm = TRUE), 1),
  sd_NI     = round(sd(nocturnal_index,     na.rm = TRUE), 1),
  min_NI    = round(min(nocturnal_index,    na.rm = TRUE), 1),
  max_NI    = round(max(nocturnal_index,    na.rm = TRUE), 1),
  pct_nachtaktiv = round(mean(nocturnal_index > 0, na.rm = TRUE) * 100, 1)
), by = tag_f][order(tag_f)]

cat("Deskriptive Statistik Nocturnal Index (Tage 1–5):\n")
print(desc_et)
cat("\n")
cat("(pct_nachtaktiv = Anteil Igel mit NI > 0, also nachts aktiver als tags)\n\n")

# Statistischer Test: Verändert sich der NI über die ersten 5 Tage?
if (nlevels(factor(et$tage_seit)) >= 2) {
  # Kruskal-Wallis (nichtparametrisch, kleine N)
  kw <- kruskal.test(nocturnal_index ~ tag_f, data = et)
  cat("Kruskal-Wallis-Test NI über Tage 1–5:\n")
  cat("  H(", kw$parameter, ") =", round(kw$statistic, 2),
      ", p =", round(kw$p.value, 3), "\n")
  cat(if (kw$p.value < 0.05)
        "  → SIGNIFIKANT: Der Nocturnal Index unterscheidet sich zwischen den ersten Tagen.\n"
      else
        "  → Nicht signifikant: Kein messbarer Unterschied zwischen Tag 1 und Tag 5.\n")

  # Paarweise Wilcoxon Tag 1 vs. Tag 5 (falls beide vorhanden)
  if (sum(et$tage_seit == 1) >= 3 && sum(et$tage_seit == 5) >= 3) {
    wt15 <- wilcox.test(
      et[tage_seit == 1, nocturnal_index],
      et[tage_seit == 5, nocturnal_index],
      exact = FALSE
    )
    cat("\nWilcoxon Tag 1 vs. Tag 5:",
        "W =", round(wt15$statistic, 1),
        ", p =", round(wt15$p.value, 3), "\n")
  }
}
cat("\n")

# ── Plot ET-1: Boxplot + Einzelpunkte pro Tag ──────────────────
p_et_box <- ggplot(et, aes(x = tag_f, y = nocturnal_index)) +
  geom_hline(yintercept = 0, color = "firebrick",
             linetype = "dashed", linewidth = 0.8) +
  geom_violin(fill = "#2C3E6B", alpha = 0.15, color = NA, width = 0.8) +
  geom_boxplot(fill = "#2C3E6B", alpha = 0.5,
               width = 0.35, outlier.shape = NA, color = "#1a2744") +
  geom_jitter(aes(color = igel), width = 0.12, size = 2.5, alpha = 0.8) +
  geom_text(data = n_pro_tag_et,
            aes(x = tag_f, y = min(et$nocturnal_index, na.rm=TRUE) - 5,
                label = paste0("n=", n_igel)),
            size = 3, color = "grey40", inherit.aes = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_color_viridis_d(option = "turbo") +
  annotate("text",
           x = 0.6, y = 3,
           label = "→ nocturnal", hjust = 0,
           color = "grey40", size = 3, fontface = "italic") +
  annotate("text",
           x = 0.6, y = -3,
           label = "→ diurnal", hjust = 0,
           color = "firebrick", size = 3, fontface = "italic") +
  labs(
    title    = "Nocturnal Index in the first 5 days after release",
    subtitle = "Each dot = one hedgehog on that day | n = number of hedgehogs\nRed line = zero line (no day/night difference)",
    x = NULL, y = "Nocturnal Index (night% − day%)",
    color = "Hedgehog"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )

# ── Plot ET-2: Spaghetti nur Tage 1–5, jeder Igel als Linie ──
p_et_spaghetti <- ggplot(et,
    aes(x = tage_seit, y = nocturnal_index,
        color = igel, group = igel)) +
  geom_hline(yintercept = 0, color = "firebrick",
             linetype = "dashed", linewidth = 0.8) +
  geom_line(alpha = 0.6, linewidth = 1.0) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "loess", span = 1,
              color = "black", linewidth = 1.5, se = TRUE,
              fill = "grey30", alpha = 0.2) +
  scale_x_continuous(breaks = 1:erste_tage_max,
                     labels = paste0("Tag ", 1:erste_tage_max)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title    = "Individual trajectories: Nocturnal Index days 1–5",
    subtitle = "Each line = one hedgehog | Black = population trend (LOESS)\nDots = observed nights",
    x = NULL, y = "Nocturnal Index (%)",
    color = "Hedgehog"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

# ── Plot ET-3: Heatmap Igel × Tag (NI als Farbe) ──────────────
# Zeigt auf einen Blick wer wann wie aktiv war
et_heat <- et[, .(nocturnal_index = mean(nocturnal_index)), by = .(igel, tage_seit)]
# Igel sortiert nach mittlerem NI über alle Tage
igel_order_ni <- et_heat[, .(mean_ni = mean(nocturnal_index)), by = igel][order(mean_ni)]$igel
et_heat[, igel_f := factor(igel, levels = igel_order_ni)]

p_et_heat <- ggplot(et_heat,
    aes(x = factor(tage_seit), y = igel_f, fill = nocturnal_index)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(nocturnal_index, 0)),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_gradient2(
    low     = "#c0392b",   # tagaktiv
    mid     = "grey90",
    high    = "#2C3E6B",   # nachtaktiv
    midpoint = 0,
    name   = "Nocturnal\nIndex (%)"
  ) +
  scale_x_discrete(labels = paste0("Day ", 1:erste_tage_max)) +
  labs(
    title    = "Heatmap: Nocturnal Index per hedgehog × day",
    subtitle = "Blue = nocturnal | Red = diurnal | Grey = no difference\nEmpty cells = no data for this day",
    x = NULL, y = "Hedgehog (sorted by mean NI)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

png(file.path(output_ordner, "erste_tage_boxplot.png"),
    width = 1400, height = 900, res = 130)
print(p_et_box)
dev.off()

png(file.path(output_ordner, "erste_tage_spaghetti.png"),
    width = 1400, height = 900, res = 130)
print(p_et_spaghetti)
dev.off()

png(file.path(output_ordner, "erste_tage_heatmap.png"),
    width = 1200, height = 900, res = 130)
print(p_et_heat)
dev.off()

cat("✓ Erste-Tage-Plots gespeichert:\n")
cat("  erste_tage_boxplot.png   — Verteilung NI pro Tag (Boxplot + Punkte)\n")
cat("  erste_tage_spaghetti.png — Individuelle Verläufe Tag 1–5\n")
cat("  erste_tage_heatmap.png   — Heatmap: wer war wann wie nachtaktiv\n\n")

# ──────────────────────────────────────────────────────────────
# NOCTURNAL INDEX & FRÜPHASEN-ANALYSE
# ──────────────────────────────────────────────────────────────
#
# WARUM?
# Das Standard-GAMM (Modell A) testet die absolute Nachtaktivität
# über den gesamten Zeitraum. Dabei können kurzfristige Anpassungen
# in der ersten Woche "verschwinden", weil:
#   (1) der Smooth über viele Wochen geglättet wird
#   (2) individuelle Unterschiede im Aktivitätsniveau den Trend verdecken
#
# Besser für die Frage "Gibt es eine Anpassungszeit?":
#   → Nocturnal Index = Nacht% - Tag%  (relativ, kontrolliert Gesamtaktivität)
#   → log(tage_seit + 1) als Zeitachse (streckt frühe Tage auseinander)
#   → Modell nur auf Tage 1–21 (fokussiert auf Anpassungsphase)
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("Nocturnal Index & Frühphasen-Analyse\n")
cat("════════════════════════════════════\n\n")

# ── Nocturnal Index — bereits in Erste-Tage-Analyse berechnet ──
# ni_alle wurde oben als vollständiger NI-Datensatz erstellt,
# hier weiter verwenden (alle Tage, nicht nur 1–5)
ni_wide <- ni_alle

if (any(!is.na(ni_wide$nocturnal_index))) {
  cat("Nocturnal Index (alle Tage, n =", nrow(ni_wide[!is.na(nocturnal_index)]), "Igel×Tage)\n")
  cat("  Bereich:", round(min(ni_wide$nocturnal_index, na.rm=TRUE),1),
      "bis", round(max(ni_wide$nocturnal_index, na.rm=TRUE),1), "\n")
  cat("  Mittelwert:", round(mean(ni_wide$nocturnal_index, na.rm=TRUE),1),
      "| > 0 = nachtaktiv\n\n")
} else {
  cat("HINWEIS: Nocturnal Index konnte nicht berechnet werden.\n\n")
}

# ── Plot NI-1: Nocturnal Index Spaghetti ──
if (any(!is.na(ni_wide$nocturnal_index))) {
  ni_ok <- ni_wide[!is.na(nocturnal_index)]
  ni_ok[, log_tage := log(tage_seit + 1)]

  p_ni_spaghetti <- ggplot(ni_ok,
      aes(x = tage_seit, y = nocturnal_index, color = igel, group = igel)) +
    annotate("rect",
             xmin = heatmap_vis_cutoff + 0.5,
             xmax = max(ni_ok$tage_seit, na.rm = TRUE) + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = "grey80", alpha = 0.35) +
    annotate("text",
             x = heatmap_vis_cutoff + 1,
             y = max(ni_ok$nocturnal_index, na.rm = TRUE) * 0.9,
             label = paste0("N \u2264 10\n(ab Tag ", heatmap_vis_cutoff + 1, ")"),
             hjust = 0, vjust = 1, size = 2.5, color = "grey40") +
    geom_vline(xintercept = heatmap_vis_cutoff + 0.5,
               linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "firebrick", linewidth = 0.8, linetype = "dashed") +
    geom_line(alpha = 0.4, linewidth = 0.8) +
    geom_smooth(aes(group = 1), method = "loess", span = 0.5,
                color = "black", linewidth = 1.5, se = TRUE,
                fill = "grey30", alpha = 0.2) +
    annotate("text", x = min(ni_ok$tage_seit, na.rm=TRUE),
             y = 2, label = "\u2192 nocturnal", hjust = 0,
             color = "grey30", size = 3, fontface = "italic") +
    annotate("text", x = min(ni_ok$tage_seit, na.rm=TRUE),
             y = -2, label = "\u2192 diurnal", hjust = 0,
             color = "firebrick", size = 3, fontface = "italic") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title    = "Nocturnal Index over time (night% − day%)",
      subtitle = paste0("0 = equal | Positive = nocturnal | Negative = diurnal | ",
                        "Grey = N \u2264 10 (from day ", heatmap_vis_cutoff + 1, ")"),
      x = "Days since release", y = "Nocturnal Index (%)", color = "Hedgehog"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8, color = "grey40"))

  # ── Plot NI-2: Frühphase LOESS auf log-Zeitskala ──
  # Die logarithmische x-Achse streckt Tag 1–7 auseinander
  # und zeigt besser ob sich was in der ersten Woche tut
  p_ni_log <- ggplot(ni_ok[tage_seit >= 1],
      aes(x = tage_seit, y = nocturnal_index, color = igel, group = igel)) +
    geom_hline(yintercept = 0, color = "firebrick", linewidth = 0.8, linetype = "dashed") +
    geom_line(alpha = 0.35, linewidth = 0.7) +
    geom_smooth(aes(group = 1), method = "loess", span = 0.6,
                color = "black", linewidth = 1.5, se = TRUE,
                fill = "grey30", alpha = 0.2) +
    scale_x_log10(breaks = c(1, 2, 3, 5, 7, 14, 21, 42),
                  labels = c("T1","T2","T3","T5","T7","T14","T21","T42")) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title    = "Nocturnal Index — log time axis (focus: early days)",
      subtitle = "Logarithmic x-axis: days 1–7 are stretched apart\nBlack = population trend | Red = zero line (no preference)",
      x = "Days since release (log scale)", y = "Nocturnal Index (%)", color = "Hedgehog"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  # ── Plot NI-3: Facetten pro Igel mit log-Achse ──
  p_ni_facetten <- ggplot(ni_ok[tage_seit >= 1],
      aes(x = tage_seit, y = nocturnal_index)) +
    geom_hline(yintercept = 0, color = "firebrick", linewidth = 0.6, linetype = "dashed") +
    geom_point(color = "#2C3E6B", alpha = 0.5, size = 1.3) +
    geom_smooth(method = "loess", span = 0.8,
                color = "#E8A87C", fill = "#E8A87C", alpha = 0.3, linewidth = 1) +
    facet_wrap(~ igel, scales = "free_x") +
    scale_x_log10(breaks = c(1, 2, 3, 5, 7, 14, 21, 42),
                  labels = c("T1","T2","T3","T5","T7","T14","T21","T42")) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title    = "Nocturnal Index per hedgehog (log time scale)",
      subtitle = "Orange line = individual LOESS trend | Red = zero line",
      x = "Days since release (log scale)", y = "Nocturnal Index (%)"
    ) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "#EEF2F7", color = NA))

  png(file.path(output_ordner, "ni_spaghetti.png"),
      width = 1600, height = 900, res = 130)
  print(p_ni_spaghetti)
  dev.off()

  png(file.path(output_ordner, "ni_log_zeitachse.png"),
      width = 1600, height = 900, res = 130)
  print(p_ni_log)
  dev.off()

  png(file.path(output_ordner, "ni_facetten_log.png"),
      width = 2000, height = 1600, res = 130)
  print(p_ni_facetten)
  dev.off()

  cat("✓ Nocturnal-Index-Plots gespeichert:\n")
  cat("  ni_spaghetti.png      — alle Igel, lineare Zeitachse\n")
  cat("  ni_log_zeitachse.png  — alle Igel, LOG-Zeitachse (frühe Tage sichtbar)\n")
  cat("  ni_facetten_log.png   — jeder Igel einzeln, LOG-Zeitachse\n\n")
}

# ── Frühphasen-GAMM (Tage 1–21, log-Zeit) ─────────────────────
cat("── Frühphasen-GAMM (Anpassungszeit, Tage 1–21) ──\n")
cat("Zeitskala: log(tage_seit + 1) — sensitiver für frühe Veränderungen\n\n")

frueh_tage   <- 21   # Analysefenster für Frühphase
ni_frueh     <- ni_wide[!is.na(nocturnal_index) & tage_seit <= frueh_tage & tage_seit >= 1]
ni_frueh[, log_tage := log(tage_seit + 1)]
ni_frueh[, igel      := factor(igel)]

cat("Datenpunkte Frühphase:", nrow(ni_frueh), "\n")
cat("Igel in Frühphase:   ", nlevels(ni_frueh$igel), "\n\n")

if (nrow(ni_frueh) >= 20 && nlevels(ni_frueh$igel) >= 3) {

  # Transformation für Beta-Regression
  n_tot <- ni_frueh[, .N, by = igel][, mean(N)]
  # Nocturnal Index liegt in (-100, 100) → auf (0,1) skalieren
  ni_frueh[, ni_01       := (nocturnal_index + 100) / 200]
  ni_frueh[, ni_transf   := (ni_01 * (n_tot - 1) + 0.5) / n_tot]
  ni_frueh    <- ni_frueh[ni_transf > 0 & ni_transf < 1]

  set.seed(42)
  modell_frueh <- tryCatch(
    bam(
      ni_transf ~
        s(log_tage, bs = "tp", k = 6) +
        s(log_tage, igel, bs = "fs", k = 4, m = 1) +
        s(igel, bs = "re"),
      family = betar(link = "logit"),
      data   = ni_frueh,
      method = "fREML"
    ),
    error = function(e) {
      cat("HINWEIS: Frühphasen-GAMM konnte nicht gefittet werden:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(modell_frueh)) {
    cat("── Frühphasen-Modell: Zusammenfassung ──\n")
    print(summary(modell_frueh))

    # Plot: vorhergesagter Trend in der Frühphase (zurück in Original-Skala)
    nd_frueh <- data.table(
      log_tage = seq(log(1 + 1), log(frueh_tage + 1), length.out = 200),
      igel     = factor(levels(ni_frueh$igel)[1], levels = levels(ni_frueh$igel))
    )
    pr_frueh <- predict(modell_frueh, newdata = nd_frueh,
                        exclude = c("s(log_tage,igel)", "s(igel)"),
                        se.fit = TRUE)
    nd_frueh[, fit_01 := plogis(pr_frueh$fit)]
    nd_frueh[, lwr_01 := plogis(pr_frueh$fit - 1.96 * pr_frueh$se.fit)]
    nd_frueh[, upr_01 := plogis(pr_frueh$fit + 1.96 * pr_frueh$se.fit)]
    # Zurück in NI-Skala: NI = (01 * 200) - 100
    nd_frueh[, fit := fit_01 * 200 - 100]
    nd_frueh[, lwr := lwr_01 * 200 - 100]
    nd_frueh[, upr := upr_01 * 200 - 100]
    nd_frueh[, tage := exp(log_tage) - 1]

    p_frueh_gamm <- ggplot() +
      geom_point(data = ni_frueh,
                 aes(x = exp(log_tage) - 1, y = nocturnal_index, color = igel),
                 alpha = 0.45, size = 1.8) +
      geom_hline(yintercept = 0, color = "firebrick", linewidth = 0.8, linetype = "dashed") +
      geom_ribbon(data = nd_frueh,
                  aes(x = tage, ymin = lwr, ymax = upr),
                  fill = "#2C3E6B", alpha = 0.2) +
      geom_line(data = nd_frueh,
                aes(x = tage, y = fit),
                color = "#2C3E6B", linewidth = 1.5) +
      scale_x_continuous(breaks = c(1,2,3,5,7,10,14,21),
                         labels = paste0("T", c(1,2,3,5,7,10,14,21))) +
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      scale_color_viridis_d(option = "turbo") +
      labs(
        title    = "Early-phase GAMM: Nocturnal Index days 1–21",
        subtitle = paste0(
          "Log time scale | Population trend (blue) with 95% CI\n",
          "s(log_tage): p = ",
          {
            pv <- summary(modell_frueh)$s.table
            idx <- grep("log_tage\\)", rownames(pv), fixed = TRUE)[1]
            if (!is.na(idx) && pv[idx, "p-value"] < 0.001) "< 0.001"
            else if (!is.na(idx)) as.character(round(pv[idx, "p-value"], 3))
            else "n.v."
          }
        ),
        x = "Days since release", y = "Nocturnal Index (%)", color = "Hedgehog"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))

    png(file.path(output_ordner, "gamm_fruehphase.png"),
        width = 1400, height = 900, res = 130)
    print(p_frueh_gamm)
    dev.off()
    cat("✓ Frühphasen-Plot gespeichert: gamm_fruehphase.png\n\n")
  }
} else {
  modell_frueh <- NULL
  cat("Zu wenig Daten für Frühphasen-GAMM (n =", nrow(ni_frueh),
      ", Igel =", nlevels(ni_frueh$igel), ")\n\n")
}

# ──────────────────────────────────────────────────────────────
# OPTION A: TAGESSUMMARY — NACHTAKTIVITÄT ÜBER ZEIT
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("Option A: GAMM auf Tagesebene\n")
cat("Response: Anteil aktiver Nachtminuten pro Tag\n")
cat("════════════════════════════════════\n\n")

# Nur Nachtzeilen; Beta-Regression braucht (0,1) — Werte genau 0 oder 1
# werden mit Smithson & Verkuilen (2006) Transformation angepasst
nacht_tage <- dt_tage_ok[tageszeit == "Nacht" & !is.na(pct_aktiv)]
nacht_tage[, n_total := n_min]
# Transformation: (y*(n-1) + 0.5) / n  → verschiebt 0 und 1 leicht weg vom Rand
nacht_tage[, pct_transf := (pct_aktiv * (n_total - 1) + 0.5) / n_total]
# Sicherheitscheck: alles im (0,1) Intervall
nacht_tage <- nacht_tage[pct_transf > 0 & pct_transf < 1]

# ── Beobachtungsgewichte in Modell A einbauen ──────────────────
# Igel mit kurzem Beobachtungsfenster liefern nur Daten für frühe
# Zeitpunkte → ihre Datenpunkte sind in der Gesamtstichprobe
# seltener vertreten als die Datenpunkte lang-beobachteter Igel.
# Ohne Gewichtung dominieren Langebeobachtete die späten Zeitpunkte.
# Lösung: jeder Datenpunkt eines Igels erhält das vorab berechnete
# Gewicht (kürzer beobachtet = höheres Gewicht pro Punkt).
nacht_tage <- merge(nacht_tage,
                    abdeckung[, .(igel, gewicht)],
                    by = "igel", all.x = TRUE)
nacht_tage[is.na(gewicht), gewicht := 1]  # Fallback

cat("Datenpunkte für Modell A:", nrow(nacht_tage), "\n")
cat("Tage-Bereich:            ", nacht_tage[, paste(min(tage_seit), "–", max(tage_seit))], "\n")
cat("Gewichte: min =", round(min(nacht_tage$gewicht),3),
    "| max =", round(max(nacht_tage$gewicht),3),
    "| Mittelwert =", round(mean(nacht_tage$gewicht),3), "\n\n")

# Modell A fitten
# s(tage_seit):              Populationstrend (nichtlinear, Thin-Plate-Spline)
# s(tage_seit, igel, bs="fs"): Individuelle Abweichungen mit Shrinkage
#                              (Igel mit wenig Daten → stärker Richtung Populationsmittel)
# s(igel, bs="re"):           Zufälliges Intercept je Igel
# weights = gewicht:          Korrektur für Survivorship Bias
#                             (kurz beobachtete Igel erhalten höheres Gewicht)
set.seed(42)
modell_A <- bam(
  pct_transf ~
    s(tage_seit, bs = "tp", k = 8) +
    s(tage_seit, igel, bs = "fs", k = 5, m = 1) +
    s(igel, bs = "re"),
  family  = betar(link = "logit"),
  data    = nacht_tage,
  weights = nacht_tage$gewicht,
  method  = "fREML"
)

# Vergleichsmodell OHNE Gewichte (für Sensitivitätscheck)
set.seed(42)
modell_A_ungewichtet <- bam(
  pct_transf ~
    s(tage_seit, bs = "tp", k = 8) +
    s(tage_seit, igel, bs = "fs", k = 5, m = 1) +
    s(igel, bs = "re"),
  family = betar(link = "logit"),
  data   = nacht_tage,
  method = "fREML"
)

cat("── Sensitivitätscheck: Modell A mit vs. ohne Gewichte ──\n")
cat("MIT Gewichten:    s(tage_seit) p =",
    round(summary(modell_A)$s.table[
      grep("tage_seit\\)", rownames(summary(modell_A)$s.table),
           fixed=FALSE)[1], "p-value"], 3), "\n")
cat("OHNE Gewichte:    s(tage_seit) p =",
    round(summary(modell_A_ungewichtet)$s.table[
      grep("tage_seit\\)", rownames(summary(modell_A_ungewichtet)$s.table),
           fixed=FALSE)[1], "p-value"], 3), "\n")
cat("→ Wenn beide ähnlich: Survivorship Bias hat geringen Einfluss.\n")
cat("→ Wenn verschieden: Ergebnis hängt von den Langzeitbeobachteten ab!\n\n")

cat("── Modell A: Zusammenfassung ──\n")
print(summary(modell_A))

# Autokorrelation in Residuen prüfen (Tage innerhalb Igel)
cat("\n── Autokorrelation der Residuen (ACF) ──\n")
png(file.path(output_ordner, "gamm_A_diagnostik.png"),
    width = 1400, height = 900, res = 130)
par(mfrow = c(2, 2))
gam.check(modell_A)
dev.off()
cat("✓ Diagnostikplot gespeichert: gamm_A_diagnostik.png\n")

# Konkurrenz zwischen Prädiktoren prüfen (Concurvity)
cat("\n── Concurvity (sollte < 0.8 sein) ──\n")
print(round(concurvity(modell_A, full = TRUE), 3))

# ── Plot A1: Populationstrend mit Konfidenzband ──
# Vorhersagegrid ueber vollen Modellbereich (alle Daten genutzt)
tage_seq <- data.table(
  tage_seit = seq(nacht_tage[, min(tage_seit)],
                  nacht_tage[, max(tage_seit)], length.out = 200),
  igel      = factor(levels(nacht_tage$igel)[1], levels = levels(nacht_tage$igel))
)
pred_A <- predict(modell_A, newdata = tage_seq,
                  exclude = c("s(tage_seit,igel)", "s(igel)"),
                  se.fit = TRUE)
tage_seq[, fit    := plogis(pred_A$fit)]
tage_seq[, lwr    := plogis(pred_A$fit - 1.96 * pred_A$se.fit)]
tage_seq[, upr    := plogis(pred_A$fit + 1.96 * pred_A$se.fit)]

# Rohdaten für Hintergrundpunkte
roh_punkte <- nacht_tage[, .(igel, tage_seit, pct_aktiv)]

pA1 <- ggplot() +
  # Grauer Bereich nach Tag 14 (Survivorship Bias Zone)
  annotate("rect",
           xmin = heatmap_vis_cutoff + 0.5,
           xmax = nacht_tage[, max(tage_seit)] + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.35) +
  annotate("text",
           x = heatmap_vis_cutoff + 1, y = 0.97,
           label = "N \u2264 10\nHedgehog22+23 dominate",
           hjust = 0, vjust = 1, size = 2.5, color = "grey40") +
  geom_vline(xintercept = heatmap_vis_cutoff + 0.5,
             linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_point(data = roh_punkte,
             aes(x = tage_seit, y = pct_aktiv, color = igel),
             alpha = 0.35, size = 1.5) +
  geom_ribbon(data = tage_seq,
              aes(x = tage_seit, ymin = lwr, ymax = upr),
              fill = "#2C3E6B", alpha = 0.2) +
  geom_line(data = tage_seq,
            aes(x = tage_seit, y = fit),
            color = "#2C3E6B", linewidth = 1.4) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    title    = "Nocturnal activity over time — population trend (GAMM)",
    subtitle = paste0("Line = population smooth with 95% CI | Points = individual nights per hedgehog | ",
                      "Grey = N \u2264 10 (from day ", heatmap_vis_cutoff + 1, ")"),
    x = "Days since release",
    y = "Proportion of active night minutes",
    color = "Hedgehog"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        legend.position = "right")

# ── Plot A2: Individuelle Kurven ──
# Vorhersage für jeden Igel über die Zeit
igel_lvl <- levels(nacht_tage$igel)
pred_list <- lapply(igel_lvl, function(ig) {
  tage_ig <- nacht_tage[igel == ig, .(min_t = min(tage_seit), max_t = max(tage_seit))]
  nd <- data.table(
    tage_seit = seq(tage_ig$min_t, tage_ig$max_t, length.out = 80),
    igel = factor(ig, levels = igel_lvl)
  )
  pr <- predict(modell_A, newdata = nd, se.fit = TRUE)
  nd[, fit := plogis(pr$fit)]
  nd[, lwr := plogis(pr$fit - 1.96 * pr$se.fit)]
  nd[, upr := plogis(pr$fit + 1.96 * pr$se.fit)]
  nd
})
pred_igel <- rbindlist(pred_list)

pA2 <- ggplot(pred_igel,
              aes(x = tage_seit, y = fit, color = igel, fill = igel)) +
  annotate("rect",
           xmin = heatmap_vis_cutoff + 0.5,
           xmax = max(pred_igel$tage_seit) + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.35) +
  geom_vline(xintercept = heatmap_vis_cutoff + 0.5,
             linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    title    = "Individual GAMM curves per hedgehog",
    subtitle = paste0("Dashed line = 50% threshold | Grey = N \u2264 10 (from day ",
                      heatmap_vis_cutoff + 1, ")"),
    x = "Days since release",
    y = "Proportion of active night minutes",
    color = "Hedgehog", fill = "Hedgehog"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        legend.position = "right")

# ── Plot A3: N-at-risk als Unterplot zum Populationstrend ──
# Zeigt wie viele Igel zu jedem Zeitpunkt beitragen
# → transparenz über Datengrundlage des Trends

# N pro Zeitpunkt (kontinuierlich über die tage_seq-Punkte)
n_at_seq <- sapply(tage_seq$tage_seit, function(t) {
  sum(abdeckung$t_start <= t & abdeckung$t_ende >= t)
})
tage_seq[, n_igel := n_at_seq]

pA_n <- ggplot(tage_seq, aes(x = tage_seit, y = n_igel)) +
  # Bereich nach Tag 14 schattieren
  annotate("rect",
           xmin = heatmap_vis_cutoff + 0.5,
           xmax = max(tage_seq$tage_seit) + 0.5,
           ymin = 0, ymax = Inf,
           fill = "grey80", alpha = 0.5) +
  geom_vline(xintercept = heatmap_vis_cutoff + 0.5,
             linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_area(fill = "#2C3E6B", alpha = 0.35) +
  geom_line(color = "#2C3E6B", linewidth = 0.8) +
  # Bereich mit N < 5 rot markieren
  geom_rect(data = tage_seq[n_igel < 5],
            aes(xmin = tage_seit - 0.5, xmax = tage_seit + 0.5,
                ymin = 0, ymax = n_igel),
            fill = "firebrick", alpha = 0.3, inherit.aes = FALSE) +
  scale_y_continuous(breaks = function(x) unique(floor(pretty(x))),
                     limits = c(0, NA)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    x = "Days since release",
    y = "N hedgehogs",
    caption = paste0("Dashed line = day ", heatmap_vis_cutoff,
                     " (visualisation limit) | Red = N < 5 (no population-level conclusions)")
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.caption = element_text(color = "grey50", size = 7))

# Plots kombinieren und speichern
# Drei-Panel: Populationstrend / Individuelle Kurven / N-at-risk
png(file.path(output_ordner, "gamm_A_nachtaktivitaet.png"),
    width = 1600, height = 1800, res = 130)
print(pA1 / pA2 / pA_n + patchwork::plot_layout(heights = c(3, 3, 1)))
dev.off()
cat("✓ Plot gespeichert: gamm_A_nachtaktivitaet.png (inkl. N-at-risk Unterplot)\n\n")

# ──────────────────────────────────────────────────────────────
# OPTION B: MINUTENEBENE — VERSCHIEBT SICH DER 24H-RHYTHMUS?
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("Option B: GAMM auf Minutenebene (Nacht)\n")
cat("Response: aktiv/passiv ~ s(tage_seit) + s(stunde) + ti(tage_seit,stunde)\n")
cat("Hinweis: kann bei vielen Daten einige Minuten dauern.\n")
cat("════════════════════════════════════\n\n")

cat("Datenpunkte Nacht:", format(nrow(dt_nacht_ok), big.mark = "'"), "\n\n")

# Modell B:
# s(tage_seit):                Zeittrend gesamt
# s(stunde, bs="cc"):          24h-Rhythmus (zirkulär: 0h = 24h)
# ti(tage_seit, stunde):       Interaktion: verändert sich der Tagesrhythmus über Zeit?
# s(igel, bs="re"):            Zufälliges Intercept je Igel
# s(tage_seit, igel, bs="fs"): Individueller Zeittrend (Shrinkage)

# Modell B vereinfacht: Factor Smooth entfernt (verursacht numerische Instabilität
# im Binomial-Modell durch quasi-perfekte Separation).
# s(stunde, bs="cr") statt "cc" (cyclic), da Nachtdaten nicht den vollen
# 0-24h Bereich abdecken → zirkuläre Randbedingung wäre falsch.
modell_B <- bam(
  aktiv ~
    s(tage_seit, bs = "tp", k = 8) +
    s(stunde,    bs = "cr", k = 10) +
    ti(tage_seit, stunde, bs = c("tp", "cr"), k = c(5, 6)) +
    s(igel, bs = "re"),
  family   = binomial(link = "logit"),
  data     = dt_nacht_ok,
  method   = "fREML",
  discrete = TRUE
)

cat("── Modell B: Zusammenfassung ──\n")
print(summary(modell_B))

png(file.path(output_ordner, "gamm_B_diagnostik.png"),
    width = 1400, height = 900, res = 130)
par(mfrow = c(2, 2))
gam.check(modell_B)
dev.off()
cat("✓ Diagnostikplot gespeichert: gamm_B_diagnostik.png\n")

# ── Plot B1: 2D-Heatmap Tage × Stunde (Populationseffekt) ──
# Modell wird auf ALLEN verfügbaren Daten geschätzt (bessere Schätzung der
# Zufallseffekte). Die Visualisierung ist auf Tag 1–heatmap_vis_cutoff begrenzt
# (definiert im Einstellungs-Block oben), weil ab Tag ~15 nur noch ≤10 Tiere
# im Datensatz sind (Survivorship Bias) und Igel22+Igel23 dominieren.

grid_B <- expand.grid(
  tage_seit = seq(dt_nacht_ok[, min(tage_seit)],
                  min(dt_nacht_ok[, max(tage_seit)], heatmap_vis_cutoff),
                  length.out = 60),
  stunde    = seq(0, 23.5, by = 0.5)
)
setDT(grid_B)
grid_B[, igel := factor(levels(dt_nacht_ok$igel)[1],
                         levels = levels(dt_nacht_ok$igel))]
pred_B <- predict(modell_B, newdata = grid_B,
                  exclude = "s(igel)",
                  type = "response")
grid_B[, p_aktiv := pred_B]

# N-at-risk-Label für Heatmap-Untertitel
n_t1  <- n_pro_tag[tage_seit == 1,  n_igel] |> head(1)
n_t7  <- n_pro_tag[tage_seit == 7,  n_igel] |> head(1)
n_t14 <- n_pro_tag[tage_seit == 14, n_igel] |> head(1)
n_label <- paste0("N-at-risk: Tag 1 = ", n_t1 %||% "?",
                   " | Tag 7 = ", n_t7 %||% "?",
                   " | Tag 14 = ", n_t14 %||% "?")

pB1 <- ggplot(grid_B, aes(x = tage_seit, y = stunde, fill = p_aktiv)) +
  geom_tile() +
  scale_fill_gradientn(
    colours  = c("#0d1b2a", "#1e3a5f", "#e9c46a", "#e76f51"),
    limits   = c(0, 1),
    labels   = percent_format(accuracy = 1),
    name     = "P(active)"
  ) +
  scale_y_continuous(breaks = seq(0, 23, 3),
                     labels = sprintf("%02d:00", seq(0, 23, 3))) +
  scale_x_continuous(breaks = seq(1, heatmap_vis_cutoff, 1)) +
  labs(
    title    = paste0("Activity probability — days 1–", heatmap_vis_cutoff,
                      " after release"),
    subtitle = paste0("Population effect (model on all available data) | ",
                      n_label, "\n",
                      "Note: From day 15, Hedgehog22 (blindness) & Hedgehog23 ",
                      "dominate the estimate (survivorship bias)"),
    x = "Days since release",
    y = "Time of day (CET/CEST)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"))

# ── Plot B2: Zeittrend (marginaler Effekt, Stunde rausgemittelt) ──
# Voller Modellbereich wird gezeigt; ab Tag heatmap_vis_cutoff grau schattiert
# (gleiche Logik wie Spaghetti/A1/A2: Daten sichtbar, Unsicherheitsbereich markiert)
grid_zeit <- data.table(
  tage_seit = seq(dt_nacht_ok[, min(tage_seit)],
                  dt_nacht_ok[, max(tage_seit)], length.out = 200),
  stunde    = mean(dt_nacht_ok$stunde),
  igel      = factor(levels(dt_nacht_ok$igel)[1],
                     levels = levels(dt_nacht_ok$igel))
)
pred_zeit <- predict(modell_B, newdata = grid_zeit,
                     exclude = c("s(igel)",
                                 "s(stunde)", "ti(tage_seit,stunde)"),
                     se.fit = TRUE, type = "link")
grid_zeit[, fit := plogis(pred_zeit$fit)]
grid_zeit[, lwr := plogis(pred_zeit$fit - 1.96 * pred_zeit$se.fit)]
grid_zeit[, upr := plogis(pred_zeit$fit + 1.96 * pred_zeit$se.fit)]

pB2 <- ggplot(grid_zeit, aes(x = tage_seit)) +
  annotate("rect",
           xmin = heatmap_vis_cutoff + 0.5,
           xmax = dt_nacht_ok[, max(tage_seit)] + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.35) +
  annotate("text",
           x = heatmap_vis_cutoff + 1,
           y = max(grid_zeit$upr, na.rm = TRUE) * 0.97,
           label = paste0("N \u2264 10\nHedgehog22+23\ndominate"),
           hjust = 0, vjust = 1, size = 2.5, color = "grey40") +
  geom_vline(xintercept = heatmap_vis_cutoff + 0.5,
             linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#2C3E6B", alpha = 0.2) +
  geom_line(aes(y = fit), color = "#2C3E6B", linewidth = 1.3) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    title    = "Marginal time trend — s(tage_seit) from Model B",
    subtitle = paste0("Hour averaged out | Grey = N \u2264 10 (from day ",
                      heatmap_vis_cutoff + 1, ") — curve shown for completeness"),
    x = "Days since release",
    y = "P(active) — population mean"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"))

png(file.path(output_ordner, "gamm_B_rhythmus.png"),
    width = 1600, height = 1400, res = 130)
print(pB1 / pB2)
dev.off()
cat("✓ Plot gespeichert: gamm_B_rhythmus.png\n\n")

# ──────────────────────────────────────────────────────────────
# BOUT-ANALYSE
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("Bout-Analyse\n")
cat("════════════════════════════════════\n\n")

# ── Bout Criterion via Log-Survivor-Plot ──
# Pausenlängen zwischen zwei aktiven Minuten berechnen
cat("Berechne Inter-Bout-Intervalle (IBI)...\n")

bout_ibis <- dt_nacht_ok[order(igel, datum, stunde), {
  # Aktive Minuten identifizieren und Pausen berechnen
  aktiv_idx  <- which(aktiv == 1)
  if (length(aktiv_idx) < 2) {
    .(ibi = numeric(0))
  } else {
    pausen <- diff(aktiv_idx)   # Abstand in Minuten zwischen aktiven Minuten
    .(ibi = pausen[pausen > 1]) # Pausen > 1 Min = echte Unterbrechung
  }
}, by = .(igel, datum)]

# Log-Survivor-Plot der Pausenlängen
ibi_vals <- bout_ibis$ibi
ibi_sort <- sort(ibi_vals)
survivor <- (length(ibi_sort):1) / length(ibi_sort)

ibi_dt <- data.table(ibi = ibi_sort, survivor = survivor)
# Auf max. 60 Minuten Pause begrenzen (längere = sicher verschiedene Bouts)
ibi_plot <- ibi_dt[ibi <= 60]

p_logsurv <- ggplot(ibi_plot, aes(x = ibi, y = log(survivor))) +
  geom_line(color = "#2C3E6B", linewidth = 1) +
  geom_vline(xintercept = c(5, 10, 15, 20), linetype = "dashed",
             color = "grey60", alpha = 0.7) +
  labs(
    title    = "Log-survivor plot of inter-bout intervals",
    subtitle = "Inflection point = bout criterion (typically 5–20 minutes)\nVertical lines: 5, 10, 15, 20 min as possible criteria",
    x = "Pause between active minutes (min)",
    y = "log(survival probability)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

png(file.path(output_ordner, "bout_logsurvivor.png"),
    width = 1200, height = 700, res = 130)
print(p_logsurv)
dev.off()
cat("✓ Log-Survivor-Plot gespeichert: bout_logsurvivor.png\n")
cat("  → Schaue dir den Plot an und wähle das Bout-Kriterium!\n")
cat("  → Dann unten 'bout_kriterium' anpassen und ab Zeile 'BOUT-METRIKEN' erneut ausführen.\n\n")

# ──────────────────────────────────────────────────────────────
# BOUT-METRIKEN BERECHNEN
# (nach Sichtung des Log-Survivor-Plots anpassen)
# ──────────────────────────────────────────────────────────────

# !! HIER ANPASSEN nach Log-Survivor-Plot !!
bout_kriterium <- 10   # Minimale Pause in Minuten = zwei verschiedene Bouts

cat("Bout-Kriterium:", bout_kriterium, "Minuten\n")

# Bouts identifizieren und Metriken pro Nacht berechnen
berechne_bouts <- function(aktiv_vec, bout_krit) {
  n  <- length(aktiv_vec)
  if (n == 0 || sum(aktiv_vec) == 0) {
    # WICHTIG: alle Werte konsistent als numeric (nicht integer/NA-Mischung)
    return(list(n_bouts    = 0L,
                mean_dauer = NA_real_,
                max_dauer  = NA_real_,
                mean_ibi   = NA_real_))
  }
  # Run-Length-Encoding für aktive Blöcke
  rle_res  <- rle(aktiv_vec)
  laengen  <- as.numeric(rle_res$lengths)   # explizit numeric → verhindert integer/double-Konflikt
  werte    <- rle_res$values

  # Kurze Pausen (< bout_krit) mit benachbarten Bouts zusammenführen
  for (iter in 1:3) {
    passiv_idx <- which(werte == 0)
    zu_mergen  <- passiv_idx[laengen[passiv_idx] < bout_krit]
    if (length(zu_mergen) == 0) break
    j <- zu_mergen[1]
    if (j > 1 && j < length(werte)) {
      neue_laenge <- laengen[j-1] + laengen[j] + laengen[j+1]
      laengen <- c(laengen[seq_len(j-2)], neue_laenge, laengen[(j+2):length(laengen)])
      werte   <- c(werte[seq_len(j-2)], 1L, werte[(j+2):length(werte)])
    }
  }

  aktiv_bouts <- laengen[werte == 1]
  passiv_ibi  <- laengen[werte == 0]
  passiv_ibi  <- passiv_ibi[passiv_ibi >= bout_krit]

  list(
    n_bouts    = length(aktiv_bouts),                                         # integer OK
    mean_dauer = as.numeric(mean(aktiv_bouts)),                               # immer double
    max_dauer  = as.numeric(max(aktiv_bouts)),                                # immer double
    mean_ibi   = if (length(passiv_ibi) > 0) as.numeric(mean(passiv_ibi))    # immer double
                 else NA_real_
  )
}

cat("Berechne Bout-Metriken pro Nacht...\n")

bout_metriken <- dt_nacht_ok[order(igel, datum, stunde), {
  bm <- berechne_bouts(aktiv, bout_kriterium)
  .(
    n_bouts    = bm$n_bouts,
    mean_dauer = bm$mean_dauer,
    max_dauer  = bm$max_dauer,
    mean_ibi   = bm$mean_ibi,
    tage_seit  = tage_seit[1]
  )
}, by = .(igel, datum)]

cat("Bout-Metriken berechnet:", nrow(bout_metriken), "Nächte\n\n")

# Überblick
cat("── Bout-Metriken Übersicht ──\n")
print(bout_metriken[, .(
  n_bouts_med    = median(n_bouts, na.rm = TRUE),
  dauer_med_min  = median(mean_dauer, na.rm = TRUE),
  ibi_med_min    = median(mean_ibi, na.rm = TRUE)
)])

# ── GAMM: Verändern sich Bout-Metriken über Zeit? ──
cat("\n── GAMM auf Bout-Metriken ──\n")

# Anzahl Bouts pro Nacht (Poisson)
modell_bouts_n <- bam(
  n_bouts ~
    s(tage_seit, bs = "tp", k = 7) +
    s(tage_seit, igel, bs = "fs", k = 4, m = 1) +
    s(igel, bs = "re"),
  family = poisson(link = "log"),
  data   = bout_metriken[!is.na(n_bouts)],
  method = "fREML"
)

# Mittlere Boutdauer pro Nacht (Gamma — positiv und rechtsschief)
modell_bouts_dauer <- bam(
  mean_dauer ~
    s(tage_seit, bs = "tp", k = 7) +
    s(tage_seit, igel, bs = "fs", k = 4, m = 1) +
    s(igel, bs = "re"),
  family = Gamma(link = "log"),
  data   = bout_metriken[!is.na(mean_dauer) & mean_dauer > 0],
  method = "fREML"
)

cat("\nBout-Anzahl Modell:\n"); print(summary(modell_bouts_n))
cat("\nBout-Dauer Modell:\n");  print(summary(modell_bouts_dauer))

# ── Plots Bout-Metriken ──
plot_bout_smooth <- function(modell, bout_dt, y_var, titel, y_lab, transform = exp) {
  rng <- bout_dt[!is.na(get(y_var)), .(mn = min(tage_seit), mx = max(tage_seit))]
  nd  <- data.table(
    tage_seit = seq(rng$mn, rng$mx, length.out = 200),
    igel = factor(levels(bout_dt$igel)[1], levels = levels(bout_dt$igel))
  )
  pr <- predict(modell, newdata = nd,
                exclude = c("s(igel)", "s(tage_seit,igel)"),
                se.fit = TRUE, type = "link")
  nd[, fit := transform(pr$fit)]
  nd[, lwr := transform(pr$fit - 1.96 * pr$se.fit)]
  nd[, upr := transform(pr$fit + 1.96 * pr$se.fit)]

  ggplot() +
    geom_point(data = bout_dt[!is.na(get(y_var))],
               aes(x = tage_seit, y = get(y_var), color = igel),
               alpha = 0.3, size = 1.5) +
    geom_ribbon(data = nd,
                aes(x = tage_seit, ymin = lwr, ymax = upr),
                fill = "#2C3E6B", alpha = 0.2) +
    geom_line(data = nd,
              aes(x = tage_seit, y = fit),
              color = "#2C3E6B", linewidth = 1.3) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title = titel,
         x = "Days since release", y = y_lab, color = "Hedgehog") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
}

p_n_bouts <- plot_bout_smooth(
  modell_bouts_n, bout_metriken, "n_bouts",
  "Number of activity bouts per night",
  "Number of bouts")

p_dauer_bouts <- plot_bout_smooth(
  modell_bouts_dauer, bout_metriken[mean_dauer > 0], "mean_dauer",
  "Mean bout duration per night",
  "Duration (minutes)")

png(file.path(output_ordner, "gamm_bout_metriken.png"),
    width = 1600, height = 1400, res = 130)
print(p_n_bouts / p_dauer_bouts)
dev.off()
cat("✓ Plot gespeichert: gamm_bout_metriken.png\n\n")

# ──────────────────────────────────────────────────────────────
# ERGEBNISTABELLE SPEICHERN
# ──────────────────────────────────────────────────────────────

cat("════════════════════════════════════\n")
cat("Ergebnisse\n")
cat("════════════════════════════════════\n\n")

# Modellvergleich über R² (approx.) und EDF
ergebnisse <- data.table(
  Modell      = c("A: Nachtanteil (Beta)", "B: Aktiv/passiv Nacht (Binomial)",
                  "Bout N (Poisson)", "Bout Dauer (Gamma)"),
  Formel      = c("pct_aktiv_nacht ~ s(tage) + s(tage,igel,fs) + s(igel,re)",
                  "aktiv ~ s(tage) + s(stunde,cc) + ti(tage,stunde) + s(igel,re)",
                  "n_bouts ~ s(tage) + s(tage,igel,fs) + s(igel,re)",
                  "mean_dauer ~ s(tage) + s(tage,igel,fs) + s(igel,re)"),
  Deviance_erkl = c(
    round(summary(modell_A)$dev.expl * 100, 1),
    round(summary(modell_B)$dev.expl * 100, 1),
    round(summary(modell_bouts_n)$dev.expl * 100, 1),
    round(summary(modell_bouts_dauer)$dev.expl * 100, 1)
  ),
  AIC = round(c(AIC(modell_A), AIC(modell_B),
                AIC(modell_bouts_n), AIC(modell_bouts_dauer)), 1)
)

cat("Modellgüte:\n")
print(ergebnisse)

# Als RDS speichern (für spätere Interpretation)
saveRDS(list(
  modell_A           = modell_A,
  modell_B           = modell_B,
  modell_bouts_n     = modell_bouts_n,
  modell_bouts_dauer = modell_bouts_dauer,
  bout_metriken      = bout_metriken,
  igel_einschluss    = igel_einschluss,
  igel_ausschluss    = igel_ausschluss,
  bout_kriterium     = bout_kriterium
), file.path(output_ordner, "gamm_modelle.rds"))

cat("\n✓ Alle Modelle gespeichert: gamm_modelle.rds\n")

# ──────────────────────────────────────────────────────────────
# AUTOMATISCHER ERGEBNISBERICHT (Word)
# ──────────────────────────────────────────────────────────────

cat("\nErstelle Word-Ergebnisbericht...\n")

if (!requireNamespace("officer",   quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
library(officer)
library(flextable)

# ── Hilfsfunktionen für officer ──

# Schöne Smooth-Term-Tabelle aus einem GAMM
smooth_tabelle <- function(mod, titel) {
  s  <- summary(mod)
  st <- as.data.frame(s$s.table)
  st <- cbind(Term = rownames(st), round(st, 3))
  rownames(st) <- NULL
  names(st)[names(st) == "p-value"] <- "p-Wert"
  names(st)[names(st) == "edf"]     <- "EDF"
  # p-Wert-Spalte als Text formatieren
  st[["p-Wert"]] <- ifelse(st[["p-Wert"]] < 0.001, "< 0.001",
                    ifelse(st[["p-Wert"]] < 0.01,  "< 0.01",
                    ifelse(st[["p-Wert"]] < 0.05,  "< 0.05",
                           as.character(round(st[["p-Wert"]], 3)))))
  ft <- flextable(st) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#2C3E6B") |>
    color(part = "header", color = "white") |>
    bg(i = seq(1, nrow(st), 2), bg = "#F0F4F8") |>
    border_outer(part = "all", border = fp_border(color = "#CCCCCC", width = 1)) |>
    fontsize(size = 9, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    autofit()
  ft
}

# Kennzahlen-Box als Tabelle
kennzahlen_tabelle <- function(mod, family_name) {
  s   <- summary(mod)
  dev <- round(s$dev.expl * 100, 1)
  n   <- nrow(mod$model)
  rsq <- if (!is.null(s$r.sq)) round(s$r.sq, 3) else NA
  df_tab <- data.frame(
    Kennzahl = c("Stichprobengröße (n)", "Erklärte Devianz", "R² (adj.)", "AIC",
                 "Verteilung", "Schätzmethode"),
    Wert     = c(format(n, big.mark = "'"),
                 paste0(dev, " %"),
                 ifelse(is.na(rsq), "—", rsq),
                 round(AIC(mod), 1),
                 family_name,
                 "fREML")
  )
  flextable(df_tab) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#2C3E6B") |>
    color(part = "header", color = "white") |>
    bg(i = seq(1, nrow(df_tab), 2), bg = "#F0F4F8") |>
    fontsize(size = 9, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    autofit()
}

# Bild einfügen (mit Fehlertoleranz)
add_plot_safe <- function(doc, pfad, breite = 15, hoehe = 11) {
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

# Abschnittsüberschrift mit blauer Trennlinie
add_section <- function(doc, titel, level = 1) {
  style_name <- if (level == 1) "heading 1" else "heading 2"
  body_add_par(doc, titel, style = style_name)
}

# ── Dokument aufbauen ──

doc <- read_docx()

# Seitenränder
doc <- doc |>
  body_set_default_section(
    prop_section(
      page_margins = page_mar(top = 2, bottom = 2, left = 3, right = 2.5),
      page_size    = page_size(width = 21 / 2.54, height = 29.7 / 2.54)
    )
  )

# ── TITELSEITE ──
doc <- doc |>
  body_add_par("Igelbesenderung Niedersachsen", style = "heading 1") |>
  body_add_par("GAMM-Ergebnisbericht — Aktivitätsanalyse nach Auswilderung",
               style = "heading 2") |>
  body_add_par(paste("Erstellt:", format(Sys.time(), "%d.%m.%Y %H:%M")),
               style = "Normal") |>
  body_add_par(paste("Bout-Kriterium:", bout_kriterium, "Minuten"),
               style = "Normal") |>
  body_add_break()

# ── 1. STICHPROBE ──
doc <- add_section(doc, "1  Stichprobe und Beobachtungsabdeckung")

doc <- body_add_par(doc, paste0(
  "ZEITFENSTER-CUTOFF (N-at-risk-Schwelle): Alle Modelle (A, B, Bout-Analyse) ",
  "wurden auf den Zeitraum beschränkt, in dem mindestens ", n_min_igel,
  " Igel Daten liefern. Der automatisch berechnete Cutoff liegt bei Tag ",
  cutoff_tag, ". Danach ist kein sinnvoller Populationsschluss mehr möglich, ",
  "da späte Zeitpunkte von nur noch wenigen Langzeitbeobachteten dominiert ",
  "werden (Survivorship Bias). Das effektive Analysefenster umfasst Tag ",
  dt_tage_ok[, min(tage_seit)], " bis Tag ", cutoff_tag, " mit ",
  nlevels(dt_tage_ok$igel), " Individuen. ",
  "Zusätzlich wurden Beobachtungsgewichte eingesetzt ",
  "(kürzer beobachtete Igel = höheres Gewicht pro Datenpunkt)."),
  style = "Normal")

doc <- add_section(doc, "Gantt-Chart: Beobachtungsfenster je Igel", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "abdeckung_gantt.png"),
                     breite = 16, hoehe = 10)

doc <- add_section(doc, "N-at-risk pro Woche", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "abdeckung_n_risk.png"),
                     breite = 14, hoehe = 7)

# N-at-risk Tabelle
tbl_risk <- as.data.frame(n_risk_dt)
tbl_risk$tage_seit <- paste0("Tag ", tbl_risk$tage_seit)
names(tbl_risk) <- c("Zeitpunkt", "N Igel")
ft_risk <- flextable(tbl_risk) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#2C3E6B") |>
  color(part = "header", color = "white") |>
  bg(i = which(tbl_risk$`N Igel` < 5), bg = "#FFE4E1") |>  # rot = < 5 Igel
  fontsize(size = 9, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  autofit()
doc <- body_add_par(doc, "Rot markiert = weniger als 5 Igel → Populationsaussagen eingeschränkt.", style = "Normal")
doc <- body_add_flextable(doc, ft_risk)

# Eingeschlossene Igel
doc <- body_add_par(doc,
  paste0("Für die GAMM-Analyse wurden ", nrow(igel_einschluss),
         " von ", nrow(igel_einschluss) + nrow(igel_ausschluss),
         " Igeln eingeschlossen. Ausschlusskriterien: ",
         "< ", min_naechte, " Nächte, < ", min_nacht_min,
         " absolute Nachtminuten oder < 1% Nachtaktivität."),
  style = "Normal")

# Einschluss-Tabelle
tbl_ein <- igel_einschluss[, .(
  igel,
  `Nächte`           = n_naechte,
  `Nachtminuten`     = n_nacht_min,
  `% Nachtdaten`     = round(pct_nacht_dat, 1),
  `% aktiv (Nacht)`  = round(mean_aktiv_n,  1)
)]
ft_ein <- flextable(as.data.frame(tbl_ein)) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#2C3E6B") |>
  color(part = "header", color = "white") |>
  bg(i = seq(1, nrow(tbl_ein), 2), bg = "#F0F4F8") |>
  fontsize(size = 9, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  autofit()
doc <- body_add_flextable(doc, ft_ein)

if (nrow(igel_ausschluss) > 0) {
  doc <- body_add_par(doc, "Ausgeschlossene Individuen:", style = "Normal")
  tbl_aus <- igel_ausschluss[, .(
    igel,
    Nächte         = n_naechte,
    Nachtminuten   = n_nacht_min,
    `% Nachtdaten` = round(pct_nacht_dat, 1),
    Grund          = ausschluss_grund)]
  ft_aus <- flextable(as.data.frame(tbl_aus)) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#8B0000") |>
    color(part = "header", color = "white") |>
    fontsize(size = 9, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    autofit()
  doc <- body_add_flextable(doc, ft_aus)
}

doc <- add_section(doc, "Explorative Rohdaten", level = 2)
doc <- body_add_par(doc,
  paste0("Die folgenden Plots zeigen die Nachtaktivität der einzelnen Igel ",
         "vor jeder Modellierung — so lässt sich direkt ablesen, ob ein ",
         "Trend sichtbar ist oder ob die Tiere von Anfang an stabil nachtaktiv sind."),
  style = "Normal")

# Deskriptive Tabelle
tbl_desc <- as.data.frame(desc_igel)
ft_desc <- flextable(tbl_desc) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#2C3E6B") |>
  color(part = "header", color = "white") |>
  bg(i = seq(1, nrow(tbl_desc), 2), bg = "#F0F4F8") |>
  fontsize(size = 9, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  autofit()
doc <- body_add_flextable(doc, ft_desc)

# Früh-/Spätphase Vergleich
if (all(c("Früh (≤14d)", "Spät (>14d)") %in% frueh_spaet$phase)) {
  wt <- wilcox.test(
    nacht_roh[tage_seit <= grenze,  pct_aktiv_pct],
    nacht_roh[tage_seit  > grenze,  pct_aktiv_pct],
    exact = FALSE
  )
  doc <- body_add_par(doc,
    paste0("Früh- vs. Spätphase (Wilcoxon-Test): W = ",
           round(wt$statistic, 1), ", p = ", round(wt$p.value, 3),
           if (wt$p.value < 0.05)
             " → Signifikanter Unterschied zwischen erster und späterer Phase."
           else
             " → Kein signifikanter Unterschied — stabile Nachtaktivität von Beginn an."),
    style = "Normal")
}

doc <- add_section(doc, "Spaghetti-Plot: Nachtaktivität aller Igel", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "explorativ_spaghetti.png"),
                     breite = 16, hoehe = 9)
doc <- add_section(doc, "Facetten-Plot: Jeder Igel einzeln", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "explorativ_facetten.png"),
                     breite = 18, hoehe = 14)

doc <- body_add_break(doc)

# ── 1a. ERSTE-TAGE-ANALYSE ──
doc <- add_section(doc, "1a  Erste Tage nach Auswilderung (Tag 1–5)")
doc <- body_add_par(doc, paste0(
  "Die isolierte Betrachtung der ersten ", erste_tage_max, " Tage hat einen ",
  "entscheidenden methodischen Vorteil: In diesem Zeitfenster sind fast alle ",
  "Igel noch im Datensatz — es gibt nahezu keinen Survivorship Bias. ",
  "Unterschiede die hier sichtbar sind, spiegeln echte individuelle Variation ",
  "und eventuelle Anpassungseffekte wider, nicht Selektionseffekte durch ",
  "unterschiedliche Beobachtungsdauer. ",
  "Als Maß wurde der Nocturnal Index (NI = Nachtaktivität% − Tagaktivität%) ",
  "verwendet: positiv = nachtaktiver als tagaktiv."),
  style = "Normal")

# Deskriptive Tabelle
ft_et_desc <- flextable(as.data.frame(desc_et)) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#2C3E6B") |>
  color(part = "header", color = "white") |>
  bg(i = seq(1, nrow(desc_et), 2), bg = "#F0F4F8") |>
  fontsize(size = 9, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  autofit()
doc <- body_add_flextable(doc, ft_et_desc)

# Kruskal-Wallis-Ergebnis
if (exists("kw")) {
  doc <- body_add_par(doc, paste0(
    "Kruskal-Wallis-Test (Tag 1–", erste_tage_max, "): ",
    "H(", kw$parameter, ") = ", round(kw$statistic, 2),
    ", p = ", round(kw$p.value, 3),
    if (kw$p.value < 0.05)
      " → Signifikant: Der Nocturnal Index unterscheidet sich zwischen den ersten Tagen."
    else
      " → Nicht signifikant: Kein messbarer Unterschied zwischen Tag 1 und Tag 5 — die Tiere sind von Beginn an nachtaktiv."),
    style = "Normal")
}

doc <- add_section(doc, "Heatmap: Nocturnal Index je Igel und Tag", level = 2)
doc <- body_add_par(doc,
  "Jede Zelle zeigt den Nocturnal Index eines Igels an einem bestimmten Tag. Blau = nachtaktiv, Rot = tagaktiv. Leere Felder = kein Datenpunkt für diesen Tag.",
  style = "Normal")
doc <- add_plot_safe(doc, file.path(output_ordner, "erste_tage_heatmap.png"),
                     breite = 14, hoehe = 10)

doc <- add_section(doc, "Boxplot: Verteilung des Nocturnal Index pro Tag", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "erste_tage_boxplot.png"),
                     breite = 14, hoehe = 9)

doc <- add_section(doc, "Individuelle Verläufe: Tag 1–5", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "erste_tage_spaghetti.png"),
                     breite = 14, hoehe = 9)

doc <- body_add_break(doc)

# ── 1b. NOCTURNAL INDEX & FRÜHPHASEN-ANALYSE ──
doc <- add_section(doc, "1b  Nocturnal Index und Frühphasen-Analyse")
doc <- body_add_par(doc, paste0(
  "Um die Frage nach einer Anpassungszeit feiner zu beantworten, wurde zusätzlich ",
  "der Nocturnal Index (NI = Nachtaktivität% − Tagaktivität%) als Maß für die ",
  "relative Nocturnalität berechnet. Ein positiver NI bedeutet: der Igel ist ",
  "nachtaktiver als tagaktiver. Diese Größe kontrolliert für individuelle ",
  "Unterschiede im Gesamtaktivitätsniveau und ist unabhängig von der Tageslänge. ",
  "Zusätzlich wurde ein gesondertes GAMM auf die ersten 21 Tage (Frühphase) ",
  "mit logarithmischer Zeitachse gefittet — diese Skala streckt die ersten Tage ",
  "auseinander und ist sensitiver für kurzfristige Anpassungsprozesse."),
  style = "Normal")

doc <- add_section(doc, "Nocturnal Index — Populationstrend", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "ni_spaghetti.png"),
                     breite = 16, hoehe = 9)

doc <- add_section(doc, "Nocturnal Index — Log-Zeitachse (Frühe Tage)", level = 2)
doc <- body_add_par(doc, paste0(
  "Die logarithmische Zeitachse streckt die ersten Tage nach der Auswilderung ",
  "auseinander. So wird sichtbar, ob sich der Nocturnal Index direkt nach ",
  "der Auswilderung verändert — ein Effekt, der auf der linearen Zeitachse ",
  "oft nicht erkennbar wäre."),
  style = "Normal")
doc <- add_plot_safe(doc, file.path(output_ordner, "ni_log_zeitachse.png"),
                     breite = 16, hoehe = 9)

doc <- add_section(doc, "Nocturnal Index — Individuelle Verläufe (Log-Zeitachse)", level = 2)
doc <- add_plot_safe(doc, file.path(output_ordner, "ni_facetten_log.png"),
                     breite = 18, hoehe = 14)

if (!is.null(modell_frueh)) {
  doc <- add_section(doc, "Frühphasen-GAMM (Tage 1–21, log-Zeit)", level = 2)
  doc <- body_add_par(doc, paste0(
    "Beta-GAMM auf Nocturnal Index, nur Tage 1–21, mit log(tage_seit + 1) als ",
    "Prädiktor. Dieser Ansatz ist gezielt sensitiv für kurzfristige Adaptation ",
    "in den ersten Wochen nach der Auswilderung."),
    style = "Normal")
  doc <- body_add_flextable(doc, kennzahlen_tabelle(modell_frueh, "Beta (logit), NI Tage 1–21"))
  doc <- body_add_flextable(doc, smooth_tabelle(modell_frueh, "Frühphasen-Modell"))
  doc <- add_plot_safe(doc, file.path(output_ordner, "gamm_fruehphase.png"),
                       breite = 16, hoehe = 9)
}

doc <- body_add_break(doc)

# ── 2. MODELL A — NACHTAKTIVITÄT ÜBER ZEIT ──
doc <- add_section(doc, "2  Modell A — Nachtaktivität über die Zeit")
doc <- body_add_par(doc,
  paste0("Beta-GAMM (Tagesebene). Response: transformierter Anteil aktiver Nachtminuten. ",
         "n = ", nrow(nacht_tage), " Beobachtungen (Igel × Nacht)."),
  style = "Normal")

doc <- add_section(doc, "Modellgüte", level = 2)
doc <- body_add_flextable(doc, kennzahlen_tabelle(modell_A, "Beta (logit), gewichtet"))

doc <- add_section(doc, "Smooth Terms (Signifikanz nichtlinearer Effekte)", level = 2)
doc <- body_add_par(doc,
  paste0("EDF = Estimated Degrees of Freedom. EDF ≈ 1 = näherungsweise linear; ",
         "EDF > 1 = nichtlinear. Ein signifikanter s(tage_seit)-Term bedeutet, ",
         "dass sich die Nachtaktivität systematisch über die Zeit verändert."),
  style = "Normal")
doc <- body_add_flextable(doc, smooth_tabelle(modell_A, "Modell A (gewichtet)"))

doc <- add_section(doc, "Sensitivitätscheck: Survivorship-Bias-Korrektur", level = 2)
{
  p_gew    <- summary(modell_A)$s.table
  p_ungew  <- summary(modell_A_ungewichtet)$s.table
  idx_g    <- grep("tage_seit\\)", rownames(p_gew),  fixed = FALSE)[1]
  idx_u    <- grep("tage_seit\\)", rownames(p_ungew), fixed = FALSE)[1]
  pval_g   <- if(!is.na(idx_g))  p_gew[idx_g,  "p-value"] else NA
  pval_u   <- if(!is.na(idx_u))  p_ungew[idx_u, "p-value"] else NA

  sens_text <- paste0(
    "Das Modell wurde mit und ohne Beobachtungsgewichte gefittet. ",
    "Gewichtet (Survivorship-Korrektur): s(tage_seit) p = ",
    ifelse(is.na(pval_g), "n.v.", ifelse(pval_g < 0.001, "< 0.001", round(pval_g, 3))),
    ". Ungewichtet: p = ",
    ifelse(is.na(pval_u), "n.v.", ifelse(pval_u < 0.001, "< 0.001", round(pval_u, 3))),
    ". ",
    if (!is.na(pval_g) && !is.na(pval_u) &&
        ((pval_g < 0.05) == (pval_u < 0.05))) {
      "Beide Modelle kommen zum gleichen Signifikanzergebnis — der Survivorship Bias hat keinen qualitativ entscheidenden Einfluss auf die Schlussfolgerung."
    } else {
      "ACHTUNG: Die Modelle kommen zu verschiedenen Signifikanzergebnissen! Der Survivorship Bias beeinflusst das Ergebnis substanziell — Interpretation mit Vorsicht."
    }
  )
  doc <- body_add_par(doc, sens_text, style = "Normal")
}

doc <- add_section(doc, "Abbildung: Populationstrend und individuelle Kurven", level = 2)
doc <- add_plot_safe(doc,
  file.path(output_ordner, "gamm_A_nachtaktivitaet.png"),
  breite = 16, hoehe = 14)

doc <- add_section(doc, "Modelldiagnostik (gam.check)", level = 2)
doc <- add_plot_safe(doc,
  file.path(output_ordner, "gamm_A_diagnostik.png"),
  breite = 16, hoehe = 10)

doc <- body_add_break(doc)

# ── 3. MODELL B — 24H-RHYTHMUS ──
doc <- add_section(doc, "3  Modell B — 24h-Rhythmus und Verschiebung über Zeit")
doc <- body_add_par(doc,
  paste0("Binomial-GAMM (Minutenebene, nur Nacht). Response: aktiv (1) / passiv (0). ",
         "n = ", format(nrow(dt_nacht_ok), big.mark = "'"), " Nachtminuten. ",
         "Der Tensor-Interaktionsterm ti(tage_seit, stunde) zeigt, ob sich der ",
         "Aktivitätszeitpunkt innerhalb der Nacht über die Zeit verschiebt."),
  style = "Normal")

doc <- add_section(doc, "Modellgüte", level = 2)
doc <- body_add_flextable(doc, kennzahlen_tabelle(modell_B, "Binomial (logit)"))

doc <- add_section(doc, "Smooth Terms", level = 2)
doc <- body_add_flextable(doc, smooth_tabelle(modell_B, "Modell B"))
doc <- body_add_par(doc,
  paste0("Ein signifikanter ti(tage_seit,stunde)-Term bedeutet: Der 24h-Rhythmus ",
         "verändert sich signifikant über die Zeit nach der Auswilderung."),
  style = "Normal")

doc <- add_section(doc, "Abbildung: 24h-Heatmap und Zeittrend", level = 2)
doc <- add_plot_safe(doc,
  file.path(output_ordner, "gamm_B_rhythmus.png"),
  breite = 16, hoehe = 14)

doc <- body_add_break(doc)

# ── 4. BOUT-ANALYSE ──
doc <- add_section(doc, "4  Bout-Analyse")
doc <- body_add_par(doc,
  paste0("Bout-Kriterium (aus Log-Survivor-Plot): ", bout_kriterium, " Minuten. ",
         "Berechnet über ", nrow(bout_metriken), " Nächte (", nlevels(bout_metriken$igel),
         " Igel). Mittlere Boutanzahl pro Nacht: ",
         round(mean(bout_metriken$n_bouts, na.rm = TRUE), 1),
         " | Mittlere Boutdauer: ",
         round(mean(bout_metriken$mean_dauer, na.rm = TRUE), 1), " Min."),
  style = "Normal")

doc <- add_section(doc, "Log-Survivor-Plot (Bout-Kriterium)", level = 2)
doc <- add_plot_safe(doc,
  file.path(output_ordner, "bout_logsurvivor.png"),
  breite = 14, hoehe = 8)

doc <- add_section(doc, "Bout-Anzahl über Zeit (Poisson-GAMM)", level = 2)
doc <- body_add_flextable(doc, kennzahlen_tabelle(modell_bouts_n, "Poisson (log)"))
doc <- body_add_flextable(doc, smooth_tabelle(modell_bouts_n, "Bout-Anzahl"))

doc <- add_section(doc, "Boutdauer über Zeit (Gamma-GAMM)", level = 2)
doc <- body_add_flextable(doc, kennzahlen_tabelle(modell_bouts_dauer, "Gamma (log)"))
doc <- body_add_flextable(doc, smooth_tabelle(modell_bouts_dauer, "Bout-Dauer"))

doc <- add_section(doc, "Abbildung: Bout-Metriken über Zeit", level = 2)
doc <- add_plot_safe(doc,
  file.path(output_ordner, "gamm_bout_metriken.png"),
  breite = 16, hoehe = 14)

doc <- body_add_break(doc)

# ── 5. MODELLVERGLEICH ──
doc <- add_section(doc, "5  Modellvergleich")
ft_erg <- flextable(as.data.frame(ergebnisse)) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#2C3E6B") |>
  color(part = "header", color = "white") |>
  bg(i = seq(1, nrow(ergebnisse), 2), bg = "#F0F4F8") |>
  fontsize(size = 9, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  autofit()
doc <- body_add_flextable(doc, ft_erg)

# ── 6. BIOLOGISCHE INTERPRETATION ──
doc <- body_add_break(doc)
doc <- add_section(doc, "6  Biologische Interpretation der Ergebnisse")

# Hilfsfunktion: p-Wert lesbar aus Smooth-Tabelle extrahieren
get_p <- function(mod, term) {
  st <- summary(mod)$s.table
  rn <- rownames(st)
  idx <- grep(term, rn, fixed = TRUE)
  if (length(idx) == 0) return(NA_real_)
  st[idx[1], "p-value"]
}
get_edf <- function(mod, term) {
  st <- summary(mod)$s.table
  rn <- rownames(st)
  idx <- grep(term, rn, fixed = TRUE)
  if (length(idx) == 0) return(NA_real_)
  round(st[idx[1], "edf"], 1)
}

# ── Modell A Interpretation
p_A_zeit  <- get_p(modell_A,  "tage_seit)")   # Populationstrend
p_A_igel  <- get_p(modell_A,  "igel)")        # Zufälliger Effekt Igel
sig_A_zeit <- !is.na(p_A_zeit) && p_A_zeit < 0.05
sig_A_igel <- !is.na(p_A_igel) && p_A_igel < 0.05

interp_A <- if (sig_A_zeit) {
  paste0(
    "Modell A zeigt einen SIGNIFIKANTEN Zeittrend (p = ",
    ifelse(p_A_zeit < 0.001, "< 0.001", round(p_A_zeit, 3)),
    "): Die Nachtaktivität der Igel verändert sich systematisch ",
    "in den Wochen nach der Auswilderung. ")
} else {
  paste0(
    "Modell A zeigt KEINEN signifikanten Zeittrend (p = ",
    ifelse(is.na(p_A_zeit), "n.v.", round(p_A_zeit, 3)),
    "). Das bedeutet: Die Nachtaktivität der Igel bleibt ",
    "über die gesamte Beobachtungszeit stabil — es gibt keine ",
    "systematische Zu- oder Abnahme nach der Auswilderung. ",
    "Dieses Ergebnis ist biologisch bedeutsam: Die Tiere zeigen ",
    "sofort nach der Auswilderung ihr artgemäßes nächtliches ",
    "Aktivitätsmuster und passen sich nicht erst ein. ")
}
interp_A <- paste0(interp_A,
  if (sig_A_igel) {
    paste0("Allerdings unterscheiden sich die Individuen ",
           "signifikant in ihrem Aktivitätsniveau (p < 0.001): ",
           "Manche Igel sind generell aktiver als andere, ",
           "unabhängig von der Zeit nach der Auswilderung.")
  } else {
    "Auch zwischen den Individuen gibt es keine signifikanten Unterschiede."
  })

doc <- add_section(doc, "Modell A: Nachtaktivität über Zeit", level = 2)
doc <- body_add_par(doc, interp_A, style = "Normal")

# ── Modell B Interpretation
p_B_tage  <- get_p(modell_B, "tage_seit)")
p_B_stunde <- get_p(modell_B, "stunde)")
p_B_inter <- get_p(modell_B, "ti(tage_seit,stunde)")
sig_B_inter <- !is.na(p_B_inter) && p_B_inter < 0.05
sig_B_stunde <- !is.na(p_B_stunde) && p_B_stunde < 0.05

interp_B1 <- if (sig_B_stunde) {
  "Modell B bestätigt: Innerhalb der Nacht gibt es einen signifikanten zeitlichen Rhythmus — die Igel sind nicht gleichmäßig die ganze Nacht aktiv, sondern zeigen Präferenzen für bestimmte Nachtphasen (z.B. Aktivitätsspitzen kurz nach Einbruch der Dunkelheit)."
} else {
  "Modell B zeigt keinen klar strukturierten Aktivitätszeitpunkt innerhalb der Nacht — die Aktivität ist über die Nachtstunden relativ gleichmäßig verteilt."
}

interp_B2 <- if (sig_B_inter) {
  paste0(
    "Der Interaktionsterm ti(tage_seit, stunde) ist SIGNIFIKANT (p = ",
    ifelse(p_B_inter < 0.001, "< 0.001", round(p_B_inter, 3)),
    "): Der zeitliche Rhythmus innerhalb der Nacht VERÄNDERT sich über die ",
    "Wochen nach der Auswilderung. Die Igel verschieben also den Zeitpunkt ",
    "ihrer Nachtaktivität — ein Zeichen für Habitateingewöhnung oder ",
    "saisonale Anpassung.")
} else {
  paste0(
    "Der Interaktionsterm ti(tage_seit, stunde) ist NICHT signifikant (p = ",
    ifelse(is.na(p_B_inter), "n.v.", round(p_B_inter, 3)),
    "): Der Aktivitätszeitpunkt innerhalb der Nacht bleibt stabil — ",
    "die Igel zeigen direkt nach der Auswilderung dasselbe Aktivitätsmuster ",
    "innerhalb der Nacht wie nach Wochen. Kein Eingewöhnungseffekt nachweisbar.")
}

doc <- add_section(doc, "Modell B: 24h-Rhythmus und zeitliche Verschiebung", level = 2)
doc <- body_add_par(doc, interp_B1, style = "Normal")
doc <- body_add_par(doc, interp_B2, style = "Normal")

# ── Bout-Analyse Interpretation
p_bn  <- get_p(modell_bouts_n,     "tage_seit)")
p_bd  <- get_p(modell_bouts_dauer, "tage_seit)")
edf_bn <- get_edf(modell_bouts_n,  "s(tage_seit)")
sig_bn <- !is.na(p_bn) && p_bn < 0.05
sig_bd <- !is.na(p_bd) && p_bd < 0.05

med_n  <- round(median(bout_metriken$n_bouts,    na.rm = TRUE), 1)
med_d  <- round(median(bout_metriken$mean_dauer, na.rm = TRUE), 1)

interp_bout <- paste0(
  "Im Median zeigen die Igel ", med_n, " Aktivitätsbouts pro Nacht ",
  "mit einer mittleren Boutdauer von ", med_d, " Minuten. ")

interp_bout <- paste0(interp_bout,
  if (sig_bn) {
    paste0(
      "Die ANZAHL der Bouts pro Nacht verändert sich signifikant über die Zeit ",
      "(p = ", ifelse(p_bn < 0.001, "< 0.001", round(p_bn, 3)),
      ", EDF = ", edf_bn, "). Ein EDF > 1 deutet auf einen nicht-linearen Verlauf hin: ",
      "Die Igel verändern nicht einfach gradlinig ihre Boutanzahl, sondern es gibt ",
      "eine komplexere zeitliche Dynamik — z.B. eine anfängliche Anpassungsphase. ",
      "Dies ist der auffälligste Befund der Bout-Analyse.")
  } else {
    paste0("Die Anzahl der Bouts pro Nacht bleibt über die Zeit stabil ",
           "(p = ", ifelse(is.na(p_bn), "n.v.", round(p_bn, 3)), ").")
  })

interp_bout <- paste0(interp_bout, " ")

interp_bout <- paste0(interp_bout,
  if (sig_bd) {
    paste0("Auch die DAUER der Bouts verändert sich signifikant über die Zeit ",
           "(p = ", ifelse(p_bd < 0.001, "< 0.001", round(p_bd, 3)), ").")
  } else {
    paste0("Die Boutdauer bleibt dagegen stabil (p = ",
           ifelse(is.na(p_bd), "n.v.", round(p_bd, 3)),
           "): Die einzelnen Aktivitätsphasen werden weder länger noch kürzer.")
  })

doc <- add_section(doc, "Bout-Analyse: Struktur der Aktivitätsphasen", level = 2)
doc <- body_add_par(doc, interp_bout, style = "Normal")

# ── Gesamtfazit
fazit_teile <- c(
  "GESAMTFAZIT:",
  if (!sig_A_zeit && sig_A_igel)
    "Die Nachtaktivität der Igel ist direkt nach der Auswilderung stabil und verändert sich nicht systematisch über die Zeit — ein Hinweis auf sofortige Habituation. Individuelle Unterschiede im Aktivitätsniveau sind jedoch ausgeprägt."
  else if (sig_A_zeit)
    "Die Nachtaktivität verändert sich über die Zeit — weitere Interpretation aus Plot A ablesen."
  else
    "Keine signifikanten Effekte in Modell A nachgewiesen.",
  if (sig_B_inter)
    "Die Lage der Aktivitätsphasen innerhalb der Nacht verschiebt sich signifikant."
  else
    "Der 24h-Rhythmus bleibt innerhalb der Nacht stabil.",
  if (sig_bn)
    "Die Anzahl der Aktivitätsbouts pro Nacht verändert sich nonlinear — dies ist die deutlichste zeitliche Dynamik im Datensatz."
  else
    "Keine zeitliche Veränderung in der Boutstruktur nachgewiesen."
)
fazit_text <- paste(fazit_teile, collapse = " ")

doc <- add_section(doc, "Gesamtfazit", level = 2)
doc <- body_add_par(doc, fazit_text, style = "Normal")

# ── Speichern ──
bericht_pfad <- file.path(output_ordner, "Block2_GAMM_Ergebnisbericht.docx")
print(doc, target = bericht_pfad)
cat("\n✓ Word-Ergebnisbericht gespeichert:", basename(bericht_pfad), "\n")

cat("\n── Erzeugte Dateien ──\n")
cat("  gamm_A_nachtaktivitaet.png  — Populationstrend + individuelle Kurven\n")
cat("  gamm_A_diagnostik.png       — gam.check() für Modell A\n")
cat("  gamm_B_rhythmus.png         — 2D-Heatmap Tage × Uhrzeit\n")
cat("  gamm_B_diagnostik.png       — gam.check() für Modell B\n")
cat("  bout_logsurvivor.png        — Bout-Kriterium ablesen\n")
cat("  gamm_bout_metriken.png      — Bout-Anzahl und Boutdauer über Zeit\n")
cat("  gamm_modelle.rds            — alle Modellobjekte\n")
cat("  GAMM_Ergebnisbericht.docx   — vollständiger Word-Bericht\n")
cat("\n✓ GAMM-Analyse abgeschlossen!\n")



# Modellzusammenfassung speichern:
sink(file.path(output_ordner, "gamm_kennzahlen.txt"))
summary(modell_A)
summary(modell_B)
sink()
# ══════════════════════════════════════════════════════════════
# EXCEL OVERVIEW  (Block2_Uebersicht.xlsx)
# ══════════════════════════════════════════════════════════════

cat("\n════════════════════════════════════\n")
cat("Excel overview (Block2_Uebersicht.xlsx)\n")
cat("════════════════════════════════════\n\n")

# ── Colour palette ─────────────────────────────────────────────
XC_DARK  <- "#1a2e4a"
XC_MID   <- "#2C5F8A"
XC_LIGHT <- "#D6E4F0"
XC_LBLUE <- "#EBF3FB"
XC_GREY  <- "#F5F5F5"
XC_WHITE <- "#FFFFFF"
XC_YELL  <- "#FFF9C4"
XC_RED   <- "#FFE0E0"
XC_GREEN <- "#E8F5E9"

# ── Style helpers ──────────────────────────────────────────────
xl_hdr <- createStyle(fontName="Arial", fontSize=10, fontColour=XC_WHITE,
  fgFill=XC_MID, halign="CENTER", valign="center",
  textDecoration="bold", wrapText=TRUE,
  border="TopBottomLeftRight", borderColour="#AABFD4")

xl_title <- createStyle(fontName="Arial", fontSize=14, fontColour=XC_WHITE,
  fgFill=XC_DARK, textDecoration="bold", halign="left", valign="center", indent=1)

xl_sub <- createStyle(fontName="Arial", fontSize=10, fontColour="#444444",
  fgFill=XC_LBLUE, halign="left", valign="center",
  textDecoration="italic", indent=1)

xl_sec <- createStyle(fontName="Arial", fontSize=11, fontColour=XC_WHITE,
  fgFill=XC_MID, textDecoration="bold", halign="left", valign="center", indent=1)

xl_cell <- function(bg=XC_WHITE, ha="center") {
  createStyle(fontName="Arial", fontSize=10, fgFill=bg,
    halign=ha, valign="center",
    border="TopBottomLeftRight", borderColour="#AABFD4", wrapText=TRUE)
}

xl_safe <- function(v) {
  v <- as.numeric(v[!is.na(v)])
  if (length(v) == 0) return(c(N=0,Min=NA,Max=NA,Median=NA,Mean=NA,SD=NA))
  c(N=length(v), Min=round(min(v),1), Max=round(max(v),1),
    Median=round(median(v),1), Mean=round(mean(v),1), SD=round(sd(v),1))
}

wb <- createWorkbook()

# ────────────────────────────────────────────────────────────────
# Sheet 1: Overview (Dashboard)
# ────────────────────────────────────────────────────────────────
addWorksheet(wb, "Overview", gridLines=FALSE)
ws <- "Overview"

# Title
mergeCells(wb, ws, cols=1:8, rows=1)
writeData(wb, ws,
  paste0("Block 2 — GAMM Analysis Overview  |  VHF Hedgehog Telemetry  |  Updated: ",
         format(Sys.Date(), "%d.%m.%Y")),
  startRow=1, startCol=1)
addStyle(wb, ws, xl_title, rows=1, cols=1:8, gridExpand=TRUE)
setRowHeights(wb, ws, rows=1, heights=28)

# Subtitle
mergeCells(wb, ws, cols=1:8, rows=2)
writeData(wb, ws,
  paste0("N = ", nlevels(dt_tage_ok$igel), " hedgehogs (GAMM model)  |  ",
         "Observation window: day ", dt_tage_ok[, min(tage_seit)],
         " – day ", cutoff_tag,
         "  |  Bout criterion: ", bout_kriterium, " min"),
  startRow=2, startCol=1)
addStyle(wb, ws, xl_sub, rows=2, cols=1:8, gridExpand=TRUE)
setRowHeights(wb, ws, rows=2, heights=18)
setRowHeights(wb, ws, rows=3, heights=6)

# ── Section: Inclusion / Exclusion ────────────────────────────
cur <- 4
mergeCells(wb, ws, cols=1:8, rows=cur)
writeData(wb, ws, "Sample inclusion / exclusion", startRow=cur, startCol=1)
addStyle(wb, ws, xl_sec, rows=cur, cols=1:8, gridExpand=TRUE)
setRowHeights(wb, ws, rows=cur, heights=20); cur <- cur+1

inc_hdrs <- c("","Included","Excluded","Threshold: nights",
              "Threshold: night-min","Mean obs. duration (d)")
writeData(wb, ws, as.data.frame(t(inc_hdrs)), startRow=cur, startCol=1, colNames=FALSE)
addStyle(wb, ws, xl_hdr, rows=cur, cols=1:6, gridExpand=TRUE)
setRowHeights(wb, ws, rows=cur, heights=18); cur <- cur+1

inc_row <- data.frame(
  Label    = "Hedgehogs",
  Included = nrow(igel_einschluss),
  Excluded = nrow(igel_ausschluss),
  Thr_N    = min_naechte,
  Thr_Min  = min_nacht_min,
  MeanObs  = round(mean(abdeckung$beob_dauer), 1)
)
writeData(wb, ws, inc_row, startRow=cur, startCol=1, colNames=FALSE)
addStyle(wb, ws, xl_cell(XC_WHITE,"center"), rows=cur, cols=1:6, gridExpand=TRUE)
setRowHeights(wb, ws, rows=cur, heights=16); cur <- cur+2

# ── Section: Model performance ─────────────────────────────────
mergeCells(wb, ws, cols=1:8, rows=cur)
writeData(wb, ws, "Model performance", startRow=cur, startCol=1)
addStyle(wb, ws, xl_sec, rows=cur, cols=1:8, gridExpand=TRUE)
setRowHeights(wb, ws, rows=cur, heights=20); cur <- cur+1

mod_hdrs <- c("Model","Response","Family","n","Dev. expl. (%)","AIC","p(time smooth)")
writeData(wb, ws, as.data.frame(t(mod_hdrs)), startRow=cur, startCol=1, colNames=FALSE)
addStyle(wb, ws, xl_hdr, rows=cur, cols=1:7, gridExpand=TRUE)
setRowHeights(wb, ws, rows=cur, heights=18); cur <- cur+1

get_p_xl <- function(mod, term) {
  st <- summary(mod)$s.table
  idx <- grep(term, rownames(st), fixed=TRUE)
  if (length(idx)==0) return(NA_real_)
  round(st[idx[1], "p-value"], 4)
}

mod_perf <- list(
  list("A: Night proportion (Beta GAMM)", "pct_aktiv_nacht", "Beta (logit), weighted",
       nrow(modell_A$model), round(summary(modell_A)$dev.expl*100,1), round(AIC(modell_A),1),
       get_p_xl(modell_A, "tage_seit)")),
  list("B: Active/passive minutes (Binomial)", "aktiv (0/1)", "Binomial (logit)",
       nrow(modell_B$model), round(summary(modell_B)$dev.expl*100,1), round(AIC(modell_B),1),
       get_p_xl(modell_B, "tage_seit)")),
  list("Bout count (Poisson GAMM)", "n_bouts", "Poisson (log)",
       nrow(modell_bouts_n$model), round(summary(modell_bouts_n)$dev.expl*100,1),
       round(AIC(modell_bouts_n),1), get_p_xl(modell_bouts_n, "tage_seit)")),
  list("Bout duration (Gamma GAMM)", "mean_dauer", "Gamma (log)",
       nrow(modell_bouts_dauer$model), round(summary(modell_bouts_dauer)$dev.expl*100,1),
       round(AIC(modell_bouts_dauer),1), get_p_xl(modell_bouts_dauer, "tage_seit)"))
)
if (!is.null(modell_frueh)) {
  mod_perf[[5]] <- list("Early-phase GAMM (Beta, days 1–21)", "Nocturnal Index (0–1)",
    "Beta (logit)", nrow(modell_frueh$model),
    round(summary(modell_frueh)$dev.expl*100,1), round(AIC(modell_frueh),1),
    get_p_xl(modell_frueh, "log_tage)"))
}

for (i in seq_along(mod_perf)) {
  m   <- mod_perf[[i]]
  bg  <- if (i %% 2 == 0) XC_GREY else XC_WHITE
  pv  <- m[[7]]
  pv_txt <- if (is.na(pv)) "—" else if (pv < 0.001) "< 0.001" else as.character(pv)
  rv <- data.frame(Model=m[[1]], Resp=m[[2]], Fam=m[[3]],
                   n=m[[4]], Dev=m[[5]], AIC=m[[6]], P=pv_txt)
  writeData(wb, ws, rv, startRow=cur, startCol=1, colNames=FALSE)
  # highlight if significant
  sig_bg <- if (!is.na(pv) && pv < 0.05) XC_GREEN else bg
  addStyle(wb, ws, xl_cell(bg,"left"),    rows=cur, cols=1:3, gridExpand=TRUE)
  addStyle(wb, ws, xl_cell(bg,"center"),  rows=cur, cols=4:6, gridExpand=TRUE)
  addStyle(wb, ws, xl_cell(sig_bg,"center"), rows=cur, cols=7, gridExpand=TRUE)
  setRowHeights(wb, ws, rows=cur, heights=16); cur <- cur+1
}
setRowHeights(wb, ws, rows=cur, heights=6); cur <- cur+1

# ── Section: Descriptive statistics (nocturnal activity) ───────
mergeCells(wb, ws, cols=1:8, rows=cur)
writeData(wb, ws, "Descriptive statistics — nocturnal activity (% active night minutes)",
          startRow=cur, startCol=1)
addStyle(wb, ws, xl_sec, rows=cur, cols=1:8, gridExpand=TRUE)
setRowHeights(wb, ws, rows=cur, heights=20); cur <- cur+1

stat_hdrs <- c("Variable","N","Min.","Max.","Median","Mean","SD","Unit")
writeData(wb, ws, as.data.frame(t(stat_hdrs)), startRow=cur, startCol=1, colNames=FALSE)
addStyle(wb, ws, xl_hdr, rows=cur, cols=1:8, gridExpand=TRUE)
setRowHeights(wb, ws, rows=cur, heights=18); cur <- cur+1

stat_vars_xl <- list(
  list("Nocturnal activity (% active night min)", nacht_roh$pct_aktiv_pct, "%"),
  list("Nocturnal Index",   ni_wide$nocturnal_index[!is.na(ni_wide$nocturnal_index)], "%"),
  list("Number of bouts / night", bout_metriken$n_bouts, "n"),
  list("Mean bout duration / night", bout_metriken$mean_dauer[!is.na(bout_metriken$mean_dauer)], "min"),
  list("Observation duration", abdeckung$beob_dauer, "days")
)
for (i in seq_along(stat_vars_xl)) {
  item <- stat_vars_xl[[i]]; s <- xl_safe(item[[2]])
  bg   <- if (i %% 2 == 0) XC_GREY else XC_WHITE
  rv   <- data.frame(Var=item[[1]], N=s["N"],
    Min=ifelse(is.na(s["Min"]),"—",as.character(s["Min"])),
    Max=ifelse(is.na(s["Max"]),"—",as.character(s["Max"])),
    Med=ifelse(is.na(s["Median"]),"—",as.character(s["Median"])),
    Mean=ifelse(is.na(s["Mean"]),"—",as.character(s["Mean"])),
    SD=ifelse(is.na(s["SD"]),"—",as.character(s["SD"])),
    Unit=item[[3]])
  writeData(wb, ws, rv, startRow=cur, startCol=1, colNames=FALSE)
  addStyle(wb, ws, xl_cell(bg,"left"),   rows=cur, cols=1, gridExpand=TRUE)
  addStyle(wb, ws, xl_cell(bg,"center"), rows=cur, cols=2:8, gridExpand=TRUE)
  setRowHeights(wb, ws, rows=cur, heights=16); cur <- cur+1
}
setColWidths(wb, ws, cols=1:8, widths=c(38,8,10,10,12,12,10,8))

# ────────────────────────────────────────────────────────────────
# Sheet 2: Included animals
# ────────────────────────────────────────────────────────────────
addWorksheet(wb, "Included animals", gridLines=FALSE)
ws2 <- "Included animals"
freezePane(wb, ws2, firstRow=TRUE)

# Merge abdeckung with igel_qc for full info
incl_tbl <- merge(igel_einschluss, abdeckung[, .(igel, t_start, t_ende, beob_dauer)],
                  by="igel", all.x=TRUE)
incl_tbl <- incl_tbl[order(beob_dauer, decreasing=TRUE)]

# Add mean nocturnal activity from nacht_roh
act_sum <- nacht_roh[, .(
  mean_pct_aktiv = round(mean(pct_aktiv_pct), 1),
  median_ni      = round(median(
    ni_wide$nocturnal_index[ni_wide$igel == igel], na.rm=TRUE
  ), 1)
), by=igel]
# safer merge via desc_igel
act_sum2 <- desc_igel[, .(igel, mean_aktiv=Mittelwert, median_aktiv=Median)]
incl_tbl <- merge(incl_tbl, act_sum2, by="igel", all.x=TRUE)

out_incl <- incl_tbl[, .(
  `Hedgehog`           = igel,
  `Nights (N)`         = n_naechte,
  `Night minutes`      = n_nacht_min,
  `% night data`       = round(pct_nacht_dat, 1),
  `Mean active % (night)` = round(mean_aktiv_n, 1),
  `Mean nocturnal act.`= mean_aktiv,
  `Obs. start (day)`   = t_start,
  `Obs. end (day)`     = t_ende,
  `Obs. duration (d)`  = beob_dauer
)]

writeData(wb, ws2, out_incl, startRow=1, startCol=1, headerStyle=xl_hdr)
setRowHeights(wb, ws2, rows=1, heights=22)
for (i in seq_len(nrow(out_incl))) {
  bg <- if (i %% 2 == 0) XC_GREY else XC_WHITE
  addStyle(wb, ws2, xl_cell(bg,"center"), rows=i+1, cols=1:9, gridExpand=TRUE)
  setRowHeights(wb, ws2, rows=i+1, heights=16)
}
setColWidths(wb, ws2, cols=1:9, widths=c(12,10,14,12,20,20,16,14,16))

# ────────────────────────────────────────────────────────────────
# Sheet 3: Excluded animals
# ────────────────────────────────────────────────────────────────
addWorksheet(wb, "Excluded animals", gridLines=FALSE)
ws3 <- "Excluded animals"

if (nrow(igel_ausschluss) > 0) {
  freezePane(wb, ws3, firstRow=TRUE)
  out_excl <- igel_ausschluss[, .(
    `Hedgehog`       = igel,
    `Nights (N)`     = n_naechte,
    `Night minutes`  = n_nacht_min,
    `% night data`   = round(pct_nacht_dat, 1),
    `Mean active %`  = round(mean_aktiv_n,  1),
    `Exclusion reason` = ausschluss_grund
  )]
  excl_hdr <- createStyle(fontName="Arial", fontSize=10, fontColour=XC_WHITE,
    fgFill="#8B2020", halign="CENTER", valign="center",
    textDecoration="bold", wrapText=TRUE,
    border="TopBottomLeftRight", borderColour="#AABFD4")
  writeData(wb, ws3, out_excl, startRow=1, startCol=1, headerStyle=excl_hdr)
  setRowHeights(wb, ws3, rows=1, heights=22)
  for (i in seq_len(nrow(out_excl))) {
    bg <- if (i %% 2 == 0) XC_GREY else XC_WHITE
    addStyle(wb, ws3, xl_cell(bg,"center"), rows=i+1, cols=1:5, gridExpand=TRUE)
    addStyle(wb, ws3, xl_cell(bg,"left"),   rows=i+1, cols=6,   gridExpand=TRUE)
    setRowHeights(wb, ws3, rows=i+1, heights=16)
  }
  setColWidths(wb, ws3, cols=1:6, widths=c(12,10,14,12,14,38))
} else {
  writeData(wb, ws3, data.frame(Note="No hedgehogs excluded."), startRow=1, startCol=1)
}

# ────────────────────────────────────────────────────────────────
# Sheet 4: Model results (smooth terms)
# ────────────────────────────────────────────────────────────────
addWorksheet(wb, "Model results", gridLines=FALSE)
ws4 <- "Model results"

mod_list_xl <- list(
  list(modell_A,           "Model A — Night proportion (Beta GAMM, weighted)"),
  list(modell_B,           "Model B — Active/passive minutes (Binomial GAMM)"),
  list(modell_bouts_n,     "Bout count model (Poisson GAMM)"),
  list(modell_bouts_dauer, "Bout duration model (Gamma GAMM)")
)
if (!is.null(modell_frueh))
  mod_list_xl[[5]] <- list(modell_frueh, "Early-phase GAMM (Beta, days 1–21, log scale)")

cur4 <- 1
for (ml in mod_list_xl) {
  mod <- ml[[1]]; lbl <- ml[[2]]
  # Section header
  mergeCells(wb, ws4, cols=1:6, rows=cur4)
  writeData(wb, ws4, lbl, startRow=cur4, startCol=1)
  addStyle(wb, ws4, xl_sec, rows=cur4, cols=1:6, gridExpand=TRUE)
  setRowHeights(wb, ws4, rows=cur4, heights=20); cur4 <- cur4+1
  # Subtitle: deviance + AIC
  s   <- summary(mod)
  dev <- round(s$dev.expl*100, 1)
  rsq <- if (!is.null(s$r.sq)) round(s$r.sq, 3) else NA
  mergeCells(wb, ws4, cols=1:6, rows=cur4)
  writeData(wb, ws4,
    paste0("Deviance explained: ", dev, "%  |  R² adj.: ",
           ifelse(is.na(rsq),"—",rsq), "  |  AIC: ", round(AIC(mod),1),
           "  |  n = ", format(nrow(mod$model), big.mark=",")),
    startRow=cur4, startCol=1)
  addStyle(wb, ws4, xl_sub, rows=cur4, cols=1:6, gridExpand=TRUE)
  setRowHeights(wb, ws4, rows=cur4, heights=16); cur4 <- cur4+1
  # Smooth term table
  # Beta/Binomial GAMMs → "F"; Poisson/Gamma GAMMs → "Chi.sq"
  st <- as.data.frame(s$s.table)
  stat_col   <- if ("F"      %in% names(st)) "F"      else
                if ("Chi.sq" %in% names(st)) "Chi.sq" else NA_character_
  stat_label <- if (identical(stat_col, "F"))      "F"      else
                if (identical(stat_col, "Chi.sq")) "Chi.sq" else "Stat."
  sm_hdrs <- c("Smooth term","EDF","Ref.df", stat_label, "p-value","Significance")
  writeData(wb, ws4, as.data.frame(t(sm_hdrs)), startRow=cur4, startCol=1, colNames=FALSE)
  addStyle(wb, ws4, xl_hdr, rows=cur4, cols=1:6, gridExpand=TRUE)
  setRowHeights(wb, ws4, rows=cur4, heights=18); cur4 <- cur4+1
  for (j in seq_len(nrow(st))) {
    pv   <- st[j,"p-value"]
    sig  <- if (is.na(pv)) "—" else if (pv<0.001) "***" else if (pv<0.01) "**" else if (pv<0.05) "*" else "n.s."
    pv_s <- if (is.na(pv)) "—" else if (pv<0.001) "< 0.001" else as.character(round(pv,4))
    bg   <- if (j %% 2 == 0) XC_GREY else XC_WHITE
    sig_bg <- if (!is.na(pv) && pv < 0.05) XC_GREEN else bg
    stat_val <- if (!is.na(stat_col) && stat_col %in% names(st)) {
      v <- st[j, stat_col]
      if (is.numeric(v)) round(v, 2) else as.character(v)
    } else NA_real_
    rv <- data.frame(
      Term  = rownames(st)[j],
      EDF   = round(st[j,"edf"],   2),
      RefDF = round(st[j,"Ref.df"],2),
      F_val = stat_val,
      P     = pv_s,
      Sig   = sig
    )
    writeData(wb, ws4, rv, startRow=cur4, startCol=1, colNames=FALSE)
    addStyle(wb, ws4, xl_cell(bg,"left"),      rows=cur4, cols=1,   gridExpand=TRUE)
    addStyle(wb, ws4, xl_cell(bg,"center"),    rows=cur4, cols=2:5, gridExpand=TRUE)
    addStyle(wb, ws4, xl_cell(sig_bg,"center"),rows=cur4, cols=6,   gridExpand=TRUE)
    setRowHeights(wb, ws4, rows=cur4, heights=16); cur4 <- cur4+1
  }
  cur4 <- cur4+2  # gap between models
}
setColWidths(wb, ws4, cols=1:6, widths=c(36,10,10,10,12,14))

# ── Save Excel ──────────────────────────────────────────────────
xl_pfad <- file.path(output_ordner, "Block2_Uebersicht.xlsx")
saveWorkbook(wb, xl_pfad, overwrite=TRUE)
cat("✓ Excel overview saved:", basename(xl_pfad), "\n\n")

# ══════════════════════════════════════════════════════════════
# METHODENBERICHT  (Block2_Methodenbericht.docx, Deutsch)
# ══════════════════════════════════════════════════════════════

cat("Erstelle Methodenbericht (Word, Deutsch)...\n")

doc_m <- read_docx()

# Pipe-Helfer (officer nutzt Pipe-Syntax, aber wir arbeiten mit Funktionen)
h1m <- function(d, t) body_add_par(d, t, style = "heading 1")
h2m <- function(d, t) body_add_par(d, t, style = "heading 2")
h3m <- function(d, t) body_add_par(d, t, style = "heading 3")
pm  <- function(d, t) body_add_par(d, t, style = "Normal")
brm <- function(d)    body_add_break(d)

# ── Titelseite ─────────────────────────────────────────────────
doc_m <- h1m(doc_m, "Block 2 — Methodenbericht")
doc_m <- pm(doc_m,  "GAMM-Aktivitaetsanalyse | VHF Igelbesenderung Niedersachsen")
doc_m <- pm(doc_m,  "Projekt: Wildtierstation Sachsenhagen / TiHo Hannover")
doc_m <- pm(doc_m,  paste0("Erstellt: ", format(Sys.time(), "%d.%m.%Y, %H:%M"), " Uhr"))
doc_m <- pm(doc_m,  paste0("R-Version: ", R.version$major, ".", R.version$minor))
doc_m <- brm(doc_m)

# ── Vorbemerkung ───────────────────────────────────────────────
doc_m <- h1m(doc_m, "Vorbemerkung: Forschungsfrage und Skriptaufbau")
doc_m <- pm(doc_m, paste0(
  "Block 2 (Block2_GAMM.R) beantwortet die zentrale Frage: Veraendert sich die ",
  "Nachtaktivitaet von Igeln nach der Auswilderung systematisch ueber die Zeit? ",
  "Gibt es individuelle Unterschiede, die unabhaengig vom Zeittrend bestehen? ",
  "Das Skript verwendet Generalized Additive Mixed Models (GAMMs) via mgcv::bam(), ",
  "da diese Methode nichtlineare Zeittrends, individuelle Zufallseffekte und ",
  "Nicht-Normalverteilung der Response in einem einzigen Modell vereint. ",
  "Voraussetzung: Block0_Datenpipeline.R muss vorab ausgefuehrt worden sein und ",
  "gamm_tagesdaten.rds sowie gamm_nachtminuten.rds im Output-Ordner abgelegt haben."
))
doc_m <- brm(doc_m)

# ══════════════════════════════════════════════════════════════
doc_m <- h1m(doc_m, "1  Datengrundlage")
# ══════════════════════════════════════════════════════════════

doc_m <- h2m(doc_m, "1.1  Quelldateien")
doc_m <- pm(doc_m, paste0(
  "Block 2 laedt zwei RDS-Dateien aus Block 0: gamm_tagesdaten.rds enthaelt ",
  "fuer jedes Tier und jeden Tag einen aggregierten Datensatz mit dem Anteil ",
  "aktiver Nachtminuten (pct_aktiv), der Gesamtzahl beobachteter Minuten (n_min) ",
  "und der Anzahl Tage seit Auswilderung (tage_seit). gamm_nachtminuten.rds enthaelt ",
  "eine Zeile pro Minute der Nacht, also ob das Tier in dieser Minute aktiv (1) ",
  "oder passiv (0) war. Diese Minutendaten werden fuer Modell B (Binomial-GAMM) ",
  "und die Bout-Analyse benoetigt."
))
doc_m <- pm(doc_m, paste0(
  "Warum zwei Datenebenen: Die Tagesebene (Modell A) ist effizienter zu rechnen ",
  "und beantwortet direkt die Frage nach dem Zeittrend der Nachtaktivitaet. ",
  "Die Minutenebene (Modell B) erfordert deutlich mehr Rechenaufwand, erlaubt aber ",
  "zusaetzlich zu testen, ob sich der interne Rhythmus innerhalb der Nacht ",
  "(Aktivitaetszeitpunkt) ueber die Zeit nach Auswilderung verschiebt."
))

doc_m <- h2m(doc_m, "1.2  Einschlusskriterien")
doc_m <- pm(doc_m, paste0(
  "Nicht alle Tiere koennen in das GAMM eingeschlossen werden. Drei Kriterien ",
  "werden geprueft: (1) Mindestanzahl Naechte mit Nachtdaten (Standard: ",
  min_naechte, " Naechte). (2) Mindestanzahl absoluter Nachtminuten gesamt ",
  "(Standard: ", min_nacht_min, " Minuten). (3) Mittlere Nachtaktivitaet >= 1% ",
  "— Tiere darunter zeigen ein statisches Signal, das auf einen toten Sender oder ",
  "ein verungluecktes Tier hindeutet."
))
doc_m <- pm(doc_m, paste0(
  "WICHTIG: Es wird KEIN prozentualer Schwellenwert fuer den Anteil der Nachtdaten ",
  "an den Gesamtdaten verwendet. Der Grund: Im Hochsommer (Juni/Juli, ~52 Grad N) ",
  "dauert die Nacht biologisch nur ca. 8 Stunden — maximal ~33 Prozent aller ",
  "Tagesminuten fallen dann in die Nacht. Ein prozentualer Schwellenwert (z.B. 30 Prozent) ",
  "wuerde Sommer-Auswilderungen systematisch benachteiligen, obwohl die Datenqualitaet ",
  "dieselbe ist. Der absolute Schwellenwert (Anzahl Nachtminuten) ist tageslangen-neutral."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Die Einschlusskriterien sind konservativ gewaehlt, koennen aber trotzdem ",
  "Tiere ausschliessen, die biologisch interessant waeren (z.B. sehr kurz beobachtete ",
  "Tiere mit klarem Aktivitaetsmuster). Eine Anpassung der Schwellenwerte ist im ",
  "Einstellungs-Block am Skriptanfang moeglich und sollte je nach Datenlage geprueft werden."
))

doc_m <- h2m(doc_m, "1.3  Dynamischer Zeitcutoff (N-at-risk-Schwelle)")
doc_m <- pm(doc_m, paste0(
  "Ein zentrales methodisches Problem in Laengsschnittdaten mit variierender ",
  "Beobachtungsdauer ist der Survivorship Bias: Spaete Zeitpunkte werden nur noch ",
  "von wenigen Langzeitbeobachteten benoetigt — ihr individuelles Aktivitaetsmuster ",
  "bestimmt dann den 'Populationstrend', obwohl er eigentlich nur diese wenigen Tiere ",
  "beschreibt. Beispiel: Ab Tag 15 sind in dieser Stichprobe nur noch <= 10 Igel ",
  "im Datensatz, und ab Tag ~20 dominieren Igel22 und Igel23 die Schätzung."
))
doc_m <- pm(doc_m, paste0(
  "Loesung: Ein automatischer Cutoff-Mechanismus berechnet fuer jeden Tag, wie viele ",
  "eingeschlossene Tiere noch Daten liefern (N-at-risk). Der letzte Tag, an dem noch ",
  "mindestens n_min_igel = ", n_min_igel, " Tiere Daten haben, wird als Zeitcutoff ",
  "gesetzt. Alle Modelldaten werden auf dieses Fenster beschraenkt. Im aktuellen ",
  "Datensatz liegt der Cutoff bei Tag ", cutoff_tag, "."
))
doc_m <- pm(doc_m, paste0(
  "Wichtig: Die Visualisierungen zeigen aber einen weiter gehenden grauen Bereich ",
  "(ab Tag ", heatmap_vis_cutoff + 1, "), der als 'N <= 10' markiert ist. Dies soll ",
  "den Leser darueber informieren, dass die Kurven in diesem Bereich von sehr wenigen ",
  "Tieren gestuetzt werden und Populationsaussagen hier eingeschraenkt sind. Die ",
  "Modelle nutzen weiterhin alle verfuegbaren Daten bis zum Cutoff fuer eine ",
  "stabilere Schätzung der Zufallseffekte."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Der Cutoff-Wert n_min_igel ist eine willkuerliche Setzung. Bei sehr ",
  "kleinen Stichproben koennte er auf 3-4 Tiere reduziert werden. Die Entscheidung ",
  "sollte biologisch begruendet sein: Wie viele Individuen sind das Mindeste, um ",
  "eine Aussage ueber die Population zu treffen?"
))

doc_m <- h2m(doc_m, "1.4  Beobachtungsgewichte (Survivorship Bias Korrektur)")
doc_m <- pm(doc_m, paste0(
  "Auch innerhalb des Analysefensters sind nicht alle Tiere gleich stark repraesentiert: ",
  "Ein Tier mit 60 Beobachtungstagen liefert 60 Datenpunkte, eines mit 10 Tagen nur 10. ",
  "Ohne Korrektur dominieren die langzeitbeobachteten Tiere die spaeten Zeitpunkte und ",
  "beeinflussen dadurch den gesaetzten Trend disproportional."
))
doc_m <- pm(doc_m, paste0(
  "Loesung: Jedem Datenpunkt eines Tieres wird ein Gewicht zugewiesen nach der Formel ",
  "Gewicht = 1 / sqrt(Beobachtungsdauer_in_Tagen), normiert auf Mittelwert = 1. ",
  "Die Quadratwurzel daempft extreme Gewichte (ein Tier mit 5 Tagen Beobachtung wuerde ",
  "sonst ein 12-fach hoeheres Gewicht bekommen als eines mit 60 Tagen — mit Wurzel nur ",
  "das 3.5-fache). Sensitivitaet: Modell A wird mit und ohne Gewichte gefittet. Wenn ",
  "beide Versionen dasselbe Signifikanzergebnis liefern, hat der Bias keinen ",
  "qualitativ entscheidenden Einfluss."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Die Gewichtung ist eine Heuristik, kein theoretisch exaktes Verfahren. ",
  "Alternativ koennten mixed models mit einem zeitabhaengigen Sampling-Term modelliert ",
  "werden (komplexer). Bei gleichmaessig verteilten Beobachtungsdauern ist die Korrektur ",
  "nicht notwendig."
))
doc_m <- brm(doc_m)

# ══════════════════════════════════════════════════════════════
doc_m <- h1m(doc_m, "2  Explorative Analysen")
# ══════════════════════════════════════════════════════════════

doc_m <- h2m(doc_m, "2.1  Rohdaten-Spaghetti-Plot und Facetten-Plot")
doc_m <- pm(doc_m, paste0(
  "Vor jedem Modell werden die Rohdaten visualisiert, um direkt abzulesen, ob ",
  "ein Trend sichtbar ist. Der Spaghetti-Plot zeigt jedes Tier als einzelne Linie ",
  "plus einem LOESS-Populationstrend (schwarz). Der Facetten-Plot zeigt jedes Tier ",
  "in einem eigenen Panel mit individuellem LOESS-Trend (orange). Beide Plots ",
  "verwenden als x-Achse die Tage seit Auswilderung und als y-Achse den Anteil ",
  "aktiver Nachtminuten in Prozent (pct_aktiv_pct = pct_aktiv * 100)."
))
doc_m <- pm(doc_m, paste0(
  "Begruendung: Explorative Visualisierungen vor der Modellierung sind essenziell. ",
  "Sie zeigen sofort, ob ein visuell offensichtlicher Trend vorhanden ist, ob ",
  "einzelne Tiere Ausreisser sind, und ob die Modellannahmen plausibel erscheinen. ",
  "Eine starke visuelle Evidenz ergibt und unterstuetzt statistisch signifikante Befunde."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: LOESS-Trends glaetten die Daten und koennen bei wenigen Tieren ",
  "sehr breit werden oder stark durch Ausreisser beeinflusst sein. Die Konfidenzbaender ",
  "des LOESS sind nur explorative Hilfslinien, keine statistischen Inferenzaussagen."
))

doc_m <- h2m(doc_m, "2.2  Frueh- vs. Spaetphasen-Vergleich (Wilcoxon)")
doc_m <- pm(doc_m, paste0(
  "Als erster statistischer Test wird geprueft, ob sich die mittlere Nachtaktivitaet ",
  "in der fruehen Phase (Tage 1–", grenze, ") von der spaeten Phase (Tage >", grenze,
  ") unterscheidet. Hierfuer wird der Wilcoxon-Rangsummentest eingesetzt ",
  "(nicht-parametrisch, da Normalverteilung der Nachtaktivitaetswerte nicht vorausgesetzt)."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Dieser einfache Test ignoriert die hierarchische Struktur der Daten ",
  "(mehrere Beobachtungen pro Tier), behandelt alle Beobachtungspunkte als unabhaengig, ",
  "und hat wenig Power bei kleinen Stichproben. Er ist als grober Indikator gedacht, ",
  "nicht als Ersatz fuer das GAMM."
))

doc_m <- h2m(doc_m, "2.3  Erste-Tage-Analyse (Tage 1–5, Nocturnal Index)")
doc_m <- pm(doc_m, paste0(
  "Die Tage 1 bis 5 nach der Auswilderung werden isoliert analysiert, weil: ",
  "(1) In diesem Zeitfenster sind fast alle Tiere noch im Datensatz — Survivorship Bias ",
  "ist minimal. (2) Wenn es eine Eingewoehnungszeit gibt, wuerde sie sich in den ersten ",
  "Tagen als veraenderter Nocturnal Index zeigen. ",
  "Als Mass wird der Nocturnal Index (NI = Nachtaktivitaet% minus Tagaktivitaet%) ",
  "verwendet. NI > 0 bedeutet: Das Tier ist nachtaktiver als tagaktiv. NI = 0 ",
  "bedeutet gleichmaessige Aktivitaet ueber den Tag. NI < 0 bedeutet ueberraschenderweise ",
  "mehr Tagaktivitaet. Statistisch: Kruskal-Wallis-Test ueber Tage 1–5, ",
  "sowie paarweiser Wilcoxon-Test Tag 1 vs. Tag 5."
))
doc_m <- pm(doc_m, paste0(
  "Begruendung fuer Nocturnal Index statt absolutem Nachtanteil: Der NI kontrolliert ",
  "fuer individuelle Unterschiede im Gesamtaktivitaetsniveau und ist unabhaengig von ",
  "der Tageslange — sinnvoll, da Herbst-Tiere und Fruehlings-Tiere sehr verschiedene ",
  "Nacht/Tag-Verhaeltnisse haben."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Bei sehr kleinen N pro Tag (z.B. nur 5–8 Tiere) haben alle Tests ",
  "geringe statistische Power. Nicht-signifikante Ergebnisse koennen also auf ",
  "fehlende Power zurueckzufuehren sein, nicht unbedingt auf fehlende Effekte."
))
doc_m <- brm(doc_m)

# ══════════════════════════════════════════════════════════════
doc_m <- h1m(doc_m, "3  Nocturnal Index und Fruehphasen-GAMM")
# ══════════════════════════════════════════════════════════════

doc_m <- h2m(doc_m, "3.1  Nocturnal Index (alle Tage)")
doc_m <- pm(doc_m, paste0(
  "Fuer alle Tiere und alle Tage wird der Nocturnal Index aus dem breiten Format ",
  "berechnet (Tag und Nacht als separate Spalten per dcast). NI = Nachtanteil minus ",
  "Taganteil, in Prozent. Dargestellt als Spaghetti-Plot (lineare Zeitachse) und ",
  "als Plot mit log-transformierter Zeitachse. Die logarithmische Zeitachse (log10) ",
  "streckt die fruehen Tage auseinander und verkleinert spaete Zeitpunkte — ",
  "ideal um zu sehen, ob sich in der ersten Woche etwas veraendert."
))

doc_m <- h2m(doc_m, "3.2  Fruehphasen-GAMM (Tage 1–21, Beta-Regression)")
doc_m <- pm(doc_m, paste0(
  "Ein gesondertes GAMM wird nur auf die ersten ", frueh_tage, " Tage gefittet, ",
  "um gezielt die Anpassungsphase zu analysieren. Als Zeitpraediktor wird ",
  "log(tage_seit + 1) verwendet statt tage_seit — diese Transformation streckt ",
  "die fruehen Tage auseinander und erlaubt dem Smooth, feine Veraenderungen in ",
  "der ersten Woche zu erkennen, die auf linearer Skala 'verschwoemmern' wuerden."
))
doc_m <- pm(doc_m, paste0(
  "Modellformel: ni_transf ~ s(log_tage, bs='tp', k=6) + ",
  "s(log_tage, igel, bs='fs', k=4, m=1) + s(igel, bs='re'). ",
  "Response: der auf (0,1) transformierte und mit Smithson-Verkuilen-Methode ",
  "gestauchte Nocturnal Index (noetig fuer Beta-Regression, die keine exakten 0 oder 1 erlaubt). ",
  "Die Familie betar(link='logit') passt zu diesem stetigen Anteilsdatum."
))
doc_m <- pm(doc_m, paste0(
  "Begruendung: Ein GAMM auf nur 21 Tage mit logarithmischer Zeitskala ist die ",
  "statistisch saerste Methode, um kurzfristige, moeglicherweis nichtlineare Anpassung ",
  "zu detektieren. Der Faktor-Smooth (fs) erlaubt dabei individuelle Kurvenformen ",
  "mit Shrinkage Richtung Populationsmittel — Tiere mit wenig Daten werden nicht ",
  "ueberparametrisiert."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Die Beta-Regression erfordert Werte strikt zwischen 0 und 1. Die ",
  "Smithson-Verkuilen-Transformation (y * (n-1) + 0.5) / n ist eine anerkannte Methode, ",
  "ist aber eine leichte Verzerrung der Rohdaten. Das Modell wird nur gefittet, wenn ",
  "mindestens 20 Datenpunkte und mindestens 3 Tiere in der Fruehphase vorhanden sind."
))
doc_m <- brm(doc_m)

# ══════════════════════════════════════════════════════════════
doc_m <- h1m(doc_m, "4  Modell A — Nachtaktivitaet ueber Zeit (Tagesebene)")
# ══════════════════════════════════════════════════════════════

doc_m <- h2m(doc_m, "4.1  Modellspezifikation")
doc_m <- pm(doc_m, paste0(
  "Modell A ist ein Beta-GAMM auf Tagesebene. Die Response ist der transformierte ",
  "Anteil aktiver Nachtminuten (pct_transf, im Intervall (0,1)). Praeditoren: ",
  "s(tage_seit, bs='tp', k=8) — nichtlinearer Populationszeittrend (Thin-Plate-Spline). ",
  "s(tage_seit, igel, bs='fs', k=5, m=1) — individuelle Abweichungen vom ",
  "Populationstrend mit Shrinkage (Tiere mit wenig Daten werden zum Mittel gezogen). ",
  "s(igel, bs='re') — zufaelliges Intercept je Tier (erklaert konstante individuelle ",
  "Niveauunterschiede in der Nachtaktivitaet)."
))
doc_m <- pm(doc_m, paste0(
  "Begruendung fuer Beta-Regression: Der Anteil aktiver Nachtminuten ist ein ",
  "stetiger Anteilswert zwischen 0 und 1. Eine Normalverteilung waere unangemessen, ",
  "da Werte nahe 0 oder 1 durch Bodensatz-/Deckeneffekte vorkommen koennen und die ",
  "Varianz bei mittleren Anteilen groesser ist als bei extremen. Die Beta-Verteilung ",
  "modelliert genau dieses Muster."
))
doc_m <- pm(doc_m, paste0(
  "Begruendung fuer bam() statt gam(): bam() (big data GAM) ist eine numerisch ",
  "effizientere Implementierung, die bei grossen Datensaetzen erheblich schneller ist. ",
  "Die fREML-Schaetzung (fast Restricted Maximum Likelihood) ist fuer Mixed Models ",
  "mit vielen Levels des Zufallseffekts (viele Igel) besonders geeignet."
))
doc_m <- pm(doc_m, paste0(
  "Smithson-Verkuilen-Transformation: Da Beta-Regression keine Werte von exakt 0 oder 1 ",
  "erlaubt, wird angewendet: pct_transf = (pct_aktiv * (n_total - 1) + 0.5) / n_total, ",
  "wobei n_total die Anzahl der Nachtminuten je Datenpunkt ist. Diese Transformation ",
  "ist abhaengig von n_total — lange beobachtete Naechte werden weniger stark verzerrt."
))

doc_m <- h2m(doc_m, "4.2  Sensitivitaetsanalyse: mit vs. ohne Gewichte")
doc_m <- pm(doc_m, paste0(
  "Modell A wird parallel ohne Beobachtungsgewichte gefittet. Wenn beide Modelle ",
  "(mit und ohne Gewichte) dasselbe Signifikanzergebnis fuer s(tage_seit) liefern, ",
  "hat der Survivorship Bias keinen qualitativen Einfluss auf die Schlussfolgerung. ",
  "Wenn die Ergebnisse divergieren, ist Vorsicht geboten und der Bias ist substanziell."
))

doc_m <- h2m(doc_m, "4.3  Diagnostik")
doc_m <- pm(doc_m, paste0(
  "Nach dem Modell wird gam.check() ausgefuehrt, das vier Diagnostik-Plots erstellt: ",
  "QQ-Plot der Residuen (Normalverteilung der Deviance-Residuen?), Residuen vs. ",
  "Fitted Values, Histogramm der Residuen, Response vs. Fitted. Ziel: Keine ",
  "starken systematischen Muster, keine extremen Ausreisser. Concurvity (Analogon ",
  "zur Multikollinearitaet bei Smooths) wird ebenfalls geprueft: Werte < 0.8 gelten ",
  "als unproblematisch."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Bei Beta-GAMMs sind die Diagnostik-Plots auf Deviance-Residuen basiert, ",
  "die nicht normalverteilt sind — der QQ-Plot ist daher kein exakter Test der ",
  "Verteilungsannahme. Die visuelle Pruefung ist dennoch sinnvoll. Autokorrelation ",
  "der Residuen innerhalb eines Tieres ueber die Zeit (zeitliche Abhaengigkeit) ist ",
  "nicht explizit modelliert, koennte aber die Standardfehler verzerren. Dies waere ",
  "mit einem AR1-Korrelationsstruktur-Term behebbar (komplexer, weniger stabil)."
))
doc_m <- brm(doc_m)

# ══════════════════════════════════════════════════════════════
doc_m <- h1m(doc_m, "5  Modell B — 24h-Rhythmus und zeitliche Verschiebung")
# ══════════════════════════════════════════════════════════════

doc_m <- h2m(doc_m, "5.1  Modellspezifikation")
doc_m <- pm(doc_m, paste0(
  "Modell B ist ein Binomial-GAMM auf Minutenebene. Jede Minute der Nacht ist ",
  "ein Datenpunkt: aktiv (1) oder passiv (0). Praeditoren: ",
  "s(tage_seit, bs='tp', k=8) — Zeittrend. ",
  "s(stunde, bs='cr', k=10) — 24h-Rhythmus innerhalb der Nacht (kubischer Spline). ",
  "ti(tage_seit, stunde, bs=c('tp','cr'), k=c(5,6)) — Tensor-Interaktionsterm: ",
  "veraendert sich der Aktivitaetszeitpunkt innerhalb der Nacht ueber die Zeit? ",
  "s(igel, bs='re') — Zufaelliges Intercept je Tier."
))
doc_m <- pm(doc_m, paste0(
  "Begruendung fuer bs='cr' statt 'cc' (cyclic): Zirkulaere Splines (cc) ",
  "erzwingen, dass der Wert am Anfang (0 Uhr) und Ende (24 Uhr) identisch sind. ",
  "Da aber nur Nachtminuten analysiert werden (nicht der volle 24h-Zyklus), ",
  "decken die Daten nicht den vollen Bereich ab — eine zirkulaere Randbedingung ",
  "waere falsch. Der kubische Spline (cr) ist hier geeigneter."
))
doc_m <- pm(doc_m, paste0(
  "Begruendung fuer Tensor-Interaktion (ti): Ein einfacher Haupteffekt-Smooth ",
  "s(stunde) wuerde nur testen, ob es einen Stunden-Rhythmus gibt, nicht ob dieser ",
  "sich veraendert. Der ti()-Term erlaubt, die Interaktion zwischen Zeitpunkt-in-der-Nacht ",
  "und Tage-nach-Auswilderung zu modellieren — genau die Forschungsfrage: ",
  "Verschiebt sich das Aktivitaetsfenster innerhalb der Nacht?"
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Das Binomial-GAMM auf Minutenebene hat einen sehr grossen Datensatz ",
  "(viele Tausend Minuten). Sehr kleine Effekte koennen dadurch statistisch signifikant ",
  "werden, obwohl sie biologisch nicht relevant sind. Die praktische Bedeutsamkeit ",
  "eines signifikanten ti()-Terms muss daher auch visuell (Heatmap) beurteilt werden. ",
  "Ausserdem: Aufeinanderfolgende Minuten desselben Tieres sind nicht unabhaengig ",
  "(Autokorrelation), was die Standardfehler in der gleichen Richtung verzerren koennte ",
  "wie bei Modell A."
))
doc_m <- brm(doc_m)

# ══════════════════════════════════════════════════════════════
doc_m <- h1m(doc_m, "6  Bout-Analyse")
# ══════════════════════════════════════════════════════════════

doc_m <- h2m(doc_m, "6.1  Was ist ein Aktivitaets-Bout?")
doc_m <- pm(doc_m, paste0(
  "Ein 'Bout' ist eine zusammenhaengende Phase mit Aktivitaet, unterbrochen von ",
  "einer Pause laenger als das Bout-Kriterium. Das Bout-Kriterium (minimum inter-bout ",
  "interval, MIBI) ist die minimale Pause zwischen zwei Aktivitaetsphasen, die sie ",
  "als separate Bouts definiert. Beispiel: Kriterium 10 Minuten bedeutet, dass eine ",
  "Pause < 10 Minuten als kurze Unterbrechung innerhalb desselben Bouts gilt."
))

doc_m <- h2m(doc_m, "6.2  Bestimmung des Bout-Kriteriums: Log-Survivor-Plot")
doc_m <- pm(doc_m, paste0(
  "Das Bout-Kriterium wird empirisch aus den Daten bestimmt. Dafuer werden alle ",
  "Inter-Bout-Intervalle (IBI = Pausen zwischen aktiven Minuten) berechnet und im ",
  "Log-Survivor-Plot (log des Anteils der IBIs laenger als x) dargestellt. ",
  "In diesem Plot erscheinen Verteilungen mit zwei Komponenten als Kurve mit einem ",
  "Knickpunkt: Kurze IBIs (innerhalb eines Bouts) bilden einen steilen Abfall, ",
  "laengere IBIs (zwischen Bouts) einen flacheren. Der Knickpunkt definiert das ",
  "Bout-Kriterium (typisch 5–20 Minuten bei Kleinsaeuger-Aktivitaetsdaten). ",
  "WICHTIG: Das Kriterium muss nach Sichtung des Plots manuell im Skript eingetragen ",
  "werden (Variable bout_kriterium, aktuell: ", bout_kriterium, " Minuten)."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Der Knickpunkt im Log-Survivor-Plot ist oft nicht eindeutig sichtbar, ",
  "besonders bei kleinen Stichproben oder wenn die Aktivitaet sehr gleichmaessig ist. ",
  "Verschiedene Bout-Kriterien fuehren zu verschiedenen Metriken — die Wahl ist eine ",
  "explizite biologische Entscheidung. Eine Sensitivitaetsanalyse mit verschiedenen ",
  "Kriterien wird empfohlen."
))

doc_m <- h2m(doc_m, "6.3  Bout-Metriken pro Nacht")
doc_m <- pm(doc_m, paste0(
  "Aus dem gewaehlen Bout-Kriterium werden pro Nacht und Tier berechnet: ",
  "n_bouts (Anzahl Aktivitaetsbouts), mean_dauer (mittlere Boutdauer in Minuten), ",
  "max_dauer (laengster Bout der Nacht), mean_ibi (mittleres Inter-Bout-Intervall). ",
  "Kurze Pausen innerhalb des Bout-Kriteriums werden mit benachbarten Bouts ",
  "zusammengefuehrt (merge-loop, max. 3 Iterationen)."
))

doc_m <- h2m(doc_m, "6.4  GAMMs auf Bout-Metriken")
doc_m <- pm(doc_m, paste0(
  "Zwei separate GAMMs testen, ob sich Boutanzahl und Boutdauer ueber die Zeit veraendern. ",
  "Boutanzahl: Poisson-GAMM (log-Link), da es sich um Zaehlwerte handelt. ",
  "Boutdauer: Gamma-GAMM (log-Link), da Dauern positiv und rechtschief sind. ",
  "Praeditoren: s(tage_seit) + s(tage_seit, igel, fs) + s(igel, re). ",
  "Ein signifikanter Zeittrend in der Boutanzahl wuerde zeigen, dass die Tiere ",
  "ihre Aktivitaet in mehr oder weniger Phasen pro Nacht aufteilen."
))
doc_m <- pm(doc_m, paste0(
  "Limitation: Bei kleinen Stichproben (wenige Naechte pro Tier) koennen Bout-Metriken ",
  "sehr variabel sein. Einzelne Naechte mit sehr vielen oder sehr wenigen Bouts ",
  "koennen die Ergebnisse stark beeinflussen. Neben dem Modell sollten die Rohdaten ",
  "immer visuell geprueft werden."
))
doc_m <- brm(doc_m)

# ══════════════════════════════════════════════════════════════
doc_m <- h1m(doc_m, "7  Modellvergleich und Ergebnisinterpretation")
# ══════════════════════════════════════════════════════════════

doc_m <- h2m(doc_m, "7.1  Modellguete")
doc_m <- pm(doc_m, paste0(
  "Zum Vergleich der Modelle werden erklärte Devianz (Deviance explained, in Prozent) ",
  "und AIC (Akaike Information Criterion) herangezogen. Die erklaerte Devianz ist das ",
  "GAMM-Aequivalent zum R² in linearen Modellen: 0 Prozent = Modell erklaert nichts, ",
  "100 Prozent = perfekte Anpassung. AIC bewertet Modellfit unter Beru cksichtigung ",
  "der Modellkomplexitaet (Anzahl Freiheitsgrade)."
))

doc_m <- h2m(doc_m, "7.2  Biologische Interpretation der Smooth Terms")
doc_m <- pm(doc_m, paste0(
  "Der EDF (Estimated Degrees of Freedom) eines Smooth-Terms gibt Auskunft ueber ",
  "seine Nichtlinearitaet: EDF ~ 1 entspricht einer linearen Beziehung, EDF > 1 ",
  "ist nichtlinear, EDF >> 1 ist komplex-kurvenfoermig. Ein signifikanter s(tage_seit) ",
  "in Modell A bedeutet, dass sich die Nachtaktivitaet systematisch veraendert. ",
  "Nicht-Signifikanz bedeutet stabile Nachtaktivitaet — biologisch interpretierbar ",
  "als sofortige Habituation ohne Anpassungszeit."
))
doc_m <- pm(doc_m, paste0(
  "Wichtige Einschraenkung: Alle p-Werte in GAMMs sind Naeherungsloesungen. Die ",
  "Freiheitsgrade der Smooth-Terme werden geschaetzt, nicht exakt berechnet. ",
  "Bei kleinen Stichproben koennen sie konservativ (zu wenig Power) oder liberal ",
  "(zu leicht signifikant) sein. Die visuelle Ueberpruefung der Smooths ergaenzt ",
  "immer die formale Signifikanzaussage."
))

doc_m <- h2m(doc_m, "7.3  Signifikanzniveau und Software")
doc_m <- pm(doc_m, paste0(
  "Alle Tests werden mit zweiseitigem Signifikanzniveau alpha = 0.05 durchgefuehrt. ",
  "Alle Modelle wurden in R (Version ", R.version$major, ".", R.version$minor, ") ",
  "mit den Paketen mgcv (GAMM-Schätzung), data.table (Datentransformation), ",
  "ggplot2 (Visualisierung), officer/flextable (Word-Bericht) und openxlsx ",
  "(Excel-Uebersicht) erstellt."
))
doc_m <- brm(doc_m)

# ══════════════════════════════════════════════════════════════
doc_m <- h1m(doc_m, "8  Output-Dateien")
# ══════════════════════════════════════════════════════════════
doc_m <- pm(doc_m, paste0(
  "Alle Output-Dateien werden im Ordner output/Block2_GAMM/ gespeichert. ",
  "Grafiken (PNG): abdeckung_gantt, abdeckung_n_risk (Beobachtungsabdeckung), ",
  "explorativ_spaghetti, explorativ_facetten (Rohdaten), ",
  "erste_tage_boxplot/spaghetti/heatmap (Tage 1–5), ",
  "ni_spaghetti/log/facetten (Nocturnal Index), gamm_fruehphase (Fruehphasen-GAMM), ",
  "gamm_A_nachtaktivitaet (Modell A: Populationstrend + individuelle Kurven + N-at-risk), ",
  "gamm_A_diagnostik (gam.check), gamm_B_rhythmus (Heatmap + Zeittrend Modell B), ",
  "gamm_B_diagnostik, bout_logsurvivor (Bout-Kriterium), gamm_bout_metriken. ",
  "Modelle als RDS: gamm_modelle.rds (alle Modellobjekte). ",
  "Word: Block2_GAMM_Ergebnisbericht.docx (Ergebnisse + Plots), ",
  "Block2_Methodenbericht.docx (dieser Methodenbericht). ",
  "Excel: Block2_Uebersicht.xlsx (Stichprobe, Modellguete, deskriptive Statistiken)."
))

# ── Speichern ──────────────────────────────────────────────────
meth_pfad <- file.path(output_ordner, "Block2_Methodenbericht.docx")
print(doc_m, target = meth_pfad)
cat("✓ Methodenbericht gespeichert:", basename(meth_pfad), "\n")

cat("\n── Output-Dateien Block 2 ──────────────────────────────\n")
cat("  Plots:  ", length(list.files(output_ordner, "\\.png$")), "PNG-Dateien\n")
cat("  Excel:  Block2_Uebersicht.xlsx\n")
cat("  Word:   Block2_GAMM_Ergebnisbericht.docx  (Ergebnisse)\n")
cat("  Word:   Block2_Methodenbericht.docx        (Methodik)\n")
cat("  RDS:    gamm_modelle.rds\n")
cat("\n✓ Block 2 komplett!\n")
