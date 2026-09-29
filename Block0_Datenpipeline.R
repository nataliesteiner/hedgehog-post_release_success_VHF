# ==============================================================
# VHF activity analysis — ALL HEDGEHOGS (batch)
# ==============================================================
# Project:  Hedgehog VHF telemetry, Lower Saxony
# Author:   Natalie Steiner
# Method:   pred_nested_loio_smoothed_wmv, minute reduction,
#           modal value, day/night via NOAA algorithm (base R)
# ==============================================================
# IMPORTANT: always run the script from line 1 (Ctrl+Alt+R in
#            RStudio, or source("Block0_Datenpipeline.R")).
#            Only then are all packages loaded correctly.
# ==============================================================

# Fehlende Pakete automatisch installieren
pakete <- c("data.table", "lubridate", "ggplot2",
            "readxl", "scales", "gridExtra", "openxlsx", "patchwork")
neu    <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) {
  cat("Installiere fehlende Pakete:", paste(neu, collapse = ", "), "\n")
  install.packages(neu)
}

library(data.table)
library(lubridate)
library(ggplot2)
library(readxl)
library(scales)
library(gridExtra)
library(openxlsx)
library(patchwork)

# ──────────────────────────────────────────────────────────────
# SONNENAUF/-UNTERGANG  (ohne suncalc — NOAA-Algorithmus in Base R)
# Genauigkeit: ±1 Minute für mitteleuropäische Breiten, ausreichend
# für Tag/Nacht-Klassifikation.
# ──────────────────────────────────────────────────────────────
.sonnenzeiten <- function(datum, lat, lon, tz = "UTC") {
  # datum: Date-Vektor, lat/lon: Skalare
  d <- as.Date(datum)
  # Julianisches Datum (Mittag UTC)
  jd <- as.numeric(d) + 2440587.5 + 0.5

  n     <- jd - 2451545.0
  L     <- (280.460 + 0.9856474 * n) %% 360
  g_deg <- (357.528 + 0.9856003 * n) %% 360
  g_rad <- g_deg * pi / 180
  lam   <- L + 1.915 * sin(g_rad) + 0.020 * sin(2 * g_rad)

  eps_deg <- 23.439 - 4e-7 * n
  eps_rad <- eps_deg * pi / 180
  lam_rad <- lam * pi / 180

  # Deklination (Grad)
  delta_rad <- asin(sin(eps_rad) * sin(lam_rad))

  # Zeitgleichung (Minuten)  — vereinfachte Version
  rl <- L * pi / 180
  rl2 <- lam_rad
  E_min <- (-1.915 * sin(g_rad) - 0.020 * sin(2 * g_rad)
            + 2.466 * sin(2 * rl2) - 0.053 * sin(4 * rl2))

  # Stundenwinkel Sonnenaufgang (cos H0)
  lat_rad <- lat * pi / 180
  # Atmosphärische Refraktion + Sonnenhalbdurchmesser: -0.8333°
  cos_h0 <- (sin(-0.8333 * pi / 180) - sin(lat_rad) * sin(delta_rad)) /
              (cos(lat_rad) * cos(delta_rad))

  # Polar-Sonderfälle: Polarnacht / Mitternachtssonne
  polar <- cos_h0 < -1 | cos_h0 > 1
  h0_deg <- ifelse(polar, NA_real_, acos(pmax(-1, pmin(1, cos_h0))) * 180 / pi)

  # Solarer Mittag (Stunden UTC)
  noon_utc <- 12 - lon / 15 - E_min / 60

  sr_utc <- noon_utc - h0_deg / 15   # Stunden UTC
  ss_utc <- noon_utc + h0_deg / 15

  # In POSIXct (UTC) umwandeln
  origin_utc <- as.POSIXct(paste0(d, " 00:00:00"), tz = "UTC")
  sunrise_utc <- origin_utc + sr_utc * 3600
  sunset_utc  <- origin_utc + ss_utc * 3600

  # Lokale Zeitzone
  sunrise_loc <- with_tz(sunrise_utc, tzone = tz)
  sunset_loc  <- with_tz(sunset_utc,  tzone = tz)

  data.table(date = d, sunrise = sunrise_loc, sunset = sunset_loc)
}

cat("✓ Alle Pakete geladen\n\n")

# ──────────────────────────────────────────────────────────────
# EINSTELLUNGEN  — nur hier anpassen
# ──────────────────────────────────────────────────────────────

# Projektwurzel — einzige Zeile die du ggf. anpassen musst
projekt_root <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"

# Ordner mit den CSV-Aktivitätsdateien
daten_ordner <- file.path(projekt_root, "data", "activity")

# Metadaten-Datei (XLSX)
meta_datei <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")

# Ausgabeordner für Plots und Tabellen (wird automatisch erstellt)
output_ordner <- file.path(projekt_root, "output", "Block0_Pipeline")

# Standort (Niedersachsen)
standort_lat <- 52.39729710523643
standort_lon <-  9.216876248766871
zeitzone     <- "Europe/Berlin"

# Klassifikationsspalte
klasse_spalte <- "pred_nested_loio_smoothed_wmv"

# Schwellenwert für "toter Sender"-Erkennung:
# Tage am Ende des Datensatzes mit weniger als X% Aktivität
# werden als statisches Signal gewertet und entfernt.
# Empfehlung: 5% (konservativ). Erhöhen = aggressiver trimmen.
schwelle_aktiv_pct <- 5

# ──────────────────────────────────────────────────────────────
# HILFSFUNKTIONEN
# ──────────────────────────────────────────────────────────────

# Modalwert (häufigster Wert in einem Vektor)
modal_val <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1L]
}

# Hauptaktivitätsphase: längster zusammenhängender Block aktiver Stunden
# (aktiv = Stunden über dem Durchschnitt der Gesamtaktivität)
# Sucht auch über Mitternacht hinweg (zirkulär)
find_hauptphase <- function(profil_dt) {
  # profil_dt: data.table mit Spalten stunde_rund (0-23) und pct_aktiv
  alle <- data.table(stunde_rund = 0:23)
  p <- merge(alle, profil_dt, by = "stunde_rund", all.x = TRUE)
  p[is.na(pct_aktiv), pct_aktiv := 0]
  setorder(p, stunde_rund)

  schwelle <- mean(p$pct_aktiv)          # Schwellenwert = Gesamtdurchschnitt
  p[, ist_aktiv := pct_aktiv >= schwelle]

  # Stunden verdoppeln für zirkuläre Suche (Aktivität über Mitternacht)
  aktiv_vec <- c(p$ist_aktiv, p$ist_aktiv)

  bester_start <- NA_integer_
  bester_len   <- 0L
  aktueller_start <- NA_integer_
  aktueller_len   <- 0L

  for (j in seq_along(aktiv_vec)) {
    if (aktiv_vec[j]) {
      if (is.na(aktueller_start)) { aktueller_start <- j; aktueller_len <- 1L }
      else aktueller_len <- aktueller_len + 1L
      # Nur Blöcke bis max. 24h (Verdoppelung verhindert Endlosschleifen)
      if (aktueller_len <= 24L && aktueller_len > bester_len) {
        bester_len   <- aktueller_len
        bester_start <- aktueller_start
      }
    } else {
      aktueller_start <- NA_integer_
      aktueller_len   <- 0L
    }
  }

  if (is.na(bester_start)) {
    return(list(start_h = NA_real_, ende_h = NA_real_, dauer_h = 0L))
  }

  start_h <- p$stunde_rund[(bester_start - 1L) %% 24L + 1L]
  ende_h  <- p$stunde_rund[(bester_start + bester_len - 2L) %% 24L + 1L]
  list(start_h = start_h, ende_h = ende_h, dauer_h = bester_len)
}

# NULL / NA Fallback
`%||%` <- function(a, b) if (length(a) > 0 && !is.na(a[1])) a[1] else b

# Farbpalette Tag/Nacht
farben <- c("Tag" = "#E8A87C", "Nacht" = "#2C3E6B")

# ──────────────────────────────────────────────────────────────
# AUSGABEORDNER ERSTELLEN
# ──────────────────────────────────────────────────────────────

if (!dir.exists(output_ordner)) {
  dir.create(output_ordner, recursive = TRUE)
  cat("Ausgabeordner erstellt:", output_ordner, "\n")
}

# ──────────────────────────────────────────────────────────────
# METADATEN LADEN  (CSV oder XLSX automatisch erkannt)
# ──────────────────────────────────────────────────────────────

# Pfad expandieren (wichtig: Leerzeichen im Pfad vertragen fread() sonst nicht)
meta_pfad <- path.expand(meta_datei)
cat("Lade Metadaten aus:", meta_pfad, "\n")

if (grepl("\\.csv$", meta_datei, ignore.case = TRUE)) {
  # CSV einlesen — file= explizit angeben, damit Leerzeichen im Pfad kein Problem sind
  meta <- fread(file = meta_pfad, encoding = "UTF-8")
} else {
  meta <- as.data.table(read_excel(meta_pfad, sheet = 1))
}

# Spaltennamen vereinheitlichen (egal ob CSV oder Excel)
setnames(meta,
  old = c("individual", "sex", "diagnosis_main", "tagging_weight",
          "soft_release_start", "date_release", "date_last_signal",
          "tagging_period"),
  new = c("igel", "geschlecht", "diagnose", "gewicht_g",
          "soft_release", "hard_release", "letztes_signal",
          "besenderungsdauer"),
  skip_absent = TRUE)

# Igel-Namen bereinigen (Leerzeichen entfernen)
meta[, igel := trimws(igel)]

# Datumsspalten parsen — erkennt beide Formate: DD.MM.YYYY und YYYY-MM-DD
parse_datum <- function(x) {
  x <- as.character(x)
  result <- suppressWarnings(as.Date(x, format = "%d.%m.%Y"))   # deutsches Format
  na_idx <- is.na(result)
  result[na_idx] <- suppressWarnings(as.Date(x[na_idx]))        # ISO-Format als Fallback
  result
}

meta[, hard_release   := parse_datum(hard_release)]
meta[, soft_release   := parse_datum(soft_release)]
meta[, letztes_signal := parse_datum(letztes_signal)]

cat("Metadaten geladen:", nrow(meta), "Igel\n\n")
print(meta[, .(igel, geschlecht, gewicht_g, hard_release,
               letztes_signal, besenderungsdauer)])

# ──────────────────────────────────────────────────────────────
# CSV-DATEIEN FINDEN
# ──────────────────────────────────────────────────────────────

alle_csvs <- list.files(daten_ordner,
                        pattern = "classification_active_passive.*\\.csv$",
                        full.names = TRUE)

cat("\nGefundene CSV-Dateien:", length(alle_csvs), "\n")

# Für jeden Igel die passende CSV suchen (flexibel — erkennt auch
# Sondernamen wie "Igel3_150.080" oder "Igel 5")
finde_csv <- function(igel_name) {
  nr <- gsub("[^0-9]", "", igel_name)          # nur Nummer extrahieren
  # Muster: "Igel" + optional Leerzeichen/Unterstrich + Nummer + Nicht-Zahl
  treffer <- grep(paste0("Igel[[:space:]_]?0*", nr, "[^0-9]"),
                  alle_csvs, value = TRUE)
  if (length(treffer) == 0) return(NA_character_)
  treffer[1]
}

meta[, csv_pfad := sapply(igel, finde_csv)]

# Überblick
cat("\n── CSV-Zuordnung ──\n")
meta[, .(igel,
         gefunden = ifelse(!is.na(csv_pfad), "✓", "✗ FEHLT"),
         datei    = basename(csv_pfad))] |> print()

# ──────────────────────────────────────────────────────────────
# HAUPTSCHLEIFE: Jeder Igel
# ──────────────────────────────────────────────────────────────

alle_summary    <- list()   # Sammelt Ergebnisse für Übersichtstabelle
fehlende_igel   <- list()   # Sammelt Igel mit Fehler/Warnung (für Qualitätscheck)
gamm_nacht_min  <- list()   # Nacht-Minutendaten für GAMM (Option B)
gamm_tage_sum   <- list()   # Tagessummary für GAMM (Option A)
chrono_min      <- list()   # Volle 24h-Minutendaten für Chronobiologie-Analyse

for (i in seq_len(nrow(meta))) {

  igel_info <- meta[i]
  igel_name <- igel_info$igel
  csv_pfad  <- igel_info$csv_pfad

  cat("\n══════════════════════════════════════\n")
  cat("Verarbeite:", igel_name, "\n")

  # ── Überspringe wenn kein CSV vorhanden ──
  if (is.na(csv_pfad)) {
    cat("  ⚠ Keine CSV-Datei gefunden — übersprungen.\n")
    fehlende_igel[[igel_name]] <- data.table(
      igel = igel_name, qc_status = "✗ FEHLER",
      qc_hinweis = "Keine CSV-Datei gefunden",
      n_tage = NA_integer_, tote_tage = NA_integer_,
      pct_aktiv_tag = NA_real_, pct_aktiv_nacht = NA_real_,
      nachtaktiv = NA_character_, pct_nacht_daten = NA_real_)
    next
  }

  hard_release <- igel_info$hard_release

  # ── 1. DATEN EINLESEN ──
  dt <- tryCatch(
    fread(csv_pfad, na.strings = c("", "NA"), showProgress = FALSE),
    error = function(e) {
      cat("  ✗ Fehler beim Einlesen:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0) next

  # Spalten umbenennen (tolerant gegenüber fehlenden Spalten)
  setnames(dt,
    old = c("Station Name", "Individual Name", "_time", "0"),
    new = c("station", "igel_col", "time_raw", "signal"),
    skip_absent = TRUE)

  cat("  Zeilen eingelesen:", format(nrow(dt), big.mark = "'"), "\n")

  # ── 2. ZEITSTEMPEL ──
  dt[, time_utc   := ymd_hms(time_raw, tz = "UTC", quiet = TRUE)]
  dt[, time_local := with_tz(time_utc, tzone = zeitzone)]
  dt[, time_min   := floor_date(time_local, unit = "minute")]

  # ── 3. KLASSIFIKATION ──
  # Verwendet: pred_nested_loio_smoothed_wmv (~20 Min Glaettung)
  # HINWEIS: Igel11 hat keine aktiven smoothed_wmv-Klassifikationen
  #   (Aktivitaetssignale zu fragmentiert/kurz fuer den 20-Min-Filter —
  #   genuines Merkmal dieses Tieres, kein technischer Fehler).
  #   Igel11 wird daher in Block 3 (Chronobiologie) ausgeschlossen.
  #   In Block 2 (GAMM) bleibt Igel11 enthalten (wmv-Klassifikation).
  if (!klasse_spalte %in% names(dt)) {
    cat("  ⚠ Spalte '", klasse_spalte, "' nicht gefunden — übersprungen.\n")
    next
  }
  dt[, klasse := fcase(
    get(klasse_spalte) == "a", "aktiv",
    get(klasse_spalte) == "p", "passiv",
    default = NA_character_
  )]

  # ── 4. MINUTENREDUKTION (Modalwert) ──
  # Pro Minute den häufigsten Wert (aktiv/passiv) nehmen
  dt_min <- dt[!is.na(klasse), .(
    klasse_modal = modal_val(klasse),
    signal_max   = max(signal, na.rm = TRUE),
    n_messungen  = .N
  ), by = time_min]

  cat("  Minuten nach Reduktion:", format(nrow(dt_min), big.mark = "'"), "\n")

  # ── 5. ZEITRAUM FILTERN ──
  dt_min[, datum := as.Date(time_min)]

  # Start: hartes Auswilderungsdatum aus Metadaten
  dt_min <- dt_min[datum >= hard_release]

  if (nrow(dt_min) == 0) {
    cat("  ⚠ Keine Daten nach dem Auswilderungsdatum — übersprungen.\n")
    fehlende_igel[[igel_name]] <- data.table(
      igel = igel_name, qc_status = "✗ FEHLER",
      qc_hinweis = "Keine Daten nach Auswilderungsdatum",
      n_tage = NA_integer_, tote_tage = NA_integer_,
      pct_aktiv_tag = NA_real_, pct_aktiv_nacht = NA_real_,
      nachtaktiv = NA_character_, pct_nacht_daten = NA_real_)
    next
  }

  # Ende: letztes Datum direkt aus der CSV-Datei ableiten
  # (Metadaten-Wert kann falsch sein — z.B. Igel9 hat 2022 als letztes Signal)
  letztes_datum_csv <- max(dt_min$datum, na.rm = TRUE)
  dt_min <- dt_min[datum <= letztes_datum_csv]   # (bereits der Fall, nur explizit)

  cat("  Zeitraum in CSV:", format(hard_release, "%d.%m.%Y"),
      "bis", format(letztes_datum_csv, "%d.%m.%Y"), "\n")
  cat("  Minuten nach Zeitraumfilter:", format(nrow(dt_min), big.mark = "'"), "\n")

  # Hinweis wenn Metadaten-Datum stark abweicht
  letztes_datum_meta <- igel_info$letztes_signal
  if (!is.na(letztes_datum_meta) &&
      abs(as.integer(letztes_datum_csv - letztes_datum_meta)) > 5) {
    cat("  ⚠ Metadaten-Enddatum (", format(letztes_datum_meta, "%d.%m.%Y"),
        ") weicht >5 Tage vom CSV-Enddatum (", format(letztes_datum_csv, "%d.%m.%Y"),
        ") ab.\n", sep = "")
  }

  # ── 5b. TOTER SENDER ERKENNEN ──
  # Strategie: "letzter sinnvoller Tag" = letzter Tag mit >= schwelle_aktiv_pct % Aktivität
  # Alle Tage danach werden als statisches Signal (toter/verlorener Sender) gewertet
  # und aus der Analyse entfernt.
  tagesdaten <- dt_min[!is.na(klasse_modal), .(
    n_min_roh  = .N,
    pct_aktiv  = mean(klasse_modal == "aktiv") * 100
  ), by = datum][order(datum)]

  aktive_tage <- tagesdaten[pct_aktiv >= schwelle_aktiv_pct, datum]

  tote_tage_entfernt <- 0L
  if (length(aktive_tage) == 0) {
    # Kein einziger Tag über Schwellenwert — alle Daten suspekt
    cat("  ⚠ Kein Tag mit >=", schwelle_aktiv_pct, "% Aktivität gefunden.",
        "Igel wird übersprungen.\n")
    fehlende_igel[[igel_name]] <- data.table(
      igel = igel_name, qc_status = "⚠⚠ KRITISCH",
      qc_hinweis = paste0("100% passiv über alle Tage — statisches Signal / toter Sender"),
      n_tage = dt_min[, uniqueN(datum)], tote_tage = NA_integer_,
      pct_aktiv_tag = NA_real_, pct_aktiv_nacht = NA_real_,
      nachtaktiv = NA_character_, pct_nacht_daten = NA_real_)
    next
  }

  letzter_aktiver_tag <- max(aktive_tage)

  if (letzter_aktiver_tag < letztes_datum_csv) {
    tote_tage_entfernt <- as.integer(letztes_datum_csv - letzter_aktiver_tag)
    cat("  ⚠ Toter-Sender-Erkennung:", tote_tage_entfernt,
        "Tag(e) mit <", schwelle_aktiv_pct, "% Aktivität am Ende entfernt\n")
    cat("    Neues Enddatum:", format(letzter_aktiver_tag, "%d.%m.%Y"), "\n")
    dt_min           <- dt_min[datum <= letzter_aktiver_tag]
    letztes_datum_csv <- letzter_aktiver_tag
  } else {
    cat("  ✓ Kein statisches Signal am Ende erkannt\n")
  }

  # ── 6. TAG / NACHT (NOAA-Algorithmus, kein suncalc) ──
  tage_vec <- as.Date(as.integer(unique(dt_min$datum)), origin = "1970-01-01")
  sun <- .sonnenzeiten(tage_vec, lat = standort_lat, lon = standort_lon, tz = zeitzone)

  dt_min <- merge(dt_min, sun[, .(date, sunrise, sunset)],
                  by.x = "datum", by.y = "date", all.x = TRUE)

  dt_min[, tageszeit := fifelse(
    time_min >= sunrise & time_min < sunset, "Tag", "Nacht"
  )]

  dt_min[, stunde := hour(time_min) + minute(time_min) / 60]

  # Mittlere Sonnenauf/-untergangszeit (für Peak-Berechnungen und Plots)
  sunrise_h <- mean(hour(dt_min$sunrise) + minute(dt_min$sunrise) / 60, na.rm = TRUE)
  sunset_h  <- mean(hour(dt_min$sunset)  + minute(dt_min$sunset)  / 60, na.rm = TRUE)

  # ── 7. TAGE NACH AUSWILDERUNG ──
  dt_min[, tage_seit := as.integer(datum - hard_release) + 1L]

  # ── 8. SUMMARY ──

  # 8a. Pro Tag und Tageszeit (für Plots)
  summary_igel <- dt_min[!is.na(klasse_modal), .(
    n_min     = .N,
    pct_aktiv = round(mean(klasse_modal == "aktiv") * 100, 1)
  ), by = .(datum, tageszeit, tage_seit)]

  # 8b. Gesamtdatenmenge Tag vs. Nacht
  daten_tn <- dt_min[!is.na(klasse_modal), .(n_min = .N), by = tageszeit]
  n_tag   <- daten_tn[tageszeit == "Tag",   n_min] %||% 0
  n_nacht <- daten_tn[tageszeit == "Nacht", n_min] %||% 0
  n_total <- n_tag + n_nacht
  pct_nacht_daten <- round(n_nacht / n_total * 100, 1)

  # 8c. Aktivitätsanteil Tag und Nacht getrennt
  aktiv_tn <- dt_min[!is.na(klasse_modal), .(
    pct_aktiv = round(mean(klasse_modal == "aktiv") * 100, 1)
  ), by = tageszeit]
  pct_aktiv_tag   <- aktiv_tn[tageszeit == "Tag",   pct_aktiv] %||% NA_real_
  pct_aktiv_nacht <- aktiv_tn[tageszeit == "Nacht", pct_aktiv] %||% NA_real_

  # 8d. 24h-Profil und Peak-Stunde
  profil_24h_igel <- dt_min[!is.na(klasse_modal), .(
    pct_aktiv = mean(klasse_modal == "aktiv") * 100
  ), by = .(stunde_rund = round(stunde) %% 24L)]

  # Peak gesamt (Stunde mit höchster Aktivität über 24h)
  peak_gesamt_h <- profil_24h_igel[which.max(pct_aktiv), stunde_rund]

  # Peak nachts (nur Nachtstunden)
  nacht_stunden <- profil_24h_igel[
    stunde_rund < round(sunrise_h) | stunde_rund >= round(sunset_h)]
  peak_nacht_h <- if (nrow(nacht_stunden) > 0)
    nacht_stunden[which.max(pct_aktiv), stunde_rund] else NA_integer_

  # Peak tagsüber
  tag_stunden <- profil_24h_igel[
    stunde_rund >= round(sunrise_h) & stunde_rund < round(sunset_h)]
  peak_tag_h <- if (nrow(tag_stunden) > 0)
    tag_stunden[which.max(pct_aktiv), stunde_rund] else NA_integer_

  # 8e. Hauptaktivitätsphase (längster aktiver Block über 24h)
  phase <- find_hauptphase(profil_24h_igel)

  # 8f. Flag: wenig Nachtdaten (< 30% der Gesamtdaten)
  wenig_nacht <- pct_nacht_daten < 30

  # 8g. Alles in eine Zeile für die Gesamttabelle
  sg <- data.table(
    igel                = igel_name,
    geschlecht          = igel_info$geschlecht %||% NA_character_,
    gewicht_g           = igel_info$gewicht_g  %||% NA_real_,
    diagnose            = igel_info$diagnose   %||% NA_character_,
    hard_release             = hard_release,
    letztes_signal_csv       = letztes_datum_csv,
    tote_sender_tage_entfernt = tote_tage_entfernt,
    n_tage                   = dt_min[, uniqueN(datum)],
    # ── Datenmenge ──
    n_min_tag           = n_tag,
    n_min_nacht         = n_nacht,
    pct_nacht_daten     = pct_nacht_daten,
    wenig_nachtdaten    = ifelse(wenig_nacht, "⚠ JA", "ok"),
    # ── Aktivität ──
    pct_aktiv_tag       = pct_aktiv_tag,
    pct_aktiv_nacht     = pct_aktiv_nacht,
    nachtaktiv          = fcase(
      pct_aktiv_nacht > pct_aktiv_tag, "ja",
      pct_aktiv_tag > pct_aktiv_nacht, "nein",
      default = "unklar"),
    # ── Peak-Stunden ──
    peak_uhrzeit_gesamt = sprintf("%02d:00", peak_gesamt_h),
    peak_uhrzeit_nacht  = sprintf("%02d:00", peak_nacht_h),
    peak_uhrzeit_tag    = sprintf("%02d:00", peak_tag_h),
    # ── Hauptaktivitätsphase ──
    hauptphase_beginn   = sprintf("%02d:00", phase$start_h),
    hauptphase_ende     = sprintf("%02d:00", phase$ende_h),
    hauptphase_dauer_h  = phase$dauer_h
  )
  alle_summary[[igel_name]] <- sg

  # ── Daten für GAMM-Analyse akkumulieren ──

  # Tagessummary (Option A): eine Zeile pro Igel × Datum × Tageszeit
  tage_tn <- dt_min[!is.na(klasse_modal), .(
    n_min     = .N,
    pct_aktiv = mean(klasse_modal == "aktiv"),   # Anteil 0–1 (für Beta-Regression)
    n_aktiv   = sum(klasse_modal == "aktiv")
  ), by = .(datum, tageszeit, tage_seit)]
  tage_tn[, igel := igel_name]
  gamm_tage_sum[[igel_name]] <- tage_tn

  # Minutendaten Nacht (Option B): eine Zeile pro Nachtminute
  nacht_min <- dt_min[tageszeit == "Nacht" & !is.na(klasse_modal), .(
    igel     = igel_name,
    datum,
    tage_seit,
    stunde   = stunde,                            # kontinuierlich 0–24
    aktiv    = as.integer(klasse_modal == "aktiv")
  )]
  gamm_nacht_min[[igel_name]] <- nacht_min

  # Minutendaten 24h (Chronobiologie): alle Tagesstunden für Actogramme & Cosinor
  alle_min_24h <- dt_min[!is.na(klasse_modal), .(
    igel      = igel_name,
    datum,
    tage_seit,
    tageszeit,
    stunde    = stunde,
    aktiv     = as.integer(klasse_modal == "aktiv")
  )]
  chrono_min[[igel_name]] <- alle_min_24h

  # ── QC-Eintrag für diesen Igel ──
  qc_hinweis_liste <- character(0)
  if (tote_tage_entfernt > 0)
    qc_hinweis_liste <- c(qc_hinweis_liste,
      paste0(tote_tage_entfernt, " Tage (stat. Sender) entfernt"))
  if (sg$n_tage < 7)
    qc_hinweis_liste <- c(qc_hinweis_liste, "< 7 Tage Daten")
  if (wenig_nacht)
    qc_hinweis_liste <- c(qc_hinweis_liste, "< 30% Nachtdaten")
  if (!is.na(pct_aktiv_tag) && !is.na(pct_aktiv_nacht) &&
      pct_aktiv_tag > pct_aktiv_nacht + 10)
    qc_hinweis_liste <- c(qc_hinweis_liste, "Scheinbar tagaktiv")

  qc_status_igel <- if (length(qc_hinweis_liste) == 0) {
    "✓ OK"
  } else if (any(grepl("tagaktiv", qc_hinweis_liste))) {
    "⚠⚠ KRITISCH"
  } else {
    "⚠ WARNUNG"
  }

  fehlende_igel[[igel_name]] <- data.table(
    igel          = igel_name,
    qc_status     = qc_status_igel,
    qc_hinweis    = if (length(qc_hinweis_liste) == 0) "—" else paste(qc_hinweis_liste, collapse = " | "),
    n_tage        = sg$n_tage,
    tote_tage     = tote_tage_entfernt,
    pct_aktiv_tag   = pct_aktiv_tag,
    pct_aktiv_nacht = pct_aktiv_nacht,
    nachtaktiv    = sg$nachtaktiv,
    pct_nacht_daten = pct_nacht_daten)

  # ── 9. PLOTS ──
  # sunrise_h und sunset_h bereits nach .sonnenzeiten-Merge berechnet

  sub_txt <- sprintf("%s | %s g | %s | Auswilderung: %s | Daten bis: %s (%d Tage)",
    igel_info$geschlecht %||% "?",
    igel_info$gewicht_g  %||% "?",
    igel_info$diagnose   %||% "?",
    format(hard_release, "%d.%m.%Y"),
    format(letztes_datum_csv, "%d.%m.%Y"),
    as.integer(letztes_datum_csv - hard_release) + 1L)

  # Plot A: Aktivitätsanteil Tag vs. Nacht über Zeit
  pA <- ggplot(summary_igel,
               aes(x = tage_seit, y = pct_aktiv,
                   color = tageszeit, group = tageszeit)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = farben) +
    scale_y_continuous(limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title    = paste(igel_name, "— Aktivitätsanteil Tag vs. Nacht"),
         subtitle = sub_txt,
         x = "Tage nach Auswilderung",
         y = "Anteil aktiver Minuten (%)",
         color = "Tageszeit") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8, color = "grey50"))

  # Plot B: 24h-Aktivitätsprofil
  profil_24h <- dt_min[!is.na(klasse_modal), .(
    pct_aktiv = mean(klasse_modal == "aktiv") * 100
  ), by = .(stunde_rund = round(stunde) %% 24L)]

  pB <- ggplot(profil_24h,
               aes(x    = stunde_rund, y = pct_aktiv,
                   fill = ifelse(stunde_rund < sunrise_h |
                                 stunde_rund >= sunset_h, "Nacht", "Tag"))) +
    annotate("rect", xmin = -0.5, xmax = sunrise_h,
             ymin = -Inf, ymax = Inf, fill = "#2C3E6B", alpha = 0.08) +
    annotate("rect", xmin = sunset_h, xmax = 23.5,
             ymin = -Inf, ymax = Inf, fill = "#2C3E6B", alpha = 0.08) +
    geom_col(alpha = 0.85, show.legend = TRUE) +
    scale_fill_manual(values = farben, name = "Tageszeit") +
    scale_x_continuous(breaks = seq(0, 23, 3),
                       labels = sprintf("%02d:00", seq(0, 23, 3))) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = paste(igel_name, "— 24h-Aktivitätsprofil"),
         x = "Uhrzeit (MEZ/CEST)", y = "Anteil aktiver Minuten (%)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          plot.title   = element_text(face = "bold"),
          legend.position = "bottom")

  # Plot C: Heatmap (Datum × Stunde)
  hm <- dt_min[!is.na(klasse_modal), .(
    pct_aktiv = mean(klasse_modal == "aktiv") * 100
  ), by = .(datum, stunde_rund = round(stunde) %% 24L)]

  pC <- ggplot(hm, aes(x = stunde_rund, y = as.factor(datum),
                        fill = pct_aktiv)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_vline(xintercept = sunrise_h, color = "gold",
               linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = sunset_h,  color = "darkorange",
               linewidth = 1, linetype = "dashed") +
    scale_fill_gradient2(low = "#264653", mid = "#e9c46a",
                         high = "#e76f51", midpoint = 50,
                         limits = c(0, 100), name = "% aktiv") +
    scale_x_continuous(breaks = seq(0, 23, 3),
                       labels = sprintf("%02d:00", seq(0, 23, 3))) +
    labs(title = paste(igel_name, "— Aktivitäts-Heatmap"),
         x = "Uhrzeit (MEZ/CEST)", y = "Datum") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold"))

  # Alle drei Plots speichern
  plot_datei <- file.path(output_ordner,
                           paste0(igel_name, "_aktivitaet.png"))
  png(plot_datei, width = 1600, height = 1900, res = 130)
  grid.arrange(pA, pB, pC, ncol = 1,
    top = grid::textGrob(
      paste("VHF Aktivitätsanalyse —", igel_name),
      gp = grid::gpar(fontsize = 14, fontface = "bold")))
  dev.off()

  cat("  ✓ Plot gespeichert:", basename(plot_datei), "\n")
}

# ──────────────────────────────────────────────────────────────
# GESAMTÜBERSICHT ALLER IGEL
# ──────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════\n")
cat("Erstelle Gesamtübersicht...\n")

summary_alle <- rbindlist(alle_summary, fill = TRUE)

cat("\n=== GESAMTÜBERSICHT ALLE IGEL ===\n")
print(summary_alle[order(hard_release)])

# ── Als formatierte Excel-Tabelle speichern ──
tbl <- summary_alle[order(hard_release)]

wb_out <- createWorkbook()
addWorksheet(wb_out, "Igel Übersicht")

# Spaltentitel auf Deutsch
col_labels <- c(
  "Igel", "Geschlecht", "Gewicht (g)", "Diagnose",
  "Auswilderung", "Letztes Signal (CSV)", "Entfernte Tage (stat. Sender)",
  "Anzahl Tage",
  "Minuten Tag", "Minuten Nacht", "% Nachtdaten",
  "Wenig Nachtdaten (<30%)",
  "% aktiv Tag", "% aktiv Nacht", "Nachtaktiv",
  "Peak Uhrzeit (gesamt)", "Peak Uhrzeit (Nacht)", "Peak Uhrzeit (Tag)",
  "Hauptphase Beginn", "Hauptphase Ende", "Hauptphase Dauer (h)"
)

# Styles
style_header <- createStyle(
  fontName = "Arial", fontSize = 10, fontColour = "white",
  fgFill = "#2C3E6B", halign = "center", valign = "center",
  textDecoration = "bold", wrapText = TRUE, border = "Bottom",
  borderColour = "white"
)
style_warnung <- createStyle(
  fontName = "Arial", fontSize = 10,
  fgFill = "#FDECEA", fontColour = "#C0392B", textDecoration = "bold"
)
style_ok <- createStyle(
  fontName = "Arial", fontSize = 10, fgFill = "#EAF4EA"
)
style_zahl <- createStyle(
  fontName = "Arial", fontSize = 10, numFmt = "0.0", halign = "center"
)
style_pct <- createStyle(
  fontName = "Arial", fontSize = 10, numFmt = "0.0\"%\"", halign = "center"
)
style_datum <- createStyle(
  fontName = "Arial", fontSize = 10, numFmt = "DD.MM.YYYY", halign = "center"
)
style_text_center <- createStyle(
  fontName = "Arial", fontSize = 10, halign = "center"
)
style_normal <- createStyle(fontName = "Arial", fontSize = 10)

# Header schreiben
writeData(wb_out, "Igel Übersicht", as.data.table(t(col_labels)),
          startRow = 1, startCol = 1, colNames = FALSE)
addStyle(wb_out, "Igel Übersicht", style_header,
         rows = 1, cols = seq_along(col_labels), gridExpand = TRUE)
setRowHeights(wb_out, "Igel Übersicht", rows = 1, heights = 38)

# Daten schreiben
writeData(wb_out, "Igel Übersicht", tbl,
          startRow = 2, startCol = 1, colNames = FALSE)

n_zeilen <- nrow(tbl)
daten_rows <- 2:(n_zeilen + 1)

# Grundformatierung alle Zellen
addStyle(wb_out, "Igel Übersicht", style_normal,
         rows = daten_rows, cols = 1:ncol(tbl), gridExpand = TRUE)

# Zahlenformat: Gewicht, Minuten, Tage
addStyle(wb_out, "Igel Übersicht", style_zahl,
         rows = daten_rows, cols = c(3, 7, 8, 9, 10, 21), gridExpand = TRUE)

# Prozentformat: % Nachtdaten, % aktiv Tag/Nacht
addStyle(wb_out, "Igel Übersicht", style_pct,
         rows = daten_rows, cols = c(11, 13, 14), gridExpand = TRUE)

# Datum-Format: Auswilderung, letztes Signal
addStyle(wb_out, "Igel Übersicht", style_datum,
         rows = daten_rows, cols = c(5, 6), gridExpand = TRUE)

# Zentriert: Uhrzeiten, nachtaktiv, Geschlecht
addStyle(wb_out, "Igel Übersicht", style_text_center,
         rows = daten_rows, cols = c(2, 15, 16, 17, 18, 19, 20), gridExpand = TRUE)

# Bedingte Formatierung: Warnung-Zeilen einfärben
for (z in seq_len(n_zeilen)) {
  excel_row <- z + 1L
  if (!is.na(tbl$wenig_nachtdaten[z]) && tbl$wenig_nachtdaten[z] == "⚠ JA") {
    addStyle(wb_out, "Igel Übersicht", style_warnung,
             rows = excel_row, cols = 12, gridExpand = FALSE)
  } else {
    addStyle(wb_out, "Igel Übersicht", style_ok,
             rows = excel_row, cols = 12, gridExpand = FALSE)
  }
  # Nachtaktiv einfärben
  if (!is.na(tbl$nachtaktiv[z])) {
    farbe_nacht <- switch(tbl$nachtaktiv[z],
      "ja"     = createStyle(fontName="Arial", fontSize=10,
                             fgFill="#D5E8D4", fontColour="#2D6A2D",
                             halign="center", textDecoration="bold"),
      "nein"   = createStyle(fontName="Arial", fontSize=10,
                             fgFill="#FFE6CC", fontColour="#8B4500",
                             halign="center", textDecoration="bold"),
      "unklar" = createStyle(fontName="Arial", fontSize=10,
                             fgFill="#F5F5F5", fontColour="#666666",
                             halign="center"),
      style_text_center
    )
    addStyle(wb_out, "Igel Übersicht", farbe_nacht,
             rows = excel_row, cols = 15, gridExpand = FALSE)
  }
}

# Spaltenbreiten
setColWidths(wb_out, "Igel Übersicht", cols = 1:ncol(tbl),
             widths = c(8, 10, 10, 22, 14, 16, 18, 10,
                        11, 12, 11, 14, 12, 13, 10,
                        14, 14, 13, 14, 13, 13))

# Erste Zeile einfrieren
freezePane(wb_out, "Igel Übersicht", firstRow = TRUE)

# Als Tabelle formatieren
addFilter(wb_out, "Igel Übersicht", row = 1, cols = 1:ncol(tbl))

# ──────────────────────────────────────────────────────────────
# QUALITÄTSCHECK-SHEET direkt zu wb_out hinzufügen
# (vor dem einzigen saveWorkbook-Aufruf)
# ──────────────────────────────────────────────────────────────

cat("\nErstelle Qualitätscheck-Sheet...\n")

# Alle Einträge zusammenführen (erfolgreiche + fehlgeschlagene Igel)
qc_alle <- rbindlist(fehlende_igel, fill = TRUE)
# Nach Igelnummer sortieren (Igel1, Igel2, …)
qc_alle[, nr := as.integer(gsub("[^0-9]", "", igel))]
setorder(qc_alle, nr)
qc_alle[, nr := NULL]

# Spaltenbezeichnungen für QC
qc_col_labels <- c(
  "Igel", "Status", "Hinweise",
  "Anzahl Tage", "Entfernte Sender-Tage",
  "% aktiv Tag", "% aktiv Nacht", "Nachtaktiv", "% Nachtdaten"
)

# Sheets direkt in wb_out ergänzen — kein erneutes Laden nötig
addWorksheet(wb_out, "Qualitätscheck")
addWorksheet(wb_out, "Legende")

# ── Styles für QC ──
qc_style_header <- createStyle(
  fontName = "Arial", fontSize = 10, fontColour = "white",
  fgFill = "#1a1a2e", halign = "center", valign = "center",
  textDecoration = "bold", wrapText = TRUE
)
qc_style_ok <- createStyle(
  fontName = "Arial", fontSize = 10, fgFill = "#D5F5D5", fontColour = "#1a5c1a"
)
qc_style_warn <- createStyle(
  fontName = "Arial", fontSize = 10, fgFill = "#FFF3CD", fontColour = "#7d5a00"
)
qc_style_kritisch <- createStyle(
  fontName = "Arial", fontSize = 10, fgFill = "#FDECEA", fontColour = "#8B0000",
  textDecoration = "bold"
)
qc_style_fehler <- createStyle(
  fontName = "Arial", fontSize = 10, fgFill = "#F2DEDE", fontColour = "#C0392B"
)
qc_style_pct <- createStyle(
  fontName = "Arial", fontSize = 10, numFmt = "0.0\"%\"", halign = "center"
)
qc_style_center <- createStyle(
  fontName = "Arial", fontSize = 10, halign = "center"
)
qc_style_normal <- createStyle(fontName = "Arial", fontSize = 10)

# Header schreiben
writeData(wb_out, "Qualitätscheck", as.data.table(t(qc_col_labels)),
          startRow = 1, startCol = 1, colNames = FALSE)
addStyle(wb_out, "Qualitätscheck", qc_style_header,
         rows = 1, cols = seq_along(qc_col_labels), gridExpand = TRUE)
setRowHeights(wb_out, "Qualitätscheck", rows = 1, heights = 36)

# Daten schreiben
writeData(wb_out, "Qualitätscheck", qc_alle,
          startRow = 2, startCol = 1, colNames = FALSE)

# Zeilenweise einfärben je nach Status
n_qc <- nrow(qc_alle)
for (z in seq_len(n_qc)) {
  excel_row <- z + 1L
  status_z  <- qc_alle$qc_status[z]

  zeil_style <- switch(status_z,
    "✓ OK"        = qc_style_ok,
    "⚠ WARNUNG"   = qc_style_warn,
    "⚠⚠ KRITISCH" = qc_style_kritisch,
    "✗ FEHLER"    = qc_style_fehler,
    qc_style_normal
  )
  addStyle(wb_out, "Qualitätscheck", zeil_style,
           rows = excel_row, cols = 1:length(qc_col_labels), gridExpand = TRUE)
  # Prozentzahlen zentriert
  addStyle(wb_out, "Qualitätscheck", qc_style_pct,
           rows = excel_row, cols = c(6, 7, 9), gridExpand = TRUE)
  addStyle(wb_out, "Qualitätscheck", qc_style_center,
           rows = excel_row, cols = c(4, 5, 8), gridExpand = TRUE)
}

# Spaltenbreiten QC
setColWidths(wb_out, "Qualitätscheck",
             cols = 1:length(qc_col_labels),
             widths = c(10, 13, 45, 12, 18, 12, 13, 11, 13))
freezePane(wb_out, "Qualitätscheck", firstRow = TRUE)
addFilter(wb_out, "Qualitätscheck", row = 1, cols = 1:length(qc_col_labels))

# ── Legende ──
legende_daten <- data.table(
  Status = c("✓ OK", "⚠ WARNUNG", "⚠⚠ KRITISCH", "✗ FEHLER"),
  Bedeutung = c(
    "Keine Auffälligkeiten — Daten können direkt interpretiert werden",
    "Kleinere Einschränkungen (z.B. kurze Überwachung, tote Sendertage entfernt, wenig Nachtdaten) — mit Vorsicht interpretieren",
    "Starke Auffälligkeit — z.B. scheinbar tagaktiv oder 100% passiv — Einzelfallprüfung notwendig",
    "Kritischer Fehler — Igel konnte nicht analysiert werden (fehlende CSV oder keine Daten nach Auswilderung)"
  )
)

legende_header_style <- createStyle(
  fontName = "Arial", fontSize = 11, textDecoration = "bold",
  fgFill = "#1a1a2e", fontColour = "white"
)
writeData(wb_out, "Legende", legende_daten, startRow = 1, startCol = 1, colNames = TRUE)
addStyle(wb_out, "Legende", legende_header_style,
         rows = 1, cols = 1:2, gridExpand = TRUE)

legende_farben <- list(
  "✓ OK"        = qc_style_ok,
  "⚠ WARNUNG"   = qc_style_warn,
  "⚠⚠ KRITISCH" = qc_style_kritisch,
  "✗ FEHLER"    = qc_style_fehler
)
for (z in 1:4) {
  addStyle(wb_out, "Legende", legende_farben[[z]],
           rows = z + 1L, cols = 1:2, gridExpand = TRUE)
}
setColWidths(wb_out, "Legende", cols = 1:2, widths = c(16, 90))
setRowHeights(wb_out, "Legende", rows = 2:5, heights = 36)

# ──────────────────────────────────────────────────────────────
# DASHBOARD-SHEET  (Kennzahlen auf einen Blick)
# ──────────────────────────────────────────────────────────────

cat("\nErstelle Dashboard-Sheet...\n")

addWorksheet(wb_out, "Dashboard", gridLines = FALSE)

# Hilfsstyle-Funktion (lokal)
ds_style <- function(fg = "FFFFFF", fc = "000000", bold = FALSE,
                     halign = "center", size = 10, italic = FALSE,
                     border = FALSE, wrapText = FALSE) {
  add_hash <- function(x) ifelse(startsWith(x, "#"), x, paste0("#", x))
  fg <- add_hash(fg); fc <- add_hash(fc)
  td <- c(if (bold) "bold", if (italic) "italic")
  if (length(td) == 0) td <- NULL
  createStyle(fontName = "Arial", fontSize = size, fontColour = fc,
              fgFill = fg, textDecoration = td,
              halign = halign, valign = "center",
              wrapText = wrapText,
              border = if (border) "TopBottomLeftRight" else NULL,
              borderColour = if (border) "#AABFD4" else NULL)
}

# Zeilen-/Spaltenbreiten
setColWidths(wb_out, "Dashboard",
             cols = 1:10, widths = c(22, 13, 13, 13, 13, 13, 13, 13, 13, 5))

# ── Titelzeile ─────────────────────────────────────────────────
mergeCells(wb_out, "Dashboard", cols = 1:9, rows = 1)
writeData(wb_out, "Dashboard",
  paste0("Block 0 — Aktivitätsdaten-Übersicht  |  VHF Igelbesenderung  |  Stand: ",
         format(Sys.Date(), "%d.%m.%Y")),
  startRow = 1, startCol = 1)
addStyle(wb_out, "Dashboard", ds_style("1a2e4a", "FFFFFF", bold = TRUE,
         halign = "left", size = 14), rows = 1, cols = 1:9, gridExpand = TRUE)
setRowHeights(wb_out, "Dashboard", rows = 1, heights = 28)

# ── Kennzahlen-Zeile ───────────────────────────────────────────
n_proc    <- nrow(summary_alle)
n_ok      <- sum(qc_alle$qc_status == "✓ OK",        na.rm = TRUE)
n_warn    <- sum(qc_alle$qc_status == "⚠ WARNUNG",   na.rm = TRUE)
n_krit    <- sum(qc_alle$qc_status == "⚠⚠ KRITISCH", na.rm = TRUE)
n_fehler  <- sum(qc_alle$qc_status == "✗ FEHLER",    na.rm = TRUE)
n_nacht_j <- sum(summary_alle$nachtaktiv == "ja",     na.rm = TRUE)
n_nacht_n <- sum(summary_alle$nachtaktiv == "nein",   na.rm = TRUE)

mergeCells(wb_out, "Dashboard", cols = 1:9, rows = 2)
writeData(wb_out, "Dashboard",
  paste0("N = ", n_proc, " Tiere verarbeitet  |  ",
         n_ok, " OK  |  ", n_warn, " Warnungen  |  ",
         n_fehler + n_krit, " Fehler/Kritisch  |  ",
         n_nacht_j, " nachtaktiv  |  ", n_nacht_n, " tagaktiv/unklar"),
  startRow = 2, startCol = 1)
addStyle(wb_out, "Dashboard", ds_style("EBF3FB", "444444", italic = TRUE,
         halign = "left", size = 10), rows = 2, cols = 1:9, gridExpand = TRUE)
setRowHeights(wb_out, "Dashboard", rows = 2, heights = 16)
setRowHeights(wb_out, "Dashboard", rows = 3, heights = 6)

# ── Abschnitt 1: Stichprobe nach Geschlecht ───────────────────
setRowHeights(wb_out, "Dashboard", rows = 4, heights = 20)
mergeCells(wb_out, "Dashboard", cols = 1:9, rows = 4)
writeData(wb_out, "Dashboard", "Stichprobe nach Geschlecht", startRow = 4, startCol = 1)
addStyle(wb_out, "Dashboard", ds_style("2C5F8A", "FFFFFF", bold = TRUE,
         halign = "left", size = 11), rows = 4, cols = 1:9, gridExpand = TRUE)

# Header
s1_hdrs <- c("Gruppe", "N", "Ø Tracking (d)", "Min. (d)", "Max. (d)",
              "Ø % aktiv Tag", "Ø % aktiv Nacht", "N nachtaktiv", "N warnung")
writeData(wb_out, "Dashboard", as.data.frame(t(s1_hdrs)), startRow = 5, startCol = 1, colNames = FALSE)
addStyle(wb_out, "Dashboard", ds_style("2C5F8A", "FFFFFF", bold = TRUE, border = TRUE),
         rows = 5, cols = 1:9, gridExpand = TRUE)
setRowHeights(wb_out, "Dashboard", rows = 5, heights = 18)

ds_row <- 6
for (grp in c("Male", "Female", "GESAMT")) {
  if (grp == "GESAMT") {
    sub <- summary_alle
    qc_sub <- qc_alle
    lbl <- "Gesamt"
    bg  <- "D6E4F0"
  } else {
    sub <- summary_alle[geschlecht == grp]
    qc_sub <- qc_alle[igel %in% sub$igel]
    lbl <- if (grp == "Male") "Männchen" else "Weibchen"
    bg  <- if (grp == "Male") "E3F0FB" else "FCE8EC"
  }
  if (nrow(sub) == 0) next
  tp_v <- sub$n_tage[!is.na(sub$n_tage)]
  row_data <- data.frame(
    Gruppe   = lbl,
    N        = nrow(sub),
    Tracking = round(mean(tp_v), 1),
    Min      = if (length(tp_v)) min(tp_v) else NA,
    Max      = if (length(tp_v)) max(tp_v) else NA,
    PctTag   = round(mean(sub$pct_aktiv_tag[!is.na(sub$pct_aktiv_tag)]), 1),
    PctNacht = round(mean(sub$pct_aktiv_nacht[!is.na(sub$pct_aktiv_nacht)]), 1),
    Nachtaktiv = sum(sub$nachtaktiv == "ja", na.rm = TRUE),
    Warnung  = sum(qc_sub$qc_status %in% c("⚠ WARNUNG","⚠⚠ KRITISCH"), na.rm = TRUE)
  )
  writeData(wb_out, "Dashboard", row_data, startRow = ds_row, startCol = 1, colNames = FALSE)
  addStyle(wb_out, "Dashboard", ds_style(bg, "222222", halign = "left",  border = TRUE),
           rows = ds_row, cols = 1, gridExpand = TRUE)
  addStyle(wb_out, "Dashboard", ds_style(bg, "222222", halign = "center", border = TRUE),
           rows = ds_row, cols = 2:9, gridExpand = TRUE)
  if (grp == "GESAMT")
    addStyle(wb_out, "Dashboard", ds_style(bg, "222222", bold = TRUE, halign = "center", border = TRUE),
             rows = ds_row, cols = 1:9, gridExpand = TRUE)
  setRowHeights(wb_out, "Dashboard", rows = ds_row, heights = 16)
  ds_row <- ds_row + 1
}

setRowHeights(wb_out, "Dashboard", rows = ds_row, heights = 6); ds_row <- ds_row + 1

# ── Abschnitt 2: Aktivitätskennzahlen ─────────────────────────
mergeCells(wb_out, "Dashboard", cols = 1:9, rows = ds_row)
writeData(wb_out, "Dashboard", "Aktivitätskennzahlen (alle Tiere mit Daten)",
          startRow = ds_row, startCol = 1)
addStyle(wb_out, "Dashboard", ds_style("2C5F8A", "FFFFFF", bold = TRUE,
         halign = "left", size = 11), rows = ds_row, cols = 1:9, gridExpand = TRUE)
setRowHeights(wb_out, "Dashboard", rows = ds_row, heights = 20)
ds_row <- ds_row + 1

# Aktivitäts-Header
ak_hdrs <- c("Kennzahl", "Mittelwert", "Median", "Min.", "Max.", "SD",
              "N Tiere", "", "")
writeData(wb_out, "Dashboard", as.data.frame(t(ak_hdrs)), startRow = ds_row, startCol = 1, colNames = FALSE)
addStyle(wb_out, "Dashboard", ds_style("2C5F8A", "FFFFFF", bold = TRUE, border = TRUE),
         rows = ds_row, cols = 1:9, gridExpand = TRUE)
setRowHeights(wb_out, "Dashboard", rows = ds_row, heights = 18)
ds_row <- ds_row + 1

ak_vars <- list(
  list("% aktiv Tag",          summary_alle$pct_aktiv_tag),
  list("% aktiv Nacht",        summary_alle$pct_aktiv_nacht),
  list("% Nachtdaten",         summary_alle$pct_nacht_daten),
  list("Tracking-Dauer (d)",   summary_alle$n_tage),
  list("Hauptphase Dauer (h)", summary_alle$hauptphase_dauer_h)
)
for (i in seq_along(ak_vars)) {
  item <- ak_vars[[i]]
  v    <- as.numeric(item[[2]][!is.na(item[[2]])])
  bg   <- if (i %% 2 == 0) "F5F5F5" else "FFFFFF"
  row_d <- data.frame(
    Lbl  = item[[1]],
    Mean = if (length(v)) round(mean(v), 1) else NA,
    Med  = if (length(v)) round(median(v), 1) else NA,
    Min  = if (length(v)) round(min(v), 1) else NA,
    Max  = if (length(v)) round(max(v), 1) else NA,
    SD   = if (length(v)) round(sd(v), 1) else NA,
    N    = length(v),
    E1   = "", E2 = ""
  )
  writeData(wb_out, "Dashboard", row_d, startRow = ds_row, startCol = 1, colNames = FALSE)
  addStyle(wb_out, "Dashboard", ds_style(bg, "222222", halign = "left",  border = TRUE),
           rows = ds_row, cols = 1, gridExpand = TRUE)
  addStyle(wb_out, "Dashboard", ds_style(bg, "222222", halign = "center", border = TRUE),
           rows = ds_row, cols = 2:9, gridExpand = TRUE)
  setRowHeights(wb_out, "Dashboard", rows = ds_row, heights = 16)
  ds_row <- ds_row + 1
}

setRowHeights(wb_out, "Dashboard", rows = ds_row, heights = 6); ds_row <- ds_row + 1

# ── Abschnitt 3: QC-Status Übersicht ──────────────────────────
mergeCells(wb_out, "Dashboard", cols = 1:9, rows = ds_row)
writeData(wb_out, "Dashboard", "Qualitätsstatus", startRow = ds_row, startCol = 1)
addStyle(wb_out, "Dashboard", ds_style("2C5F8A", "FFFFFF", bold = TRUE,
         halign = "left", size = 11), rows = ds_row, cols = 1:9, gridExpand = TRUE)
setRowHeights(wb_out, "Dashboard", rows = ds_row, heights = 20)
ds_row <- ds_row + 1

qc_stat_hdrs <- c("Status", "N", "% aller Tiere", "Bedeutung", "", "", "", "", "")
writeData(wb_out, "Dashboard", as.data.frame(t(qc_stat_hdrs)),
          startRow = ds_row, startCol = 1, colNames = FALSE)
addStyle(wb_out, "Dashboard", ds_style("2C5F8A", "FFFFFF", bold = TRUE, border = TRUE),
         rows = ds_row, cols = 1:9, gridExpand = TRUE)
setRowHeights(wb_out, "Dashboard", rows = ds_row, heights = 18)
ds_row <- ds_row + 1

qc_stat_items <- list(
  list("✓ OK",        n_ok,     "D5F5D5", "1a5c1a",
       "Keine Auffälligkeiten"),
  list("⚠ WARNUNG",   n_warn,   "FFF3CD", "7d5a00",
       "Kurze Überwachung, tote Sendertage entfernt oder wenig Nachtdaten"),
  list("⚠⚠ KRITISCH", n_krit,   "FDECEA", "8B0000",
       "Scheinbar tagaktiv oder 100% passiv — Einzelfallprüfung"),
  list("✗ FEHLER",    n_fehler, "F2DEDE", "C0392B",
       "Keine CSV-Datei oder keine Daten nach Auswilderung")
)
n_all_qc <- nrow(qc_alle)
for (item in qc_stat_items) {
  pct <- if (n_all_qc > 0) round(item[[2]] / n_all_qc * 100, 1) else 0
  row_d <- data.frame(Status=item[[1]], N=item[[2]],
                      Pct=paste0(pct, "%"), Bed=item[[5]],
                      E1="",E2="",E3="",E4="",E5="")
  writeData(wb_out, "Dashboard", row_d, startRow=ds_row, startCol=1, colNames=FALSE)
  st <- ds_style(item[[3]], item[[4]], border=TRUE)
  addStyle(wb_out, "Dashboard", st, rows=ds_row, cols=1:9, gridExpand=TRUE)
  addStyle(wb_out, "Dashboard",
           ds_style(item[[3]], item[[4]], halign="left", border=TRUE),
           rows=ds_row, cols=c(1,4), gridExpand=TRUE)
  setRowHeights(wb_out, "Dashboard", rows=ds_row, heights=16)
  ds_row <- ds_row + 1
}

setRowHeights(wb_out, "Dashboard", rows=ds_row, heights=6); ds_row <- ds_row + 1

# ── Abschnitt 4: Kompakte Tier-Tabelle ────────────────────────
mergeCells(wb_out, "Dashboard", cols=1:9, rows=ds_row)
writeData(wb_out, "Dashboard", "Kompakte Tier-Übersicht", startRow=ds_row, startCol=1)
addStyle(wb_out, "Dashboard", ds_style("2C5F8A", "FFFFFF", bold=TRUE,
         halign="left", size=11), rows=ds_row, cols=1:9, gridExpand=TRUE)
setRowHeights(wb_out, "Dashboard", rows=ds_row, heights=20)
ds_row <- ds_row + 1

tier_hdrs <- c("Tier-ID", "Geschlecht", "Tracking (d)", "% aktiv Tag",
               "% aktiv Nacht", "Nachtaktiv", "Peak Nacht", "QC-Status", "Hinweis")
writeData(wb_out, "Dashboard", as.data.frame(t(tier_hdrs)),
          startRow=ds_row, startCol=1, colNames=FALSE)
addStyle(wb_out, "Dashboard", ds_style("2C5F8A","FFFFFF",bold=TRUE,border=TRUE),
         rows=ds_row, cols=1:9, gridExpand=TRUE)
setRowHeights(wb_out, "Dashboard", rows=ds_row, heights=18)
ds_row <- ds_row + 1

tbl_kompakt <- merge(
  summary_alle[, .(igel, geschlecht, n_tage, pct_aktiv_tag, pct_aktiv_nacht,
                    nachtaktiv, peak_uhrzeit_nacht)],
  qc_alle[, .(igel, qc_status, qc_hinweis)],
  by = "igel", all.x = TRUE
)
tbl_kompakt[, nr := as.integer(gsub("[^0-9]","", igel))]
setorder(tbl_kompakt, nr); tbl_kompakt[, nr := NULL]

for (i in seq_len(nrow(tbl_kompakt))) {
  r <- tbl_kompakt[i]
  bg <- switch(r$qc_status %||% "—",
    "✓ OK"        = "FFFFFF",
    "⚠ WARNUNG"   = "FFF9E6",
    "⚠⚠ KRITISCH" = "FFF0EE",
    "✗ FEHLER"    = "FFE8E8",
    if (i %% 2 == 0) "F5F5F5" else "FFFFFF"
  )
  na_bg <- if (i %% 2 == 0) "F5F5F5" else "FFFFFF"
  row_d <- data.frame(
    Igel    = r$igel,
    Sex     = r$geschlecht %||% "—",
    Tage    = r$n_tage,
    PctTag  = r$pct_aktiv_tag,
    PctNacht= r$pct_aktiv_nacht,
    Nacht   = r$nachtaktiv %||% "—",
    Peak    = r$peak_uhrzeit_nacht %||% "—",
    Status  = r$qc_status %||% "—",
    Hinweis = r$qc_hinweis %||% "—"
  )
  writeData(wb_out, "Dashboard", row_d, startRow=ds_row, startCol=1, colNames=FALSE)
  addStyle(wb_out, "Dashboard", ds_style(bg, "222222", border=TRUE, halign="center"),
           rows=ds_row, cols=1:9, gridExpand=TRUE)
  addStyle(wb_out, "Dashboard", ds_style(bg, "222222", border=TRUE, halign="left"),
           rows=ds_row, cols=c(1,9), gridExpand=TRUE)
  # Nachtaktiv-Farbe
  nacht_bg <- switch(r$nachtaktiv %||% "—",
    "ja"="D5E8D4", "nein"="FFE6CC", "EEF0F0")
  nacht_fc <- switch(r$nachtaktiv %||% "—",
    "ja"="2D6A2D", "nein"="8B4500", "666666")
  addStyle(wb_out, "Dashboard",
           ds_style(nacht_bg, nacht_fc, bold=TRUE, border=TRUE),
           rows=ds_row, cols=6, gridExpand=TRUE)
  setRowHeights(wb_out, "Dashboard", rows=ds_row, heights=15)
  ds_row <- ds_row + 1
}

cat("✓ Dashboard-Sheet erstellt\n")

# ── Einmaliges Speichern aller Sheets ──
xlsx_pfad <- file.path(output_ordner, "igel_alle_summary.xlsx")
saveWorkbook(wb_out, xlsx_pfad, overwrite = TRUE)
cat("\n✓ Excel gespeichert:", basename(xlsx_pfad),
    "(Sheets: Igel Übersicht | Qualitätscheck | Legende)\n")
cat("  Übersicht-Spalten:", paste(col_labels, collapse = " | "), "\n")

# ──────────────────────────────────────────────────────────────
# ÜBERSICHTSPLOT: ALLE IGEL NEBENEINANDER
# ──────────────────────────────────────────────────────────────

# Übersichtsplot: Tag vs. Nacht für alle Igel — ins Long-Format bringen
uebersicht_long <- melt(
  summary_alle[, .(igel, hard_release, pct_aktiv_tag, pct_aktiv_nacht,
                   wenig_nachtdaten)],
  id.vars       = c("igel", "hard_release", "wenig_nachtdaten"),
  measure.vars  = c("pct_aktiv_tag", "pct_aktiv_nacht"),
  variable.name = "tageszeit",
  value.name    = "pct_aktiv"
)
uebersicht_long[, tageszeit := fifelse(tageszeit == "pct_aktiv_tag", "Tag", "Nacht")]

p_alle <- ggplot(
    uebersicht_long[!is.na(pct_aktiv)],
    aes(x    = reorder(igel, hard_release),
        y    = pct_aktiv,
        fill = tageszeit,
        alpha = wenig_nachtdaten == "⚠ JA")) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = farben, name = "Tageszeit") +
  scale_alpha_manual(values = c("TRUE" = 0.45, "FALSE" = 0.88),
                     name = "Wenig Nachtdaten", labels = c("nein", "ja (< 30%)")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  labs(
    title    = "Alle Igel — Aktivitätsanteil Tag vs. Nacht",
    subtitle = "Transparente Balken = weniger als 30% Nachtdaten (Vorsicht bei Interpretation)",
    x = NULL, y = "Anteil aktiver Minuten (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey50"),
        legend.position = "bottom")

ggsave(file.path(output_ordner, "igel_ALLE_uebersicht.png"),
       p_alle, width = 14, height = 6, dpi = 150)
cat("✓ Übersichtsplot gespeichert: igel_ALLE_uebersicht.png\n")

# ──────────────────────────────────────────────────────────────
# DESKRIPTIVE PLOTS  (6-Panel-Übersicht)
# ──────────────────────────────────────────────────────────────

cat("\nErstelle deskriptive Übersichtsplots...\n")

# Datenbasis: summary_alle + Metadaten zusammenführen
# Diagnosegruppen vereinfachen (wie Block1)
meta_plot <- copy(meta)
meta_plot[, diag_gruppe := fcase(
  grepl("orphan",  tolower(diagnose)),                          "Waise",
  grepl("trauma",  tolower(diagnose)) & !grepl("blind", tolower(diagnose)), "Trauma",
  grepl("blind",   tolower(diagnose)),                          "Blindheit",
  grepl("fungal",  tolower(diagnose)),                          "Pilz",
  grepl("ecto|endo|parasit", tolower(diagnose)),                "Parasiten",
  default = "Sonstige"
)]
meta_plot[, saison := fcase(
  month(hard_release) %in% 3:5,  "Frühling",
  month(hard_release) %in% 6:8,  "Sommer",
  month(hard_release) %in% 9:11, "Herbst",
  default = "Winter"
)]
meta_plot[, alter := ifelse(
  grepl("orphan", tolower(diagnose)) |
    (!is.na(gewicht_g) & as.numeric(gewicht_g) < 300),
  "Jungtier", "Adult"
)]

# Merge
sa <- merge(summary_alle,
            meta_plot[, .(igel, diag_gruppe, saison, alter)],
            by = "igel", all.x = TRUE)

# Gemeinsame Farbpaletten
pal_sex   <- c("Male" = "#2C7BB6", "Female" = "#D7191C")
pal_diag  <- c("Parasiten" = "#74C476", "Trauma" = "#FC8D59",
               "Waise"     = "#9E9AC8", "Blindheit" = "#BDBDBD",
               "Pilz"      = "#FDAE61", "Sonstige"  = "#EEEEEE")
pal_nacht <- c("ja" = "#2D6A2D", "nein" = "#C0392B", "unklar" = "#888888")
pal_saison <- c("Frühling" = "#74C476", "Sommer" = "#E8A838",
                "Herbst"   = "#FC8D59", "Winter" = "#4292C6")

theme_desk <- theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, color = "grey50"),
        panel.grid.minor = element_blank())

# ── Plot D1: Beobachtungsdauer je Tier (sortiert) ─────────────
sa_sort <- sa[!is.na(n_tage)][order(n_tage)]
sa_sort[, igel_f := factor(igel, levels = igel)]

p_d1 <- ggplot(sa_sort, aes(x = n_tage, y = igel_f, fill = geschlecht)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_vline(xintercept = 7, linetype = "dashed",
             color = "#C0392B", linewidth = 0.8) +
  annotate("text", x = 7.5, y = 1.5, label = "7 d\n(IS-Minimum)",
           color = "#C0392B", size = 2.8, hjust = 0) +
  scale_fill_manual(values = pal_sex, name = "Geschlecht",
                    labels = c("Male" = "Männchen", "Female" = "Weibchen")) +
  scale_x_continuous(breaks = c(0, 7, 14, 21, 30, 46)) +
  labs(title    = "Beobachtungsdauer nach Auswilderung",
       subtitle = "Rote Linie = Mindestdauer für IS-Analyse (7 Tage)",
       x = "Tage", y = NULL) +
  theme_desk +
  theme(axis.text.y = element_text(size = 8))

# ── Plot D2: Aktivität Tag vs. Nacht (Scatter) ────────────────
sa_act <- sa[!is.na(pct_aktiv_tag) & !is.na(pct_aktiv_nacht)]

p_d2 <- ggplot(sa_act, aes(x = pct_aktiv_tag, y = pct_aktiv_nacht,
                            color = nachtaktiv, label = igel)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.7) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text(vjust = -0.6, size = 2.5, alpha = 0.7) +
  scale_color_manual(values = pal_nacht, name = "Nachtaktiv",
                     labels = c("ja"="Ja","nein"="Nein","unklar"="Unklar")) +
  labs(title    = "Aktivität: Tag vs. Nacht",
       subtitle = "Über der Diagonalen = mehr nachts aktiv",
       x = "% aktiv tagsüber", y = "% aktiv nachts") +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  theme_desk

# ── Plot D3: % Nachtaktivität nach Geschlecht (Boxplot) ───────
p_d3 <- ggplot(sa[!is.na(pct_aktiv_nacht) & !is.na(geschlecht)],
               aes(x = geschlecht, y = pct_aktiv_nacht, fill = geschlecht)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.8, aes(color = nachtaktiv)) +
  scale_fill_manual(values  = pal_sex, guide = "none") +
  scale_color_manual(values = pal_nacht, name = "Nachtaktiv",
                     labels = c("ja"="Ja","nein"="Nein","unklar"="Unklar")) +
  scale_x_discrete(labels = c("Male"="Männchen","Female"="Weibchen")) +
  labs(title    = "Nachtaktivität nach Geschlecht",
       subtitle = "Anteil aktiver Minuten in der Nacht",
       x = NULL, y = "% aktiv (Nacht)") +
  theme_desk

# ── Plot D4: % Nachtaktivität nach Diagnose (Boxplot) ─────────
p_d4 <- ggplot(sa[!is.na(pct_aktiv_nacht) & !is.na(diag_gruppe)],
               aes(x = reorder(diag_gruppe, pct_aktiv_nacht, median, na.rm=TRUE),
                   y = pct_aktiv_nacht, fill = diag_gruppe)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = pal_diag, guide = "none") +
  labs(title    = "Nachtaktivität nach Diagnosegruppe",
       subtitle = "Geordnet nach Median",
       x = NULL, y = "% aktiv (Nacht)") +
  theme_desk +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ── Plot D5: Peak-Aktivitätszeit (Histogramm) ─────────────────
sa_peak <- sa[!is.na(peak_uhrzeit_nacht) & peak_uhrzeit_nacht != "NA:00"]
sa_peak[, peak_h := as.integer(sub(":.*", "", peak_uhrzeit_nacht))]
sa_peak <- sa_peak[!is.na(peak_h)]

p_d5 <- ggplot(sa_peak, aes(x = peak_h, fill = geschlecht)) +
  geom_histogram(binwidth = 2, alpha = 0.85, color = "white", position = "stack") +
  annotate("rect", xmin = -0.5, xmax = 6, ymin = -Inf, ymax = Inf,
           fill = "#2C3E6B", alpha = 0.06) +
  annotate("rect", xmin = 20, xmax = 24.5, ymin = -Inf, ymax = Inf,
           fill = "#2C3E6B", alpha = 0.06) +
  scale_fill_manual(values = pal_sex, name = "Geschlecht",
                    labels = c("Male"="Männchen","Female"="Weibchen")) +
  scale_x_continuous(breaks = seq(0, 23, 3),
                     labels = sprintf("%02d:00", seq(0, 23, 3))) +
  labs(title    = "Peak-Aktivität Uhrzeit (Nacht)",
       subtitle = "Blau = typische Nachtstunden",
       x = "Uhrzeit", y = "Anzahl Tiere") +
  theme_desk +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ── Plot D6: Nachtaktiv-Status & Saison ───────────────────────
sa_saison <- sa[!is.na(saison) & !is.na(nachtaktiv)]
sa_saison[, nachtaktiv_label := fcase(
  nachtaktiv == "ja",    "Nachtaktiv",
  nachtaktiv == "nein",  "Tag- / unklar",
  default                = "Unklar"
)]

p_d6 <- ggplot(sa_saison,
               aes(x = saison, fill = nachtaktiv_label)) +
  geom_bar(position = "fill", alpha = 0.85, width = 0.6) +
  geom_text(stat = "count",
            aes(label = after_stat(count)),
            position = position_fill(vjust = 0.5),
            color = "white", size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Nachtaktiv"="#2D6A2D",
                                "Tag- / unklar"="#C0392B",
                                "Unklar"="#888888"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(title    = "Nachtaktivität nach Saison",
       subtitle = "Anteil nachtaktiver Tiere",
       x = NULL, y = "Anteil") +
  theme_desk

# ── Zusammenführen & Speichern ────────────────────────────────
p_desk <- (p_d1 | p_d2) / (p_d3 | p_d4) / (p_d5 | p_d6) +
  plot_annotation(
    title    = "Block 0 — Deskriptive Übersicht der VHF-Aktivitätsdaten",
    subtitle = paste0("N = ", nrow(summary_alle), " Tiere verarbeitet  |  ",
                      n_nacht_j, " nachtaktiv  |  Stand: ",
                      format(Sys.Date(), "%d.%m.%Y")),
    theme = theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(output_ordner, "deskriptiv_uebersicht.png"),
       p_desk, width = 16, height = 18, dpi = 150)
cat("✓ Deskriptive Übersicht gespeichert: deskriptiv_uebersicht.png\n")

# ── Zusatzplot: Tracking-Gantt (Timeline aller Tiere) ─────────
gantt_dt <- sa[!is.na(hard_release) & !is.na(n_tage)]
gantt_dt[, t_ende := hard_release + n_tage]
gantt_dt[, igel_f := factor(igel,
  levels = gantt_dt[order(hard_release), igel])]

p_gantt <- ggplot(gantt_dt,
                  aes(y = igel_f, color = saison)) +
  geom_segment(aes(x = hard_release, xend = t_ende, yend = igel_f),
               linewidth = 5, alpha = 0.75) +
  geom_point(aes(x = hard_release), shape = 21, fill = "white",
             size = 3, stroke = 1.2) +
  scale_color_manual(values = pal_saison, name = "Saison") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  labs(
    title    = "Beobachtungszeiträume aller Igel",
    subtitle = "Punkt = Auswilderung | Balken = aktive Beobachtungsperiode | Farbe = Saison",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        axis.text.x   = element_text(angle = 30, hjust = 1),
        axis.text.y   = element_text(size = 9),
        panel.grid.major.y = element_blank())

ggsave(file.path(output_ordner, "deskriptiv_gantt.png"),
       p_gantt, width = 14, height = 10, dpi = 150)
cat("✓ Gantt-Zeitstrahl gespeichert: deskriptiv_gantt.png\n")

cat("\n✓ Fertig! Alle Dateien im Ordner:", normalizePath(output_ordner), "\n")

# ──────────────────────────────────────────────────────────────
# DATEN FÜR GAMM-ANALYSE SPEICHERN
# ──────────────────────────────────────────────────────────────

cat("\nSpeichere GAMM-Datensätze...\n")

gamm_tage_alle  <- rbindlist(gamm_tage_sum,  fill = TRUE)
gamm_nacht_alle <- rbindlist(gamm_nacht_min, fill = TRUE)

# Igel als geordneter Faktor (Igel1, Igel2, …)
igel_reihenfolge <- paste0("Igel", sort(as.integer(
  gsub("[^0-9]", "", unique(gamm_tage_alle$igel)))))
gamm_tage_alle[,  igel := factor(igel, levels = igel_reihenfolge)]
gamm_nacht_alle[, igel := factor(igel, levels = igel_reihenfolge)]

saveRDS(gamm_tage_alle,
        file.path(output_ordner, "gamm_tagesdaten.rds"))
saveRDS(gamm_nacht_alle,
        file.path(output_ordner, "gamm_nachtminuten.rds"))

# 24h-Minutendaten für Chronobiologie (Actogramme, Cosinor, Rayleigh)
chrono_alle <- rbindlist(chrono_min, fill = TRUE)
chrono_alle[, igel := factor(igel, levels = igel_reihenfolge)]
saveRDS(chrono_alle,
        file.path(output_ordner, "chrono_minuten24h.rds"))

cat("✓ gamm_tagesdaten.rds    →", nrow(gamm_tage_alle),
    "Zeilen (Igel × Tag × Tageszeit)\n")
cat("✓ gamm_nachtminuten.rds  →", format(nrow(gamm_nacht_alle), big.mark="'"),
    "Nachtminuten\n")
cat("✓ chrono_minuten24h.rds  →", format(nrow(chrono_alle), big.mark="'"),
    "Minuten (24h, alle Igel) für Chronobiologie\n")
