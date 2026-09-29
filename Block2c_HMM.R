# ==============================================================
# Block 2c — HMM activity classification
# ==============================================================
# Project:  Hedgehog VHF telemetry, Lower Saxony
# Author:   Natalie Steiner
#
# Method:
#   2-state gamma hidden Markov model on rolling_var (rolling
#   signal variance from raw CSVs) as a data-driven alternative
#   to the wmv classification.
#
# Question:
#   Do HMM and wmv yield consistent activity patterns?
#   Which bout criterion emerges empirically from the HMM?
#
# Comparisons:
#   1. Emission distributions (rolling_var per HMM state)
#   2. HMM vs. wmv agreement rate per hedgehog
#   3. Log-survivor plot (HMM IBIs vs. wmv IBIs)
#   4. Bout-metric boxplots (HMM vs. wmv)
#   5. 24h activity profile (HMM vs. wmv)
#   6. Night-activity time course (HMM vs. wmv)

# ──────────────────────────────────────────────────────────────
# PAKETE
# ──────────────────────────────────────────────────────────────
pakete <- c("data.table", "lubridate", "ggplot2", "patchwork",
            "scales", "suncalc", "readxl", "depmixS4",
            "officer", "flextable", "openxlsx")
neu <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) {
  cat("Installiere fehlende Pakete:", paste(neu, collapse = ", "), "\n")
  install.packages(neu)
}
invisible(lapply(pakete, library, character.only = TRUE))
cat("✓ Alle Pakete geladen\n\n")

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a[1] else b

# ──────────────────────────────────────────────────────────────
# EINSTELLUNGEN — hier anpassen
# ──────────────────────────────────────────────────────────────
projekt_root  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
daten_ordner  <- file.path(projekt_root, "data", "activity")
meta_datei    <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
# nacht_rds enthält smoothed_wmv (Block0, Zeile 57) — wird in Block2c NICHT
# für den wmv-Vergleich verwendet. Referenz kommt direkt aus dt_hmm$aktiv_wmv.
# nacht_rds     <- file.path(projekt_root, "output", "Block0_Pipeline", "gamm_nachtminuten.rds")
output_ordner <- file.path(projekt_root, "output", "Block2c_HMM")

dir.create(output_ordner, showWarnings = FALSE, recursive = TRUE)

# Standort
standort_lat <- 52.39729710523643
standort_lon <-  9.216876248766871
zeitzone     <- "Europe/Berlin"

# Ausschlusskriterien (wie Block2b)
min_naechte        <- 3L
min_nacht_min      <- 60L
min_pct_aktiv      <- 1.0
schwelle_aktiv_pct <- 5.0   # Toter-Sender-Schwelle

# ── Zeitfenster: Mindest-Stichprobengröße ──────────────────────
# Analysen nur für Tage mit ≥ n_min_igel Tieren (Survivorship-Bias-Schutz).
# Basierend auf den Daten: n ≥ 5 Igel = bis Tag 27.
# Wird automatisch aus den Daten ermittelt — nicht manuell anpassen.
n_min_igel <- 5L

# Bout-Kriterium für wmv (aus Block2b, für Vergleich)
bout_kriterium_wmv <- 10L   # Minuten
# Bout-Kriterium für HMM wird automatisch aus Log-Survivor abgeleitet;
# du kannst es hier manuell überschreiben falls gewünscht:
bout_kriterium_hmm_manuell <- NA  # NA = automatisch aus Log-Survivor

# HMM-Einstellungen
hmm_n_starts  <- 5L          # Anzahl Zufallsstarts (mehr = stabiler, langsamer)
hmm_maxit     <- 500L        # Maximale EM-Iterationen
set.seed(42)

# Farben
farbe_hmm  <- "#7B2D8B"   # Lila — HMM
farbe_wmv  <- "#E07B39"   # Orange — wmv (wie Block2b)

# ──────────────────────────────────────────────────────────────
# HILFSFUNKTIONEN (identisch Block2b)
# ──────────────────────────────────────────────────────────────

parse_datum <- function(x) {
  x      <- as.character(x)
  result <- suppressWarnings(as.Date(x, format = "%d.%m.%Y"))
  na_idx <- is.na(result)
  result[na_idx] <- suppressWarnings(as.Date(x[na_idx]))
  result
}

modal_val_chr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1L]
}

berechne_ibis <- function(aktiv_vec) {
  aktiv_idx <- which(aktiv_vec == 1L)
  if (length(aktiv_idx) < 2L) return(numeric(0))
  pausen <- diff(aktiv_idx)
  pausen[pausen > 1L]
}

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

logsurv_dt <- function(aktiv_liste, stunde_col = "stunde") {
  ibis <- aktiv_liste[order(igel, datum, get(stunde_col)), {
    .(ibi = berechne_ibis(aktiv))
  }, by = .(igel, datum)]
  ibi_vals <- ibis$ibi
  if (length(ibi_vals) == 0L) return(NULL)
  ibi_sort <- sort(ibi_vals)
  survivor <- (length(ibi_sort):1L) / length(ibi_sort)
  data.table(ibi = ibi_sort, survivor = survivor)
}

nacht_shift  <- function(x) ifelse(x < 12, x + 24, x)
nacht_breaks <- c(18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
nacht_labels <- sprintf("%02d:00", nacht_breaks %% 24L)
nacht_limits <- c(18.5, 30.5)

# ──────────────────────────────────────────────────────────────
# METADATEN LADEN
# ──────────────────────────────────────────────────────────────
cat("── Metadaten laden ──\n")
meta <- as.data.table(read_excel(path.expand(meta_datei), sheet = 1))
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
# A) ROHDATEN LADEN & AUF MINUTENEBENE AGGREGIEREN
# ──────────────────────────────────────────────────────────────
cat("════════════════════════════════════════════════════════\n")
cat("A) ROHDATEN LADEN (rolling_var aus CSV)\n")
cat("════════════════════════════════════════════════════════\n\n")

alle_csvs <- list.files(path.expand(daten_ordner),
                        pattern = "classification_active_passive.*\\.csv$",
                        full.names = TRUE)
cat("  Gefundene CSVs:", length(alle_csvs), "\n\n")

finde_csv <- function(igel_name) {
  nr      <- gsub("[^0-9]", "", igel_name)
  treffer <- grep(paste0("Igel[[:space:]_]?0*", nr, "[^0-9]"),
                  alle_csvs, value = TRUE)
  if (length(treffer) == 0L) return(NA_character_)
  treffer[1L]
}
meta[, csv_pfad := sapply(igel, finde_csv)]

hmm_rohdaten_liste <- list()

for (i in seq_len(nrow(meta))) {

  igel_info      <- meta[i]
  igel_name      <- igel_info$igel
  csv_pfad       <- igel_info$csv_pfad
  hard_release_i <- igel_info$hard_release

  if (is.na(csv_pfad)) {
    cat("  ⚠", igel_name, "— keine CSV, übersprungen\n")
    next
  }

  # 1. Einlesen
  dt <- tryCatch(
    fread(csv_pfad, na.strings = c("", "NA"), showProgress = FALSE),
    error = function(e) {
      cat("  ✗", igel_name, "Lesefehler:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0L) next

  # rolling_var vorhanden?
  if (!"rolling_var" %in% names(dt)) {
    cat("  ⚠", igel_name, "— Spalte 'rolling_var' fehlt, übersprungen\n")
    next
  }

  # Spalten umbenennen
  setnames(dt,
    old = c("_time"),
    new = c("time_raw"),
    skip_absent = TRUE)

  # 2. Zeitstempel
  dt[, time_utc   := ymd_hms(time_raw, tz = "UTC", quiet = TRUE)]
  dt[, time_local := with_tz(time_utc, tzone = zeitzone)]
  dt[, time_min   := floor_date(time_local, unit = "minute")]
  dt <- dt[!is.na(time_min) & !is.na(rolling_var)]

  if (nrow(dt) == 0L) next

  # 3. Minutenaggregation — mittlere rolling_var + wmv-Klassifikation (für Vergleich)
  # Fallback-Kette: wmv → mv → loio — gilt generell für alle Igel
  wmv_hat  <- "pred_nested_loio_wmv" %in% names(dt) &&
               any(dt[["pred_nested_loio_wmv"]] %in% c("a", "p"))
  mv_hat   <- "pred_nested_loio_mv"  %in% names(dt) &&
               any(dt[["pred_nested_loio_mv"]]  %in% c("a", "p"))
  loio_hat <- "pred_nested_loio"     %in% names(dt) &&
               any(dt[["pred_nested_loio"]]     %in% c("a", "p"))

  klasse_spalte_i <- if (wmv_hat)  "pred_nested_loio_wmv" else
                     if (mv_hat)   "pred_nested_loio_mv"  else
                     if (loio_hat) "pred_nested_loio"     else NA_character_

  if (!is.na(klasse_spalte_i) && klasse_spalte_i != "pred_nested_loio_wmv") {
    cat("  [", klasse_spalte_i, " → Fallback]", igel_name,
        "— wmv leer\n", sep = "")
  }

  dt_min <- dt[, .(
    rolling_var_mean = mean(rolling_var, na.rm = TRUE),
    n_det            = .N
  ), by = time_min]

  if (!is.na(klasse_spalte_i)) {
    klasse_raw <- dt[[klasse_spalte_i]]
    dt[, klasse_ref := fcase(
      klasse_raw == "a", "aktiv",
      klasse_raw == "p", "passiv",
      default = NA_character_
    )]
    klasse_agg <- dt[!is.na(klasse_ref), .(
      klasse_modal = modal_val_chr(klasse_ref)
    ), by = time_min]
    dt_min <- merge(dt_min, klasse_agg, by = "time_min", all.x = TRUE)
    dt_min[, aktiv_wmv := as.integer(klasse_modal == "aktiv")]
  } else {
    dt_min[, klasse_modal := NA_character_]
    dt_min[, aktiv_wmv    := NA_integer_]
  }

  dt_min[, datum := as.Date(time_min, tz = zeitzone)]

  # 4. Hard-Release-Filter
  if (!is.na(hard_release_i)) {
    dt_min <- dt_min[datum >= hard_release_i]
  }
  if (nrow(dt_min) == 0L) {
    cat("  ⚠", igel_name, "— keine Daten nach Auswilderungsdatum\n")
    next
  }

  # 5. Toter-Sender-Erkennung
  if (!is.na(klasse_spalte_i)) {
    tagesdaten <- dt_min[!is.na(aktiv_wmv), .(
      pct_aktiv = mean(aktiv_wmv) * 100
    ), by = datum][order(datum)]
    aktive_tage <- tagesdaten[pct_aktiv >= schwelle_aktiv_pct, datum]
    if (length(aktive_tage) == 0L) {
      cat("  ⚠", igel_name, "— kein aktiver Tag erkennbar, übersprungen\n")
      next
    }
    letzter_aktiver_tag <- max(aktive_tage)
    dt_min <- dt_min[datum <= letzter_aktiver_tag]
  }

  # 6. Tag/Nacht via suncalc
  tage_df <- data.frame(date = unique(dt_min$datum),
                        lat  = standort_lat, lon = standort_lon)
  sun <- as.data.table(getSunlightTimes(data = tage_df, tz = zeitzone,
                                         keep = c("sunrise", "sunset")))
  dt_min <- merge(dt_min, sun[, .(date, sunrise, sunset)],
                  by.x = "datum", by.y = "date", all.x = TRUE)
  dt_min[, tageszeit := fifelse(
    time_min >= sunrise & time_min < sunset, "Tag", "Nacht"
  )]
  dt_min[, stunde    := hour(time_min) + minute(time_min) / 60]
  dt_min[, tage_seit := as.integer(datum - hard_release_i) + 1L]

  # 7. Nur Nacht behalten
  nacht_i <- dt_min[
    tageszeit == "Nacht" & !is.na(rolling_var_mean),
    .(igel       = igel_name,
      datum,
      tage_seit,
      stunde,
      time_min,
      rolling_var_mean,
      aktiv_wmv,
      quelle_wmv  = klasse_spalte_i)
  ]

  if (nrow(nacht_i) == 0L) {
    cat("  ⚠", igel_name, "— keine Nachtminuten\n")
    next
  }

  hmm_rohdaten_liste[[igel_name]] <- nacht_i
  cat("  ✓", igel_name, "—", format(nrow(nacht_i), big.mark = "'"),
      "Nachtminuten | rolling_var: Median =",
      round(median(nacht_i$rolling_var_mean, na.rm = TRUE), 2), "\n")
}

dt_hmm_raw <- rbindlist(hmm_rohdaten_liste)
cat("\n  Gesamt:", format(nrow(dt_hmm_raw), big.mark = "'"),
    "Nachtminuten | Igel:", uniqueN(dt_hmm_raw$igel), "\n\n")

# Qualitätskontrolle
qc_hmm <- dt_hmm_raw[, .(
  n_naechte = uniqueN(datum),
  n_min     = .N
), by = igel]
qc_hmm[, einschluss := n_naechte >= min_naechte & n_min >= min_nacht_min]
igel_hmm <- qc_hmm[einschluss == TRUE, igel]
dt_hmm   <- dt_hmm_raw[igel %in% igel_hmm]

cat("  Nach QC:", length(igel_hmm), "Igel eingeschlossen\n")
if (any(!qc_hmm$einschluss)) {
  cat("  Ausgeschlossen:\n")
  print(qc_hmm[einschluss == FALSE])
}
cat("\n")

# ──────────────────────────────────────────────────────────────
# DYNAMISCHER ZEITCUTOFF — N-at-risk-Schwelle
# ──────────────────────────────────────────────────────────────
# Späte Zeitpunkte (ab Tag ~28) haben zu wenig Tiere für eine
# belastbare Populationsaussage. Das HMM wird nur auf dem Fenster
# Tag 1–cutoff_tag gefittet → konsistent mit Block2_GAMM.R.
# ──────────────────────────────────────────────────────────────

n_pro_tag_hmm <- dt_hmm[, .(n_igel = uniqueN(igel)), by = tage_seit]
cutoff_tag_2c <- n_pro_tag_hmm[n_igel >= n_min_igel, max(tage_seit)]

cat("══════════════════════════════════════════════\n")
cat("DYNAMISCHER ZEITCUTOFF (N-at-risk-Schwelle)\n")
cat("══════════════════════════════════════════════\n")
cat(sprintf("  Mindest-N:   %d Igel pro Zeitpunkt\n", n_min_igel))
cat(sprintf("  Cutoff Tag:  %d\n", cutoff_tag_2c))
cat(sprintf("  Zeitraum:    Tag 1 – %d nach Auswilderung\n", cutoff_tag_2c))

# N-at-risk-Tabelle um Cutoff
cat("\n  N-at-risk um den Cutoff:\n")
print(n_pro_tag_hmm[tage_seit >= max(1, cutoff_tag_2c - 5) &
                     tage_seit <= cutoff_tag_2c + 3][order(tage_seit)])
cat("══════════════════════════════════════════════\n\n")

# HMM-Datensatz auf Zeitfenster filtern
dt_hmm   <- dt_hmm[tage_seit <= cutoff_tag_2c]
igel_hmm <- dt_hmm[, unique(igel)]
cat(sprintf("  Igel im HMM nach Cutoff: %d\n\n", length(igel_hmm)))

# ──────────────────────────────────────────────────────────────
# B) WMV FÜR VERGLEICH — direkt aus CSV-Pipeline (korrekt)
# ──────────────────────────────────────────────────────────────
# aktiv_wmv in dt_hmm stammt aus Abschnitt A: pred_nested_loio_wmv
# (Fallback: pred_nested_loio_mv) — identisch mit Block2b-Pipeline.
#
# gamm_nachtminuten.rds wird NICHT verwendet: enthält laut Block0
# (Zeile 57) pred_nested_loio_smoothed_wmv, also smoothed_wmv.
# Bout-Kriterium 10 min ist korrekt für wmv (~5 min Auflösung).
cat("── B) wmv-Referenz aus CSV-Pipeline (dt_hmm$aktiv_wmv) ──\n")

# dt_hmm ist bereits auf cutoff_tag_2c gefiltert (Abschnitt A).
dt_wmv_ref <- dt_hmm[!is.na(aktiv_wmv),
                      .(igel, datum, tage_seit, stunde, aktiv = aktiv_wmv)]

if (nrow(dt_wmv_ref) == 0L) {
  cat("  ⚠ Keine wmv-Klassifikation in den CSVs gefunden — Vergleich entfällt\n\n")
  dt_wmv_ref <- NULL
} else {
  cat("  wmv-Nachtminuten:", format(nrow(dt_wmv_ref), big.mark = "'"),
      "| Igel:", uniqueN(dt_wmv_ref$igel), "\n\n")
}

igel_beide <- if (!is.null(dt_wmv_ref))
  intersect(igel_hmm, unique(dt_wmv_ref$igel)) else igel_hmm

cat(sprintf("  In HMM und wmv (Fenster Tag 1–%d): %d Igel\n\n",
            cutoff_tag_2c, length(igel_beide)))

# ──────────────────────────────────────────────────────────────
# C) HMM FITTEN — 2-Zustands-Gamma-Modell
# ──────────────────────────────────────────────────────────────
cat("════════════════════════════════════════════════════════\n")
cat("C) HMM FITTEN (2 Zustände, Gamma, pooled)\n")
cat("════════════════════════════════════════════════════════\n\n")
cat("  Hinweis: Je nach Datenmenge kann das Fitting einige\n")
cat("  Minuten dauern. Bitte warten...\n\n")

# Epsilon hinzufügen: Gamma erfordert strikt positive Werte
# Minima ~ 1e-3 ist sicher unter dem realen Wertebereich
dt_hmm[, rv_pos := pmax(rolling_var_mean, 1e-3)]

# Daten MÜSSEN exakt sortiert sein — ntimes muss mit Reihenfolge übereinstimmen
setorder(dt_hmm, igel, datum, stunde)

# Sequenzlängen: jede Nacht pro Igel = eine unabhängige Sequenz
# (verhindert Zustandsübergänge über Tageslücken hinweg)
seq_info <- dt_hmm[, .N, by = .(igel, datum)]
setorder(seq_info, igel, datum)
ntimes_vec <- seq_info$N

stopifnot(sum(ntimes_vec) == nrow(dt_hmm))  # Sanity-Check

# Startparameter aus Daten schätzen
# Naive Trennung: unteres Quartil = passiv, oberes Quartil = aktiv
mu_passiv_start <- median(dt_hmm[rv_pos < quantile(rv_pos, 0.3), rv_pos])
mu_aktiv_start  <- median(dt_hmm[rv_pos > quantile(rv_pos, 0.7), rv_pos])
cat(sprintf("  Startparameter: mu_passiv = %.3f | mu_aktiv = %.2f\n\n",
            mu_passiv_start, mu_aktiv_start))

# Modell-Vorlage
mod_template <- depmix(
  response = list(rv_pos ~ 1),
  data     = as.data.frame(dt_hmm),
  nstates  = 2L,
  family   = list(Gamma(link = "log")),
  ntimes   = ntimes_vec
)

# Mehrere Zufallsstarts → bestes Modell nach Log-Likelihood wählen
# (verhindert lokale Optima im EM-Algorithmus)
cat("  Starte", hmm_n_starts, "EM-Läufe (random starts)...\n")
fit_liste  <- vector("list", hmm_n_starts)
ll_werte   <- rep(-Inf, hmm_n_starts)

for (k in seq_len(hmm_n_starts)) {
  cat(sprintf("    Lauf %d/%d ... ", k, hmm_n_starts))
  tryCatch({
    fit_k <- fit(
      mod_template,
      verbose    = FALSE,
      emcontrol  = em.control(maxit = hmm_maxit, tol = 1e-8,
                               random.start = TRUE)
    )
    ll_werte[k]  <- logLik(fit_k)
    fit_liste[[k]] <- fit_k
    cat(sprintf("LogLik = %.1f\n", ll_werte[k]))
  }, error = function(e) cat("fehlgeschlagen\n"))
}

best_idx <- which.max(ll_werte)
if (ll_werte[best_idx] == -Inf) stop("Alle HMM-Läufe fehlgeschlagen.")

fit_hmm <- fit_liste[[best_idx]]
cat(sprintf("\n  ✓ Bestes Modell: Lauf %d | LogLik = %.1f\n\n",
            best_idx, ll_werte[best_idx]))

# Modellzusammenfassung
cat("── HMM-Parameter ──\n")
print(summary(fit_hmm))
cat("\n")

# ──────────────────────────────────────────────────────────────
# D) VITERBI-DEKODIERUNG → aktiv/passiv je Minute
# ──────────────────────────────────────────────────────────────
cat("── D) Viterbi-Dekodierung ──\n")

post <- posterior(fit_hmm)
dt_hmm[, state_raw := post$state]

# Zustand 1 = passiv (niedrige rv_pos), Zustand 2 = aktiv (hohe rv_pos)
mu_s1 <- dt_hmm[state_raw == 1L, mean(rv_pos)]
mu_s2 <- dt_hmm[state_raw == 2L, mean(rv_pos)]

if (mu_s1 > mu_s2) {
  # Zustände vertauscht — korrigieren
  dt_hmm[, state_raw := ifelse(state_raw == 1L, 2L, 1L)]
  cat("  Zustandsbezeichnung korrigiert (Zustand 1 war aktiver).\n")
}

dt_hmm[, aktiv_hmm := as.integer(state_raw == 2L)]

pct_aktiv_hmm <- dt_hmm[, mean(aktiv_hmm) * 100]
pct_aktiv_wmv_raw <- if (!is.null(dt_wmv_ref))
  dt_wmv_ref[igel %in% igel_beide, mean(aktiv) * 100] else NA

cat(sprintf("  Nachtaktivität HMM:  %.1f%%\n", pct_aktiv_hmm))
if (!is.na(pct_aktiv_wmv_raw))
  cat(sprintf("  Nachtaktivität wmv:  %.1f%%\n", pct_aktiv_wmv_raw))
cat("\n")

# ──────────────────────────────────────────────────────────────
# E) BOUT-KRITERIUM FÜR HMM (Log-Survivor)
# ──────────────────────────────────────────────────────────────
cat("── E) Log-Survivor → Bout-Kriterium HMM ──\n")

ls_hmm <- logsurv_dt(dt_hmm[igel %in% igel_beide,
                              .(igel, datum, stunde, aktiv = aktiv_hmm)])
ls_wmv_v <- if (!is.null(dt_wmv_ref)) {
  dt_wmv_v <- dt_wmv_ref[igel %in% igel_beide]
  logsurv_dt(dt_wmv_v[, .(igel, datum, stunde, aktiv)])
} else NULL

# Bout-Kriterium für HMM: manuell oder aus Log-Survivor
if (!is.na(bout_kriterium_hmm_manuell)) {
  bout_kriterium_hmm <- as.integer(bout_kriterium_hmm_manuell)
  cat("  Bout-Kriterium HMM: manuell gesetzt =", bout_kriterium_hmm, "min\n")
} else if (!is.null(ls_hmm) && nrow(ls_hmm) > 10) {
  # Automatisch: Knick via zweite Ableitung auf Log-Skala (0–60 min)
  ls_sub <- ls_hmm[ibi <= 60 & ibi >= 1]
  if (nrow(ls_sub) > 5) {
    logsurv_smooth <- loess(log(survivor) ~ ibi, data = ls_sub, span = 0.3)
    x_seq <- seq(1, 60, by = 0.5)
    y_pred <- predict(logsurv_smooth, newdata = data.frame(ibi = x_seq))
    d2     <- diff(diff(y_pred))
    knick  <- x_seq[which.max(d2) + 1L]
    bout_kriterium_hmm <- max(5L, min(30L, round(knick)))
    cat(sprintf("  Bout-Kriterium HMM: automatisch aus Log-Survivor = %d min\n",
                bout_kriterium_hmm))
    cat("  → Bitte Log-Survivor-Plot (01_log_survivor.png) prüfen\n")
    cat("    und 'bout_kriterium_hmm_manuell' oben anpassen falls nötig.\n")
  } else {
    bout_kriterium_hmm <- bout_kriterium_wmv
    cat("  ⚠ Zu wenige IBIs für Auto-Kriterium — verwende wmv-Kriterium:", bout_kriterium_hmm, "min\n")
  }
} else {
  bout_kriterium_hmm <- bout_kriterium_wmv
  cat("  ⚠ kein Log-Survivor verfügbar — verwende wmv-Kriterium:", bout_kriterium_hmm, "min\n")
}
cat("\n")

# ──────────────────────────────────────────────────────────────
# F) BOUT-METRIKEN
# ──────────────────────────────────────────────────────────────
cat("── F) Bout-Metriken berechnen ──\n")

dt_hmm_v <- dt_hmm[igel %in% igel_beide]

bm_hmm <- dt_hmm_v[order(igel, datum, stunde), {
  bm <- berechne_bouts(aktiv_hmm, bout_kriterium_hmm)
  .(n_bouts    = bm$n_bouts,
    mean_dauer = bm$mean_dauer,
    max_dauer  = bm$max_dauer,
    mean_ibi   = bm$mean_ibi,
    tage_seit  = tage_seit[1L])
}, by = .(igel, datum)]
bm_hmm[, Variante := "HMM"]

bm_wmv <- if (!is.null(dt_wmv_ref)) {
  dt_wmv_v <- dt_wmv_ref[igel %in% igel_beide]
  dt_wmv_v[order(igel, datum, stunde), {
    bm <- berechne_bouts(aktiv, bout_kriterium_wmv)
    .(n_bouts    = bm$n_bouts,
      mean_dauer = bm$mean_dauer,
      max_dauer  = bm$max_dauer,
      mean_ibi   = bm$mean_ibi,
      tage_seit  = tage_seit[1L])
  }, by = .(igel, datum)][, Variante := "wmv"]
} else NULL

cat("  Bout-Metriken HMM:  ", nrow(bm_hmm), "Nächte\n")
if (!is.null(bm_wmv)) cat("  Bout-Metriken wmv:  ", nrow(bm_wmv), "Nächte\n")
cat("\n")

# ──────────────────────────────────────────────────────────────
# G) ÜBEREINSTIMMUNGSANALYSE HMM vs. wmv
# ──────────────────────────────────────────────────────────────
cat("── G) Übereinstimmungsanalyse HMM vs. wmv ──\n")

# Merge auf Minutenebene: HMM-Klassifikation vs. wmv aus den CSVs
vergl <- dt_hmm[!is.na(aktiv_wmv) & igel %in% igel_beide,
                .(igel, datum, stunde, aktiv_hmm, aktiv_wmv)]

uebereinstimmung <- vergl[, .(
  n_min       = .N,
  pct_agree   = mean(aktiv_hmm == aktiv_wmv, na.rm = TRUE) * 100,
  pct_hmm_akt = mean(aktiv_hmm) * 100,
  pct_wmv_akt = mean(aktiv_wmv, na.rm = TRUE) * 100
), by = igel]

cat("  Mittlere Übereinstimmung HMM vs. wmv: ",
    round(mean(uebereinstimmung$pct_agree), 1), "%\n\n")
print(uebereinstimmung[order(pct_agree)][,
  .(igel, pct_agree = round(pct_agree, 1),
    pct_hmm_akt = round(pct_hmm_akt, 1),
    pct_wmv_akt = round(pct_wmv_akt, 1))])
cat("\n")

# ──────────────────────────────────────────────────────────────
# H) PLOTS
# ──────────────────────────────────────────────────────────────
cat("════════════════════════════════════════════════════════\n")
cat("H) PLOTS\n")
cat("════════════════════════════════════════════════════════\n\n")

# Hilfsfunktion: Plottheme
theme_igel <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(size = 12, face = "bold", color = "#1F4E79"),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      panel.grid.minor = element_blank()
    )
}

# ── Plot 1: Emissionsverteilungen (rolling_var je HMM-Zustand) ──
cat("  Plot 1: Emissionsverteilungen...\n")

dt_hmm[, Zustand := ifelse(aktiv_hmm == 1L, "aktiv (HMM)", "passiv (HMM)")]

# Log-Skala für bessere Lesbarkeit (rv_pos ist stark rechtsskew)
p_emission <- ggplot(dt_hmm[rv_pos < quantile(rv_pos, 0.99)],
                     aes(x = log1p(rv_pos), fill = Zustand)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80,
                 alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("aktiv (HMM)" = farbe_hmm,
                                "passiv (HMM)" = "#AAAAAA")) +
  labs(
    title    = "HMM emission distributions",
    subtitle = "log(1 + rolling_var) per decoded state  |  validates state separation",
    x        = "log(1 + rolling_var)",
    y        = "Density",
    fill     = NULL
  ) +
  theme_igel() +
  theme(legend.position = "top")

ggsave(file.path(output_ordner, "01_emissionsverteilungen.png"),
       p_emission, width = 10, height = 5, dpi = 200)
cat("  ✓ 01_emissionsverteilungen.png\n")

# ── Plot 2: Übereinstimmungsrate je Igel ──
cat("  Plot 2: Übereinstimmungsrate...\n")

p_agree <- ggplot(uebereinstimmung,
                  aes(x = reorder(igel, pct_agree), y = pct_agree)) +
  geom_col(aes(fill = pct_agree), width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = mean(uebereinstimmung$pct_agree),
             linetype = "dashed", color = "grey40") +
  scale_fill_gradient(low = "#F4A261", high = "#264653") +
  scale_y_continuous(limits = c(0, 100), labels = label_percent(scale = 1)) +
  coord_flip() +
  labs(
    title    = "Minute-level agreement: HMM vs. wmv",
    subtitle = paste0("Dashed = mean (", round(mean(uebereinstimmung$pct_agree), 1), "%)  ",
                      "| n = ", nrow(vergl[!is.na(aktiv_wmv)]), " min total"),
    x        = NULL,
    y        = "Agreement (%)"
  ) +
  theme_igel()

ggsave(file.path(output_ordner, "02_uebereinstimmung.png"),
       p_agree, width = 8, height = 6, dpi = 200)
cat("  ✓ 02_uebereinstimmung.png\n")

# ── Plot 3: Log-Survivor overlaid (HMM vs. wmv) ──
cat("  Plot 3: Log-Survivor-Plot...\n")

if (!is.null(ls_hmm) && !is.null(ls_wmv_v)) {
  lbl_hmm <- paste0("HMM  (Kriterium ", bout_kriterium_hmm, " min)")
  lbl_wmv <- paste0("wmv  (Kriterium ", bout_kriterium_wmv, " min)")

  ls_hmm[, Variante := lbl_hmm]
  ls_wmv_v[, Variante := lbl_wmv]
  ls_all <- rbindlist(list(ls_hmm, ls_wmv_v))

  farb_ls <- setNames(c(farbe_hmm, farbe_wmv), c(lbl_hmm, lbl_wmv))

  p_ls <- ggplot(ls_all[ibi <= 60], aes(x = ibi, y = log(survivor),
                                         color = Variante)) +
    geom_line(linewidth = 1.1) +
    geom_vline(xintercept = bout_kriterium_hmm,
               linetype = "dashed", color = farbe_hmm, alpha = 0.8) +
    geom_vline(xintercept = bout_kriterium_wmv,
               linetype = "dashed", color = farbe_wmv, alpha = 0.8) +
    scale_color_manual(values = farb_ls) +
    labs(
      title    = "Log-survivor plot of inter-bout intervals",
      subtitle = paste0("HMM n = ", nrow(ls_hmm), " IBIs  |  ",
                        "wmv n = ", nrow(ls_wmv_v), " IBIs  |  ",
                        "n = ", length(igel_beide), " Igel"),
      x        = "Pause between active minutes (min)",
      y        = "log(survival probability)",
      color    = NULL
    ) +
    theme_igel() +
    theme(legend.position = "top")

  ggsave(file.path(output_ordner, "03_log_survivor.png"),
         p_ls, width = 10, height = 6, dpi = 200)
  cat("  ✓ 03_log_survivor.png\n")
} else {
  cat("  ⚠ Log-Survivor nicht verfügbar (zu wenige IBIs)\n")
}

# ── Plot 4: Bout-Metriken Boxplots ──
cat("  Plot 4: Bout-Metriken...\n")

if (!is.null(bm_wmv)) {
  bm_all <- rbindlist(list(bm_hmm, bm_wmv), fill = TRUE)
  bm_all[, Variante := factor(Variante, levels = c("HMM", "wmv"))]

  farb_bm <- c(HMM = farbe_hmm, wmv = farbe_wmv)

  make_box <- function(yvar, titel, ylab) {
    ggplot(bm_all[!is.na(get(yvar))],
           aes(x = Variante, y = get(yvar), fill = Variante)) +
      geom_boxplot(outlier.shape = 21, outlier.size = 1.5,
                   outlier.alpha = 0.4, width = 0.5) +
      stat_summary(fun = median, geom = "text",
                   aes(label = round(after_stat(y), 1)),
                   vjust = -0.6, size = 3.5, color = "grey20") +
      scale_fill_manual(values = farb_bm, guide = "none") +
      labs(title = titel, x = NULL, y = ylab) +
      theme_igel()
  }

  p_bm <- (make_box("n_bouts",    "Bouts per night",    "n") +
           make_box("mean_dauer", "Mean bout duration", "min") +
           make_box("max_dauer",  "Max bout duration",  "min") +
           make_box("mean_ibi",   "Mean IBI",           "min")) +
    plot_layout(nrow = 1) +
    plot_annotation(
      title    = "Bout metrics: HMM vs. wmv",
      subtitle = paste0("Bout criterion: HMM = ", bout_kriterium_hmm,
                        " min | wmv = ", bout_kriterium_wmv,
                        " min | n = ", length(igel_beide), " Igel"),
      theme = theme(plot.title    = element_text(size = 12, face = "bold",
                                                  color = "#1F4E79"),
                    plot.subtitle = element_text(size = 10, color = "grey40"))
    )

  ggsave(file.path(output_ordner, "04_bout_metriken.png"),
         p_bm, width = 14, height = 6, dpi = 200)
  cat("  ✓ 04_bout_metriken.png\n")
}

# ── Plot 5: 24h-Aktivitätsprofil ──
cat("  Plot 5: 24h-Aktivitätsprofil...\n")

profil_hmm <- dt_hmm_v[, .(
  pct_aktiv  = mean(aktiv_hmm) * 100,
  Variante   = "HMM"
), by = .(stunde_plot = nacht_shift(floor(stunde)))]

if (!is.null(dt_wmv_ref)) {
  profil_wmv <- dt_wmv_ref[igel %in% igel_beide, .(
    pct_aktiv = mean(aktiv) * 100,
    Variante  = "wmv"
  ), by = .(stunde_plot = nacht_shift(floor(stunde)))]
  profil_all <- rbindlist(list(profil_hmm, profil_wmv))
} else {
  profil_all <- profil_hmm
}

farb_prof <- c(HMM = farbe_hmm, wmv = farbe_wmv)

p_profil <- ggplot(profil_all,
                   aes(x = stunde_plot, y = pct_aktiv,
                       color = Variante, group = Variante)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(breaks = nacht_breaks, labels = nacht_labels,
                     limits = nacht_limits) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, NA)) +
  scale_color_manual(values = farb_prof) +
  labs(
    title    = "Nocturnal activity profile: HMM vs. wmv",
    subtitle = paste0("Mean % active minutes per hour  |  n = ",
                      length(igel_beide), " Igel"),
    x        = "Time of night",
    y        = "% active minutes",
    color    = NULL
  ) +
  theme_igel() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_ordner, "05_aktivitaetsprofil.png"),
       p_profil, width = 10, height = 6, dpi = 200)
cat("  ✓ 05_aktivitaetsprofil.png\n")

# ── Plot 6: Zeitverlauf ──
cat("  Plot 6: Zeitverlauf...\n")

nacht_hmm_ts <- dt_hmm_v[, .(
  pct_aktiv = mean(aktiv_hmm) * 100,
  Variante  = "HMM"
), by = .(tage_seit)]

if (!is.null(dt_wmv_ref)) {
  nacht_wmv_ts <- dt_wmv_ref[igel %in% igel_beide, .(
    pct_aktiv = mean(aktiv) * 100,
    Variante  = "wmv"
  ), by = tage_seit]
  ts_all <- rbindlist(list(nacht_hmm_ts, nacht_wmv_ts))
} else {
  ts_all <- nacht_hmm_ts
}

p_ts <- ggplot(ts_all[tage_seit >= 1],
               aes(x = tage_seit, y = pct_aktiv,
                   color = Variante, group = Variante)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "loess", span = 0.35, se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = farb_prof) +
  scale_x_continuous(limits = c(1, cutoff_tag_2c),
                     breaks = seq(0, cutoff_tag_2c, by = 7)) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(
    title    = "Time trend: nocturnal activity after release",
    subtitle = paste0("LOESS smoothing  |  n = ", length(igel_beide),
                      " Igel  |  Analysis window: Day 1–", cutoff_tag_2c,
                      " (n ≥ ", n_min_igel, " animals/day)"),
    x        = "Days since release",
    y        = "% active night minutes",
    color    = NULL
  ) +
  theme_igel() +
  theme(legend.position = "top")

ggsave(file.path(output_ordner, "06_zeitverlauf.png"),
       p_ts, width = 12, height = 6, dpi = 200)
cat("  ✓ 06_zeitverlauf.png\n\n")

# ──────────────────────────────────────────────────────────────
# I) ZUSAMMENFASSUNG
# ──────────────────────────────────────────────────────────────
cat("════════════════════════════════════════════════════════\n")
cat("ZUSAMMENFASSUNG — HMM vs. wmv\n")
cat("════════════════════════════════════════════════════════\n\n")

cat(sprintf("  Igel im HMM-Datensatz:       %d\n", length(igel_hmm)))
cat(sprintf("  Igel im Vergleich (beide):   %d\n", length(igel_beide)))
cat(sprintf("  HMM LogLik (best run):       %.1f\n", ll_werte[best_idx]))
cat(sprintf("  Bout-Kriterium HMM:          %d min\n", bout_kriterium_hmm))
cat(sprintf("  Bout-Kriterium wmv:          %d min\n", bout_kriterium_wmv))
cat(sprintf("  Ø Übereinstimmung HMM/wmv:   %.1f%%\n\n",
            mean(uebereinstimmung$pct_agree)))

cat("── Bout-Metriken Median ──\n")
cat(sprintf("  %-30s  %6s  %6s\n", "", "HMM", "wmv"))
cat(sprintf("  %-30s  %6.1f  %6s\n", "Bouts pro Nacht:",
  median(bm_hmm$n_bouts, na.rm = TRUE),
  if (!is.null(bm_wmv)) round(median(bm_wmv$n_bouts, na.rm = TRUE), 1) else "—"))
cat(sprintf("  %-30s  %6.1f  %6s\n", "Mittl. Boutdauer (min):",
  median(bm_hmm$mean_dauer, na.rm = TRUE),
  if (!is.null(bm_wmv)) round(median(bm_wmv$mean_dauer, na.rm = TRUE), 1) else "—"))
cat(sprintf("  %-30s  %6.1f  %6s\n", "Mittl. IBI (min):",
  median(bm_hmm$mean_ibi, na.rm = TRUE),
  if (!is.null(bm_wmv)) round(median(bm_wmv$mean_ibi, na.rm = TRUE), 1) else "—"))
cat("\n")

cat("── Interpretation ──\n")
mean_agree <- mean(uebereinstimmung$pct_agree)
if (mean_agree >= 85) {
  cat("  ✓ Hohe Übereinstimmung (≥ 85%) → HMM bestätigt wmv-Klassifikation.\n")
  cat("    Die wmv-Klassifikation ist für die GAMM-Analyse geeignet.\n")
} else if (mean_agree >= 70) {
  cat("  ~ Moderate Übereinstimmung (70–85%) → HMM und wmv zeigen ähnliche,\n")
  cat("    aber nicht identische Muster. Beide Klassifikationen für die\n")
  cat("    GAMM-Analyse vergleichen (Sensitivitätsanalyse empfohlen).\n")
} else {
  cat("  ⚠ Geringe Übereinstimmung (< 70%) → HMM und wmv weichen deutlich ab.\n")
  cat("    Log-Survivor-Plot und Emissionsverteilungen prüfen.\n")
  cat("    Möglicherweise HMM-Klassifikation für GAMM verwenden.\n")
}
cat("\n")

cat("── Ausgabedateien ──\n")
cat("  Ordner:", output_ordner, "\n")
cat("   01_emissionsverteilungen.png  — rolling_var je HMM-Zustand\n")
cat("   02_uebereinstimmung.png       — % Übereinstimmung je Igel\n")
cat("   03_log_survivor.png           — Log-Survivor HMM vs. wmv\n")
cat("   04_bout_metriken.png          — Bout-Metriken Boxplots\n")
cat("   05_aktivitaetsprofil.png      — 24h-Profil HMM vs. wmv\n")
cat("   06_zeitverlauf.png            — Zeitverlauf Nachtaktivität\n")
cat("════════════════════════════════════════════════════════\n\n")

cat("→ Nächster Schritt:\n")
cat("  1. Plot 01 prüfen: sind die zwei Zustände klar getrennt?\n")
cat("  2. Plot 03 prüfen: bestätigt der Log-Survivor-Knick das Kriterium?\n")
cat("  3. Falls Übereinstimmung ≥ 85%: wmv in Block2_GAMM.R verwenden.\n")
cat("  4. Falls Übereinstimmung < 85%: HMM-Ausgabe als Sensitivitätscheck\n")
cat("     in Block2_GAMM.R einbauen.\n\n")

# HMM-dekodierte Klassifikation als RDS speichern
# (kann in Block2_GAMM.R als Alternative zu wmv geladen werden)
hmm_export <- dt_hmm[, .(igel, datum, tage_seit, stunde,
                           aktiv_hmm, rolling_var_mean)]
saveRDS(hmm_export, file.path(output_ordner, "hmm_nachtminuten.rds"))
cat("  ✓ hmm_nachtminuten.rds gespeichert (für Block2_GAMM.R als Alternative)\n\n")

# ──────────────────────────────────────────────────────────────
# J) EXCEL-ÜBERBLICK
# ──────────────────────────────────────────────────────────────
cat("════════════════════════════════════════════════════════\n")
cat("J) EXCEL & BERICHT erstellen...\n")
cat("════════════════════════════════════════════════════════\n\n")

library(openxlsx)

# ── Styles ──
stil_h_blau <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white",
                            fgFill = "#2C3E6B", halign = "center", valign = "center",
                            textDecoration = "bold", border = "TopBottomLeftRight",
                            borderColour = "#CCCCCC", wrapText = TRUE)
stil_h_lila <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white",
                            fgFill = "#7B2D8B", halign = "center", valign = "center",
                            textDecoration = "bold", border = "TopBottomLeftRight",
                            borderColour = "#CCCCCC", wrapText = TRUE)
stil_h_org  <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white",
                            fgFill = "#E07B39", halign = "center", valign = "center",
                            textDecoration = "bold", border = "TopBottomLeftRight",
                            borderColour = "#CCCCCC", wrapText = TRUE)
stil_normal <- createStyle(fontName = "Arial", fontSize = 9,
                            border = "TopBottomLeftRight", borderColour = "#CCCCCC")
stil_zahl1  <- createStyle(fontName = "Arial", fontSize = 9, numFmt = "0.0",
                            border = "TopBottomLeftRight", borderColour = "#CCCCCC")
stil_pct    <- createStyle(fontName = "Arial", fontSize = 9, numFmt = "0.0",
                            border = "TopBottomLeftRight", borderColour = "#CCCCCC")
stil_zebra  <- createStyle(fontName = "Arial", fontSize = 9, fgFill = "#F0F4F8",
                            border = "TopBottomLeftRight", borderColour = "#CCCCCC")
stil_zebra1 <- createStyle(fontName = "Arial", fontSize = 9, fgFill = "#F0F4F8",
                            numFmt = "0.0", border = "TopBottomLeftRight",
                            borderColour = "#CCCCCC")

wb <- createWorkbook()

# ── Sheet 1: Pro-Igel Übersicht ──
addWorksheet(wb, "Uebersicht je Igel")

# Quellinformation je Igel zusammenstellen
quelle_info <- as.data.frame(dt_hmm_raw[, .(
  quelle_wmv = unique(quelle_wmv)[1]
), by = igel])

# Alle Teile als base-R data.frames mergen — verhindert data.table-Namens-Bugs
df_qc <- as.data.frame(qc_hmm[, .(igel, n_naechte, n_min, einschluss)])

df_agree_ex <- if (nrow(uebereinstimmung) > 0) {
  as.data.frame(uebereinstimmung[, .(
    igel,
    pct_agree   = round(pct_agree,   1),
    pct_hmm_akt = round(pct_hmm_akt, 1),
    pct_wmv_akt = round(pct_wmv_akt, 1)
  )])
} else {
  data.frame(igel = character(0), pct_agree = numeric(0),
             pct_hmm_akt = numeric(0), pct_wmv_akt = numeric(0))
}

# Bout-Metriken je Igel
bm_hmm_igel <- as.data.frame(bm_hmm[, .(
  hmm_bouts_med  = round(median(n_bouts,    na.rm = TRUE), 1),
  hmm_dauer_med  = round(median(mean_dauer, na.rm = TRUE), 1),
  hmm_ibi_med    = round(median(mean_ibi,   na.rm = TRUE), 1)
), by = igel])

# Schrittweise mergen
df_igel <- merge(df_qc, df_agree_ex,  by = "igel", all.x = TRUE)
df_igel <- merge(df_igel, quelle_info, by = "igel", all.x = TRUE)
df_igel <- merge(df_igel, bm_hmm_igel, by = "igel", all.x = TRUE)

if (!is.null(bm_wmv)) {
  bm_wmv_igel <- as.data.frame(bm_wmv[, .(
    wmv_bouts_med = round(median(n_bouts,    na.rm = TRUE), 1),
    wmv_dauer_med = round(median(mean_dauer, na.rm = TRUE), 1),
    wmv_ibi_med   = round(median(mean_ibi,   na.rm = TRUE), 1)
  ), by = igel])
  df_igel <- merge(df_igel, bm_wmv_igel, by = "igel", all.x = TRUE)
}

df_igel <- df_igel[order(df_igel$igel), ]

# Spaltennamen setzen — mit expliziter Längenprüfung
neue_namen <- c("Igel", "Nächte", "Minuten (Nacht)", "Einschluss HMM",
                "Übereinstimmung HMM/wmv (%)", "% aktiv (HMM)", "% aktiv (wmv-Ref)",
                "wmv/mv-Quelle (Vergleich)",
                "HMM Bouts/Nacht", "HMM Dauer (min)", "HMM IBI (min)")
if (!is.null(bm_wmv))
  neue_namen <- c(neue_namen, "wmv Bouts/Nacht", "wmv Dauer (min)", "wmv IBI (min)")

if (length(neue_namen) == ncol(df_igel)) {
  names(df_igel) <- neue_namen
} else {
  cat(sprintf("  ⚠ Spaltenanzahl: df_igel hat %d, Namen haben %d — Original behalten\n",
              ncol(df_igel), length(neue_namen)))
  cat("  Aktuelle Spalten:", paste(names(df_igel), collapse = ", "), "\n")
}

nc1 <- ncol(df_igel)
writeData(wb, 1, df_igel, startRow = 1, startCol = 1, colNames = TRUE)
addStyle(wb, 1, stil_h_blau, rows = 1, cols = 1:min(4, nc1), gridExpand = TRUE)
if (nc1 > 4)
  addStyle(wb, 1, stil_h_lila, rows = 1, cols = 5:nc1, gridExpand = TRUE)
for (i in seq_len(nrow(df_igel))) {
  st  <- if (i %% 2 == 0) stil_zebra  else stil_normal
  st1 <- if (i %% 2 == 0) stil_zebra1 else stil_zahl1
  addStyle(wb, 1, st,  rows = i + 1, cols = 1:min(4, nc1), gridExpand = TRUE)
  if (nc1 > 4)
    addStyle(wb, 1, st1, rows = i + 1, cols = 5:nc1, gridExpand = TRUE)
}
col_widths <- c(10, 10, 16, 16, 28, 16, 18, 24, rep(16, max(nc1 - 8, 0)))
setColWidths(wb, 1, cols = seq_len(nc1), widths = col_widths[seq_len(nc1)])
setRowHeights(wb, 1, rows = 1, heights = 42)

# ── Sheet 2: HMM-Modellparameter ──
addWorksheet(wb, "HMM-Parameter")

# Transition matrix
trans_mat <- getpars(fit_hmm)
# Extract emission parameters from summary
hmm_sum <- capture.output(summary(fit_hmm))

# Emissionsparameter: beobachtete Mittelwerte je dekodiertem Zustand
# (zuverlässiger als coef() auf depmixS4-Objekten, biologisch interpretierbar)
mu_s1_fit <- round(dt_hmm[state_raw == 1L, mean(rv_pos, na.rm = TRUE)], 4)
mu_s2_fit <- round(dt_hmm[state_raw == 2L, mean(rv_pos, na.rm = TRUE)], 2)
cv_s1_fit <- round(dt_hmm[state_raw == 1L, sd(rv_pos, na.rm = TRUE)] /
                   max(mu_s1_fit, 1e-6), 3)
cv_s2_fit <- round(dt_hmm[state_raw == 2L, sd(rv_pos, na.rm = TRUE)] /
                   max(mu_s2_fit, 1e-6), 3)

cat(sprintf("  Zustand 1 (passiv): µ = %.4f | CV = %.3f\n", mu_s1_fit, cv_s1_fit))
cat(sprintf("  Zustand 2 (aktiv):  µ = %.2f  | CV = %.3f\n", mu_s2_fit, cv_s2_fit))

# Verweildauer: erst über trprob(), dann empirisch aus Viterbi-Sequenz
pmat <- tryCatch(
  trprob(fit_hmm),
  error = function(e) matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)
)
verweildauer_passiv <- round(1 / max(1 - pmat[1, 1], 0.01), 1)
verweildauer_aktiv  <- round(1 / max(1 - pmat[2, 2], 0.01), 1)

# Empirische Verweildauer als Fallback (aus RLE der Viterbi-Zustände)
dwell_raw <- dt_hmm[order(igel, datum, stunde), {
  r <- rle(state_raw)
  .(dauer = r$lengths, state = r$values)
}, by = .(igel, datum)]
dwell_passiv_emp <- round(mean(dwell_raw[state == 1L, dauer], na.rm = TRUE), 1)
dwell_aktiv_emp  <- round(mean(dwell_raw[state == 2L, dauer], na.rm = TRUE), 1)

# Falls trprob-Fallback (beide = 10) → empirische Werte verwenden
if (verweildauer_passiv == 10 && verweildauer_aktiv == 10) {
  verweildauer_passiv <- dwell_passiv_emp
  verweildauer_aktiv  <- dwell_aktiv_emp
  cat("  Verweildauer: empirisch aus Viterbi-Sequenz\n")
} else {
  cat(sprintf("  Verweildauer (Modell): passiv %.1f min | aktiv %.1f min\n",
              verweildauer_passiv, verweildauer_aktiv))
}

df_params <- data.frame(
  Parameter = c(
    "Anzahl Zustände",
    "Verteilungsfamilie",
    "Observationsvariable",
    "Optimierungsverfahren",
    "Anzahl Zufallsstarts",
    "Maximale EM-Iterationen",
    "Log-Likelihood (bester Lauf)",
    "",
    "Zustand 1 (passiv) — mittl. rolling_var",
    "Zustand 1 (passiv) — Variationskoeffizient",
    "Zustand 2 (aktiv)  — mittl. rolling_var",
    "Zustand 2 (aktiv)  — Variationskoeffizient",
    "",
    "Mittl. Verweildauer passiv (min)",
    "Mittl. Verweildauer aktiv (min)",
    "",
    "Bout-Kriterium HMM (Log-Survivor)",
    "Bout-Kriterium wmv (Block2b)"
  ),
  Wert = c(
    "2",
    "Gamma (log-link)",
    "rolling_var (mittl. Signalvarianz je Minute)",
    "EM-Algorithmus (depmixS4)",
    as.character(hmm_n_starts),
    as.character(hmm_maxit),
    round(ll_werte[best_idx], 1),
    "",
    round(mu_s1_fit, 3),
    round(cv_s1_fit, 3),
    round(mu_s2_fit, 2),
    round(cv_s2_fit, 3),
    "",
    verweildauer_passiv,
    verweildauer_aktiv,
    "",
    paste0(bout_kriterium_hmm, " min"),
    paste0(bout_kriterium_wmv, " min")
  ),
  stringsAsFactors = FALSE
)

writeData(wb, 2, df_params, startRow = 1, startCol = 1, colNames = TRUE)
addStyle(wb, 2, stil_h_blau, rows = 1, cols = 1:2, gridExpand = TRUE)
for (i in seq_len(nrow(df_params))) {
  st <- if (i %% 2 == 0) stil_zebra else stil_normal
  addStyle(wb, 2, st, rows = i + 1, cols = 1:2, gridExpand = TRUE)
}
setColWidths(wb, 2, cols = 1:2, widths = c(46, 36))

# ── Sheet 3: Gesamtstatistik ──
addWorksheet(wb, "Gesamtstatistik")

# Fallback-Info
fallback_hmm <- dt_hmm_raw[quelle_wmv == "pred_nested_loio_mv", unique(igel)]

# Gesamtstatistik: als Vektoren aufbauen, dann data.frame — verhindert NULL-in-c()-Bugs
ges_kennzahl <- c(
  "Igel mit Rohdaten (CSV vorhanden)",
  "Igel im HMM-Datensatz (nach QC)",
  "Igel im Vergleich HMM/wmv",
  "Igel mit mv-Fallback (wmv leer)",
  "",
  "Nachtminuten gesamt (HMM-Datensatz)",
  "Davon aktiv (HMM)",
  "Davon aktiv (wmv-Referenz)",
  "",
  "Ø Übereinstimmung HMM vs. wmv (%)",
  "Min. Übereinstimmung (schlechtester Igel)",
  "Max. Übereinstimmung (bester Igel)",
  "",
  "Median Bouts/Nacht — HMM",
  "Median Boutdauer (min) — HMM",
  "Median IBI (min) — HMM"
)
ges_wert <- c(
  as.character(length(hmm_rohdaten_liste)),
  as.character(length(igel_hmm)),
  as.character(length(igel_beide)),
  if (length(fallback_hmm) > 0) paste(sort(fallback_hmm), collapse = ", ") else "keine",
  "",
  format(nrow(dt_hmm), big.mark = "'"),
  paste0(round(mean(dt_hmm$aktiv_hmm) * 100, 1), " %"),
  if (!is.null(dt_wmv_ref))
    paste0(round(dt_wmv_ref[igel %in% igel_beide, mean(aktiv)] * 100, 1), " %")
  else "—",
  "",
  as.character(round(mean(uebereinstimmung$pct_agree), 1)),
  as.character(round(min(uebereinstimmung$pct_agree),  1)),
  as.character(round(max(uebereinstimmung$pct_agree),  1)),
  "",
  as.character(round(median(bm_hmm$n_bouts,    na.rm = TRUE), 1)),
  as.character(round(median(bm_hmm$mean_dauer, na.rm = TRUE), 1)),
  as.character(round(median(bm_hmm$mean_ibi,   na.rm = TRUE), 1))
)
if (!is.null(bm_wmv)) {
  ges_kennzahl <- c(ges_kennzahl,
    "Median Bouts/Nacht — wmv",
    "Median Boutdauer (min) — wmv",
    "Median IBI (min) — wmv")
  ges_wert <- c(ges_wert,
    as.character(round(median(bm_wmv$n_bouts,    na.rm = TRUE), 1)),
    as.character(round(median(bm_wmv$mean_dauer, na.rm = TRUE), 1)),
    as.character(round(median(bm_wmv$mean_ibi,   na.rm = TRUE), 1)))
}
df_ges <- data.frame(Kennzahl = ges_kennzahl, Wert = ges_wert,
                     stringsAsFactors = FALSE)

writeData(wb, 3, df_ges, startRow = 1, startCol = 1, colNames = TRUE)
addStyle(wb, 3, stil_h_blau, rows = 1, cols = 1:2, gridExpand = TRUE)
for (i in seq_len(nrow(df_ges))) {
  st <- if (i %% 2 == 0) stil_zebra else stil_normal
  addStyle(wb, 3, st, rows = i + 1, cols = 1:2, gridExpand = TRUE)
}
setColWidths(wb, 3, cols = 1:2, widths = c(44, 30))

excel_pfad <- file.path(output_ordner, "Block2c_HMM_Uebersicht.xlsx")
saveWorkbook(wb, excel_pfad, overwrite = TRUE)
cat("  ✓ Excel gespeichert:", excel_pfad, "\n\n")

# ──────────────────────────────────────────────────────────────
# K) WORD-METHODENBERICHT
# ──────────────────────────────────────────────────────────────
library(officer)
library(flextable)

add_plot_safe <- function(doc, pfad, breite = 16, hoehe = 10) {
  if (file.exists(pfad)) {
    doc <- body_add_img(doc, pfad,
                        width  = breite / 2.54,
                        height = hoehe  / 2.54)
  } else {
    doc <- body_add_par(doc,
      paste0("[Abbildung nicht gefunden: ", basename(pfad), "]"),
      style = "Normal")
  }
  doc
}

add_section <- function(doc, titel, level = 1) {
  body_add_par(doc, titel,
               style = if (level == 1) "heading 1" else "heading 2")
}

make_ft <- function(df, header_farbe = "#2C3E6B", zebra = TRUE) {
  nr <- nrow(as.data.frame(df))
  ft <- flextable(as.data.frame(df)) |>
    bold(part = "header") |>
    bg(part = "header", bg = header_farbe) |>
    color(part = "header", color = "white") |>
    border_outer(part = "all",
                 border = fp_border(color = "#CCCCCC", width = 1)) |>
    border_inner_h(part = "body",
                   border = fp_border(color = "#E8E8E8", width = 0.5)) |>
    fontsize(size = 9, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    autofit()
  if (zebra && nr >= 2)
    ft <- bg(ft, i = seq(2, nr, 2), bg = "#F0F4F8")
  ft
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
  body_add_par("Block 2c — HMM-Aktivitätsklassifikation",
               style = "heading 2") |>
  body_add_par("VHF-Igelbesenderung — Europäischer Igel (Erinaceus europaeus)",
               style = "Normal")

df_meta <- data.frame(
  Feld = c("Autorin", "Institution", "Standort", "Analysedatum",
           "Methode", "Bout-Kriterium HMM", "Bout-Kriterium wmv (Referenz)"),
  Wert = c("Natalie Steiner",
           "Stiftung Tierärztliche Hochschule Hannover",
           "Wildtierstation Sachsenhagen, Niedersachsen",
           format(Sys.Date(), "%d. %B %Y"),
           "2-Zustands-Gamma-HMM auf rolling_var (depmixS4)",
           paste0(bout_kriterium_hmm, " Minuten (Log-Survivor)"),
           paste0(bout_kriterium_wmv, " Minuten (Block 2b)")),
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
  width(j = 1, width = 6   / 2.54) |>
  width(j = 2, width = 11.5 / 2.54)
doc <- body_add_flextable(doc, ft_meta)
doc <- body_add_break(doc)

# ── 1. FRAGESTELLUNG ──
doc <- add_section(doc, "1  Fragestellung")

doc <- body_add_par(doc, paste0(
  "Ziel dieses Analysebausteins ist die Validierung der wmv-Klassifikation ",
  "(pred_nested_loio_wmv; 300-s-Rollfenster, ~5 min effektive Auflösung) ",
  "durch eine vollständig datengetriebene Alternative: ein ",
  "2-Zustands-Hidden-Markov-Modell (HMM) auf der kontinuierlichen Signalvarianz ",
  "(rolling_var) aus den Roh-CSV-Dateien. ",
  "Das HMM klassifiziert jede Minute als 'aktiv' oder 'passiv', ohne ",
  "vorab ein Zeitfenster oder einen Schwellenwert vorzugeben. ",
  "Eine hohe Übereinstimmung zwischen HMM und wmv bestätigt, dass die ",
  "wmv-Klassifikation das zugrundeliegende Aktivitätssignal korrekt abbildet."),
  style = "Normal")

# ── 2. METHODIK ──
doc <- add_section(doc, "2  Methodik")

doc <- add_section(doc, "2.1  Datenbasis und Vorverarbeitung", level = 2)
doc <- body_add_par(doc, paste0(
  "Grundlage sind die Roh-CSV-Dateien aus dem tRackIT-System ",
  "(je eine Datei pro Igel, Namensschema: ",
  "classification_active_passive_IgelX_300s_window_plus_filter.csv). ",
  "Die Spalte rolling_var enthält die rollende Varianz des VHF-Rohsignals ",
  "über ein 300-s-Fenster (berechnete Detektionsebene, ~1.2 s Intervall). ",
  "Ein aktiver Igel erzeugt durch Körperbewegungen Signalschwankungen ",
  "(rolling_var >> 0), ein ruhender Igel ein stabiles Signal (rolling_var ≈ 0). ",
  "Für die HMM-Analyse wurde rolling_var auf Minutenebene gemittelt. ",
  "Anschliessend wurden dieselben Filter wie in Block 0 und Block 2b angewendet: ",
  "Hard-Release-Filter (nur Daten ab Auswilderungsdatum), ",
  "Toter-Sender-Erkennung (Ausschluss nach letztem aktivem Tag), ",
  "Nacht-Filter via suncalc (Sonnenuntergang bis Sonnenaufgang)."),
  style = "Normal")

doc <- body_add_par(doc, paste0(
  "Für Igel 1–6 war die Spalte pred_nested_loio_wmv im tRackIT-Export leer ",
  "(älteres Exportformat). Als Vergleichsreferenz wurde automatisch auf ",
  "pred_nested_loio_mv zurückgegriffen. ",
  "Der HMM selbst ist davon unabhängig: rolling_var ist für alle Tiere vollständig vorhanden. ",
  if (length(fallback_hmm) > 0)
    paste0("Folgende Igel nutzen den mv-Fallback: ",
           paste(sort(fallback_hmm), collapse = ", "), ".")
  else "Kein Fallback notwendig."),
  style = "Normal")

doc <- add_section(doc, "2.2  Hidden Markov Modell", level = 2)
doc <- body_add_par(doc, paste0(
  "Ein Hidden Markov Modell (HMM) modelliert eine beobachtbare Zeitreihe ",
  "(hier: mittlere rolling_var je Minute) als Funktion eines verborgenen ",
  "Zustands (hier: 'aktiv' oder 'passiv'). Die Übergangswahrscheinlichkeiten ",
  "zwischen den Zuständen bestimmen die typische Verweildauer und damit ",
  "implizit die Bout-Zeitskala — ohne dass ein Fenster oder Schwellenwert ",
  "vorab festgelegt werden muss."),
  style = "Normal")

doc <- body_add_par(doc, paste0(
  "Modellspezifikation: 2 verborgene Zustände (aktiv / passiv), ",
  "Gamma-Verteilung für die Emissionen (log-link), ",
  "jede Nacht pro Individuum als eigenständige Sequenz (keine Übergänge ",
  "über Tageslücken). Fitting via EM-Algorithmus (R-Paket depmixS4; ",
  "Visser & Speekenbrink 2010) mit ", hmm_n_starts, " Zufallsstarts und ",
  "Selektion nach maximaler Log-Likelihood. ",
  "Zustandsfolge per Viterbi-Algorithmus dekodiert."),
  style = "Normal")

# Modellparameter-Tabelle
df_params_doc <- df_params[df_params$Parameter != "", ]
doc <- body_add_flextable(doc, make_ft(df_params_doc))

# ── 3. VALIDIERUNG: EMISSIONSVERTEILUNGEN ──
doc <- add_section(doc, "3  Modellvalidierung: Emissionsverteilungen")

doc <- body_add_par(doc, paste0(
  "Die Emissionsverteilungen zeigen die Verteilung von log(1 + rolling_var) ",
  "getrennt für die beiden dekodierten Zustände. Eine gute Trennung der ",
  "Verteilungen bestätigt, dass das HMM zwei biologisch distinkte Zustände ",
  "identifiziert hat. ",
  sprintf("Zustand 'passiv' hat eine mittlere rolling_var von %.3f, ",
          mu_s1_fit),
  sprintf("Zustand 'aktiv' von %.2f — Faktor %.0fx höher.",
          mu_s2_fit, mu_s2_fit / max(mu_s1_fit, 0.001))),
  style = "Normal")

doc <- add_plot_safe(doc,
  file.path(output_ordner, "01_emissionsverteilungen.png"),
  breite = 16, hoehe = 8)

# ── 4. ÜBEREINSTIMMUNG HMM vs. WMV ──
doc <- add_section(doc, "4  Übereinstimmung HMM vs. wmv")

doc <- body_add_par(doc, paste0(
  "Die minutenweise Übereinstimmung zwischen HMM-Klassifikation und ",
  "wmv-Referenzklassifikation (pred_nested_loio_wmv bzw. mv-Fallback) ",
  "beträgt im Mittel ",
  round(mean(uebereinstimmung$pct_agree), 1), "% ",
  "(Spanne: ", round(min(uebereinstimmung$pct_agree), 1), "–",
  round(max(uebereinstimmung$pct_agree), 1), "%). ",
  if (mean(uebereinstimmung$pct_agree) >= 85)
    paste0("Dies entspricht einer hohen Übereinstimmung und bestätigt, ",
           "dass die wmv-Klassifikation das Aktivitätssignal korrekt erfasst.")
  else if (mean(uebereinstimmung$pct_agree) >= 70)
    paste0("Dies entspricht einer moderaten Übereinstimmung. ",
           "Eine Sensitivitätsanalyse mit beiden Klassifikationen wird empfohlen.")
  else
    paste0("Die Übereinstimmung ist gering. ",
           "Die HMM-Klassifikation sollte als Alternative zur wmv-Methode ",
           "in Betracht gezogen werden.")),
  style = "Normal")

# Übereinstimmungstabelle je Igel
df_agree_doc <- uebereinstimmung[order(-pct_agree), .(
  Igel                    = igel,
  `Übereinstimmung (%)`   = round(pct_agree,   1),
  `% aktiv (HMM)`         = round(pct_hmm_akt, 1),
  `% aktiv (wmv/mv-Ref.)` = round(pct_wmv_akt, 1),
  `Minuten gesamt`        = n_min
)]
doc <- body_add_flextable(doc, make_ft(df_agree_doc))
doc <- add_plot_safe(doc,
  file.path(output_ordner, "02_uebereinstimmung.png"),
  breite = 12, hoehe = 8)

# ── 5. LOG-SURVIVOR ──
doc <- add_section(doc, "5  Log-Survivor-Plot und Bout-Kriterium")

doc <- body_add_par(doc, paste0(
  "Der Log-Survivor-Plot der Inter-Bout-Intervalle (IBI) zeigt, ob das HMM ",
  "eine biologisch interpretierbare Bout-Struktur produziert. ",
  "Ein Knick in der Kurve markiert den Übergang zwischen kurzen ",
  "intra-bout-Pausen und längeren inter-bout-Intervallen und liefert ",
  "das empirische Bout-Kriterium. ",
  sprintf("Für die HMM-Klassifikation wurde ein Kriterium von %d Minuten ",
          bout_kriterium_hmm),
  "abgeleitet",
  if (!is.na(bout_kriterium_hmm_manuell))
    " (manuell gesetzt)."
  else
    " (automatisch via zweiter Ableitung der Log-Survivor-Kurve)."),
  style = "Normal")

if (!is.null(ls_hmm) && !is.null(ls_wmv_v)) {
  df_ls_tab <- data.frame(
    Variante = c("HMM", "wmv"),
    `IBIs gesamt` = c(nrow(ls_hmm), nrow(ls_wmv_v)),
    `Bout-Kriterium (min)` = c(bout_kriterium_hmm, bout_kriterium_wmv),
    check.names = FALSE
  )
  doc <- body_add_flextable(doc, make_ft(df_ls_tab))
}

doc <- add_plot_safe(doc,
  file.path(output_ordner, "03_log_survivor.png"),
  breite = 16, hoehe = 8)

# ── 6. BOUT-METRIKEN ──
doc <- add_section(doc, "6  Bout-Metriken")

doc <- body_add_par(doc, paste0(
  "Die Bout-Metriken (Bouts pro Nacht, Boutdauer, IBI) zeigen, wie ähnlich ",
  "die biologischen Muster zwischen HMM- und wmv-Klassifikation sind. ",
  "Stimmen die Mediane überein, spricht dies für die Robustheit der Ergebnisse ",
  "gegenüber der Klassifikationsmethode."),
  style = "Normal")

df_bm_doc <- data.frame(
  Metrik = c("Bouts pro Nacht (Median)",
             "Mittlere Boutdauer — Median (min)",
             "Maximale Boutdauer — Median (min)",
             "Mittleres IBI — Median (min)"),
  HMM = c(
    round(median(bm_hmm$n_bouts,    na.rm = TRUE), 1),
    round(median(bm_hmm$mean_dauer, na.rm = TRUE), 1),
    round(median(bm_hmm$max_dauer,  na.rm = TRUE), 1),
    round(median(bm_hmm$mean_ibi,   na.rm = TRUE), 1)
  ),
  wmv = if (!is.null(bm_wmv)) c(
    round(median(bm_wmv$n_bouts,    na.rm = TRUE), 1),
    round(median(bm_wmv$mean_dauer, na.rm = TRUE), 1),
    round(median(bm_wmv$max_dauer,  na.rm = TRUE), 1),
    round(median(bm_wmv$mean_ibi,   na.rm = TRUE), 1)
  ) else rep("—", 4),
  stringsAsFactors = FALSE
)
doc <- body_add_flextable(doc, make_ft(df_bm_doc))

if (file.exists(file.path(output_ordner, "04_bout_metriken.png")))
  doc <- add_plot_safe(doc,
    file.path(output_ordner, "04_bout_metriken.png"),
    breite = 16, hoehe = 8)

# ── 7. AKTIVITÄTSPROFIL & ZEITVERLAUF ──
doc <- add_section(doc, "7  Aktivitätsprofil und Zeitverlauf")

doc <- body_add_par(doc, paste0(
  "Das 24-h-Aktivitätsprofil und der Zeitverlauf dienen als abschliessende ",
  "Plausibilitätsprüfung: Wenn HMM und wmv dieselbe glockenförmige Nachtaktivität ",
  "und denselben zeitlichen Trend zeigen, ist die Klassifikationsmethode ",
  "für die Kernaussagen der Arbeit irrelevant — das biologische Signal ist robust."),
  style = "Normal")

doc <- add_section(doc, "24-h-Aktivitätsprofil", level = 2)
doc <- add_plot_safe(doc,
  file.path(output_ordner, "05_aktivitaetsprofil.png"),
  breite = 16, hoehe = 8)

doc <- add_section(doc, "Zeitverlauf nach Auswilderung", level = 2)
doc <- add_plot_safe(doc,
  file.path(output_ordner, "06_zeitverlauf.png"),
  breite = 16, hoehe = 8)

# ── 8. SCHLUSSFOLGERUNG ──
doc <- add_section(doc, "8  Schlussfolgerung und Empfehlung")

mean_agree_doc <- mean(uebereinstimmung$pct_agree)
empfehlung <- if (mean_agree_doc >= 85) {
  paste0(
    "Die mittlere minutenweise Übereinstimmung von ",
    round(mean_agree_doc, 1),
    "% bestätigt, dass die wmv-Klassifikation (pred_nested_loio_wmv / mv-Fallback) ",
    "das Aktivitätssignal korrekt und konsistent mit der vollständig ",
    "datengetriebenen HMM-Methode erfasst. ",
    "Die wmv-Klassifikation wird daher für die GAMM-Analyse in Block 2 empfohlen. ",
    "Das HMM liefert mit einem Bout-Kriterium von ", bout_kriterium_hmm,
    " Minuten ein unabhängig validiertes Referenzkriterium.")
} else if (mean_agree_doc >= 70) {
  paste0(
    "Die mittlere Übereinstimmung von ",
    round(mean_agree_doc, 1),
    "% ist moderat. Es wird empfohlen, die GAMM-Analyse sowohl mit der ",
    "wmv- als auch mit der HMM-Klassifikation durchzuführen (Sensitivitätsanalyse) ",
    "und die Robustheit der Hauptergebnisse zu prüfen.")
} else {
  paste0(
    "Die mittlere Übereinstimmung von ",
    round(mean_agree_doc, 1),
    "% ist gering. Die HMM-Klassifikation (gespeichert als ",
    "hmm_nachtminuten.rds) sollte als primäre Klassifikation für die ",
    "GAMM-Analyse in Betracht gezogen werden. ",
    "Die Emissionsverteilungen und den Log-Survivor-Plot sorgfältig prüfen.")
}

doc <- body_add_par(doc, empfehlung, style = "Normal")

df_empf <- data.frame(
  Kriterium = c(
    "Datengetriebene Klassifikation ohne Vorannahmen",
    "Übereinstimmung mit wmv/mv-Referenz",
    "Log-Survivor-Knick vorhanden",
    "Mittl. Verweildauer aktiv (min)",
    "Bout-Kriterium (empirisch)",
    "Empfehlung für Block 2 GAMM"
  ),
  HMM = c(
    "Ja — rolling_var direkt modelliert",
    paste0(round(mean_agree_doc, 1), "%"),
    if (!is.null(ls_hmm) && nrow(ls_hmm) > 10) "Prüfen (Plot 03)" else "zu wenige IBIs",
    paste0(verweildauer_aktiv, " min"),
    paste0(bout_kriterium_hmm, " min"),
    if (mean_agree_doc >= 85) "Validiert wmv" else "Sensitivitätsanalyse"
  ),
  `wmv/mv` = c(
    "Nein — Fenster vorgeschrieben (300 s)",
    "—",
    "Ja (Block 2b, Knick ~5–8 min)",
    "—",
    paste0(bout_kriterium_wmv, " min"),
    if (mean_agree_doc >= 85) "Empfohlen (primär)" else "Als Alternative prüfen"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ft_empf <- flextable(df_empf) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#2C3E6B") |>
  color(part = "header", color = "white") |>
  bg(i = nrow(df_empf), bg = "#E8F4E8") |>
  bg(i = seq(2, nrow(df_empf) - 1, 2), bg = "#F0F4F8") |>
  bold(i = nrow(df_empf)) |>
  border_outer(part = "all", border = fp_border(color = "#CCCCCC", width = 1)) |>
  border_inner_h(part = "body", border = fp_border(color = "#E8E8E8", width = 0.5)) |>
  fontsize(size = 9, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  width(j = 1, width = 6.5 / 2.54) |>
  width(j = 2, width = 6   / 2.54) |>
  width(j = 3, width = 6   / 2.54)
doc <- body_add_flextable(doc, ft_empf)

# ── Dokument speichern ──
docx_pfad <- file.path(output_ordner, "Block2c_HMM_Bericht.docx")
print(doc, target = docx_pfad)
cat("  ✓ Word-Bericht gespeichert:", docx_pfad, "\n\n")

cat("════════════════════════════════════════════════════════\n")
cat("FERTIG — alle Ausgabedateien erstellt:\n")
cat("  →", excel_pfad, "\n")
cat("  →", docx_pfad, "\n")
cat("  →", file.path(output_ordner, "hmm_nachtminuten.rds"), "\n")
cat("════════════════════════════════════════════════════════\n\n")
