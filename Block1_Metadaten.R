# ==============================================================
# Block 1 — Metadata Analysis
# VHF Hedgehog Telemetry, Lower Saxony
# ==============================================================
# Project:  Hedgehog VHF Telemetry, Wildtierstation Sachsenhagen / TiHo Hannover
# Author:   Natalie Steiner
#
# Content:
#   1. Load & clean metadata
#   2. Descriptive visualisations (individual plots, English)
#   3. German Word methods report
#   4. Excel overview (auto-generated from current data)
#
# Requirements:
#   - data/excel_files/data_igel.xlsx present
#
# Note: Chronobiology / rhythm analysis is handled exclusively in Block 3.
#       The linking of chrono results with metadata is also done there.
# ==============================================================
# Run from line 1 (Ctrl+Alt+R in RStudio).
# ==============================================================

pakete <- c("data.table", "ggplot2", "readxl", "scales", "patchwork",
            "officer", "flextable", "lubridate", "openxlsx")
neu <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(neu) > 0) install.packages(neu)

library(data.table)
library(ggplot2)
library(readxl)
library(scales)
library(patchwork)
library(officer)
library(flextable)
library(lubridate)
library(openxlsx)

cat("✓ Packages loaded\n\n")

# ──────────────────────────────────────────────────────────────
# SETTINGS
# ──────────────────────────────────────────────────────────────

projekt_root  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
output_ordner <- file.path(projekt_root, "output", "Block1_Metadaten")
dir.create(output_ordner, showWarnings = FALSE, recursive = TRUE)

# Colour palettes (English labels in plots)
pal_sex   <- c("Male" = "#2C7BB6", "Female" = "#D7191C")
pal_age   <- c("Juvenile" = "#FEC44F", "Adult" = "#2C7BB6")
pal_seas  <- c("Autumn" = "#E8A838", "Spring" = "#74C476", "Summer" = "#2CA25F")
pal_diag  <- c("Orphan"    = "#9E9AC8",
               "Parasites" = "#74C476",
               "Trauma"    = "#FC8D59",
               "Fungal"    = "#FDAE61",
               "Blindness" = "#D9D9D9",
               "Other"     = "#BDBDBD")

# Shared theme
theme_b1 <- theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey45"),
        panel.grid.minor = element_blank())

# Helper: save individual plot
save_plot <- function(p, name, w = 1600, h = 1200, res = 130) {
  pfad <- file.path(output_ordner, name)
  png(pfad, width = w, height = h, res = res)
  print(p)
  dev.off()
  cat("  ✓", name, "\n")
  invisible(pfad)
}

# ──────────────────────────────────────────────────────────────
# 1. LOAD & CLEAN METADATA
# ──────────────────────────────────────────────────────────────

cat("Loading metadata from data_igel.xlsx...\n")

meta_raw <- as.data.table(read_excel(
  file.path(projekt_root, "data", "excel_files", "data_igel.xlsx"), sheet = 1))
setnames(meta_raw, names(meta_raw), tolower(gsub(" ", "_", names(meta_raw))))

# Fix known typo in original column name
if ("weigth_entry" %in% names(meta_raw))
  setnames(meta_raw, "weigth_entry", "weight_entry")

meta_raw[, individual := trimws(individual)]
meta <- meta_raw[!is.na(id)]
cat("  Animals loaded:", nrow(meta), "\n")

# ── Age class ────────────────────────────────────────────────
# Juveniles = admitted as orphan OR admission weight < 300 g
# (300 g is the threshold below which hedgehogs are not winter-ready)
meta[, age_class := ifelse(
  grepl("orphan", tolower(diagnosis_main)) |
    (!is.na(weight_entry) & weight_entry < 300),
  "Juvenile", "Adult"
)]

# ── Diagnosis groups ─────────────────────────────────────────
meta[, diag_group := fcase(
  grepl("orphan",            tolower(diagnosis_main)),                     "Orphan",
  grepl("trauma",            tolower(diagnosis_main)) &
    !grepl("blind",          tolower(diagnosis_main)),                     "Trauma",
  grepl("blind",             tolower(diagnosis_main)),                     "Blindness",
  grepl("fungal",            tolower(diagnosis_main)),                     "Fungal",
  grepl("ecto|endo|parasit", tolower(diagnosis_main)),                     "Parasites",
  default = "Other"
)]

# ── Release season ───────────────────────────────────────────
meta[, release_month := month(date_release)]
meta[, season := fcase(
  release_month %in% c(9, 10, 11), "Autumn",
  release_month %in% c(3,  4,  5), "Spring",
  release_month %in% c(6,  7,  8), "Summer",
  default = "Winter"
)]

# Keep German label for Excel (backwards compatibility)
meta[, saison_auswild := fcase(
  release_month %in% c(9, 10, 11), "Herbst",
  release_month %in% c(3,  4,  5), "Frühling",
  release_month %in% c(6,  7,  8), "Sommer",
  default = "Winter"
)]

# ── Numeric columns ──────────────────────────────────────────
for (col in c("weight_entry", "tagging_weight", "weight_gain", "time_reha")) {
  meta[, (col) := suppressWarnings(as.numeric(get(col)))]
}

meta[, igel := individual]

# ── Console summary ──────────────────────────────────────────
cat("\n── Sample overview ──\n")
cat("  Total:       ", nrow(meta), "hedgehogs\n")
cat("  Males:       ", meta[sex == "Male",    .N], "\n")
cat("  Females:     ", meta[sex == "Female",  .N], "\n")
cat("  Juveniles:   ", meta[age_class == "Juvenile", .N], "\n")
cat("  Adults:      ", meta[age_class == "Adult",    .N], "\n")
reha_v <- meta[!is.na(time_reha), time_reha]
cat(sprintf("  Rehab (d):   n=%d | Median=%.0f | Mean=%.1f | SD=%.1f | Min=%.0f | Max=%.0f\n",
    length(reha_v), median(reha_v), mean(reha_v), sd(reha_v), min(reha_v), max(reha_v)))
cat("\n── Diagnosis groups ──\n")
print(meta[, .N, by = diag_group][order(-N)])

# ══════════════════════════════════════════════════════════════
# 2. DESCRIPTIVE PLOTS  (individual files, English)
# ══════════════════════════════════════════════════════════════

cat("\n════════════════════════════════════\n")
cat("2. Descriptive plots\n")
cat("════════════════════════════════════\n")

# ── Plot 01: Diagnosis groups ─────────────────────────────────
p01 <- ggplot(meta, aes(x = reorder(diag_group, -table(diag_group)[diag_group]),
                         fill = diag_group)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = pal_diag, guide = "none") +
  labs(title    = "Diagnosis groups",
       subtitle  = paste0("N = ", nrow(meta), " hedgehogs"),
       x = NULL, y = "Number of hedgehogs") +
  theme_b1 +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p01, "meta_01_diagnosis.png")

# ── Plot 02: Rehabilitation duration by diagnosis ─────────────
p02 <- ggplot(meta[!is.na(time_reha)],
              aes(x = reorder(diag_group, time_reha, median, na.rm = TRUE),
                  y = time_reha, fill = diag_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = pal_diag, guide = "none") +
  labs(title    = "Rehabilitation duration by diagnosis",
       subtitle  = "Ordered by median; points = individual animals",
       x = NULL, y = "Days in rehabilitation") +
  theme_b1 +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p02, "meta_02_reha_duration.png")

# ── Plot 03: Weight development ───────────────────────────────
meta_wt <- melt(
  meta[!is.na(weight_entry) & !is.na(tagging_weight),
       .(igel, diag_group, age_class,
         `Admission` = weight_entry,
         `Release`   = tagging_weight)],
  id.vars = c("igel", "diag_group", "age_class"),
  variable.name = "Timepoint", value.name = "Weight_g"
)
p03 <- ggplot(meta_wt, aes(x = Timepoint, y = Weight_g,
                             group = igel, color = age_class)) +
  geom_line(alpha = 0.45, linewidth = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = pal_age, name = "Age class") +
  labs(title    = "Body weight: admission vs. release",
       subtitle  = "Lines connect the same individual",
       x = NULL, y = "Body weight (g)") +
  theme_b1
save_plot(p03, "meta_03_weight_development.png")

# ── Plot 04: Release month ────────────────────────────────────
p04 <- ggplot(meta, aes(x = factor(release_month, levels = 1:12,
                                    labels = month.abb),
                         fill = season)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = pal_seas, name = "Season") +
  labs(title    = "Release month",
       subtitle  = "Colour = season of release",
       x = "Month", y = "Number of hedgehogs") +
  theme_b1
save_plot(p04, "meta_04_release_month.png")

# ── Plot 05: Sex distribution ─────────────────────────────────
sex_tbl <- meta[, .N, by = .(sex, age_class)]
p05 <- ggplot(sex_tbl, aes(x = sex, y = N, fill = age_class)) +
  geom_col(position = "stack", width = 0.5, alpha = 0.85) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5),
            size = 4, color = "white", fontface = "bold") +
  scale_fill_manual(values = pal_age, name = "Age class") +
  scale_x_discrete(labels = c("Male" = "Male", "Female" = "Female")) +
  labs(title    = "Sex and age class distribution",
       subtitle  = paste0("N = ", nrow(meta)),
       x = NULL, y = "Number of hedgehogs") +
  theme_b1
save_plot(p05, "meta_05_sex_age.png")

# ── Plot 06: Tracking duration (Gantt) ───────────────────────
meta_gantt <- meta[!is.na(date_release) & !is.na(tagging_period)][order(date_release)]
meta_gantt[, igel_f  := factor(igel, levels = rev(igel))]
meta_gantt[, t_start := date_release]
meta_gantt[, t_end   := date_release + tagging_period]

p06 <- ggplot(meta_gantt, aes(y = igel_f, color = season)) +
  geom_segment(aes(x = t_start, xend = t_end, yend = igel_f),
               linewidth = 5, alpha = 0.75) +
  geom_point(aes(x = t_start), shape = 21, fill = "white",
             size = 3, stroke = 1.2) +
  scale_color_manual(values = pal_seas, name = "Release season") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  labs(title    = "Monitoring periods — all hedgehogs",
       subtitle  = "Dot = release date | Bar = monitoring period | Colour = season",
       x = NULL, y = NULL) +
  theme_b1 +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major.y = element_blank())
save_plot(p06, "meta_06_gantt.png", w = 1800, h = 1400)

# ══════════════════════════════════════════════════════════════
# 3. METHODENBERICHT (Word, Deutsch)
# ══════════════════════════════════════════════════════════════

cat("\n════════════════════════════════════\n")
cat("3. Methodenbericht (Word)\n")
cat("════════════════════════════════════\n")

# ── Word-Dokument ──────────────────────────────────────────────
doc <- read_docx()

h1 <- function(doc, txt) body_add_par(doc, txt, style = "heading 1")
h2 <- function(doc, txt) body_add_par(doc, txt, style = "heading 2")
h3 <- function(doc, txt) body_add_par(doc, txt, style = "heading 3")
p  <- function(doc, txt) body_add_par(doc, txt, style = "Normal")
br <- function(doc)      body_add_break(doc)

# ── Titelseite ────────────────────────────────────────────────
doc <- h1(doc, "Block 1 — Methodenbericht")
doc <- p(doc, "Metadatenanalyse | VHF Igelbesenderung Niedersachsen")
doc <- p(doc, "Projekt: Wildtierstation Sachsenhagen / TiHo Hannover")
doc <- p(doc, paste0("Erstellt: ", format(Sys.time(), "%d.%m.%Y, %H:%M"), " Uhr"))
doc <- p(doc, paste0("R-Version: ", R.version$major, ".", R.version$minor))
doc <- br(doc)

# ── Vorbemerkung ──────────────────────────────────────────────
doc <- h1(doc, "Vorbemerkung: Abgrenzung Block 1 zu Block 3")
doc <- p(doc, paste0(
  "Block 1 (Block1_Metadaten.R) ist ausschliesslich ein Metadaten-Analyseskript. ",
  "Es fuehrt KEINE Rhythmusanalysen und KEINE Gruppenvergleiche chronobiologischer ",
  "Parameter durch. Die Berechnung von Cosinor-Amplitude, Akrophase und Rayleigh-rho ",
  "sowie alle statistischen Vergleiche dieser Parameter (z.B. nach Geschlecht, ",
  "Altersklasse, Diagnose oder Saison) erfolgen vollstaendig in Block 3 ",
  "(Block3_Chronobiologie.R). Die Verknuepfung von Metadaten mit Rhythmusergebnissen ",
  "findet ebenfalls erst in Block 3 statt, damit alle chronobiologischen Analysen ",
  "an einem Ort zentralisiert sind und Block 1 uebersichtlich und eigenstaendig bleibt."
))
doc <- br(doc)

# ══════════════════════════════════════════════════════════════
doc <- h1(doc, "1  Datengrundlage")
# ══════════════════════════════════════════════════════════════

doc <- h2(doc, "1.1  Quelldatei")
doc <- p(doc, paste0(
  "Alle Metadaten stammen aus der Excel-Datei data_igel.xlsx ",
  "(Pfad: data/excel_files/data_igel.xlsx). Diese Datei wird mit read_excel() ",
  "aus dem Paket readxl eingelesen. Das erste Arbeitsblatt (sheet = 1) wird ",
  "verwendet. Die Datei enthaelt fuer jedes besenderte Tier eine Zeile mit ",
  "klinischen und administrativen Daten aus der Pflegestation."
))
doc <- p(doc, paste0(
  "Warum diese Quelle: Die Klinikdaten der Wildtierstation sind die einzige ",
  "vollstaendige Informationsquelle ueber Einlieferungsdiagnose, Gewichtsverlauf, ",
  "Reha-Dauer und Auswilderungszeitpunkt. Eine direkte Verknuepfung mit den VHF-Daten ",
  "ist ueber die Tier-ID (individual / igel) moeglich."
))

doc <- h2(doc, "1.2  Datenbereinigung")
doc <- p(doc, paste0(
  "Nach dem Einlesen werden alle Spaltennamen in Kleinbuchstaben konvertiert und ",
  "Leerzeichen durch Unterstriche ersetzt (tolower + gsub). Ein bekannter Tippfehler ",
  "in der Originalquelle (Spalte 'weigth_entry') wird automatisch korrigiert zu ",
  "'weight_entry'. Igelnamen werden von fuehrenden und nachgestellten Leerzeichen ",
  "befreit (trimws). Zeilen ohne ID-Eintrag werden entfernt (!is.na(id))."
))
doc <- p(doc, paste0(
  "Die numerischen Spalten weight_entry, tagging_weight, weight_gain und time_reha ",
  "werden explizit in numerische Werte konvertiert (as.numeric), da Excel-Importe ",
  "diese Felder gelegentlich als Character einlesen. Warnungen bei nicht konvertierbaren ",
  "Werten werden unterdrueckt (suppressWarnings), betroffene Felder erhalten NA."
))
doc <- p(doc, paste0(
  "Limitation: Die Metadaten koennen manuelle Eingabefehler enthalten (z.B. falsch ",
  "erfasste Gewichte oder Daten). Es findet keine automatische Plausibilitaetspruefung ",
  "statt. Ungewoehnliche Werte (z.B. negatives Gewicht oder Reha-Dauer = 0) werden ",
  "nicht herausgefiltert und koennten statistische Analysen verzerren."
))

# ══════════════════════════════════════════════════════════════
doc <- h1(doc, "2  Ableitung von Analysevariablen")
# ══════════════════════════════════════════════════════════════

doc <- h2(doc, "2.1  Altersklasse (age_class)")
doc <- p(doc, paste0(
  "Jedes Tier wird einer von zwei Altersklassen zugewiesen: Juvenile oder Adult. ",
  "Ein Tier gilt als Jungtier (Juvenile), wenn mindestens eine der folgenden ",
  "Bedingungen erfuellt ist: (a) Die Einlieferungsdiagnose enthaelt das Schluesselwort ",
  "'orphan' (Waise), oder (b) das Einlieferungsgewicht betraegt weniger als 300 g. ",
  "Tiere, die keine dieser Bedingungen erfuellen, werden als Adult klassifiziert."
))
doc <- p(doc, paste0(
  "Begruendung: Der 300-g-Schwellenwert ist ein in der Igelrehabilitation ",
  "etablierter Richtwert. Tiere unter 300 g gelten im Herbst als nicht winterfit, ",
  "da sie nicht genuegend Fettreserven fuer den Winterschlaf aufgebaut haben. ",
  "Waisen werden unabhaengig vom Gewicht als Jungtiere klassifiziert, da sie ohne ",
  "Mutter aufgezogen wurden und ihr Entwicklungsstand unsicher ist."
))
doc <- p(doc, paste0(
  "Limitation: Eine direkte Altersbestimmung (z.B. anhand von Geburtsdatum oder ",
  "Zahnstatus) ist nicht moeglich. Der Gewichtsschwellenwert ist ein Proxy, ",
  "kein exaktes Altermass. Tiere, die knapp ueber 300 g lagen, koennten trotzdem ",
  "Jungtiere gewesen sein. Ausserdem ist weight_entry bei manchen Tieren nicht ",
  "vorhanden (NA), was die Klassifizierung unsicherer macht."
))

doc <- h2(doc, "2.2  Diagnosegruppen (diag_group)")
doc <- p(doc, paste0(
  "Die klinischen Einlieferungsdiagnosen aus der Datenbank sind sehr heterogen ",
  "und enthalten viele Einzelkategorien, die fuer statistische Analysen zu kleinen ",
  "Gruppen fuehren wuerden. Daher werden sie in sechs vereinfachte Gruppen ",
  "zusammengefasst, basierend auf Schluesselwoertern in der Spalte diagnosis_main:"
))
doc <- p(doc, paste0(
  "'Orphan': Diagnose enthaelt 'orphan'. ",
  "'Trauma': Diagnose enthaelt 'trauma', aber nicht 'blind'. ",
  "'Blindness': Diagnose enthaelt 'blind'. ",
  "'Fungal': Diagnose enthaelt 'fungal'. ",
  "'Parasites': Diagnose enthaelt 'ecto', 'endo' oder 'parasit'. ",
  "'Other': Alle nicht zugeordneten Diagnosen. ",
  "Die Zuordnung erfolgt hierarchisch: Eine Diagnose, die sowohl 'trauma' als auch ",
  "'blind' enthaelt, wird als 'Blindness' klassifiziert."
))
doc <- p(doc, paste0(
  "Begruendung: Die Vereinfachung erhoet die statistisch verwertbaren Gruppengroessen ",
  "und macht inhaltlich aehnliche Faelle vergleichbar. Tiere mit aehnlichen ",
  "Grunderkrankungen sind am ehesten in biologisch sinnvolle Vergleichsgruppen einteilbar."
))
doc <- p(doc, paste0(
  "Limitation: Die Vereinfachung kann zu Informationsverlust fuehren. Tiere mit ",
  "Mehrfachdiagnosen werden nur einer Gruppe zugeordnet. Die Gruppe 'Other' ist ",
  "inhaltlich heterogen und schwer interpretierbar. Die Schluesselwort-Suche ist ",
  "abhaengig von konsistenter Schreibweise in der Quelldatei."
))

doc <- h2(doc, "2.3  Auswilderungssaison (season)")
doc <- p(doc, paste0(
  "Aus dem Auswilderungsdatum (date_release) wird der Kalendermonat extrahiert ",
  "und einer von vier Jahreszeiten zugeordnet: Fruehling (Maerz bis Mai), ",
  "Sommer (Juni bis August), Herbst (September bis November), Winter (Dezember bis Februar). ",
  "Winter ist in dieser Stichprobe nicht vertreten."
))
doc <- p(doc, paste0(
  "Begruendung: Die Jahreszeit der Auswilderung beeinflusst Umgebungsbedingungen ",
  "(Temperatur, Taglaenge, Nahrungsangebot) und damit moeglicherweise die ",
  "Uebertragbarkeit von Stationserfahrungen auf die Wildnis. Herbst-Auswilderungen ",
  "sind fuer Igel besonders kritisch, da kurze Zeit vor dem Winterschlaf bleibt."
))
doc <- p(doc, paste0(
  "Limitation: Die astronomischen Jahreszeiten (Sonnenwende/Tagundnachtgleiche) ",
  "wurden hier durch feste Monatsgrenzen approximiert. Da keine Winter-Auswilderungen ",
  "in der Stichprobe sind, ist ein vollstaendiger Saisonvergleich nicht moeglich."
))

# ══════════════════════════════════════════════════════════════
doc <- h1(doc, "3  Deskriptive Visualisierungen")
# ══════════════════════════════════════════════════════════════

doc <- p(doc, paste0(
  "Block 1 erstellt sechs deskriptive Grafiken zur Stichprobenbeschreibung. ",
  "Alle Grafiken werden als einzelne PNG-Dateien im Ausgabeordner gespeichert ",
  "(meta_01_diagnosis.png bis meta_06_gantt.png). Die Grafiken sind auf Englisch ",
  "beschriftet, um direkte Verwendung in wissenschaftlichen Publikationen zu ",
  "ermoeglichen. Alle verwenden ein einheitliches minimalistisches Theme ",
  "(theme_minimal aus ggplot2)."
))

doc <- h2(doc, "meta_01_diagnosis.png — Diagnosegruppen")
doc <- p(doc, paste0(
  "Balkendiagramm der Haeufigkeit jeder Diagnosegruppe. Die Balken sind nach ",
  "Haeufigkeit absteigend sortiert. Absolute Zahlen werden ueber den Balken ",
  "angezeigt. Jede Diagnosegruppe hat eine eigene Farbe (farbliche Palette pal_diag). ",
  "Zweck: Schneller Ueberblick ueber die Zusammensetzung der Stichprobe nach ",
  "Einlieferungsgrund. Hilft einzuschaetzen, ob bestimmte Diagnosegruppen ",
  "statistisch analysierfaehig sind (Mindestgruppengroesse)."
))

doc <- h2(doc, "meta_02_reha_duration.png — Reha-Dauer nach Diagnose")
doc <- p(doc, paste0(
  "Boxplot der Aufenthaltsdauer (days in rehabilitation) fuer jede Diagnosegruppe, ",
  "geordnet nach Median. Einzelne Tiere werden als Punkte (jitter) ueber den Boxplots ",
  "angezeigt. Tiere ohne bekannte Reha-Dauer (NA) werden ausgeschlossen. ",
  "Zweck: Zeigt, ob bestimmte Diagnosen mit laengeren oder kuerzeren Klinikaufenthalten ",
  "verbunden sind. Waisen sind typischerweise laenger in der Station als Traumatiere."
))

doc <- h2(doc, "meta_03_weight_development.png — Gewichtsentwicklung")
doc <- p(doc, paste0(
  "Verbundenes Punkt-Linien-Diagramm (connected dots): Fuer jedes Tier werden ",
  "Einlieferungsgewicht und Auswilderungsgewicht als Punkte dargestellt, verbunden ",
  "durch eine Linie. Die Farbe kodiert die Altersklasse (Juvenile / Adult). ",
  "Nur Tiere mit beiden Gewichtswerten werden dargestellt. ",
  "Zweck: Visualisiert den Gewichtsverlauf waehrend der Rehabilitation auf ",
  "Einzeltierebene. Abwaertslinien zeigen Gewichtsverlust, was als Qualitaetsindikator ",
  "dienen kann. Grosse Streuung im Einlieferungsgewicht reflektiert die heterogene Stichprobe."
))

doc <- h2(doc, "meta_04_release_month.png — Auswilderungsmonat")
doc <- p(doc, paste0(
  "Balkendiagramm der Auswilderungen nach Kalendermonat (Januar bis Dezember). ",
  "Die Balkenfarbe kodiert die Jahreszeit (Spring / Summer / Autumn). ",
  "Zweck: Zeigt die zeitliche Verteilung der Auswilderungen. Haeufungen im Herbst ",
  "reflektieren den typischen Saisonverlauf der Igelrehabilitation (viele ",
  "Jungtier-Aufnahmen im Spaetsommer/Herbst)."
))

doc <- h2(doc, "meta_05_sex_age.png — Geschlecht und Altersklasse")
doc <- p(doc, paste0(
  "Gestapeltes Balkendiagramm: Fuer Maennchen und Weibchen wird jeweils die ",
  "Anzahl Jungtiere (Juvenile) und Adulte (Adult) als gestapelte Balken dargestellt. ",
  "Absolute Zahlen werden innerhalb der Segmente angezeigt. ",
  "Zweck: Kompakte Uebersicht ueber die demographische Struktur der Stichprobe."
))

doc <- h2(doc, "meta_06_gantt.png — Beobachtungszeitraum (Gantt)")
doc <- p(doc, paste0(
  "Gantt-Diagramm: Fuer jedes Tier wird ein horizontaler Balken vom Auswilderungsdatum ",
  "(Punkt) bis zum letzten VHF-Signal dargestellt. Die Tiere sind chronologisch nach ",
  "Auswilderungsdatum sortiert. Die Balkenfarbe kodiert die Auswilderungssaison. ",
  "Nur Tiere mit bekannter Tracking-Dauer werden dargestellt. ",
  "Zweck: Gibt auf einen Blick Auskunft ueber Beobachtungsdauer, zeitliche Ueberlappungen ",
  "und saisonale Verteilung. Sehr kurze Balken (<7 Tage) signalisieren Tiere, die fuer ",
  "Aktivitaetsanalysen ausgeschlossen werden muessen."
))

# ══════════════════════════════════════════════════════════════
doc <- h1(doc, "4  Excel-Uebersicht (Block1_Uebersicht.xlsx)")
# ══════════════════════════════════════════════════════════════

doc <- p(doc, paste0(
  "Die Excel-Datei wird vollstaendig automatisch aus den aktuellen Metadaten ",
  "generiert und bei jedem Skriptdurchlauf neu erstellt (overwrite = TRUE). ",
  "Sie enthaelt vier Arbeitsblaetter:"
))

doc <- h2(doc, "Sheet 'Overview' — Dashboard")
doc <- p(doc, paste0(
  "Zusammenfassung der wichtigsten Kenngzahlen der Stichprobe: ",
  "Aufteilung nach Auswilderungssaison (N, Geschlecht, mittlere Tracking-Dauer, ",
  "mittlere Reha-Dauer, mittlere Gewichtszunahme, Tier-IDs) und nach ",
  "Diagnosegruppe (N, Geschlecht, Altersklasse, mittlere Reha-Dauer, Gewichtszunahme). ",
  "Ausserdem deskriptive Statistiken (N, Min, Max, Median, Mittelwert, SD) fuer ",
  "Einlieferungsgewicht, Auswilderungsgewicht, Gewichtszunahme, Reha-Dauer ",
  "und Tracking-Dauer."
))

doc <- h2(doc, "Sheet 'All animals' — Alle Tiere")
doc <- p(doc, paste0(
  "Eine Zeile pro Tier mit allen relevanten Metadaten: Tier-ID, Geschlecht, ",
  "Altersklasse, Diagnosegruppe, Originaldiagnose, Einlieferungsdatum, Gewichte, ",
  "Gewichtszunahme, Reha-Dauer, Auswilderungsdatum, Saison, Tracking-Dauer. ",
  "Jungtiere werden gelb hinterlegt, Tiere mit negativer Gewichtszunahme in Rot."
))

doc <- h2(doc, "Sheet 'By season' — Nach Saison")
doc <- p(doc, paste0(
  "Tiere gruppiert nach Auswilderungssaison. Jede Saison hat einen farbigen ",
  "Gruppenheader, gefolgt von einer Tabelle der Tiere dieser Saison und einer ",
  "Mittelwertzeile. Ermoeglicht direkten Vergleich einzelner Tiere innerhalb ",
  "einer Saison."
))

doc <- h2(doc, "Sheet 'By diagnosis' — Nach Diagnose")
doc <- p(doc, paste0(
  "Analoge Struktur zu 'By season', aber gruppiert nach Diagnosegruppe. ",
  "Jede Diagnosegruppe hat einen eigenen farbigen Block. Erleichtert ",
  "den Vergleich von Tieren mit aehnlicher Einlieferungsursache."
))

doc <- br(doc)

# ── Speichern ─────────────────────────────────────────────────
bericht_pfad <- file.path(output_ordner, "Block1_Methodenbericht.docx")
print(doc, target = bericht_pfad)
cat("✓ Methodenbericht gespeichert:", basename(bericht_pfad), "\n\n")

# ══════════════════════════════════════════════════════════════
# 4. EXCEL OVERVIEW  (auto-generated from current data)
# ══════════════════════════════════════════════════════════════

cat("════════════════════════════════════\n")
cat("4. Excel overview\n")
cat("════════════════════════════════════\n\n")

# Summary counts for Overview sheet
n_total  <- nrow(meta)
n_male   <- meta[sex == "Male",          .N]
n_female <- meta[sex == "Female",        .N]
n_juv    <- meta[age_class == "Juvenile",.N]
n_adult  <- meta[age_class == "Adult",   .N]

# ── Helper functions ───────────────────────────────────────────
xls_stats <- function(v) {
  v <- as.numeric(v[!is.na(v)])
  if (length(v) == 0) return(c(N=0, Min=NA, Max=NA, Median=NA, Mean=NA, SD=NA))
  c(N = length(v), Min = round(min(v), 1), Max = round(max(v), 1),
    Median = round(median(v), 1), Mean = round(mean(v), 1), SD = round(sd(v), 1))
}
safe_mean <- function(v) {
  v <- as.numeric(v[!is.na(v)])
  if (length(v) == 0) return(NA_real_)
  round(mean(v), 1)
}

# ── Colour constants ───────────────────────────────────────────
C_DARK  <- "#1a2e4a"
C_MID   <- "#2C5F8A"
C_LIGHT <- "#D6E4F0"
C_LBLUE <- "#EBF3FB"
C_GREY  <- "#F5F5F5"
C_WHITE <- "#FFFFFF"
C_YELL  <- "#FFF9C4"
C_RED   <- "#FFE0E0"

hdr_style <- createStyle(fontName="Arial", fontSize=10, fontColour=C_WHITE,
  fgFill=C_MID, halign="CENTER", valign="center",
  textDecoration="bold", wrapText=TRUE,
  border="TopBottomLeftRight", borderColour="#AABFD4")

title_style <- createStyle(fontName="Arial", fontSize=14, fontColour=C_WHITE,
  fgFill=C_DARK, textDecoration="bold", halign="left", valign="center", indent=1)

sub_style <- createStyle(fontName="Arial", fontSize=10, fontColour="#666666",
  fgFill=C_LBLUE, halign="left", valign="center",
  textDecoration="italic", indent=1)

section_style <- createStyle(fontName="Arial", fontSize=11, fontColour=C_WHITE,
  fgFill=C_MID, textDecoration="bold", halign="left", valign="center", indent=1)

cell_style <- function(bg=C_WHITE, halign="center") {
  createStyle(fontName="Arial", fontSize=10, fgFill=bg,
    halign=halign, valign="center",
    border="TopBottomLeftRight", borderColour="#AABFD4", wrapText=TRUE)
}

wb_xl <- createWorkbook()

# ──────────────────────────────────────────────────────────────
# Sheet 1: Overview (Dashboard)
# ──────────────────────────────────────────────────────────────
addWorksheet(wb_xl, "Overview", gridLines=FALSE)
ws_u <- "Overview"

mergeCells(wb_xl, ws_u, cols=1:9, rows=1)
writeData(wb_xl, ws_u,
  paste0("Block 1 — Metadata Overview  |  VHF Hedgehog Telemetry  |  Updated: ",
         format(Sys.Date(), "%d.%m.%Y")),
  startRow=1, startCol=1)
addStyle(wb_xl, ws_u, title_style, rows=1, cols=1:9, gridExpand=TRUE)
setRowHeights(wb_xl, ws_u, rows=1, heights=28)

n_vhf <- nrow(meta)
n_is  <- sum(!is.na(meta$tagging_period) & meta$tagging_period >= 7)
mergeCells(wb_xl, ws_u, cols=1:9, rows=2)
writeData(wb_xl, ws_u,
  paste0("N = ", n_total, " animals  |  ",
         n_male, " males, ", n_female, " females  |  ",
         n_juv,  " juveniles, ", n_adult, " adults  |  ",
         n_is, " animals with ≥ 7 days tracking (IS-eligible)"),
  startRow=2, startCol=1)
addStyle(wb_xl, ws_u, sub_style, rows=2, cols=1:9, gridExpand=TRUE)
setRowHeights(wb_xl, ws_u, rows=2, heights=18)
setRowHeights(wb_xl, ws_u, rows=3, heights=6)

# Season section
cur_row <- 4
mergeCells(wb_xl, ws_u, cols=1:9, rows=cur_row)
writeData(wb_xl, ws_u, "Season of release", startRow=cur_row, startCol=1)
addStyle(wb_xl, ws_u, section_style, rows=cur_row, cols=1:9, gridExpand=TRUE)
setRowHeights(wb_xl, ws_u, rows=cur_row, heights=20); cur_row <- cur_row + 1

sai_hdrs <- c("Season","N","Males","Females",
              "Mean tracking (d)","Mean rehab (d)","Mean weight gain (g)","Animal IDs")
writeData(wb_xl, ws_u, as.data.frame(t(sai_hdrs)), startRow=cur_row, startCol=1, colNames=FALSE)
addStyle(wb_xl, ws_u, hdr_style, rows=cur_row, cols=1:8, gridExpand=TRUE)
setRowHeights(wb_xl, ws_u, rows=cur_row, heights=18); cur_row <- cur_row + 1

for (sai in c("Spring","Summer","Autumn")) {
  sub <- meta[season == sai]
  if (nrow(sub) == 0) next
  bg <- if (cur_row %% 2 == 0) C_GREY else C_WHITE
  rv  <- data.frame(Season=sai, N=nrow(sub),
    Males=sub[sex=="Male",.N], Females=sub[sex=="Female",.N],
    Track=safe_mean(sub$tagging_period), Reha=safe_mean(sub$time_reha),
    WtGain=safe_mean(sub$weight_gain), IDs=paste(sub$igel, collapse=", "))
  writeData(wb_xl, ws_u, rv, startRow=cur_row, startCol=1, colNames=FALSE)
  addStyle(wb_xl, ws_u, cell_style(bg,"center"), rows=cur_row, cols=1:7, gridExpand=TRUE)
  addStyle(wb_xl, ws_u, cell_style(bg,"left"),   rows=cur_row, cols=8, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_u, rows=cur_row, heights=16); cur_row <- cur_row + 1
}
setRowHeights(wb_xl, ws_u, rows=cur_row, heights=6); cur_row <- cur_row + 1

# Diagnosis section
mergeCells(wb_xl, ws_u, cols=1:9, rows=cur_row)
writeData(wb_xl, ws_u, "Diagnosis groups", startRow=cur_row, startCol=1)
addStyle(wb_xl, ws_u, section_style, rows=cur_row, cols=1:9, gridExpand=TRUE)
setRowHeights(wb_xl, ws_u, rows=cur_row, heights=20); cur_row <- cur_row + 1

diag_hdrs <- c("Diagnosis","N","Males","Females","Juvenile","Adult",
               "Mean rehab (d)","Mean weight gain (g)")
writeData(wb_xl, ws_u, as.data.frame(t(diag_hdrs)), startRow=cur_row, startCol=1, colNames=FALSE)
addStyle(wb_xl, ws_u, hdr_style, rows=cur_row, cols=1:8, gridExpand=TRUE)
setRowHeights(wb_xl, ws_u, rows=cur_row, heights=18); cur_row <- cur_row + 1

for (dg in c("Parasites","Trauma","Orphan","Blindness","Fungal","Other")) {
  sub <- meta[diag_group == dg]
  if (nrow(sub) == 0) next
  bg <- if (cur_row %% 2 == 0) C_GREY else C_WHITE
  rv <- data.frame(Diag=dg, N=nrow(sub),
    Males=sub[sex=="Male",.N], Females=sub[sex=="Female",.N],
    Juv=sub[age_class=="Juvenile",.N], Adult=sub[age_class=="Adult",.N],
    Reha=safe_mean(sub$time_reha), WtGain=safe_mean(sub$weight_gain))
  writeData(wb_xl, ws_u, rv, startRow=cur_row, startCol=1, colNames=FALSE)
  addStyle(wb_xl, ws_u, cell_style(bg,"left"),   rows=cur_row, cols=1, gridExpand=TRUE)
  addStyle(wb_xl, ws_u, cell_style(bg,"center"), rows=cur_row, cols=2:8, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_u, rows=cur_row, heights=16); cur_row <- cur_row + 1
}
setRowHeights(wb_xl, ws_u, rows=cur_row, heights=6); cur_row <- cur_row + 1

# Descriptive statistics section
mergeCells(wb_xl, ws_u, cols=1:9, rows=cur_row)
writeData(wb_xl, ws_u, "Descriptive statistics", startRow=cur_row, startCol=1)
addStyle(wb_xl, ws_u, section_style, rows=cur_row, cols=1:9, gridExpand=TRUE)
setRowHeights(wb_xl, ws_u, rows=cur_row, heights=20); cur_row <- cur_row + 1

stat_hdrs <- c("Variable","N","Min.","Max.","Median","Mean","SD","Unit")
writeData(wb_xl, ws_u, as.data.frame(t(stat_hdrs)), startRow=cur_row, startCol=1, colNames=FALSE)
addStyle(wb_xl, ws_u, hdr_style, rows=cur_row, cols=1:8, gridExpand=TRUE)
setRowHeights(wb_xl, ws_u, rows=cur_row, heights=18); cur_row <- cur_row + 1

stat_vars <- list(
  list("Admission weight",  meta$weight_entry,    "g"),
  list("Release weight",    meta$tagging_weight,  "g"),
  list("Weight gain",       meta$weight_gain,     "g"),
  list("Rehabilitation",    meta$time_reha,       "days"),
  list("Tracking duration", meta$tagging_period,  "days")
)
for (i in seq_along(stat_vars)) {
  item <- stat_vars[[i]]; s <- xls_stats(item[[2]])
  bg <- if (i %% 2 == 0) C_GREY else C_WHITE
  rv <- data.frame(Var=item[[1]], N=s["N"],
    Min=ifelse(is.na(s["Min"]),"—",as.character(s["Min"])),
    Max=ifelse(is.na(s["Max"]),"—",as.character(s["Max"])),
    Med=ifelse(is.na(s["Median"]),"—",as.character(s["Median"])),
    Mean=ifelse(is.na(s["Mean"]),"—",as.character(s["Mean"])),
    SD=ifelse(is.na(s["SD"]),"—",as.character(s["SD"])),
    Unit=item[[3]])
  writeData(wb_xl, ws_u, rv, startRow=cur_row, startCol=1, colNames=FALSE)
  addStyle(wb_xl, ws_u, cell_style(bg,"left"),   rows=cur_row, cols=1, gridExpand=TRUE)
  addStyle(wb_xl, ws_u, cell_style(bg,"center"), rows=cur_row, cols=2:8, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_u, rows=cur_row, heights=16); cur_row <- cur_row + 1
}
setColWidths(wb_xl, ws_u, cols=1:9, widths=c(22,8,10,10,12,12,16,32,5))

# ──────────────────────────────────────────────────────────────
# Sheet 2: All animals
# ──────────────────────────────────────────────────────────────
addWorksheet(wb_xl, "All animals", gridLines=FALSE)
ws_t <- "All animals"
freezePane(wb_xl, ws_t, firstRow=TRUE)

tbl_alle <- meta[order(date_release), .(
  `Animal ID`           = igel,
  Sex                   = sex,
  `Age class`           = age_class,
  `Diagnosis group`     = diag_group,
  `Diagnosis (original)`= diagnosis_main,
  `Admission date`      = format(date_entry,   "%d.%m.%Y"),
  `Adm. weight (g)`     = weight_entry,
  `Rel. weight (g)`     = tagging_weight,
  `Weight gain (g)`     = weight_gain,
  `Rehab (d)`           = time_reha,
  `Release date`        = format(date_release, "%d.%m.%Y"),
  Season                = season,
  `Tracking (d)`        = tagging_period,
  `VHF data`            = "Yes"
)]

writeData(wb_xl, ws_t, tbl_alle, startRow=1, startCol=1, headerStyle=hdr_style)
setRowHeights(wb_xl, ws_t, rows=1, heights=22)

for (i in seq_len(nrow(tbl_alle))) {
  rn <- i + 1
  is_juv <- tbl_alle[i, `Age class`] == "Juvenile"
  bg <- if (is_juv) C_YELL else if (i %% 2 == 0) C_GREY else C_WHITE
  addStyle(wb_xl, ws_t, cell_style(bg,"center"), rows=rn, cols=1:14, gridExpand=TRUE)
  wg <- tbl_alle[i, `Weight gain (g)`]
  if (!is.na(wg) && is.numeric(wg) && wg < 0)
    addStyle(wb_xl, ws_t, cell_style(C_RED,"center"), rows=rn, cols=9, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_t, rows=rn, heights=16)
}
setColWidths(wb_xl, ws_t, cols=1:14,
  widths=c(9,10,10,15,24,12,13,13,13,10,12,10,11,10))

# ──────────────────────────────────────────────────────────────
# Sheet 3: By season
# ──────────────────────────────────────────────────────────────
addWorksheet(wb_xl, "By season", gridLines=FALSE)
ws_s <- "By season"

seas_bg <- c(Spring="#C8E6C9", Summer="#FFF9C4", Autumn="#FFE0B2")
sai_col_hdrs <- c("Animal ID","Sex","Age","Diagnosis","Rehab (d)",
                  "Adm. wt (g)","Rel. wt (g)","Wt gain (g)","Tracking (d)")
cur_row <- 1
for (sai in c("Spring","Summer","Autumn")) {
  sub <- meta[season == sai]
  if (nrow(sub) == 0) next
  sc <- seas_bg[sai]
  mergeCells(wb_xl, ws_s, cols=1:9, rows=cur_row)
  writeData(wb_xl, ws_s, paste0(sai, "  (N = ", nrow(sub), ")"),
            startRow=cur_row, startCol=1)
  addStyle(wb_xl, ws_s,
    createStyle(fontName="Arial", fontSize=12, fontColour=C_DARK,
                fgFill=sc, textDecoration="bold",
                halign="left", valign="center", indent=1),
    rows=cur_row, cols=1:9, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_s, rows=cur_row, heights=22); cur_row <- cur_row + 1
  writeData(wb_xl, ws_s, as.data.frame(t(sai_col_hdrs)),
            startRow=cur_row, startCol=1, colNames=FALSE)
  addStyle(wb_xl, ws_s, hdr_style, rows=cur_row, cols=1:9, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_s, rows=cur_row, heights=18); cur_row <- cur_row + 1
  tbl_s <- sub[order(date_release), .(
    igel, sex, age_class, diag_group,
    time_reha, weight_entry, tagging_weight, weight_gain, tagging_period)]
  for (j in seq_len(nrow(tbl_s))) {
    bg <- if (j %% 2 == 0) C_GREY else C_WHITE
    writeData(wb_xl, ws_s, tbl_s[j], startRow=cur_row, startCol=1, colNames=FALSE)
    addStyle(wb_xl, ws_s, cell_style(bg,"center"), rows=cur_row, cols=1:9, gridExpand=TRUE)
    setRowHeights(wb_xl, ws_s, rows=cur_row, heights=15); cur_row <- cur_row + 1
  }
  oe <- data.frame(Lbl=paste0("Mean ", sai), S2="—", S3="—", S4="—",
    S5=safe_mean(sub$time_reha), S6=safe_mean(sub$weight_entry),
    S7=safe_mean(sub$tagging_weight), S8=safe_mean(sub$weight_gain),
    S9=safe_mean(sub$tagging_period))
  writeData(wb_xl, ws_s, oe, startRow=cur_row, startCol=1, colNames=FALSE)
  addStyle(wb_xl, ws_s,
    createStyle(fontName="Arial", fontSize=9, fgFill=sc,
                textDecoration="bold", halign="center", valign="center",
                border="TopBottomLeftRight", borderColour="#AABFD4"),
    rows=cur_row, cols=1:9, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_s, rows=cur_row, heights=16); cur_row <- cur_row + 2
}
setColWidths(wb_xl, ws_s, cols=1:9, widths=c(10,10,10,14,9,13,13,14,11))

# ──────────────────────────────────────────────────────────────
# Sheet 4: By diagnosis
# ──────────────────────────────────────────────────────────────
addWorksheet(wb_xl, "By diagnosis", gridLines=FALSE)
ws_d <- "By diagnosis"

diag_bg <- c(Parasites="#E3F2FD", Trauma="#FCE4EC", Orphan="#F3E5F5",
             Blindness="#EEEEEE", Fungal="#E8F5E9", Other="#FAFAFA")
diag_col_hdrs <- c("Animal ID","Sex","Age","Season",
                   "Rehab (d)","Adm. wt (g)","Rel. wt (g)","Wt gain (g)","VHF data")
cur_row_d <- 1
for (dg in c("Parasites","Trauma","Orphan","Blindness","Fungal","Other")) {
  sub <- meta[diag_group == dg]
  if (nrow(sub) == 0) next
  dc <- diag_bg[dg]
  mergeCells(wb_xl, ws_d, cols=1:9, rows=cur_row_d)
  writeData(wb_xl, ws_d, paste0(dg, "  (N = ", nrow(sub), ")"),
            startRow=cur_row_d, startCol=1)
  addStyle(wb_xl, ws_d,
    createStyle(fontName="Arial", fontSize=12, fontColour=C_DARK,
                fgFill=dc, textDecoration="bold",
                halign="left", valign="center", indent=1),
    rows=cur_row_d, cols=1:9, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_d, rows=cur_row_d, heights=22); cur_row_d <- cur_row_d + 1
  writeData(wb_xl, ws_d, as.data.frame(t(diag_col_hdrs)),
            startRow=cur_row_d, startCol=1, colNames=FALSE)
  addStyle(wb_xl, ws_d, hdr_style, rows=cur_row_d, cols=1:9, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_d, rows=cur_row_d, heights=18); cur_row_d <- cur_row_d + 1
  tbl_d <- sub[order(date_release), .(
    igel, sex, age_class, season,
    time_reha, weight_entry, tagging_weight, weight_gain, vhf="Yes")]
  for (j in seq_len(nrow(tbl_d))) {
    bg <- if (j %% 2 == 0) C_GREY else C_WHITE
    writeData(wb_xl, ws_d, tbl_d[j], startRow=cur_row_d, startCol=1, colNames=FALSE)
    addStyle(wb_xl, ws_d, cell_style(bg,"center"), rows=cur_row_d, cols=1:9, gridExpand=TRUE)
    setRowHeights(wb_xl, ws_d, rows=cur_row_d, heights=15); cur_row_d <- cur_row_d + 1
  }
  oe_d <- data.frame(Lbl=paste0("Mean ", dg), S2="—", S3="—", S4="—",
    S5=safe_mean(sub$time_reha), S6=safe_mean(sub$weight_entry),
    S7=safe_mean(sub$tagging_weight), S8=safe_mean(sub$weight_gain), S9="—")
  writeData(wb_xl, ws_d, oe_d, startRow=cur_row_d, startCol=1, colNames=FALSE)
  addStyle(wb_xl, ws_d,
    createStyle(fontName="Arial", fontSize=9, fgFill=dc,
                textDecoration="bold", halign="center", valign="center",
                border="TopBottomLeftRight", borderColour="#AABFD4"),
    rows=cur_row_d, cols=1:9, gridExpand=TRUE)
  setRowHeights(wb_xl, ws_d, rows=cur_row_d, heights=16); cur_row_d <- cur_row_d + 2
}
setColWidths(wb_xl, ws_d, cols=1:9, widths=c(12,10,10,10,9,13,13,14,10))

# ── Save Excel ─────────────────────────────────────────────────
xl_pfad <- file.path(output_ordner, "Block1_Uebersicht.xlsx")
saveWorkbook(wb_xl, xl_pfad, overwrite=TRUE)
cat("✓ Excel overview saved:", basename(xl_pfad), "\n\n")

cat("── Output files ──────────────────────────────────\n")
cat("  Plots (", length(list.files(output_ordner, "meta_.*\\.png")),
    " files, meta_01 – meta_06): ", output_ordner, "\n", sep="")
cat("  Word:  Block1_Methodenbericht.docx\n")
cat("  Excel: Block1_Uebersicht.xlsx\n")
cat("\n✓ Block 1 complete!\n")
