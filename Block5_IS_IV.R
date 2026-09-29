# ==============================================================================
# Interdaily Stability (IS) and Intradaily Variability (IV)
# Non-parametric circadian-rhythm metrics after van Someren et al. 1999
#
# Question: How well is the hedgehogs' activity rhythm synchronised with the
#           24h environmental cycle? Does this synchronisation change over
#           time after release?
#
# Method:  IS = Interdaily Stability (0-1; 1 = identical pattern every day)
#          IV = Intradaily Variability (0-inf; low = more stable, long rhythm)
#          RA = Relative Amplitude (M10/L5-based; 1 = maximum contrast)
#
# Data:    chrono_minuten24h.rds — full 24h minute data of all hedgehogs
# ==============================================================================

library(data.table)
library(ggplot2)
library(lubridate)
library(scales)
library(patchwork)

# ── Einstellungen ─────────────────────────────────────────────────────────────

# Projektwurzel — einzige Zeile die du ggf. anpassen musst
projekt_root    <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"

output_ordner   <- file.path(projekt_root, "output", "Block5_IS_IV")
dir.create(output_ordner, showWarnings = FALSE, recursive = TRUE)
chrono_rds      <- file.path(projekt_root, "output", "Block0_Pipeline", "chrono_minuten24h.rds")
meta_pfad       <- file.path(projekt_root, "data", "excel_files", "data_igel.xlsx")
epoch_min       <- 60L    # Epochenlänge in Minuten (Standard: 60 = stündlich)

# Mindest-Beobachtungstage für IS/IV-Berechnung
# Begründung: IS/IV sind nichtparametrische Maße die ausreichend viele
# Tages-Wiederholungen benötigen um stabil zu sein. Van Someren et al. (1999)
# empfehlen mindestens 7 Tage. Empirischer Vergleich (Block0_wmv_Vergleich.R)
# zeigte instabile Schätzungen bei n_tage < 7. Wert muss mit MIN_TAGE_IS_IV
# in Block8_RehabSuccess.R übereinstimmen.
min_tage_is     <- 7L     # Mindest-Beobachtungstage für IS/IV-Berechnung (van Someren 1999)

is_fenster_tage <- 7L     # Fenstergröße für gleitende IS/IV (Tage)
vis_cutoff      <- 14L    # Grau-Shading ab diesem Tag (analog GAMM)
n_min_plot      <- 5L     # Mindestanzahl Tiere für Plot-Cutoff (tage_seit_max)

cat("══════════════════════════════════════════════════════\n")
cat("  IS / IV Analyse — Nicht-parametrische Chronometrik\n")
cat("══════════════════════════════════════════════════════\n\n")

# ── Daten laden ───────────────────────────────────────────────────────────────
if (!file.exists(chrono_rds)) {
  stop("chrono_minuten24h.rds nicht gefunden. Bitte Block0_Datenpipeline.R ausführen.")
}

dt <- readRDS(chrono_rds)
setDT(dt)
cat("Geladene Minuten:", format(nrow(dt), big.mark = "'"), "\n")
cat("Igel:", dt[, uniqueN(igel)], "\n")
cat("Spalten:", paste(names(dt), collapse = ", "), "\n\n")

# ── Methodische Hinweise zur Klassifikationsvariable ─────────────────────────
# Die Aktivitätsdaten basieren auf pred_nested_loio_smoothed_wmv (tRackIT).
# Diese Variable enthält zwei Glättungsstufen:
#   1. 300s Rolling Window (Weighted Majority Voting)
#   2. ~20-minütige zusätzliche Glättungsstufe
#
# Konsequenz für IV: Die effektive zeitliche Auflösung liegt bei ~20 Minuten.
# IV misst kurzfristige Aktivitätswechsel — diese werden durch die Glättung
# artifiziell unterdrückt. IV-Werte sind daher systematisch unterschätzt.
#
# Entscheidung: IV wird beibehalten, aber NUR als relatives Vergleichsmaß
# zwischen Individuen interpretiert (alle Tiere gleich betroffen).
# Die Verwendung von loio_mwv (ungeglättet) wurde verworfen, da sie
# methodologische Inkonsistenz mit Block2–Block8 erzeugt hätte.
#
# IS ist weniger betroffen, da IS auf stündlichen Mittelwerten basiert und
# die Glättung auf dieser Aggregationsebene kaum ins Gewicht fällt.
#
# Entscheidungsgrundlage für smoothed_wmv: siehe Block0_wmv_Vergleich.R
# und Block0_Methodenbericht.docx Abschnitt 3.1.1 und 6.8.
# ─────────────────────────────────────────────────────────────────────────────

# Metadaten laden für Gruppenvergleiche
if (file.exists(meta_pfad)) {
  meta_raw <- as.data.table(readxl::read_excel(meta_pfad))
  setnames(meta_raw, names(meta_raw), tolower(trimws(gsub("\\s+", "_", names(meta_raw)))))

  meta_ok <- data.table(
    igel          = trimws(as.character(meta_raw[["individual"]])),
    sex           = meta_raw[["sex"]],
    weight_entry  = suppressWarnings(as.numeric(meta_raw[["weigth_entry"]])),
    release_date  = suppressWarnings(as.Date(meta_raw[["date_release"]]))
  )

  meta_ok[!grepl("^Igel", igel, ignore.case = TRUE),
          igel := paste0("Igel", igel)]

  meta_ok[, alter := ifelse(!is.na(weight_entry) & weight_entry < 300,
                             "Jungtier", "Adult")]

  meta_ok[, monat := month(release_date)]
  meta_ok[, saison_auswild := fcase(
    monat %in% c(3, 4, 5),   "Frühling",
    monat %in% c(6, 7, 8),   "Sommer",
    monat %in% c(9, 10, 11), "Herbst",
    default = NA_character_
  )]
  meta_ok[, c("weight_entry", "release_date", "monat") := NULL]

  cat("Metadaten geladen:", nrow(meta_ok), "Zeilen\n")
  cat("Igel-Namen:", paste(head(meta_ok$igel, 5), collapse = ", "), "\n")
  cat("Alter-Verteilung:", paste(meta_ok[, .N, alter][order(alter)]$alter,
                                  meta_ok[, .N, alter][order(alter)]$N,
                                  sep = "=", collapse = ", "), "\n\n")
} else {
  meta_ok <- NULL
  cat("⚠ Metadaten nicht gefunden — Gruppenvergleiche entfallen.\n\n")
}

# ── Hilfsfunktionen ───────────────────────────────────────────────────────────

# Aggregiert Minutendaten zu stündlichen Epochen (Anteil aktiver Minuten 0–1)
aggregiere_epochen <- function(dt_igel, epoch_min = 60L) {
  dt_igel[, epoch_h := floor(stunde * 60 / epoch_min) * (epoch_min / 60)]
  ep <- dt_igel[, .(x = mean(aktiv, na.rm = TRUE)), by = .(datum, epoch_h)]
  setorder(ep, datum, epoch_h)
  ep
}

# IS: Interdaily Stability
# Misst wie ähnlich das 24h-Profil von Tag zu Tag ist.
# IS ~ 1: jeder Tag identisch; IS ~ 0: zufällig
berechne_IS <- function(ep, min_tage = min_tage_is, epoch_min_l = epoch_min) {
  p   <- 24L / (epoch_min_l / 60L)   # Epochen pro Tag
  n   <- nrow(ep)
  # Mindestanforderung nach van Someren (1997): >= min_tage BEOBACHTUNGSTAGE.
  # NICHT min_tage*24 befüllte Zellen — VHF-Ortungen decken nie jede Stunde
  # jedes Tages ab, sonst wäre IS (v.a. im gleitenden 7-Tage-Fenster) fast
  # immer NA. Die Tage-Prüfung ist das biologisch gemeinte Kriterium.
  if (uniqueN(ep$datum) < min_tage) return(NA_real_)
  if (n < 2L) return(NA_real_)
  x_mean  <- mean(ep$x, na.rm = TRUE)
  ep[, epoch_idx := .GRP, by = epoch_h]
  xh_mean <- ep[, mean(x, na.rm = TRUE), by = epoch_idx]$V1
  zaehler <- n * sum((xh_mean - x_mean)^2, na.rm = TRUE)
  nenner  <- p * sum((ep$x  - x_mean)^2, na.rm = TRUE)
  if (nenner == 0) return(NA_real_)
  min(zaehler / nenner, 1)
}

# IV: Intradaily Variability
# Misst wie fragmentiert/unstetig der Rhythmus ist.
# IV ~ 0: lange, stabile Aktivitätsphasen; IV > 1: sehr fragmentiert
# HINWEIS: Nur als relatives Vergleichsmaß zwischen Individuen interpretieren
# (systematische Unterschätzung durch ~20-Min-Glättung in smoothed_wmv)
berechne_IV <- function(ep) {
  n      <- nrow(ep)
  if (n < 4L) return(NA_real_)
  x_mean <- mean(ep$x, na.rm = TRUE)
  nenner <- sum((ep$x - x_mean)^2, na.rm = TRUE)
  if (nenner == 0) return(NA_real_)
  diff_sq <- diff(ep$x)^2
  (n * sum(diff_sq, na.rm = TRUE)) / ((n - 1L) * nenner)
}

# RA: Relative Amplitude (M10 vs L5)
# M10 = mittlere Aktivität der aktivsten 10h; L5 = ruhigste 5h
# RA = (M10 - L5) / (M10 + L5)
berechne_RA <- function(ep, epoch_min_l = epoch_min) {
  p <- 24L / (epoch_min_l / 60L)
  if (nrow(ep) < p) return(NA_real_)
  profil <- ep[, mean(x, na.rm = TRUE), by = epoch_h][order(epoch_h)]$V1
  if (length(profil) < p) return(NA_real_)
  m10_n  <- as.integer(10L / (epoch_min_l / 60L))
  l5_n   <- as.integer(5L  / (epoch_min_l / 60L))
  profil2 <- c(profil, profil)
  m10 <- max(sapply(1:length(profil),
                    function(i) mean(profil2[i:(i + m10_n - 1L)])))
  l5  <- min(sapply(1:length(profil),
                    function(i) mean(profil2[i:(i + l5_n  - 1L)])))
  if ((m10 + l5) == 0) return(NA_real_)
  (m10 - l5) / (m10 + l5)
}

# ── Gesamtwert pro Tier ───────────────────────────────────────────────────────
cat("════════════════════════════════════\n")
cat("1. Gesamt-IS/IV pro Tier\n")
cat("════════════════════════════════════\n\n")

igel_liste <- dt[, levels(droplevels(factor(igel)))]

is_iv_gesamt <- rbindlist(lapply(igel_liste, function(ig) {
  sub <- dt[igel == ig]
  n_tage <- sub[, uniqueN(datum)]
  if (n_tage < min_tage_is) {
    return(data.table(igel = ig, n_tage = n_tage,
                      IS = NA_real_, IV = NA_real_, RA = NA_real_))
  }
  ep <- aggregiere_epochen(sub, epoch_min)
  data.table(
    igel   = ig,
    n_tage = n_tage,
    IS     = round(berechne_IS(ep), 4),
    IV     = round(berechne_IV(ep), 4),
    RA     = round(berechne_RA(ep), 4)
  )
}))

if (!is.null(meta_ok)) {
  is_iv_gesamt <- merge(is_iv_gesamt, meta_ok, by = "igel", all.x = TRUE)
}

n_mit_is    <- nrow(is_iv_gesamt[!is.na(IS)])
n_zu_kurz   <- nrow(is_iv_gesamt[is.na(IS) & !is.na(n_tage) & n_tage > 0])
tiere_kurz  <- is_iv_gesamt[is.na(IS), igel]

cat(sprintf("Ergebnisse: %d von %d Tieren mit gültigem IS/IV\n",
            n_mit_is, nrow(is_iv_gesamt)))
if (n_zu_kurz > 0)
  cat(sprintf("  Ausgeschlossen (< %d Tage): %s\n",
              min_tage_is, paste(tiere_kurz, collapse = ", ")))
cat("\n")
print(is_iv_gesamt[order(igel)])
cat("\n")

# Plot-Cutoff aus Rohdaten berechnen (vor Zeitverlaufs-Loop)
# Letzter Tag mit >= n_min_plot Tieren — konsistent mit Block2 GAMM
n_roh_tag     <- dt[!is.na(aktiv), .(n_tiere = uniqueN(igel)), by = tage_seit]
tage_seit_max <- n_roh_tag[n_tiere >= n_min_plot, max(tage_seit, na.rm = TRUE)]
cat(sprintf("Plot-Cutoff (aus Rohdaten): Tag %d (letzter Tag mit N >= %d Tieren)\n\n",
            as.integer(tage_seit_max), n_min_plot))

cat("── Deskriptive Statistik ──\n")
for (var in c("IS", "IV", "RA")) {
  vals <- is_iv_gesamt[[var]]
  vals <- vals[!is.na(vals)]
  cat(sprintf("  %-3s  N=%d  Median=%.3f  Mean=%.3f  SD=%.3f  Min=%.3f  Max=%.3f\n",
              var, length(vals), median(vals), mean(vals), sd(vals),
              min(vals), max(vals)))
}
cat("  ⚠ IV: nur relatives Vergleichsmaß — Absolutwerte durch smoothed_wmv unterschätzt\n\n")

# ── Gleitendes IS/IV über Zeit (pro Tier) ────────────────────────────────────
cat("════════════════════════════════════\n")
cat("2. Zeitverlauf IS/IV (gleitend, Fenster =", is_fenster_tage, "Tage)\n")
cat("════════════════════════════════════\n\n")

is_iv_zeit <- rbindlist(lapply(igel_liste, function(ig) {
  sub  <- dt[igel == ig & !is.na(aktiv)]
  if (sub[, uniqueN(tage_seit)] < min_tage_is) return(NULL)
  tage <- sort(unique(sub$tage_seit[sub$tage_seit <= tage_seit_max]))

  rbindlist(lapply(tage, function(t) {
    t_min <- t - floor(is_fenster_tage / 2)
    t_max <- t + floor(is_fenster_tage / 2)
    sub_f <- sub[tage_seit >= t_min & tage_seit <= t_max]
    n_t   <- sub_f[, uniqueN(datum)]
    if (n_t < min_tage_is) return(NULL)
    ep <- aggregiere_epochen(copy(sub_f), epoch_min)
    data.table(
      igel      = ig,
      tage_seit = t,
      n_tage    = n_t,
      IS        = round(berechne_IS(ep), 4),
      IV        = round(berechne_IV(ep), 4)
    )
  }), fill = TRUE)
}), fill = TRUE)

if (!is.null(meta_ok) && nrow(is_iv_zeit) > 0) {
  is_iv_zeit <- merge(is_iv_zeit, meta_ok, by = "igel", all.x = TRUE)
}

cat("Zeitverlaufs-Datenpunkte:", nrow(is_iv_zeit[!is.na(IS)]), "\n")

# Cutoff bereits oben berechnet — Zeitverlauf wurde auf tage_seit <= tage_seit_max begrenzt

# ── Gruppenvergleiche ─────────────────────────────────────────────────────────
if (!is.null(meta_ok) && "alter" %in% names(is_iv_gesamt)) {
  cat("════════════════════════════════════\n")
  cat("3. Gruppenvergleiche\n")
  cat("════════════════════════════════════\n\n")

  vergleiche <- list()

  for (var in c("IS", "IV", "RA")) {
    # Alter
    if ("alter" %in% names(is_iv_gesamt)) {
      grp <- split(is_iv_gesamt[[var]][!is.na(is_iv_gesamt[[var]])],
                   is_iv_gesamt$alter[!is.na(is_iv_gesamt[[var]])])
      if (length(grp) == 2 && all(sapply(grp, length) >= 3)) {
        wt <- wilcox.test(grp[[1]], grp[[2]])
        vergleiche[[paste0(var, "_alter")]] <- data.table(
          test   = paste0("Wilcoxon (Alter, ", var, ")"),
          W      = round(wt$statistic, 2),
          p_wert = round(wt$p.value, 4),
          signif = ifelse(wt$p.value < 0.05, "ja *", "nein")
        )
      }
    }
    # Saison
    if ("saison_auswild" %in% names(is_iv_gesamt)) {
      vals_ok <- is_iv_gesamt[!is.na(get(var)) & !is.na(saison_auswild)]
      if (vals_ok[, uniqueN(saison_auswild)] >= 2) {
        kt <- kruskal.test(as.formula(paste(var, "~ saison_auswild")),
                           data = vals_ok)
        vergleiche[[paste0(var, "_saison")]] <- data.table(
          test   = paste0("Kruskal-Wallis (Saison, ", var, ")"),
          W      = round(kt$statistic, 3),
          p_wert = round(kt$p.value, 4),
          signif = ifelse(kt$p.value < 0.05, "ja *", "nein")
        )
      }
    }
  }

  if (length(vergleiche) > 0) {
    erg_vgl <- rbindlist(vergleiche)
    print(erg_vgl)
    cat("\n")
  }
}

# Korrelation IS/IV/RA ~ Beobachtungsdauer (n_tage)
cat("── Korrelation IS/IV/RA ~ Beobachtungsdauer ──\n")
for (var in c("IS", "IV", "RA")) {
  sub_cor <- is_iv_gesamt[!is.na(get(var)) & !is.na(n_tage)]
  if (nrow(sub_cor) >= 5) {
    ct <- cor.test(sub_cor[[var]], sub_cor$n_tage, method = "spearman", exact = FALSE)
    cat(sprintf("  %-3s ~ n_tage  rho=%.3f  p=%.3f  n=%d\n",
                var, ct$estimate, ct$p.value, nrow(sub_cor)))
  }
}
cat("\n")

# ── Plots ─────────────────────────────────────────────────────────────────────
cat("════════════════════════════════════\n")
cat("4. Plots\n")
cat("════════════════════════════════════\n\n")

farben_alter <- c("Jungtier" = "#E76F51", "Adult" = "#2A9D8F")

# ── Plot 1: IS und IV Zeitverlauf (Spaghetti) ─────────────────────────────────
# Auf Tage mit >= n_min_plot Tieren beschränken (verhindert dass ein Ausreißer
# wie Igel5 die x-Achse bis Tag 220 zieht)
is_iv_zeit_plot <- is_iv_zeit[!is.na(IS) & !is.na(igel) & tage_seit <= tage_seit_max]

p_is_zeit <- ggplot(is_iv_zeit_plot,
                    aes(x = tage_seit, y = IS, group = igel, color = igel)) +
  geom_line(alpha = 0.55, linewidth = 0.6) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE,
              color = "black", fill = "grey30", alpha = 0.15, linewidth = 1.2) +
  annotate("rect",
           xmin = vis_cutoff + 0.5, xmax = max(is_iv_zeit_plot$tage_seit, na.rm = TRUE) + 0.5,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.35) +
  geom_vline(xintercept = vis_cutoff + 0.5,
             linetype = "dashed", color = "grey50", linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.01)) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  labs(title = "Interdaily Stability (IS) über Zeit",
       subtitle = paste0("IS = Synchronisation mit dem 24h-Umweltrhythmus | ",
                         "Fenster: ±", floor(is_fenster_tage/2), " Tage | ",
                         "Grau = N reduziert (ab Tag ", vis_cutoff + 1, ")"),
       x = "Tage nach Auswilderung", y = "IS (0–1)") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"))

p_iv_zeit <- ggplot(is_iv_zeit_plot,
                    aes(x = tage_seit, y = IV, group = igel, color = igel)) +
  geom_line(alpha = 0.55, linewidth = 0.6) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE,
              color = "black", fill = "grey30", alpha = 0.15, linewidth = 1.2) +
  annotate("rect",
           xmin = vis_cutoff + 0.5, xmax = max(is_iv_zeit_plot$tage_seit, na.rm = TRUE) + 0.5,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.35) +
  geom_vline(xintercept = vis_cutoff + 0.5,
             linetype = "dashed", color = "grey50", linewidth = 0.8) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  labs(title = "Intradaily Variability (IV) über Zeit — nur relativ interpretieren",
       subtitle = paste0("IV = Fragmentierung des Rhythmus (niedriger = stabiler) | ",
                         "⚠ Absolutwerte durch ~20-Min-Glättung unterschätzt"),
       x = "Tage nach Auswilderung", y = "IV (relativ)") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#C55A11"))

png(file.path(output_ordner, "is_iv_zeitverlauf.png"),
    width = 1600, height = 1200, res = 150)
print(p_is_zeit / p_iv_zeit)
dev.off()
cat("✓ is_iv_zeitverlauf.png\n")

# ── Plot 2: IS/IV/RA Übersicht pro Tier (Barplot) ─────────────────────────────
is_iv_long <- melt(is_iv_gesamt[!is.na(IS)],
                   id.vars    = c("igel", "n_tage",
                                  intersect(c("alter", "sex", "saison_auswild"),
                                            names(is_iv_gesamt))),
                   measure.vars = c("IS", "IV", "RA"),
                   variable.name = "metrik",
                   value.name    = "wert")

is_iv_long <- is_iv_long[!is.na(wert)]
is_iv_long[, igel_n := paste0(igel, " (n=", n_tage, ")")]

color_var <- if ("alter" %in% names(is_iv_long) &&
                  any(!is.na(is_iv_long$alter))) "alter" else NULL

p_uebersicht <- ggplot(is_iv_long,
                       aes(x = reorder(igel, wert),
                           y = wert,
                           fill = if (!is.null(color_var)) get(color_var) else igel)) +
  geom_col(alpha = 0.85, width = 0.75) +
  facet_wrap(~ metrik, scales = "free_y", ncol = 1,
             labeller = labeller(metrik = c(
               IS = "IS — Interdaily Stability (Synchronisation)",
               IV = "IV — Intradaily Variability (relativ; Glättungsartefakt beachten)",
               RA = "RA — Relative Amplitude (Tag/Nacht-Kontrast)"
             ))) +
  coord_flip() +
  { if (!is.null(color_var)) scale_fill_manual(values = farben_alter, name = "Alter")
    else scale_fill_viridis_d(option = "turbo", guide = "none") } +
  labs(title = "IS / IV / RA je Igel — Gesamt-Rhythmusqualität",
       subtitle = paste0("Berechnet auf allen verfügbaren 24h-Epochen | ",
                         "Epochenlänge: ", epoch_min, " Minuten | ",
                         "IV nur relativ zwischen Individuen vergleichen"),
       x = NULL, y = "Wert") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        strip.text    = element_text(face = "bold", size = 9))

png(file.path(output_ordner, "is_iv_uebersicht.png"),
    width = 1400, height = 1600, res = 150)
print(p_uebersicht)
dev.off()
cat("✓ is_iv_uebersicht.png\n")

# ── Plot 3: Gruppenvergleich IS nach Alter und Saison ────────────────────────
if (!is.null(meta_ok)) {
  plot_list <- list()

  for (grp_var in intersect(c("alter", "saison_auswild"), names(is_iv_gesamt))) {
    sub <- is_iv_gesamt[!is.na(IS) & !is.na(get(grp_var))]
    if (nrow(sub) < 4) next

    medians <- sub[, .(med = round(median(IS, na.rm = TRUE), 3),
                       n   = .N), by = grp_var]
    setnames(medians, grp_var, "grp")

    p <- ggplot(sub, aes(x = .data[[grp_var]], y = IS, fill = .data[[grp_var]])) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.4) +
      geom_jitter(aes(color = .data[[grp_var]]), width = 0.1, size = 2.5, alpha = 0.8) +
      geom_text(data = medians,
                aes(x = grp, y = med,
                    label = paste0("Md=", med, "\nn=", n)),
                vjust = -0.5, size = 3, fontface = "bold",
                inherit.aes = FALSE) +
      { if (grp_var == "alter")
          scale_fill_manual(values = farben_alter, guide = "none")
        else scale_fill_viridis_d(guide = "none") } +
      { if (grp_var == "alter")
          scale_color_manual(values = farben_alter, guide = "none")
        else scale_color_viridis_d(guide = "none") } +
      labs(title = paste0("IS nach ", grp_var),
           x = NULL, y = "IS (0–1)") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))

    plot_list[[grp_var]] <- p
  }

  if (length(plot_list) >= 1) {
    png(file.path(output_ordner, "is_iv_gruppen.png"),
        width = 1400, height = 700, res = 150)
    if (length(plot_list) == 2) {
      print(plot_list[[1]] | plot_list[[2]])
    } else {
      print(plot_list[[1]])
    }
    dev.off()
    cat("✓ is_iv_gruppen.png\n")
  }
}

# ── Ergebnisse speichern ──────────────────────────────────────────────────────
cat("\n")
saveRDS(list(gesamt = is_iv_gesamt, zeitverlauf = is_iv_zeit),
        file.path(output_ordner, "is_iv_ergebnisse.rds"))
cat("✓ is_iv_ergebnisse.rds gespeichert\n")

sink(file.path(output_ordner, "is_iv_kennzahlen.txt"))
cat("IS / IV / RA Ergebnisse\n")
cat("══════════════════════════════════\n\n")
cat("── Gesamt-Kennzahlen ──\n")
for (var in c("IS", "IV", "RA")) {
  vals <- is_iv_gesamt[[var]]
  vals <- vals[!is.na(vals)]
  cat(sprintf("%-3s  N=%d  Median=%.3f  Mean=%.3f  SD=%.3f  Min=%.3f  Max=%.3f\n",
              var, length(vals), median(vals), mean(vals), sd(vals),
              min(vals), max(vals)))
}
cat("\nIV-Hinweis: Absolutwerte durch ~20-Min-Glättung (smoothed_wmv) unterschätzt.\n")
cat("IV nur als relatives Vergleichsmaß zwischen Individuen interpretieren.\n")
cat("\n── Werte pro Tier ──\n")
print(is_iv_gesamt[, c("igel", "n_tage", "IS", "IV", "RA",
                        intersect(c("alter", "sex", "saison_auswild"),
                                  names(is_iv_gesamt))), with = FALSE])
sink()
cat("✓ is_iv_kennzahlen.txt gespeichert\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Excel-Übersicht
# ══════════════════════════════════════════════════════════════════════════════
cat("════════════════════════════════════\n")
cat("5. Excel-Übersicht\n")
cat("════════════════════════════════════\n\n")

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)

  wb <- createWorkbook()

  # Stile
  hs       <- createStyle(fontName = "Arial", fontSize = 11, fontColour = "white",
                           fgFill = "#1F4E79", textDecoration = "bold",
                           halign = "center", valign = "center",
                           border = "TopBottomLeftRight", borderColour = "#BFBFBF")
  cs       <- createStyle(fontName = "Arial", fontSize = 10,
                           border = "TopBottomLeftRight", borderColour = "#D9D9D9")
  cs_warn  <- createStyle(fontName = "Arial", fontSize = 10, fgFill = "#FCE4D6",
                           border = "TopBottomLeftRight", borderColour = "#D9D9D9")
  cs_sig   <- createStyle(fontName = "Arial", fontSize = 10, fgFill = "#E2EFDA",
                           border = "TopBottomLeftRight", borderColour = "#D9D9D9")
  cs_na    <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "#999999",
                           border = "TopBottomLeftRight", borderColour = "#D9D9D9")
  cs_note  <- createStyle(fontName = "Arial", fontSize = 9, fontColour = "#595959",
                           wrapText = TRUE)
  cs_title <- createStyle(fontName = "Arial", fontSize = 13, textDecoration = "bold",
                           fontColour = "#1F4E79")
  cs_head2 <- createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                           fontColour = "#2E75B6")

  addWorksheet(wb, "Übersicht")
  addWorksheet(wb, "Individuen")
  addWorksheet(wb, "Gruppenvergleiche")
  addWorksheet(wb, "Zeitverlauf")

  # ── Sheet 1: Übersicht ──────────────────────────────────────────────────────
  sh1 <- "Übersicht"
  writeData(wb, sh1, "IS / IV / RA — Block5 Ergebnisübersicht", startRow = 1)
  addStyle(wb, sh1, cs_title, rows = 1, cols = 1)
  writeData(wb, sh1, paste0("Erstellt: ", format(Sys.time(), "%d.%m.%Y %H:%M")), startRow = 2)
  writeData(wb, sh1,
    paste0("Epochenlänge: ", epoch_min, " Min | Min. Beobachtungstage: ",
           min_tage_is, " | Zeitverlauf-Fenster: ±", floor(is_fenster_tage/2), " Tage"),
    startRow = 3)

  writeData(wb, sh1, "Deskriptive Statistik", startRow = 5)
  addStyle(wb, sh1, cs_head2, rows = 5, cols = 1)

  uebers_df <- data.frame(
    Metrik       = c("IS", "IV (relativ)", "RA"),
    Beschreibung = c(
      "Interdaily Stability: Synchronisation mit 24h-Rhythmus (0=zufällig, 1=identisch)",
      "Intradaily Variability: Fragmentierung — NUR relativ interpretieren (Glättungsartefakt)",
      "Relative Amplitude: Tag/Nacht-Kontrast M10/L5 (0=kein Kontrast, 1=maximal)"
    ),
    N       = sapply(c("IS","IV","RA"), function(v) sum(!is.na(is_iv_gesamt[[v]]))),
    Median  = sapply(c("IS","IV","RA"), function(v) round(median(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    Mittelw = sapply(c("IS","IV","RA"), function(v) round(mean(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    SD      = sapply(c("IS","IV","RA"), function(v) round(sd(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    Min     = sapply(c("IS","IV","RA"), function(v) round(min(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    Max     = sapply(c("IS","IV","RA"), function(v) round(max(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    stringsAsFactors = FALSE
  )
  writeData(wb, sh1, uebers_df, startRow = 6, headerStyle = hs)
  addStyle(wb, sh1, cs, rows = 7:9, cols = 1:8, gridExpand = TRUE)
  addStyle(wb, sh1, cs_warn, rows = 8, cols = 1:8, gridExpand = TRUE)   # IV orange

  writeData(wb, sh1,
    paste0("Hinweis IV: Die Variable smoothed_wmv enthält eine ~20-Min-Glättung, die ",
           "kurzfristige Aktivitätswechsel artifiziell unterdrückt. IV-Absolutwerte sind ",
           "systematisch unterschätzt und nur als relatives Vergleichsmaß zwischen Individuen ",
           "interpretierbar. Verwendung von loio_mwv wurde verworfen (Inkonsistenz mit Block2-8)."),
    startRow = 11)
  addStyle(wb, sh1, cs_note, rows = 11, cols = 1)
  mergeCells(wb, sh1, cols = 1:8, rows = 11)
  setColWidths(wb, sh1, cols = 1:8, widths = c(15, 58, 5, 8, 8, 8, 8, 8))
  setRowHeights(wb, sh1, rows = 11, heights = 45)

  # ── Sheet 2: Individuen ─────────────────────────────────────────────────────
  sh2 <- "Individuen"
  writeData(wb, sh2, "IS / IV / RA pro Individuum", startRow = 1)
  addStyle(wb, sh2, cs_title, rows = 1, cols = 1)
  writeData(wb, sh2,
    paste0("N = ", nrow(is_iv_gesamt), " Tiere | ",
           sum(!is.na(is_iv_gesamt$IS)), " mit gültigem IS/IV (>= ", min_tage_is, " Tage)"),
    startRow = 2)

  cols_sh2 <- intersect(c("igel","n_tage","IS","IV","RA","alter","sex","saison_auswild"),
                        names(is_iv_gesamt))
  ind_df <- as.data.frame(is_iv_gesamt[order(igel), ..cols_sh2])
  writeData(wb, sh2, ind_df, startRow = 4, headerStyle = hs)
  n_ind <- nrow(ind_df)
  addStyle(wb, sh2, cs, rows = 5:(4 + n_ind), cols = 1:ncol(ind_df), gridExpand = TRUE)

  # NA-Zeilen grau
  na_rows <- which(is.na(ind_df$IS))
  if (length(na_rows) > 0)
    addStyle(wb, sh2, cs_na, rows = 4 + na_rows, cols = 1:ncol(ind_df), gridExpand = TRUE)

  # Zahlenformat
  num_cols <- which(names(ind_df) %in% c("IS","IV","RA"))
  if (length(num_cols) > 0) {
    num_style <- createStyle(numFmt = "0.000")
    addStyle(wb, sh2, num_style, rows = 5:(4 + n_ind), cols = num_cols, gridExpand = TRUE)
  }
  setColWidths(wb, sh2, cols = 1:ncol(ind_df), widths = c(12,10,8,8,8,12,8,14))

  # ── Sheet 3: Gruppenvergleiche ──────────────────────────────────────────────
  sh3 <- "Gruppenvergleiche"
  writeData(wb, sh3, "Statistische Tests — IS / IV / RA", startRow = 1)
  addStyle(wb, sh3, cs_title, rows = 1, cols = 1)
  writeData(wb, sh3, "Grün = signifikant (p < 0,05)", startRow = 2)

  test_rows <- list()
  for (var in c("IS", "IV", "RA")) {
    # Wilcoxon Alter
    if ("alter" %in% names(is_iv_gesamt)) {
      vals <- is_iv_gesamt[!is.na(get(var))]
      grp  <- split(vals[[var]], vals$alter)
      if (length(grp) == 2 && all(sapply(grp, length) >= 3)) {
        wt <- wilcox.test(grp[[1]], grp[[2]])
        test_rows[[length(test_rows)+1]] <- data.frame(
          Metrik = var, Faktor = "Alter (Adult vs. Jungtier)",
          Test = "Wilcoxon", Statistik = round(wt$statistic, 2),
          p_Wert = round(wt$p.value, 4),
          Signifikant = ifelse(wt$p.value < 0.05, "ja *", "nein"),
          stringsAsFactors = FALSE)
      }
    }
    # Kruskal-Wallis Saison
    if ("saison_auswild" %in% names(is_iv_gesamt)) {
      vals <- is_iv_gesamt[!is.na(get(var)) & !is.na(saison_auswild)]
      if (vals[, uniqueN(saison_auswild)] >= 2) {
        kt <- kruskal.test(as.formula(paste(var, "~ saison_auswild")), data = vals)
        test_rows[[length(test_rows)+1]] <- data.frame(
          Metrik = var, Faktor = "Saison der Auswilderung",
          Test = "Kruskal-Wallis", Statistik = round(kt$statistic, 3),
          p_Wert = round(kt$p.value, 4),
          Signifikant = ifelse(kt$p.value < 0.05, "ja *", "nein"),
          stringsAsFactors = FALSE)
      }
    }
    # Spearman IS/IV/RA ~ n_tage
    sub_cor <- is_iv_gesamt[!is.na(get(var)) & !is.na(n_tage)]
    if (nrow(sub_cor) >= 5) {
      ct <- cor.test(sub_cor[[var]], sub_cor$n_tage, method = "spearman", exact = FALSE)
      test_rows[[length(test_rows)+1]] <- data.frame(
        Metrik = var, Faktor = "Beobachtungsdauer (n_tage)",
        Test = paste0("Spearman"), Statistik = round(ct$estimate, 3),
        p_Wert = round(ct$p.value, 4),
        Signifikant = ifelse(ct$p.value < 0.05, "ja *", "nein"),
        stringsAsFactors = FALSE)
    }
  }

  if (length(test_rows) > 0) {
    test_df <- do.call(rbind, test_rows)
    writeData(wb, sh3, test_df, startRow = 4, headerStyle = hs)
    n_tests <- nrow(test_df)
    addStyle(wb, sh3, cs, rows = 5:(4 + n_tests), cols = 1:6, gridExpand = TRUE)
    sig_rows_xl <- which(test_df$Signifikant == "ja *")
    if (length(sig_rows_xl) > 0)
      addStyle(wb, sh3, cs_sig, rows = 4 + sig_rows_xl, cols = 1:6, gridExpand = TRUE)
    setColWidths(wb, sh3, cols = 1:6, widths = c(8, 26, 15, 12, 10, 14))
  } else {
    writeData(wb, sh3, "Keine Gruppenvergleiche möglich (Metadaten fehlen).", startRow = 4)
    test_df <- NULL
  }

  # ── Sheet 4: Zeitverlauf ────────────────────────────────────────────────────
  sh4 <- "Zeitverlauf"
  writeData(wb, sh4, paste0("Gleitender IS/IV-Zeitverlauf (7-Tage-Fenster)"), startRow = 1)
  addStyle(wb, sh4, cs_title, rows = 1, cols = 1)
  writeData(wb, sh4,
    paste0("N Datenpunkte: ", nrow(is_iv_zeit[!is.na(IS)]), " | IV nur relativ interpretieren"),
    startRow = 2)

  zeit_cols <- intersect(c("igel","tage_seit","n_tage","IS","IV","alter","saison_auswild"),
                         names(is_iv_zeit))
  zeit_df <- as.data.frame(is_iv_zeit[!is.na(IS), ..zeit_cols][order(igel, tage_seit)])
  writeData(wb, sh4, zeit_df, startRow = 4, headerStyle = hs)
  if (nrow(zeit_df) > 0) {
    addStyle(wb, sh4, cs, rows = 5:(4 + nrow(zeit_df)), cols = 1:ncol(zeit_df),
             gridExpand = TRUE)
    num_cols_z <- which(names(zeit_df) %in% c("IS","IV"))
    if (length(num_cols_z) > 0)
      addStyle(wb, sh4, createStyle(numFmt = "0.0000"),
               rows = 5:(4 + nrow(zeit_df)), cols = num_cols_z, gridExpand = TRUE)
  }
  setColWidths(wb, sh4, cols = 1:ncol(zeit_df), widths = c(12,10,10,10,10,12,14))

  excel_pfad <- file.path(output_ordner, "Block5_IS_IV_Uebersicht.xlsx")
  saveWorkbook(wb, excel_pfad, overwrite = TRUE)
  cat("✓ Block5_IS_IV_Uebersicht.xlsx gespeichert\n")
} else {
  cat("  openxlsx nicht verfügbar — Excel-Export übersprungen\n")
  test_df <- NULL
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Methodenbericht (Word)
# ══════════════════════════════════════════════════════════════════════════════
cat("════════════════════════════════════\n")
cat("6. Methodenbericht (Word)\n")
cat("════════════════════════════════════\n\n")

if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")

if (requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)) {
  library(officer)
  library(flextable)

  make_ft <- function(df, hl_rows = NULL) {
    brd <- fp_border(color = "#BFBFBF", width = 0.5)
    ft <- flextable(df) |>
      font(fontname = "Arial", part = "all") |>
      fontsize(size = 10, part = "all") |>
      bold(part = "header") |>
      color(part = "header", color = "white") |>
      bg(part = "header", bg = "#1F4E79") |>
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

  doc <- doc |>
    body_add_par("Block 5: IS / IV / RA — Chronobiologische Rhythmusqualität",
                 style = "heading 1") |>
    body_add_par("Igelbesenderung Niedersachsen — TiHo Hannover | Natalie Steiner",
                 style = "Normal") |>
    body_add_par(paste0("Erstellt: ", format(Sys.Date(), "%d. %B %Y")), style = "Normal") |>
    body_add_par("", style = "Normal")

  # 1. Fragestellung
  doc <- doc |>
    body_add_par("1. Fragestellung", style = "heading 2") |>
    body_add_par(paste0(
      "Wie gut synchronisiert ist der zirkadiane Aktivitätsrhythmus ausgewilderter Igel mit ",
      "dem 24-Stunden-Umweltrhythmus? Verändert sich diese Synchronisation im Verlauf der Zeit ",
      "nach Auswilderung? Unterscheidet er sich zwischen Altersgruppen oder Auswilderungssaisons? ",
      "Dazu werden drei nicht-parametrische Metriken nach van Someren et al. (1997) berechnet: ",
      "Interdaily Stability (IS), Intradaily Variability (IV) und Relative Amplitude (RA)."
    ), style = "Normal") |>
    body_add_par("", style = "Normal")

  # 2. Methodik
  doc <- doc |>
    body_add_par("2. Methodik", style = "heading 2") |>
    body_add_par("2.1 Datenbasis", style = "heading 3") |>
    body_add_par(paste0(
      "Datenquelle: chrono_minuten24h.rds (Block0_Datenpipeline.R). ",
      "Aktivitätsklassifikation: pred_nested_loio_smoothed_wmv (tRackIT). ",
      "Epochenlänge: ", epoch_min, " Minuten (stündliche Aggregation auf Anteil aktiver Minuten). ",
      "Mindestanzahl Beobachtungstage: ", min_tage_is, " (van Someren et al. 1997). ",
      "Gleitendes Fenster für den Zeitverlauf: ±", floor(is_fenster_tage/2), " Tage (= ",
      is_fenster_tage, " Tage gesamt). Visualisierungs-Cutoff: Tag ", vis_cutoff, "."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("2.2 Metriken", style = "heading 3") |>
    body_add_par(paste0(
      "IS (Interdaily Stability) misst die Ähnlichkeit des 24h-Aktivitätsprofils ",
      "zwischen aufeinanderfolgenden Tagen. IS = 1: identisches Muster jeden Tag; IS = 0: ",
      "zufällig. Berechnung: Verhältnis der Varianz des mittleren Stundenprofils zur ",
      "Gesamtvarianz aller Stundenwerte. IS ist wenig anfällig für Glättungsartefakte, ",
      "da es auf stündlichen Mittelwerten basiert."
    ), style = "Normal") |>
    body_add_par(paste0(
      "IV (Intradaily Variability) misst die Fragmentierung des Rhythmus — also wie häufig ",
      "und abrupt Aktivitätswechsel auftreten. Niedrige IV = lange, stabile Phasen; ",
      "hohe IV = viele kurze Wechsel. Berechnung: Verhältnis der quadrierten ",
      "Stunden-zu-Stunden-Differenzen zur Gesamtvarianz."
    ), style = "Normal") |>
    body_add_par(paste0(
      "RA (Relative Amplitude) quantifiziert den Kontrast zwischen der aktivsten ",
      "10-Stunden-Phase (M10) und der ruhigsten 5-Stunden-Phase (L5): ",
      "RA = (M10 - L5) / (M10 + L5). RA = 1: maximaler Tag/Nacht-Kontrast; ",
      "RA = 0: kein Kontrast."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("2.3 Methodische Einschränkung: IV-Interpretation", style = "heading 3") |>
    body_add_par(paste0(
      "Die verwendete Aktivitätsvariable (smoothed_wmv) enthält neben dem 300-Sekunden-",
      "Rolling-Window-WMV eine zusätzliche ~20-minütige Glättungsstufe. Diese unterdrückt ",
      "kurzfristige Aktivitätswechsel artifiziell und führt zu einer systematischen ",
      "Unterschätzung der IV-Absolutwerte."
    ), style = "Normal") |>
    body_add_par(paste0(
      "Alternative: Die ungeglättete Variable loio_mwv wurde als Datenbasis für IV erwogen, ",
      "aber verworfen — sie erzeugt methodologische Inkonsistenz mit allen anderen ",
      "Analyseblöcken (Block2–Block8), die ausnahmslos auf smoothed_wmv basieren. ",
      "Da die Glättung alle Individuen in gleicher Weise betrifft, bleiben IV-Werte als ",
      "relatives Vergleichsmaß zwischen Individuen interpretierbar. In Publikationen sollte ",
      "dieser Vorbehalt explizit genannt werden."
    ), style = "Normal") |>
    body_add_par("", style = "Normal")

  # 3. Ergebnisse
  doc <- doc |>
    body_add_par("3. Ergebnisse", style = "heading 2") |>
    body_add_par("3.1 Deskriptive Kennzahlen", style = "heading 3")

  desc_df <- data.frame(
    Metrik  = c("IS", "IV (relativ)", "RA"),
    N       = sapply(c("IS","IV","RA"), function(v) sum(!is.na(is_iv_gesamt[[v]]))),
    Median  = sapply(c("IS","IV","RA"), function(v) round(median(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    Mittelw = sapply(c("IS","IV","RA"), function(v) round(mean(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    SD      = sapply(c("IS","IV","RA"), function(v) round(sd(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    Min     = sapply(c("IS","IV","RA"), function(v) round(min(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    Max     = sapply(c("IS","IV","RA"), function(v) round(max(is_iv_gesamt[[v]], na.rm=TRUE), 3)),
    stringsAsFactors = FALSE
  )
  doc <- doc |>
    body_add_flextable(make_ft(desc_df)) |>
    body_add_par("", style = "Normal")

  # Individuen-Tabelle
  doc <- doc |>
    body_add_par("3.2 Werte pro Individuum", style = "heading 3")

  cols_ind <- intersect(c("igel","n_tage","IS","IV","RA","alter","sex","saison_auswild"),
                        names(is_iv_gesamt))
  ind_tab  <- as.data.frame(is_iv_gesamt[order(igel), ..cols_ind])
  names(ind_tab)[names(ind_tab) == "n_tage"]        <- "N Tage"
  names(ind_tab)[names(ind_tab) == "saison_auswild"] <- "Saison"

  doc <- doc |>
    body_add_flextable(make_ft(ind_tab)) |>
    body_add_par("", style = "Normal")

  # 4. Statistische Tests
  doc <- doc |>
    body_add_par("4. Statistische Tests", style = "heading 2") |>
    body_add_par(paste0(
      "Gruppenvergleiche: Wilcoxon-Rangsummentest (Alter: Adult vs. Jungtier) und ",
      "Kruskal-Wallis-Test (Saison der Auswilderung). ",
      "Zusammenhangsanalyse: Spearman-Rangkorrelation (IS/IV/RA ~ Beobachtungsdauer). ",
      "Signifikanzniveau: α = 0,05."
    ), style = "Normal") |>
    body_add_par("", style = "Normal")

  if (exists("test_df") && !is.null(test_df) && nrow(test_df) > 0) {
    sig_rows_w <- which(test_df$Signifikant == "ja *")
    doc <- doc |>
      body_add_flextable(make_ft(test_df,
                                  hl_rows = if (length(sig_rows_w) > 0) sig_rows_w else NULL)) |>
      body_add_par("", style = "Normal")
  } else {
    doc <- doc |>
      body_add_par("Keine Gruppenvergleiche möglich (Metadaten nicht verfügbar).",
                   style = "Normal") |>
      body_add_par("", style = "Normal")
  }

  # 5. Interpretation
  doc <- doc |>
    body_add_par("5. Interpretation", style = "heading 2") |>
    body_add_par(paste0(
      "IS-Werte > 0,5 sprechen für eine gute zirkadiane Synchronisation der Igel mit dem ",
      "Umweltrhythmus nach der Auswilderung. Hohe RA-Werte (nahe 1,0) bestätigen einen ",
      "deutlichen Tag/Nacht-Kontrast der Aktivität, der für nachtaktive Arten erwartet wird. ",
      "Der gleitende IS-Zeitverlauf erlaubt Rückschlüsse darauf, ob die circadiane ",
      "Organisation sich in den Wochen nach der Auswilderung stabilisiert oder verändert."
    ), style = "Normal") |>
    body_add_par(paste0(
      "IV-Werte sollten ausschließlich relativ zwischen Individuen verglichen werden. ",
      "Tiere mit höherer IV zeigen — innerhalb der verfügbaren zeitlichen Auflösung — ",
      "fragmentiertere Aktivitätsmuster, was auf schlechtere Kondition, höhere ",
      "Störungsexposition oder Orientierungsprobleme hinweisen kann."
    ), style = "Normal") |>
    body_add_par("", style = "Normal") |>
    body_add_par("6. Referenzen", style = "heading 2") |>
    body_add_par(paste0(
      "van Someren, E.J.W., Lijzenga, C., Mirmiran, M., Swaab, D.F. (1997). Long-term fitness ",
      "training improves sleep quality of elderly women. Journal of Sleep Research, 6(4), 223–229."
    ), style = "Normal") |>
    body_add_par(paste0(
      "van Someren, E.J.W., Swaab, D.F., Colenda, C.C., Cohen, W., McCall, W.V., Rosenquist, P.B. ",
      "(1999). Bright light therapy: improved sensitivity to its effects on rest-activity rhythms in ",
      "Alzheimer patients by application of nonparametric methods. Chronobiology International, ",
      "16(4), 505–518."
    ), style = "Normal")

  bericht_pfad <- file.path(output_ordner, "Block5_IS_IV_Methodenbericht.docx")
  print(doc, target = bericht_pfad)
  cat("✓ Block5_IS_IV_Methodenbericht.docx gespeichert\n")
} else {
  cat("  officer/flextable nicht verfügbar — Word-Bericht übersprungen\n")
}

# ── Abschluss ─────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════\n")
cat("Block5 IS/IV Analyse abgeschlossen!\n")
cat("Erzeugte Dateien:\n")
cat("  is_iv_zeitverlauf.png          — IS + IV Spaghetti über Zeit\n")
cat("  is_iv_uebersicht.png           — IS/IV/RA Barplot je Tier\n")
cat("  is_iv_gruppen.png              — Gruppenvergleiche (Alter, Saison)\n")
cat("  is_iv_ergebnisse.rds           — Ergebnisobjekte (gesamt + zeitverlauf)\n")
cat("  is_iv_kennzahlen.txt           — Kennzahlen-Textdatei\n")
cat("  Block5_IS_IV_Uebersicht.xlsx   — Excel-Übersicht (4 Sheets)\n")
cat("  Block5_IS_IV_Methodenbericht.docx — Word-Methodenbericht\n")
cat("══════════════════════════════════\n")
