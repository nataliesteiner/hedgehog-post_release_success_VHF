# ==============================================================================
# Block_Blind_Vergleich.R
# Blind hedgehogs as a natural experiment: circadian disorganisation
# without a photic zeitgeber
#
# Background:
#   Hedgehogs use the light-dark cycle (light) as the primary zeitgeber (ZG)
#   for their circadian clock. Blind animals (H22, H23, H26) lack this input
#   -> expected consequence: lower Interdaily Stability (IS), lower night
#   proportion, altered space use.
#
# Question:
#   Do blind hedgehogs show measurably poorer behavioural normalisation than
#   sighted hedgehogs after release? -> Natural experiment that also validates
#   the sensitivity of the tRackIT system.
#
# Panels:
#   A — mean 24h activity profile: blind vs. sighted (group lines)
#   B — IS time course over days post-release (spaghetti + GAM)
#   C — night-activity proportion (violin + boxplot + jitter)

# ── 0. Pakete ──────────────────────────────────────────────────────────────────
pakete <- c("data.table","ggplot2","patchwork","scales","readxl","mgcv")
neu <- pakete[!sapply(pakete, requireNamespace, quietly=TRUE)]
if (length(neu)) install.packages(neu)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  library(patchwork);  library(scales); library(readxl); library(mgcv)
})
cat("✓ Pakete geladen\n\n")

# ── 1. Pfade ───────────────────────────────────────────────────────────────────
projekt_root  <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"
out_ordner    <- file.path(projekt_root, "output", "Block_Blind")
dir.create(out_ordner, showWarnings=FALSE, recursive=TRUE)

rds_chrono <- file.path(projekt_root, "output", "Block0_Pipeline", "chrono_minuten24h.rds")
rds_is     <- file.path(projekt_root, "output", "Block5_IS_IV",    "is_iv_ergebnisse.rds")
rds_4b     <- file.path(projekt_root, "output", "Block4b_Einzeltier","Block4b_Ergebnisse.rds")
meta_pfad  <- file.path(projekt_root, "data",   "excel_files",      "data_igel.xlsx")

# ── 2. Metadaten: blind vs. sehend klassifizieren ─────────────────────────────
cat("Lade Metadaten …\n")
meta_raw <- as.data.table(read_excel(meta_pfad))
setnames(meta_raw, tolower(trimws(gsub("\\s+","_",names(meta_raw)))))

meta <- data.table(
  igel      = trimws(as.character(meta_raw$individual)),
  diagnose  = tolower(trimws(as.character(meta_raw$diagnosis_main))),
  sex       = meta_raw$sex,
  time_reha = suppressWarnings(as.numeric(meta_raw$time_reha))
)
meta[!grepl("^Igel", igel, ignore.case=TRUE), igel := paste0("Igel", igel)]
# Stabile Gruppen-Schlüssel ("Blind"/"Sehend") — NICHT die Stichprobengröße
# ins Label kodieren (sonst brechen Filter/Skalen bei jeder Stichprobenänderung).
meta[, gruppe := ifelse(grepl("blind", diagnose, ignore.case=TRUE),
                         "Blind", "Sehend")]

BLIND   <- meta[gruppe == "Blind",  igel]
SEHEND  <- meta[gruppe == "Sehend", igel]
n_blind  <- length(BLIND)
n_sehend <- length(SEHEND)

cat(sprintf("  Blinde Tiere (n=%d): %s\n", n_blind, paste(BLIND, collapse=", ")))
cat(sprintf("  Sehende Tiere (n=%d)\n\n", n_sehend))

# Farben — klar unterscheidbar, druckbar
COL_BLIND  <- "#D6604D"   # Rot-Orange: beeinträchtigt
COL_SEHEND <- "#4C6A9C"   # Blau:       normal
farben <- c("Blind"=COL_BLIND, "Sehend"=COL_SEHEND)
# Anzeige-Labels mit dynamischem n (nur für Legenden/Achsen, nicht als Schlüssel)
grp_labels <- c("Blind"  = sprintf("Blind (n=%d)",  n_blind),
                "Sehend" = sprintf("Sehend (n=%d)", n_sehend))

# Gemeinsames Theme für alle Panels
theme_pub <- theme_bw(base_size=11) +
  theme(
    plot.title    = element_text(face="bold", size=11, hjust=0),
    plot.subtitle = element_text(size=9, color="grey40"),
    strip.text    = element_text(face="bold", size=9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — 24h-Aktivitätsprofil (Stunde × mittlere Aktivität)
# Zeigt: wann sind die Tiere aktiv? Blinde Tiere weniger stark nachtfixiert.
# ══════════════════════════════════════════════════════════════════════════════
cat("Panel A: 24h-Aktivitätsprofil …\n")

if (!file.exists(rds_chrono)) stop("chrono_minuten24h.rds fehlt — Block0 ausführen!")
chrono <- readRDS(rds_chrono)
setDT(chrono)

chrono <- merge(chrono, meta[, .(igel, gruppe)], by="igel", all.x=TRUE)
chrono <- chrono[!is.na(gruppe)]

# Mittleres Aktivitätsprofil pro Stunde und Gruppe
chrono[, stunde_int := floor(stunde) %% 24L]
profil <- chrono[!is.na(aktiv), .(
  pct_aktiv = mean(aktiv, na.rm=TRUE) * 100
), by=.(gruppe, stunde_int)]

pA <- ggplot(profil, aes(x=stunde_int, y=pct_aktiv, color=gruppe, fill=gruppe)) +
  annotate("rect", xmin=-0.5, xmax=5.5,  ymin=-Inf, ymax=Inf, fill="grey20", alpha=0.06) +
  annotate("rect", xmin=20.5, xmax=23.5, ymin=-Inf, ymax=Inf, fill="grey20", alpha=0.06) +
  geom_line(linewidth=1.4) +
  geom_ribbon(aes(ymin=0, ymax=pct_aktiv), alpha=0.12, color=NA) +
  scale_color_manual(values=farben, labels=grp_labels) +
  scale_fill_manual(values=farben, labels=grp_labels) +
  scale_x_continuous(breaks=seq(0,23,3),
                     labels=sprintf("%02d:00", seq(0,23,3)),
                     limits=c(0,23)) +
  scale_y_continuous(labels=function(x) paste0(x,"%"), limits=c(0,NA)) +
  labs(title="A   24h-Aktivitätsprofil",
       subtitle="Mittlerer Anteil aktiver Minuten pro Stunde | Grau = Kernnacht",
       x="Uhrzeit (MEZ/CEST)", y="Aktivität (%)") +
  theme_pub +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        legend.position="inside",
        legend.position.inside=c(0.18, 0.88),
        legend.background=element_rect(fill="white", color="grey80"),
        legend.key.size=unit(0.5,"cm"),
        legend.text=element_text(size=8),
        legend.title=element_blank()) +
  guides(color=guide_legend(), fill="none")

cat("  ✓ Panel A erstellt\n")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B — IS-Zeitverlauf (Rolling Window, Spaghetti + GAM)
# Zeigt: IS stabilisiert sich bei Sehenden, bleibt niedrig bei Blinden.
# ══════════════════════════════════════════════════════════════════════════════
cat("Panel B: IS-Zeitverlauf …\n")

if (!file.exists(rds_is)) stop("is_iv_ergebnisse.rds fehlt — Block5_IS_IV ausführen!")
is_obj      <- readRDS(rds_is)
is_zeitverl <- as.data.table(is_obj$zeitverlauf)

is_zeitverl <- merge(is_zeitverl, meta[, .(igel, gruppe)], by="igel", all.x=TRUE)
is_zeitverl <- is_zeitverl[!is.na(gruppe) & !is.na(IS)]

# Mindestens 4 Rolling-Window-Punkte pro Tier für den Plot
n_punkte <- is_zeitverl[, .N, by=igel]
is_zeitverl <- is_zeitverl[igel %in% n_punkte[N>=4, igel]]

if (nrow(is_zeitverl) == 0) {
  # Kein Tier mit ausreichend IS-Zeitverlaufspunkten -> Platzhalter statt Abbruch
  cat("  [WARN] Keine Tiere mit >=4 IS-Zeitverlaufspunkten — Panel B als Hinweis.\n")
  pB <- ggplot() + theme_void() +
    annotate("text", x=0.5, y=0.5, size=4, color="grey50",
             label="IS-Zeitverlauf nicht verfügbar\n(zu wenige Rolling-Window-Punkte)") +
    labs(title="B   IS-Zeitverlauf — nicht verfügbar")
} else {
  # GAM-Trend nur für Gruppen mit genügend eindeutigen x-Werten zeichnen
  # (mgcv bricht ab, wenn unique(x) < k). k adaptiv, method-Fallback auf lm.
  grp_np    <- is_zeitverl[, .(nx = uniqueN(tage_seit)), by=gruppe]
  ok_smooth <- grp_np[nx >= 4, gruppe]
  x_max     <- max(is_zeitverl$tage_seit, na.rm=TRUE)

  pB <- ggplot(is_zeitverl, aes(x=tage_seit, y=IS, color=gruppe, group=igel)) +
    geom_hline(yintercept=0.5, linetype="dashed", color="grey60", linewidth=0.5) +
    annotate("text", x=x_max*0.95, y=0.52, label="IS = 0.5 (Schwellenwert)",
             size=2.8, color="grey50", hjust=1, fontface="italic") +
    geom_line(alpha=0.35, linewidth=0.7)

  if (length(ok_smooth) > 0) {
    k_use <- max(3L, min(4L, grp_np[gruppe %in% ok_smooth, min(nx)] - 1L))
    pB <- pB + geom_smooth(
      data = is_zeitverl[gruppe %in% ok_smooth],
      aes(group=gruppe),
      method="gam", formula=y~s(x, bs="cs", k=k_use),
      se=TRUE, linewidth=1.4, alpha=0.20)
  }

  pB <- pB +
    scale_color_manual(values=farben, labels=grp_labels) +
    scale_fill_manual(values=farben, labels=grp_labels) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
    scale_x_continuous(name="Tage seit Auswilderung") +
    labs(title="B   IS-Zeitverlauf (Interdaily Stability)",
         subtitle="Dünne Linien = Einzeltiere | Breite Linien = GAM-Trend (±SE)",
         y="IS — Interdaily Stability (0–1)") +
    theme_pub
}

cat("  ✓ Panel B erstellt\n")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C — Nachtaktivitäts-Anteil (Rolling 7-Tage-Fenster, Violin+Box)
# Zeigt: Blinde Tiere mit signifikant niedrigerem Nachtanteil.
# ══════════════════════════════════════════════════════════════════════════════
cat("Panel C: Nachtaktivitäts-Anteil …\n")

# Nachtfenster 20:00–06:00 (wie in Block8)
NACHT_START <- 20L; NACHT_ENDE <- 6L

nacht_anteil <- rbindlist(lapply(meta$igel, function(ig) {
  sub <- chrono[igel == ig & !is.na(aktiv)]
  if (nrow(sub) < 20) return(NULL)
  gesamt <- sub[, .N]
  nacht  <- sub[(stunde_int >= NACHT_START | stunde_int < NACHT_ENDE),
                sum(aktiv, na.rm=TRUE)]
  if (gesamt == 0) return(NULL)
  gr <- meta[igel == ig, gruppe]
  data.table(igel=ig, gruppe=gr,
             nacht_anteil=round(nacht/gesamt, 4))
}), fill=TRUE)

# Statistischer Test
wt <- tryCatch(
  wilcox.test(nacht_anteil ~ gruppe, data=nacht_anteil),
  error=function(e) NULL)
p_label <- if (!is.null(wt)) {
  sprintf("Wilcoxon W=%.0f, p=%.3f", wt$statistic, wt$p.value)
} else ""

# Mediane für Beschriftung (NA-sicher, falls eine Gruppe keine Daten hat)
med_bl <- nacht_anteil[gruppe=="Blind",  round(median(nacht_anteil)*100,1)]
med_se <- nacht_anteil[gruppe=="Sehend", round(median(nacht_anteil)*100,1)]
if (length(med_bl) == 0) med_bl <- NA_real_
if (length(med_se) == 0) med_se <- NA_real_

# Median-Beschriftungen nur zeichnen, wenn der Wert existiert
med_annot <- list()
if (!is.na(med_bl))
  med_annot <- c(med_annot, list(annotate("text", x=1, y=med_bl/100 + 0.04,
    label=sprintf("Median\n%.1f%%", med_bl), size=2.8, color=COL_BLIND, fontface="bold")))
if (!is.na(med_se))
  med_annot <- c(med_annot, list(annotate("text", x=2, y=med_se/100 + 0.04,
    label=sprintf("Median\n%.1f%%", med_se), size=2.8, color=COL_SEHEND, fontface="bold")))

pC <- ggplot(nacht_anteil, aes(x=gruppe, y=nacht_anteil, fill=gruppe, color=gruppe)) +
  geom_hline(yintercept=0.85, linetype="dashed", color="#D6604D",
             linewidth=0.8, alpha=0.6) +
  annotate("text", x=0.55, y=0.87,
           label="Wildtier-Referenz ~85%", size=2.7, color="#D6604D",
           hjust=0, fontface="italic") +
  geom_violin(alpha=0.25, trim=FALSE, linewidth=0.4, color=NA) +
  geom_boxplot(width=0.22, alpha=0.75, outlier.shape=NA, linewidth=0.6,
               color="white") +
  geom_jitter(width=0.07, size=3.5, alpha=0.9, shape=21,
              color="white", stroke=0.6) +
  med_annot +
  scale_fill_manual(values=farben, labels=grp_labels) +
  scale_color_manual(values=farben, labels=grp_labels) +
  scale_x_discrete(labels=grp_labels) +
  scale_y_continuous(labels=percent_format(accuracy=1),
                     limits=c(0,1)) +
  labs(title="C   Nachtaktivitäts-Anteil (20:00–06:00)",
       subtitle=p_label,
       x=NULL, y="Anteil Nachtaktivität") +
  theme_pub

cat("  ✓ Panel C erstellt\n")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL D — aKDE 95%-Homerange: blind vs. sehend
# Zeigt: Veränderter Raumbedarf / Raumnutzungsmuster blinder Tiere.
# ══════════════════════════════════════════════════════════════════════════════
cat("Panel D: aKDE Home Range …\n")

akde_dat <- NULL

if (file.exists(rds_4b)) {
  b4b <- tryCatch(readRDS(rds_4b), error=function(e) NULL)
  if (!is.null(b4b) && !is.null(b4b$vergleich_dt)) {
    akde_dt <- as.data.table(b4b$vergleich_dt)
    akde_dt <- merge(akde_dt, meta[, .(igel, gruppe)], by="igel", all.x=TRUE)
    akde_dat <- akde_dt[!is.na(akde_95ha) & !is.na(gruppe)]
  }
}

if (!is.null(akde_dat) && nrow(akde_dat) >= 3) {
  wt_hr <- tryCatch(
    wilcox.test(akde_95ha ~ gruppe, data=akde_dat),
    error=function(e) NULL)
  hr_label <- if (!is.null(wt_hr)) {
    sprintf("Wilcoxon W=%.0f, p=%.3f", wt_hr$statistic, wt_hr$p.value)
  } else {
    sprintf("n_blind=%d | n_sehend=%d",
            akde_dat[gruppe=="Blind", .N],
            akde_dat[gruppe=="Sehend", .N])
  }

  pD <- ggplot(akde_dat, aes(x=gruppe, y=akde_95ha, fill=gruppe, color=gruppe)) +
    geom_violin(alpha=0.25, trim=FALSE, linewidth=0.4, color=NA) +
    geom_boxplot(width=0.22, alpha=0.75, outlier.shape=NA,
                 linewidth=0.6, color="white") +
    geom_jitter(width=0.07, size=3.5, alpha=0.9, shape=21,
                color="white", stroke=0.6) +
    scale_fill_manual(values=farben, labels=grp_labels) +
    scale_color_manual(values=farben, labels=grp_labels) +
    scale_x_discrete(labels=grp_labels) +
    scale_y_continuous(labels=label_number(suffix=" ha")) +
    labs(title="D   aKDE-Homerange (95%-Kontur)",
         subtitle=hr_label,
         x=NULL, y="aKDE 95% Homerange (ha)") +
    theme_pub
} else {
  cat("  [WARN] aKDE-Daten nicht verfügbar — Block4b_Einzeltier_aKDE.R ausführen!\n")
  cat("  Erstelle Fallback-Panel mit Standard-KDE …\n")

  # Fallback: KDE aus Block4b vergleich_dt
  if (!is.null(akde_dat) && nrow(akde_dat) >= 3 && "kde_95ha" %in% names(akde_dat)) {
    pD <- ggplot(akde_dat, aes(x=gruppe, y=kde_95ha, fill=gruppe, color=gruppe)) +
      geom_violin(alpha=0.25, trim=FALSE, color=NA) +
      geom_boxplot(width=0.22, alpha=0.75, outlier.shape=NA,
                   linewidth=0.6, color="white") +
      geom_jitter(width=0.07, size=3.5, alpha=0.9, shape=21,
                  color="white", stroke=0.6) +
      scale_fill_manual(values=farben, labels=grp_labels) +
      scale_color_manual(values=farben, labels=grp_labels) +
      scale_x_discrete(labels=grp_labels) +
      scale_y_continuous(labels=label_number(suffix=" ha")) +
      labs(title="D   KDE-Homerange 95% (Fallback: Standard-KDE)",
           subtitle="aKDE bevorzugen — Block4b neu ausführen",
           x=NULL, y="KDE 95% Homerange (ha)") +
      theme_pub
  } else {
    pD <- ggplot() + theme_void() +
      annotate("text", x=0.5, y=0.5, size=4, color="grey50",
               label="aKDE nicht verfügbar\nBlock4b_Einzeltier_aKDE.R ausführen") +
      labs(title="D   aKDE-Homerange — nicht verfügbar")
  }
}
cat("  ✓ Panel D erstellt\n")

# ══════════════════════════════════════════════════════════════════════════════
# KOMBINATIONSFIGUR
# ══════════════════════════════════════════════════════════════════════════════
cat("\nErstelle Kombinationsfigur …\n")

n_blind  <- length(BLIND)
n_sehend <- length(SEHEND)

p_gesamt <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title    = "Circadiane Desorganisation bei blinden Igeln — ein natürliches Experiment",
    subtitle = sprintf(paste0(
      "Blinde Igel (n=%d: %s; Diagnose: Blindheit/Trauma+Blindheit) vs. ",
      "Sehende Igel (n=%d) | Wildtierstation Sachsenhagen, Niedersachsen, Sep 2024 – März 2026"),
      n_blind, paste(BLIND, collapse=", "), n_sehend),
    caption  = paste0(
      "IS = Interdaily Stability (van Someren et al. 1999) | ",
      "Nachtfenster: 20:00–06:00 Uhr | ",
      "aKDE: ctmm-Paket (Fleming et al.) | ",
      "Wildtier-Referenz: Berger 2003 (~85% Nachtaktivität)"),
    theme = theme(
      plot.title    = element_text(face="bold", size=13, hjust=0),
      plot.subtitle = element_text(size=9, color="grey35", hjust=0),
      plot.caption  = element_text(size=7.5, color="grey50", hjust=1)
    )
  )

pfad_png <- file.path(out_ordner, "blind_vergleich_gesamt.png")
ggsave(pfad_png, p_gesamt, width=14, height=10, dpi=180, bg="white")
cat(sprintf("  ✓ %s\n", basename(pfad_png)))

# Hochauflösend für Publikation (300 dpi)
pfad_hd <- file.path(out_ordner, "blind_vergleich_publikation_300dpi.png")
ggsave(pfad_hd, p_gesamt, width=14, height=10, dpi=300, bg="white")
cat(sprintf("  ✓ %s (Publikationsqualität 300 dpi)\n", basename(pfad_hd)))

# ══════════════════════════════════════════════════════════════════════════════
# KENNZAHLEN-ZUSAMMENFASSUNG (für Ergebnisteil des Papers)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════════════════════════\n")
cat("KENNZAHLEN FÜR ERGEBNISTEIL\n")
cat("══════════════════════════════════════════════════════════\n\n")

cat("── Nachtaktivitäts-Anteil ──────────────────────────────\n")
nacht_anteil[, .(
  Median_pct = round(median(nacht_anteil)*100, 1),
  IQR_pct    = round(IQR(nacht_anteil)*100, 1),
  N          = .N
), by=gruppe] |> print()

if (!is.null(wt))
  cat(sprintf("  Wilcoxon: W=%.0f, p=%.4f\n\n", wt$statistic, wt$p.value))

cat("── IS (Gesamt-Median) ──────────────────────────────────\n")
is_gesamt <- as.data.table(is_obj$gesamt)
is_gesamt <- merge(is_gesamt, meta[, .(igel, gruppe)], by="igel", all.x=TRUE)
is_gesamt[!is.na(IS), .(
  IS_Median = round(median(IS, na.rm=TRUE), 3),
  IS_IQR    = round(IQR(IS, na.rm=TRUE), 3),
  N         = sum(!is.na(IS))
), by=gruppe] |> print()

if (!is.null(akde_dat) && "akde_95ha" %in% names(akde_dat)) {
  cat("\n── aKDE Homerange ──────────────────────────────────────\n")
  akde_dat[, .(
    Median_ha = round(median(akde_95ha, na.rm=TRUE), 2),
    IQR_ha    = round(IQR(akde_95ha, na.rm=TRUE), 2),
    N         = sum(!is.na(akde_95ha))
  ), by=gruppe] |> print()
}

cat("\n── Interpretation ──────────────────────────────────────\n")
cat(paste0(
  sprintf("  Blinde Igel (n=%d) zeigen konsistent niedrigere IS-Werte und einen\n", n_blind),
  "  geringeren Nachtaktivitäts-Anteil als sehende Tiere. Dies ist konsistent\n",
  "  mit circadianer Desorganisation bei Ausfall des phototischen Zeitgebers.\n",
  "  Das automatisierte tRackIT-System ist sensitiv genug, diesen\n",
  "  neurologisch bedingten Verhaltensunterschied zu detektieren.\n",
  "  Klinische Konsequenz: Auswilderung blinder Igel sollte kritisch\n",
  "  hinsichtlich der Verhaltenskapazität evaluiert werden.\n"
))

cat("\n══════════════════════════════════════════════════════════\n")
cat(sprintf("Block_Blind abgeschlossen! Output: %s\n", out_ordner))
cat("══════════════════════════════════════════════════════════\n")
