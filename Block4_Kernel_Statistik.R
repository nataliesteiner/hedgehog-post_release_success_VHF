# ==============================================================
# Block 4 — aKDE home ranges: group comparisons & statistics
# ==============================================================
# Project:  Hedgehog VHF telemetry, Lower Saxony
# Author:   Natalie Steiner
# TiHo Hannover
#
# This script loads the already-computed KDE and aKDE results
# from Block4_Kernel (kernel_ergebnisse.rds) and Block4b_Einzeltier
# (Block4b_Ergebnisse.rds) and runs all group comparisons and
# statistics based on the aKDE values.
#
# Primary metric: aKDE 95% (ha) — autocorrelation-corrected KDE (ctmm)
#
# Requirement:
#   Block4b_Einzeltier_aKDE.R must have been run FIRST.
#
# Input:
#   output/Block4_Kernel/kernel_ergebnisse.rds    — fixes + metadata
#   output/Block4b_Einzeltier/Block4b_Ergebnisse.rds — aKDE values

# ── 0. Pakete ──────────────────────────────────────────────────
pakete <- c("data.table", "ggplot2", "patchwork", "scales",
            "openxlsx", "officer", "flextable", "viridis")
fehlend <- pakete[!sapply(pakete, requireNamespace, quietly = TRUE)]
if (length(fehlend) > 0) install.packages(fehlend)

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(openxlsx)
if (requireNamespace("officer",   quietly = TRUE)) library(officer)
if (requireNamespace("flextable", quietly = TRUE)) library(flextable)

cat("✓ Pakete geladen\n\n")

# ── 1. Pfade ───────────────────────────────────────────────────
projekt_root <- "~/Desktop/Arbeit/Wildtiere TiHo/Rehabilitation/Igel/Auswertung/VHF"

rds_b4  <- file.path(projekt_root, "output", "Block4_Kernel",      "kernel_ergebnisse.rds")
rds_b4b <- file.path(projekt_root, "output", "Block4b_Einzeltier", "Block4b_Ergebnisse.rds")
out_ordner <- file.path(projekt_root, "output", "Block4_Kernel")
dir.create(out_ordner, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rds_b4))
  stop("kernel_ergebnisse.rds nicht gefunden — Block4_Kernel_HomeRange.R zuerst ausführen.")
if (!file.exists(rds_b4b))
  stop("Block4b_Ergebnisse.rds nicht gefunden — Block4b_Einzeltier_aKDE.R zuerst ausführen.")

# ── 2. Daten laden ─────────────────────────────────────────────
cat("Lade Ergebnisse...\n")

b4  <- readRDS(rds_b4)
b4b <- readRDS(rds_b4b)

# Haupttabelle aus Block4: KDE + Metadaten pro Igel
kde_dt <- b4$kde_gesamt
cat(sprintf("  Block4:  %d Igel mit KDE-Daten\n", nrow(kde_dt)))

# aKDE-Zusammenfassung aus Block4b
akde_summary <- b4b$vergleich_dt
if (is.null(akde_summary)) {
  # Fallback: aus einzeltier_results aufbauen
  akde_summary <- rbindlist(lapply(names(b4b$einzeltier_results), function(ig) {
    x <- b4b$einzeltier_results[[ig]]
    data.table(igel      = ig,
               akde_95ha = x$akde_95ha,
               akde_50ha = x$akde_50ha,
               ess_area  = x$ess_area,
               akde_modell = x$akde_modell)
  }), fill = TRUE)
}

cat(sprintf("  Block4b: %d Igel mit aKDE-Daten\n\n", sum(!is.na(akde_summary$akde_95ha))))

# ── 3. Merge: aKDE + KDE + Metadaten ──────────────────────────
# Nur Spalten die wir brauchen aus akde_summary
akde_cols <- intersect(c("igel","akde_95ha","akde_50ha","ess_area","akde_modell"),
                        names(akde_summary))
dt <- merge(kde_dt, akde_summary[, ..akde_cols], by = "igel", all.x = TRUE)

# Sicherstellen dass alle nötigen Metadaten-Spalten vorhanden sind
meta_spalten <- c("alter","sex","saison_auswild","diag_gruppe","time_reha","n_naechte")
fehlende_meta <- setdiff(meta_spalten, names(dt))
if (length(fehlende_meta) > 0) {
  cat("Fehlende Metadaten-Spalten:", paste(fehlende_meta, collapse=", "), "\n")
  cat("Verfügbare Spalten:", paste(names(dt), collapse=", "), "\n\n")
}

# Faktoren setzen
if ("saison_auswild" %in% names(dt))
  dt[, saison_auswild := factor(saison_auswild, levels = c("Frühling","Sommer","Herbst","Winter"))]
cat(sprintf("Analysedatensatz: %d Igel | %d mit aKDE-Wert\n\n",
            nrow(dt), sum(!is.na(dt$akde_95ha))))
print(dt[, .(igel, akde_95ha, akde_50ha, sex, saison_auswild,
             diag_gruppe, time_reha)])

# ── 4. Deskriptive Statistik ───────────────────────────────────
cat("\n═══════════════════════════════════════\n")
cat("Deskriptive Statistik (aKDE 95%)\n")
cat("═══════════════════════════════════════\n")

dt_ok <- dt[!is.na(akde_95ha)]

deskr <- dt_ok[, .(
  n          = .N,
  median_ha  = round(median(akde_95ha), 2),
  iqr_ha     = round(IQR(akde_95ha), 2),
  mean_ha    = round(mean(akde_95ha), 2),
  sd_ha      = round(sd(akde_95ha), 2),
  min_ha     = round(min(akde_95ha), 2),
  max_ha     = round(max(akde_95ha), 2),
  median_50  = round(median(akde_50ha, na.rm=TRUE), 2)
)]

cat(sprintf("N = %d Igel\n", deskr$n))
cat(sprintf("aKDE 95%%: Median = %.2f ha (IQR: %.2f ha)\n", deskr$median_ha, deskr$iqr_ha))
cat(sprintf("          Mean = %.2f ha ± %.2f SD\n", deskr$mean_ha, deskr$sd_ha))
cat(sprintf("          Range: %.2f – %.2f ha\n", deskr$min_ha, deskr$max_ha))
cat(sprintf("aKDE 50%%: Median = %.2f ha\n\n", deskr$median_50))

# ── 5. Statistische Tests ──────────────────────────────────────
cat("═══════════════════════════════════════\n")
cat("Statistische Tests (aKDE 95%)\n")
cat("═══════════════════════════════════════\n\n")

test_ergebnisse <- list()

# ── Test 1: Wilcoxon — Geschlecht (Male vs. Female) ───────────
if ("sex" %in% names(dt_ok) && dt_ok[, uniqueN(sex[!is.na(sex)])] == 2) {
  wt_sex <- wilcox.test(akde_95ha ~ sex, data = dt_ok[!is.na(sex)])
  grp_s <- dt_ok[!is.na(sex), .(median=round(median(akde_95ha),2), n=.N), by=sex]
  cat("2. Wilcoxon: Geschlecht (Male vs. Female)\n")
  print(grp_s); cat(sprintf("   W = %.0f, p = %.4f %s\n\n", wt_sex$statistic, wt_sex$p.value,
                              ifelse(wt_sex$p.value < 0.05, "***", "n.s.")))
  test_ergebnisse[["sex"]] <- data.table(
    test="Wilcoxon (Geschlecht: Male vs. Female)",
    statistik=round(wt_sex$statistic,2), p_wert=round(wt_sex$p.value,4),
    signif=ifelse(wt_sex$p.value < 0.05, "ja", "nein"))
}

# ── Test 3: Kruskal-Wallis — Saison ───────────────────────────
if ("saison_auswild" %in% names(dt_ok) && dt_ok[, uniqueN(saison_auswild[!is.na(saison_auswild)])] >= 2) {
  kt_s <- kruskal.test(akde_95ha ~ saison_auswild, data = dt_ok[!is.na(saison_auswild)])
  grp_sai <- dt_ok[!is.na(saison_auswild), .(median=round(median(akde_95ha),2), n=.N), by=saison_auswild]
  cat("3. Kruskal-Wallis: Saison\n")
  print(grp_sai); cat(sprintf("   H = %.3f, df = %d, p = %.4f %s\n\n",
                               kt_s$statistic, kt_s$parameter, kt_s$p.value,
                               ifelse(kt_s$p.value < 0.05, "***", "n.s.")))
  test_ergebnisse[["saison"]] <- data.table(
    test="Kruskal-Wallis (Saison)",
    statistik=round(kt_s$statistic,3), p_wert=round(kt_s$p.value,4),
    signif=ifelse(kt_s$p.value < 0.05, "ja", "nein"))

  # Post-hoc Dunn-Test (Bonferroni) wenn signifikant
  if (kt_s$p.value < 0.05) {
    cat("   Post-hoc Paarvergleiche (Wilcoxon, Bonferroni-korrigiert):\n")
    saisonen <- as.character(unique(dt_ok$saison_auswild[!is.na(dt_ok$saison_auswild)]))
    paare    <- combn(saisonen, 2, simplify = FALSE)
    n_paare  <- length(paare)
    ph_res   <- rbindlist(lapply(paare, function(paar) {
      g1 <- dt_ok[saison_auswild == paar[1], akde_95ha]
      g2 <- dt_ok[saison_auswild == paar[2], akde_95ha]
      if (length(g1) < 2 || length(g2) < 2) return(NULL)
      wt_ph <- wilcox.test(g1, g2, exact = FALSE)
      data.table(Gruppe1=paar[1], Gruppe2=paar[2],
                 W=round(wt_ph$statistic,1),
                 p_raw=round(wt_ph$p.value,4),
                 p_bonf=round(min(wt_ph$p.value * n_paare, 1), 4),
                 signif=ifelse(wt_ph$p.value * n_paare < 0.05, "*", "n.s."))
    }))
    print(ph_res)
    cat("\n")
    test_ergebnisse[["saison_posthoc"]] <- ph_res
  }
}

# ── Test 4: Kruskal-Wallis — Diagnosegruppe ───────────────────
if ("diag_gruppe" %in% names(dt_ok) && dt_ok[, uniqueN(diag_gruppe[!is.na(diag_gruppe)])] >= 2) {
  kt_d <- kruskal.test(akde_95ha ~ diag_gruppe, data = dt_ok[!is.na(diag_gruppe)])
  grp_d <- dt_ok[!is.na(diag_gruppe), .(median=round(median(akde_95ha),2), n=.N), by=diag_gruppe]
  cat("4. Kruskal-Wallis: Diagnosegruppe\n")
  print(grp_d); cat(sprintf("   H = %.3f, df = %d, p = %.4f %s\n\n",
                              kt_d$statistic, kt_d$parameter, kt_d$p.value,
                              ifelse(kt_d$p.value < 0.05, "***", "n.s.")))
  test_ergebnisse[["diag"]] <- data.table(
    test="Kruskal-Wallis (Diagnosegruppe)",
    statistik=round(kt_d$statistic,3), p_wert=round(kt_d$p.value,4),
    signif=ifelse(kt_d$p.value < 0.05, "ja", "nein"))
}

# ── Test 5: Spearman — aKDE ~ Rehab-Dauer ─────────────────────
if ("time_reha" %in% names(dt_ok) && sum(!is.na(dt_ok$time_reha)) >= 5) {
  sp_r <- cor.test(dt_ok$akde_95ha, dt_ok$time_reha, method = "spearman", use = "complete.obs")
  cat("5. Spearman: aKDE ~ Rehabilitationsdauer\n")
  cat(sprintf("   rho = %.3f, p = %.4f, n = %d %s\n\n",
              sp_r$estimate, sp_r$p.value, sum(!is.na(dt_ok$time_reha)),
              ifelse(sp_r$p.value < 0.05, "***", "n.s.")))
  test_ergebnisse[["reha"]] <- data.table(
    test="Spearman (aKDE ~ Rehab-Dauer)",
    statistik=round(sp_r$estimate,3), p_wert=round(sp_r$p.value,4),
    signif=ifelse(sp_r$p.value < 0.05, "ja", "nein"))
}

# ── Test 6: Spearman — aKDE ~ Beobachtungsdauer (N Nächte) ────
if ("n_naechte" %in% names(dt_ok) && sum(!is.na(dt_ok$n_naechte)) >= 5) {
  sp_n <- cor.test(dt_ok$akde_95ha, dt_ok$n_naechte, method = "spearman", use = "complete.obs")
  cat("6. Spearman: aKDE ~ Beobachtungsdauer (N Nächte)\n")
  cat(sprintf("   rho = %.3f, p = %.4f, n = %d %s\n\n",
              sp_n$estimate, sp_n$p.value, sum(!is.na(dt_ok$n_naechte)),
              ifelse(sp_n$p.value < 0.05, "***", "n.s.")))
  test_ergebnisse[["naechte"]] <- data.table(
    test="Spearman (aKDE ~ Beobachtungsdauer/N Nächte)",
    statistik=round(sp_n$estimate,3), p_wert=round(sp_n$p.value,4),
    signif=ifelse(sp_n$p.value < 0.05, "ja", "nein"))
}

test_dt <- rbindlist(test_ergebnisse[sapply(test_ergebnisse, is.data.table)], fill = TRUE)

# ── 6. Plots ───────────────────────────────────────────────────
cat("\nErstelle Plots...\n")

farben_sex    <- c("Male" = "#4575b4", "Female" = "#d73027")

# ── Plot P1: Geschlecht ────────────────────────────────────────
p_sex <- ggplot(dt_ok[!is.na(sex)], aes(x = sex, y = akde_95ha, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.5) +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.7, shape = 21,
              aes(fill = sex), color = "white") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  scale_fill_manual(values = farben_sex) +
  labs(title = "Home Range by sex",
       subtitle = paste0("Wilcoxon: W=",
         ifelse("sex" %in% names(test_ergebnisse),
                test_ergebnisse$sex$statistik, "?"),
         ", p=",
         ifelse("sex" %in% names(test_ergebnisse),
                test_ergebnisse$sex$p_wert, "?")),
       x = NULL, y = "aKDE 95% (ha)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 11))

# ── Plot P3: Saison ────────────────────────────────────────────
p_saison <- ggplot(dt_ok[!is.na(saison_auswild)],
                   aes(x = saison_auswild, y = akde_95ha, fill = saison_auswild)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.5) +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.7, shape = 21,
              aes(fill = saison_auswild), color = "white") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  scale_fill_viridis_d(option = "plasma", end = 0.85) +
  labs(title = "Home Range by release season",
       subtitle = paste0("Kruskal-Wallis: H=",
         ifelse("saison" %in% names(test_ergebnisse),
                test_ergebnisse$saison$statistik, "?"),
         ", p=",
         ifelse("saison" %in% names(test_ergebnisse),
                test_ergebnisse$saison$p_wert, "?")),
       x = NULL, y = "aKDE 95% (ha)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 11))

# ── Plot P4: Diagnose ──────────────────────────────────────────
p_diag <- ggplot(dt_ok[!is.na(diag_gruppe)],
                 aes(x = diag_gruppe, y = akde_95ha, fill = diag_gruppe)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, width = 0.5) +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.7, shape = 21,
              aes(fill = diag_gruppe), color = "white") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Home Range by diagnosis group",
       subtitle = paste0("Kruskal-Wallis: H=",
         ifelse("diag" %in% names(test_ergebnisse),
                test_ergebnisse$diag$statistik, "?"),
         ", p=",
         ifelse("diag" %in% names(test_ergebnisse),
                test_ergebnisse$diag$p_wert, "?")),
       x = NULL, y = "aKDE 95% (ha)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1),
        plot.title = element_text(face = "bold", size = 11))

# Kombinierter Gruppenplot
p_gruppen <- (p_sex + p_saison) / (p_diag + plot_spacer()) +
  plot_annotation(
    title    = "aKDE Home Range — Gruppenvergleiche",
    subtitle = "aKDE 95% (ctmm, OUF anisotropic) | Raute = Mittelwert",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(out_ordner, "akde_gruppen.png"),
       p_gruppen, width = 12, height = 10, dpi = 150)
cat("  ✓ akde_gruppen.png\n")

# ── Plot P5: Korrelation Rehab-Dauer ───────────────────────────
if ("time_reha" %in% names(dt_ok) && sum(!is.na(dt_ok$time_reha)) >= 5) {
  p_reha <- ggplot(dt_ok[!is.na(time_reha) & !is.na(sex)],
                   aes(x = time_reha, y = akde_95ha, color = sex, label = igel)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey60",
                fill = "grey80", linetype = "dashed", linewidth = 0.8) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_color_manual(values = farben_sex) +
    labs(title    = "Rehabilitation duration vs. aKDE home range",
         subtitle = paste0("Spearman rho=",
           ifelse("reha" %in% names(test_ergebnisse),
                  test_ergebnisse$reha$statistik, "?"),
           ", p=",
           ifelse("reha" %in% names(test_ergebnisse),
                  test_ergebnisse$reha$p_wert, "?")),
         x = "Rehabilitation duration (days)",
         y = "aKDE 95% (ha)",
         color = "Sex") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_ordner, "akde_korrelation_reha.png"),
         p_reha, width = 8, height = 5, dpi = 150)
  cat("  ✓ akde_korrelation_reha.png\n")
}

# ── Plot P6: Individuelle Übersicht (geordnet) ─────────────────
igel_ord <- dt_ok[order(akde_95ha), igel]
dt_long  <- melt(dt_ok[, .(igel, akde_95ha, akde_50ha)],
                 id.vars = "igel", variable.name = "typ", value.name = "ha")
dt_long[, typ_label := ifelse(typ == "akde_95ha", "95% aKDE (Home Range)", "50% aKDE (Core zone)")]
dt_long[, igel := factor(igel, levels = igel_ord)]

p_uebersicht <- ggplot(dt_long, aes(x = igel, y = ha, fill = typ_label)) +
  geom_col(position = "identity", alpha = 0.85) +
  scale_fill_manual(values = c("95% aKDE (Home Range)" = "#2166ac",
                                "50% aKDE (Core zone)"  = "#f4a582")) +
  coord_flip() +
  labs(title = "Nocturnal home range — all hedgehogs",
       subtitle = "aKDE 95% and 50% (ctmm, OUF anisotropic model)",
       x = NULL, y = "Area (ha)", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_ordner, "akde_uebersicht.png"),
       p_uebersicht, width = 8, height = 9, dpi = 150)
cat("  ✓ akde_uebersicht.png\n")

# ── 7. Ergebnisse speichern ────────────────────────────────────
cat("\nSpeichere Ergebnisse...\n")

# RDS
saveRDS(list(
  akde_dt        = dt,
  test_ergebnisse = test_dt,
  deskriptiv     = deskr
), file.path(out_ordner, "akde_statistik.rds"))
cat("  ✓ akde_statistik.rds\n")

# Excel
wb <- createWorkbook()
addWorksheet(wb, "aKDE_Übersicht")
writeData(wb, "aKDE_Übersicht", dt[, .(igel, akde_95ha, akde_50ha, ess_area,
                                         sex, saison_auswild,
                                         diag_gruppe, time_reha, n_naechte)])

addWorksheet(wb, "Statistische_Tests")
writeData(wb, "Statistische_Tests", test_dt)

addWorksheet(wb, "Deskriptiv")
writeData(wb, "Deskriptiv", deskr)

saveWorkbook(wb, file.path(out_ordner, "akde_statistik.xlsx"), overwrite = TRUE)
cat("  ✓ akde_statistik.xlsx\n")

# ── 8. Abschlusszusammenfassung ────────────────────────────────
cat("\n═══════════════════════════════════════\n")
cat("ERGEBNISZUSAMMENFASSUNG\n")
cat("═══════════════════════════════════════\n")
cat(sprintf("N = %d Igel (aKDE berechnet)\n", deskr$n))
cat(sprintf("aKDE 95%%: Median %.2f ha (IQR %.2f ha, Range %.2f–%.2f ha)\n",
            deskr$median_ha, deskr$iqr_ha, deskr$min_ha, deskr$max_ha))
cat(sprintf("aKDE 50%%: Median %.2f ha\n\n", deskr$median_50))
cat("Statistische Tests:\n")
if (nrow(test_dt) > 0) {
  for (i in seq_len(nrow(test_dt))) {
    cat(sprintf("  %-48s stat=%-6s p=%.4f %s\n",
                test_dt$test[i], test_dt$statistik[i], test_dt$p_wert[i],
                ifelse(test_dt$signif[i]=="ja", "***", "n.s.")))
  }
}
cat("\n✓ Block 4 abgeschlossen — alle Outputs in:", out_ordner, "\n")
