# Post-release behaviour and rehabilitation success in European hedgehogs

Analysis code for a study on the post-release behaviour and rehabilitation success of
26 European hedgehogs (*Erinaceus europaeus*) released at a single study site in
Sachsenhagen, Lower Saxony, Germany, and tracked with an automated VHF telemetry system
(reverse-GPS / multilateration).

**Author:** Natalie Steiner, University of Veterinary Medicine Hannover (TiHo) / Wildtierstation Sachsenhagen
**Associated paper:** *(add citation / DOI once available)*

---

## What this repository contains

R scripts for the full analysis pipeline. The workflow is organised in numbered "Blocks";
run them in order. Helper and figure scripts can be run once the corresponding Block has produced its outputs.

| Script | Purpose |
|---|---|
| `Block0_Datenpipeline.R` | Activity data pipeline (all animals, batch) |
| `Block0_wmv_Vergleich.R` | Comparison of tRackIT classifications (wmv vs. smoothed) |
| `Block1_Metadaten.R` | Metadata / descriptive analysis of study animals |
| `Block2_GAMM.R` | Hierarchical additive mixed model of night activity over time |
| `Block2b_BoutVergleich.R` | Activity-bout comparison |
| `Block2c_HMM.R` | Hidden Markov model activity classification |
| `Block2d_ReleaseEffect.R` | Post-release effect analysis |
| `Block3_Chronobiologie.R` | Chronobiological analysis of activity after release |
| `Block4_Kernel_HomeRange.R` | Space use: fixes and spatial metrics |
| `Block4_Kernel_Statistik.R` | aKDE home ranges: group comparisons and statistics |
| `Block4b_Einzeltier_aKDE.R` | Individual-level aKDE comparison |
| `Block4d_TemporaleAnalyse.R` | Temporal movement analysis |
| `Block4e_Tagesschlafplaetze.R` | Daytime nest/resting-site analysis |
| `Block5_IS_IV.R` | Interdaily Stability & Intradaily Variability (circadian metrics) |
| `Block5_Umweltdaten.R` | Weather × home range & habitat analysis |
| `Block6_Dispersal.R` | Dispersal distance and speed from the release site |
| `Block7_DetektionsValidierung.R` | Detection validation |
| `Block8_RehabSuccess.R` | Rehabilitation success analysis |
| `Block_Blind_Vergleich.R` | Comparison for blind vs. sighted animals |


## Expected folder layout

The scripts use a project root and expect this structure:

```
<projekt_root>/
├── scripts/            # this repository
├── data/
│   ├── activity/           # activity classification CSVs (tRackIT output)
│   ├── excel_files/        # data_igel.xlsx (animal metadata) and habitat layers
│   └── kernel_files/       # multilateration .gpkg files (localizations)
└── output/                 # created by the scripts
```

At the top of the Block scripts, set `projekt_root` to the folder that contains `data/`
and `output/` on your machine.

## Data availability

The raw tracking data are **not** included in this repository (size and access reasons). The author can be contacted for request.

## Requirements

R (>= 4.x). Core packages used across the pipeline:

`data.table`, `dplyr`, `tidyr`, `lubridate`, `readxl`, `openxlsx`,
`sf`, `sp`, `adehabitatHR`, `ctmm`, `raster`, `tidyterra`, `maptiles`, `osmdata`,
`rnaturalearth`, `rnaturalearthdata`, `ggspatial`,
`mgcv`, `gratia`, `lme4`, `lmerTest`, `circular`, `suncalc`, `rdwd`,
`ggplot2`, `patchwork`, `scales`, `ggrepel`, `ggdist`, `viridis`,
`officer`, `flextable`.

Run `sessionInfo_report.R` after a full analysis run to record the exact package
versions used (produces `R_package_versions.csv`, `sessionInfo.txt`, `R_citations.bib`).

## How to cite

If you use this code, please cite the associated paper and this repository

## License

MIT License

Copyright (c) 2026 Natalie Steiner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Contact

Natalie Steiner
natalie.steiner@tiho-hannover.de (institutional) · natalie-steiner@hotmail.com (permanent)
