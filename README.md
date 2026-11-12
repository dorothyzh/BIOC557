# BIOC557 multi-omics lab

[Interactive lab](https://dorothyzh.github.io/BIOC557/) · [R setup guide](https://dorothyzh.github.io/BIOC557/R_START_HERE.html)

## One-click default report in RStudio

1. Install R (4.3 or newer) and RStudio Desktop.
2. Download a teaching ZIP from the links below and extract it completely.
3. Double-click `BIOC557.Rproj`, open `START_HERE.R`, and click **Source**.

Missing packages install into the project folder automatically. The first run needs internet access. The complete HTML report opens after the default analysis finishes. No parameter or path edits are needed. For later runs, open `BIOC557_Default_Lab.Rmd` and click **Knit**.

| Default study | Data + R Markdown | Completed report |
| --- | --- | --- |
| Honeybee head: protein + RNA, MOFA K=5 | [Download](https://dorothyzh.github.io/BIOC557/downloads/BIOC557_bee-own-2025_R.zip) | [Preview](https://dorothyzh.github.io/BIOC557/downloads/BIOC557_Bee_Default_Example.html) |
| CCLE: protein + RNA + metabolites, MOFA K=5 | [Download](https://dorothyzh.github.io/BIOC557/downloads/BIOC557_Three_Omics_R_Kit.zip) | [Preview](https://dorothyzh.github.io/BIOC557/downloads/BIOC557_Three_Omics_R_Example.html) |
| Human pancreatic islets: three views, PCA | [Download](https://dorothyzh.github.io/BIOC557/downloads/BIOC557_human-islets-2025_R.zip) | [Preview](https://dorothyzh.github.io/BIOC557/downloads/BIOC557_Islets_Default_Example.html) |

The repository root is also a runnable RStudio project with the honeybee-head default. The data, functions, R Markdown and startup script are included. Each report exports all its figures as PNG/PDF/SVG, model matrices, aligned gene IDs, result tables and software versions under `outputs/`.

The interactive website supports additional parameter experiments. The beginner reports run a single fixed default workflow. Data provenance and analysis methods are documented in [Data & methods](https://dorothyzh.github.io/BIOC557/sources.html).
