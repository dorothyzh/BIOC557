# BIOC557 · Multi-omics Lab

**Student website: https://dorothyzh.github.io/BIOC557/**

Interactive teaching lab for RNA, proteomics and metabolomics, prepared for BIOC557 by Huan Zhong, UBC.

## Additional 2025 demos

- **Human pancreatic islets:** six RRID-matched donors (3 non-diabetic, 3 T2D), with RNA, proteomics and untargeted metabolite features. ACS Omega (2025), https://doi.org/10.1021/acsomega.4c10637 . Technical metabolomics injections are collapsed per donor. This small cohort illustrates methods and confounding, not clinical prediction.
- **Honeybee eggs:** 81 exact queen-ID matches (18 May, 63 June), with lipidomics and metabolomics. Chapman et al., Scientific Reports (2025), https://doi.org/10.1038/s41598-025-94670-2 . Two views are used; no transcriptome or proteome is fabricated.

Each new demo has source-file hashes, processing steps, a dataset download and runnable R code. Human GO/KEGG catalogs are not applied to bee features or unidentified human m/z–retention-time signals.

## Presentation and figure exports

The animated landing visual illustrates molecular layers and shared factors; it is conceptual, can be paused, and respects reduced-motion preferences. Text is larger throughout. Quantitative plots use white backgrounds, no grid, labeled axes and SVG/high-resolution PNG/source-data downloads. Network coordinates are illustrative layouts, not measured quantities. Factor interpretation adds descriptive metadata associations and cross-view top-feature correlations, with observed-pair counts and no causal claim.

## Use the lab

Open the student website in a desktop browser. Start with quality control, inspect each molecular layer, fit joint PCA and MOFA, compare factor counts and preparation choices, explore networks and GO/pathways, then export your workspace. Each browser performs its own calculations in an independent Web Worker. The first visit downloads the scientific runtime; results remain in that browser tab until exported or reloaded. No account or local Python installation is required.

The R teaching kit and executed R walkthrough are linked from the site. Students can perform the same methods using R/MOFA2.

## Public teaching data

61 uniquely matched CCLE/DepMap cell-line identities: 28 breast and 33 skin lines. Inputs contain 1,200 protein features, 1,200 RNA genes, and 225 metabolite features; the default model selects 100 variable features per view. These are different public assays linked by cell-line identity, not a common aliquot or culture batch.

- Nusinow et al. *Quantitative Proteomics of the Cancer Cell Line Encyclopedia*. Cell, 2020. https://doi.org/10.1016/j.cell.2019.12.023
- Ghandi et al. *Next-generation characterization of the Cancer Cell Line Encyclopedia*. Nature, 2019. https://doi.org/10.1038/s41586-019-1186-3 (RNA release: DepMap public_20q3).
- Li et al. *The landscape of cancer cell line metabolism*. Nature Medicine, 2019. https://doi.org/10.1038/s41591-019-0404-8

Source URLs, processing, and original file hashes are recorded in `data/provenance.json`. Annotation versions and mappings are included in `data/enrichment_annotations.json`. See the site's Data & methods page and `THIRD_PARTY_NOTICES.txt` for provenance and upstream terms.

## Hosting and updates

GitHub Pages publishes the `main` branch from `/`. `.nojekyll` preserves the static assets. All application/data/download URLs are relative to this project path. The browser calculation engine uses the official Pyodide distribution and bundled mofapy2; no classroom API server is needed for this hosted edition.

To preview locally: `python3 -m http.server 8000`, then open `http://localhost:8000/`.
