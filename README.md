# bakanposs

Decomposing the satellite–tower ET divergence between MSG-based
satellite products (LSA SAF Meteosat ETv3) and two eddy-covariance
flux towers in Mediterranean Spain:

- **Oran** — rainfed winter cereal (2018–2020)
- **Tarazona** — drip-irrigated almond orchard (2020–2024)

---

## Quick start

```bash
git clone <repo-url> && cd bakanposs
pip install -r requirements.txt

# Reproduce the main poster figures
python3 analysis_A_v31.py         # Fig 4   — recovery τ comparison
python3 analysis_A_v32.py         # Fig 4b  — Tarazona blind-spot decomposition
python3 analysis_B_v6_driver_attribution.py   # Fig — driver attribution
```

Outputs land under `output/analysis_A/v31/`, `output/analysis_A/v32/`,
`output/analysis_B/v6/`.

---

## Headline result

| | Oran (rainfed) | Tarazona (drip) |
|---|---|---|
| Pearson r (EC vs ETv3) | **0.82** | **0.07** |
| Bias | +1 % | **−70 %** |
| KGE | +0.61 | −0.21 |
| Annual ET peak shift | +6 d | **−44 d** (Sat earlier) |

At rainfed Oran, Meteosat ETv3 reproduces EC tower ET, NDVI, and GPP
faithfully.  At drip-irrigated Tarazona, NDVI and GPP still agree, but
ET breaks down — the satellite peaks ~6 weeks earlier than the tower
and captures none of the irrigation amplitude.

The blind spot is structural, not a calibration issue: the 1-ha orchard
is invisible to the 5-km satellite chain because of three combined
resolution limits (spatial dilution, surface-only microwave SM, and an
SVAT model calibrated for natural vegetation).

---

## Repository layout

```
bakanposs/
├── README.md                    this file
├── CLAUDE.md                    AI-assistant navigation hub
├── requirements.txt
├── data_loaders.py              common data loaders
├── run_analysis_C.py            Analysis-C runner
│
├── analysis_A_v31.py            CURRENT — recovery τ comparison (poster Fig 4)
├── analysis_A_v32.py            CURRENT — Tarazona blind-spot (poster Fig 4b)
├── analysis_B_v3_bias_tau.py    CURRENT — bias-pool builder (input to v32)
├── analysis_B_v6_driver_attribution.py  CURRENT — driver bars
├── analysis_C_v2_ndvi_phenology.py      CURRENT — NDVI phenology
│
├── docs/                        all written documentation
│   ├── RESEARCH_OVERVIEW.md
│   ├── ANALYSIS_A_FINAL.md      Analysis-A consolidated results
│   ├── ANALYSIS_A_FAQ.md        common questions
│   ├── ANALYSIS_B_PLAN.md       Analysis-B design
│   ├── ANALYSIS_C_PLAN.md       Analysis-C design
│   ├── SATELLITE_ET_NOTES.md    satellite ET caveats and references
│   ├── REPORT_to_site_collaborators.md  progress report for site PIs
│   ├── paper_outline.md / paper_methods_results.md
│   ├── RUN_ANALYSIS_C.md        Analysis-C runner notes
│   └── sessions/                cumulative session history
│
├── data/                        input data (master_full_v2.csv etc.)
├── output/                      all script outputs
│   ├── analysis_A/v15–v32/
│   ├── analysis_B/v1–v6/
│   ├── analysis_C/v1, v2, last/
│   └── bias_stats/
│
├── archive/scripts/             legacy script versions (history preserved)
├── scripts/                     helper scripts (GEE, bias stats, loaders)
├── poster/                      A0-portrait poster template
└── reports/                     auxiliary reports
```

---

## Document index

| Topic | File |
|---|---|
| Research overview | [docs/RESEARCH_OVERVIEW.md](docs/RESEARCH_OVERVIEW.md) |
| Analysis A — recovery τ | [docs/ANALYSIS_A_FINAL.md](docs/ANALYSIS_A_FINAL.md) |
| Analysis B — satellite ET | [docs/ANALYSIS_B_PLAN.md](docs/ANALYSIS_B_PLAN.md) |
| Analysis C — NDVI phenology | [docs/ANALYSIS_C_PLAN.md](docs/ANALYSIS_C_PLAN.md) |
| Satellite ET caveats | [docs/SATELLITE_ET_NOTES.md](docs/SATELLITE_ET_NOTES.md) |
| Progress report to site PIs | [docs/REPORT_to_site_collaborators.md](docs/REPORT_to_site_collaborators.md) |
| Paper draft | [docs/paper_methods_results.md](docs/paper_methods_results.md) |
| Session history | [docs/sessions/](docs/sessions/) |

---

## Status

- ✅ Analysis A — τ comparison, MDE, amplitude scaling (v31)
- ✅ Analysis B — bias recovery + driver attribution (v3, v6)
- 🟡 Analysis C — NDVI phenology (v2 in progress)
- 🟡 Poster — A0-portrait draft ready, narrative refinement ongoing
- ⬜ H-SAF SM products (H141 / H142 / H26 / H28) — to be fetched and
  substituted into the Analysis B driver-attribution regression
- ⬜ Manuscript — outline drafted, target *Agricultural and Forest
  Meteorology* or *Remote Sensing of Environment*

---

## Citation / contact

Shion Nagamine, Kazuhito Ichii
Center for Environmental Remote Sensing (CEReS)
Chiba University, Japan

Data kindly provided by the Tarazona and Oran flux tower teams (see
[docs/REPORT_to_site_collaborators.md](docs/REPORT_to_site_collaborators.md)
for acknowledgments).
