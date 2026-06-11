# Running Analysis C

This document describes how to run Analysis C (v1 legacy and v2 phenology) with
the package-based layout introduced in 2026-06.

## Install (one time)

```bash
# clone + create venv
python3 -m venv .venv
source .venv/bin/activate

# editable install — picks up `bakanposs/` and all submodules
pip install -e .
```

## Configure data paths

Site / data paths are defined in `bakanposs/analysis_a.py::PATHS`. If your data
lives elsewhere, edit that dict OR (recommended) override at runtime via
environment variables — see `.env.example` for the variable names supported by
`run_analysis_C.py`:

```bash
cp .env.example .env
$EDITOR .env
set -a; source .env; set +a
```

Required inputs:
- `BAKANPOSS_ORAN_EC`   : Oran AmeriFlux half-hourly CSV
- `BAKANPOSS_TARA_EC`   : Tarazona daily summary CSV
- `BAKANPOSS_NDVI_FILE` : MOD13Q1 16-day NDVI/EVI AppEEARS CSV (covers both sites)

## Run

Three equivalent ways:

```bash
# 1) thin wrapper (recommended)
python analyses/run_analysis_C.py --version v2

# 2) python -m module form
python -m bakanposs.analysis_c.v2_phenology

# 3) console script (installed by pyproject.toml)
bakanposs-C-v2
```

Substitute `--version v1` (or the `v1_legacy` module) to run the older
multi-purpose Analysis C v1 instead.

## Output

Default output directory: `./output_analysis_C_v2/` (v2) or
`./output_analysis_C_v1/` (v1).

v2 writes:
- `v2_phenology_metrics.csv`       — per-year SOS/Peak/EOS/season
- `v2_active_period_validation.csv`— NDVI∩assumed IoU per year
- `v2_ndvi_active_tau.csv`         — τ ± CI per window (SOS-EOS, peak ±30d)
- `v2_ndvi_le_lag.csv`             — r(NDVI, LE) at lags ±20d
- `v2_{site}_merged.csv`           — EC + NDVI merge for downstream reuse
- `fig01_ndvi_seasonal_curve.png`
- `fig02_active_period_overlap.png`
- `fig03_ndvi_window_tau.png`
- `fig04_ndvi_le_lag.png`

See `docs/ANALYSIS_C_PLAN.md` for the scientific design, `reports/analysis_C_report.md`
for previous v1 results, and `docs/RESEARCH_OVERVIEW.md` for how Analysis C
relates to A and B.

## Troubleshooting

| symptom | likely cause | fix |
|---|---|---|
| `ModuleNotFoundError: bakanposs` | venv not active or package not installed | `source .venv/bin/activate && pip install -e .` |
| `FileNotFoundError` on EC / NDVI | paths in `PATHS` don't match your machine | edit `bakanposs/analysis_a.py` or set env vars |
| `[events Oran] n=0` | raw AmeriFlux precip column not picked up | inspect the `[_oran_rain_from_raw]` diagnostic line and add the column name to `_oran_rain_from_raw` |
| τ fit fails (`at-boundary`) | window too narrow / canopy not in decay regime | use `window_type=season` (default) rather than `peak` |
