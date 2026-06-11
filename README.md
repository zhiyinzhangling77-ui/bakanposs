# bakanposs

Drip irrigation × satellite ET analysis at two flux tower sites in semi-arid
southeastern Spain (Albacete province):

| Site | Location | Type | EC period |
|---|---|---|---|
| **Oran** | 38.82°N, 1.86°W | rainfed cereal (vetch / wheat / pea rotation) | 2018-01 – 2020-06 |
| **Tarazona de la Mancha (TzM)** | 39.27°N, 1.94°W | drip-irrigated almond orchard | 2020-06 – 2024-10 |

Central scientific claim:
**Satellite ET (MOD16, PML, METv3) systematically underestimates drip-irrigated
almond ET by 2.8–4.1 mm d⁻¹ in the 3 days following irrigation, then relaxes
on τ ≈ 4–6 d. This event-locked structure is fully explained by the drip
wetted-bulb mechanism (Skaggs et al., 2004; Cote et al., 2003) and is
correctable with a single exponential decay term conditioned on
days_since_irrig (RMSE −49 to −65 %).**

## Repository layout

```
bakanposs/
├── analysis_A/                  EC-only τ recovery (LE-pulse fit)
│   └── analysis_A_v9.py
│
├── analysis_B/                  Satellite τ verification (5 attempts)
│   ├── analysis_B_v1_mod16_tau.py        MOD16A2 8-day (fails: 8-day blur)
│   ├── analysis_B_v2_metv3_tau.py        METv3 daily   (fails: night-zero dilution)
│   ├── analysis_B_v3_bias_tau.py         EC−METv3 bias (SUCCESS: τ_TzM = 4.57 d)
│   ├── analysis_B_v4_metv3_30min_tau.py  METv3 30-min  (fails: diurnal swamps signal)
│   └── analysis_B_v5_oran_lst_response.py  Oran rain-LST response probe
│
├── analysis_C_v1.py             NDVI / phenology / drought-class analysis
│                                (imports analysis_A_v9 + data_loaders)
├── data_loaders.py              Shared EC daily loaders (Oran / TzM)
│
├── pipeline/                    one-shot CSV/parquet builders
│   ├── unify_ec_daily.py        EC half-hour → daily master
│   ├── aggregate_oran_30min.py  Oran 30-min recovery of LE/H/G
│   ├── add_flags.py             season / irrig bucket / drought class flags
│   ├── unify_satellite.py       GEE wide-CSVs → daily long table
│   ├── load_metv3.py            ~120 k NetCDF → daily METv3 (≈3.5 h)
│   ├── load_smap.py             SMAP 3 h CSV → daily
│   ├── merge_satellite_ec.py    EC × satellite → master_full.csv
│   ├── integrate_metv3_smap.py  METv3 + SMAP join → master_full_v2.csv
│   └── qc_master.py             QC report on master
│
├── figures/                     figure / hypothesis-test scripts
│   ├── tau_fit.py               3-product exp-decay (Δ = a·exp(−t/τ) + c)
│   ├── sds_vs_bias.py           SDS × bias scatter
│   ├── sds_v14_repro.py         in-situ SDS re-derivation
│   ├── figure_C_summer.py       summer × irrig-bucket boxplots
│   ├── hypothesis_tests.py      H1 / H4 / H6 verifications
│   ├── killer_figures.py        narrative-driving figures
│   └── first_plots.py           early-exploration plots
│
├── inspectors/                  one-off data exploration
│   ├── inspect_ec.py
│   ├── inspect_satellite.py
│   ├── inspect_metv3.py
│   └── inspect_smap.py
│
├── gee/                         Google Earth Engine JavaScript
│   ├── gee_extract.js           7-product MODIS/PML/S2/LAI/LST/ERA5/CHIRPS
│   └── gee_smap_only.js         SMAP-only (needs 6 km buffer)
│
└── docs/
    ├── paper_methods_results.md   paper Intro/Methods/Results/Discussion draft
    ├── paper_outline.md           target-journal outline
    ├── analysis_narrative.md      full narrative + statistics walkthrough
    ├── SESSION_SUMMARY.md         cross-session handoff
    └── reports/
        ├── analysis_C_report.md
        └── migration_to_data_loaders.md
```

## Data sources (NOT in repo; absolute paths in scripts)

| Dataset | Path | Cadence | Notes |
|---|---|---|---|
| EC raw | `~/Dataset/Eddy data in Spain/` | half-hour + daily | Oran/TzM CSVs |
| METv3 | `/mnt/hdd/Dataset/METv3/YYYY/MM/MMDD/*.nc` | 30 min, 0.05° | ~120 k files; `pipeline/load_metv3.py` aggregates to daily |
| SMAP L4 | `/mnt/hdd/Dataset/SMAP_OranTzM.csv` | 3 h, 9 km | `pipeline/load_smap.py`; requires **6 km buffer** in GEE |
| MOD16A2 | `/mnt/hdd/Dataset/MOD16A2_2018_2024_Oran_TzM/*.csv` | 8-day, 500 m | NASA AppEEARS export; scale already applied |
| GEE wide-CSVs | `/mnt/hdd/Dataset/Fast_OranTzM_*.csv` | 8-day, 500 m–10 m | MODIS, PML, S2, LAI, LST, ERA5, CHIRPS |

## Reproduction pipeline (cold start)

```bash
# 1. data ingestion
python3 pipeline/unify_ec_daily.py
python3 pipeline/add_flags.py
python3 pipeline/aggregate_oran_30min.py
python3 pipeline/unify_satellite.py
python3 pipeline/load_metv3.py        # 3.5 h (skip if metv3_daily_all.csv exists)
python3 pipeline/load_smap.py
python3 pipeline/merge_satellite_ec.py
python3 pipeline/integrate_metv3_smap.py
# → produces master_full_v2.csv (53 cols × 1356 site-days)

# 2. EC-only τ
python3 analysis_A/analysis_A_v9.py
# → τ_TzM = 3.36 d, τ_Oran = 2.82 d (active-pool)

# 3. Satellite τ verification (run any subset)
python3 analysis_B/analysis_B_v1_mod16_tau.py --quick
python3 analysis_B/analysis_B_v3_bias_tau.py  --quick   # the only one that succeeds

# 4. NDVI / phenology
python3 analysis_C_v1.py

# 5. Hypothesis tests + figures
python3 figures/tau_fit.py
python3 figures/sds_vs_bias.py
python3 figures/hypothesis_tests.py
```

## Branch convention

- `main` — submitted state
- `claude/compare-ec-satellite-et-ZnENi` — primary development branch
- Other `claude/*` branches — session-specific work; merge via PR when stable

## Key numerical results (current state)

| Quantity | Value | Source |
|---|---|---|
| τ (EC, TzM) | 3.36 d (SE 0.62) | `analysis_A_v9.py` v27 |
| τ (EC − METv3 bias, TzM) | 4.57 d [3.03, 13.10] | `analysis_B_v3` |
| τ (full bias model, METv3) | 6.0 d [4.6, 8.3] | `figures/tau_fit.py` |
| MBE_MOD16 (TzM) | −2.69 mm d⁻¹ | `figures/sds_vs_bias.py` |
| MBE_METv3 (TzM) | −2.34 mm d⁻¹ | `figures/sds_vs_bias.py` |
| τ-corrected RMSE reduction | −49 % to −65 % | `figures/hypothesis_tests.py` (H1) |
| ΔAIC days_since_irrig vs VPD | −66 to −153 | `figures/hypothesis_tests.py` (H4) |
| r(SWC_5cm, SMAP_rz) Oran spring | +0.80 (n=203) | `figures/hypothesis_tests.py` (H6) |
| r(SWC_5cm, SMAP_rz) TzM d0–3 | −0.19 (n=281) | `figures/hypothesis_tests.py` (H6) |

See `docs/paper_methods_results.md` for full results and `docs/analysis_narrative.md`
for methodological narrative.
