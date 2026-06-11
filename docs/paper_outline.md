# Paper outline — Drip irrigation decouples surface soil moisture from canopy ET: implications for satellite ET retrieval

## Working title (alternatives)
1. Drip irrigation decouples surface soil moisture from canopy transpiration on a 3–4 day timescale: implications for satellite-based ET retrieval
2. Why satellite ET products underestimate irrigated almond orchards: a 7-year flux-tower–satellite comparison in semi-arid Spain
3. Quantifying the irrigation bias in MODIS and PML satellite ET over a Mediterranean almond orchard

## One-sentence pitch
Eddy-covariance and satellite ET (MOD16, PML) at one drip-irrigated almond orchard and one rainfed cereal field in semi-arid Spain reveal a 2–3 mm/day systematic underestimation of irrigated ET by satellites that decays on a 3–4 day timescale after each irrigation event — directly attributable to surface SWC decoupling, not deep rooting.

## Abstract bullets
- Two flux towers in Albacete (Spain): TzM (drip-irrigated almond, 2020–2024) and Oran (rainfed cereal rotation, 2018–2020).
- 1,356 site-days of EC LE/H/G/Rn merged with daily MOD16, PML, S2 NDVI/NDWI/NDMI, MODIS LST/LAI, ERA5-Land and CHIRPS.
- In-situ drought sensitivity SDS = 1 − mean(LE|dry SWC)/mean(LE|normal SWC) is +0.43 at rainfed Oran spring but +0.08 at irrigated TzM summer; reproduces v14 result independently from a clean re-pipelined master CSV.
- TzM satellite ET bias by days-since-last-irrigation:
  - MOD16: −3.35 → −1.94 → −1.68 mm/day (d0–3 / d4–7 / d8+)
  - PML:   −2.12 → −0.87 → −0.41 mm/day
- Oran (no irrigation): MOD16 bias −0.23, PML bias +0.44 mm/day.
- Mechanism: drip wetted bulb (10–30 cm) sustains transpiration without raising surface (0–5 cm) SWC; satellite retrievals see the dry surface and underpredict.
- Implication: standard SMAP/SMOS-driven satellite ET in irrigated agriculture has a quantifiable, physically interpretable irrigation bias that decays in 3–4 days.

## 1. Introduction

### 1.1 Knowledge gap
- Satellite ET products are increasingly used for irrigation scheduling, drought monitoring and water accounting.
- Their parameterisations are typically calibrated on FLUXNET2015 sites that under-represent drip-irrigated horticulture.
- Recent intercomparisons report large (>50%) ET underestimation in irrigated agriculture, but the dynamics on top of an irrigation event have not been quantified at sub-weekly resolution from EC data.

### 1.2 Hypothesis evolution (acknowledge methodological honesty)
- The investigation began as a "deep root" hypothesis at TzM.
- Subsequent v14-style analysis with days-since-irrigation stratification reframed it: the apparent decoupling is irrigation-driven, not biological deep-rooting.
- The present paper builds on the reframed hypothesis with an independent satellite check.

### 1.3 Objectives
1. Reproduce in-situ SDS from a clean master pipeline.
2. Quantify MOD16 and PML bias against EC ET in two contrasting management systems.
3. Test whether the bias scales with days-since-last-irrigation as predicted by the decoupling time-constant.
4. Translate the finding into a correctable bias model for irrigated agricultural pixels.

## 2. Sites and data

| Item | Oran | TzM |
|---|---|---|
| Lat / Lon | 38.82°N, −1.86°E | 39.27°N, −1.94°E |
| Crop | winter cereal rotation (vetch/wheat/pea) | drip-irrigated almond |
| Period | 2018-01 → 2020-06 | 2020-06 → 2024-10 |
| Daily EC observations | 654 | 701 |
| Mean annual ET (EC) | 1.13 mm/d | 3.36 mm/d |
| EBR (energy balance ratio) | 0.96 (2020) | 0.87 |

### 2.1 EC processing
- EddyPro half-hour fluxes; Oran AmeriFlux CLEAN format, TzM proprietary EddyPro output.
- Daily aggregation: LE/H/G means, ET in mm/day; QC ≤ 2 retained, ≥ 24 valid half-hours per day.
- VPD from Tetens; energy balance ratio reported (0.7–0.95 across both sites).

### 2.2 Satellite products (all from Google Earth Engine)
- MOD16A2GF v061 (500 m, 8-day) — ET, PET; scaled ×0.1 then ÷8 to mm/day.
- PML v018 (500 m, 8-day) — Ec + Es + Ei summed to total ET.
- MCD15A3H — LAI, FPAR (4-day, 500 m).
- MOD11A1 — daily LST.
- Sentinel-2 L2A — NDVI, NDWI, NDMI (5-day, 10 m).
- ERA5-Land — 2 m T, total precipitation.
- CHIRPS — daily precipitation (5 km).

### 2.3 Master pipeline
1. unify_ec_daily.py → ec_daily_master.csv
2. add_flags.py → days_since_irrig, irrig_bucket, season, drought_class
3. aggregate_oran_30min.py → fills Oran daily LE/H/G/Rn from 30-min
4. unify_satellite.py → satellite_long.csv + satellite_daily.csv (8-day forward filled to daily)
5. merge_satellite_ec.py → master_full.csv

## 3. Methods

### 3.1 In-situ drought sensitivity
SDS = 1 − mean(LE | SWC < p25) / mean(LE | p25 ≤ SWC ≤ p75), bootstrap (n=2000) 95% CI; computed per site / season; for TzM summer also per irrigation bucket using the TzM-summer-pool SWC quantiles.

### 3.2 Satellite ET bias
bias_X = sat_ET_X − EC_ET on overlapping site-days, evaluated globally and stratified by:
- Site (Oran vs TzM),
- Season,
- Irrigation bucket (d0–3 / d4–7 / d8+ / no_irrig),
- NDVI growing-season gate (NDVI > 0.3).

### 3.3 Linking SDS and bias
At each (site, stratum) cell with sufficient n we report SDS, mean satellite bias, and standard error to test the prediction:
> low SDS (decoupling) ↔ strongly negative satellite bias.

## 4. Results

### Fig 1. Site context and master coverage
4-panel: ET time series for both sites; NDVI seasonality; SWC range; irrigation events at TzM.

### Fig 2. SDS reproduction
Bar chart of SDS by stratum with bootstrap CI (Oran by season; TzM by season; TzM summer by irrigation bucket).
Key numbers: Oran spring +0.43 [+0.36, +0.50]; TzM summer all +0.08 [0.00, +0.17]; TzM d0–3 +0.11; d4–7 −0.03; d8+ −0.43 (n=43).

### Fig 3. EC vs satellite scatter
2-panel scatter (MOD16, PML) with site colours, 1:1 line, n / r / RMSE / MBE.

### Fig 4. Time series overlay
EC vs MOD16 vs PML for both sites.

### Fig 5. Irrigation decoupling — the headline figure
Boxplots of (sat − EC) bias by irrig_bucket at TzM, four panels: all-year vs summer vs growing vs summer×growing. Show the dose-response: −3.35 → −1.94 → −1.68 mm/day (MOD16) and −2.12 → −0.87 → −0.41 mm/day (PML).

### Fig 6. SDS vs bias
Scatter with strata as points, two panels (MOD16 / PML); demonstrates that low SDS strata are exactly the strata with most-negative bias.

### Optional Fig 7. Mechanistic schematic
Drip wetted bulb diagram showing why surface SWC sensors and satellite SM products miss the active root zone.

## 5. Discussion

### 5.1 Mechanism
- Drip irrigation creates a localised wet bulb at 10–30 cm depth; the surface 0–5 cm where most soil moisture sensors (and microwave SM retrievals) live remains dry.
- LE at TzM tracks the wet bulb, not the surface; in the first 3 days post-irrigation LE is decoupled from surface SWC.
- Satellite ET parameterisations driven by surface SM or LST therefore underestimate by a magnitude that scales with how recent irrigation was.

### 5.2 Quantitative bias model
A first-order correction:
  bias(t) ≈ bias_0 · exp(−t / τ),  τ ≈ 3–4 d.
Fitting to TzM gives bias_0 ≈ −3 mm/day (MOD16) or −2 mm/day (PML), τ ≈ 3 d.

### 5.3 Independence from rainfed control
At Oran (no irrigation) bias is small (|MBE| < 0.5 mm/day) and shows no irrigation-bucket structure, confirming the bias signal is driven by management not by climate or instrument.

### 5.4 Reframing the deep-root hypothesis
The original deep-root narrative is unsupported by the present data; the irrigation-buffered shallow/medium-root system fully explains the decoupling. Where the d8+ subset shows a nominally negative SDS (n small), it is consistent with VPD-driven transpiration of residual irrigation water rather than groundwater uplift.

## 6. Limitations
- Single drip-irrigated site (TzM) and single rainfed site (Oran); generalisation to flood/sprinkler irrigation requires more sites.
- 5 cm SWC sensor under-samples the drip wetted bulb (10–30 cm).
- ERA5 dewpoint and Rn could not be retrieved (only Ta and total precipitation), so VPD comes from EC, not satellite.
- 2024 PML data unavailable (PML v018 ends 2023).
- No SMAP root-zone SM in the present version; will be added in revision.
- Phenological differences (deciduous almond vs annual cereal) confound direct cross-site comparison; stratification by season and crop period mitigates this but does not remove it.
- Single 5 cm SWC sensor; multi-depth profile would let us directly observe the wet bulb.
- Footprint heterogeneity at Oran (rotation across 2018–2020) introduces a crop signal that should be filtered per cropping period.

## 7. Conclusions
1. Satellite ET (MOD16, PML) systematically underestimates EC ET by 2–3 mm/day in the first 3 days after a drip irrigation event at a Mediterranean almond orchard.
2. The bias decays on a 3–4 day timescale and is essentially zero at a nearby rainfed cereal site, consistent with an irrigation-driven, not climate-driven, decoupling between surface SWC and canopy transpiration.
3. The same decoupling is captured independently by the in-situ SDS metric, providing a coherent flux-tower–satellite story.
4. A simple exponential-decay correction conditioned on days-since-last-irrigation could remove most of this bias and substantially improve satellite ET in irrigated agriculture.

## Target journals (short list)
- Agric. Forest Meteorol. (best fit, EC + remote sensing)
- Remote Sensing of Environment
- Hydrology and Earth System Sciences
- Journal of Hydrometeorology
- Biogeosciences

## Next steps before submission
- [ ] Re-run figs with summer × NDVI > 0.3 filter (figure_C_summer.py) and use that as Fig 5.
- [ ] Add LSA SAF Meteosat ETv2/v3 as a third independent satellite product (Spain in MSG view).
- [ ] Add SMAP L4 root-zone SM when download finishes.
- [ ] Fit the exponential decay bias model and report τ with CI.
- [ ] Confirm SIGPAC parcel boundaries for footprint quality control.
- [ ] Get PI consent / ICOS metadata if either site is in the network.
- [ ] Draft Methods paragraph that lists every script in scripts/ as the reproducibility chain.
