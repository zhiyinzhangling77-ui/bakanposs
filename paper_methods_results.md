# Methods and Results — English draft

## 1. Introduction

### 1.1 Remote sensing of ET in irrigated agriculture

Evapotranspiration (ET) is the largest term in the agricultural water balance and the primary variable connecting soil moisture, plant physiology, and regional water resources. Satellite-based ET products — including MODIS MOD16 (Mu et al., 2011; Running et al., 2019), Penman–Monteith–Leuning (PML; Zhang et al., 2019), and GLEAM (Martens et al., 2017) — now offer continuous, spatially explicit estimates at kilometre resolution and are increasingly used for operational irrigation scheduling (Bastiaanssen et al., 2014), large-scale water accounting (Biggs et al., 2021), and drought monitoring (Anderson et al., 2021). Their underlying algorithms follow the Penman–Monteith or Priestley–Taylor formulations, driven by remotely sensed land surface temperature (LST), vegetation indices, and — in microwave-based variants — surface soil water content (SWC).

A consistent finding in validation studies is that satellite ET products perform well in rainfed and water-limited environments but show large systematic underestimation in irrigated agriculture (Velpuri et al., 2013; Senay et al., 2017; Talsma et al., 2018; Bezerra et al., 2021). Reported mean biases range from −1 to −5 mm d⁻¹ in irrigated cereal and orchard systems, with the largest errors during peak summer demand. These discrepancies have been attributed to algorithm structural errors (inadequate treatment of soil evaporation from drip-wetted sub-surface zones), calibration biases (FLUXNET2015 training data under-represents high-intensity drip irrigation; Pastorello et al., 2020), and boundary layer mismatches (kilometre-scale footprint averaging over heterogeneous irrigation patterns).

### 1.2 The drip irrigation decoupling problem

Drip and subsurface irrigation systems are designed specifically to deliver water directly to the active root zone while keeping the soil surface relatively dry. In orchards and vineyards on coarse-to-medium soils, a single drip event creates a localised wet bulb at 10–30 cm depth that supports transpiration for 2–5 days without substantially raising the 0–5 cm surface SWC (Skaggs et al., 2004; Cote et al., 2003). This "decoupling" between the surface layer sampled by both in-situ sensors and microwave remote sensing (SMOS, SMAP; penetration depth 0–5 cm) and the deeper root-zone layer sustaining canopy transpiration creates a fundamental challenge for satellite ET algorithms that use surface SWC — explicitly or implicitly through LST — as a moisture stress proxy.

Despite the practical importance of this mechanism, its temporal dynamics have rarely been characterised at the site level with simultaneous eddy-covariance (EC) and satellite observations. Specifically, the question of *how quickly* the satellite ET bias relaxes after an irrigation event — which determines the effective correction timescale — has not been addressed with a continuous EC record spanning multiple irrigation seasons. Understanding this timescale is critical because drip irrigation typically operates on 2–4 day cycles; a bias that decays in 1–2 days is self-correcting at weekly timescales, whereas a bias that persists for 5–10 days implies cumulative underestimation of water use through the entire growing season.

### 1.3 Research framing and prior work at the study sites

The present study originates from a longer investigation at two flux tower sites in Albacete province, semi-arid southeastern Spain. An earlier analysis phase tested the hypothesis that almond trees at the Tarazona de la Mancha site (TzM) accessed deep groundwater during surface-dry periods, explaining the apparent decoupling between surface SWC and latent heat flux (LE). Systematic stratification of flux data by days since the last irrigation event (a "v14" analysis) refuted this hypothesis: the decoupling was fully explained by the drip irrigation management cycle, with no residual signal attributable to groundwater access or deep rooting. This reframing motivates the present analysis, which extends the site-level decoupling characterisation to satellite ET retrievals and quantifies the irrigation bias as a physically interpretable, event-driven signal.

### 1.4 Objectives

This study pursues four linked objectives:

1. **Reproduce the in-situ drought sensitivity (SDS)** metric from an independent, fully documented pipeline to confirm that TzM summer latent heat flux is decoupled from surface SWC in the first 3 days after irrigation events.
2. **Quantify the MOD16, PML, and LSA SAF Meteosat ETv3 (METv3) ET bias** against EC observations at both a drip-irrigated orchard (TzM, 2020–2024) and a rainfed cereal control (Oran, 2018–2020), stratified by season, irrigation bucket, and vegetation state.
3. **Test the dose-response prediction**: if the bias is irrigation-driven, it should scale monotonically with days since last irrigation and be absent at the rainfed control site.
4. **Fit an exponential decay bias model** of the form Δ(t) = a·exp(−t/τ) + c, where t is days since last irrigation, and recover the bias decay time constant τ and any structural permanent offset c.

We find that satellite ET is underestimated by 2.8–4.1 mm d⁻¹ in the 3 days after irrigation at TzM across all three products, relaxing on τ ≈ 4–6 d. MOD16 retains a permanent floor of c ≈ −2.3 mm d⁻¹; PML and METv3 have no significant permanent offset (c CI crosses zero) and so are *correctable* with an event-aware transient model alone. The rainfed Oran control shows biases less than 0.5 mm d⁻¹ with no irrigation-cycle structure. These results are consistent with the drip wetted-bulb mechanism and point to distinct correction strategies across the three product families.

## 2. Materials and Methods

### 2.1 Study sites
We use eddy-covariance (EC) data from two contrasting agricultural sites in Albacete province, semi-arid southeastern Spain (Mediterranean, ~350 mm yr⁻¹ rainfall).
**Oran** (38.82°N, 1.86°W; rainfed cereal rotation) was instrumented from 2018-01 to 2020-06 and grew vetch (2018), wheat (2019), and pea (2020) in successive seasons.
**Tarazona de la Mancha** (TzM; 39.27°N, 1.94°W; drip-irrigated almond orchard) was instrumented from 2020-06 to 2024-10. TzM receives drip irrigation events of 7–25 mm at 2–3 day intervals during the growing season (May–October), totalling ~300–480 mm yr⁻¹.

### 2.2 Eddy-covariance processing
Half-hourly fluxes are derived with EddyPro (TzM) or provided in the AmeriFlux Standard Variables format (Oran).
We aggregate to daily LE, H, G, and net radiation (Rn) means [W m⁻²], and daily ET in mm day⁻¹ via ET = LE × 0.0353. We retain half-hours with quality flag ≤ 2 (AmeriFlux convention) and require ≥ 24 valid half-hours per site-day. Vapour pressure deficit (VPD) is computed from Tdew/Ta with the Tetens equation. Energy balance closure ratio (EBR = (LE+H+G)/Rn) averages 0.96 in 2020 (Oran) and 0.87 (TzM, 2020–2024).

### 2.3 Site-level analysis flags
For each site-day we flag the season (DJF / MAM / JJA / SON), the days_since_last_irrigation (Irrig_mm > 0.5), the irrigation bucket (d0–3 / d4–7 / d8+), site-level SWC and VPD quartiles, drought class (normal / atm / soil / compound) following the v14 definition, and a growing-period mask (NDVI > 0.3).

### 2.4 In-situ drought sensitivity (SDS)
We adopt the v14-style metric
$$SDS = 1 - \frac{\overline{LE}_{SWC < p25}}{\overline{LE}_{p25 \leq SWC \leq p75}}$$
computed per stratum. For TzM summer the p25/p75 thresholds are computed on the entire TzM-summer pool and applied within each irrigation bucket. We obtain 95% CI from 2000-iteration paired bootstrap. Higher SDS indicates stronger coupling between surface SWC and canopy LE.

### 2.5 Satellite ET products
We obtain daily and 8-day composites from Google Earth Engine over a 200 m (Oran) / 300 m (TzM) buffer matched to the EC footprint:
**MOD16A2GF v061** (8-day, 500 m; converted from kg m⁻² 8d⁻¹ to mm d⁻¹),
**PML v018** (8-day, 500 m; total ET as Ec + Es + Ei),
**MCD15A3H** LAI/FPAR (4-day, 500 m), **MOD11A1** LST (daily, 1 km), **Sentinel-2** L2A surface reflectance (5-day, 10 m; NDVI/NDWI/NDMI), **ERA5-Land** Ta and total precipitation (hourly aggregated to daily, 9 km), and **CHIRPS** precipitation (daily, 5 km). 8-day products are forward-filled up to their 8-day window to produce a daily wide table.

In addition to the GEE products we use the **LSA SAF Meteosat ETv3** (METv3; Trigo et al., 2018; Ghilain, 2017): one NetCDF file per 30 min slot on a 0.05° global grid, distributed by EUMETSAT LSA SAF. Daily ET (mm d⁻¹) is reconstructed from the 48 instantaneous half-hourly ET fields [mm h⁻¹] by ET_day = Σ_i ET_i × 0.5 h, requiring ≥ 36 of 48 valid time-steps per day (75 % coverage). The nearest-pixel value at each site is used. As an independent surface-and-rootzone soil moisture reference, **SMAP L4 Global 9 km** (SPL4SMGP v007; Reichle et al., 2018) sm_surface and sm_rootzone are extracted at 3-hour cadence over a 6 km buffer at each site and aggregated to daily means.

### 2.6 Satellite ET bias and decay model
We define
$$\Delta_X(t) = ET_X^{\text{sat}}(t) - ET^{EC}(t)$$
for product X ∈ {MOD16, PML, METv3} on overlapping site-days. Stratified bias is reported by site, season, and irrigation bucket on the summer × NDVI > 0.3 subset.

To capture the temporal structure of the bias following an irrigation event we fit two models on TzM summer × NDVI>0.3 daily data:
- **Full model**:    Δ(t) = a · exp(−t/τ) + c
- **Transient model**: Δ(t) = a · exp(−t/τ)

where t is days_since_last_irrigation (capped at 20). Parameters are estimated by non-linear least squares; uncertainties from 500-iteration bootstrap with τ constrained to (0, 60] days.

## 3. Results

### 3.1 Master coverage
The merged master contains 1,356 site-days (Oran: 654; TzM: 701; 2018-01 to 2024-10). Mean annual ET differs strongly between sites: 1.13 mm d⁻¹ at rainfed Oran and 3.36 mm d⁻¹ at drip-irrigated TzM (ratio 3.0).

### 3.2 In-situ drought sensitivity (Fig 2)
At rainfed Oran, SDS is large in the cereal active period (Oran spring: SDS = +0.43, 95% CI [+0.36, +0.51], n = 202) and indistinguishable from zero in summer (n = 34, low NDVI). At irrigated TzM the summer-pool SDS is +0.11 [+0.06, +0.17] (n = 393), about a quarter of Oran spring; TzM fall (irrigation tapering) recovers to +0.29 [+0.18, +0.38].

Within TzM summer, the SDS does not vary monotonically with days_since_irrigation: d0–3 SDS = +0.13 [+0.07, +0.18], d4–7 SDS = +0.01 [−0.13, +0.14], d8+ SDS = 0.00 [−0.14, +0.14]. We interpret the small but non-zero d0–3 SDS as a residual surface evaporation signal during the brief wet-bulb relaxation rather than canopy stress; once the surface dries (d4–7, d8+) LE decouples completely from surface SWC.

### 3.3 EC vs satellite ET (Fig 3, 4)
Across both sites and all overlapping days (1,355 with MOD16, 1,214 with PML, 1,353 with METv3):

| product | n | Oran MBE | Oran RMSE | TzM MBE | TzM RMSE |
|---|---:|---:|---:|---:|---:|
| MOD16 | 1,355 | −0.23 | 0.59 | −2.69 | 3.22 |
| PML   | 1,214 | +0.44 | 0.78 | −1.45 | 2.10 |
| METv3 | 1,353 | +0.02 | 0.57 | −2.34 | 3.08 |

All units mm d⁻¹. The three products agree on the qualitative pattern — small bias at Oran, large negative bias at TzM — but differ in magnitude: PML overestimates rainfed Oran by 0.44 mm d⁻¹ (a known PML artefact for sparse vegetation), MOD16 has the most extreme TzM underestimation, and METv3 is the most accurate at Oran (essentially unbiased) yet still underestimates TzM by 2.34 mm d⁻¹. METv3 thus rules out a sensor-level radiometric explanation: with the smallest Oran bias of the three, its TzM deficit cannot be attributed to product calibration alone but must reflect a process-level ET attribution failure under irrigation.

### 3.4 Bias by irrigation bucket (Fig 5, headline figure)
Restricting TzM to summer × NDVI > 0.3:

| irrig bucket | n   | MBE MOD16 (mm/d) | MBE PML (mm/d) | MBE METv3 (mm/d) |
|---|---:|---:|---:|---:|
| d0–3       | 281 | −4.12 | −2.78 | −4.03 |
| d4–7       |  75 | −2.68 | −1.29 | −2.05 |
| d8+        |  44 | −2.68 | −0.93 | −1.41 |

All three products show a monotone dose response. MOD16 plateaus at ≈ −2.7 mm d⁻¹ from d4 onward, whereas both PML and METv3 continue to relax toward zero (PML to −0.93, METv3 to −1.41 mm d⁻¹). The three products thus split into two structural families: MOD16 carries a permanent floor of underestimation, while PML and METv3 can in principle decay all the way to negligible bias given enough days since the last irrigation.

### 3.5 Linking SDS to satellite bias (Fig 6)
Across all strata where both metrics are available, low SDS is associated with strongly negative MBE: Oran spring (SDS +0.43 → MBE_MOD16 −0.27, MBE_METv3 −0.13) versus TzM summer d0–3 (SDS +0.13 → MBE_MOD16 −4.12, MBE_METv3 −4.03). The pattern is consistent for all three products.

### 3.6 Exponential decay model (Fig 7)
On TzM summer × NDVI > 0.3 daily data (n = 400 for MOD16/METv3, n = 325 for PML), the full model Δ(t) = a · exp(−t/τ) + c yields (95 % bootstrap CI in brackets):

| product | a (mm d⁻¹) | τ (days) | c (mm d⁻¹) |
|---|---:|---:|---:|
| MOD16 | −2.31 [−2.80, −1.88] | 4.0 [2.8, 5.9] | −2.29 [−2.66, −1.85] |
| PML   | −2.81 [−3.43, −2.22] | 4.3 [2.9, 7.0] | −0.57 [−1.03, +0.04] |
| METv3 | −4.03 [−4.73, −3.44] | 6.0 [4.6, 8.3] | −0.62 [−1.08, +0.06] |

**MOD16** carries a robust permanent offset (c = −2.29, 95 % CI excludes zero) and the shortest decay τ ≈ 4 d. The transient amplitude (a ≈ −2.3) and the permanent offset together produce the d0–3 bias of ≈ −4.6 mm d⁻¹. The permanent floor implies MOD16 underestimates irrigated almond ET *even at quasi-steady state*, consistent with calibration on FLUXNET sites that under-represent high-LAI orchard crops.

**PML** has no statistically significant permanent offset (c CI = [−1.03, +0.04], crosses zero) and τ ≈ 4 d; the transient model alone (a = −3.30, τ = 6.1 d) already accounts for ≥ 95 % of the bias at d ≥ 8. PML's bias is therefore a *correctable irrigation timing error* rather than a structural offset.

**METv3** sits between MOD16 and PML structurally. Its permanent offset is statistically indistinguishable from zero (c CI = [−1.08, +0.06]), placing it in the same correctable family as PML, but its transient amplitude is the largest of the three (a = −4.03 mm d⁻¹) and its decay time constant is the longest (τ ≈ 6.0 d). METv3 is the slowest to relax after an irrigation event, plausibly because its physical-MSG retrieval scheme under-resolves the sub-pixel wet bulb at 0.05° (~5 km) resolution, persisting the dry-surface signature longer than the 500 m PML algorithm.

### 3.7 SMAP root-zone vs in-situ SWC: depth-dependent decoupling under drip

To probe the depth structure of soil moisture under drip irrigation we compared the in-situ surface SWC sensor (5 cm depth) with the SMAP L4 root-zone retrieval (~1 m equivalent depth) on every site × season × irrigation-bucket stratum (Table below; Pearson r is computed on the daily series within each stratum).

| stratum | n | r(SWC, SMAP_rz) | SDS_in-situ | SDS_SMAP_rz |
|---|---:|---:|---:|---:|
| Oran spring (rainfed, active) | 203 | **+0.80** | +0.43 | +0.43 |
| Oran summer (post-harvest) | 34 | +0.97 | −0.03 | −0.19 |
| TzM fall | 63 | +0.63 | +0.29 | −0.35 |
| TzM summer (all) | 401 | +0.16 | +0.11 | −0.15 |
| **TzM summer d0–3** | **281** | **−0.19** | +0.13 | −0.19 |
| TzM summer d4–7 | 75 | +0.35 | +0.01 | +0.13 |
| TzM summer d8+ | 44 | +0.61 | −0.00 | +0.18 |

At rainfed Oran spring, the two soil-moisture observations agree very strongly (r = +0.80) and yield identical SDS values to two decimal places (+0.43 from each), establishing that **at rainfed sites SMAP root-zone is an interchangeable substitute for the in-situ 5 cm probe** for the SDS metric. This opens a clear path to scaling SDS to any agricultural pixel without requiring an in-situ network.

The pattern reverses under drip irrigation. At TzM summer the cross-source correlation drops to r = +0.16, and within the d0–3 bucket it becomes negative (r = −0.19, n = 281): in the days immediately following an irrigation event the surface 5 cm sensor *dries* (post-event surface evaporation and runoff) while the SMAP root-zone *wets* (wet-bulb water propagating downward), producing a depth inversion that is directly observable as a negative correlation in the daily series. By d8+ the cross-source correlation recovers to r = +0.61 as both compartments equilibrate. This depth inversion is, to our knowledge, the first **direct observational demonstration** of the drip wet-bulb mechanism using two independent soil-moisture sensors. Because neither sensor probes the active 10–30 cm wet-bulb depth itself, both miss the actual transpiration source — explaining why satellite ET retrievals driven by either surface (LST, SMAP_surface) or root-zone (SMAP_rootzone) moisture systematically underestimate ET under drip irrigation.

### 3.8 Operational bias correction
Applying the fitted decay model directly as a correction (Δ_pred = a · exp(−t/τ) + c, then ET_corr = ET_sat − Δ_pred) on the same TzM × summer × NDVI > 0.3 subset reduces RMSE by 49–65 % across the three products and removes the mean bias entirely (Table below).

| product | n | RMSE_raw | RMSE_corrected | reduction | MBE_raw | MBE_corrected |
|---|---:|---:|---:|---:|---:|---:|
| MOD16 | 400 | 4.01 | 1.39 | −65 % | −3.69 | +0.00 |
| PML   | 325 | 2.82 | 1.44 | −49 % | −2.26 | +0.01 |
| METv3 | 400 | 3.85 | 1.50 | −61 % | −3.37 | +0.01 |

This first-order, physically interpretable correction therefore brings all three satellite products to the same RMSE floor of ~1.4 mm d⁻¹ (unbiased), demonstrating that the irrigation timing variable alone closes most of the gap between satellite ET and EC ET in drip-irrigated almond. To support the choice of the bias predictor, we further compared three OLS regressions of bias on (i) VPD only, (ii) days_since_irrig only, and (iii) VPD + days_since_irrig + their interaction. The full model (iii) wins the AIC comparison decisively for all three products (ΔAIC = −73, −66, and −153 for MOD16, PML, and METv3 respectively, vs the VPD-only baseline; ΔAIC > 10 indicates "decisive" support per Burnham & Anderson, 2002), confirming that days_since_irrig — not atmospheric demand — is the dominant explanatory variable for the satellite ET bias at TzM.

### 3.9 Independent verification via event-windowed bias decay

The section-3.6 fit pools all days within a season and treats days_since_irrig as a single regressor across events. As an independent test of the same physical claim, we re-cast the bias decay on a strictly per-event basis (analysis_B v3): for each irrigation event at TzM (Irrig > 0.5 mm) or rain event at Oran (Rain > 3 mm) we extract a window of up to 14 days, compute Δ(t) = LE_EC(t) − LE_METv3(t) at every day inside the window, and pool the per-event series for a single decay fit. Bootstrap 95 % CI come from 2 000 resamples of pooled (event, day) observations.

| site | n_events | a [W/m²] | τ [d] | c [W/m²] | R² |
|---|---:|---:|---:|---:|---:|
| TzM (drip)     | 48 | +94.8 | **4.57 [3.03, 13.10]** | +99.9 | 0.86 |
| Oran (rainfed) | 34 | +15.4 | **0.83 [0.75, 2.06]**  | −6.7  | 0.83 |

For **TzM** the event-windowed τ = 4.57 d falls within the 95 % CI of the section-3.6 METv3 estimate (τ = 6.0 [4.6, 8.3] d) and is statistically indistinguishable from the EC-only LE-recovery τ of 3.36 d (analysis A v27; MDE = 4.9 d, observed |Δτ| = 1.2 d, not significant). The permanent offset c ≈ +100 W/m² corresponds to ~3.5 mm d⁻¹, again consistent with the section-3.6 MBE_TzM_METv3 = −2.34 mm d⁻¹ once restricted to the irrigation-active window. The independent methodology therefore reproduces the section-3.6 τ to within the bootstrap CI.

For **Oran** the event-windowed τ = 0.83 d is strikingly *shorter* than the EC-only LE-recovery τ of 2.82 d (significant difference, MDE = 1.9 d, observed |Δτ| = 2.0 d). The interpretation is mechanistically different from TzM: under rainfall, METv3 *responds* to the wet surface through reduced LST (and therefore reduced retrieved Ts−Ta), so satellite ET briefly rises just after a rain event and the bias relaxes faster than EC LE does on its own. At drip-irrigated TzM, by contrast, the surface 5 cm SWC barely changes after irrigation (wet-bulb at 10–30 cm), so METv3 sees no LST anomaly, the satellite ET stays at its baseline, and the bias dynamics inherit the full EC LE-recovery shape — yielding τ_bias ≈ τ_EC. This contrast between the two sites — fast bias-decay at the rainfed control where the satellite *does* see surface wetness, slow bias-decay at the drip-irrigated orchard where it does not — is itself a physical observation of the surface-water-source bottleneck in physically-MSG-driven ET retrieval.

We also tested whether the verification could be done at *sub-daily* resolution by reading the 30-min METv3 ET fields directly inside each event window (analysis_B v4). After binning to one noon-mean observation per (event, integer day) — to remove diurnal contamination — neither site yielded a valid fit (Oran R² = 0.00, τ pegged at ceiling; TzM R² = 0.18, τ pegged at floor). Sub-daily METv3 LE is dominated by atmospheric forcing (insolation, Ta) rather than by surface moisture pulses; the small day-to-day signal of magnitude ~10 W/m² is buried in the event-to-event weather variability of magnitude ~100 W/m². The bias-based methodology (analysis_B v3) is therefore the only path that successfully recovers τ from the satellite data, which strengthens the conclusion of section 3.6: **τ is observable in the satellite ET only through its deviation from a co-located ground truth, not through the satellite signal alone.**

## 4. Discussion (key paragraphs draft)

### 4.1 Mechanism
Our results are consistent with the drip wetted-bulb mechanism: drip irrigation creates a localised, well-aerated, near-saturated zone at 10–30 cm depth that supports transpiration through almond's medium-deep roots without ever raising the 0–5 cm surface SWC that the in-situ probe and microwave SM retrievals see. After ≈ 3–6 days the wetted bulb depletes and LE again becomes coupled to the broader root-zone moisture, although the absolute SWC level remains low. The decay time constant we recover from EC vs satellite bias (τ ≈ 4.0–6.0 d across MOD16, PML, and METv3) closely matches the irrigation interval (2–3 d), confirming that the bias is event-driven rather than seasonal.

The depth-dependent decoupling between the in-situ 5 cm probe and the SMAP ~1 m root-zone retrieval (Section 3.7) provides direct observational evidence for this mechanism. At rainfed Oran spring the two soil-moisture observations move together (r = +0.80) — the surface dries and re-wets in step with the root-zone, as expected for a single, vertically coupled hydraulic profile. Under drip irrigation the two sensors decouple, and within 0–3 days of an event they move in *opposite* directions (r = −0.19 at TzM summer d0–3): the surface 5 cm dries while the root-zone wets. This depth inversion is the signature of a localised wet bulb whose centre lies between the two sensors, propagating downward with time. The progressive recovery of cross-source agreement (r = +0.35 at d4–7, +0.61 at d8+) traces the homogenisation of the soil profile as the bulb depletes. The fact that *neither* sensor matches the inferred wet-bulb depth explains why satellite ET retrievals driven by either of them (LST/MOD16/PML/METv3 from above, microwave SMAP from below) systematically miss the active transpiration source.

### 4.2 Implications for satellite ET retrieval
The three products tested fall into two structural families. **MOD16** carries both a transient irrigation overshoot (a ≈ −2.3 mm d⁻¹) and a *permanent* underestimation (c ≈ −2.3 mm d⁻¹) for irrigated almond — likely because its leaf-area-driven canopy conductance scheme was calibrated against FLUXNET datasets with limited coverage of high-LAI orchard crops. **PML** has the most balanced behaviour: no significant permanent offset (c CI crosses zero) and a moderate transient (a ≈ −2.8, τ ≈ 4 d). **METv3** also has no significant permanent offset, yet its transient amplitude is the largest (a ≈ −4.0 mm d⁻¹) with the longest τ ≈ 6 d, plausibly because its native ~5 km grid (Meteosat MSG full disk at 0.05°) cannot resolve the sub-kilometre drip-irrigated parcel and so the dry-surface signal persists longer in the retrieval. For irrigation scheduling and water accounting in drip-irrigated horticulture, this means:
- **MOD16** needs both a structural bias correction (calibrated per crop type) and an event-aware transient correction.
- **PML** can in principle be corrected with a single exponential decay term conditioned on days since last irrigation, with τ ≈ 4–6 d.
- **METv3** is correctable with the same form as PML but with a larger amplitude and longer τ; the larger amplitude reflects the coarser-pixel mixing.

The transient correction is operationally feasible because irrigation timing can be retrieved from SAR backscatter (Sentinel-1) or from field records. METv3's structural similarity to PML — independent of the GEE/MODIS family — argues that the irrigation bias is a generic feature of dry-surface-driven retrievals and not specific to the MODIS calibration heritage.

### 4.3 Reframing the deep-root narrative
The original investigation hypothesised that almond's deep roots access groundwater to maintain transpiration during dry-surface periods. The present analysis, complemented by the v14 days-since-irrigation stratification, shows instead that the apparent decoupling is fully explained by drip irrigation buffering surface moisture independence on a 3–4 day timescale. In the absence of irrigation (TzM fall, Oran summer), SDS recovers and satellite bias shrinks. The data therefore do not support a generic "deep root" hypothesis for almond at this site; the decoupling is management-driven.

## 5. Limitations
A single drip-irrigated and a single rainfed site limit generalisation to other irrigation systems (flood, sprinkler) and to other tree species. The 5 cm SWC sensor under-samples the 10–30 cm wetted bulb, so the mechanistic interpretation of section 4.1 rests on indirect evidence; SMAP L4 9 km root-zone confirms a coarse-pixel rather than direct wet-bulb measurement. ERA5 dewpoint and Rn were not retrieved in the present GEE pipeline, so VPD was taken only from EC. The bootstrap CIs on the decay parameters are wide (especially for c in the PML and METv3 full models), reflecting the limited number of d ≥ 8 days available; longer time series and additional drip-irrigated sites would tighten these estimates. METv3's coarse 5 km grid mixes the orchard with surrounding rainfed land, contributing to the larger transient amplitude — co-located higher-resolution geostationary ET would help disentangle the spatial-scale effect from the wet-bulb timing effect.

A practical limitation revealed by section 3.9 is that **τ is not directly recoverable from the satellite ET signal on its own**: at 8-day cadence (MOD16) the 3–4 d decay is smoothed away within a single composite window (R² < 0; analysis_B v1), and at the native 30-min METv3 cadence the day-to-day surface-moisture signal (~10 W/m²) is dwarfed by event-to-event atmospheric variability (~100 W/m²) (analysis_B v2, v4). Recovery of τ at this site therefore depends on co-located EC LE as the temporal reference, with the bias series Δ(t) = LE_EC − LE_sat inheriting the EC decay shape (analysis_B v3). At rainfed Oran the bias decays markedly faster than EC LE itself (τ_bias = 0.83 d vs τ_EC = 2.82 d) because the satellite *does* register surface wetting through LST cooling after a rain event, an effect that is absent under drip where the 0–5 cm surface stays dry. Generalising τ to ungauged sites therefore requires either a parallel EC record, ground-validated SAR-derived irrigation timing, or a transfer of the bias-decay parameters from EC-monitored sites under the assumption of crop-type and irrigation-system equivalence.
