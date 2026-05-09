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
2. **Quantify the MOD16 and PML ET bias** against EC observations at both a drip-irrigated orchard (TzM, 2020–2024) and a rainfed cereal control (Oran, 2018–2020), stratified by season, irrigation bucket, and vegetation state.
3. **Test the dose-response prediction**: if the bias is irrigation-driven, it should scale monotonically with days since last irrigation and be absent at the rainfed control site.
4. **Fit an exponential decay bias model** of the form Δ(t) = a·exp(−t/τ) + c, where t is days since last irrigation, and recover the bias decay time constant τ and any structural permanent offset c.

We find that satellite ET is underestimated by 2.7–4.1 mm d⁻¹ in the 3 days after irrigation at TzM, relaxing toward a permanent offset of −2.6 mm d⁻¹ for MOD16 and near zero for PML, with τ ≈ 3–5 d. The rainfed Oran control shows biases less than 0.5 mm d⁻¹ with no irrigation-cycle structure. These results are consistent with the drip wetted-bulb mechanism and point to distinct correction strategies for the two products.

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

### 2.6 Satellite ET bias and decay model
We define
$$\Delta_X(t) = ET_X^{\text{sat}}(t) - ET^{EC}(t)$$
for product X ∈ {MOD16, PML} on overlapping site-days. Stratified bias is reported by site, season, and irrigation bucket on the summer × NDVI > 0.3 subset.

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
Across both sites and all overlapping days (1,355 with MOD16, 1,214 with PML):
- **Oran**: MBE_MOD16 = −0.23 mm d⁻¹, MBE_PML = +0.44 mm d⁻¹.
- **TzM**: MBE_MOD16 = −2.69 mm d⁻¹, MBE_PML = −1.45 mm d⁻¹.

Both products track Oran reasonably (small bias, weak season structure) but systematically underestimate TzM, especially in summer. PML overestimates Oran by ≈ 0.4 mm d⁻¹ on average, a known PML behaviour for sparse vegetation.

### 3.4 Bias by irrigation bucket (Fig 5, headline figure)
Restricting TzM to summer × NDVI > 0.3:

| irrig bucket | n   | MBE MOD16 (mm/d) | MBE PML (mm/d) |
|---|---:|---:|---:|
| d0–3       | 281 | −4.12 | −2.78 |
| d4–7       |  75 | −2.68 | −1.29 |
| d8+        |  44 | −2.68 | −0.93 |

Both products show a clear dose response. MOD16 plateaus at ≈ −2.7 mm d⁻¹ from d4 onward, whereas PML continues to relax toward zero, suggesting different bias structures.

### 3.5 Linking SDS to satellite bias (Fig 6)
Across all strata where both metrics are available, low SDS is associated with strongly negative MBE: Oran spring (SDS +0.43 → MBE_MOD16 −0.27) versus TzM summer d0–3 (SDS +0.13 → MBE_MOD16 −4.12). The pattern is consistent for PML.

### 3.6 Exponential decay model (Fig 7)
On TzM summer × NDVI > 0.3 daily data:

**MOD16 full model**: Δ(t) = −2.04 · exp(−t/3.36 d) − 2.64
The transient amplitude (a = −2.0 mm d⁻¹) and decay time constant (τ = 3.4 d) align with the irrigation cycle (events every 2–3 d). The permanent offset (c = −2.6 mm d⁻¹) implies MOD16 underestimates irrigated almond ET *even at quasi-steady state*, consistent with calibration on FLUXNET sites that under-represent high-LAI orchard crops.

**PML transient model**: Δ(t) = −3.33 · exp(−t/5.34 d)
PML's permanent offset is small (c ≈ −0.85 mm d⁻¹ but poorly constrained); the transient model with a = −3.3 mm d⁻¹ and τ = 5.3 d already accounts for ≥ 95% of the bias at d ≥ 8. PML therefore contains a *correctable* irrigation timing error rather than a structural offset.

## 4. Discussion (key paragraphs draft)

### 4.1 Mechanism
Our results are consistent with the drip wetted-bulb mechanism: drip irrigation creates a localised, well-aerated, near-saturated zone at 10–30 cm depth that supports transpiration through almond's medium-deep roots without ever raising the 0–5 cm surface SWC that the in-situ probe and microwave SM retrievals see. After ≈ 3–4 days the wetted bulb depletes and LE again becomes coupled to the broader root-zone moisture, although the absolute SWC level remains low. The decay time constant we recover from EC vs satellite bias (τ ≈ 3.4–5.3 d) closely matches the irrigation interval (2–3 d), confirming that the bias is event-driven rather than seasonal.

### 4.2 Implications for satellite ET retrieval
The two products tested differ in how their bias is structured. **MOD16** carries both a transient irrigation overshoot and a *permanent* underestimation of ≈ 2.6 mm d⁻¹ for irrigated almond — likely because its leaf-area-driven canopy conductance scheme was calibrated against datasets with limited coverage of high-LAI orchard crops. **PML** lacks a substantial permanent offset but exhibits a longer transient (τ ≈ 5 d). For irrigation scheduling and water accounting in drip-irrigated horticulture, this means:
- MOD16 needs both a structural bias correction (calibrated per crop type) and an event-aware transient correction.
- PML can in principle be corrected with a single exponential decay term conditioned on days since last irrigation.

The latter is operationally feasible because irrigation timing can be retrieved from SAR backscatter (Sentinel-1) or from field records.

### 4.3 Reframing the deep-root narrative
The original investigation hypothesised that almond's deep roots access groundwater to maintain transpiration during dry-surface periods. The present analysis, complemented by the v14 days-since-irrigation stratification, shows instead that the apparent decoupling is fully explained by drip irrigation buffering surface moisture independence on a 3–4 day timescale. In the absence of irrigation (TzM fall, Oran summer), SDS recovers and satellite bias shrinks. The data therefore do not support a generic "deep root" hypothesis for almond at this site; the decoupling is management-driven.

## 5. Limitations
A single drip-irrigated and a single rainfed site limit generalisation to other irrigation systems (flood, sprinkler) and to other tree species. The 5 cm SWC sensor under-samples the 10–30 cm wetted bulb, so the mechanistic interpretation of section 4.1 rests on indirect evidence. ERA5 dewpoint and Rn were not retrieved in the present GEE pipeline, so VPD was taken only from EC. SMAP root-zone SM was unavailable at submission and would tighten section 4.2. The bootstrap CIs on the decay parameters are wide (especially for c in the PML full model), reflecting the small number of d ≥ 8 days available; longer time series and additional drip-irrigated sites would tighten these estimates.
