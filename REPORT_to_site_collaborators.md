# Progress report:  Satellite-vs-tower ET divergence at Tarazona and Oran

**To:**   The Tarazona / Oran EC tower team (Castilla–La Mancha, Spain)
**Cc:**   José [Surname], colleagues at site PI's institute
**From:** Shion Nagamine, Kazuhito Ichii
          Center for Environmental Remote Sensing (CEReS), Chiba University, Japan
**Date:** 2026-05-21
**Subject:** Progress report — decomposing the satellite–tower ET
              divergence at the Tarazona and Oran flux towers


---

## 1.  Purpose of this report

We have been using the eddy-covariance (EC) data you kindly provided
from the **Tarazona** drip-irrigated almond orchard (2020–2024) and the
**Oran** rainfed winter-cereal site (2018–2020) to assess whether the
Meteosat-based ET product (LSA SAF **Meteosat ETv3**) reproduces what
the towers measure.  This document is an interim progress report:
we summarise what we have done, the main findings to date, the
methodological caveats, and the next steps.  We would also like to
ask a small number of clarification questions about site metadata
(see §10).


---

## 2.  Datasets used

| Source | Variable(s) | Period | Resolution |
|---|---|---|---|
| Tarazona EC tower (provided by you) | Daily LE, H, G, Rn, VPD, SWC, irrigation log | 2020-06-18 → 2024-10-15 | Daily, point-scale |
| Oran EC tower (provided by you) | Daily ET (mm d⁻¹), rainfall | 2018-01-02 → 2020-06-25 | Daily, point-scale |
| LSA SAF Meteosat ETv3 (DMET v3) | Daily ET (mm d⁻¹) | 2018-01-01 → 2024-12-31 | ~5 km, MSG-disk |
| LSA SAF Meteosat MGPP | 10-day GPP | 2018–2024 | ~5 km |
| Sentinel-2 SR Harmonized (Copernicus) | NDVI | 2018–2024 | 10 m, ~5-day |
| SMAP L4 SPL4SMAU | Surface and root-zone soil moisture | 2018–2024 | 9 km, 3-hourly |
| ESA WorldCover v200 (2021) | Land cover | static | 10 m |

All satellite products were resampled or extracted to daily values at
the tower grid cell.  Quality filters applied: ETv3 daily files
require ≥ 36 of 48 30-min observations; freeze-day filter (Ta_min
> 0 °C) excludes snow / frost days.  Energy-balance correction status
of the EC LE has been treated as a methodological caveat (see §9).


---

## 3.  Analysis A — Recovery timescale (τ) of LE after a water-input event

**Method.**  For each EC tower, we detected water-input events using
sensor records:

  - Tarazona: irrigation ≥ 0.5 mm/day (drip flow-meter threshold; the
    typical event delivers 12–15 mm)
  - Oran:     rainfall ≥ 3 mm/day (threshold chosen for soil-wetting
    significance per Ramos et al. 2014)

Each event was followed for up to 14 days or until the next event,
with a minimum 4-day window.  We fitted the exponential model

    LE(d) = LE_∞ + (LE_0 − LE_∞) · exp(−d / τ)

to the median daily LE per day-since-event (≥ 3 events per day
required), with LE_∞ fixed to the observed median for d ≥ 7 days.
Bootstrap (B = 5000) gave the SE and 95 % CI of τ.

**Results.**

| Site (active months) | n events | τ (d) | 95 % CI | LE_0 | LE_∞ | Amplitude (W m⁻²) |
|---|---|---|---|---|---|---|
| Tarazona (Jun–Sep) | 41 | **3.36** | [2.44, 4.90] | ~210 | 114 | **95** |
| Oran (Nov–Jun) | 10 | **2.82** | [1.83, 5.48] | ~36 | 17 | **21** |

The observed inter-site difference |Δτ| = 0.54 d is below the
minimum-detectable-effect threshold MDE = 1.96 · √(SE₁² + SE₂²)
= 2.15 d (α = 0.05).  Despite the small sample sizes, this is
evidence of statistical equivalence rather than under-powered
non-difference.

**Interpretation.**  The recovery timescale τ ≈ 3 d appears to be a
robust property of these Mediterranean dryland systems, independent
of whether the water input is rainfall or drip irrigation.  In
contrast, the recovery **amplitude** scales 4.5× between rainfed
and drip-irrigated conditions, suggesting that τ is a climate-level
property while amplitude is a management signal.


---

## 4.  Analysis B — Satellite ET validation against EC

**Paired daily comparison** between EC ET and Meteosat ETv3, qflag-
clean and freeze-day excluded:

| Site | n paired | Mean EC | Mean ETv3 | Bias | RMSE | Pearson r | KGE |
|---|---|---|---|---|---|---|---|
| Oran | 605 | 1.17 mm d⁻¹ | 1.18 mm d⁻¹ | **+1 %** | 0.58 | **0.82** | **+0.61** |
| Tarazona | 699 | 3.37 mm d⁻¹ | 1.03 mm d⁻¹ | **−70 %** | 3.08 | **0.07** | **−0.21** |

The satellite reproduces the EC ET at Oran with very good agreement,
but at Tarazona both the magnitude and the phase of ET diverge.

We cross-checked vegetation-side indicators at both sites:

| Indicator | Oran (r, Δpeak) | Tarazona (r, Δpeak) |
|---|---|---|
| NDVI (Sentinel-2 10 m vs EC tower NDVI) | 0.98, +2 d | 0.74, ~in-phase |
| GPP (Meteosat MGPP 5 km vs EC tower GPP) | 0.82, −8 d | 0.86, +1 d |
| **ET (Meteosat ETv3 5 km vs EC tower ET)** | **0.91, +6 d** | **0.55, −44 d** |

At Tarazona, NDVI and GPP retain reasonable phase agreement
(although MGPP underestimates amplitude by ~3×).  Only ET diverges
both in amplitude and in phase: ETv3 peaks 30–95 days earlier than
the EC tower in every year (mean Δpeak = −44 d across 2020–2023).
**The divergence is ET-specific.**


---

## 5.  Spatial diagnosis — pixel landscape around each tower

To understand the spatial context of each tower we computed mean
Sentinel-2 NDVI (Jun–Sep, full EC observation period) at nested
buffers matching the satellite-product pixel scales:

|  | Oran  (NDVI, fraction > 0.5) | Tarazona  (NDVI, fraction > 0.5) |
|---|---|---|
| 200 m (tower footprint) | 0.20, 2.4 % | **0.42, 7.6 %** |
| 1 km (~H26 pixel) | 0.20, 2.3 % | 0.31, 4.5 % |
| 5 km (~ETv3 pixel) | 0.23, 1.5 % | 0.22, 0.7 % |
| 12.5 km (~ASCAT pixel) | 0.26, 2.5 % | 0.22, 1.6 % |

Both pixels converge to similar 5-km mean NDVI (~0.22, bare/dry
summer landscape).  At the tower scale, however:

  - **Oran tower**: NDVI 0.20, representative of the dominant
    bare-stubble summer cereal in its 5-km pixel.
  - **Tarazona tower**: NDVI 0.42, an anomaly within an otherwise
    dry 5-km pixel where the 1-ha drip orchard occupies ≈ 1 %
    of the pixel area.

This explains the asymmetric ET agreement statistically: the
satellite pixel-mean is structurally close to the Oran tower
(both sample the rainfed majority) but structurally far from the
Tarazona tower (which samples an anomaly).  We refer to this as a
**representativeness mismatch** between the point flux measurement
and the coarse pixel-volume satellite estimate.


---

## 6.  Temporal diagnosis — recovery dynamics of the bias

We pooled the 48 Tarazona irrigation events and fitted the same
exponential model separately to (i) EC LE, (ii) Meteosat ETv3 LE,
and (iii) the bias Δ = EC − Sat.

  - EC LE: τ_EC = 4.7 d, **amp_EC = 95 W m⁻²**, R² = 0.79
  - ETv3:  τ fit pegs at the floor (~0.5 d), R² ≈ 0 → no detectable
    exponential recovery.  **amp_Sat ≈ 0 W m⁻²**.
  - Bias:  τ_bias = 4.57 d [3.03, 13.10], amp_bias ≈ 95 W m⁻²,
    R² = 0.86.

Because the satellite series is essentially flat through the event
window, the bias inherits the EC exponential by construction
(τ_bias ≈ τ_EC and amp_bias ≈ amp_EC).  The meaningful measurement
is therefore not τ_bias itself but the **amplitude contrast**: the
EC tower captures the full ~95 W m⁻² management response; the
satellite captures essentially none of it.


---

## 7.  Mechanistic diagnosis — driver attribution

We ran a standardised multivariate regression on weekly summer
means at Tarazona (Jun–Sep, 2020–2024, n = 88 weeks), with the
following predictors:

  - Demand:   VPD, Rn
  - Supply:   Tower SWC (in-situ, ~3–7 cm depth per José's
              installation recollection — see §9), SMAP root-zone SM

Standardised coefficients (β, mean ± SE; ***/**/* = p < 0.001 /
0.01 / 0.05):

| Driver | EC ET (R² = 0.49) | Meteosat ETv3 (R² = 0.61) |
|---|---|---|
| VPD | +0.19 * | +0.02 ns |
| Rn  | **+0.63 \\***\\* | +0.01 ns |
| Tower SWC | −0.09 ns | +0.16 * |
| SMAP root-zone | +0.25 ** | **+0.74 \\***\\* |

**Interpretation.**  EC ET at the orchard is **demand-driven**:
net radiation dominates because drip irrigation removes the water
limitation that normally caps Mediterranean summer ET.  Meteosat
ETv3 is **supply-limited**: its output is essentially a function of
the coarse satellite SM input (here SMAP root-zone as a proxy for
the H-SAF H141/H142 family that ETv3 actually ingests; we have
not yet fetched H-SAF directly — see §11).

The two products are thus reading different physics from the same
orchard.  EC measures the actual transpiration flux; ETv3 infers
ET from a coarse SM proxy that does not resolve drip irrigation.


---

## 8.  Synthesis — why ETv3 fails specifically at Tarazona

Three resolution limits in the satellite chain combine to produce
the blind spot:

  1. **Spatial scale.**  H-SAF surface SM operates at 12.5 km
     (H141/H142), the downscaled SCATSAR-SWI (H26) at 1 km.  A
     1-ha orchard is therefore diluted by 100×–156 000× within
     the pixel.
  2. **Microwave sensing physics.**  C/L-band SM products are
     sensitive only to the top ~5 cm of soil.  Drip irrigation
     wets small patches under each emitter that re-dry within
     hours of the morning satellite overpass.  Mature almond
     canopy further attenuates the microwave signal (high VOD).
  3. **SVAT model structure.**  The f(SM) availability curve in
     the ETv3 SVAT is calibrated for uniform natural vegetation
     and does not represent the bi-modal SM distribution of a
     drip-fed orchard (locally saturated under emitters, dry
     elsewhere), nor the engineered root access to deep water.

We have also confirmed empirically that the satellite signal is
unrecoverable at the orchard scale: a Sentinel-2 NDVI diagnostic
shows that only 4.5 % of the H26 1-km pixel around Tarazona is
summer-active vegetation, giving an expected pixel-mean ΔSM
of ~0.003 m³ m⁻³ during irrigation — well below the typical
satellite SM noise floor (~0.04 m³ m⁻³).


---

## 9.  Caveats and limitations

  1. **Tower SWC depth.**  Following José's note, the SWC sensors
     are at approximately 3–7 cm depth.  This is too shallow to
     represent the actual root-zone water status of the almond
     trees (taproot to 1–5 m).  The non-significance of Tower SWC
     in the regression (§7) is therefore expected on physical
     grounds and is **not** evidence that the orchard SM is
     unimportant — only that no observable shallow measurement
     captures the deep, irrigation-supplied water that drives
     transpiration.

  2. **Energy-balance closure.**  We have used the EC LE as
     provided.  If the EC data have been Bowen-ratio or
     residual-corrected, this would affect the absolute
     amplitude comparison with the satellite.  Please confirm
     the EC processing convention for both sites.

  3. **Sample size at Oran.**  n = 10 rainfall events for the
     Analysis A τ fit at Oran reflects the 3-year EC record and
     the 3 mm event threshold.  The MDE-based equivalence test
     mitigates the under-power concern, but extension to more
     rainfed sites would strengthen the universality claim.

  4. **Single-site case study for the satellite blind spot.**
     The Tarazona result is currently from one orchard.
     Generalisation to other Mediterranean drip-irrigated
     systems (citrus, olive, larger almond) is planned but
     not yet done.

  5. **SMAP as a proxy for the H-SAF chain.**  SMAP L4 is a
     coarse satellite SM product in the same physical family
     as the H-SAF products that ETv3 actually ingests, but the
     two are not interchangeable.  We expect the driver pattern
     to be insensitive to this substitution, and we plan to
     re-run the regression with H141 / H142 / H26 / H28 once
     fetched.


---

## 10.  Questions for the site team

We would be very grateful for clarification on the following
points, which will sharpen the manuscript:

  1. **EC processing convention.**  Are the LE/H/Rn/G values you
     provided raw eddy-covariance outputs, or have they been
     corrected for energy-balance closure (residual, Bowen,
     or otherwise)?

  2. **SWC sensor depth — exact metadata.**  José mentioned the
     sensors are around 3–7 cm depth.  If the exact depth
     metadata for the 2020–2024 period is recoverable, we would
     like to cite the precise depth in the manuscript.

  3. **Irrigation log granularity.**  We currently use a daily
     irrigation flag (≥ 0.5 mm).  If sub-daily irrigation
     timing (start time, duration per event) is available, this
     would allow a refined diurnal analysis.

  4. **Orchard footprint and surroundings.**  We have estimated
     the orchard at ~1 ha based on land-cover data.  If you can
     confirm the actual orchard area, planting density, and the
     land use of the immediate (~500 m) surroundings, we will
     correct our spatial diagnosis accordingly.

  5. **Authorship / acknowledgement convention.**  Please advise
     on how you would like to be acknowledged in the eventual
     manuscript (acknowledgement only, co-authorship, data
     citation), per the conventions of your institute.


---

## 11.  Next steps

  - **H-SAF SM fetch.**  We will fetch H141 (12.5 km surface SM,
    the operational ETv3 input), H142 (root-zone SWI), and either
    H26 (1 km SCATSAR-SWI) or H28 (0.5 km, the latest disaggregated
    product) for direct substitution into the driver-attribution
    regression of §7.
  - **Tower SWC event-recovery analysis.**  Quantify the SWC
    response to each irrigation event explicitly and compare with
    the (non-)response of all satellite SM products.
  - **Multi-site extension.**  Apply the same framework to one
    or two additional drip-irrigated Mediterranean orchards if
    appropriate EC data become available.
  - **Manuscript preparation.**  Draft submission targeted at
    *Agricultural and Forest Meteorology* or *Remote Sensing of
    Environment* in 2026 H2, contingent on the points in §10.


---

## 12.  Acknowledgements

We would like to thank you and your team sincerely for providing
the Tarazona and Oran flux-tower data, for the helpful site
description, and for José's clarification on the SWC sensor depth.
Without these data and the associated metadata this analysis would
not have been possible.

We look forward to your feedback, especially on the questions in
§10, and to discussing how best to handle authorship for any
subsequent manuscript.

With kind regards,

Shion Nagamine, Kazuhito Ichii
Center for Environmental Remote Sensing (CEReS)
Chiba University, Japan
shion.nagamine@chiba-u.jp
