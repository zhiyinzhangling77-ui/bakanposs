"""Strict v14 SDS reproduction: summer-only filter + NDVI gate + global SWC quantiles.

SDS (Soil-moisture Drought Sensitivity) = 1 - mean(LE | dry) / mean(LE | normal)
Higher SDS means stronger drought sensitivity.

Key choices to match v14:
  * NDVI > 0.3 to restrict to vegetated growing periods (drops bare-soil days
    where ET is dominated by evaporation, not transpiration).
  * For TzM summer, SWC thresholds (p25, p75) are computed on the WHOLE TzM
    summer pool, then applied within each irrig_bucket. This makes the
    irrig-bucket comparisons apples-to-apples; v14's bucket SDS is reported
    against the same SWC reference, not bucket-internal quantiles.
  * For Oran by season, thresholds are computed per season (no irrigation).

Reads:  ec_daily_master_flagged.csv
Writes: sds_v14_results.csv  (table of SDS estimates with 95% CI)
"""
from pathlib import Path
import numpy as np
import pandas as pd

CSV_IN = Path("/home/shion-nagamine/bakanposs/ec_daily_master_flagged.csv")
CSV_OUT = Path("/home/shion-nagamine/bakanposs/sds_v14_results.csv")

VARS = {"LE": "LE_Wm2", "ET": "ET_mm"}
NDVI_GATE = 0.3
N_BOOT = 2000
RNG = np.random.default_rng(42)


def sds_from_arrays(values: np.ndarray, swc: np.ndarray,
                    p25: float, p75: float) -> float:
    dry = swc < p25
    norm = (swc >= p25) & (swc <= p75)
    if dry.sum() < 5 or norm.sum() < 5:
        return np.nan
    m_norm = np.nanmean(values[norm])
    m_dry = np.nanmean(values[dry])
    if not np.isfinite(m_norm) or m_norm <= 0:
        return np.nan
    return 1.0 - m_dry / m_norm


def bootstrap_ci(values: np.ndarray, swc: np.ndarray, p25: float, p75: float,
                 n: int = N_BOOT) -> tuple[float, float, float]:
    idx = np.arange(len(values))
    samples = []
    for _ in range(n):
        b = RNG.choice(idx, size=len(idx), replace=True)
        samples.append(sds_from_arrays(values[b], swc[b], p25, p75))
    samples = np.array(samples)
    samples = samples[np.isfinite(samples)]
    if len(samples) < 50:
        return np.nan, np.nan, np.nan
    return (float(np.median(samples)),
            float(np.quantile(samples, 0.025)),
            float(np.quantile(samples, 0.975)))


def add_ef(df: pd.DataFrame) -> pd.DataFrame:
    if "LE_Wm2" in df and "Rn_Wm2" in df:
        ef = df["LE_Wm2"] / df["Rn_Wm2"].replace(0, np.nan)
        df = df.assign(EF=ef.where((ef > -0.2) & (ef < 1.5)))
    return df


def gate(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["NDVI"].fillna(0) > NDVI_GATE].copy()


def compute_block(df: pd.DataFrame, label: str,
                  p25: float | None = None,
                  p75: float | None = None) -> list[dict]:
    rows = []
    if len(df) < 20:
        return rows
    swc = df["SWC"].to_numpy()
    if p25 is None:
        p25 = float(np.nanquantile(swc, 0.25))
    if p75 is None:
        p75 = float(np.nanquantile(swc, 0.75))

    df_ef = add_ef(df)
    var_cols = {"LE": "LE_Wm2", "ET": "ET_mm", "EF": "EF"}
    for var_name, col in var_cols.items():
        if col not in df_ef:
            continue
        vals = df_ef[col].to_numpy()
        ok = np.isfinite(vals) & np.isfinite(swc)
        if ok.sum() < 30:
            continue
        v, s = vals[ok], swc[ok]
        med, lo, hi = bootstrap_ci(v, s, p25, p75)
        rows.append({
            "stratum": label, "var": var_name,
            "n_total": int(ok.sum()),
            "n_dry": int((s < p25).sum()),
            "n_norm": int(((s >= p25) & (s <= p75)).sum()),
            "swc_p25": round(p25, 4),
            "swc_p75": round(p75, 4),
            "SDS": med, "CI95_lo": lo, "CI95_hi": hi,
        })
    return rows


def main() -> None:
    df = pd.read_csv(CSV_IN, parse_dates=["date"])
    df = gate(df)
    rows = []

    # Oran: per-season thresholds, no irrig stratification
    for season in ["winter", "spring", "summer", "fall"]:
        sub = df[(df.site == "Oran") & (df.season == season)]
        rows.extend(compute_block(sub, f"Oran_{season}"))

    # TzM summer: compute thresholds on the whole TzM summer pool,
    # then apply the same thresholds within each irrig_bucket.
    tzm_summer = df[(df.site == "TzM") & (df.season == "summer")]
    swc_summer = tzm_summer["SWC"].to_numpy()
    p25_s = float(np.nanquantile(swc_summer, 0.25))
    p75_s = float(np.nanquantile(swc_summer, 0.75))

    rows.extend(compute_block(tzm_summer, "TzM_summer_all", p25_s, p75_s))
    for bk in ["irrig_d0-3", "irrig_d4-7", "irrig_d8plus"]:
        sub = tzm_summer[tzm_summer.irrig_bucket == bk]
        rows.extend(compute_block(sub, f"TzM_summer_{bk}", p25_s, p75_s))

    # TzM shoulder seasons for context
    for season in ["spring", "fall"]:
        sub = df[(df.site == "TzM") & (df.season == season)]
        rows.extend(compute_block(sub, f"TzM_{season}"))

    out = pd.DataFrame(rows)
    out.to_csv(CSV_OUT, index=False)
    pd.set_option("display.float_format", "{:+.3f}".format)
    print(f"\nwrote {CSV_OUT}\n")
    print(f"NDVI gate: > {NDVI_GATE}")
    print(f"TzM summer SWC thresholds: p25={p25_s:.4f}  p75={p75_s:.4f}")
    print()
    print(out.to_string(index=False))

    print("\n--- v14 reference values ---")
    print("  Oran spring        SDS ~ +0.69 [+0.62, +0.72]")
    print("  TzM summer all     SDS ~ +0.27 [+0.14, +0.35]")
    print("  TzM summer d0-3    SDS ~  0.00")
    print("  TzM summer d4-7    SDS ~ +0.38")
    print("  TzM summer d8+     SDS ~ +0.28")


if __name__ == "__main__":
    main()
