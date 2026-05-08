"""Strict v14 SDS reproduction: summer-only filter + bootstrap CI.

SDS (Soil-moisture Drought Sensitivity) = 1 - mean(LE | SWC dry) / mean(LE | SWC normal)
where dry = SWC < p25 within site-season, normal = p25 <= SWC <= p75.
Higher SDS means stronger drought sensitivity.

Reproduces the v14 numbers:
  Oran spring SDS ~ +0.69
  TzM  summer SDS by irrig_bucket: d0-3 ~ 0,  d4-7 ~ +0.38,  d8+ ~ +0.28

Reads:  ec_daily_master_flagged.csv
Writes: sds_v14_results.csv  (table of SDS estimates with 95% CI)
Prints: results table to stdout
"""
from pathlib import Path
import numpy as np
import pandas as pd

CSV_IN = Path("/home/shion-nagamine/bakanposs/ec_daily_master_flagged.csv")
CSV_OUT = Path("/home/shion-nagamine/bakanposs/sds_v14_results.csv")

VARS = {"LE": "LE_Wm2", "ET": "ET_mm"}  # EF computed inline if Rn available
N_BOOT = 2000
RNG = np.random.default_rng(42)


def sds(values: np.ndarray, mask_dry: np.ndarray, mask_norm: np.ndarray) -> float:
    if mask_norm.sum() < 5 or mask_dry.sum() < 5:
        return np.nan
    m_norm = np.nanmean(values[mask_norm])
    m_dry = np.nanmean(values[mask_dry])
    if m_norm <= 0 or not np.isfinite(m_norm):
        return np.nan
    return 1.0 - m_dry / m_norm


def bootstrap_ci(values: np.ndarray, swc: np.ndarray, p25: float, p75: float,
                 n: int = N_BOOT) -> tuple[float, float, float]:
    idx_all = np.arange(len(values))
    samples = []
    for _ in range(n):
        boot = RNG.choice(idx_all, size=len(idx_all), replace=True)
        v = values[boot]
        s = swc[boot]
        dry = s < p25
        norm = (s >= p25) & (s <= p75)
        samples.append(sds(v, dry, norm))
    samples = np.array(samples)
    samples = samples[np.isfinite(samples)]
    if len(samples) < 50:
        return np.nan, np.nan, np.nan
    return (np.median(samples),
            np.quantile(samples, 0.025),
            np.quantile(samples, 0.975))


def compute_for(df: pd.DataFrame, label: str) -> list[dict]:
    """Compute SDS for every variable on a sub-DataFrame."""
    rows = []
    if len(df) < 20:
        return rows
    swc = df["SWC"].to_numpy()
    p25 = np.nanquantile(swc, 0.25)
    p75 = np.nanquantile(swc, 0.75)
    for var_name, col in VARS.items():
        if col not in df:
            continue
        vals = df[col].to_numpy()
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
            "SDS": med, "CI95_lo": lo, "CI95_hi": hi,
        })
    # EF separately if Rn available
    if "Rn_Wm2" in df and df["Rn_Wm2"].notna().any() and df["LE_Wm2"].notna().any():
        ef = df["LE_Wm2"] / df["Rn_Wm2"].replace(0, np.nan)
        ef = ef.where((ef > -0.2) & (ef < 1.5))
        df_ef = df.assign(EF=ef)
        ok = np.isfinite(df_ef["EF"]) & np.isfinite(swc)
        if ok.sum() >= 30:
            v = df_ef["EF"].to_numpy()[ok]
            s = swc[ok]
            med, lo, hi = bootstrap_ci(v, s, p25, p75)
            rows.append({
                "stratum": label, "var": "EF",
                "n_total": int(ok.sum()),
                "n_dry": int((s < p25).sum()),
                "n_norm": int(((s >= p25) & (s <= p75)).sum()),
                "SDS": med, "CI95_lo": lo, "CI95_hi": hi,
            })
    return rows


def main() -> None:
    df = pd.read_csv(CSV_IN, parse_dates=["date"])
    rows = []

    # Oran: by season, no irrigation buckets
    for season in ["winter", "spring", "summer", "fall"]:
        sub = df[(df.site == "Oran") & (df.season == season)]
        rows.extend(compute_for(sub, f"Oran_{season}"))

    # TzM: summer-only, then split by irrig_bucket
    tzm_summer = df[(df.site == "TzM") & (df.season == "summer")]
    rows.extend(compute_for(tzm_summer, "TzM_summer_all"))
    for bk in ["irrig_d0-3", "irrig_d4-7", "irrig_d8plus"]:
        sub = tzm_summer[tzm_summer.irrig_bucket == bk]
        rows.extend(compute_for(sub, f"TzM_summer_{bk}"))

    # TzM: shoulder seasons (spring/fall) for context
    for season in ["spring", "fall"]:
        sub = df[(df.site == "TzM") & (df.season == season)]
        rows.extend(compute_for(sub, f"TzM_{season}"))

    out = pd.DataFrame(rows)
    out.to_csv(CSV_OUT, index=False)

    # pretty print
    print(f"\nwrote {CSV_OUT}\n")
    pd.set_option("display.float_format", "{:+.3f}".format)
    print(out.to_string(index=False))

    print("\n--- v14 reference values ---")
    print("  Oran spring        SDS ~ +0.69 [+0.62, +0.72]")
    print("  TzM summer all     SDS ~ +0.27 [+0.14, +0.35]")
    print("  TzM summer d0-3    SDS ~  0.00")
    print("  TzM summer d4-7    SDS ~ +0.38")
    print("  TzM summer d8+     SDS ~ +0.28")


if __name__ == "__main__":
    main()
