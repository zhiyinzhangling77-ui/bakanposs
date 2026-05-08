"""Unify the two daily EC files (TzM almond + Oran cereal) into one master CSV.

Inputs (paths assume original layout on shion-nagamine machine):
  - Daily_Summary_Filtered_forPred_ActEne26.csv      (TzM, almond)
  - Oran_EddyDaily_MASTER_2018_2020_correct.csv      (Oran, cereal)

Output:
  - ec_daily_master.csv                              (both sites concatenated)

Unified schema:
  date, site, crop, year, doy,
  ET_mm, LE_Wm2, H_Wm2, G_Wm2, Rn_Wm2,
  Ta_mean, Ta_max, Ta_min,
  VPD_kPa_mean, VPD_kPa_max,
  SWC, P_mm, Irrig_mm, NDVI, GPP_gC_m2_d, n_obs_le
"""
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("/home/shion-nagamine/Dataset/Eddy data in Spain")
OUT = Path("/home/shion-nagamine/bakanposs/ec_daily_master.csv")

UNIFIED_COLS = [
    "date", "site", "crop", "year", "doy",
    "ET_mm", "LE_Wm2", "H_Wm2", "G_Wm2", "Rn_Wm2",
    "Ta_mean", "Ta_max", "Ta_min",
    "VPD_kPa_mean", "VPD_kPa_max",
    "SWC", "P_mm", "Irrig_mm",
    "NDVI", "GPP_gC_m2_d", "n_obs_le",
]


def load_tzm() -> pd.DataFrame:
    """TzM almond daily (file 1).

    Date format: YYYY/MM/DD. ET_sum is mm/day. LE/H/G/Rn _avg are W/m^2 means.
    """
    p = BASE / "Daily_Summary_Filtered_forPred_ActEne26.csv"
    df = pd.read_csv(p)
    df["date"] = pd.to_datetime(df["date"], format="%Y/%m/%d", errors="coerce")
    out = pd.DataFrame({
        "date": df["date"],
        "site": "TzM",
        "crop": "almond",
        "year": df["Year"],
        "doy": df["DOY"],
        "ET_mm": df["ET_sum"],
        "LE_Wm2": df["LE_avg"],
        "H_Wm2": df["H_avg"],
        "G_Wm2": df["G_avg"],
        "Rn_Wm2": df["NetRad_avg"],
        "Ta_mean": df["air_temperatureC_mean"],
        "Ta_max": df["air_temperatureC_max"],
        "Ta_min": df["air_temperatureC_min"],
        # File1's VPD_mean appears to be in Pa (~1858), VPD_kPa is kPa.
        "VPD_kPa_mean": df["VPD_mean"] / 1000.0,
        "VPD_kPa_max": df["VPD_max"] / 1000.0,
        "SWC": df["SWC_avg"],
        "P_mm": df["Rain_mm"].astype(float),
        "Irrig_mm": df["Irrig_mm"],
        "NDVI": df["NDVI_interp"],
        "GPP_gC_m2_d": df["GPP_sum"],
        "n_obs_le": np.nan,
    })
    return out


def load_oran() -> pd.DataFrame:
    """Oran cereal daily (file 4).

    Date format: DD/MM/YYYY. ET is Eddy_ET_mm_d. Daily LE/H/G not present
    (only NETRAD aggregates); leave as NaN, recompute later from 30-min if
    needed.
    """
    p = BASE / "Oran_EddyDaily_MASTER_2018_2020_correct.csv"
    df = pd.read_csv(p)
    df["date"] = pd.to_datetime(df["Date"], format="%d/%m/%Y", errors="coerce")
    out = pd.DataFrame({
        "date": df["date"],
        "site": "Oran",
        "crop": df["crop"].fillna("cereal"),
        "year": df["Year"],
        "doy": df["DOY"],
        "ET_mm": df["Eddy_ET_mm_d"],
        "LE_Wm2": np.nan,
        "H_Wm2": np.nan,
        "G_Wm2": np.nan,
        "Rn_Wm2": df["NETRAD_mean_Wm2"],
        "Ta_mean": df["Eddy_Tair_mean"],
        "Ta_max": df["Eddy_Tair_max"],
        "Ta_min": df["Eddy_Tair_min"],
        "VPD_kPa_mean": df["Eddy_VPD_mean_kPa"],
        "VPD_kPa_max": df["Eddy_VPD_max_kPa"],
        "SWC": df["SWC_avg_mean"],
        "P_mm": df["Rain_mm"].astype(float),
        "Irrig_mm": 0.0,
        "NDVI": df["NDVI_int"],
        "GPP_gC_m2_d": df["Eddy_GPP_gC_m2_d"],
        "n_obs_le": df["n_obs_le"],
    })
    return out


def main() -> None:
    parts = [load_tzm(), load_oran()]
    df = pd.concat(parts, ignore_index=True)
    df = df.dropna(subset=["date"]).sort_values(["site", "date"]).reset_index(drop=True)
    df = df[UNIFIED_COLS]

    df.to_csv(OUT, index=False)

    print(f"wrote {OUT}  rows={len(df)}")
    print("\nrange per site:")
    print(df.groupby("site")["date"].agg(["min", "max", "count"]))
    print("\nmissing per column:")
    print(df.isna().sum().to_string())
    print("\nhead:")
    print(df.head().to_string())


if __name__ == "__main__":
    main()
