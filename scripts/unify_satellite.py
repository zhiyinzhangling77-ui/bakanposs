"""Unify GEE-downloaded satellite CSVs into a long-format daily table.

Most Fast_*.csv files are wide: columns are "YYYY-MM-DD_<variable>" plus
buf/name/.geo, with one row per site. Two files are already long:
  - OranTzM_S2_NDVI.csv  (name, date, product, B8/B4/B11/B3)
  - PML_V2_TimeSeries_Oran_TzM.csv  (name, date, GPP/ET/Ec/Es/Ei) — 2023 only

We melt the wide files, apply unit conversions, harmonise the long files,
concat into one big long table, then pivot to (date x site) x variable.

Reads:  the 8 GEE files in ~/Downloads/drive-download-.../
Writes: satellite_long.csv      raw long table for inspection
        satellite_daily.csv     wide-by-variable, one row per (date, site),
                                 8-day products forward-filled to daily
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

BASE = Path("/mnt/hdd/Dataset")
REPO = Path(__file__).parent.parent
OUT_LONG = REPO / "satellite_long.csv"
OUT_DAILY = REPO / "satellite_daily.csv"

DATE_VAR = re.compile(r"^(\d{4}-\d{2}-\d{2})_(.+)$")

WIDE_FILES = {
    "CHIRPS": "Fast_OranTzM_CHIRPS.csv",
    "ERA5":   "Fast_OranTzM_ERA5_Daily.csv",
    "LAI":    "Fast_OranTzM_LAI.csv",
    "LST":    "Fast_OranTzM_LST.csv",
    "MOD16":  "Fast_OranTzM_MOD16.csv",
    "PML":    "Fast_OranTzM_PML.csv",
}


def melt_wide(p: Path, source: str) -> pd.DataFrame:
    df = pd.read_csv(p, low_memory=False)
    if "name" not in df.columns:
        raise ValueError(f"{p.name}: missing 'name' column")

    date_cols = [c for c in df.columns if DATE_VAR.match(c)]
    if not date_cols:
        return pd.DataFrame()

    long = df.melt(id_vars=["name"], value_vars=date_cols,
                   var_name="key", value_name="value")
    parts = long["key"].str.extract(DATE_VAR.pattern)
    long["date"] = pd.to_datetime(parts[0])
    long["variable"] = parts[1]
    long["source"] = source
    long["value"] = pd.to_numeric(long["value"], errors="coerce")
    return long[["date", "name", "source", "variable", "value"]]


def load_s2() -> pd.DataFrame:
    p = BASE / "OranTzM_S2_NDVI.csv"
    df = pd.read_csv(p)
    df["date"] = pd.to_datetime(df["date"])
    for c in ["B8", "B4", "B11", "B3"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["NDVI"] = (df["B8"] - df["B4"]) / (df["B8"] + df["B4"])
    df["NDWI"] = (df["B3"] - df["B8"]) / (df["B3"] + df["B8"])
    df["NDMI"] = (df["B8"] - df["B11"]) / (df["B8"] + df["B11"])
    rows = []
    for var in ["NDVI", "NDWI", "NDMI"]:
        sub = df[["date", "name", var]].rename(columns={var: "value"})
        sub["source"] = "S2"
        sub["variable"] = var
        rows.append(sub[["date", "name", "source", "variable", "value"]])
    return pd.concat(rows)


def load_pml_long() -> pd.DataFrame:
    p = BASE / "PML_V2_TimeSeries_Oran_TzM.csv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p)
    df["date"] = pd.to_datetime(df["date"])
    rows = []
    for var in ["GPP", "ET", "Ec", "Es", "Ei"]:
        sub = df[["date", "name", var]].rename(columns={var: "value"})
        sub["source"] = "PML_V2_long"
        sub["variable"] = var
        rows.append(sub[["date", "name", "source", "variable", "value"]])
    return pd.concat(rows)


def apply_units(long: pd.DataFrame) -> pd.DataFrame:
    """Convert to physical units: ET in mm/day, T in C, precip in mm, Rn in W/m²."""
    out = long.copy()

    # MOD16: stored as kg/m^2/8day with scale 0.1 → mm/day
    mod_mask = (out.source == "MOD16") & out.variable.isin(["ET", "PET"])
    out.loc[mod_mask, "value"] = out.loc[mod_mask, "value"] * 0.1 / 8.0

    # LST: scale 0.02 to Kelvin, convert to C
    lst_mask = (out.source == "LST") & out.variable.str.contains("LST")
    out.loc[lst_mask, "value"] = out.loc[lst_mask, "value"] * 0.02 - 273.15

    # LAI: scale 0.1 (MCD15A3H raw 0-100 → 0-10)
    lai_mask = (out.source == "LAI") & out.variable.isin(["Lai"])
    out.loc[lai_mask, "value"] = out.loc[lai_mask, "value"] * 0.1
    fpar_mask = (out.source == "LAI") & out.variable.isin(["Fpar"])
    out.loc[fpar_mask, "value"] = out.loc[fpar_mask, "value"] * 0.01

    # ERA5: temperatures K → C, precip m → mm, Rn J/m² (daily) → W/m²
    era_temp = (out.source == "ERA5") & out.variable.str.contains("temperature_2m")
    out.loc[era_temp, "value"] = out.loc[era_temp, "value"] - 273.15
    era_p = (out.source == "ERA5") & (out.variable == "total_precipitation")
    out.loc[era_p, "value"] = out.loc[era_p, "value"] * 1000.0
    era_rn = (out.source == "ERA5") & (out.variable == "surface_net_solar_radiation")
    out.loc[era_rn, "value"] = out.loc[era_rn, "value"] / 86400.0

    return out


def derive_extras(long: pd.DataFrame) -> pd.DataFrame:
    """Compute MOD16 ET as the only ET, PML total ET = Ec+Es+Ei, ERA5 VPD."""
    extras = []

    # PML wide-file total ET
    pml = long[long.source == "PML"].pivot_table(
        index=["date", "name"], columns="variable", values="value"
    ).reset_index()
    if {"Ec", "Es", "Ei"}.issubset(pml.columns):
        pml["value"] = pml[["Ec", "Es", "Ei"]].sum(axis=1, min_count=1)
        pml["source"] = "PML"
        pml["variable"] = "ET"
        extras.append(pml[["date", "name", "source", "variable", "value"]])

    # ERA5 VPD from Ta and Tdp
    era = long[long.source == "ERA5"].pivot_table(
        index=["date", "name"], columns="variable", values="value"
    ).reset_index()
    if {"temperature_2m", "dewpoint_temperature_2m"}.issubset(era.columns):
        ta = era["temperature_2m"]
        td = era["dewpoint_temperature_2m"]
        es = 0.6108 * np.exp(17.27 * ta / (ta + 237.3))   # kPa
        ea = 0.6108 * np.exp(17.27 * td / (td + 237.3))
        era["value"] = (es - ea).clip(lower=0)
        era["source"] = "ERA5"
        era["variable"] = "VPD_kPa"
        extras.append(era[["date", "name", "source", "variable", "value"]])

    if extras:
        return pd.concat([long] + extras, ignore_index=True)
    return long


def to_daily_wide(long: pd.DataFrame) -> pd.DataFrame:
    """Pivot long → wide by source/variable; forward-fill 8-day series to daily."""
    long["col"] = long["source"] + "_" + long["variable"]
    wide = long.pivot_table(index=["date", "name"], columns="col",
                            values="value", aggfunc="mean").reset_index()

    # build a full daily index per site
    full = []
    for site in wide["name"].unique():
        s = wide[wide["name"] == site].set_index("date").sort_index()
        idx = pd.date_range(s.index.min(), s.index.max(), freq="D")
        s = s.reindex(idx).rename_axis("date")
        s["name"] = site
        # forward-fill 8-day MOD16/PML/LAI/LST etc up to 8 days
        ffill_cols = [c for c in s.columns
                      if c.startswith(("MOD16_", "PML_", "LAI_", "LST_", "S2_"))]
        s[ffill_cols] = s[ffill_cols].ffill(limit=8)
        full.append(s.reset_index())
    out = pd.concat(full, ignore_index=True)
    out = out.rename(columns={"name": "site"})
    return out


def main() -> None:
    longs = []
    for source, fname in WIDE_FILES.items():
        p = BASE / fname
        if not p.exists():
            print(f"[skip] {fname} missing")
            continue
        m = melt_wide(p, source)
        if len(m):
            longs.append(m)
            print(f"  {source:7s}  {len(m):>6d} rows  "
                  f"vars={sorted(m.variable.unique())}")

    longs.append(load_s2())
    print(f"  {'S2':7s}  {len(longs[-1]):>6d} rows  "
          f"vars={sorted(longs[-1].variable.unique())}")

    pml_long = load_pml_long()
    if len(pml_long):
        longs.append(pml_long)
        print(f"  {'PML_V2':7s}  {len(pml_long):>6d} rows (2023 only, "
              f"flagged separately)")

    long = pd.concat(longs, ignore_index=True)
    long = long.dropna(subset=["value"])
    long = apply_units(long)
    long = derive_extras(long)

    long.to_csv(OUT_LONG, index=False)
    print(f"\nwrote {OUT_LONG}  ({len(long):,} rows)")

    daily = to_daily_wide(long)
    daily.to_csv(OUT_DAILY, index=False)
    print(f"wrote {OUT_DAILY}  ({len(daily):,} daily rows × "
          f"{len(daily.columns) - 2} vars)")

    et_cols = [c for c in daily.columns if "_ET" in c]
    if et_cols:
        print("\nET column summary (mm/day):")
        print(daily.groupby("site")[et_cols].agg(["count", "mean"]).round(2))


if __name__ == "__main__":
    main()
