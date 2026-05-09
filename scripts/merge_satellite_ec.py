"""Merge EC daily master with the unified satellite daily table.

Inputs:
  ec_daily_master_flagged.csv      EC + flags (irrig_bucket, season, etc.)
  satellite_daily.csv              date×site wide table from unify_satellite.py

Output:
  master_full.csv                  one row per (date, site) with EC + satellite
"""
from pathlib import Path
import pandas as pd

EC = Path("/home/shion-nagamine/bakanposs/ec_daily_master_flagged.csv")
SAT = Path("/home/shion-nagamine/bakanposs/satellite_daily.csv")
OUT = Path("/home/shion-nagamine/bakanposs/master_full.csv")


def main() -> None:
    ec = pd.read_csv(EC, parse_dates=["date"])
    sat = pd.read_csv(SAT, parse_dates=["date"])

    print(f"EC  rows: {len(ec):,}  dates {ec.date.min().date()}-{ec.date.max().date()}")
    print(f"SAT rows: {len(sat):,}  dates {sat.date.min().date()}-{sat.date.max().date()}")

    df = ec.merge(sat, on=["date", "site"], how="left")
    df.to_csv(OUT, index=False)

    print(f"\nwrote {OUT}  rows={len(df)}")

    et_cols = [c for c in ["ET_mm", "MOD16_ET", "PML_ET"] if c in df.columns]
    print("\nET coverage (non-null) per site:")
    print(df.groupby("site")[et_cols].count())
    print("\nET mean per site:")
    print(df.groupby("site")[et_cols].mean().round(3))

    # Bias preview
    if {"MOD16_ET", "ET_mm"}.issubset(df.columns):
        df["bias_MOD16"] = df["MOD16_ET"] - df["ET_mm"]
    if {"PML_ET", "ET_mm"}.issubset(df.columns):
        df["bias_PML"] = df["PML_ET"] - df["ET_mm"]
    bias_cols = [c for c in df.columns if c.startswith("bias_")]
    if bias_cols:
        print("\nsat - EC bias (mm/day) by site:")
        print(df.groupby("site")[bias_cols].agg(["count", "mean", "median"])
              .round(2))
        print("\nsat - EC bias by irrig_bucket (TzM only):")
        print(df[df.site == "TzM"].groupby("irrig_bucket")[bias_cols]
              .agg(["count", "mean", "median"]).round(2))


if __name__ == "__main__":
    main()
