"""Integrate METv3 and SMAP into master_full.csv.

Run AFTER:
  - load_metv3.py  → metv3_daily_all.csv
  - load_smap.py   → smap_daily.csv
  - merge_satellite_ec.py  → master_full.csv (MOD16 + PML already merged)

This script left-joins METv3 and SMAP onto master_full.csv, reports
bias statistics for all three satellite ET products, and writes
master_full_v2.csv.

Usage:
    python3 scripts/integrate_metv3_smap.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).parent.parent
MASTER = REPO / "master_full.csv"
METV3 = REPO / "metv3_daily_all.csv"
SMAP  = REPO / "smap_daily.csv"
OUT   = REPO / "master_full_v2.csv"


# ── helpers ──────────────────────────────────────────────────────────────────

def bias_table(df: pd.DataFrame, sat_col: str, ec_col: str,
               label: str) -> None:
    d = df.dropna(subset=[sat_col, ec_col]).copy()
    d["delta"] = d[sat_col] - d[ec_col]
    print(f"\n{'─'*60}")
    print(f"{label}  (n={len(d):,})")
    print(f"  overall MBE = {d['delta'].mean():.3f} mm/d  "
          f"RMSE = {np.sqrt((d['delta']**2).mean()):.3f}")
    for site in ["Oran", "TzM"]:
        s = d[d["site"] == site]
        if len(s) < 5:
            continue
        print(f"  {site} (n={len(s):,})  MBE={s['delta'].mean():.3f}  "
              f"RMSE={np.sqrt((s['delta']**2).mean()):.3f}")

    # irrigation bucket
    if "irrig_bucket" in d.columns:
        tzm = d[d["site"] == "TzM"]
        is_sum = tzm.get("season", pd.Series(dtype=str)) == "JJA"
        is_ndvi = tzm.get("is_growing", pd.Series(dtype=bool))
        sub = tzm[is_sum & is_ndvi] if (is_sum & is_ndvi).any() else tzm
        if len(sub) > 10:
            print(f"  TzM summer×NDVI>0.3 by irrigation bucket:")
            for bkt in ["d0-3", "d4-7", "d8+"]:
                g = sub[sub["irrig_bucket"] == bkt]
                if len(g) < 3:
                    continue
                print(f"    {bkt}  n={len(g):>3d}  "
                      f"MBE={g['delta'].mean():.3f}  "
                      f"std={g['delta'].std():.3f}")


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    # load master
    if not MASTER.exists():
        raise FileNotFoundError(f"{MASTER} not found — run merge_satellite_ec.py first")
    master = pd.read_csv(MASTER, parse_dates=["date"])
    print(f"master_full  : {len(master):,} rows  cols={list(master.columns)[:8]}...")

    # ── METv3 ──────────────────────────────────────────────────────────────
    if METV3.exists():
        metv3 = pd.read_csv(METV3, parse_dates=["date"])
        print(f"\nMETv3 all    : {len(metv3):,} site-days")
        print(f"  date range : {metv3['date'].min()} → {metv3['date'].max()}")
        print(f"  sites      : {sorted(metv3['site'].unique())}")
        print(f"  non-null ET: {metv3['ET_mm'].notna().sum():,}  "
              f"({100*metv3['ET_mm'].notna().mean():.1f}%)")
        # site name harmonisation
        metv3 = metv3.rename(columns={"ET_mm": "ET_metv3_mm",
                                       "n_obs": "metv3_n_obs",
                                       "qflag_mean": "metv3_qflag"})
        master = master.merge(
            metv3[["date", "site", "ET_metv3_mm", "metv3_n_obs", "metv3_qflag"]],
            on=["date", "site"], how="left"
        )
        print(f"  joined     : {master['ET_metv3_mm'].notna().sum():,} non-null rows")
    else:
        print(f"\n[SKIP] {METV3} not found — run load_metv3.py first")

    # ── SMAP ───────────────────────────────────────────────────────────────
    if SMAP.exists():
        smap = pd.read_csv(SMAP, parse_dates=["date"])
        print(f"\nSMAP daily   : {len(smap):,} site-days")
        print(f"  date range : {smap['date'].min()} → {smap['date'].max()}")
        print(f"  non-null surface : {smap['sm_surface_mean'].notna().sum():,}  "
              f"({100*smap['sm_surface_mean'].notna().mean():.1f}%)")
        print(f"  non-null rootzone: {smap['sm_rootzone_mean'].notna().sum():,}  "
              f"({100*smap['sm_rootzone_mean'].notna().mean():.1f}%)")
        master = master.merge(
            smap[["date", "site", "sm_surface_mean", "sm_rootzone_mean",
                   "n_obs_surface", "n_obs_rootzone"]].rename(
                columns={"sm_surface_mean": "smap_surface",
                         "sm_rootzone_mean": "smap_rootzone",
                         "n_obs_surface": "smap_n_surface",
                         "n_obs_rootzone": "smap_n_rootzone"}),
            on=["date", "site"], how="left"
        )
        print(f"  joined SMAP surface : {master['smap_surface'].notna().sum():,} non-null")
    else:
        print(f"\n[SKIP] {SMAP} not found — run load_smap.py first")

    # ── bias summary ────────────────────────────────────────────────────────
    ec_col = "ET_mm" if "ET_mm" in master.columns else "ET_EC_mm"
    # try to find EC ET column
    for c in ["ET_mm", "ET_EC_mm", "et_ec"]:
        if c in master.columns and master[c].notna().sum() > 100:
            ec_col = c
            break

    for sat, col in [("MOD16", "mod16_et"), ("PML", "pml_et"),
                     ("METv3", "ET_metv3_mm")]:
        if col in master.columns:
            bias_table(master, col, ec_col, f"{sat} vs EC")

    # ── write output ────────────────────────────────────────────────────────
    master.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}  ({len(master):,} rows, {master.shape[1]} cols)")

    # column inventory
    print("\ncolumn inventory:")
    for c in sorted(master.columns):
        nn = master[c].notna().sum()
        print(f"  {c:<30s}  {nn:>6,} non-null  ({100*nn/len(master):.0f}%)")


if __name__ == "__main__":
    main()
