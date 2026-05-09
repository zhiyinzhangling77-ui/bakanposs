"""Load LSA SAF Meteosat ETv3 at Oran and TzM, aggregate to daily mm/day.

Each NetCDF file holds one 30-min instantaneous ET field [mm/h] on a
global 0.05 deg grid. Daily ET (mm/day) = sum_over_48(ET_i * 0.5 h).
We require >=36 of 48 timesteps per day (75% coverage).

Resumable: writes per-year CSVs into ./metv3_daily/ and the combined
metv3_daily_all.csv at the end. Re-running skips years already done
unless --force is set.

Usage:
    python3 scripts/load_metv3.py            # process all years
    python3 scripts/load_metv3.py 2020 2021  # specific years
    python3 scripts/load_metv3.py --force    # re-do everything
"""
from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd
import xarray as xr

BASE = Path("/home/shion-nagamine/Dataset/METv3")
OUT_DIR = Path("/home/shion-nagamine/bakanposs/metv3_daily")
OUT_DIR.mkdir(exist_ok=True)
OUT_ALL = Path("/home/shion-nagamine/bakanposs/metv3_daily_all.csv")

SITES = {
    "Oran": (38.82, -1.86),
    "TzM":  (39.2660, -1.9397),
}
MIN_OBS_PER_DAY = 36   # of 48 half-hours
DT_HOURS = 0.5


def find_nearest_indices(ds: xr.Dataset) -> dict:
    """Return {site: (i_lat, i_lon)} from the first dataset."""
    lats = ds["lat"].values
    lons = ds["lon"].values
    out = {}
    for name, (lat, lon) in SITES.items():
        i_lat = int(np.argmin(np.abs(lats - lat)))
        i_lon = int(np.argmin(np.abs(lons - lon)))
        out[name] = (i_lat, i_lon)
        print(f"  {name}: nearest grid {lats[i_lat]:.3f}, {lons[i_lon]:.3f} "
              f"(target {lat:.3f}, {lon:.3f})")
    return out


def parse_timestamp(p: Path) -> pd.Timestamp | None:
    """Parse YYYYMMDD_HHMM from filename."""
    stem = p.stem  # 20180101_0000
    try:
        d, t = stem.split("_")
        return pd.Timestamp(f"{d[:4]}-{d[4:6]}-{d[6:8]} {t[:2]}:{t[2:4]}")
    except Exception:
        return None


def process_year(year: int, indices: dict | None) -> tuple[pd.DataFrame, dict]:
    year_dir = BASE / str(year)
    files = sorted(year_dir.rglob("*.nc"))
    if not files:
        print(f"  {year}: no files")
        return pd.DataFrame(), indices or {}

    print(f"  {year}: {len(files):,} files")

    # establish nearest-pixel indices once
    if indices is None:
        with xr.open_dataset(files[0], engine="h5netcdf") as ds0:
            indices = find_nearest_indices(ds0)

    rows = []
    t0 = time.time()
    for k, p in enumerate(files):
        ts = parse_timestamp(p)
        if ts is None:
            continue
        try:
            with xr.open_dataset(p, engine="h5netcdf") as ds:
                et = ds["ET"].values
                qf = ds["quality_flag"].values if "quality_flag" in ds else None
                for site, (i_lat, i_lon) in indices.items():
                    v = float(et[0, i_lat, i_lon])
                    q = float(qf[0, i_lat, i_lon]) if qf is not None else np.nan
                    rows.append({"ts": ts, "site": site,
                                 "ET_mmh": v, "qflag": q})
        except Exception as e:
            print(f"\n    [skip] {p.name}: {e}")
            continue

        if (k + 1) % 500 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (k + 1) * (len(files) - k - 1)
            print(f"    [{k+1:>6,}/{len(files):,}]  "
                  f"elapsed {elapsed/60:.1f} min  ETA {eta/60:.1f} min")

    inst = pd.DataFrame(rows)
    if inst.empty:
        return inst, indices

    # daily aggregation
    inst["date"] = inst["ts"].dt.normalize()
    inst["ET_mm_per_30min"] = inst["ET_mmh"].clip(lower=0) * DT_HOURS
    daily = inst.groupby(["date", "site"]).agg(
        ET_mm=("ET_mm_per_30min", "sum"),
        n_obs=("ET_mm_per_30min", "count"),
        qflag_mean=("qflag", "mean"),
    ).reset_index()
    daily.loc[daily["n_obs"] < MIN_OBS_PER_DAY, "ET_mm"] = np.nan
    return daily, indices


def main(args: list[str]) -> None:
    force = "--force" in args
    args = [a for a in args if a != "--force"]

    if args:
        years = [int(a) for a in args]
    else:
        years = sorted(int(p.name) for p in BASE.iterdir()
                       if p.is_dir() and p.name.isdigit())
    print(f"years to process: {years}  (force={force})")

    indices = None
    for y in years:
        out_csv = OUT_DIR / f"metv3_daily_{y}.csv"
        if out_csv.exists() and not force:
            print(f"\n[{y}] already done -> {out_csv}")
            continue
        print(f"\n[{y}] processing ...")
        daily, indices = process_year(y, indices)
        if not daily.empty:
            daily.to_csv(out_csv, index=False)
            print(f"  wrote {out_csv}  ({len(daily)} site-days)")

    # combine all years
    parts = [pd.read_csv(p, parse_dates=["date"])
             for p in sorted(OUT_DIR.glob("metv3_daily_*.csv"))]
    if parts:
        combined = pd.concat(parts, ignore_index=True)
        combined.to_csv(OUT_ALL, index=False)
        print(f"\nwrote {OUT_ALL}  ({len(combined):,} site-days)")
        print("\nmean daily ET (mm/d) per site, per year:")
        combined["year"] = combined["date"].dt.year
        print(combined.pivot_table(index="year", columns="site",
                                     values="ET_mm", aggfunc="mean").round(2))


if __name__ == "__main__":
    main(sys.argv[1:])
