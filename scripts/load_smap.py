"""Load SMAP 3-hour soil moisture, aggregate to daily means at Oran/TzM.

Accepts either:
  (A) the new single-file GEE export at /mnt/hdd/Dataset/SMAP_OranTzM.csv
      (columns: name, date_time, sm_surface, sm_rootzone)
  (B) the old per-year files at /mnt/hdd/Dataset/SMAP_3hr/SMAP_3hr_YYYY.csv

Cadence: 8 obs/day (01:30, 04:30, ..., 22:30).

Writes: smap_daily.csv  (date, site, sm_surface_mean, sm_rootzone_mean,
                         n_obs_surface, n_obs_rootzone)
"""
from pathlib import Path
import pandas as pd

NEW_FILE = Path("/mnt/hdd/Dataset/SMAP_OranTzM.csv")
OLD_DIR  = Path("/mnt/hdd/Dataset/SMAP_3hr")
REPO = Path(__file__).parent.parent
OUT = REPO / "smap_daily.csv"
MIN_OBS_PER_DAY = 4   # of 8 (50%)


def load_raw() -> pd.DataFrame:
    if NEW_FILE.exists():
        print(f"reading single GEE export: {NEW_FILE}")
        df = pd.read_csv(NEW_FILE, parse_dates=["date_time"])
        df["site"] = df["name"]
        return df

    files = sorted(OLD_DIR.glob("SMAP_3hr_*.csv"))
    if not files:
        raise FileNotFoundError("no SMAP files found at either expected path")

    print(f"reading {len(files)} per-year SMAP files from {OLD_DIR}")
    raws = []
    for f in files:
        df = pd.read_csv(f, parse_dates=["date_time"])
        df["site"] = df["name"]
        raws.append(df)
    return pd.concat(raws, ignore_index=True)


def main() -> None:
    raw = load_raw()
    print(f"raw rows: {len(raw):,}")
    print(f"date range: {raw['date_time'].min()} - {raw['date_time'].max()}")
    print(f"non-null surface : {raw['sm_surface'].notna().sum():,}  "
          f"({100 * raw['sm_surface'].notna().mean():.1f}%)")
    print(f"non-null rootzone: {raw['sm_rootzone'].notna().sum():,}  "
          f"({100 * raw['sm_rootzone'].notna().mean():.1f}%)")

    raw["date"] = raw["date_time"].dt.normalize()
    g = raw.groupby(["date", "site"])
    daily = pd.DataFrame({
        "sm_surface_mean":  g["sm_surface"].mean(),
        "sm_rootzone_mean": g["sm_rootzone"].mean(),
        "n_obs_surface":    g["sm_surface"].count(),
        "n_obs_rootzone":   g["sm_rootzone"].count(),
    }).reset_index()

    bad_s = daily["n_obs_surface"] < MIN_OBS_PER_DAY
    bad_r = daily["n_obs_rootzone"] < MIN_OBS_PER_DAY
    daily.loc[bad_s, "sm_surface_mean"] = pd.NA
    daily.loc[bad_r, "sm_rootzone_mean"] = pd.NA

    daily.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}  ({len(daily):,} site-days)")

    print("\nannual mean per site (m^3/m^3):")
    daily["year"] = pd.to_datetime(daily["date"]).dt.year
    print(daily.pivot_table(
        index="year", columns="site",
        values=["sm_surface_mean", "sm_rootzone_mean"], aggfunc="mean"
    ).round(3))


if __name__ == "__main__":
    main()
