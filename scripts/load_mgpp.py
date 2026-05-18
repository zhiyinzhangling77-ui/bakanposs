"""Load LSA SAF Meteosat MGPP (10-day composite GPP) at Oran and TzM.

The MGPP product is a 10-day composite (decadal) over the MSG disk. Each
NetCDF file's timestamp marks the start of the 10-day window. We extract
the nearest-pixel value at each EC tower and write a per-decade CSV plus
a daily-forward-filled CSV for plotting alongside EC.

Usage:
    python3 scripts/load_mgpp.py            # process all years found
    python3 scripts/load_mgpp.py 2018 2019  # specific years
    python3 scripts/load_mgpp.py --force    # re-do everything

After it finishes:
    git add data/mgpp_decadal_all.csv data/mgpp_daily_all.csv
    git commit -m "add Meteosat MGPP extracts"
    git push
"""
from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd
import xarray as xr

BASE = Path("/home/shion-nagamine/Dataset/MGPP")
REPO = Path(__file__).parent.parent
OUT_DIR = REPO / "data" / "mgpp_decadal"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_DECADAL = REPO / "data" / "mgpp_decadal_all.csv"
OUT_DAILY   = REPO / "data" / "mgpp_daily_all.csv"

SITES = {
    "Oran": (38.82, -1.86),
    "TzM":  (39.2660, -1.9397),
}

# Probable variable names inside MSG MGPP NetCDFs (LSA SAF naming varies).
# We try them in order. Override here if your files use a different name.
GPP_VAR_CANDIDATES = ["GPP 10D", "MGPP", "GPP", "gpp"]
QFLAG_CANDIDATES   = ["quality_flag", "qflag", "QF", "QC"]


def parse_timestamp(p: Path) -> pd.Timestamp | None:
    """Parse YYYYMMDDHHMM at the tail of the filename stem.

    Example: NETCDF4_LSASAF_MSG_MGPP_MSG-Disk_201801010000  ->  2018-01-01.
    """
    stem = p.stem
    tail = stem.split("_")[-1]
    if len(tail) >= 8 and tail[:8].isdigit():
        d = tail[:8]
        return pd.Timestamp(f"{d[:4]}-{d[4:6]}-{d[6:8]}")
    return None


def find_nearest_indices(ds: xr.Dataset) -> dict:
    """Return {site: (i_y, i_x)} using whichever coord names the file uses."""
    if "lat" in ds.coords and "lon" in ds.coords:
        lats = ds["lat"].values
        lons = ds["lon"].values
        if lats.ndim == 1 and lons.ndim == 1:
            out = {}
            for name, (lat, lon) in SITES.items():
                i_y = int(np.argmin(np.abs(lats - lat)))
                i_x = int(np.argmin(np.abs(lons - lon)))
                out[name] = (i_y, i_x)
                print(f"  {name}: grid ({lats[i_y]:.3f}, {lons[i_x]:.3f})"
                      f"  target ({lat:.3f}, {lon:.3f})")
            return out
    # fall back: 2-D lat/lon arrays (typical for MSG geostationary projection)
    lats2 = ds["lat"].values if "lat" in ds else ds["latitude"].values
    lons2 = ds["lon"].values if "lon" in ds else ds["longitude"].values
    out = {}
    for name, (lat, lon) in SITES.items():
        d2 = (lats2 - lat) ** 2 + (lons2 - lon) ** 2
        i_y, i_x = np.unravel_index(np.nanargmin(d2), d2.shape)
        out[name] = (int(i_y), int(i_x))
        print(f"  {name}: grid ({lats2[i_y, i_x]:.3f}, {lons2[i_y, i_x]:.3f})"
              f"  target ({lat:.3f}, {lon:.3f})")
    return out


def pick_var(ds: xr.Dataset, candidates: list[str]) -> str | None:
    for v in candidates:
        if v in ds.variables:
            return v
    return None


def process_year(year: int, indices: dict | None) -> tuple[pd.DataFrame, dict]:
    year_dir = BASE / str(year)
    files = sorted(year_dir.rglob("*.nc"))
    if not files:
        print(f"  {year}: no files")
        return pd.DataFrame(), indices or {}
    print(f"  {year}: {len(files):,} files")

    if indices is None:
        with xr.open_dataset(files[0], engine="h5netcdf") as ds0:
            indices = find_nearest_indices(ds0)
            var_g = pick_var(ds0, GPP_VAR_CANDIDATES)
            print(f"  GPP variable: {var_g!r}"
                  f"  (available: {list(ds0.data_vars)})")
            if var_g is None:
                raise RuntimeError(
                    f"None of {GPP_VAR_CANDIDATES} found in {files[0].name}. "
                    f"Edit GPP_VAR_CANDIDATES in this script.")
            attrs = ds0[var_g].attrs
            print(f"  attrs: units={attrs.get('units')!r} "
                  f"scale_factor={attrs.get('scale_factor')} "
                  f"add_offset={attrs.get('add_offset')} "
                  f"_FillValue={attrs.get('_FillValue')} "
                  f"missing_value={attrs.get('missing_value')}")

    rows = []
    t0 = time.time()
    for k, p in enumerate(files):
        ts = parse_timestamp(p)
        if ts is None:
            continue
        try:
            with xr.open_dataset(p, engine="h5netcdf") as ds:
                var_g = pick_var(ds, GPP_VAR_CANDIDATES)
                var_q = pick_var(ds, QFLAG_CANDIDATES)
                if var_g is None:
                    continue
                gpp = ds[var_g].values
                qf  = ds[var_q].values if var_q else None
                # squeeze leading singleton time dim if present
                if gpp.ndim == 3:
                    gpp = gpp[0]
                if qf is not None and qf.ndim == 3:
                    qf = qf[0]
                for site, (i_y, i_x) in indices.items():
                    v = float(gpp[i_y, i_x])
                    q = float(qf[i_y, i_x]) if qf is not None else np.nan
                    rows.append({"date": ts, "site": site,
                                 "GPP_mgpp": v, "qflag": q})
        except Exception as e:
            print(f"\n    [skip] {p.name}: {e}")
            continue
        if (k + 1) % 50 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (k + 1) * (len(files) - k - 1)
            print(f"    [{k+1:>4}/{len(files)}]  elapsed {elapsed:.1f}s "
                  f"ETA {eta:.1f}s")

    return pd.DataFrame(rows), indices


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
        out_csv = OUT_DIR / f"mgpp_decadal_{y}.csv"
        if out_csv.exists() and not force:
            print(f"\n[{y}] already done -> {out_csv}")
            continue
        print(f"\n[{y}] processing ...")
        decadal, indices = process_year(y, indices)
        if not decadal.empty:
            decadal.to_csv(out_csv, index=False)
            print(f"  wrote {out_csv}  ({len(decadal)} rows)")

    # combine all years
    parts = [pd.read_csv(p, parse_dates=["date"])
             for p in sorted(OUT_DIR.glob("mgpp_decadal_*.csv"))]
    if not parts:
        print("no decadal CSVs produced; nothing to combine.")
        return
    decadal_all = pd.concat(parts, ignore_index=True).sort_values(
        ["site", "date"]).reset_index(drop=True)
    decadal_all.to_csv(OUT_DECADAL, index=False)
    print(f"\nwrote {OUT_DECADAL}  ({len(decadal_all):,} rows)")

    # expand 10-day composites to daily (forward-fill the value across
    # the 10-day window) for easy joining with daily EC time series.
    daily_parts = []
    for site, g in decadal_all.groupby("site"):
        g = g.sort_values("date").reset_index(drop=True)
        if g.empty:
            continue
        idx = pd.date_range(g["date"].min(),
                             g["date"].max() + pd.Timedelta(days=9),
                             freq="D")
        s = (g.set_index("date")["GPP_mgpp"].reindex(idx).ffill(limit=9))
        q = (g.set_index("date")["qflag"].reindex(idx).ffill(limit=9))
        daily_parts.append(pd.DataFrame({
            "date": idx, "site": site,
            "GPP_mgpp": s.values, "qflag": q.values,
        }))
    if daily_parts:
        daily_all = pd.concat(daily_parts, ignore_index=True)
        daily_all.to_csv(OUT_DAILY, index=False)
        print(f"wrote {OUT_DAILY}  ({len(daily_all):,} site-days)")

    print("\nannual mean MGPP (g C m-2 d-1 if units already daily):")
    decadal_all["year"] = decadal_all["date"].dt.year
    print(decadal_all.pivot_table(index="year", columns="site",
                                   values="GPP_mgpp", aggfunc="mean").round(3))
    print("\nNOTE: verify MGPP units (kg C m-2 s-1 vs g C m-2 d-1) in the "
          "NetCDF metadata before plotting. Apply scale factor if needed.")


if __name__ == "__main__":
    main(sys.argv[1:])
