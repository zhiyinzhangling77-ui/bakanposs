"""Inspect the LSA SAF Meteosat ETv3 NetCDF files to understand structure."""
from pathlib import Path
import xarray as xr

BASE = Path("/home/shion-nagamine/Dataset/METv3")

SAMPLE_PATHS = [
    BASE / "2018/01/0101/20180101_0000.nc",
    BASE / "2018/06/0615/20180615_1200.nc",
    BASE / "2020/07/0715/20200715_1200.nc",
]


def main() -> None:
    print(f"Looking under: {BASE}")
    print(f"Top-level entries: {sorted(p.name for p in BASE.iterdir())[:20]}")

    # Count files in 2018 to confirm cadence
    y = BASE / "2018"
    if y.exists():
        ncs = list(y.rglob("*.nc"))
        print(f"\n2018: {len(ncs):,} .nc files")
        if ncs:
            print(f"first: {ncs[0]}")
            print(f"last : {ncs[-1]}")

    for p in SAMPLE_PATHS:
        if not p.exists():
            print(f"\n[MISSING] {p}")
            continue
        print("\n" + "=" * 78)
        print(p)
        print("=" * 78)
        try:
            ds = xr.open_dataset(p)
        except Exception as e:
            print(f"open failed: {e}")
            try:
                ds = xr.open_dataset(p, engine="h5netcdf")
            except Exception as e2:
                print(f"h5netcdf failed too: {e2}")
                continue
        print(ds)
        for v in ds.data_vars:
            arr = ds[v]
            print(f"\nvar {v}: shape={arr.shape}, dtype={arr.dtype}, "
                  f"attrs={dict(arr.attrs)}")
        print(f"\ncoords: {list(ds.coords)}")
        for c in ds.coords:
            arr = ds[c]
            try:
                print(f"  {c}: range [{arr.min().values}, {arr.max().values}], "
                      f"shape={arr.shape}")
            except Exception:
                print(f"  {c}: shape={arr.shape}")


if __name__ == "__main__":
    main()
