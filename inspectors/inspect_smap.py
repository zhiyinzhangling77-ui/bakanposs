"""Inspect SMAP 3-hour CSV files."""
from pathlib import Path
import pandas as pd

BASE = Path("/mnt/hdd/Dataset/SMAP_3hr")


def main() -> None:
    print(f"Looking under: {BASE}")
    if not BASE.exists():
        print("missing")
        return
    files = sorted(BASE.glob("*.csv"))
    print(f"\nfound {len(files)} CSVs:")
    for f in files:
        print(f"  {f.name}  size={f.stat().st_size / 1e6:.1f} MB")

    for f in files[:2]:  # only first two for inspection
        print("\n" + "=" * 78)
        print(f.name)
        print("=" * 78)
        df = pd.read_csv(f, nrows=5)
        print("shape (head only):", df.shape)
        print("columns:", list(df.columns))
        print(df.head(5).to_string())

        full = pd.read_csv(f, low_memory=False)
        print(f"\nfull shape: {full.shape}")
        for c in full.columns:
            if any(k in c.lower() for k in ["date", "time", "ts"]):
                vals = full[c].dropna().astype(str)
                if len(vals):
                    print(f"  {c}: first={vals.iloc[0]}  last={vals.iloc[-1]}")
        # Site identifiers
        for c in ["name", "site", "lat", "lon", "latitude", "longitude"]:
            if c in full.columns:
                print(f"  {c} unique: {sorted(full[c].dropna().unique().tolist())[:10]}")


if __name__ == "__main__":
    main()
