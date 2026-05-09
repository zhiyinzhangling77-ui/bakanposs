"""Inspect satellite CSV files downloaded from GEE."""
from pathlib import Path
import pandas as pd

BASE = Path("/home/shion-nagamine/Downloads/drive-download-20260509T065316Z-3-001")

FILES = [
    "Fast_OranTzM_CHIRPS.csv",
    "Fast_OranTzM_ERA5_Daily.csv",
    "Fast_OranTzM_LAI.csv",
    "Fast_OranTzM_LST.csv",
    "Fast_OranTzM_MOD16.csv",
    "Fast_OranTzM_PML.csv",
    "OranTzM_S2_NDVI.csv",
    "PML_V2_TimeSeries_Oran_TzM.csv",
]


def main() -> None:
    for f in FILES:
        p = BASE / f
        print("\n" + "=" * 78)
        print(p.name)
        print("=" * 78)
        if not p.exists():
            print("MISSING")
            continue
        try:
            df = pd.read_csv(p, nrows=5)
        except Exception:
            df = pd.read_csv(p, nrows=5, sep=None, engine="python")
        print("shape (head only):", df.shape)
        print("columns:", list(df.columns))
        print(df.head(3).to_string())

        # Try to find date column and infer total row count
        full = pd.read_csv(p, low_memory=False)
        print(f"\nfull shape: {full.shape}")
        for c in full.columns:
            if any(k in c.lower() for k in ["date", "time", "system:index"]):
                vals = full[c].dropna().astype(str)
                if len(vals):
                    print(f"  {c}: first={vals.iloc[0]} last={vals.iloc[-1]}")


if __name__ == "__main__":
    main()
