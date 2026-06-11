"""Inspect EC data files: print shape, columns, head and dtypes for each."""
from pathlib import Path
import pandas as pd

BASE = Path("/home/shion-nagamine/Dataset/Eddy data in Spain")

FILES = [
    "Daily_Summary_Filtered_forPred_ActEne26.csv",
    "EddyAlmond_Raw5years_withG.csv",
    "Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv",
    "Oran_EddyDaily_MASTER_2018_2020_correct.csv",
]


def resolve(name: str) -> Path:
    p = BASE / name
    if p.exists():
        return p
    stem = name.replace(".csv", "")
    for cand in BASE.glob(stem + "*"):
        return cand
    return p


def read_head(p: Path, n: int = 5) -> pd.DataFrame:
    try:
        return pd.read_csv(p, nrows=n)
    except Exception:
        return pd.read_csv(p, nrows=n, sep=None, engine="python")


def main() -> None:
    for name in FILES:
        p = resolve(name)
        print("\n" + "=" * 78)
        print(p.name)
        print("=" * 78)
        if not p.exists():
            print("MISSING:", p)
            continue
        df = read_head(p)
        print("shape (head only):", df.shape)
        print("columns:", list(df.columns))
        print("dtypes:")
        print(df.dtypes.to_string())
        print("first 3 rows:")
        print(df.head(3).to_string())


if __name__ == "__main__":
    main()
