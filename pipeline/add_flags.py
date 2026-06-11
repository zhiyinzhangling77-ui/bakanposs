"""Add analysis flags to ec_daily_master.csv.

Builds on the reframed hypothesis: irrigation buffers SWC-ET coupling on a
~3-4 day timescale. The added columns let us stratify any downstream
EC-vs-satellite comparison by the same buckets used in the v14 analysis.

Adds:
  days_since_irrig    days since last irrigation event (Irrig_mm > 0.5).
                      NaN if no prior event in record. 0 on the event itself.
  irrig_bucket        categorical: irrig_d0-3 / d4-7 / d8plus / no_irrig
  season              winter (DJF) / spring (MAM) / summer (JJA) / fall (SON)
  swc_q               SWC quartile within site (dry/low/mid/wet)
  vpd_q               VPD quartile within site (lo/mid/hi/xhi)
  drought_class       normal / atm / soil / compound  (per-site quantile def)
                      atm:      VPD > p75 & SWC > p25
                      soil:     SWC < p25 & VPD < p75
                      compound: VPD > p75 & SWC < p25
                      normal:   else
  is_growing          True if NDVI > 0.3 (rough growth-period mask)

Reads: ec_daily_master.csv
Writes: ec_daily_master_flagged.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).parent.parent
CSV_IN = REPO / "ec_daily_master.csv"
CSV_OUT = REPO / "ec_daily_master_flagged.csv"

IRRIG_THRESHOLD = 0.5  # mm; below this treated as zero


def add_days_since_irrig(df: pd.DataFrame) -> pd.DataFrame:
    """Per-site running counter of days since last irrigation event."""
    out = []
    for site, g in df.groupby("site", sort=False):
        g = g.sort_values("date").copy()
        is_event = (g["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD).to_numpy()
        days = np.full(len(g), np.nan)
        last = -1
        for i, ev in enumerate(is_event):
            if ev:
                last = i
                days[i] = 0
            elif last >= 0:
                days[i] = i - last
        g["days_since_irrig"] = days
        out.append(g)
    return pd.concat(out).sort_index()


def bucket_irrig(d: float) -> str:
    if np.isnan(d):
        return "no_irrig"
    if d <= 3:
        return "irrig_d0-3"
    if d <= 7:
        return "irrig_d4-7"
    return "irrig_d8plus"


def add_season(df: pd.DataFrame) -> pd.DataFrame:
    m = df["date"].dt.month
    season = np.select(
        [m.isin([12, 1, 2]), m.isin([3, 4, 5]),
         m.isin([6, 7, 8]), m.isin([9, 10, 11])],
        ["winter", "spring", "summer", "fall"],
        default="other",
    )
    df["season"] = season
    return df


def add_quantile_flags(df: pd.DataFrame) -> pd.DataFrame:
    df["swc_q"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df["vpd_q"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df["drought_class"] = pd.Series("normal", index=df.index, dtype="object")
    for site, g in df.groupby("site", sort=False):
        idx = g.index
        swc = g["SWC"]
        vpd = g["VPD_kPa_mean"]
        df.loc[idx, "swc_q"] = pd.qcut(
            swc, 4, labels=["dry", "low", "mid", "wet"], duplicates="drop"
        ).astype("object")
        df.loc[idx, "vpd_q"] = pd.qcut(
            vpd, 4, labels=["lo", "mid", "hi", "xhi"], duplicates="drop"
        ).astype("object")
        swc_p25 = swc.quantile(0.25)
        vpd_p75 = vpd.quantile(0.75)
        atm  = (vpd >  vpd_p75) & (swc >= swc_p25)
        soil = (vpd <= vpd_p75) & (swc <  swc_p25)
        comp = (vpd >  vpd_p75) & (swc <  swc_p25)
        cls = pd.Series("normal", index=g.index, dtype="object")
        cls[atm] = "atm"
        cls[soil] = "soil"
        cls[comp] = "compound"
        df.loc[idx, "drought_class"] = cls
    return df


def add_growing_mask(df: pd.DataFrame) -> pd.DataFrame:
    df["is_growing"] = df["NDVI"] > 0.3
    return df


def summarize(df: pd.DataFrame) -> None:
    print("\nirrigation buckets per site:")
    print(df.groupby(["site", "irrig_bucket"]).size().unstack(fill_value=0))

    print("\nseason × site:")
    print(df.groupby(["site", "season"]).size().unstack(fill_value=0))

    print("\ndrought_class × site:")
    print(df.groupby(["site", "drought_class"]).size().unstack(fill_value=0))

    print("\nTzM ET by irrig_bucket (the reframed result):")
    t = df[df.site == "TzM"].dropna(subset=["ET_mm", "SWC"])
    print(t.groupby("irrig_bucket")["ET_mm"].agg(["count", "mean", "std"]))

    print("\nTzM corr(SWC, ET_mm) by irrig_bucket:")
    for bk, g in t.groupby("irrig_bucket"):
        if len(g) > 10:
            print(f"  {bk:14s}  n={len(g):4d}  r={g.SWC.corr(g.ET_mm):+.3f}")

    print("\nOran ET by season (no irrigation):")
    o = df[df.site == "Oran"].dropna(subset=["ET_mm", "SWC"])
    print(o.groupby("season")["ET_mm"].agg(["count", "mean", "std"]))


def main() -> None:
    df = pd.read_csv(CSV_IN, parse_dates=["date"])
    df = add_days_since_irrig(df)
    df["irrig_bucket"] = df["days_since_irrig"].apply(bucket_irrig)
    df = add_season(df)
    df = add_quantile_flags(df)
    df = add_growing_mask(df)
    df = df.sort_values(["site", "date"]).reset_index(drop=True)
    df.to_csv(CSV_OUT, index=False)
    print(f"wrote {CSV_OUT}  rows={len(df)}")
    summarize(df)


if __name__ == "__main__":
    main()
