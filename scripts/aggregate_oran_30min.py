"""Aggregate Oran 30-min EC fluxes to daily means and merge into the master.

The Oran daily MASTER file lacks LE/H/G daily means (only NETRAD aggregates),
so reading the 30-min CLEAN file ourselves is the only way to get the energy
balance closure (EBR) for Oran.

QC convention in this file (AmeriFlux): LE_QC, H_QC -> 1=best, 2=ok, 3=bad.
We keep QC <= 2 and require >=24 valid half-hours per day (50% coverage).

Reads:
  Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv  (30-min)
  ec_daily_master_flagged.csv                    (current master)

Writes:
  oran_daily_fluxes.csv             (Oran daily LE/H/G/Rn/EBR/qc_frac)
  ec_daily_master_complete.csv      (master with Oran fluxes filled in)
"""
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("/home/shion-nagamine/Dataset/Eddy data in Spain")
RAW = BASE / "Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv"
MASTER_IN = Path("/home/shion-nagamine/bakanposs/ec_daily_master_flagged.csv")
DAILY_OUT = Path("/home/shion-nagamine/bakanposs/oran_daily_fluxes.csv")
MASTER_OUT = Path("/home/shion-nagamine/bakanposs/ec_daily_master_complete.csv")

QC_OK = {1, 2}
MIN_OBS_PER_DAY = 24  # half-hours; 24 = 50% coverage


def load_30min() -> pd.DataFrame:
    """Robust loader: build timestamp from year + Julian + Time_hours.

    The DateTime / TIMESTAMP string columns mix "YYYY/MM/DD" and
    "YYYY/MM/DD HH:MM:SS" so pandas locks onto the date-only format and
    loses 98% of rows. year + Julian + Time_hours are clean numerics that
    are guaranteed present for every half-hour.
    """
    cols = ["year", "Julian", "Time_hours",
            "LE", "LE_QC", "H", "H_QC", "NETRAD", "G",
            "VPD", "TA_1_1_1", "SWC_1_1_1", "P", "energy_closure"]
    df = pd.read_csv(RAW, usecols=cols, low_memory=False)
    print(f"  read {len(df):,} raw rows")

    for col in cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    bad_idx = df[["year", "Julian", "Time_hours"]].isna().any(axis=1)
    if bad_idx.any():
        print(f"  {bad_idx.sum():,} rows had bad year/Julian/Time_hours (dropped)")
    df = df[~bad_idx].copy()

    df["ts"] = (
        pd.to_datetime(df["year"].astype(int).astype(str), format="%Y")
        + pd.to_timedelta(df["Julian"].astype(int) - 1, unit="D")
        + pd.to_timedelta(df["Time_hours"], unit="h")
    )

    df["LE_QC"] = df["LE_QC"].astype("Int64")
    df["H_QC"] = df["H_QC"].astype("Int64")
    df["date"] = df["ts"].dt.normalize()
    print(f"  kept {len(df):,} half-hour rows  "
          f"({df['ts'].min().date()} → {df['ts'].max().date()})")
    return df


def mask_qc(s: pd.Series, qc: pd.Series) -> pd.Series:
    """Keep s where qc is in QC_OK, else NaN."""
    qc_int = pd.to_numeric(qc, errors="coerce")
    return s.where(qc_int.isin(QC_OK))


def aggregate_daily(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["LE_q"] = mask_qc(df.LE, df.LE_QC)
    df["H_q"] = mask_qc(df.H, df.H_QC)

    g = df.groupby("date")
    daily = pd.DataFrame({
        "LE_Wm2": g["LE_q"].mean(),
        "H_Wm2": g["H_q"].mean(),
        "G_Wm2": g["G"].mean(),
        "Rn_Wm2": g["NETRAD"].mean(),
        "Ta_30min_mean": g["TA_1_1_1"].mean(),
        "VPD_30min_mean_kPa": g["VPD"].mean() / 10.0,  # source is hPa
        "SWC_30min_mean": g["SWC_1_1_1"].mean(),
        "n_LE_ok": g["LE_q"].count(),
        "n_H_ok": g["H_q"].count(),
        "n_total": g.size(),
    }).reset_index()

    # require enough observations
    bad = daily["n_LE_ok"] < MIN_OBS_PER_DAY
    daily.loc[bad, ["LE_Wm2", "H_Wm2", "G_Wm2"]] = np.nan

    # energy balance ratio
    turb = daily["LE_Wm2"] + daily["H_Wm2"] + daily["G_Wm2"]
    daily["EBR"] = turb / daily["Rn_Wm2"].replace(0, np.nan)
    daily["qc_frac"] = daily["n_LE_ok"] / daily["n_total"]
    return daily


def merge_into_master(master: pd.DataFrame, oran_daily: pd.DataFrame) -> pd.DataFrame:
    out = master.copy()
    od = oran_daily.copy()
    od["site"] = "Oran"
    od = od.rename(columns={
        "LE_Wm2": "LE_Wm2_new",
        "H_Wm2": "H_Wm2_new",
        "G_Wm2": "G_Wm2_new",
        "Rn_Wm2": "Rn_Wm2_new",
    })
    out["date"] = pd.to_datetime(out["date"])
    od["date"] = pd.to_datetime(od["date"])
    merged = out.merge(
        od[["date", "site", "LE_Wm2_new", "H_Wm2_new", "G_Wm2_new",
            "Rn_Wm2_new", "EBR", "qc_frac"]],
        on=["date", "site"], how="left",
    )
    # Fill Oran flux columns from 30-min aggregation; preserve TzM as-is
    for col in ["LE_Wm2", "H_Wm2", "G_Wm2"]:
        new = merged[col + "_new"]
        merged[col] = merged[col].where(merged[col].notna(), new)
        merged.drop(columns=[col + "_new"], inplace=True)
    # Rn: prefer existing; only fill where missing
    merged["Rn_Wm2"] = merged["Rn_Wm2"].where(merged["Rn_Wm2"].notna(),
                                              merged["Rn_Wm2_new"])
    merged.drop(columns=["Rn_Wm2_new"], inplace=True)
    return merged


def main() -> None:
    print("loading 30-min ...")
    raw = load_30min()
    print(f"  {len(raw):,} rows, {raw['date'].min().date()} → {raw['date'].max().date()}")

    print("\naggregating to daily ...")
    daily = aggregate_daily(raw)
    daily.to_csv(DAILY_OUT, index=False)
    print(f"wrote {DAILY_OUT}  rows={len(daily)}")
    print("\nEBR per year (target ~0.7-0.9):")
    daily["year"] = daily["date"].dt.year
    print(daily.groupby("year")["EBR"].agg(["count", "mean", "median"]))

    print("\nflux means per year (W/m²):")
    print(daily.groupby("year")[["LE_Wm2", "H_Wm2", "G_Wm2", "Rn_Wm2"]].mean()
          .round(1))

    print("\nmerging into master ...")
    master = pd.read_csv(MASTER_IN, parse_dates=["date"])
    print(f"master before: Oran LE non-null = {master[master.site=='Oran']['LE_Wm2'].notna().sum()}")
    merged = merge_into_master(master, daily)
    print(f"master after : Oran LE non-null = {merged[merged.site=='Oran']['LE_Wm2'].notna().sum()}")
    merged.to_csv(MASTER_OUT, index=False)
    print(f"\nwrote {MASTER_OUT}  rows={len(merged)}")


if __name__ == "__main__":
    main()
