"""Sanity checks on ec_daily_master.csv.

Verifies the assumptions baked into unify_ec_daily.py:
  1. VPD unit: TzM VPD_mean was assumed Pa, divided by 1000 to get kPa
  2. ET magnitude: ET_mm should be 0-10 mm/day at these sites
  3. Crop rotation at Oran: list the unique crops and the date span of each
  4. Period coverage and missing-data per site/year
  5. Energy balance numbers (Rn vs LE+H+G) where available
  6. Irrigation events at TzM
  7. SWC vs ET decoupling preview (drought signal sniff)

No figures, only printed tables. Run AFTER unify_ec_daily.py.
"""
from pathlib import Path
import numpy as np
import pandas as pd

CSV = Path("/home/shion-nagamine/bakanposs/ec_daily_master.csv")
RAW_TZM = Path(
    "/home/shion-nagamine/Dataset/Eddy data in Spain/"
    "Daily_Summary_Filtered_forPred_ActEne26.csv"
)


def hr(title: str) -> None:
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78)


def check_vpd_unit() -> None:
    """Re-open the raw TzM file and compare VPD_mean vs VPD_kPa directly."""
    hr("1. TzM VPD unit check (Pa vs kPa)")
    raw = pd.read_csv(RAW_TZM, usecols=["VPD_mean", "VPD_kPa", "VPD_max"])
    print("VPD_mean : ", raw.VPD_mean.describe()[["min", "50%", "max"]].to_dict())
    print("VPD_kPa  : ", raw.VPD_kPa.describe()[["min", "50%", "max"]].to_dict())
    ratio = (raw.VPD_mean / raw.VPD_kPa).replace([np.inf, -np.inf], np.nan).dropna()
    print(f"ratio VPD_mean / VPD_kPa : median={ratio.median():.3f}  "
          f"mean={ratio.mean():.3f}")
    if 900 < ratio.median() < 1100:
        verdict = "VPD_mean is in Pa (divide by 1000 → kPa). Current unify is correct."
    elif 0.9 < ratio.median() < 1.1:
        verdict = ("VPD_mean is ALREADY in kPa! Remove the /1000 in load_tzm. "
                   "Current unify gives values 1000x too small.")
    else:
        verdict = f"unexpected ratio {ratio.median():.3f}; inspect manually."
    print("VERDICT:", verdict)


def check_et_range(df: pd.DataFrame) -> None:
    hr("2. ET range per site")
    print(df.groupby("site")["ET_mm"].describe()[["count", "mean", "min", "max"]])
    high = df[df.ET_mm > 12]
    if len(high):
        print(f"\n!! {len(high)} rows with ET>12 mm/d (likely outlier):")
        print(high[["date", "site", "ET_mm", "Rn_Wm2", "VPD_kPa_mean"]].head())
    neg = df[df.ET_mm < 0]
    if len(neg):
        print(f"\n!! {len(neg)} rows with ET<0:")
        print(neg[["date", "site", "ET_mm"]].head())


def check_oran_crops(df: pd.DataFrame) -> None:
    hr("3. Oran crop rotation")
    o = df[df.site == "Oran"].copy()
    print("unique crops:", o.crop.unique())
    spans = o.groupby("crop")["date"].agg(["min", "max", "count"])
    print(spans)
    print("\nyearly crop counts:")
    print(o.groupby([o.date.dt.year, "crop"]).size().unstack(fill_value=0))


def check_coverage(df: pd.DataFrame) -> None:
    hr("4. Period coverage")
    cov = df.groupby(["site", df.date.dt.year]).agg(
        n_days=("date", "count"),
        ET_n=("ET_mm", "count"),
        Rn_n=("Rn_Wm2", "count"),
        SWC_n=("SWC", "count"),
        NDVI_n=("NDVI", "count"),
    )
    print(cov)
    print("\nmissing per site per column:")
    miss = df.groupby("site").apply(lambda g: g.isna().sum())
    print(miss[["ET_mm", "LE_Wm2", "H_Wm2", "G_Wm2", "Rn_Wm2",
                "VPD_kPa_mean", "SWC", "NDVI"]])


def check_energy_balance(df: pd.DataFrame) -> None:
    """Available LE+H+G vs Rn for TzM (Oran daily lacks LE/H/G)."""
    hr("5. Energy balance closure (where flux daily means available)")
    s = df.dropna(subset=["LE_Wm2", "H_Wm2", "G_Wm2", "Rn_Wm2"]).copy()
    if len(s) == 0:
        print("no daily flux means available")
        return
    s["turb"] = s.LE_Wm2 + s.H_Wm2 + s.G_Wm2
    s["EBR"] = s.turb / s.Rn_Wm2.replace(0, np.nan)
    print(s.groupby("site")["EBR"].describe()[["count", "mean", "50%"]])
    print("EBR ~0.7-0.9 is typical; <0.5 → check sign/unit.")


def check_irrigation(df: pd.DataFrame) -> None:
    hr("6. Irrigation events at TzM")
    t = df[(df.site == "TzM") & (df.Irrig_mm > 0)]
    print(f"#irrigation days = {len(t)}")
    if len(t):
        print(t.Irrig_mm.describe()[["count", "mean", "max"]].to_dict())
        print("monthly irrigation totals:")
        print(t.groupby(t.date.dt.to_period("M"))["Irrig_mm"].sum().head(20))


def check_swc_et_decoupling(df: pd.DataFrame) -> None:
    hr("7. SWC-ET decoupling preview (drought signal)")
    for site in ["Oran", "TzM"]:
        s = df[df.site == site].dropna(subset=["SWC", "ET_mm"])
        if len(s) < 30:
            continue
        # split by SWC quantile
        s["swc_q"] = pd.qcut(s.SWC, 4, labels=["dry", "low", "mid", "wet"],
                             duplicates="drop")
        agg = s.groupby("swc_q", observed=True)["ET_mm"].agg(["count", "mean", "max"])
        print(f"\n{site}")
        print(agg)
        print(f"corr(SWC, ET) = {s.SWC.corr(s.ET_mm):+.3f}")


def main() -> None:
    check_vpd_unit()
    df = pd.read_csv(CSV, parse_dates=["date"])
    check_et_range(df)
    check_oran_crops(df)
    check_coverage(df)
    check_energy_balance(df)
    check_irrigation(df)
    check_swc_et_decoupling(df)


if __name__ == "__main__":
    main()
