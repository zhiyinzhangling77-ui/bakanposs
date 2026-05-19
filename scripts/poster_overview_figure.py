"""Poster overview figure: NDVI, GPP, ET dynamics at two flux towers.

Each indicator is loaded from its INDEPENDENT source so that satellite
lines remain plotted on EC-missing days. master_full_v2.csv (EC-joined)
was discarding satellite values on EC-gap days; we bypass that here.

Sources:
  EC tower fluxes/NDVI      : data/master_full_v2.csv      (EC days only)
  Sentinel-2 NDVI           : data/OranTzM_S2_NDVI.csv     (B4/B8 -> NDVI;
                              full S2 revisit, independent of EC)
  Meteosat MGPP (10-day GPP): data/mgpp_decadal_all.csv    (independent)
  Meteosat ETv3 (daily ET)  : data/metv3_daily_all.csv     (independent)

Layout (3 rows x 2 cols):
  Row 1: NDVI (EC tower vs Sentinel-2)
  Row 2: GPP  (EC tower vs Meteosat MGPP)
  Row 3: ET   (EC tower vs Meteosat ETv3)

Spatial / temporal scales are annotated in each legend entry.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

REPO = Path(__file__).parent.parent
EC_CSV    = REPO / "data" / "master_full_v2.csv"
METV3_CSV = REPO / "data" / "metv3_daily_all.csv"
MGPP_CSV  = REPO / "data" / "mgpp_decadal_all.csv"
S2_CSV    = REPO / "data" / "OranTzM_S2_NDVI.csv"
OUT_DIR = REPO / "figures" / "poster"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SITES = ["Oran", "TzM"]
SITE_LABEL = {
    "Oran": "Oran  ·  rainfed vetch  ·  2018–2020",
    "TzM":  "Tarazona  ·  drip-irrig. almond  ·  2020–2024",
}

C_EC      = "#1a1a1a"
C_EC_RAW  = "#c8c8c8"
C_METV3   = "#d1493c"
C_MGPP    = "#2a8c63"
C_S2      = "#1e6fb3"
SHADE     = "#f6e7b8"
GRID      = "#dddddd"


def smooth(s: pd.Series, w: int = 7) -> pd.Series:
    return s.rolling(w, center=True, min_periods=max(3, w // 2)).mean()


def smooth_break(s: pd.Series, gap_days: int = 7, w: int = 7) -> pd.Series:
    """Centered 7-day mean, but NaN inside any run of >gap_days missing days."""
    out = s.rolling(w, center=True, min_periods=max(3, w // 2)).mean()
    miss = s.isna()
    streak = miss.groupby((~miss).cumsum()).cumsum()
    return out.where(~(streak > gap_days))


def break_gaps(dates, values, max_gap_days: int = 15):
    """Return (dates, values) with NaN inserted where consecutive valid
    observations are more than max_gap_days apart, so matplotlib breaks
    the line rather than bridging gaps."""
    d = pd.to_datetime(pd.Series(dates))
    v = pd.Series(values, index=d.index)
    keep = v.notna() & d.notna()
    d = d[keep].reset_index(drop=True)
    v = v[keep].reset_index(drop=True)
    order = d.argsort()
    d = d.iloc[order].reset_index(drop=True)
    v = v.iloc[order].reset_index(drop=True)
    if len(d) == 0:
        return d.to_numpy(), v.to_numpy()
    gaps = d.diff().dt.days
    out_d, out_v = [], []
    for i in range(len(d)):
        if i > 0 and gaps.iloc[i] > max_gap_days:
            out_d.append(d.iloc[i - 1] + pd.Timedelta(days=1))
            out_v.append(np.nan)
        out_d.append(d.iloc[i])
        out_v.append(v.iloc[i])
    return (pd.to_datetime(out_d).to_numpy(),
            np.asarray(out_v, dtype=float))


def reindex_daily(df: pd.DataFrame, d0, d1) -> pd.DataFrame:
    full = pd.date_range(d0, d1, freq="D")
    return df.set_index("date").reindex(full).rename_axis("date").reset_index()


def shade_springs(ax, date_min, date_max):
    for y in range(date_min.year, date_max.year + 1):
        a = pd.Timestamp(f"{y}-03-01")
        b = pd.Timestamp(f"{y}-06-30")
        if b < date_min or a > date_max:
            continue
        ax.axvspan(max(a, date_min), min(b, date_max),
                   color=SHADE, alpha=0.55, zorder=0)


def style_time_ax(ax, d0, d1):
    ax.set_xlim(d0, d1)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.grid(which="major", color=GRID, lw=0.8, zorder=0)


def load_s2_ndvi() -> pd.DataFrame:
    """Compute NDVI from B4/B8 and aggregate to one value per date+site."""
    raw = pd.read_csv(S2_CSV, parse_dates=["date"])
    raw = raw.rename(columns={"name": "site"})
    raw["NDVI_s2"] = (raw["B8"] - raw["B4"]) / (raw["B8"] + raw["B4"])
    raw = raw.dropna(subset=["NDVI_s2"])
    return (raw.groupby(["site", "date"], as_index=False)["NDVI_s2"]
               .mean())


def main() -> None:
    ec_full    = pd.read_csv(EC_CSV,    parse_dates=["date"])
    metv3_full = pd.read_csv(METV3_CSV, parse_dates=["date"])
    mgpp_full  = pd.read_csv(MGPP_CSV,  parse_dates=["date"])
    mgpp_full  = mgpp_full[(mgpp_full["GPP_mgpp"] > 0) &
                           (mgpp_full["GPP_mgpp"] < 30)].copy()
    s2_full    = load_s2_ndvi()

    fig, axes = plt.subplots(
        3, 2, figsize=(13.5, 9.0),
        gridspec_kw={"hspace": 0.32, "wspace": 0.18},
    )
    panel_labels = [["a", "b"], ["c", "d"], ["e", "f"]]

    for j, site in enumerate(SITES):
        ec_site    = ec_full   [ec_full   ["site"] == site].sort_values("date")
        metv3_site = metv3_full[metv3_full["site"] == site].sort_values("date")
        mgpp_site  = mgpp_full [mgpp_full ["site"] == site].sort_values("date")
        s2_site    = s2_full   [s2_full   ["site"] == site].sort_values("date")

        # x-range = EC observation window per site
        d0, d1 = ec_site["date"].min(), ec_site["date"].max()

        # reindex each source to its OWN daily continuous frame
        ec_d    = reindex_daily(ec_site,    d0, d1)
        metv3_d = reindex_daily(metv3_site[(metv3_site["date"] >= d0) &
                                            (metv3_site["date"] <= d1)],
                                 d0, d1)
        mgpp_d  = mgpp_site[(mgpp_site["date"] >= d0) &
                            (mgpp_site["date"] <= d1)].copy()
        s2_d    = s2_site[(s2_site["date"] >= d0) &
                          (s2_site["date"] <= d1)].copy()

        # ---- Row 1: NDVI ----
        ax = axes[0, j]
        shade_springs(ax, d0, d1)
        s2_x, s2_y = break_gaps(s2_d["date"], s2_d["NDVI_s2"],
                                 max_gap_days=20)
        ax.plot(s2_x, s2_y,
                "-o", color=C_S2, lw=1.0, ms=3.2, mec="white", mew=0.4,
                alpha=0.9, zorder=2,
                label="Sentinel-2 NDVI  ·  10 m  ·  ~5-day")
        ax.plot(ec_d["date"], smooth_break(ec_d["NDVI"]),
                color=C_EC, lw=1.7, zorder=3,
                label="EC tower NDVI  ·  ~100 m  ·  daily (7-d)")
        ax.set_ylabel("NDVI (–)", fontsize=10.5)
        ax.set_ylim(-0.05, 0.85)
        style_time_ax(ax, d0, d1)
        ax.set_title(SITE_LABEL[site], fontsize=12.0, pad=8)
        if j == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=9.5)
        ax.text(-0.085, 1.07, panel_labels[0][j], transform=ax.transAxes,
                fontsize=14, fontweight="bold")

        # ---- Row 2: GPP ----
        ax = axes[1, j]
        shade_springs(ax, d0, d1)
        ax.plot(ec_d["date"], ec_d["GPP_gC_m2_d"], color=C_EC_RAW,
                lw=0.6, alpha=0.9, zorder=1)
        m_d, m_v = break_gaps(mgpp_d["date"], mgpp_d["GPP_mgpp"],
                              max_gap_days=15)
        ax.plot(m_d, m_v,
                "-o", color=C_MGPP, lw=1.2, ms=3.6, mec="white", mew=0.5,
                alpha=0.95, zorder=2,
                label="Meteosat MGPP  ·  0.05° (~5 km)  ·  10-day")
        ax.plot(ec_d["date"], smooth_break(ec_d["GPP_gC_m2_d"]),
                color=C_EC, lw=1.7, zorder=3,
                label="EC tower GPP  ·  ~100 m  ·  daily (7-d)")
        ax.set_ylabel(r"GPP (gC m$^{-2}$ d$^{-1}$)", fontsize=10.5)
        ymax = max(16, np.nanpercentile(ec_d["GPP_gC_m2_d"], 99) + 2)
        ax.set_ylim(-1, ymax)
        style_time_ax(ax, d0, d1)
        if j == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=9.5)
        ax.text(-0.085, 1.07, panel_labels[1][j], transform=ax.transAxes,
                fontsize=14, fontweight="bold")

        # ---- Row 3: ET (uses independent Meteosat source!) ----
        ax = axes[2, j]
        shade_springs(ax, d0, d1)
        ax.plot(metv3_d["date"], smooth_break(metv3_d["ET_mm"]),
                color=C_METV3, lw=1.4, alpha=0.95, zorder=2,
                label="Meteosat ETv3  ·  0.05° (~5 km)  ·  daily (7-d)")
        ax.plot(ec_d["date"], smooth_break(ec_d["ET_mm"]),
                color=C_EC, lw=1.7, zorder=3,
                label="EC tower ET  ·  ~100 m  ·  daily (7-d)")
        ax.set_ylabel("ET (mm day$^{-1}$)", fontsize=10.5)
        ax.set_ylim(-0.3, 8.5)
        style_time_ax(ax, d0, d1)
        if j == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=9.5)
        ax.text(-0.085, 1.07, panel_labels[2][j], transform=ax.transAxes,
                fontsize=14, fontweight="bold")

    fig.text(0.5, 0.005,
             "Shaded band: Mar–Jun (Mediterranean growing season). "
             "Each legend entry: product · spatial pixel · temporal "
             "cadence. EC = 7-day centered mean. Meteosat ETv3 and MGPP "
             "are loaded from their independent CSVs so satellite lines "
             "remain plotted on EC-missing days.",
             ha="center", fontsize=8.5, style="italic", color="#555555")

    out_png = OUT_DIR / "overview.png"
    out_pdf = OUT_DIR / "overview.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
