"""Poster overview figure: NDVI, GPP, ET dynamics at two flux towers.

Layout (3 rows x 2 cols):
  Row 1: NDVI (EC tower vs Sentinel-2 S2_NDVI)
  Row 2: GPP  (EC tower vs Meteosat MGPP 10-day)
  Row 3: ET   (EC tower vs Meteosat ETv3 daily)

All rows shade Mar-Jun (peak Mediterranean growing season) per year.

Spatial / temporal scales annotated in each legend entry:
  EC tower         : ~100 m footprint    , daily aggregate (30-min raw)
  Sentinel-2 NDVI  : 10 m pixel          , ~5-day revisit
  Meteosat MGPP    : 0.05° (~5 km)       , 10-day composite
  Meteosat ETv3    : 0.05° (~5 km)       , daily aggregate (30-min raw)

Inputs:  data/master_full_v2.csv, data/mgpp_decadal_all.csv
Outputs: figures/poster/overview.png  (+ .pdf)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

REPO = Path(__file__).parent.parent
CSV = REPO / "data" / "master_full_v2.csv"
MGPP_CSV = REPO / "data" / "mgpp_decadal_all.csv"
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
C_S2      = "#1e6fb3"   # Sentinel-2 NDVI
SHADE     = "#f6e7b8"
GRID      = "#dddddd"


def smooth(s: pd.Series, w: int = 7) -> pd.Series:
    return s.rolling(w, center=True, min_periods=max(3, w // 2)).mean()


def break_gaps(dates, values, max_gap_days: int = 15):
    """Return (dates, values) with NaN inserted where consecutive points
    are more than ``max_gap_days`` apart, so matplotlib breaks the line
    rather than drawing a long bridge across missing data.
    Input arrays are first sorted by date and dropped of any NaN values.
    """
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


def reindex_daily(sub: pd.DataFrame) -> pd.DataFrame:
    """Reindex to continuous daily dates so missing days become explicit
    NaN rows; smoothed lines then naturally break across long gaps."""
    full = pd.date_range(sub["date"].min(), sub["date"].max(), freq="D")
    return sub.set_index("date").reindex(full).rename_axis("date").reset_index()


def smooth_break(s: pd.Series, gap_days: int = 7, w: int = 7) -> pd.Series:
    """Centered 7-day mean, then force NaN inside any data gap of more
    than ``gap_days`` consecutive missing input values."""
    out = s.rolling(w, center=True, min_periods=max(3, w // 2)).mean()
    miss = s.isna()
    streak = miss.groupby((~miss).cumsum()).cumsum()
    # mask out smoothed points that fall inside a gap > gap_days
    long_gap = streak > gap_days
    out = out.where(~long_gap)
    return out


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


def main() -> None:
    df = pd.read_csv(CSV, parse_dates=["date"])
    mgpp = pd.read_csv(MGPP_CSV, parse_dates=["date"])
    mgpp = mgpp[(mgpp["GPP_mgpp"] > 0) & (mgpp["GPP_mgpp"] < 30)].copy()

    fig, axes = plt.subplots(
        3, 2, figsize=(13.5, 9.0),
        gridspec_kw={"hspace": 0.32, "wspace": 0.18},
    )

    panel_labels = [["a", "b"], ["c", "d"], ["e", "f"]]

    for j, site in enumerate(SITES):
        sub = df[df["site"] == site].sort_values("date").copy()
        sub = reindex_daily(sub)
        msub = mgpp[mgpp["site"] == site].sort_values("date").copy()
        d0, d1 = sub["date"].min(), sub["date"].max()
        msub = msub[(msub["date"] >= d0) & (msub["date"] <= d1)]

        # ---- Row 1: NDVI (EC vs Sentinel-2) ----
        ax = axes[0, j]
        shade_springs(ax, d0, d1)
        s2_d, s2_v = break_gaps(sub["date"], sub["S2_NDVI"], max_gap_days=15)
        ax.plot(s2_d, s2_v,
                "-o", color=C_S2, lw=1.0, ms=3.2, mec="white", mew=0.4,
                alpha=0.9, zorder=2,
                label="Sentinel-2 NDVI  ·  10 m  ·  ~5-day")
        ax.plot(sub["date"], smooth_break(sub["NDVI"]),
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

        # ---- Row 2: GPP (EC vs Meteosat MGPP) ----
        ax = axes[1, j]
        shade_springs(ax, d0, d1)
        ax.plot(sub["date"], sub["GPP_gC_m2_d"], color=C_EC_RAW,
                lw=0.6, alpha=0.9, zorder=1)
        m_d, m_v = break_gaps(msub["date"], msub["GPP_mgpp"], max_gap_days=15)
        ax.plot(m_d, m_v,
                "-o", color=C_MGPP, lw=1.2, ms=3.6, mec="white", mew=0.5,
                alpha=0.95, zorder=2,
                label="Meteosat MGPP  ·  0.05° (~5 km)  ·  10-day")
        ax.plot(sub["date"], smooth_break(sub["GPP_gC_m2_d"]),
                color=C_EC, lw=1.7, zorder=3,
                label="EC tower GPP  ·  ~100 m  ·  daily (7-d)")
        ax.set_ylabel(r"GPP (gC m$^{-2}$ d$^{-1}$)", fontsize=10.5)
        ymax = max(16, np.nanpercentile(sub["GPP_gC_m2_d"], 99) + 2)
        ax.set_ylim(-1, ymax)
        style_time_ax(ax, d0, d1)
        if j == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=9.5)
        ax.text(-0.085, 1.07, panel_labels[1][j], transform=ax.transAxes,
                fontsize=14, fontweight="bold")

        # ---- Row 3: ET (EC vs Meteosat) ----
        ax = axes[2, j]
        shade_springs(ax, d0, d1)
        ax.plot(sub["date"], smooth_break(sub["ET_metv3_mm"]),
                color=C_METV3, lw=1.4, alpha=0.95, zorder=2,
                label="Meteosat ETv3  ·  0.05° (~5 km)  ·  daily (7-d)")
        ax.plot(sub["date"], smooth_break(sub["ET_mm"]),
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
             "Each legend entry annotates product name · spatial pixel · "
             "temporal cadence. EC values are smoothed with a 7-day "
             "centered mean; lines break across data gaps >7 d (daily) "
             "or >15 d (sparse).",
             ha="center", fontsize=8.5, style="italic", color="#555555")

    out_png = OUT_DIR / "overview.png"
    out_pdf = OUT_DIR / "overview.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
