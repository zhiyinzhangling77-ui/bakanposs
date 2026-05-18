"""Poster overview figure: NDVI, GPP, ET dynamics at two flux towers.

Layout (3 rows x 2 cols):
  Row 1: NDVI (EC tower vs Sentinel-2 S2_NDVI)
  Row 2: GPP  (EC tower vs Meteosat MGPP 10-day)
  Row 3: ET   (EC tower vs Meteosat ETv3 daily)

All rows shade Mar-Jun (peak Mediterranean growing season) per year.

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
        msub = mgpp[mgpp["site"] == site].sort_values("date").copy()
        d0, d1 = sub["date"].min(), sub["date"].max()
        msub = msub[(msub["date"] >= d0) & (msub["date"] <= d1)]

        # ---- Row 1: NDVI (EC vs Sentinel-2) ----
        ax = axes[0, j]
        shade_springs(ax, d0, d1)
        s2 = sub.dropna(subset=["S2_NDVI"])
        ax.plot(s2["date"], s2["S2_NDVI"],
                "-o", color=C_S2, lw=1.0, ms=3.2, mec="white", mew=0.4,
                alpha=0.9, zorder=2, label="Sentinel-2 NDVI")
        ax.plot(sub["date"], smooth(sub["NDVI"]),
                color=C_EC, lw=1.7, zorder=3, label="EC tower NDVI (7-d)")
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
        ax.plot(msub["date"], msub["GPP_mgpp"],
                "-o", color=C_MGPP, lw=1.2, ms=3.6, mec="white", mew=0.5,
                alpha=0.95, zorder=2, label="Meteosat MGPP (10-d)")
        ax.plot(sub["date"], smooth(sub["GPP_gC_m2_d"]),
                color=C_EC, lw=1.7, zorder=3, label="EC tower (7-d)")
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
        ax.plot(sub["date"], smooth(sub["ET_metv3_mm"]),
                color=C_METV3, lw=1.4, alpha=0.95, zorder=2,
                label="Meteosat ETv3 (7-d)")
        ax.plot(sub["date"], smooth(sub["ET_mm"]),
                color=C_EC, lw=1.7, zorder=3, label="EC tower (7-d)")
        ax.set_ylabel("ET (mm day$^{-1}$)", fontsize=10.5)
        ax.set_ylim(-0.3, 8.5)
        style_time_ax(ax, d0, d1)
        if j == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=9.5)
        ax.text(-0.085, 1.07, panel_labels[2][j], transform=ax.transAxes,
                fontsize=14, fontweight="bold")

    fig.text(0.5, 0.005,
             "Shaded band: Mar–Jun (typical Mediterranean growing season). "
             "Sentinel-2 NDVI (10 m) ; Meteosat MGPP = LSA SAF 10-day GPP "
             "composite (0.05°) ; Meteosat ETv3 = LSA SAF DMET v3 (0.05°).",
             ha="center", fontsize=8.8, style="italic", color="#555555")

    out_png = OUT_DIR / "overview.png"
    out_pdf = OUT_DIR / "overview.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
