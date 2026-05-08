"""First-look figures from ec_daily_master.csv.

Run AFTER unify_ec_daily.py.

Outputs PNGs in ./figs/:
  fig1_et_timeseries.png   ET daily time series, both sites overlaid
  fig2_ndvi_et.png         NDVI vs ET seasonal cycle, by site
  fig3_vpd_et.png          ET vs VPD scatter, colored by season
  fig4_swc_et.png          ET vs surface SWC, colored by VPD
"""
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

CSV = Path("/home/shion-nagamine/bakanposs/ec_daily_master.csv")
FIGS = Path("/home/shion-nagamine/bakanposs/figs")
FIGS.mkdir(exist_ok=True)

SITE_COLORS = {"Oran": "#1f77b4", "TzM": "#d62728"}


def load() -> pd.DataFrame:
    df = pd.read_csv(CSV, parse_dates=["date"])
    df["month"] = df["date"].dt.month
    return df


def fig_et_timeseries(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    for ax, site in zip(axes, ["Oran", "TzM"]):
        s = df[df.site == site]
        ax.plot(s.date, s.ET_mm, ".", ms=2, color=SITE_COLORS[site])
        ax.set_ylabel(f"{site} ET [mm/d]")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("date")
    fig.suptitle("Daily EC ET")
    fig.tight_layout()
    fig.savefig(FIGS / "fig1_et_timeseries.png", dpi=150)
    plt.close(fig)


def fig_ndvi_et(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, site in zip(axes, ["Oran", "TzM"]):
        s = df[df.site == site].dropna(subset=["NDVI", "ET_mm"])
        sc = ax.scatter(s.NDVI, s.ET_mm, c=s.month, cmap="twilight", s=8)
        ax.set_xlabel("NDVI")
        ax.set_title(site)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("ET [mm/d]")
    fig.colorbar(sc, ax=axes, label="month")
    fig.savefig(FIGS / "fig2_ndvi_et.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_vpd_et(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, site in zip(axes, ["Oran", "TzM"]):
        s = df[df.site == site].dropna(subset=["VPD_kPa_mean", "ET_mm"])
        sc = ax.scatter(s.VPD_kPa_mean, s.ET_mm, c=s.month, cmap="twilight", s=8)
        ax.set_xlabel("VPD [kPa]")
        ax.set_title(site)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("ET [mm/d]")
    fig.colorbar(sc, ax=axes, label="month")
    fig.savefig(FIGS / "fig3_vpd_et.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_swc_et(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, site in zip(axes, ["Oran", "TzM"]):
        s = df[df.site == site].dropna(subset=["SWC", "ET_mm", "VPD_kPa_mean"])
        sc = ax.scatter(s.SWC, s.ET_mm, c=s.VPD_kPa_mean, cmap="plasma", s=8,
                        vmin=0, vmax=4)
        ax.set_xlabel("surface SWC [-]")
        ax.set_title(site)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("ET [mm/d]")
    fig.colorbar(sc, ax=axes, label="VPD [kPa]")
    fig.savefig(FIGS / "fig4_swc_et.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    df = load()
    print(df.groupby("site").agg(
        n=("date", "count"),
        ET_mean=("ET_mm", "mean"),
        ET_max=("ET_mm", "max"),
    ))
    fig_et_timeseries(df)
    fig_ndvi_et(df)
    fig_vpd_et(df)
    fig_swc_et(df)
    print(f"figures written to {FIGS}")


if __name__ == "__main__":
    main()
