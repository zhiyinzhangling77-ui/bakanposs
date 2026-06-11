"""Killer figures for the EC vs satellite ET paper.

Reads master_full.csv (EC + satellite merged) and produces:

  fig_A_scatter_ET.png       sat_ET vs EC_ET, site-coloured, 1:1 line, RMSE/r
  fig_B_timeseries_ET.png    EC vs MOD16 vs PML time series, by site
  fig_C_bias_vs_irrig.png    (sat_ET - EC_ET) vs days_since_irrig boxplots
                              the headline figure for the irrigation
                              decoupling hypothesis
  fig_D_seasonal_bias.png    monthly mean bias by site
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO = Path(__file__).parent.parent
CSV = REPO / "master_full.csv"
if not CSV.exists():
    CSV = REPO / "master_full_v2.csv"
FIGS = REPO / "figs"
FIGS.mkdir(exist_ok=True)

SITE_COLOR = {"Oran": "#1f77b4", "TzM": "#d62728"}
PRODUCTS = [("MOD16_ET", "MOD16", "tab:purple"),
            ("PML_ET",   "PML",   "tab:green")]


def metrics(x: pd.Series, y: pd.Series) -> tuple[int, float, float, float, float]:
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 5:
        return n, np.nan, np.nan, np.nan, np.nan
    x, y = x[m], y[m]
    r = float(np.corrcoef(x, y)[0, 1])
    rmse = float(np.sqrt(np.mean((y - x) ** 2)))
    mbe = float(np.mean(y - x))
    slope = float(np.polyfit(x, y, 1)[0])
    return n, r, rmse, mbe, slope


def fig_A_scatter(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    for ax, (col, label, _c) in zip(axes, PRODUCTS):
        if col not in df:
            continue
        for site, g in df.groupby("site"):
            ax.scatter(g["ET_mm"], g[col], s=10, alpha=0.5,
                       color=SITE_COLOR[site], label=site)
            n, r, rmse, mbe, slope = metrics(g["ET_mm"], g[col])
            ax.text(0.04, 0.95 - 0.07 * (0 if site == "Oran" else 1),
                    f"{site}: n={n} r={r:+.2f} RMSE={rmse:.2f} MBE={mbe:+.2f}",
                    transform=ax.transAxes, fontsize=9,
                    color=SITE_COLOR[site])
        lim = [0, max(df["ET_mm"].max(), df[col].max()) * 1.05]
        ax.plot(lim, lim, "k--", lw=0.8)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("EC ET [mm/d]")
        ax.set_ylabel(f"{label} ET [mm/d]")
        ax.set_title(label)
        ax.grid(alpha=0.3)
        ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_A_scatter_ET.png", dpi=150)
    plt.close(fig)


def fig_B_timeseries(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(13, 6), sharex=True)
    for ax, site in zip(axes, ["Oran", "TzM"]):
        s = df[df.site == site].sort_values("date")
        ax.plot(s["date"], s["ET_mm"], "k.", ms=2, label="EC", alpha=0.6)
        for col, label, c in PRODUCTS:
            if col in s:
                ax.plot(s["date"], s[col], "-", lw=1, color=c,
                        alpha=0.7, label=label)
        ax.set_ylabel(f"{site} ET [mm/d]")
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("date")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_B_timeseries_ET.png", dpi=150)
    plt.close(fig)


def fig_C_bias_vs_irrig(df: pd.DataFrame) -> None:
    """The headline figure: (sat - EC) by irrigation bucket, TzM only."""
    t = df[df.site == "TzM"].copy()
    bucket_order = ["irrig_d0-3", "irrig_d4-7", "irrig_d8plus", "no_irrig"]
    t["irrig_bucket"] = pd.Categorical(t["irrig_bucket"],
                                       categories=bucket_order, ordered=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, (col, label, c) in zip(axes, PRODUCTS):
        if col not in t:
            continue
        bias = t[col] - t["ET_mm"]
        data = []
        for bk in bucket_order:
            v = bias[t.irrig_bucket == bk].dropna()
            data.append(v.to_numpy())
        bp = ax.boxplot(data, labels=bucket_order, showfliers=False,
                        patch_artist=True)
        for box in bp["boxes"]:
            box.set(facecolor=c, alpha=0.4)
        ax.axhline(0, color="k", lw=0.7, ls="--")
        for i, v in enumerate(data):
            ax.scatter(np.full_like(v, i + 1, dtype=float), v,
                       s=4, alpha=0.3, color=c)
            ax.text(i + 1, ax.get_ylim()[1] * 0.95, f"n={len(v)}",
                    ha="center", fontsize=8)
        ax.set_title(f"{label} - EC ET")
        ax.set_ylabel("bias [mm/d]   sat - EC")
        ax.set_xlabel("days since last irrigation")
        ax.grid(alpha=0.3)
    fig.suptitle("TzM (almond): satellite ET underestimates EC right after irrigation",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_C_bias_vs_irrig.png", dpi=150)
    plt.close(fig)


def fig_D_seasonal_bias(df: pd.DataFrame) -> None:
    df = df.copy()
    df["month"] = df["date"].dt.month

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, (col, label, c) in zip(axes, PRODUCTS):
        if col not in df:
            continue
        df[f"bias_{label}"] = df[col] - df["ET_mm"]
        for site, g in df.groupby("site"):
            mb = g.groupby("month")[f"bias_{label}"].agg(["mean", "count"])
            mb = mb[mb["count"] >= 5]
            ax.plot(mb.index, mb["mean"], "o-", color=SITE_COLOR[site],
                    label=site)
        ax.axhline(0, color="k", lw=0.7, ls="--")
        ax.set_title(f"{label} - EC ET")
        ax.set_xlabel("month")
        ax.set_ylabel("monthly mean bias [mm/d]")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_D_seasonal_bias.png", dpi=150)
    plt.close(fig)


def main() -> None:
    df = pd.read_csv(CSV, parse_dates=["date"])
    print(f"loaded {len(df):,} rows")
    print("ET overlap counts (EC + each satellite product):")
    for col, label, _ in PRODUCTS:
        if col in df:
            n = (df["ET_mm"].notna() & df[col].notna()).sum()
            print(f"  EC × {label:6s}: {n}")

    fig_A_scatter(df)
    fig_B_timeseries(df)
    fig_C_bias_vs_irrig(df)
    fig_D_seasonal_bias(df)
    print(f"\nfigures written to {FIGS}")


if __name__ == "__main__":
    main()
