"""Refined headline figure: summer-only and growing-season-filtered bias plots.

The all-year bias-by-irrig_bucket plot mixes winter days with summer ones.
Filtering to summer or to growing season (NDVI > 0.3) sharpens the
irrigation decoupling signal because winter d8+ days also have near-zero
EC ET, which artificially shrinks the bias.

Reads: master_full.csv
Writes:
  figs/fig_C2_bias_by_irrig_seasonal.png    4 panels: all / summer / growing
                                            for both MOD16 and PML
  fig_C_summer_summary.csv                  bucket statistics
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CSV = Path("/home/shion-nagamine/bakanposs/master_full.csv")
FIGS = Path("/home/shion-nagamine/bakanposs/figs")
FIGS.mkdir(exist_ok=True)
OUT_CSV = Path("/home/shion-nagamine/bakanposs/fig_C_summer_summary.csv")

PRODUCTS = [("MOD16_ET", "MOD16", "tab:purple"),
            ("PML_ET",   "PML",   "tab:green")]
BUCKETS = ["irrig_d0-3", "irrig_d4-7", "irrig_d8plus"]


def bucket_stats(df: pd.DataFrame, sat_col: str) -> pd.DataFrame:
    bias = df[sat_col] - df["ET_mm"]
    rows = []
    for bk in BUCKETS:
        sub = bias[df.irrig_bucket == bk].dropna()
        if len(sub) < 5:
            continue
        rows.append({
            "irrig_bucket": bk,
            "n": len(sub),
            "mean_bias": float(sub.mean()),
            "median_bias": float(sub.median()),
            "std_bias": float(sub.std()),
            "ci_lo": float(sub.quantile(0.025)),
            "ci_hi": float(sub.quantile(0.975)),
        })
    return pd.DataFrame(rows)


def panel(ax, df: pd.DataFrame, sat_col: str, label: str, c: str,
          title: str) -> None:
    bias = df[sat_col] - df["ET_mm"]
    data = []
    n_each = []
    for bk in BUCKETS:
        v = bias[df.irrig_bucket == bk].dropna().to_numpy()
        data.append(v)
        n_each.append(len(v))
    bp = ax.boxplot(data, tick_labels=BUCKETS, showfliers=False,
                    patch_artist=True)
    for box in bp["boxes"]:
        box.set(facecolor=c, alpha=0.4)
    ax.axhline(0, color="k", lw=0.7, ls="--")
    for i, v in enumerate(data):
        ax.scatter(np.full_like(v, i + 1, dtype=float), v,
                   s=4, alpha=0.3, color=c)
    ymax = ax.get_ylim()[1]
    for i, n in enumerate(n_each):
        ax.text(i + 1, ymax * 0.95, f"n={n}", ha="center", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("bias [mm/d]   sat - EC")
    ax.set_xlabel("days since last irrigation")
    ax.grid(alpha=0.3)


def main() -> None:
    df = pd.read_csv(CSV, parse_dates=["date"])
    t = df[df.site == "TzM"].copy()
    t["irrig_bucket"] = pd.Categorical(
        t["irrig_bucket"], categories=BUCKETS + ["no_irrig"], ordered=True)

    masks = {
        "all year":     pd.Series(True, index=t.index),
        "summer only":  t.season == "summer",
        "growing (NDVI>0.3)": t.NDVI.fillna(0) > 0.3,
        "summer × NDVI>0.3":  (t.season == "summer") & (t.NDVI.fillna(0) > 0.3),
    }

    # write a single CSV summary for paper table
    rows = []
    for filter_name, m in masks.items():
        sub = t[m]
        for col, label, _ in PRODUCTS:
            stats = bucket_stats(sub, col)
            stats["filter"] = filter_name
            stats["product"] = label
            rows.append(stats)
    summary = pd.concat(rows, ignore_index=True)
    summary = summary[["filter", "product", "irrig_bucket",
                       "n", "mean_bias", "median_bias",
                       "std_bias", "ci_lo", "ci_hi"]]
    summary.to_csv(OUT_CSV, index=False)
    pd.set_option("display.float_format", "{:+.2f}".format)
    print(summary.to_string(index=False))
    print(f"\nwrote {OUT_CSV}")

    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharey="row")
    for r, (col, label, c) in enumerate(PRODUCTS):
        for j, (fname, m) in enumerate(masks.items()):
            panel(axes[r, j], t[m], col, label, c,
                  f"{label} | {fname}")
    fig.suptitle("TzM (almond): satellite ET bias by irrigation bucket and "
                 "season filter", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_C2_bias_by_irrig_seasonal.png", dpi=150)
    plt.close(fig)
    print(f"wrote {FIGS / 'fig_C2_bias_by_irrig_seasonal.png'}")


if __name__ == "__main__":
    main()
