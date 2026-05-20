"""
================================================================
analysis_B_v6_driver_attribution.py
  STEP 4 (Fig 4) for poster: driver attribution at drip-irrigated
  Tarazona.  Tests whether EC ET is demand-driven and Meteosat ETv3
  is supply-limited, using a standardised linear regression.
================================================================

Hypothesis
----------
At Tarazona (1-ha drip-irrigated almond orchard):

  EC tower ET    → demand-driven
                   Irrigation + 1–5 m taproot keep water non-limiting,
                   so the orchard transpires as fast as atmospheric
                   demand allows (Rn, VPD).

  Meteosat ETv3  → supply-limited
                   The SVAT models the 5-km pixel as rainfed
                   Mediterranean.  Its f(SM) availability function
                   scales ET with the coarse satellite SM input.
                   Drip irrigation is invisible at this scale, so
                   ETv3 returns the rainfed surroundings' ET.

If the hypothesis is correct, the standardised regression
coefficients should split cleanly:

  β_demand (EC)  ≫  β_supply (EC)
  β_demand (Sat) ≪  β_supply (Sat)

Period selection: Jun–Sep (Tarazona active months)
--------------------------------------------------
- Almond active transpiration : May (leaf-out) → Oct (harvest).
- Peak irrigation             : Jun–Sep (≥ 12 mm/event, 7–47 events/mo).
- Peak EC–Sat divergence      : Jul–Aug.
- Analysis A v31 active window: TARA_ACTIVE_MONTHS = [6, 7, 8, 9].

We use **Jun–Sep** (4 months) rather than meteorological JJA (3 months)
because September still carries substantial transpiration and irrigation,
and matching Analysis A's window keeps the τ baseline (3.36 d, amp 95
W m⁻²) directly comparable.

Predictors
----------
Atmospheric demand:
    VPD_kPa_mean   mean daily vapour-pressure deficit (kPa)
    Rn_Wm2         net radiation (W m⁻²)
Water supply:
    SWC            tower in-situ SWC, ~3–7 cm depth (shallow)
    smap_rootzone  SMAP L4 root-zone SM (9 km)

Note on satellite SM
--------------------
Meteosat ETv3 ingests H-SAF products (H141 12.5 km surface, H142
12.5 km SWI, H26 1 km downscaled SWI).  At time of writing those are
**not yet fetched**.  We use SMAP L4 root-zone SM as a coarse-satellite
SM **proxy** — same physical limitation (coarse pixel, surface-derived
microwave SM, blind to drip irrigation).  The driver pattern is
expected to be insensitive to which coarse SM product is used; this
script will be re-run with H142/H26 substituting smap_rootzone once
those data are available.

Output
------
output_analysis_B_v6/
  v6_driver_attribution_coefficients.csv   raw β + SE + t + p per
                                           predictor for both targets
  v6_driver_attribution_data.csv           weekly means used in the fit
  fig01_driver_attribution_bars.png/pdf    standardised β bar chart
                                           (poster Fig 4)
"""

from __future__ import annotations
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False
warnings.filterwarnings("ignore")


# ================================================================
# CONFIG
# ================================================================

SITE_LABEL_IN_MASTER = "TzM"             # column value in master_full_v2
ACTIVE_MONTHS = [6, 7, 8, 9]             # Jun–Sep (Tarazona active period)

# Predictors and their pretty (English) display names
PREDICTORS = [
    "VPD_kPa_mean",
    "Rn_Wm2",
    "SWC",
    "smap_rootzone",
]
PRETTY = {
    "VPD_kPa_mean":  "VPD",
    "Rn_Wm2":        "Net radiation",
    "SWC":           "Tower SWC\n(in-situ, 3–7 cm)",
    "smap_rootzone": "SMAP root-zone\n(satellite, 9 km)",
}
GROUP = {
    "VPD_kPa_mean":  "demand",
    "Rn_Wm2":        "demand",
    "SWC":           "supply",
    "smap_rootzone": "supply",
}

# Targets
TARGET_EC  = "ET_mm"
TARGET_SAT = "ET_metv3_mm"

# Colors
EC_COL       = "#1D9E75"   # EC tower (green, consistent with v31/v32)
SAT_COL      = "#D55E00"   # Meteosat ETv3 (orange)
DEMAND_COL   = "#3878C6"   # group label colour (blue)
SUPPLY_COL   = "#9C27B0"   # group label colour (purple)


def parse_args():
    p = argparse.ArgumentParser(description="Driver attribution v6")
    p.add_argument("--master", default="data/master_full_v2.csv")
    p.add_argument("--out", default="output_analysis_B_v6")
    p.add_argument("--site", default=SITE_LABEL_IN_MASTER)
    return p.parse_args()


# ================================================================
# Data
# ================================================================

def load_weekly(master_csv: Path, site_label: str) -> pd.DataFrame:
    df = pd.read_csv(master_csv)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    sub = df[df["site"] == site_label].copy()
    sub = sub[sub["date"].dt.month.isin(ACTIVE_MONTHS)]

    needed = [TARGET_EC, TARGET_SAT] + PREDICTORS
    missing = [c for c in needed if c not in sub.columns]
    if missing:
        raise SystemExit(f"[FATAL] master file missing columns: {missing}")

    weekly = (sub.set_index("date")[needed]
                 .resample("W")
                 .mean()
                 .dropna())
    return weekly


# ================================================================
# Regression
# ================================================================

def standardised_regression(y: pd.Series, X: pd.DataFrame) -> dict:
    """Z-score y and each column of X, then fit OLS with an intercept.
    Returns the standardised coefficients, SEs, t-values, p-values
    and R² (the intercept becomes ~0 by construction)."""
    y_z = StandardScaler().fit_transform(y.values.reshape(-1, 1)).ravel()
    X_z = StandardScaler().fit_transform(X.values)
    X_z = sm.add_constant(X_z)
    fit = sm.OLS(y_z, X_z).fit()
    return dict(
        predictors=list(X.columns),
        coef=fit.params[1:].tolist(),
        se=fit.bse[1:].tolist(),
        tvalues=fit.tvalues[1:].tolist(),
        pvalues=fit.pvalues[1:].tolist(),
        r2=float(fit.rsquared),
        n=int(len(y_z)),
    )


def sig_label(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# ================================================================
# Figure
# ================================================================

def draw_bar_comparison(out_dir: Path, res_ec: dict, res_sat: dict,
                          period_str: str):
    preds = res_ec["predictors"]
    n = len(preds)
    x = np.arange(n)
    width = 0.38

    fig, ax = plt.subplots(figsize=(11, 6.5), facecolor="white",
                             constrained_layout=True)

    # bars
    ax.bar(x - width/2, res_ec["coef"], width,
              yerr=res_ec["se"],
              color=EC_COL, edgecolor="black", lw=0.6,
              label=f"EC tower  (R² = {res_ec['r2']:.2f},  "
                     f"n = {res_ec['n']} weeks)",
              capsize=4, error_kw=dict(ecolor="black", lw=1.0))
    ax.bar(x + width/2, res_sat["coef"], width,
              yerr=res_sat["se"],
              color=SAT_COL, edgecolor="black", lw=0.6,
              label=f"Meteosat ETv3  (R² = {res_sat['r2']:.2f},  "
                     f"n = {res_sat['n']} weeks)",
              capsize=4, error_kw=dict(ecolor="black", lw=1.0))

    # significance stars — placed above (or below) the error-bar cap,
    # not on top of the bar itself, so they never overlap the whisker.
    all_coefs = np.array(res_ec["coef"] + res_sat["coef"])
    all_ses   = np.array(res_ec["se"]   + res_sat["se"])
    span = float(np.max(np.abs(all_coefs) + all_ses))
    pad  = span * 0.08         # gap above the error-bar cap

    def _star_y(c: float, se: float) -> float:
        """Star y position: just above the upper cap for c >= 0,
        just below the lower cap for c < 0."""
        return (c + se + pad) if c >= 0 else (c - se - pad)

    for i, (c, se, p) in enumerate(zip(res_ec["coef"],
                                              res_ec["se"],
                                              res_ec["pvalues"])):
        ax.text(x[i] - width/2, _star_y(c, se), sig_label(p),
                  ha="center", va="bottom" if c >= 0 else "top",
                  fontsize=13, weight="bold")
    for i, (c, se, p) in enumerate(zip(res_sat["coef"],
                                              res_sat["se"],
                                              res_sat["pvalues"])):
        ax.text(x[i] + width/2, _star_y(c, se), sig_label(p),
                  ha="center", va="bottom" if c >= 0 else "top",
                  fontsize=13, weight="bold")

    # zero line
    ax.axhline(0, color="black", lw=0.8)

    # group separator
    n_demand = sum(1 for p_ in preds if GROUP[p_] == "demand")
    ax.axvline(n_demand - 0.5, color="#AAAAAA", lw=1.0, ls=":", zorder=0)

    # y-axis — enlarged so star labels above the upper CI cap clear the
    # group-label banner and don't crowd the bars.
    y_max =  span * 1.70
    y_min = -span * 0.95
    ax.set_ylim(y_min, y_max)

    # group annotations at top
    y_label = y_max * 0.94
    ax.text((n_demand - 1) / 2, y_label, "ATMOSPHERIC DEMAND",
              ha="center", va="top", fontsize=13, weight="bold",
              color=DEMAND_COL,
              bbox=dict(boxstyle="round,pad=0.45", fc="white",
                          ec=DEMAND_COL, lw=1.4))
    ax.text((n_demand + n - 1) / 2, y_label, "WATER SUPPLY",
              ha="center", va="top", fontsize=13, weight="bold",
              color=SUPPLY_COL,
              bbox=dict(boxstyle="round,pad=0.45", fc="white",
                          ec=SUPPLY_COL, lw=1.4))

    # x-axis
    ax.set_xticks(x)
    ax.set_xticklabels([PRETTY[p] for p in preds], fontsize=12)
    ax.tick_params(axis="y", labelsize=11)

    ax.set_ylabel("Standardised regression coefficient  (β)", fontsize=14)
    ax.set_title(
        "Driver attribution:  EC tower vs Meteosat ETv3 at drip-irrigated Tarazona\n"
        + period_str,
        fontsize=14, weight="bold")

    # Significance legend (lower-left, larger so it is poster-readable)
    ax.text(0.01, 0.02,
              "*** p < 0.001    ** p < 0.01    * p < 0.05    "
              "ns: not significant",
              transform=ax.transAxes, fontsize=12, color="#333",
              va="bottom", ha="left",
              bbox=dict(boxstyle="round,pad=0.3", fc="white",
                          ec="#BBBBBB", lw=0.8, alpha=0.9))

    # Main colour-coded legend (EC vs Sat) — larger
    ax.legend(loc="lower right", fontsize=13, framealpha=0.95,
                borderpad=0.7, labelspacing=0.6)
    ax.grid(alpha=0.25, axis="y")

    png = out_dir / "fig01_driver_attribution_bars.png"
    pdf = out_dir / "fig01_driver_attribution_bars.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  [save] {png}")
    print(f"  [save] {pdf}")


# ================================================================
# Main
# ================================================================

def main():
    args = parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    print("=" * 64)
    print("解析B v6 — Driver attribution (Tarazona, Jun–Sep, weekly)")
    print("=" * 64)

    weekly = load_weekly(Path(args.master), args.site)
    date_min = weekly.index.min().date()
    date_max = weekly.index.max().date()
    n = len(weekly)
    period_str = (f"(weekly means, Jun–Sep {date_min.year}–{date_max.year}, "
                     f"n = {n} weeks)")
    print(f"  weekly rows after Jun–Sep + non-null filter: n = {n}")
    print(f"  period: {date_min} → {date_max}")

    y_ec  = weekly[TARGET_EC]
    y_sat = weekly[TARGET_SAT]
    X     = weekly[PREDICTORS]

    print("\n--- Standardised regression: EC tower ET ---")
    res_ec = standardised_regression(y_ec, X)
    print(f"  R² = {res_ec['r2']:.3f}")
    for p_name, c, se, p in zip(res_ec["predictors"], res_ec["coef"],
                                       res_ec["se"], res_ec["pvalues"]):
        print(f"  {p_name:18s}  β = {c:+.3f}  SE = {se:.3f}  "
                f"p = {p:.3g}  {sig_label(p)}")

    print("\n--- Standardised regression: Meteosat ETv3 ---")
    res_sat = standardised_regression(y_sat, X)
    print(f"  R² = {res_sat['r2']:.3f}")
    for p_name, c, se, p in zip(res_sat["predictors"], res_sat["coef"],
                                       res_sat["se"], res_sat["pvalues"]):
        print(f"  {p_name:18s}  β = {c:+.3f}  SE = {se:.3f}  "
                f"p = {p:.3g}  {sig_label(p)}")

    # Save tables
    rows = []
    for tgt_name, res in [("EC_ET", res_ec), ("Meteosat_ETv3", res_sat)]:
        for p_name, c, se, t, pval in zip(res["predictors"], res["coef"],
                                                  res["se"], res["tvalues"],
                                                  res["pvalues"]):
            rows.append(dict(
                target=tgt_name, predictor=p_name,
                group=GROUP[p_name],
                coef=c, se=se, tvalue=t, pvalue=pval,
                significance=sig_label(pval), r2=res["r2"], n=res["n"]
            ))
    pd.DataFrame(rows).to_csv(
        out / "v6_driver_attribution_coefficients.csv", index=False)
    print(f"\n  [save] {out}/v6_driver_attribution_coefficients.csv")

    weekly.to_csv(out / "v6_driver_attribution_data.csv")
    print(f"  [save] {out}/v6_driver_attribution_data.csv")

    print("\n--- Figure (v6) ---")
    draw_bar_comparison(out, res_ec, res_sat, period_str)

    print(f"\n[done] {out}/")


if __name__ == "__main__":
    main()
