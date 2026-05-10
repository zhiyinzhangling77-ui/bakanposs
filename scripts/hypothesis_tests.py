"""Test the testable hypotheses (H1, H4, H6) on master_full_v2.csv.

Hypotheses tested here:
  H1  tau-based bias correction reduces TzM summer RMSE substantially.
      For each product, predict bias(t) = a*exp(-t/tau) + c with the
      tau_fit parameters; subtract from raw satellite ET; compare RMSE
      and MBE on TzM summer x NDVI>0.3.

  H4  Adding days_since_irrig to a VPD-only bias model improves AIC.
      Model 1: bias ~ VPD (atmospheric-demand-only baseline)
      Model 2: bias ~ days_since_irrig (irrigation-cycle-only)
      Model 3: bias ~ VPD + days_since_irrig + VPD x days_since_irrig
      Lower AIC = better. AIC penalises extra parameters so an
      improvement in fit must outweigh a 2-per-parameter penalty.

  H6  SMAP L4 root-zone correlates with in-situ SWC at Oran (rainfed)
      and produces a comparable SDS metric.
      Computed at Oran summer x NDVI>0.3 only (no irrigation noise).

Untestable here (need external data):
  H2  drip vs flood vs sprinkler -> need other FLUXNET sites
  H3  crop root depth -> need olive / vegetable / vine sites
  H5  METv3 5km pixel mixing -> need SIGPAC parcels or sub-pixel landcover
  H7  SDS as paper index across Mediterranean -> EuroFLUX inventory
  H8  regional ET correction validation -> regional ET maps

Reads:  master_full_v2.csv, tau_fit_summary.csv
Writes: hypothesis_tests_summary.csv
        figs/fig_H1_correction.png
        figs/fig_H4_aic.png
        figs/fig_H6_smap_oran.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

REPO = Path(__file__).parent.parent
CSV = REPO / "master_full_v2.csv"
TAU = REPO / "tau_fit_summary.csv"
FIGS = REPO / "figs"
FIGS.mkdir(exist_ok=True)
OUT = REPO / "hypothesis_tests_summary.csv"

PRODUCTS = [("MOD16_ET",    "MOD16", "tab:purple"),
            ("PML_ET",      "PML",   "tab:green"),
            ("ET_metv3_mm", "METv3", "tab:orange")]


def model_full(t, a, tau, c):
    return a * np.exp(-t / tau) + c


# ── H1 ───────────────────────────────────────────────────────────────────────

def h1_correction(df: pd.DataFrame, tau_fit: pd.DataFrame) -> list[dict]:
    """Apply tau correction, compare pre/post RMSE on TzM summer."""
    tzm = df[(df.site == "TzM") & (df.season == "summer") &
             (df.NDVI.fillna(0) > 0.3)].copy()
    print(f"  TzM summer x NDVI>0.3: n={len(tzm)}")

    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, (col, label, c) in zip(axes, PRODUCTS):
        if col not in tzm.columns:
            continue
        sub = tau_fit[(tau_fit["product"] == label) &
                       (tau_fit["model"] == "full")]
        if len(sub) == 0:
            continue
        a   = float(sub["a_med"].iloc[0])
        tau = float(sub["tau_med"].iloc[0])
        cc  = float(sub["c_med"].iloc[0])

        d = tzm.dropna(subset=[col, "ET_mm", "days_since_irrig"]).copy()
        t = d["days_since_irrig"].clip(upper=20).to_numpy()
        d["pred_bias"] = model_full(t, a, tau, cc)
        d["sat_corr"] = d[col] - d["pred_bias"]

        bias_raw  = (d[col] - d["ET_mm"]).to_numpy()
        bias_corr = (d["sat_corr"] - d["ET_mm"]).to_numpy()
        rmse_raw  = float(np.sqrt(np.mean(bias_raw ** 2)))
        rmse_corr = float(np.sqrt(np.mean(bias_corr ** 2)))
        rows.append({
            "test": "H1", "product": label, "n": len(d),
            "RMSE_raw": rmse_raw, "RMSE_corr": rmse_corr,
            "RMSE_reduction_pct": (rmse_raw - rmse_corr) / rmse_raw * 100,
            "MBE_raw":  float(bias_raw.mean()),
            "MBE_corr": float(bias_corr.mean()),
            "tau_used": tau, "a_used": a, "c_used": cc,
        })
        print(f"  {label}: RMSE {rmse_raw:.2f} -> {rmse_corr:.2f}  "
              f"({(rmse_raw-rmse_corr)/rmse_raw*100:+.1f}%)  "
              f"MBE {bias_raw.mean():+.2f} -> {bias_corr.mean():+.2f}")

        ax.scatter(d["ET_mm"], d[col], s=10, alpha=0.4, color=c,
                    label=f"raw  RMSE={rmse_raw:.2f}")
        ax.scatter(d["ET_mm"], d["sat_corr"], s=10, alpha=0.5, color="k",
                    label=f"corr RMSE={rmse_corr:.2f}")
        m = max(d["ET_mm"].max(), d[col].max(), 10)
        ax.plot([0, m], [0, m], "k-", lw=0.7)
        ax.set_xlabel("EC ET [mm/d]")
        ax.set_ylabel(f"{label} ET [mm/d]")
        ax.set_title(f"{label}: tau correction at TzM summer")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

    fig.suptitle("H1: exponential decay correction reduces TzM summer error",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_H1_correction.png", dpi=150)
    plt.close(fig)
    return rows


# ── H4 ───────────────────────────────────────────────────────────────────────

def fit_ols(X: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Return (AIC, SSE) for OLS y ~ X."""
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    sse = float(np.sum((y - yhat) ** 2))
    sigma2 = sse / n
    if sigma2 <= 0:
        return np.inf, sse
    ll = -n / 2 * np.log(2 * np.pi * sigma2) - sse / (2 * sigma2)
    aic = 2 * k - 2 * ll
    return float(aic), sse


def h4_aic(df: pd.DataFrame) -> list[dict]:
    """Compare three bias models with AIC."""
    tzm = df[(df.site == "TzM") & (df.season == "summer") &
             (df.NDVI.fillna(0) > 0.3)].copy()
    tzm["d"] = tzm["days_since_irrig"].clip(upper=20)

    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, (col, label, c) in zip(axes, PRODUCTS):
        if col not in tzm.columns:
            continue
        d = tzm.dropna(subset=[col, "ET_mm", "VPD_kPa_mean", "d"]).copy()
        d["bias"] = d[col] - d["ET_mm"]
        if len(d) < 30:
            continue
        n = len(d)
        y = d["bias"].to_numpy()

        X1 = np.column_stack([np.ones(n), d["VPD_kPa_mean"]])
        X2 = np.column_stack([np.ones(n), d["d"]])
        X3 = np.column_stack([np.ones(n), d["VPD_kPa_mean"], d["d"],
                               d["VPD_kPa_mean"] * d["d"]])
        aic_v, sse_v   = fit_ols(X1, y)
        aic_d, sse_d   = fit_ols(X2, y)
        aic_vd, sse_vd = fit_ols(X3, y)

        best = ["VPD only", "days only", "VPD+d+interaction"][
            int(np.argmin([aic_v, aic_d, aic_vd]))
        ]
        delta = min(aic_v, aic_d, aic_vd) - aic_v   # vs VPD baseline
        rows.append({
            "test": "H4", "product": label, "n": n,
            "AIC_VPD": aic_v, "AIC_d": aic_d, "AIC_VPDd_int": aic_vd,
            "best_model": best, "deltaAIC_best_vs_VPD": delta,
            "SSE_VPD": sse_v, "SSE_d": sse_d, "SSE_VPDd_int": sse_vd,
        })
        print(f"  {label} (n={n}): AIC VPD={aic_v:.0f}  d={aic_d:.0f}  "
              f"VPD+d={aic_vd:.0f}  -> best={best}")

        labels = ["VPD\nonly", "days_irrig\nonly", "VPD+d+\ninteraction"]
        bars = ax.bar(labels, [aic_v, aic_d, aic_vd],
                       color=c, alpha=0.7, edgecolor="k")
        for bar, a in zip(bars, [aic_v, aic_d, aic_vd]):
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height(), f"{a:.0f}",
                    ha="center", va="bottom", fontsize=10)
        ax.set_ylabel("AIC (lower is better)")
        ax.set_title(f"{label}: bias model AIC")
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle("H4: days_since_irrig captures bias structure better than VPD alone",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_H4_aic.png", dpi=150)
    plt.close(fig)
    return rows


# ── H6 ───────────────────────────────────────────────────────────────────────

def compute_sds(arr: np.ndarray, sm: np.ndarray) -> tuple[float, int, int]:
    p25, p75 = np.nanquantile(sm, [0.25, 0.75])
    dry = sm < p25
    norm = (sm >= p25) & (sm <= p75)
    if dry.sum() < 5 or norm.sum() < 5:
        return np.nan, int(dry.sum()), int(norm.sum())
    sds = 1.0 - np.nanmean(arr[dry]) / np.nanmean(arr[norm])
    return float(sds), int(dry.sum()), int(norm.sum())


def h6_smap_multistratum(df: pd.DataFrame) -> list[dict]:
    """SMAP root-zone vs in-situ SWC across multiple strata.

    The original test on Oran summer (n=34, post-harvest) was limited.
    Here we evaluate the SMAP-as-substitute hypothesis on every stratum
    where SDS is meaningful: Oran spring (active growth, n>>200 expected),
    Oran summer, TzM spring/summer/fall, and TzM summer by irrig bucket.
    """
    NDVI_GATE = 0.3
    strata = []

    # site x season
    for site in ["Oran", "TzM"]:
        for season in ["spring", "summer", "fall"]:
            sub = df[(df.site == site) & (df.season == season) &
                     (df.NDVI.fillna(0) > NDVI_GATE)]
            if len(sub) >= 20:
                strata.append((f"{site}_{season}", site, sub))

    # TzM summer x irrig bucket
    tzm_sum = df[(df.site == "TzM") & (df.season == "summer") &
                 (df.NDVI.fillna(0) > NDVI_GATE)]
    for bk in ["irrig_d0-3", "irrig_d4-7", "irrig_d8plus"]:
        sub = tzm_sum[tzm_sum.irrig_bucket == bk]
        if len(sub) >= 20:
            strata.append((f"TzM_summer_{bk}", "TzM", sub))

    rows = []
    print(f"  evaluating {len(strata)} strata")
    for label, site, sub in strata:
        # 1. correlation
        both = sub.dropna(subset=["SWC", "smap_rootzone"])
        r, p = (np.nan, np.nan)
        if len(both) >= 10:
            r, p = stats.pearsonr(both["SWC"], both["smap_rootzone"])

        # 2. SDS for each SM source.  Oran lacks LE_Wm2, so use ET_mm.
        et_col = "LE_Wm2" if (sub["LE_Wm2"].notna().sum() >= 30) else "ET_mm"
        sds_results = {}
        for sm_col, key in [("SWC", "in_situ"),
                            ("smap_rootzone", "smap_rz"),
                            ("smap_surface", "smap_surf")]:
            d = sub.dropna(subset=[et_col, sm_col])
            if len(d) < 20:
                sds_results[key] = (np.nan, len(d))
                continue
            sds, _, _ = compute_sds(d[et_col].to_numpy(), d[sm_col].to_numpy())
            sds_results[key] = (sds, len(d))

        rows.append({
            "test": "H6", "stratum": label, "site": site,
            "n_total": len(sub),
            "et_col_used": et_col,
            "r_swc_smap_rz": r, "p_swc_smap_rz": p, "n_corr": len(both),
            "SDS_in_situ":  sds_results["in_situ"][0],
            "SDS_smap_rz":  sds_results["smap_rz"][0],
            "SDS_smap_surf": sds_results["smap_surf"][0],
            "n_in_situ":  sds_results["in_situ"][1],
            "n_smap_rz":  sds_results["smap_rz"][1],
            "n_smap_surf": sds_results["smap_surf"][1],
        })
        print(f"  {label:<28s}  n={len(sub):>3d}  "
              f"r={r:+.2f}  "
              f"SDS_insitu={sds_results['in_situ'][0]:+.2f}  "
              f"SDS_smap={sds_results['smap_rz'][0]:+.2f}")

    # ── figure ──────────────────────────────────────────────────────────────
    summary = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    # Panel 1: SDS comparison (in-situ vs SMAP rz) per stratum
    ax = axes[0]
    s = summary.dropna(subset=["SDS_in_situ", "SDS_smap_rz"])
    for site, marker in [("Oran", "o"), ("TzM", "s")]:
        ss = s[s.site == site]
        ax.scatter(ss["SDS_in_situ"], ss["SDS_smap_rz"],
                   s=ss["n_total"] * 0.5, alpha=0.7,
                   marker=marker, label=site,
                   edgecolor="k", linewidth=0.5)
        for _, r_ in ss.iterrows():
            ax.annotate(r_["stratum"].replace("_", "\n"),
                        xy=(r_["SDS_in_situ"], r_["SDS_smap_rz"]),
                        xytext=(5, 5), textcoords="offset points",
                        fontsize=7)
    lims = [-0.6, 0.6]
    ax.plot(lims, lims, "k--", lw=0.7)
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("SDS using in-situ SWC (5 cm)")
    ax.set_ylabel("SDS using SMAP root-zone")
    ax.set_title("SDS substitution: 1:1 line = perfect proxy")
    ax.legend(); ax.grid(alpha=0.3)

    # Panel 2: correlation r per stratum
    ax = axes[1]
    s2 = summary.sort_values("r_swc_smap_rz", ascending=False)
    colors = ["tab:blue" if x == "Oran" else "tab:red" for x in s2.site]
    bars = ax.barh(range(len(s2)), s2["r_swc_smap_rz"],
                    color=colors, alpha=0.8, edgecolor="k")
    ax.set_yticks(range(len(s2)))
    ax.set_yticklabels(s2["stratum"], fontsize=8)
    for i, (v, n) in enumerate(zip(s2["r_swc_smap_rz"], s2["n_corr"])):
        if np.isfinite(v):
            ax.text(v, i, f"  r={v:+.2f}  n={n}",
                    va="center", fontsize=8)
    ax.axvline(0.5, color="green", lw=0.7, ls="--",
               label="r=0.5 acceptable")
    ax.set_xlabel("Pearson r(in-situ SWC, SMAP root-zone)")
    ax.set_title("Cross-source agreement per stratum")
    ax.legend(); ax.grid(axis="x", alpha=0.3)
    ax.set_xlim(-0.2, 1.05)

    # Panel 3: focal Oran spring scatter
    ax = axes[2]
    oran_spring = df[(df.site == "Oran") & (df.season == "spring") &
                      (df.NDVI.fillna(0) > NDVI_GATE)]
    both = oran_spring.dropna(subset=["SWC", "smap_rootzone"])
    if len(both) >= 10:
        r_os, p_os = stats.pearsonr(both["SWC"], both["smap_rootzone"])
        ax.scatter(both["SWC"], both["smap_rootzone"],
                   s=15, alpha=0.5, color="tab:blue")
        m, b = np.polyfit(both["SWC"], both["smap_rootzone"], 1)
        xx = np.linspace(both["SWC"].min(), both["SWC"].max(), 50)
        ax.plot(xx, m * xx + b, "k--", lw=1)
        ax.set_title(f"Oran spring (focal): r={r_os:+.2f} "
                     f"p={p_os:.1e} n={len(both)}")
        ax.set_xlabel("in-situ SWC (5 cm) [vol %]")
        ax.set_ylabel("SMAP root-zone [m^3/m^3]")
        ax.grid(alpha=0.3)

    fig.suptitle("H6 multi-stratum: SMAP root-zone vs in-situ SWC across "
                 "all meaningful strata", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_H6_smap_multistratum.png", dpi=150)
    plt.close(fig)
    return rows


# Backward-compat wrapper so existing main() call still works
def h6_smap_oran(df: pd.DataFrame) -> list[dict]:
    return h6_smap_multistratum(df)


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    df  = pd.read_csv(CSV, parse_dates=["date"])
    tau = pd.read_csv(TAU)

    rows = []
    print("=" * 60)
    print("H1: tau-based bias correction effectiveness")
    print("=" * 60)
    rows.extend(h1_correction(df, tau))

    print("\n" + "=" * 60)
    print("H4: AIC of bias models (VPD vs days_since_irrig)")
    print("=" * 60)
    rows.extend(h4_aic(df))

    print("\n" + "=" * 60)
    print("H6: SMAP root-zone as drought predictor at Oran")
    print("=" * 60)
    rows.extend(h6_smap_oran(df))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    pd.set_option("display.float_format", "{:+.3f}".format)
    print("\n" + "=" * 60)
    print(out.to_string(index=False))
    print(f"\nwrote {OUT}")
    print(f"figures in {FIGS}/fig_H1_*, fig_H4_*, fig_H6_*")


if __name__ == "__main__":
    main()
