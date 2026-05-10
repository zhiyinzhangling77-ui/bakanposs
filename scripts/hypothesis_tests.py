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


def h6_smap_oran(df: pd.DataFrame) -> list[dict]:
    """SMAP root-zone vs in-situ SWC at rainfed Oran summer."""
    oran = df[(df.site == "Oran") & (df.season == "summer") &
              (df.NDVI.fillna(0) > 0.3)].copy()
    print(f"  Oran summer x NDVI>0.3: n={len(oran)}")

    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # ── correlation between SWC and SMAP rootzone ──
    both = oran.dropna(subset=["SWC", "smap_rootzone"])
    if len(both) >= 10:
        r, p = stats.pearsonr(both["SWC"], both["smap_rootzone"])
    else:
        r, p = np.nan, np.nan

    ax = axes[0]
    ax.scatter(both["SWC"], both["smap_rootzone"], s=15, alpha=0.5,
               color="tab:blue")
    if len(both) >= 2:
        m, b = np.polyfit(both["SWC"], both["smap_rootzone"], 1)
        xx = np.linspace(both["SWC"].min(), both["SWC"].max(), 50)
        ax.plot(xx, m*xx + b, "k--", lw=1)
    ax.set_xlabel("in-situ SWC (5 cm)  [vol %]")
    ax.set_ylabel("SMAP root-zone (~1 m)  [m^3/m^3]")
    ax.set_title(f"Oran summer: r={r:.2f}  p={p:.1e}  n={len(both)}")
    ax.grid(alpha=0.3)

    # ── SDS using each soil-moisture source ──
    d_swc = oran.dropna(subset=["ET_mm", "SWC"])
    sds_swc, n_dry_s, n_norm_s = compute_sds(
        d_swc["ET_mm"].to_numpy(), d_swc["SWC"].to_numpy())
    d_smap = oran.dropna(subset=["ET_mm", "smap_rootzone"])
    sds_smap, n_dry_m, n_norm_m = compute_sds(
        d_smap["ET_mm"].to_numpy(), d_smap["smap_rootzone"].to_numpy())
    d_sms = oran.dropna(subset=["ET_mm", "smap_surface"])
    sds_sms, n_dry_ms, n_norm_ms = compute_sds(
        d_sms["ET_mm"].to_numpy(), d_sms["smap_surface"].to_numpy())

    rows.append({
        "test": "H6", "stratum": "Oran summer NDVI>0.3",
        "r_swc_smap_rz": r, "p_swc_smap_rz": p,
        "SDS_in_situ_SWC": sds_swc,
        "SDS_smap_rootzone": sds_smap,
        "SDS_smap_surface": sds_sms,
        "n_swc": len(d_swc), "n_smap_rz": len(d_smap),
        "n_smap_surf": len(d_sms),
    })
    print(f"  SDS_in-situ      = {sds_swc:+.3f}  (n={len(d_swc)})")
    print(f"  SDS_SMAP rootzone= {sds_smap:+.3f}  (n={len(d_smap)})")
    print(f"  SDS_SMAP surface = {sds_sms:+.3f}  (n={len(d_sms)})")
    print(f"  r(SWC, SMAP_rz)  = {r:+.3f}  p={p:.2e}")

    ax = axes[1]
    bars = ax.bar(["in-situ\nSWC", "SMAP\nroot-zone", "SMAP\nsurface"],
                   [sds_swc, sds_smap, sds_sms],
                   color=["tab:blue", "tab:orange", "tab:red"],
                   alpha=0.8)
    for bar, v in zip(bars, [sds_swc, sds_smap, sds_sms]):
        if np.isfinite(v):
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height(), f"{v:+.2f}",
                    ha="center", va="bottom", fontsize=11)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_ylabel("SDS")
    ax.set_title("Oran summer SDS by soil moisture source")
    ax.grid(axis="y", alpha=0.3)

    # ── time series of SWC and SMAP rootzone ──
    ax = axes[2]
    o = oran.sort_values("date")
    if "SWC" in o:
        ax2 = ax.twinx()
        ax.plot(o["date"], o["SWC"], color="tab:blue",
                lw=1, label="SWC (5 cm)")
        ax.set_ylabel("SWC (5 cm) [vol %]", color="tab:blue")
        ax2.plot(o["date"], o["smap_rootzone"], color="tab:orange",
                 lw=1, label="SMAP rootzone")
        ax2.set_ylabel("SMAP rootzone [m^3/m^3]", color="tab:orange")
    ax.set_xlabel("date")
    ax.set_title("Oran: SWC vs SMAP rootzone time series")
    ax.grid(alpha=0.3)

    fig.suptitle("H6: SMAP root-zone as alternative drought predictor at Oran",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_H6_smap_oran.png", dpi=150)
    plt.close(fig)
    return rows


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
