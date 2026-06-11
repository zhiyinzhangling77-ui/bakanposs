"""Fit an exponential decay model for satellite ET bias vs days_since_irrig.

Two model structures are tested for each satellite product:

  Model A (transient + permanent):
    bias(t) = a * exp(-t / tau) + c

  Model B (transient only, c = 0):
    bias(t) = a * exp(-t / tau)

The MOD16 result already showed a plateau at ~-2.7 mm/day, so model A
should fit it better; PML decayed toward 0, so model B should suffice.
We fit per-day-bin medians (TzM summer, NDVI > 0.3) and report tau and a
with bootstrap CIs.

Reads: master_full.csv
Writes:
  tau_fit_summary.csv
  figs/fig_F_tau_fit.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

REPO = Path(__file__).parent.parent
# Use master_full_v2 if it exists (includes METv3); otherwise fall back
CSV = REPO / "master_full_v2.csv"
if not CSV.exists():
    CSV = REPO / "master_full.csv"
OUT_CSV = REPO / "tau_fit_summary.csv"
FIGS = REPO / "figs"
FIGS.mkdir(exist_ok=True)

MIN_N = 8
N_BOOT = 500
RNG = np.random.default_rng(42)


def model_full(t, a, tau, c):
    return a * np.exp(-t / tau) + c


def model_transient(t, a, tau):
    return a * np.exp(-t / tau)


def fit_one(t: np.ndarray, y: np.ndarray, w: np.ndarray, fn, p0):
    try:
        popt, pcov = curve_fit(fn, t, y, p0=p0, sigma=1.0 / np.sqrt(w),
                                maxfev=5000)
        return popt
    except Exception:
        return None


def bootstrap_params(t: np.ndarray, y: np.ndarray, w: np.ndarray,
                     fn, p0, n: int = N_BOOT) -> dict:
    samples = []
    n_pts = len(t)
    for _ in range(n):
        idx = RNG.choice(n_pts, size=n_pts, replace=True)
        if len(np.unique(idx)) < len(p0) + 1:
            continue
        popt = fit_one(t[idx], y[idx], w[idx], fn, p0)
        if popt is None:
            continue
        samples.append(popt)
    if len(samples) < 50:
        return {}
    arr = np.array(samples)
    out = {}
    for i, name in enumerate(["a", "tau", "c"][:arr.shape[1]]):
        out[f"{name}_med"] = float(np.median(arr[:, i]))
        out[f"{name}_lo"] = float(np.quantile(arr[:, i], 0.025))
        out[f"{name}_hi"] = float(np.quantile(arr[:, i], 0.975))
    return out


def prepare(df: pd.DataFrame, sat_col: str) -> pd.DataFrame:
    """Per-day-since-irrigation bin: median bias and weight (sqrt n)."""
    bias = df[sat_col] - df["ET_mm"]
    sub = pd.DataFrame({"d": df["days_since_irrig"], "bias": bias})
    sub = sub.dropna()
    sub = sub[sub["d"] >= 0]
    sub["d_bin"] = sub["d"].astype(int).clip(upper=20)
    bins = sub.groupby("d_bin")["bias"].agg(["count", "median", "mean"])
    bins = bins[bins["count"] >= MIN_N].reset_index()
    return bins


def fit_raw(t: np.ndarray, y: np.ndarray, fn, p0):
    """Fit on raw daily data (more stable than binned medians)."""
    try:
        popt, pcov = curve_fit(fn, t, y, p0=p0, maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except Exception:
        return None, None


def bootstrap_raw(t: np.ndarray, y: np.ndarray, fn, p0,
                  n: int = N_BOOT) -> dict:
    samples = []
    n_pts = len(t)
    for _ in range(n):
        idx = RNG.choice(n_pts, size=n_pts, replace=True)
        popt, _ = fit_raw(t[idx], y[idx], fn, p0)
        if popt is None or np.any(np.abs(popt) > 1e3):
            continue
        if popt[1] <= 0 or popt[1] > 60:
            continue
        samples.append(popt)
    if len(samples) < 50:
        return {}
    arr = np.array(samples)
    out = {}
    for i, name in enumerate(["a", "tau", "c"][:arr.shape[1]]):
        out[f"{name}_med"] = float(np.median(arr[:, i]))
        out[f"{name}_lo"] = float(np.quantile(arr[:, i], 0.025))
        out[f"{name}_hi"] = float(np.quantile(arr[:, i], 0.975))
    return out


def plot_one(ax, bins: pd.DataFrame, label: str, color: str,
             popt_full, popt_trans):
    ax.scatter(bins["d_bin"], bins["median"],
               s=bins["count"] * 0.5, color=color, alpha=0.7, edgecolors="k",
               label=f"observed median (size ∝ n)")
    t_grid = np.linspace(0, 20, 200)
    if popt_full is not None:
        a, tau, c = popt_full
        ax.plot(t_grid, model_full(t_grid, a, tau, c),
                color="k", lw=2, ls="-",
                label=f"full: a={a:+.2f}, τ={tau:.1f} d, c={c:+.2f}")
    if popt_trans is not None:
        a, tau = popt_trans
        ax.plot(t_grid, model_transient(t_grid, a, tau),
                color="k", lw=1.5, ls="--",
                label=f"transient: a={a:+.2f}, τ={tau:.1f} d")
    ax.axhline(0, color="gray", lw=0.7)
    ax.set_xlabel("days since last irrigation")
    ax.set_ylabel(f"({label} − EC) ET   median  [mm/d]")
    ax.set_title(label)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="lower right")


def main() -> None:
    df = pd.read_csv(CSV, parse_dates=["date"])
    t = df[(df.site == "TzM")
           & (df.season == "summer")
           & (df["NDVI"].fillna(0) > 0.3)].copy()
    print(f"TzM summer × NDVI>0.3: n={len(t)}")

    rows = []

    products = [("MOD16_ET", "MOD16", "tab:purple"),
                ("PML_ET",   "PML",   "tab:green")]
    if "ET_metv3_mm" in t.columns and t["ET_metv3_mm"].notna().sum() > 30:
        products.append(("ET_metv3_mm", "METv3", "tab:orange"))

    n_prod = len(products)
    fig, axes = plt.subplots(1, n_prod, figsize=(6 * n_prod, 5), sharex=True)
    if n_prod == 1:
        axes = [axes]

    for ax, (col, label, color) in zip(axes, products):
        if col not in t:
            continue
        bias = (t[col] - t["ET_mm"])
        d = t["days_since_irrig"]
        ok = d.notna() & bias.notna() & (d >= 0)
        ts_raw = d[ok].astype(float).clip(upper=20).to_numpy()
        ys_raw = bias[ok].astype(float).to_numpy()

        bins = prepare(t, col)
        if len(bins) < 4:
            print(f"  {label}: not enough day bins")
            continue

        popt_full, perr_full = fit_raw(ts_raw, ys_raw, model_full,
                                         p0=[-2.0, 3.0, -1.0])
        popt_trans, perr_trans = fit_raw(ts_raw, ys_raw, model_transient,
                                           p0=[-2.0, 3.0])
        boot_full = bootstrap_raw(ts_raw, ys_raw, model_full,
                                    p0=[-2.0, 3.0, -1.0])
        boot_trans = bootstrap_raw(ts_raw, ys_raw, model_transient,
                                     p0=[-2.0, 3.0])

        rows.append({"product": label, "model": "full",
                     "n_obs": len(ts_raw),
                     **(boot_full if boot_full else {})})
        rows.append({"product": label, "model": "transient",
                     "n_obs": len(ts_raw),
                     **(boot_trans if boot_trans else {})})

        plot_one(ax, bins, label, color, popt_full, popt_trans)
        if popt_full is not None:
            sd = perr_full if perr_full is not None else [np.nan] * 3
            print(f"  {label} full:      a={popt_full[0]:+.2f}±{sd[0]:.2f}  "
                  f"τ={popt_full[1]:.2f}±{sd[1]:.2f}  "
                  f"c={popt_full[2]:+.2f}±{sd[2]:.2f}")
        if popt_trans is not None:
            sd = perr_trans if perr_trans is not None else [np.nan] * 2
            print(f"  {label} transient: a={popt_trans[0]:+.2f}±{sd[0]:.2f}  "
                  f"τ={popt_trans[1]:.2f}±{sd[1]:.2f}")

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_CSV, index=False)
    print(f"\nwrote {OUT_CSV}")
    pd.set_option("display.float_format", "{:+.3f}".format)
    print(summary.to_string(index=False))

    fig.suptitle("Satellite ET bias decay after irrigation at TzM "
                 "(summer × NDVI>0.3)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_F_tau_fit.png", dpi=150)
    plt.close(fig)
    print(f"wrote {FIGS / 'fig_F_tau_fit.png'}")


if __name__ == "__main__":
    main()
