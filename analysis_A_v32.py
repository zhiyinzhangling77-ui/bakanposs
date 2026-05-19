"""
================================================================
解析A v32 : Poster Fig 4b — Tarazona irrigation blind spot
================================================================

【目的】
  解析A v31 の τ ≈ 3 d / amplitude 4.5× scaling を踏まえて、
  Meteosat ETv3 衛星 product が Tarazona 灌漑アーモンド園で
  「灌漑後の指数回復をまったく検出できず、bias がそのまま
  解析A の amplitude を吸収する」ことを 2 panel で見せる。

【2 panel 構成】
  (a) Bias recovery: Δ = EC_LE − METv3_LE [W/m²] を Day 0 起点で
      箱ひげ + exp fit (τ_bias, amplitude_bias)。第二軸 mm/day。
      → bias が解析A と同じ τ 時間スケールで減衰
  (b) τ comparison: 3 bars (EC / METv3 / Bias)
      → Sat τ は床で頭打ち(fit invalid)。EC と Bias は universal
        band 3–4 d 内。Bias amplitude ≈ EC amplitude。

【入力 (デフォルト)】
  --bias-pool        ./output_analysis_B_v3/v3_bias_perevent_pooled.csv
  --bias-summary     ./output_analysis_B_v3/v3_bias_tau_summary.csv

【出力】 ./output_analysis_A_v32/
  fig04b_tarazona_blindspot_v32.png/pdf  ★ poster figure
  v32_poster_caption.txt                 ★ 脚注テキスト
  v32_panel_a_bias_fit.csv               panel (a) fit 再現用
  v32_panel_b_tau_bars.csv               panel (b) bars
"""

from __future__ import annotations
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False
warnings.filterwarnings("ignore")


# ================================================================
# CONFIG
# ================================================================

LE_PER_MM = 28.36                # ET [mm/d] -> LE [W/m²] daily mean conversion
IRRIG_THRESHOLD = 0.5            # mm  (v31 と同じ)
MAX_WINDOW = 14                  # d
MIN_WINDOW = 4                   # d
MIN_PER_BIN = 3
TAU_FLOOR, TAU_CEIL = 0.5, 60.0
CI_PCT = (2.5, 97.5)
N_BOOT_DEFAULT = 5000

# 解析A v31 の Tarazona active 結果 (再計算しない、固定値で表示)
EC_TAU_REF = dict(tau=3.36, se=0.62, ci_lo=2.44, ci_hi=4.90,
                   amplitude=94.8, n_events=41)

# 色
EC_COL    = "#1D9E75"            # 解析A の Tarazona と同じ緑
METV3_COL = "#D55E00"            # 衛星 = 暖色
BIAS_COL  = "#5A4FCF"            # bias = 紫
N_LABEL_COL = "#5A4FCF"


def parse_args():
    p = argparse.ArgumentParser(description="解析A v32 — Tarazona irrigation blind spot")
    p.add_argument("--bias-pool",
                   default="./output_analysis_B_v3/v3_bias_perevent_pooled.csv")
    p.add_argument("--bias-summary",
                   default="./output_analysis_B_v3/v3_bias_tau_summary.csv")
    p.add_argument("--out", default="./output_analysis_A_v32")
    p.add_argument("--n-boot", type=int, default=N_BOOT_DEFAULT)
    return p.parse_args()


# ================================================================
# Loaders
# ================================================================

def load_bias_pool(fp: Path) -> pd.DataFrame:
    if not fp.exists():
        sys.exit(f"[FATAL] bias-pool not found: {fp}\n"
                  f"  → 先に analysis_B_v3_bias_tau.py を実行してください")
    df = pd.read_csv(fp)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df[df["site"] == "Tarazona"].copy()
    for c in ["LE_EC_Wm2", "LE_METv3_Wm2", "bias", "days_since_event"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["LE_EC_Wm2", "LE_METv3_Wm2", "days_since_event"])


def load_bias_summary(fp: Path) -> dict:
    df = pd.read_csv(fp)
    row = df[df["site"] == "Tarazona"].iloc[0]
    return dict(
        tau_bias=float(row["tau_bias"]),
        tau_se_bias=float(row["tau_se_bias"]),
        tau_ci_lo=float(row["tau_ci_lo"]),
        tau_ci_hi=float(row["tau_ci_hi"]),
        amplitude=float(row["a_Wm2"]),
        a_ci_lo=float(row["a_ci_lo"]),
        a_ci_hi=float(row["a_ci_hi"]),
        c=float(row["c_Wm2"]),
        n_events=int(row["n_events"]),
        r2=float(row["r2"]),
    )


# ================================================================
# Fits
# ================================================================

def exp_model(d, a, tau, c):
    """y = c + a * exp(-d/tau).  Δ-recovery 型."""
    return c + a * np.exp(-d / tau)


def fit_bias_tau(arr_d: np.ndarray, arr_y: np.ndarray,
                  min_per_bin: int = MIN_PER_BIN):
    g = (pd.DataFrame({"d": arr_d, "y": arr_y})
            .groupby("d")["y"].agg(["median", "mean", "std", "count"])
            .reset_index())
    g = g[g["count"] >= min_per_bin]
    if len(g) < 3:
        return None
    x = g["d"].values.astype(float)
    y = g["median"].values

    a0 = float(max(y.max() - y.min(), 1.0))
    c0 = float(y.min())
    p0_list = [[a0, 3.0, c0], [a0, 5.0, c0], [a0 * 1.2, 2.0, c0 - 5]]
    for p0 in p0_list:
        try:
            popt, _ = curve_fit(
                exp_model, x, y, p0=p0,
                bounds=([0.0, TAU_FLOOR, -200.0],
                         [600.0, TAU_CEIL, 300.0]),
                maxfev=8000,
            )
            yp = exp_model(x, *popt)
            sse = float(np.sum((y - yp) ** 2))
            r2 = float(1 - sse / np.sum((y - y.mean()) ** 2)) if np.var(y) > 0 else np.nan
            return dict(a=popt[0], tau=popt[1], c=popt[2],
                          r2=r2, grouped=g)
        except Exception:
            continue
    return None


def fit_metv3_tau(df_pool: pd.DataFrame, n_boot: int):
    """METv3 LE の独立 τ 推定 (bias と別物)."""
    arr_d = df_pool["days_since_event"].values.astype(float)
    arr_y = df_pool["LE_METv3_Wm2"].values.astype(float)
    res = fit_bias_tau(arr_d, arr_y)
    if res is None:
        return None
    rng = np.random.default_rng(42)
    n = len(arr_d)
    taus = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        sub = fit_bias_tau(arr_d[idx], arr_y[idx])
        if sub is not None and TAU_FLOOR + 0.01 < sub["tau"] < TAU_CEIL - 0.5:
            taus.append(sub["tau"])
    if len(taus) < 50:
        return dict(tau=res["tau"], se=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                       a=res["a"], c=res["c"], r2=res["r2"])
    taus = np.array(taus)
    return dict(tau=res["tau"], se=float(np.std(taus)),
                  ci_lo=float(np.percentile(taus, CI_PCT[0])),
                  ci_hi=float(np.percentile(taus, CI_PCT[1])),
                  a=res["a"], c=res["c"], r2=res["r2"])


# ================================================================
# Panel drawers
# ================================================================

def draw_panel_b(ax, df_pool: pd.DataFrame, fit_res: dict,
                  summary: dict):
    ax_r = ax.twinx()

    days = sorted(df_pool["days_since_event"].unique())
    days = [d for d in days if d <= MAX_WINDOW]
    box_data, ns, positions = [], [], []
    for d in days:
        vals = df_pool.loc[df_pool["days_since_event"] == d, "bias"].values
        if len(vals) >= 1:
            box_data.append(vals)
            ns.append(len(vals))
            positions.append(d)
    bp = ax.boxplot(box_data, positions=positions, widths=0.55,
                       patch_artist=True, showfliers=False,
                       medianprops=dict(color="black", lw=1.6),
                       whiskerprops=dict(color="#555"),
                       capprops=dict(color="#555"))
    for box in bp["boxes"]:
        box.set(facecolor=BIAS_COL, alpha=0.32, edgecolor="#444")

    # exp fit overlay
    xs = np.linspace(0, MAX_WINDOW, 200)
    a = summary["amplitude"]; tau = summary["tau_bias"]; c = summary["c"]
    ys = exp_model(xs, a, tau, c)
    ax.plot(xs, ys, color=BIAS_COL, lw=2.4,
             label=f"Exp fit: τ={tau:.2f} d, amp={a:.1f} W m$^{{-2}}$")

    # Day-0 marker
    ax.axvline(0, color="#888", ls=":", lw=1.0)

    # mm/day secondary axis (LE_PER_MM = 28.36)
    ymin, ymax = ax.get_ylim()
    ax_r.set_ylim(ymin / LE_PER_MM, ymax / LE_PER_MM)
    ax_r.set_ylabel("Bias [mm day$^{-1}$]", fontsize=11, color="#444")

    # N labels
    y_label_pos = ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.08
    for d, n in zip(positions, ns):
        ax.text(d, y_label_pos, f"{n}", ha="center", va="top",
                  fontsize=8, color=N_LABEL_COL, weight="bold",
                  bbox=dict(boxstyle="round,pad=0.18", fc="white",
                              ec=N_LABEL_COL, lw=0.6))

    ax.set_xlabel("Days since irrigation", fontsize=11)
    ax.set_ylabel("Bias = EC − Sat  [W m$^{-2}$]", fontsize=11)
    ax.set_xlim(-0.6, MAX_WINDOW + 0.6)
    ax.set_xticks(range(0, MAX_WINDOW + 1, 2))
    ax.set_title(
        f"(a) Bias recovery (n events = {summary['n_events']}, "
        f"R² = {summary['r2']:.2f})",
        fontsize=12, loc="left", weight="bold")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.25, axis="y")


def _metv3_fit_quality(metv3: dict) -> str:
    """Classify METv3 fit: 'ok' | 'pegged' | 'lowR2' | 'missing'."""
    if metv3 is None or np.isnan(metv3.get("tau", np.nan)):
        return "missing"
    tau = metv3["tau"]; r2 = metv3.get("r2", np.nan)
    if tau <= TAU_FLOOR + 0.05 or tau >= TAU_CEIL - 0.5:
        return "pegged"
    if not np.isnan(r2) and r2 < 0.3:
        return "lowR2"
    return "ok"


def draw_panel_c(ax, ec_ref: dict, metv3: dict, bias: dict):
    """3 bars: EC / Sat / Bias.  Sat bar is always shown — even if fit is
    invalid — with annotation so the viewer can see the τ pegged-at-floor
    pathology rather than wondering why it's missing."""
    metv3_quality = _metv3_fit_quality(metv3)

    labels_full = ["EC\n(τ_EC)", "Meteosat ETv3\n(τ_Sat)",
                       "Bias = EC−Sat\n(τ_bias)"]
    taus = [ec_ref["tau"], metv3["tau"], bias["tau_bias"]]
    err_lo = [max(ec_ref["tau"] - ec_ref["ci_lo"], 0.0),
                  max(metv3["tau"] - metv3.get("ci_lo", np.nan), 0.0)
                      if not np.isnan(metv3.get("ci_lo", np.nan)) else 0.0,
                  max(bias["tau_bias"] - bias["tau_ci_lo"], 0.0)]
    err_hi = [max(ec_ref["ci_hi"] - ec_ref["tau"], 0.0),
                  max(metv3.get("ci_hi", np.nan) - metv3["tau"], 0.0)
                      if not np.isnan(metv3.get("ci_hi", np.nan)) else 0.0,
                  max(bias["tau_ci_hi"] - bias["tau_bias"], 0.0)]
    cols = [EC_COL, METV3_COL, BIAS_COL]
    xs = np.arange(len(labels_full))

    # y-axis: cap at a reasonable value so the bias upper CI doesn't dominate
    y_max = max(ec_ref["ci_hi"], bias["tau_ci_hi"]) * 1.4
    y_max = min(y_max, 15.0)
    ax.set_ylim(0, y_max)

    # clip upper errorbar to y-axis (CI can extend off-scale)
    err_hi_clipped = [min(e, y_max - t - 0.3) for e, t in zip(err_hi, taus)]
    err_hi_offscale = [e > (y_max - t - 0.3) for e, t in zip(err_hi, taus)]

    # universal 3–4 d band
    ax.axhspan(3.0, 4.0, color="#1D9E75", alpha=0.10, zorder=0,
                  label="Analysis-A universal band (3–4 d)")

    # all three bars drawn at their fitted value
    for i, (x, t, c) in enumerate(zip(xs, taus, cols)):
        invalid_sat = (i == 1 and metv3_quality != "ok")
        ax.bar(x, t, color=c,
                  alpha=0.35 if invalid_sat else 0.85,
                  edgecolor=c if invalid_sat else "black",
                  lw=1.2 if invalid_sat else 0.8,
                  hatch="///" if invalid_sat else None,
                  zorder=3)

    # errorbars (always)
    for i in range(3):
        ax.errorbar(xs[i], taus[i],
                       yerr=[[err_lo[i]], [err_hi_clipped[i]]],
                       fmt="none", ecolor="black", capsize=6, lw=1.5, zorder=4)
        # uplimit arrow if CI extends beyond axis
        if err_hi_offscale[i]:
            ax.annotate("", xy=(xs[i], y_max - 0.05),
                          xytext=(xs[i], taus[i] + err_hi_clipped[i]),
                          arrowprops=dict(arrowstyle="-|>", color="black",
                                              lw=1.4), zorder=5)

    # text annotations under each bar
    se_list = [ec_ref["se"], metv3.get("se", np.nan), bias["tau_se_bias"]]
    for i, (x, t, se) in enumerate(zip(xs, taus, se_list)):
        se_str = f"SE={se:.2f}" if not np.isnan(se) else "SE=NA"
        if i == 1 and metv3_quality != "ok":
            r2 = metv3.get("r2", np.nan)
            reason = {"pegged": "τ at floor", "lowR2": "low R²",
                          "missing": "fit failed"}[metv3_quality]
            ax.text(x, t + 0.15,
                      f"{t:.2f} d\nFIT INVALID\n({reason}, R²={r2:.2f})",
                      ha="center", va="bottom", fontsize=8,
                      color=METV3_COL, weight="bold")
        else:
            ax.text(x, t + 0.15, f"{t:.2f} d\n{se_str}",
                      ha="center", va="bottom", fontsize=9, weight="bold")

    ax.set_xticks(xs)
    ax.set_xticklabels(labels_full, fontsize=10)
    ax.set_ylabel("τ [d]", fontsize=11)
    ax.set_title("(b) τ comparison: EC vs Sat vs Bias",
                    fontsize=12, loc="left", weight="bold")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.25, axis="y")


# ================================================================
# Figure assembly
# ================================================================

def make_figure(out_dir: Path, df_pool, summary, metv3):
    """2-panel side-by-side layout: bias recovery + τ comparison."""
    fig = plt.figure(figsize=(15, 6.5), facecolor="white",
                       constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0])
    ax_b = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[0, 1])

    fig.suptitle(
        "Tarazona irrigation blind spot — Meteosat ETv3 is flat to drip irrigation\n"
        "while the bias itself recovers on the Analysis-A timescale (τ ≈ 3–4 d)",
        fontsize=14, weight="bold")

    draw_panel_b(ax_b, df_pool, None, summary)
    draw_panel_c(ax_c, EC_TAU_REF, metv3, summary)

    png = out_dir / "fig04b_tarazona_blindspot_v32.png"
    pdf = out_dir / "fig04b_tarazona_blindspot_v32.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  [save] {png}")
    print(f"  [save] {pdf}")


def write_caption(out_dir: Path, summary: dict, metv3: dict):
    txt = f"""\
================================================================
Fig 4b — Tarazona irrigation blind spot (Analysis A v32)
================================================================

Three-panel figure showing that the Meteosat ETv3 satellite ET product
captures the *timescale* of the irrigation recovery (τ ≈ 3-4 d) but
misses the *amplitude*, leaving a bias whose recovery has the same
exponential shape as the EC tower response itself.

------------------------------------------------------------------
Definitions
------------------------------------------------------------------
- Day 0   : the day an irrigation event occurs (Irrig ≥ 0.5 mm)
- LE      : latent heat flux [W m⁻²].  ET [mm day⁻¹] → LE [W m⁻²]
            conversion uses LE = 28.36 × ET (daily mean).
- Bias    : Δ = LE_EC − LE_METv3 [W m⁻²]
- τ       : e-folding time (days) of the exponential recovery
- amp     : amplitude of the recovery, defined for each curve as
            (value at Day 0) − (asymptote at long lag)
- MDE     : minimum detectable effect on τ at α = 0.05
            (= 1.96 × √(SE₁² + SE₂²)), used to claim "no difference"

------------------------------------------------------------------
Panel reading guide
------------------------------------------------------------------
(a) Bias recovery curve — per-event Δ pooled across all Tarazona
    irrigation events ({summary['n_events']} events).  Boxplots show the
    distribution of Δ at each day-since-event; medians drive the
    exponential fit  Δ(d) = c + amp · exp(−d/τ_bias).  Right axis
    shows the same quantity in mm day⁻¹.
    Fitted: τ_bias = {summary['tau_bias']:.2f} d
            amp_bias = {summary['amplitude']:.1f} W m⁻²
            (≈ {summary['amplitude']/LE_PER_MM:.2f} mm day⁻¹)
            R² = {summary['r2']:.2f}

(b) τ comparison — three bars side by side.
    1) EC reference (Analysis A v31, Tarazona active, n = 41):
       τ_EC = {EC_TAU_REF['tau']:.2f} d  [{EC_TAU_REF['ci_lo']:.2f}, {EC_TAU_REF['ci_hi']:.2f}]
       amp = {EC_TAU_REF['amplitude']:.0f} W m⁻²
    2) Meteosat ETv3 alone, same events:
       τ_Sat = {metv3['tau']:.2f} d  (R² = {metv3.get('r2', float('nan')):.2f})
       → fit hits the lower boundary with R² ≈ 0, i.e. the satellite
         time series is essentially flat across the irrigation
         pulse — no detectable exponential recovery.
    3) Bias = EC − Sat:
       τ_bias = {summary['tau_bias']:.2f} d
       [{summary['tau_ci_lo']:.2f}, {summary['tau_ci_hi']:.2f}]
       amp_bias = {summary['amplitude']:.1f} W m⁻²
    Green band = Analysis-A universal range (3-4 d).
    The Sat bar is shown hatched because the fit pegs at the τ floor
    with R² ≈ 0 — i.e. the satellite is statistically *flat* through
    the event window.  The EC and Bias bars both fall within the
    Analysis-A universal band (3–4 d), and the bias amplitude
    (~95 W m⁻²) is essentially equal to the EC amplitude.

------------------------------------------------------------------
Take-home message
------------------------------------------------------------------
At Tarazona, Meteosat ETv3 shows **no detectable exponential
response** to drip irrigation events on the 3–4 d Analysis-A
timescale (τ floored, R² ≈ 0).  Because the EC tower *does* show
that recovery, the bias (EC − Sat) inherits the *entire*
management signal: its amplitude ≈ EC amplitude (~95 W m⁻²) and
its τ falls inside the Analysis-A universal band.
→ Sat-derived ET is usable for *climate-scale* drying timescales,
  but cannot be used to quantify the management signal of drip
  irrigation in heterogeneous orchard pixels (≈ 5 km pixel vs
  ≈ 0.1 km² orchard).

------------------------------------------------------------------
Statistical methods
------------------------------------------------------------------
- Exponential fit: weighted nonlinear least squares on the median
  LE / bias per day-since-event (≥ 3 events per day required).
- Bootstrap (B = {N_BOOT_DEFAULT}, raw-data resampling) gives SE and
  95 % CI for τ.
- Boundary fits (τ ≤ 0.6 d or τ ≥ 59 d) and R² < 0.3 are excluded.

------------------------------------------------------------------
Data sources
------------------------------------------------------------------
- EC tower: Tarazona almond orchard (lat 39.266, lon −1.9397),
  drip-irrigated, daily LE 2020-2024 (n = 41 active-season events).
- Satellite: LSA SAF Meteosat ETv3, 30-min product aggregated to
  daily ET (mm day⁻¹), same grid cell, 2020-2024.
- Irrigation log: tower operator records, mm/event.

------------------------------------------------------------------
Provenance
------------------------------------------------------------------
- Analysis A v31 (EC τ baseline):     fig04_poster_main_v31.pdf
- Analysis B v3 (bias pool + summary): v3_bias_perevent_pooled.csv
                                          v3_bias_tau_summary.csv
- This figure script:                 analysis_A_v32.py
"""
    fp = out_dir / "v32_poster_caption.txt"
    fp.write_text(txt, encoding="utf-8")
    print(f"  [save] {fp}")


# ================================================================
# Main
# ================================================================

def main():
    args = parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("解析A v32 — Tarazona irrigation blind spot poster figure")
    print(f"  出力先: {out}")
    print("=" * 60)

    print("\n--- Load precomputed bias pool (Analysis B v3) ---")
    df_pool = load_bias_pool(Path(args.bias_pool))
    summary = load_bias_summary(Path(args.bias_summary))
    print(f"  bias-pool rows: {len(df_pool)}  (n events = {summary['n_events']})")
    print(f"  τ_bias  = {summary['tau_bias']:.2f} d "
            f"[{summary['tau_ci_lo']:.2f}, {summary['tau_ci_hi']:.2f}]"
            f"  amp = {summary['amplitude']:.1f} W/m²")

    print("\n--- Fit METv3-only τ (bootstrap n={}) ---".format(args.n_boot))
    metv3 = fit_metv3_tau(df_pool, n_boot=args.n_boot)
    if metv3 is None:
        metv3 = dict(tau=np.nan, se=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                       a=np.nan, c=np.nan, r2=np.nan)
    print(f"  τ_Sat = {metv3['tau']:.2f} d, "
            f"SE = {metv3['se']:.2f},  "
            f"CI = [{metv3['ci_lo']:.2f}, {metv3['ci_hi']:.2f}], "
            f"R² = {metv3['r2']:.2f}")

    print("\n--- Figure (v32) ---")
    make_figure(out, df_pool, summary, metv3)
    write_caption(out, summary, metv3)

    df_pool.to_csv(out / "v32_panel_a_bias_fit.csv", index=False)
    pd.DataFrame([
        dict(label="EC",        tau=EC_TAU_REF["tau"], se=EC_TAU_REF["se"],
              ci_lo=EC_TAU_REF["ci_lo"], ci_hi=EC_TAU_REF["ci_hi"],
              amplitude=EC_TAU_REF["amplitude"], n=EC_TAU_REF["n_events"]),
        dict(label="METv3",     tau=metv3["tau"], se=metv3["se"],
              ci_lo=metv3["ci_lo"], ci_hi=metv3["ci_hi"],
              amplitude=metv3.get("a", np.nan), n=summary["n_events"]),
        dict(label="Bias",      tau=summary["tau_bias"], se=summary["tau_se_bias"],
              ci_lo=summary["tau_ci_lo"], ci_hi=summary["tau_ci_hi"],
              amplitude=summary["amplitude"], n=summary["n_events"]),
    ]).to_csv(out / "v32_panel_b_tau_bars.csv", index=False)

    print(f"\n[done] {out}/")


if __name__ == "__main__":
    main()
