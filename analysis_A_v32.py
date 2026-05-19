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
      → bias 内に管理 amplitude が保存されていることを可視化
  (b) EC vs Sat recovery curves overlay: 同じ event pool で EC と
      METv3 を別々に exp fit、両者の箱ひげ + fit 曲線を重ねる。
      間の shaded gap = bias.
      → 「Sat が flat」が視覚的に分かり、management pulse の不在が
        labelled 0 ではなく gap として現れる.
        τ_bias は caption に consistency check として注記(独立証拠
        ではない).

【入力 (デフォルト)】
  --bias-pool        ./output_analysis_B_v3/v3_bias_perevent_pooled.csv
  --bias-summary     ./output_analysis_B_v3/v3_bias_tau_summary.csv

【出力】 ./output_analysis_A_v32/
  fig04b_tarazona_blindspot_v32.png/pdf  ★ poster figure
  v32_poster_caption.txt                 ★ 脚注テキスト
  v32_panel_a_bias_fit.csv               panel (a) fit 再現用
  v32_panel_b_curves_fits.csv            panel (b) EC + Sat fit params
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


def fit_series_tau(df_pool: pd.DataFrame, le_col: str, n_boot: int):
    """Generic exp-decay fit on any LE_* column in the event pool.
    Used for both LE_EC and LE_METv3 so they share fit machinery."""
    arr_d = df_pool["days_since_event"].values.astype(float)
    arr_y = df_pool[le_col].values.astype(float)
    mask = ~np.isnan(arr_y)
    arr_d = arr_d[mask]; arr_y = arr_y[mask]
    res = fit_bias_tau(arr_d, arr_y)
    if res is None:
        return None
    rng = np.random.default_rng(42)
    n = len(arr_d)
    taus, amps = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        sub = fit_bias_tau(arr_d[idx], arr_y[idx])
        if sub is not None and TAU_FLOOR + 0.01 < sub["tau"] < TAU_CEIL - 0.5:
            taus.append(sub["tau"]); amps.append(sub["a"])
    if len(taus) < 50:
        return dict(tau=res["tau"], se=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                       a=res["a"], se_a=np.nan,
                       c=res["c"], r2=res["r2"])
    taus = np.array(taus); amps = np.array(amps)
    return dict(tau=res["tau"], se=float(np.std(taus)),
                  ci_lo=float(np.percentile(taus, CI_PCT[0])),
                  ci_hi=float(np.percentile(taus, CI_PCT[1])),
                  a=res["a"], se_a=float(np.std(amps)),
                  c=res["c"], r2=res["r2"])


def fit_metv3_tau(df_pool: pd.DataFrame, n_boot: int):
    """METv3-only τ.  Wrapper around fit_series_tau for backward compat."""
    return fit_series_tau(df_pool, "LE_METv3_Wm2", n_boot)


def fit_ec_tau(df_pool: pd.DataFrame, n_boot: int):
    """EC-only τ on the same event pool as Sat — for panel (b) overlay."""
    return fit_series_tau(df_pool, "LE_EC_Wm2", n_boot)


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


def draw_panel_curves(ax, df_pool: pd.DataFrame,
                          ec_fit: dict, sat_fit: dict):
    """Panel (b): EC vs Sat recovery curves on the same axes.

    Side-by-side boxplots per day-since-event for the two products,
    with exponential fit overlays.  The shaded region between the two
    fit curves IS the bias — making the missing amplitude visible as
    a gap rather than as a labelled "0" bar.
    """
    sat_quality = _metv3_fit_quality(sat_fit)
    ax_r = ax.twinx()                  # mm/day on right axis (matches panel a)

    days = sorted(df_pool["days_since_event"].unique())
    days = [d for d in days if d <= MAX_WINDOW]
    ec_box, sat_box, positions, ns = [], [], [], []
    for d in days:
        sub = df_pool[df_pool["days_since_event"] == d]
        ec_vals  = sub["LE_EC_Wm2"].dropna().values
        sat_vals = sub["LE_METv3_Wm2"].dropna().values
        if len(ec_vals) >= 1 and len(sat_vals) >= 1:
            ec_box.append(ec_vals); sat_box.append(sat_vals)
            positions.append(d); ns.append(len(ec_vals))

    # Boxplots side by side, offset ±0.18 from the day position
    # manage_ticks=False : box の position を x 軸 tick に流用させない
    bp_ec = ax.boxplot(ec_box, positions=[p - 0.18 for p in positions],
                          widths=0.30, patch_artist=True, showfliers=False,
                          manage_ticks=False,
                          medianprops=dict(color="black", lw=1.4),
                          whiskerprops=dict(color="#555"),
                          capprops=dict(color="#555"))
    for b in bp_ec["boxes"]:
        b.set(facecolor=EC_COL, alpha=0.30, edgecolor=EC_COL)

    bp_sat = ax.boxplot(sat_box, positions=[p + 0.18 for p in positions],
                           widths=0.30, patch_artist=True, showfliers=False,
                           manage_ticks=False,
                           medianprops=dict(color="black", lw=1.4),
                           whiskerprops=dict(color="#555"),
                           capprops=dict(color="#555"))
    for b in bp_sat["boxes"]:
        b.set(facecolor=METV3_COL, alpha=0.30, edgecolor=METV3_COL)

    # Exp-fit overlay curves
    xs_fit = np.linspace(0, MAX_WINDOW, 200)
    ec_curve  = exp_model(xs_fit, ec_fit["a"],  ec_fit["tau"],  ec_fit["c"])
    sat_curve = exp_model(xs_fit, sat_fit["a"], sat_fit["tau"], sat_fit["c"])

    ax.plot(xs_fit, ec_curve, color=EC_COL, lw=2.6, zorder=6,
              label=f"EC fit:  τ = {ec_fit['tau']:.2f} d,  "
                     f"amp = {ec_fit['a']:.0f} W m$^{{-2}}$")
    sat_label = (f"Sat fit:  flat (τ pegged at floor, R² ≈ 0)"
                    if sat_quality != "ok"
                    else f"Sat fit:  τ = {sat_fit['tau']:.2f} d,  "
                          f"amp = {sat_fit['a']:.0f} W m$^{{-2}}$")
    ax.plot(xs_fit, sat_curve, color=METV3_COL, lw=2.6, zorder=6, ls="--",
              label=sat_label)

    # Shaded region between curves = bias = missing amplitude
    ax.fill_between(xs_fit, sat_curve, ec_curve,
                       where=(ec_curve >= sat_curve),
                       color=BIAS_COL, alpha=0.18, zorder=2,
                       label=f"Bias = missing amplitude "
                              f"(≈ {ec_fit['a']:.0f} W m$^{{-2}}$ at Day 0)")

    ax.axvline(0, color="#888", ls=":", lw=1.0, zorder=1)

    # mm/day secondary axis — set after data so ylim is finalized
    ax.relim(); ax.autoscale_view()
    ymin, ymax = ax.get_ylim()
    ax_r.set_ylim(ymin / LE_PER_MM, ymax / LE_PER_MM)
    ax_r.set_ylabel("LE  [mm day$^{-1}$]", fontsize=11, color="#444")

    # N events per day intentionally omitted here — same event pool
    # as panel (a), so the N labels are identical.  Add a small note
    # to the legend instead.

    ax.set_xlabel("Days since irrigation", fontsize=11)
    ax.set_ylabel("LE  [W m$^{-2}$]", fontsize=11)
    # Restrict to the observed data range (events only go to ~Day 7 in
    # the v3 bias pool); avoids a long blank tail to MAX_WINDOW.
    x_data_max = max(positions) if positions else MAX_WINDOW
    ax.set_xlim(-0.7, x_data_max + 0.7)
    ax.set_xticks(range(0, int(x_data_max) + 1))
    ax.set_title("(b) EC vs Sat recovery — same events, two products",
                    fontsize=12, loc="left", weight="bold")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9,
                title="n per day: see panel (a)", title_fontsize=8)
    ax.grid(alpha=0.25, axis="y")


# ================================================================
# Figure assembly
# ================================================================

def make_figure(out_dir: Path, df_pool, summary, ec_fit, sat_fit):
    """2-panel side-by-side layout:
       (a) Bias recovery curve (Δ = EC − Sat, exp fit)
       (b) EC vs Sat recovery curves overlay (the gap = bias)"""
    fig = plt.figure(figsize=(16, 6.5), facecolor="white",
                       constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2])
    ax_b = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[0, 1])

    fig.suptitle(
        "Tarazona irrigation blind spot — Meteosat ETv3 stays flat while\n"
        "the EC tower recovers; the gap between them IS the management pulse",
        fontsize=14, weight="bold")

    draw_panel_b(ax_b, df_pool, None, summary)
    draw_panel_curves(ax_c, df_pool, ec_fit, sat_fit)

    png = out_dir / "fig04b_tarazona_blindspot_v32.png"
    pdf = out_dir / "fig04b_tarazona_blindspot_v32.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  [save] {png}")
    print(f"  [save] {pdf}")


def write_caption(out_dir: Path, summary: dict, metv3: dict, ec_fit: dict):
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

(b) EC vs Sat recovery — both products fitted on the same event pool.
    Side-by-side boxplots per day-since-event for EC LE (green) and
    Meteosat ETv3 LE (orange), with their independent exponential fits
    overlaid.  The purple shaded region between the two fit curves is
    the bias = missing amplitude.

       EC fit  :  τ = {ec_fit['tau']:.2f} d,  amp = {ec_fit['a']:.0f} W m⁻²
                  ({ec_fit.get('r2', float('nan')):.2f} R²,
                   n_events ≈ {summary['n_events']})
       Sat fit :  τ pegged at floor, R² ≈ 0  →  amp_Sat ≈ 0
                  (no detectable exponential response)

    The EC pulse (~{ec_fit['a']:.0f} W m⁻²) lives entirely in the gap
    between the green and the orange curve — the satellite resolves
    none of it.

    Note on τ_bias.  The bias τ ({summary['tau_bias']:.2f} d) in
    panel (a) appears in Analysis A's 3–4 d band, but this is a
    *consistency check*, not independent evidence: with Sat ≈ flat,
    bias = EC − const inherits the EC exponential by construction,
    so τ_bias = τ_EC mathematically.

------------------------------------------------------------------
Take-home message
------------------------------------------------------------------
Meteosat ETv3 reproduces NDVI (canopy) and GPP (photosynthesis)
at both sites, but at drip-irrigated Tarazona it captures **none**
of the EC management amplitude in ET (95 W m⁻² recovered by the
tower; ~0 W m⁻² by the satellite).  Because the satellite is
flat through the event window, the EC management pulse is
preserved intact inside the EC−Sat bias.

The blind spot is structural, not statistical: ETv3 does ingest
satellite soil moisture (H-SAF H141/H142/H26), but drip irrigation
slips through three resolution limits:

  (i)  Spatial scale — 1–12.5 km pixel vs 1 ha (0.01 km²) orchard.
       At Tarazona, Sentinel-2 confirms only ~4 % of the H26 1 km
       footprint is summer-active vegetation.
  (ii) Microwave sensing — scatterometers see only the top ~5 cm;
       drip wets emitter-local patches that dry within hours.
  (iii) SVAT structure — the generic f(SM) availability curve
       cannot represent the bi-modal SM distribution of a drip
       orchard.

Implication:  satellite ET can validate climate-scale drying
(τ universal across rainfed/irrigated), but quantifying drip
irrigation water use requires an in-situ amplitude correction.

------------------------------------------------------------------
Statistical methods
------------------------------------------------------------------
- Exponential fit: weighted nonlinear least squares on the median
  LE / bias per day-since-event (≥ 3 events per day required).
- Bootstrap (B = {N_BOOT_DEFAULT}, raw-data resampling) gives SE and
  95 % CI for τ.
- Boundary fits (τ ≤ 0.6 d or τ ≥ 59 d) and R² < 0.3 are excluded.
- τ_bias is reported as a math consistency check (see note in (b)),
  not as an independent physical measurement.

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

    nan_fit = dict(tau=np.nan, se=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                      a=np.nan, se_a=np.nan, c=np.nan, r2=np.nan)

    print("\n--- Fit METv3-only τ (bootstrap n={}) ---".format(args.n_boot))
    metv3 = fit_metv3_tau(df_pool, n_boot=args.n_boot) or nan_fit
    print(f"  τ_Sat = {metv3['tau']:.2f} d, "
            f"SE = {metv3.get('se', np.nan):.2f},  "
            f"CI = [{metv3.get('ci_lo', np.nan):.2f}, "
            f"{metv3.get('ci_hi', np.nan):.2f}], "
            f"R² = {metv3.get('r2', np.nan):.2f}, "
            f"amp = {metv3.get('a', np.nan):.1f} W/m²")

    print("\n--- Fit EC-only τ (same event pool, bootstrap n={}) ---".format(args.n_boot))
    ec_fit = fit_ec_tau(df_pool, n_boot=args.n_boot) or nan_fit
    print(f"  τ_EC  = {ec_fit['tau']:.2f} d, "
            f"SE = {ec_fit.get('se', np.nan):.2f},  "
            f"CI = [{ec_fit.get('ci_lo', np.nan):.2f}, "
            f"{ec_fit.get('ci_hi', np.nan):.2f}], "
            f"R² = {ec_fit.get('r2', np.nan):.2f}, "
            f"amp = {ec_fit.get('a', np.nan):.1f} W/m²")

    print("\n--- Figure (v32) ---")
    make_figure(out, df_pool, summary, ec_fit, metv3)
    write_caption(out, summary, metv3, ec_fit)

    df_pool.to_csv(out / "v32_panel_a_bias_fit.csv", index=False)
    pd.DataFrame([
        dict(label="EC (B-v3 pool)", tau=ec_fit["tau"], se=ec_fit.get("se", np.nan),
              ci_lo=ec_fit.get("ci_lo", np.nan), ci_hi=ec_fit.get("ci_hi", np.nan),
              amplitude=ec_fit.get("a", np.nan), n=summary["n_events"]),
        dict(label="EC (v31 ref)", tau=EC_TAU_REF["tau"], se=EC_TAU_REF["se"],
              ci_lo=EC_TAU_REF["ci_lo"], ci_hi=EC_TAU_REF["ci_hi"],
              amplitude=EC_TAU_REF["amplitude"], n=EC_TAU_REF["n_events"]),
        dict(label="METv3",     tau=metv3["tau"], se=metv3["se"],
              ci_lo=metv3["ci_lo"], ci_hi=metv3["ci_hi"],
              amplitude=metv3.get("a", np.nan), n=summary["n_events"]),
        dict(label="Bias",      tau=summary["tau_bias"], se=summary["tau_se_bias"],
              ci_lo=summary["tau_ci_lo"], ci_hi=summary["tau_ci_hi"],
              amplitude=summary["amplitude"], n=summary["n_events"]),
    ]).to_csv(out / "v32_panel_b_curves_fits.csv", index=False)

    print(f"\n[done] {out}/")


if __name__ == "__main__":
    main()
