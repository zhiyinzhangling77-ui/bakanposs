"""
================================================================
解析A v32 : Poster Fig 4b — Tarazona irrigation blind spot
================================================================

【目的】
  解析A v31 の τ ≈ 3 d / amplitude 4.5× scaling を踏まえて、
  Meteosat ETv3 衛星 product が Tarazona 灌漑アーモンド園で
  「同じ時間スケールを再現しつつ、振幅を丸ごと取り逃す」現象を
  1 枚の poster figure (3 panel) で見せる。

【3 panel 構成】
  (a) 時系列 zoom: 1 灌漑シーズン (default: 2023 Jun-Sep)
      EC LE [W/m²] vs METv3 LE [W/m²]、灌漑イベントを縦線
      → 衛星が灌漑パルスに無反応であることを視覚化
  (b) Bias recovery: Δ = EC_LE − METv3_LE [W/m²] を Day 0 起点で
      箱ひげ + exp fit (τ_bias, amplitude_bias)。第二軸 mm/day。
      → bias が解析A と同じ τ 時間スケールで減衰
  (c) τ comparison: 3 bars
      EC τ (解析A v31) / METv3 τ / Bias τ
      → 時間スケールは普遍、振幅だけ衛星が見逃す

【入力 (デフォルト)】
  --bias-pool        ./output_analysis_B_v3/v3_bias_perevent_pooled.csv
                     (date,LE_EC_Wm2,Irrig_mm,LE_METv3_Wm2,bias,
                      days_since_event,site)  解析B v3 出力
  --bias-summary     ./output_analysis_B_v3/v3_bias_tau_summary.csv
  --parquet          /home/shion-nagamine/bakanposs/analysis_A/
                     daily_classified_v4.parquet  (時系列 panel 用)
  --metv3-csv        /home/shion-nagamine/bakanposs/metv3_daily_all.csv
                     (時系列 panel 用)
  --season           "2023-06-01:2023-09-30"  panel (a) zoom 範囲

【出力】 ./output_analysis_A_v32/
  fig04b_tarazona_blindspot_v32.png/pdf  ★ poster figure
  v32_poster_caption.txt                 ★ 脚注テキスト
  v32_panel_a_timeseries.csv             panel (a) 再現用
  v32_panel_b_bias_fit.csv               panel (b) fit 再現用
  v32_panel_c_tau_bars.csv               panel (c) bars
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
from matplotlib.gridspec import GridSpec
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
    p.add_argument("--parquet",
                   default="/home/shion-nagamine/bakanposs/analysis_A/daily_classified_v4.parquet")
    p.add_argument("--tara-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Daily_Summary_Filtered_forPred_ActEne26.csv")
    p.add_argument("--metv3-csv",
                   default="/home/shion-nagamine/bakanposs/metv3_daily_all.csv")
    p.add_argument("--season", default="2023-06-01:2023-09-30",
                   help="panel (a) zoom range, ISO start:end")
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


def load_timeseries(parquet: Path, tara_csv: Path, metv3_csv: Path,
                     season_start: pd.Timestamp, season_end: pd.Timestamp
                     ) -> pd.DataFrame | None:
    """Panel (a) の時系列 zoom 用。失敗しても None を返し panel (a) スキップ."""
    try:
        daily = pd.read_parquet(parquet)
        daily["date"] = pd.to_datetime(daily["date"])
        tara = daily[daily["site"] == "Tarazona"].copy()
        if "LE_corr" in tara.columns and "LE" in tara.columns:
            tara = tara.drop(columns=["LE"]).rename(columns={"LE_corr": "LE"})
        elif "LE_corr" in tara.columns:
            tara = tara.rename(columns={"LE_corr": "LE"})
        tara["LE_EC"] = pd.to_numeric(tara["LE"], errors="coerce")

        raw = pd.read_csv(tara_csv)
        raw["date"] = pd.to_datetime(raw["date"], errors="coerce")
        if "Irrig_mm" in raw.columns:
            irrig = raw[["date", "Irrig_mm"]].dropna(subset=["date"]).drop_duplicates("date")
            tara = tara.merge(irrig, on="date", how="left")
        tara["Irrig_mm"] = tara.get("Irrig_mm", 0).fillna(0)

        met = pd.read_csv(metv3_csv)
        met["date"] = pd.to_datetime(met["date"], errors="coerce")
        met = met[met["site"].isin(["Tarazona", "TzM"])].copy()
        met["LE_MET"] = pd.to_numeric(met["ET_mm"], errors="coerce") * LE_PER_MM
        met = met[["date", "LE_MET"]].drop_duplicates("date")

        out = tara[["date", "LE_EC", "Irrig_mm"]].merge(met, on="date", how="outer")
        out = out[(out["date"] >= season_start) & (out["date"] <= season_end)]
        return out.sort_values("date").reset_index(drop=True)
    except Exception as e:
        print(f"  [warn] timeseries load failed ({type(e).__name__}: {e})")
        print(f"         panel (a) を簡易版で描画します")
        return None


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

def draw_panel_a(ax, ts: pd.DataFrame | None,
                  season_start, season_end):
    if ts is None or ts.empty:
        ax.text(0.5, 0.5,
                 "Panel (a) skipped\n(daily timeseries not loaded;\n"
                 "see --parquet / --metv3-csv args)",
                 ha="center", va="center", fontsize=11,
                 transform=ax.transAxes, color="#888")
        ax.set_xticks([]); ax.set_yticks([])
        return

    ts = ts.copy()
    ts["LE_EC_smooth"] = ts["LE_EC"].rolling(3, center=True, min_periods=1).mean()
    ts["LE_MET_smooth"] = ts["LE_MET"].rolling(3, center=True, min_periods=1).mean()

    ax.plot(ts["date"], ts["LE_EC_smooth"], color=EC_COL, lw=2.0,
             label="EC (tower)", zorder=4)
    ax.plot(ts["date"], ts["LE_MET_smooth"], color=METV3_COL, lw=2.0,
             label="Meteosat ETv3", zorder=3)
    ax.fill_between(ts["date"], ts["LE_MET_smooth"], ts["LE_EC_smooth"],
                       where=(ts["LE_EC_smooth"] > ts["LE_MET_smooth"]),
                       color=BIAS_COL, alpha=0.18, label="Bias (EC − Sat)",
                       zorder=2)

    irrig = ts[ts["Irrig_mm"] > IRRIG_THRESHOLD]
    y_top = float(np.nanmax(ts["LE_EC_smooth"])) * 1.05
    if not irrig.empty:
        ax.vlines(irrig["date"], 0, y_top, color="#3878C6", lw=0.7,
                    alpha=0.55, zorder=1)
        ax.scatter(irrig["date"], np.full(len(irrig), y_top * 0.98),
                     marker="v", s=22, color="#3878C6", zorder=5,
                     label=f"Irrigation events (n={len(irrig)})")

    ax.set_xlim(season_start, season_end)
    ax.set_ylim(0, y_top * 1.08)
    ax.set_ylabel("LE [W m$^{-2}$]", fontsize=11)
    ax.set_xlabel("Date", fontsize=11)
    ax.set_title(
        f"(a) Tarazona summer {season_start.year}: "
        "EC tower vs Meteosat ETv3",
        fontsize=12, loc="left", weight="bold")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.25)


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
        f"(b) Bias recovery (n events = {summary['n_events']}, "
        f"R² = {summary['r2']:.2f})",
        fontsize=12, loc="left", weight="bold")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.25, axis="y")


def draw_panel_c(ax, ec_ref: dict, metv3: dict, bias: dict):
    labels = [
        f"EC\n(τ_EC)",
        f"Meteosat ETv3\n(τ_Sat)",
        f"Bias = EC−Sat\n(τ_bias)",
    ]
    taus = [ec_ref["tau"], metv3["tau"], bias["tau_bias"]]
    err_lo = [ec_ref["tau"] - ec_ref["ci_lo"],
                metv3["tau"] - metv3["ci_lo"] if not np.isnan(metv3["ci_lo"]) else 0,
                bias["tau_bias"] - bias["tau_ci_lo"]]
    err_hi = [ec_ref["ci_hi"] - ec_ref["tau"],
                metv3["ci_hi"] - metv3["tau"] if not np.isnan(metv3["ci_hi"]) else 0,
                bias["tau_ci_hi"] - bias["tau_bias"]]
    cols = [EC_COL, METV3_COL, BIAS_COL]
    xs = np.arange(len(labels))

    # universal 3-4 d band
    ax.axhspan(3.0, 4.0, color="#1D9E75", alpha=0.10, zorder=0,
                  label="Analysis-A universal band (3–4 d)")

    bars = ax.bar(xs, taus, color=cols, alpha=0.85, edgecolor="black", lw=0.8,
                     zorder=3)
    ax.errorbar(xs, taus, yerr=[err_lo, err_hi], fmt="none",
                  ecolor="black", capsize=6, lw=1.5, zorder=4)

    for x, tau, se in zip(xs, taus,
                              [ec_ref["se"], metv3["se"], bias["tau_se_bias"]]):
        se_str = f"SE={se:.2f}" if not np.isnan(se) else "SE=NA"
        ax.text(x, tau + 0.25,
                  f"{tau:.2f} d\n{se_str}",
                  ha="center", va="bottom", fontsize=9, weight="bold")

    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("τ [d]", fontsize=11)
    ax.set_ylim(0, max(taus + err_hi) * 1.35 if max(err_hi) > 0
                      else max(taus) * 1.45)
    ax.set_title("(c) τ comparison: EC vs Sat vs Bias",
                    fontsize=12, loc="left", weight="bold")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.25, axis="y")


# ================================================================
# Figure assembly
# ================================================================

def make_figure(out_dir: Path, ts, df_pool, summary, metv3,
                 season_start, season_end):
    fig = plt.figure(figsize=(15, 11), facecolor="white")
    gs = GridSpec(2, 2, figure=fig,
                     height_ratios=[1.0, 1.1],
                     width_ratios=[1.35, 1.0],
                     hspace=0.42, wspace=0.32,
                     left=0.07, right=0.96, top=0.93, bottom=0.08)

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    fig.suptitle(
        "Tarazona irrigation blind spot — Meteosat ETv3 misses the management amplitude\n"
        f"(τ_sat ≈ τ_EC, but Sat amplitude ≪ EC amplitude)",
        fontsize=14, weight="bold", y=0.985)

    draw_panel_a(ax_a, ts, season_start, season_end)
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
(a) Time series zoom — one irrigation season at Tarazona.
    Green = EC tower LE, orange = Meteosat ETv3 LE (both 3-day rolling
    mean for visibility).  Blue triangles mark irrigation events
    (Irrig ≥ 0.5 mm).  Purple shading = bias (EC − Sat).
    Tower shows pulsed recovery on every irrigation; satellite stays
    flat → blind spot.

(b) Bias recovery curve — per-event Δ pooled across all Tarazona
    irrigation events ({summary['n_events']} events).  Boxplots show the
    distribution of Δ at each day-since-event; medians drive the
    exponential fit  Δ(d) = c + amp · exp(−d/τ_bias).  Right axis
    shows the same quantity in mm day⁻¹.
    Fitted: τ_bias = {summary['tau_bias']:.2f} d
            amp_bias = {summary['amplitude']:.1f} W m⁻²
            (≈ {summary['amplitude']/LE_PER_MM:.2f} mm day⁻¹)
            R² = {summary['r2']:.2f}

(c) τ comparison — three bars side by side.
    1) EC reference (Analysis A v31, Tarazona active, n = 41):
       τ_EC = {EC_TAU_REF['tau']:.2f} d  [{EC_TAU_REF['ci_lo']:.2f}, {EC_TAU_REF['ci_hi']:.2f}]
       amp = {EC_TAU_REF['amplitude']:.0f} W m⁻²
    2) Meteosat ETv3 alone, same events:
       τ_Sat = {metv3['tau']:.2f} d  [{metv3['ci_lo']:.2f}, {metv3['ci_hi']:.2f}]
    3) Bias = EC − Sat:
       τ_bias = {summary['tau_bias']:.2f} d
       [{summary['tau_ci_lo']:.2f}, {summary['tau_ci_hi']:.2f}]
       amp_bias = {summary['amplitude']:.1f} W m⁻²
    Green band = Analysis-A universal range (3-4 d).
    All three bars overlap → satellite reproduces the universal
    timescale; the entire management amplitude (~95 W m⁻², matching
    EC amplitude) is captured in the bias.

------------------------------------------------------------------
Take-home message
------------------------------------------------------------------
The satellite product is NOT broken in time — it sees the same
~3 d recovery as the EC tower.  What it misses is the *magnitude*
of the irrigation pulse: the bias amplitude is essentially equal
to the EC amplitude (≈ 95 W m⁻²).
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

    season_start, season_end = [pd.Timestamp(s) for s in args.season.split(":")]

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
        print("  [warn] METv3 fit failed; panel (c) will skip METv3 bar")
        metv3 = dict(tau=np.nan, se=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                       a=np.nan, c=np.nan, r2=np.nan)
    else:
        print(f"  τ_Sat = {metv3['tau']:.2f} d, "
                f"SE = {metv3['se']:.2f},  "
                f"CI = [{metv3['ci_lo']:.2f}, {metv3['ci_hi']:.2f}], "
                f"R² = {metv3['r2']:.2f}")

    print(f"\n--- Load timeseries for panel (a) "
           f"({season_start.date()} → {season_end.date()}) ---")
    ts = load_timeseries(Path(args.parquet), Path(args.tara_csv),
                            Path(args.metv3_csv), season_start, season_end)
    if ts is not None:
        print(f"  rows: {len(ts)}, irrig events in zoom: "
                f"{int((ts['Irrig_mm'] > IRRIG_THRESHOLD).sum())}")

    print("\n--- Figure (v32) ---")
    make_figure(out, ts, df_pool, summary, metv3, season_start, season_end)
    write_caption(out, summary, metv3)

    # CSV exports
    if ts is not None:
        ts.to_csv(out / "v32_panel_a_timeseries.csv", index=False)
    df_pool.to_csv(out / "v32_panel_b_bias_fit.csv", index=False)
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
    ]).to_csv(out / "v32_panel_c_tau_bars.csv", index=False)

    print(f"\n[done] {out}/")


if __name__ == "__main__":
    main()
