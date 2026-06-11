"""
================================================================
解析A v22 : 物理境界 + 2-param model + AIC 比較
================================================================

【v21 で発覚した問題】
  - LE_∞ Percentile CI 下限 = 0.0 (物理的不合理、bounds 緩すぎ)
  - τ Percentile CI = [2.5, 16.6] (BCa: [1.9, 9.2]) が広め
  - 14 ビンで 3 自由パラメータ → overparameterized 気味

【v22 の改善】
  改善 1. 物理境界を強化
    旧: bounds=([0, 0, 0.5], [inf, inf, 100])
    新: bounds=([50, 100, 0.5], [200, 350, 60])
    - LE_∞ ≥ 50 W/m² : Oran 雨養 cereal の baseline (絶対 floor)
    - LE_0 ≥ 100 W/m² : Tarazona 通常蒸散 (灌漑直後はそれ以上)
    - τ ≤ 60 d : 30 d 超は灌漑効果でなく季節成分

  改善 2. 2-param model (LE_∞ 固定)
    v20 Part B で d ≥ 10 の median LE = 98 W/m² が観測されている
    → これを LE_∞_obs として fix、自由 param を LE_0 と τ だけに減らす
    → 14 bin / 2 param → 自由度 12 → CI 縮小期待
    → 物理的にも妥当(外挿でなく実測値で固定)

  改善 3. AIC で 3p vs 2p モデル比較
    AIC = n*log(SSE/n) + 2k
    k=2 vs k=3 で k 増分のペナルティ < SSE 改善 か?
    AIC 小 = 良い

【出力】
  v22_comparison.csv         : 3p vs 2p の τ/CI/R²/AIC
  fig01_recovery_both.png    : 2 つの fit を重ね表示
  fig02_tau_ci_compare.png   : 3p_BCa vs 2p_BCa CI 比較
  fig03_aic_residuals.png    : AIC + 残差分析
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from scipy.stats import norm

warnings.filterwarnings("ignore")


# ================================================================
# CONFIG
# ================================================================

def parse_args():
    p = argparse.ArgumentParser(description="解析A v22 — Stable bounds + 2-param model")
    p.add_argument("--parquet", default=None)
    p.add_argument("--csv", default=None)
    p.add_argument("--out", default=None)
    p.add_argument("--n-boot", type=int, default=5000)
    p.add_argument("--min-per-bin", type=int, default=5)
    p.add_argument("--le-inf-gap", type=int, default=10,
                   help="LE_∞ 観測値の定義: days_since_irrig ≥ X の中央値 (default: 10)")
    p.add_argument("--le-inf-floor", type=float, default=50.0,
                   help="3-param fit の LE_∞ 下限 W/m² (default: 50)")
    p.add_argument("--le-0-floor", type=float, default=100.0,
                   help="LE_0 下限 W/m² (default: 100)")
    p.add_argument("--tau-ceil", type=float, default=60.0,
                   help="τ 上限 d (default: 60)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--quick", action="store_true")
    return p.parse_args()


def resolve_paths(args):
    base = Path("/home/shion-nagamine")
    parquet = Path(args.parquet) if args.parquet else \
              base / "bakanposs/analysis_A/daily_classified_v4.parquet"
    csv     = Path(args.csv) if args.csv else \
              base / "Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv"
    out     = Path(args.out) if args.out else Path("./output_analysis_A_v22")
    return parquet, csv, out


IRRIG_THRESHOLD     = 0.5
MIN_IRRIG_PER_MONTH = 2
CI_PCT              = (2.5, 97.5)
FIG_BG              = "#F8F9FA"

# 物理境界(v21 からの強化点)
LE_INF_CEIL = 200.0   # W/m², LE_∞ は灌漑時 LE_0 より明らかに低いはず
LE_0_CEIL   = 350.0   # W/m², 観測最大値の安全余裕
TAU_FLOOR   = 0.5     # d


# ================================================================
# 1. データ読込・前処理
# ================================================================

def load_and_merge(parquet, csv_path):
    if not parquet.exists(): sys.exit(f"[ERROR] {parquet} not found")
    if not csv_path.exists(): sys.exit(f"[ERROR] {csv_path} not found")
    daily = pd.read_parquet(parquet)
    daily["date"] = pd.to_datetime(daily["date"])
    raw = pd.read_csv(csv_path)
    raw["date"] = pd.to_datetime(raw["date"], errors="coerce")
    irrig_cols = [c for c in ["Irrig_mm","Rain_mm","IrrigRain_mm"] if c in raw.columns]
    extra = raw[["date"]+irrig_cols].dropna(subset=["date"]).drop_duplicates("date")
    tara = daily[daily["site"]=="Tarazona"].merge(extra, on="date", how="left")
    oran = daily[daily["site"]=="Oran"].copy()
    for c in irrig_cols: oran[c] = 0.0
    df = pd.concat([oran, tara]).sort_values(["site","date"]).reset_index(drop=True)
    return df


def add_days_since_irrig(df):
    df = df.copy()
    df["days_since_irrig"] = np.nan
    for site, idx in df.groupby("site").groups.items():
        sub = df.loc[idx].sort_values("date")
        last, out = None, []
        for _, row in sub.iterrows():
            irrig = row.get("Irrig_mm", 0) or 0
            if pd.notna(irrig) and irrig > IRRIG_THRESHOLD:
                last = row["date"]; out.append(0)
            elif last is None:
                out.append(np.nan)
            else:
                out.append((row["date"] - last).days)
        df.loc[sub.index, "days_since_irrig"] = out
    return df


def add_irrig_active_month(df):
    df = df.copy()
    df["_ym"] = df["date"].dt.to_period("M")
    df["_ev"] = df["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD
    monthly = df.groupby(["site","_ym"])["_ev"].sum().reset_index(name="n_ev")
    active = set(zip(monthly.loc[monthly["n_ev"]>=MIN_IRRIG_PER_MONTH,"site"],
                     monthly.loc[monthly["n_ev"]>=MIN_IRRIG_PER_MONTH,"_ym"]))
    df["irrig_active_month"] = df.apply(lambda r: (r["site"],r["_ym"]) in active, axis=1)
    df.drop(columns=["_ym","_ev"], inplace=True)
    return df


# ================================================================
# 2. 2 つの指数モデル
# ================================================================

def exp_model_3p(d, le_inf, le0, tau):
    """3-param: LE(d) = LE_∞ + (LE_0 − LE_∞) * exp(−d/τ)"""
    return le_inf + (le0 - le_inf) * np.exp(-d/tau)


def exp_model_2p(d, le0, tau, le_inf_fixed):
    """2-param: LE_∞ を観測値に固定"""
    return le_inf_fixed + (le0 - le_inf_fixed) * np.exp(-d/tau)


def compute_le_inf_observed(arr_d, arr_y, gap_threshold=10, min_n=10):
    """
    days_since_irrig ≥ gap_threshold の生データから median LE を計算。
    これを LE_∞_obs として 2-param model に固定する。
    """
    mask = arr_d >= gap_threshold
    n = mask.sum()
    if n < min_n:
        return np.nan, n
    return float(np.median(arr_y[mask])), int(n)


# ================================================================
# 3. fit 関数(両モデル共通)
# ================================================================

def _bin(arr_d, arr_y, min_per_bin):
    """ビン中央値生成"""
    df = pd.DataFrame({"d": arr_d, "y": arr_y})
    g = df.groupby("d")["y"].agg(["median","count"]).reset_index()
    g = g[g["count"] >= min_per_bin]
    return g


def fit_3param(arr_d, arr_y, min_per_bin, bounds):
    """3-param fit、(popt, r2, aic, grouped, sse)"""
    g = _bin(arr_d, arr_y, min_per_bin)
    if len(g) < 4:
        return None, np.nan, np.nan, g, np.nan
    x = g["d"].values.astype(float)
    y = g["median"].values
    try:
        popt, _ = curve_fit(exp_model_3p, x, y, p0=[y.min(), y.max(), 5.0],
                            bounds=bounds, maxfev=8000)
        yp = exp_model_3p(x, *popt)
        sse = np.sum((y-yp)**2)
        r2 = 1 - sse / np.sum((y-y.mean())**2)
        k, n = 3, len(y)
        aic = n*np.log(sse/n + 1e-12) + 2*k
        return popt, r2, aic, g, sse
    except Exception:
        return None, np.nan, np.nan, g, np.nan


def fit_2param(arr_d, arr_y, min_per_bin, bounds, le_inf_fixed):
    """LE_∞ 固定の 2-param fit"""
    g = _bin(arr_d, arr_y, min_per_bin)
    if len(g) < 3:  # 2 param なので 3 bin あれば一応
        return None, np.nan, np.nan, g, np.nan
    x = g["d"].values.astype(float)
    y = g["median"].values
    def model(d, le0, tau):
        return exp_model_2p(d, le0, tau, le_inf_fixed)
    try:
        popt, _ = curve_fit(model, x, y, p0=[y.max(), 5.0],
                            bounds=bounds, maxfev=8000)
        yp = model(x, *popt)
        sse = np.sum((y-yp)**2)
        r2 = 1 - sse / np.sum((y-y.mean())**2)
        k, n = 2, len(y)
        aic = n*np.log(sse/n + 1e-12) + 2*k
        return popt, r2, aic, g, sse
    except Exception:
        return None, np.nan, np.nan, g, np.nan


# ================================================================
# 4. Bootstrap (両モデル並行)
# ================================================================

def bootstrap_3param(arr_d, arr_y, n_boot, min_per_bin, bounds, seed=42):
    """Method A: 生データ復元抽出 → 3-param fit"""
    n = len(arr_d)
    rng = np.random.default_rng(seed)
    taus, le0s, le_infs = [], [], []
    n_fail = 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        popt, _, _, _, _ = fit_3param(arr_d[idx], arr_y[idx], min_per_bin, bounds)
        if popt is None:
            n_fail += 1; continue
        le_infs.append(popt[0]); le0s.append(popt[1]); taus.append(popt[2])
    return dict(taus=np.array(taus), le0s=np.array(le0s), le_infs=np.array(le_infs),
                n_success=len(taus), n_fail=n_fail, n_total=n_boot)


def bootstrap_2param(arr_d, arr_y, n_boot, min_per_bin, bounds, le_inf_fixed, seed=42):
    """Method A: 生データ復元抽出 → 2-param fit"""
    n = len(arr_d)
    rng = np.random.default_rng(seed)
    taus, le0s = [], []
    n_fail = 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        popt, _, _, _, _ = fit_2param(arr_d[idx], arr_y[idx], min_per_bin,
                                       bounds, le_inf_fixed)
        if popt is None:
            n_fail += 1; continue
        le0s.append(popt[0]); taus.append(popt[1])
    return dict(taus=np.array(taus), le0s=np.array(le0s),
                n_success=len(taus), n_fail=n_fail, n_total=n_boot)


def bca_ci(boot_stats, point_stat, jackknife_stats, alpha=0.05):
    """BCa CI"""
    boot_stats = np.asarray(boot_stats)
    if len(boot_stats) < 100: return np.nan, np.nan
    prop_less = np.mean(boot_stats < point_stat)
    if prop_less == 0 or prop_less == 1: return np.nan, np.nan
    z0 = norm.ppf(prop_less)
    j = np.asarray(jackknife_stats)
    j = j[~np.isnan(j)]
    if len(j) < 10: return np.nan, np.nan
    jm = j.mean()
    num = np.sum((jm - j)**3)
    den = 6 * (np.sum((jm - j)**2))**1.5
    a = num/den if den != 0 else 0
    z_lo = norm.ppf(alpha/2)
    z_hi = norm.ppf(1 - alpha/2)
    alpha_lo = norm.cdf(z0 + (z0+z_lo)/(1-a*(z0+z_lo)))
    alpha_hi = norm.cdf(z0 + (z0+z_hi)/(1-a*(z0+z_hi)))
    return np.percentile(boot_stats, [alpha_lo*100, alpha_hi*100])


def jackknife_tau(arr_d, arr_y, min_per_bin, bounds_3p, bounds_2p, le_inf_fixed,
                   subsample=200, seed=123):
    """tau の jackknife (BCa 用、計算負荷軽減のため subsample)"""
    n = len(arr_d)
    n_sub = min(n, subsample)
    rng = np.random.default_rng(seed)
    idx_set = rng.choice(n, size=n_sub, replace=False) if n_sub < n else np.arange(n)
    jack_3p, jack_2p = [], []
    for i in idx_set:
        mask = np.ones(n, dtype=bool); mask[i] = False
        p3, _, _, _, _ = fit_3param(arr_d[mask], arr_y[mask], min_per_bin, bounds_3p)
        p2, _, _, _, _ = fit_2param(arr_d[mask], arr_y[mask], min_per_bin,
                                     bounds_2p, le_inf_fixed)
        jack_3p.append(p3[2] if p3 is not None else np.nan)
        jack_2p.append(p2[1] if p2 is not None else np.nan)
    return np.array(jack_3p), np.array(jack_2p)


# ================================================================
# 5. 可視化
# ================================================================


def _ax(ax, title="", xl="", yl=""):
    ax.set_title(title, fontweight="bold", fontsize=11, pad=6)
    ax.set_xlabel(xl, fontsize=10); ax.set_ylabel(yl, fontsize=10)
    ax.grid(alpha=0.2, lw=0.7); ax.spines[["top","right"]].set_visible(False)


def plot_recovery_both(res, out):
    """fig01: 観測 + 3p fit + 2p fit を並列表示"""
    g = res["grouped"]
    if g is None or len(g) < 2: return
    x_obs = g["d"].values.astype(float)
    y_obs = g["median"].values
    n_per = g["count"].values

    fig, ax = plt.subplots(figsize=(12, 7))
    fig.patch.set_facecolor(FIG_BG)

    x_fit = np.linspace(0, x_obs.max()+1, 200)

    # 3-param point fit
    y_3p = exp_model_3p(x_fit, res["le_inf_3p"], res["le0_3p"], res["tau_3p"])
    ax.plot(x_fit, y_3p, color="#1A73E8", lw=2.5, ls="-",
             label=f"3-param: τ={res['tau_3p']:.2f}d, LE_∞={res['le_inf_3p']:.0f}, "
                    f"R²={res['r2_3p']:.3f}")

    # 2-param point fit (LE_∞ fixed)
    y_2p = exp_model_2p(x_fit, res["le0_2p"], res["tau_2p"], res["le_inf_obs"])
    ax.plot(x_fit, y_2p, color="#E85D04", lw=2.5, ls="--",
             label=f"2-param (LE_∞ fixed): τ={res['tau_2p']:.2f}d, "
                    f"LE_∞={res['le_inf_obs']:.0f}, R²={res['r2_2p']:.3f}")

    # Observed bin medians
    sizes = 40 + n_per * 3
    ax.scatter(x_obs, y_obs, s=sizes, c="#1D9E75", edgecolors="black", lw=1,
                zorder=10, label="Median of observed bins (size ∝ n)")
    for xi, yi, ni in zip(x_obs, y_obs, n_per):
        ax.annotate(f"n={int(ni)}", (xi, yi), xytext=(5, 5),
                     textcoords="offset points", fontsize=7.5, color="#555")

    # LE_∞ horizontal line
    ax.axhline(res["le_inf_obs"], color="#7C3AED", ls=":", lw=1.8, alpha=0.7,
                label=f"LE_∞_obs = {res['le_inf_obs']:.0f} W/m² "
                       f"(median at d≥{res['le_inf_gap']}d, n={res['le_inf_n']})")

    # AIC comparison box
    dAIC = res["aic_2p"] - res["aic_3p"]
    best = "2p" if res["aic_2p"] < res["aic_3p"] else "3p"
    aic_txt = (f"AIC comparison:\n"
               f"  3p AIC = {res['aic_3p']:.1f}\n"
               f"  2p AIC = {res['aic_2p']:.1f}\n"
               f"  ΔAIC = {dAIC:+.1f}\n"
               f"  → {best} model wins")
    ax.text(0.97, 0.97, aic_txt, transform=ax.transAxes, va="top", ha="right",
             fontsize=9.5, family="monospace",
             bbox=dict(boxstyle="round,pad=0.4", fc="white",
                       ec="#1D9E75", alpha=0.95))

    ax.set_ylim(bottom=0)
    ax.legend(loc="lower left", fontsize=9)
    _ax(ax, "fig01 — Recovery curve: 3-param vs 2-param (LE_∞ fixed)",
        "Days since last irrigation", "LE_corr [W/m²]")
    plt.tight_layout()
    fp = out/"fig01_recovery_both.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_tau_ci_compare(res, out):
    """fig02: 3p vs 2p の τ + CI 比較"""
    fig, ax = plt.subplots(figsize=(11, 6))
    fig.patch.set_facecolor(FIG_BG)

    labels = ["3-param fit\n(LE_∞ free)", "2-param fit\n(LE_∞ fixed)"]
    taus = [res["tau_3p"], res["tau_2p"]]
    lo_pct = [res["tau_3p_pctile_lo"], res["tau_2p_pctile_lo"]]
    hi_pct = [res["tau_3p_pctile_hi"], res["tau_2p_pctile_hi"]]
    lo_bca = [res["tau_3p_bca_lo"], res["tau_2p_bca_lo"]]
    hi_bca = [res["tau_3p_bca_hi"], res["tau_2p_bca_hi"]]
    cols = ["#1A73E8", "#E85D04"]

    x = np.arange(len(labels))
    w = 0.35

    # 値ラベル (個々の場合に widths が異なるためバーは細く)
    for i, (lab, t, l_p, h_p, l_b, h_b, c) in enumerate(zip(
            labels, taus, lo_pct, hi_pct, lo_bca, hi_bca, cols)):
        # Percentile CI: 細い縦線
        ax.errorbar(i - w/3, t, yerr=[[t-l_p], [h_p-t]], color=c, capsize=10,
                     lw=2.5, alpha=0.55, fmt="o", markersize=10,
                     label="Percentile CI" if i==0 else None)
        # BCa CI: 太い縦線
        if not np.isnan(l_b):
            ax.errorbar(i + w/3, t, yerr=[[t-l_b], [h_b-t]], color=c, capsize=10,
                         lw=2.5, fmt="s", markersize=10,
                         label="BCa CI" if i==0 else None)
        # 値ラベル
        ax.text(i, max(h_p, h_b if not np.isnan(h_b) else h_p) + 1.5,
                 f"τ = {t:.2f} d", ha="center", fontsize=11, fontweight="bold",
                 color=c)
        ax.text(i - w/3, l_p - 1, f"Pct\n[{l_p:.1f},{h_p:.1f}]",
                 ha="center", va="top", fontsize=8, color=c)
        if not np.isnan(l_b):
            ax.text(i + w/3, l_b - 1, f"BCa\n[{l_b:.1f},{h_b:.1f}]",
                     ha="center", va="top", fontsize=8, color=c, fontweight="bold")

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(bottom=-2)
    ax.legend(loc="upper right", fontsize=9)
    _ax(ax, "fig02 — τ + 95% CI 比較 (3-param vs 2-param)\n"
            "2-param が CI 狭い = LE_∞ 固定で不確実性減",
        "", "τ [days]")
    plt.tight_layout()
    fp = out/"fig02_tau_ci_compare.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_aic_residuals(res, out):
    """fig03: AIC 比較 + 残差プロット"""
    g = res["grouped"]
    x = g["d"].values.astype(float)
    y = g["median"].values

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.patch.set_facecolor(FIG_BG)

    # 左: AIC バー
    ax = axes[0]
    aics = [res["aic_3p"], res["aic_2p"]]
    labels = ["3-param", "2-param\n(LE_∞ fixed)"]
    cols = ["#1A73E8", "#E85D04"]
    bars = ax.bar(labels, aics, color=cols, alpha=0.85,
                   edgecolor="black", lw=1.2, width=0.55)
    for bar, a in zip(bars, aics):
        ax.text(bar.get_x()+bar.get_width()/2, a+0.5, f"{a:.1f}",
                 ha="center", va="bottom", fontsize=11, fontweight="bold")
    dAIC = res["aic_2p"] - res["aic_3p"]
    txt = (f"ΔAIC = AIC(2p) − AIC(3p) = {dAIC:+.2f}\n\n"
           f"|ΔAIC| < 2: 区別不能\n"
           f"|ΔAIC| 2-7: 弱い証拠\n"
           f"|ΔAIC| > 7: 強い証拠")
    win = "2-param" if dAIC < 0 else "3-param"
    txt += f"\n\n勝者: {win}"
    ax.text(0.5, 0.05, txt, transform=ax.transAxes, ha="center", va="bottom",
             fontsize=9.5, family="monospace",
             bbox=dict(boxstyle="round,pad=0.3", fc="white",
                       ec="black", alpha=0.85))
    _ax(ax, "AIC 比較 (低いほど良い)", "", "AIC")

    # 右: 残差 vs d
    ax = axes[1]
    y_3p = exp_model_3p(x, res["le_inf_3p"], res["le0_3p"], res["tau_3p"])
    y_2p = exp_model_2p(x, res["le0_2p"], res["tau_2p"], res["le_inf_obs"])
    resid_3p = y - y_3p
    resid_2p = y - y_2p
    ax.scatter(x, resid_3p, c="#1A73E8", s=90, alpha=0.85,
                edgecolors="black", lw=1, label=f"3p residuals (SSE={res['sse_3p']:.0f})")
    ax.scatter(x, resid_2p, c="#E85D04", marker="^", s=90, alpha=0.85,
                edgecolors="black", lw=1, label=f"2p residuals (SSE={res['sse_2p']:.0f})")
    ax.axhline(0, color="black", lw=0.8)
    _ax(ax, "残差: y_obs − y_fit", "Days since last irrigation", "Residual [W/m²]")
    ax.legend(fontsize=9)

    plt.tight_layout()
    fp = out/"fig03_aic_residuals.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 6. MAIN
# ================================================================

def main():
    args = parse_args()
    if args.quick:
        args.n_boot = 500
        print("[quick mode] n_boot = 500")

    parquet, csv_path, out_dir = resolve_paths(args)
    out_dir.mkdir(parents=True, exist_ok=True)

    bounds_3p = ([args.le_inf_floor, args.le_0_floor, TAU_FLOOR],
                  [LE_INF_CEIL, LE_0_CEIL, args.tau_ceil])
    bounds_2p = ([args.le_0_floor, TAU_FLOOR],
                  [LE_0_CEIL, args.tau_ceil])

    print("="*60)
    print("解析A v22 — 物理境界 + 2-param model + AIC 比較")
    print(f"  N_BOOT = {args.n_boot}")
    print(f"  3-param bounds: LE_∞ ∈ [{args.le_inf_floor}, {LE_INF_CEIL}], "
          f"LE_0 ∈ [{args.le_0_floor}, {LE_0_CEIL}], τ ∈ [{TAU_FLOOR}, {args.tau_ceil}]")
    print(f"  2-param bounds: LE_0 ∈ [{args.le_0_floor}, {LE_0_CEIL}], "
          f"τ ∈ [{TAU_FLOOR}, {args.tau_ceil}]")
    print(f"  LE_∞ observed: median at d ≥ {args.le_inf_gap} d")
    print("="*60)

    # ── データ読込 ──
    df = load_and_merge(parquet, csv_path)
    df = add_days_since_irrig(df)
    df = add_irrig_active_month(df)

    sub = df[(df["site"]=="Tarazona") & df["irrig_active_month"] & df["is_growing"] &
              df["LE_corr"].notna() & df["days_since_irrig"].notna()].copy()
    arr_d = sub["days_since_irrig"].values.astype(float)
    arr_y = sub["LE_corr"].values
    print(f"\n  解析対象: Tarazona 灌漑アクティブ生育日 = {len(sub)} 日")

    # ── LE_∞ 観測値 ──
    le_inf_obs, le_inf_n = compute_le_inf_observed(arr_d, arr_y, args.le_inf_gap)
    print(f"  LE_∞_obs (d ≥ {args.le_inf_gap}) = {le_inf_obs:.2f} W/m² "
          f"(n={le_inf_n})")

    # ============================================================
    # Point estimates
    # ============================================================
    print(f"\n{'='*60}\n点推定 (両モデル)\n{'='*60}")
    popt_3p, r2_3p, aic_3p, grouped, sse_3p = fit_3param(
        arr_d, arr_y, args.min_per_bin, bounds_3p)
    popt_2p, r2_2p, aic_2p, _, sse_2p = fit_2param(
        arr_d, arr_y, args.min_per_bin, bounds_2p, le_inf_obs)

    if popt_3p is None or popt_2p is None:
        sys.exit("[ERROR] 点推定 fit 失敗")

    print(f"  3-param: LE_∞={popt_3p[0]:.2f}  LE_0={popt_3p[1]:.2f}  τ={popt_3p[2]:.2f}d  "
          f"R²={r2_3p:.3f}  AIC={aic_3p:.2f}  SSE={sse_3p:.1f}")
    print(f"  2-param: LE_∞={le_inf_obs:.2f}(fixed)  LE_0={popt_2p[0]:.2f}  τ={popt_2p[1]:.2f}d  "
          f"R²={r2_2p:.3f}  AIC={aic_2p:.2f}  SSE={sse_2p:.1f}")
    dAIC = aic_2p - aic_3p
    print(f"\n  ΔAIC (2p − 3p) = {dAIC:+.2f}")
    if abs(dAIC) < 2:
        print(f"  → モデル間で区別不能 (|ΔAIC|<2)。シンプルな 2-param 推奨")
    elif aic_2p < aic_3p:
        print(f"  → 2-param が優位 (LE_∞ 固定の方が parsimonious)")
    else:
        print(f"  → 3-param が優位")

    # ============================================================
    # Bootstrap
    # ============================================================
    print(f"\n{'='*60}\nBootstrap (B={args.n_boot}, Method A)\n{'='*60}")
    print(f"  3-param bootstrap 実行中...")
    b3 = bootstrap_3param(arr_d, arr_y, args.n_boot, args.min_per_bin, bounds_3p)
    print(f"    成功 {b3['n_success']}/{b3['n_total']} (失敗 {b3['n_fail']})")
    print(f"  2-param bootstrap 実行中...")
    b2 = bootstrap_2param(arr_d, arr_y, args.n_boot, args.min_per_bin,
                           bounds_2p, le_inf_obs)
    print(f"    成功 {b2['n_success']}/{b2['n_total']} (失敗 {b2['n_fail']})")

    # Percentile CI
    p3_tau_lo, p3_tau_hi = np.percentile(b3["taus"], CI_PCT)
    p2_tau_lo, p2_tau_hi = np.percentile(b2["taus"], CI_PCT)
    p3_le_inf_lo, p3_le_inf_hi = np.percentile(b3["le_infs"], CI_PCT)

    # BCa CI (jackknife)
    print(f"\n  BCa CI 計算中 (jackknife)...")
    jack_3p, jack_2p = jackknife_tau(arr_d, arr_y, args.min_per_bin,
                                      bounds_3p, bounds_2p, le_inf_obs,
                                      subsample=200)
    b3_bca_lo, b3_bca_hi = bca_ci(b3["taus"], popt_3p[2], jack_3p)
    b2_bca_lo, b2_bca_hi = bca_ci(b2["taus"], popt_2p[1], jack_2p)

    # ── 結果統合 ──
    res = dict(
        # Point
        tau_3p=popt_3p[2], le0_3p=popt_3p[1], le_inf_3p=popt_3p[0],
        r2_3p=r2_3p, aic_3p=aic_3p, sse_3p=sse_3p,
        tau_2p=popt_2p[1], le0_2p=popt_2p[0],
        le_inf_obs=le_inf_obs, le_inf_gap=args.le_inf_gap, le_inf_n=le_inf_n,
        r2_2p=r2_2p, aic_2p=aic_2p, sse_2p=sse_2p,
        # Percentile CI
        tau_3p_pctile_lo=p3_tau_lo, tau_3p_pctile_hi=p3_tau_hi,
        tau_2p_pctile_lo=p2_tau_lo, tau_2p_pctile_hi=p2_tau_hi,
        le_inf_3p_pctile_lo=p3_le_inf_lo, le_inf_3p_pctile_hi=p3_le_inf_hi,
        # BCa CI
        tau_3p_bca_lo=b3_bca_lo, tau_3p_bca_hi=b3_bca_hi,
        tau_2p_bca_lo=b2_bca_lo, tau_2p_bca_hi=b2_bca_hi,
        # Bootstrap meta
        n_data=len(sub),
        n_3p_success=b3["n_success"], n_3p_fail=b3["n_fail"],
        n_2p_success=b2["n_success"], n_2p_fail=b2["n_fail"],
        grouped=grouped,
    )

    print(f"\n  ★ 結果まとめ:")
    print(f"    3-param: τ = {popt_3p[2]:.2f} d  "
          f"Pctile [{p3_tau_lo:.2f}, {p3_tau_hi:.2f}]  "
          f"BCa [{b3_bca_lo:.2f}, {b3_bca_hi:.2f}]")
    print(f"    2-param: τ = {popt_2p[1]:.2f} d  "
          f"Pctile [{p2_tau_lo:.2f}, {p2_tau_hi:.2f}]  "
          f"BCa [{b2_bca_lo:.2f}, {b2_bca_hi:.2f}]")
    print(f"    LE_∞ (3p): {popt_3p[0]:.1f}  Pctile [{p3_le_inf_lo:.1f}, {p3_le_inf_hi:.1f}]")
    print(f"    LE_∞ (2p): {le_inf_obs:.1f} (fixed at observed median)")

    # CI 幅
    pctile_3p_width = p3_tau_hi - p3_tau_lo
    pctile_2p_width = p2_tau_hi - p2_tau_lo
    bca_3p_width = b3_bca_hi - b3_bca_lo if not np.isnan(b3_bca_lo) else np.nan
    bca_2p_width = b2_bca_hi - b2_bca_lo if not np.isnan(b2_bca_lo) else np.nan
    print(f"\n  CI 幅比較 (狭いほど精度高):")
    print(f"    3p Pctile width = {pctile_3p_width:.2f}d")
    print(f"    2p Pctile width = {pctile_2p_width:.2f}d "
          f"(縮小率 {(1-pctile_2p_width/pctile_3p_width)*100:+.0f}%)")
    if not np.isnan(bca_3p_width):
        print(f"    3p BCa width    = {bca_3p_width:.2f}d")
    if not np.isnan(bca_2p_width):
        print(f"    2p BCa width    = {bca_2p_width:.2f}d "
              f"(縮小率 {(1-bca_2p_width/bca_3p_width)*100:+.0f}%)" if not np.isnan(bca_3p_width) else "")

    # ── CSV 保存 ──
    summary_df = pd.DataFrame([
        dict(model="3-param", tau=popt_3p[2],
             tau_pctile_lo=p3_tau_lo, tau_pctile_hi=p3_tau_hi,
             tau_bca_lo=b3_bca_lo, tau_bca_hi=b3_bca_hi,
             le0=popt_3p[1], le_inf=popt_3p[0],
             le_inf_lo=p3_le_inf_lo, le_inf_hi=p3_le_inf_hi,
             r2=r2_3p, aic=aic_3p, sse=sse_3p,
             n_boot_success=b3["n_success"]),
        dict(model="2-param (LE_inf fixed)", tau=popt_2p[1],
             tau_pctile_lo=p2_tau_lo, tau_pctile_hi=p2_tau_hi,
             tau_bca_lo=b2_bca_lo, tau_bca_hi=b2_bca_hi,
             le0=popt_2p[0], le_inf=le_inf_obs,
             le_inf_lo=le_inf_obs, le_inf_hi=le_inf_obs,
             r2=r2_2p, aic=aic_2p, sse=sse_2p,
             n_boot_success=b2["n_success"]),
    ])
    summary_df.to_csv(out_dir/"v22_comparison.csv", index=False)
    grouped.to_csv(out_dir/"v22_bin_samples.csv", index=False)
    print(f"\n  [save] v22_comparison.csv + v22_bin_samples.csv")

    # ── 可視化 ──
    if not args.no_plots:
        print(f"\n--- 可視化 ---")
        plot_recovery_both(res, out_dir)
        plot_tau_ci_compare(res, out_dir)
        plot_aic_residuals(res, out_dir)

    # ============================================================
    # 論文記述用フォーマット
    # ============================================================
    print(f"\n{'='*60}\n★ 論文記述用フォーマット\n{'='*60}")
    use_2p = (aic_2p < aic_3p) or abs(dAIC) < 2
    if use_2p:
        print(f"\n  推奨記述 (2-param model 採用):")
        bca = (b2_bca_lo, b2_bca_hi) if not np.isnan(b2_bca_lo) else (p2_tau_lo, p2_tau_hi)
        ci_method = "BCa" if not np.isnan(b2_bca_lo) else "percentile"
        print(f"""
  'We modeled LE recovery after irrigation as an exponential
   decay LE(d) = LE_∞ + (LE_0 − LE_∞) * exp(−d/τ), fixing
   LE_∞ = {le_inf_obs:.0f} W/m² to the observed median LE at
   d ≥ {args.le_inf_gap} d (n = {le_inf_n}). The fit yielded
   τ = {popt_2p[1]:.2f} d (95% {ci_method} CI: {bca[0]:.2f}–{bca[1]:.2f}),
   LE_0 = {popt_2p[0]:.0f} W/m², R² = {r2_2p:.3f}.
   Bootstrap was based on raw-data resampling (B = {args.n_boot},
   {b2['n_fail']/args.n_boot*100:.1f}% iteration failure).
   AIC favored the 2-parameter model over the 3-parameter
   alternative (ΔAIC = {dAIC:+.1f}).'""")
    else:
        print(f"\n  推奨記述 (3-param model 採用):")
        bca = (b3_bca_lo, b3_bca_hi) if not np.isnan(b3_bca_lo) else (p3_tau_lo, p3_tau_hi)
        ci_method = "BCa" if not np.isnan(b3_bca_lo) else "percentile"
        print(f"""
  'τ = {popt_3p[2]:.2f} d (95% {ci_method} CI: {bca[0]:.2f}–{bca[1]:.2f}),
   LE_∞ = {popt_3p[0]:.0f} W/m², LE_0 = {popt_3p[1]:.0f} W/m²,
   R² = {r2_3p:.3f}, B = {args.n_boot}.'""")

    print(f"\n[done] outputs in {out_dir}/")


if __name__ == "__main__":
    main()
