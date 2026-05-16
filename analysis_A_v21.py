"""
================================================================
解析A v21 : τ Bootstrap CI(Method A: 生データ復元抽出)
================================================================

【目的】
  v15-v20 で τ は点推定のみで報告していた。論文化には CI 必須。
  他のAI レビューに従い、**Method A: 生データレベルの復元抽出**で
  τ・LE_0・LE_∞ の bootstrap CI を計算する。

【Method A vs Method B】
  ❌ Method B: ビン中央値だけを復元抽出 → ビン化で圧縮された
              情報しか反映しない、CI が過小評価になる
  ✓ Method A: 元の 589 日 LE 観測値を復元抽出 → 各 iteration で
              ビン再構築 → curve_fit
              → ビン内ばらつきと ビン化情報損失の両方が CI に反映

【実装上の重要点(レビュー指摘事項)】
  1. ビンごとのサンプル数を必ず報告(透明性)
  2. curve_fit 収束失敗率を追跡し、>10% なら警告
  3. B (bootstrap 回数) は 5000 をデフォルト、最大 10000
  4. CI 計算は percentile(default)と BCa(option)両方サポート
  5. 失敗 iteration は除外し、有効サンプル数を報告

【出力】
  v21_recovery_bootstrap.csv  全期間 + 年別 τ/LE_0/LE_∞ + CI + 失敗率
  v21_bin_samples.csv         ビン毎のサンプル数(透明性)
  fig01_recovery_with_ci.png  Recovery curve + bootstrap CI band
  fig02_tau_bootstrap_dist.png τ ブートストラップ分布 histogram
  fig03_per_year_tau_ci.png   年別 τ + CI
  fig04_bin_samples.png       ビンサイズの可視化
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import norm

warnings.filterwarnings("ignore")


# ================================================================
# CONFIG
# ================================================================

def parse_args():
    p = argparse.ArgumentParser(description="解析A v21 — τ Bootstrap CI (Method A)")
    p.add_argument("--parquet", default=None)
    p.add_argument("--csv", default=None)
    p.add_argument("--out", default=None)
    p.add_argument("--years", type=int, nargs="+",
                   default=[2020, 2021, 2022, 2023, 2024])
    p.add_argument("--n-boot", type=int, default=5000,
                   help="Bootstrap 回数 (default: 5000, 推奨範囲 1000-10000)")
    p.add_argument("--min-per-bin", type=int, default=5,
                   help="ビンあたり最小サンプル数 (default: 5)")
    p.add_argument("--ci-method", choices=["percentile","bca","both"], default="percentile",
                   help="CI 計算方法 (default: percentile)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--quick", action="store_true",
                   help="クイックモード: n_boot=500 (動作確認用)")
    return p.parse_args()


def resolve_paths(args):
    base = Path("/home/shion-nagamine")
    parquet = Path(args.parquet) if args.parquet else \
              base / "bakanposs/analysis_A/daily_classified_v4.parquet"
    csv     = Path(args.csv) if args.csv else \
              base / "Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv"
    out     = Path(args.out) if args.out else Path("./output_analysis_A_v21")
    return parquet, csv, out


IRRIG_THRESHOLD     = 0.5
MIN_IRRIG_PER_MONTH = 2
CI_PCT              = (2.5, 97.5)
FIG_BG              = "#F8F9FA"
R2_MIN_TRUST        = 0.70


# ================================================================
# 1. データ読込・前処理(v20 と互換)
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
    df["year"] = df["date"].dt.year
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
# 2. 指数モデル + 単一 fit
# ================================================================

def exp_model(d, le_inf, le0, tau):
    """LE(d) = LE_∞ + (LE_0 − LE_∞) * exp(−d/τ)"""
    return le_inf + (le0 - le_inf) * np.exp(-d / tau)


def _bin_and_fit(arr_d, arr_y, min_per_bin):
    """
    生データ (days_since_irrig, LE) からビン中央値を作り curve_fit する。
    返り値: (popt[le_inf, le0, tau], r2, grouped_df) or (None, NaN, None)
    """
    df = pd.DataFrame({"d": arr_d, "y": arr_y})
    grouped = df.groupby("d")["y"].agg(["median","count"]).reset_index()
    grouped = grouped[grouped["count"] >= min_per_bin]
    if len(grouped) < 4:
        return None, np.nan, grouped
    x = grouped["d"].values.astype(float)
    y = grouped["median"].values
    try:
        popt, _ = curve_fit(exp_model, x, y,
                            p0=[y.min(), y.max(), 5.0],
                            bounds=([0, 0, 0.5], [np.inf, np.inf, 100]),
                            maxfev=8000)
        yp = exp_model(x, *popt)
        r2 = 1 - np.sum((y-yp)**2) / np.sum((y-y.mean())**2)
        return popt, r2, grouped
    except Exception:
        return None, np.nan, grouped


# ================================================================
# 3. Method A Bootstrap (生データ復元抽出)
# ================================================================

def bootstrap_recovery(arr_d, arr_y, n_boot=5000, min_per_bin=5, seed=42, verbose=True):
    """
    Method A: 生データ復元抽出 → ビン中央値再構築 → curve_fit
    各 iteration で:
        1. 生データ N 個から復元抽出 N 個
        2. days_since_irrig でグループ化、ビン中央値計算 (≥min_per_bin)
        3. curve_fit で τ・LE_0・LE_∞ 推定
    収束しない iteration は除外し、失敗率を報告。
    """
    n = len(arr_d)

    # Point estimate (元データそのまま)
    popt_pt, r2_pt, grouped_pt = _bin_and_fit(arr_d, arr_y, min_per_bin)
    if popt_pt is None:
        return dict(success=False, reason="point estimate fit failed",
                    n_data=n, grouped=grouped_pt)

    # Bootstrap
    rng = np.random.default_rng(seed)
    boot_taus, boot_le0s, boot_le_infs = [], [], []
    n_fail = 0
    n_low_bin = 0  # ビン数不足

    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        b_d, b_y = arr_d[idx], arr_y[idx]
        popt, r2, grp = _bin_and_fit(b_d, b_y, min_per_bin)
        if popt is None:
            if grp is not None and len(grp) < 4:
                n_low_bin += 1
            else:
                n_fail += 1
            continue
        boot_le_infs.append(popt[0])
        boot_le0s.append(popt[1])
        boot_taus.append(popt[2])

    n_success = len(boot_taus)
    failure_rate = (n_boot - n_success) / n_boot

    if verbose:
        print(f"    Bootstrap: {n_success}/{n_boot} 成功 "
              f"(失敗率 {failure_rate*100:.1f}%, ビン数不足 {n_low_bin})")
        if failure_rate > 0.10:
            print(f"    ⚠ 失敗率 >10% → 結果に注意 (CI が広め/不安定)")

    if n_success < n_boot * 0.5:
        return dict(success=False, reason=f"only {n_success}/{n_boot} succeeded",
                    n_data=n, grouped=grouped_pt, failure_rate=failure_rate,
                    tau=popt_pt[2], le0=popt_pt[1], le_inf=popt_pt[0], r2=r2_pt)

    # Percentile CI
    tau_lo, tau_hi       = np.percentile(boot_taus,    CI_PCT)
    le0_lo, le0_hi       = np.percentile(boot_le0s,    CI_PCT)
    le_inf_lo, le_inf_hi = np.percentile(boot_le_infs, CI_PCT)

    return dict(
        success=True,
        n_data=n,
        n_boot_success=n_success,
        n_boot_total=n_boot,
        failure_rate=failure_rate,
        n_low_bin=n_low_bin,
        # Point estimates
        tau=popt_pt[2], le0=popt_pt[1], le_inf=popt_pt[0], r2=r2_pt,
        # Percentile CIs
        tau_lo=tau_lo, tau_hi=tau_hi,
        le0_lo=le0_lo, le0_hi=le0_hi,
        le_inf_lo=le_inf_lo, le_inf_hi=le_inf_hi,
        # Raw bootstrap arrays(plot用)
        boot_taus=np.array(boot_taus),
        boot_le0s=np.array(boot_le0s),
        boot_le_infs=np.array(boot_le_infs),
        grouped=grouped_pt,
    )


# ================================================================
# 4. BCa CI(optional, より精度高い)
# ================================================================

def bca_ci_for_tau(arr_d, arr_y, boot_taus, point_tau, min_per_bin=5,
                    alpha=0.05, jackknife_n=None):
    """
    Bias-Corrected and Accelerated bootstrap CI for τ.
    Percentile に比べて分布の歪み・バイアスを補正する。
    Jackknife が計算負荷高いので、N>200 では random subsample で代用。
    """
    n = len(arr_d)
    # Bias correction
    prop_less = np.mean(boot_taus < point_tau)
    if prop_less == 0 or prop_less == 1:
        return np.nan, np.nan
    z0 = norm.ppf(prop_less)

    # Acceleration via jackknife (or subsample if N large)
    if jackknife_n is None:
        jackknife_n = min(n, 200)
    if jackknife_n < n:
        # ランダムサブサンプリングで近似
        rng = np.random.default_rng(123)
        jack_idx_set = rng.choice(n, size=jackknife_n, replace=False)
    else:
        jack_idx_set = np.arange(n)

    jack_taus = []
    for i in jack_idx_set:
        mask = np.ones(n, dtype=bool); mask[i] = False
        popt, _, _ = _bin_and_fit(arr_d[mask], arr_y[mask], min_per_bin)
        if popt is not None:
            jack_taus.append(popt[2])
    if len(jack_taus) < jackknife_n * 0.5:
        return np.nan, np.nan
    jack_taus = np.array(jack_taus)
    jack_mean = jack_taus.mean()
    num = np.sum((jack_mean - jack_taus)**3)
    den = 6 * (np.sum((jack_mean - jack_taus)**2))**1.5
    a = num/den if den != 0 else 0

    z_lo = norm.ppf(alpha/2)
    z_hi = norm.ppf(1 - alpha/2)
    alpha_lo = norm.cdf(z0 + (z0 + z_lo) / (1 - a*(z0 + z_lo)))
    alpha_hi = norm.cdf(z0 + (z0 + z_hi) / (1 - a*(z0 + z_hi)))
    return np.percentile(boot_taus, [alpha_lo*100, alpha_hi*100])


# ================================================================
# 5. 可視化
# ================================================================

def _ax(ax, title="", xl="", yl=""):
    ax.set_title(title, fontweight="bold", fontsize=11, pad=6)
    ax.set_xlabel(xl, fontsize=10); ax.set_ylabel(yl, fontsize=10)
    ax.grid(alpha=0.2, lw=0.7); ax.spines[["top","right"]].set_visible(False)


def plot_recovery_with_ci(boot_res, out):
    """fig01: Recovery curve + bootstrap CI band"""
    grouped = boot_res["grouped"]
    if grouped is None or len(grouped)<2: return
    x_obs = grouped["d"].values.astype(float)
    y_obs = grouped["median"].values
    n_per_bin = grouped["count"].values

    fig, ax = plt.subplots(figsize=(12, 7))
    fig.patch.set_facecolor(FIG_BG)

    # ── Bootstrap CI band: 各 d でモデル値の分布から ──
    x_fit = np.linspace(0, x_obs.max()+1, 200)
    boot_curves = np.array([
        exp_model(x_fit, le_inf, le0, tau)
        for le_inf, le0, tau in zip(
            boot_res["boot_le_infs"], boot_res["boot_le0s"], boot_res["boot_taus"])
    ])
    band_lo = np.percentile(boot_curves, CI_PCT[0], axis=0)
    band_hi = np.percentile(boot_curves, CI_PCT[1], axis=0)
    median_curve = np.percentile(boot_curves, 50, axis=0)

    ax.fill_between(x_fit, band_lo, band_hi, color="#1D9E75", alpha=0.18,
                     label=f"Bootstrap 95% CI band (n={boot_res['n_boot_success']})")
    ax.plot(x_fit, median_curve, "k-", lw=2.2,
             label=f"Bootstrap median curve")

    # Point estimate
    y_point = exp_model(x_fit, boot_res["le_inf"], boot_res["le0"], boot_res["tau"])
    ax.plot(x_fit, y_point, ls="--", color="#7C3AED", lw=2,
             label=f"Point fit (τ={boot_res['tau']:.1f}d, R²={boot_res['r2']:.3f})")

    # 観測中央値
    sizes = 40 + n_per_bin*3
    ax.scatter(x_obs, y_obs, s=sizes, c="#E85D04", edgecolors="black", lw=1,
                zorder=10, label=f"観測中央値 (サイズ ∝ ビン n)")
    # ビン n を点の隣に書く
    for xi, yi, ni in zip(x_obs, y_obs, n_per_bin):
        ax.annotate(f"n={int(ni)}", (xi, yi), xytext=(6, 6),
                     textcoords="offset points", fontsize=7.5, color="#555")

    # τ と CI を text で
    ci_str = f"τ = {boot_res['tau']:.2f} d  [95%CI: {boot_res['tau_lo']:.2f}, {boot_res['tau_hi']:.2f}]"
    ax.text(0.97, 0.97, ci_str, transform=ax.transAxes, va="top", ha="right",
             fontsize=11, fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#1D9E75", alpha=0.95))

    ax.set_ylim(bottom=0)
    ax.legend(loc="lower right", fontsize=9)
    _ax(ax, "fig01 — Recovery curve + Bootstrap CI band (Method A: raw data resampling)",
        "Days since last irrigation", "LE_corr [W/m²]")
    plt.tight_layout()
    fp = out/"fig01_recovery_with_ci.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_tau_bootstrap_dist(boot_res, out):
    """fig02: τ ブートストラップ分布 histogram"""
    taus = boot_res["boot_taus"]
    if len(taus) < 100: return
    fig, ax = plt.subplots(figsize=(11, 6))
    fig.patch.set_facecolor(FIG_BG)
    ax.hist(taus, bins=50, color="#1D9E75", alpha=0.75, edgecolor="black", lw=0.8)
    ax.axvline(boot_res["tau"], color="#E85D04", ls="-", lw=2.5,
                label=f"Point estimate τ = {boot_res['tau']:.2f}d")
    ax.axvline(boot_res["tau_lo"], color="black", ls="--", lw=1.5, alpha=0.7,
                label=f"95% CI lo = {boot_res['tau_lo']:.2f}d")
    ax.axvline(boot_res["tau_hi"], color="black", ls="--", lw=1.5, alpha=0.7,
                label=f"95% CI hi = {boot_res['tau_hi']:.2f}d")
    ax.axvline(np.median(taus), color="#7C3AED", ls=":", lw=2,
                label=f"Bootstrap median = {np.median(taus):.2f}d")

    info_txt = (f"N bootstrap = {boot_res['n_boot_success']}/{boot_res['n_boot_total']}\n"
                f"失敗率 = {boot_res['failure_rate']*100:.1f}%\n"
                f"ビン数不足 iter = {boot_res['n_low_bin']}")
    ax.text(0.97, 0.97, info_txt, transform=ax.transAxes, va="top", ha="right",
             fontsize=9, family="monospace",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))

    ax.legend(loc="upper left", fontsize=9)
    _ax(ax, "fig02 — τ Bootstrap distribution", "τ [days]", "Frequency")
    plt.tight_layout()
    fp = out/"fig02_tau_bootstrap_dist.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_per_year_tau_ci(year_results_df, tau_all, out):
    """fig03: 年別 τ + bootstrap CI"""
    df = year_results_df.dropna(subset=["tau"])
    if len(df) < 1: return
    fig, ax = plt.subplots(figsize=(11, 6.5))
    fig.patch.set_facecolor(FIG_BG)
    x = np.arange(len(df))
    colors = ["#1D9E75" if (not np.isnan(r2) and r2>=R2_MIN_TRUST) else "#BDBDBD"
              for r2 in df["r2"]]
    ax.bar(x, df["tau"], color=colors, alpha=0.85, edgecolor="black", lw=1.2,
            width=0.6)
    for i, (_, row) in enumerate(df.iterrows()):
        if not np.isnan(row["tau_lo"]):
            ax.errorbar(i, row["tau"],
                         yerr=[[row["tau"]-row["tau_lo"]],[row["tau_hi"]-row["tau"]]],
                         color="black", capsize=10, lw=1.6, fmt="none")
            ax.text(i, row["tau_hi"]+0.5,
                     f"{row['tau']:.1f}\n[{row['tau_lo']:.1f},{row['tau_hi']:.1f}]\n"
                     f"R²={row['r2']:.2f} n={int(row['n_data'])}",
                     ha="center", va="bottom", fontsize=8.5,
                     color="black" if (not np.isnan(row['r2']) and row['r2']>=R2_MIN_TRUST)
                     else "#888")
    ax.axhline(tau_all, color="#E85D04", ls="--", lw=2,
                label=f"全期間 τ = {tau_all:.1f}d")
    ax.axhline(5, color="gray", ls=":", lw=1, alpha=0.6, label="τ = 5d 閾値")
    ax.set_xticks(x); ax.set_xticklabels(df["year"].astype(int), fontsize=10)
    leg = [mpatches.Patch(color="#1D9E75", label="信頼年 (R²≥0.7)"),
           mpatches.Patch(color="#BDBDBD", label="低信頼年")]
    h, _ = ax.get_legend_handles_labels()
    ax.legend(handles=leg+h, fontsize=9, loc="upper right")
    _ax(ax, "fig03 — 年別 τ + Bootstrap 95% CI (Method A)", "Year", "τ [days]")
    plt.tight_layout()
    fp = out/"fig03_per_year_tau_ci.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_bin_samples(grouped, out):
    """fig04: ビンごとのサンプル数(透明性)"""
    if grouped is None or len(grouped)<2: return
    fig, ax = plt.subplots(figsize=(11, 5.5))
    fig.patch.set_facecolor(FIG_BG)
    x = grouped["d"].values.astype(int)
    n = grouped["count"].values
    bars = ax.bar(x, n, color="#1D9E75", alpha=0.85, edgecolor="black", lw=1)
    for xi, ni in zip(x, n):
        ax.text(xi, ni+max(n)*0.02, str(int(ni)),
                 ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(x)
    ax.axhline(5, color="red", ls="--", lw=1, alpha=0.6, label="min ≥5")
    ax.legend(fontsize=9)
    _ax(ax, f"fig04 — ビンごとのサンプル数 (全 {n.sum()} 日)\n"
            "ビン中央値の安定性 ∝ サンプル数",
        "Days since last irrigation", "Number of daily observations")
    plt.tight_layout()
    fp = out/"fig04_bin_samples.png"
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

    print("="*60)
    print(f"解析A v21 — τ Bootstrap CI (Method A)")
    print(f"  N_BOOT = {args.n_boot}")
    print(f"  min_per_bin = {args.min_per_bin}")
    print(f"  CI method = {args.ci_method}")
    print(f"  出力先: {out_dir}")
    print("="*60)

    # ── データ読込 ──
    df = load_and_merge(parquet, csv_path)
    df = add_days_since_irrig(df)
    df = add_irrig_active_month(df)

    # ── Tarazona 灌漑アクティブ期 + 生育期 ──
    sub = df[(df["site"]=="Tarazona") & df["irrig_active_month"] & df["is_growing"] &
              df["LE_corr"].notna() & df["days_since_irrig"].notna()].copy()
    arr_d = sub["days_since_irrig"].values.astype(float)
    arr_y = sub["LE_corr"].values
    print(f"\n  解析対象: Tarazona 灌漑アクティブ生育日 = {len(sub)} 日")

    # ============================================================
    # 全期間 τ Bootstrap
    # ============================================================
    print(f"\n{'='*60}\n全期間 τ Bootstrap (Method A)\n{'='*60}")
    boot_res = bootstrap_recovery(arr_d, arr_y,
                                   n_boot=args.n_boot,
                                   min_per_bin=args.min_per_bin)
    if not boot_res["success"]:
        print(f"  ★ Bootstrap 失敗: {boot_res.get('reason','unknown')}")
        sys.exit(1)

    print(f"\n  ★ 全期間結果")
    print(f"    τ        = {boot_res['tau']:.3f} d   "
          f"[95%CI: {boot_res['tau_lo']:.2f}, {boot_res['tau_hi']:.2f}]")
    print(f"    LE_0     = {boot_res['le0']:.1f}    "
          f"[95%CI: {boot_res['le0_lo']:.1f}, {boot_res['le0_hi']:.1f}] W/m²")
    print(f"    LE_∞     = {boot_res['le_inf']:.1f}    "
          f"[95%CI: {boot_res['le_inf_lo']:.1f}, {boot_res['le_inf_hi']:.1f}] W/m²")
    print(f"    R² (point)= {boot_res['r2']:.3f}")
    print(f"    N data   = {boot_res['n_data']}")
    print(f"    ビン数    = {len(boot_res['grouped'])}")

    # BCa CI(optional)
    if args.ci_method in ("bca", "both"):
        print(f"\n  BCa CI 計算中(jackknife)...")
        bca_lo, bca_hi = bca_ci_for_tau(
            arr_d, arr_y, boot_res["boot_taus"], boot_res["tau"],
            min_per_bin=args.min_per_bin)
        if not np.isnan(bca_lo):
            print(f"    τ BCa CI = [{bca_lo:.2f}, {bca_hi:.2f}]")
            boot_res["tau_bca_lo"] = bca_lo
            boot_res["tau_bca_hi"] = bca_hi
        else:
            print(f"    BCa CI 計算失敗(jackknife サンプル不足)")

    # ── ビン詳細を表示 ──
    print(f"\n  ビン構造(透明性のため全表示):")
    grp = boot_res["grouped"]
    print("  d  |  n  | median LE")
    print("  ---|-----|---------")
    for _, row in grp.iterrows():
        print(f"  {int(row['d']):3d}| {int(row['count']):4d}| {row['median']:7.2f}")
    grp.to_csv(out_dir/"v21_bin_samples.csv", index=False)

    # ============================================================
    # 年別 τ Bootstrap
    # ============================================================
    print(f"\n{'='*60}\n年別 τ Bootstrap\n{'='*60}")
    year_results = []
    for yr in args.years:
        sub_y = sub[sub["year"]==yr]
        arr_d_y = sub_y["days_since_irrig"].values.astype(float)
        arr_y_y = sub_y["LE_corr"].values
        print(f"\n  [{yr}] N = {len(sub_y)} 日")
        if len(sub_y) < 20:
            print(f"    サンプル不足、スキップ")
            year_results.append(dict(year=yr, n_data=len(sub_y),
                                      tau=np.nan, tau_lo=np.nan, tau_hi=np.nan,
                                      r2=np.nan, le0=np.nan, le_inf=np.nan,
                                      failure_rate=np.nan,
                                      n_boot_success=0))
            continue
        res = bootstrap_recovery(arr_d_y, arr_y_y,
                                  n_boot=args.n_boot,
                                  min_per_bin=args.min_per_bin,
                                  verbose=False)
        if not res["success"]:
            print(f"    Bootstrap 失敗: {res.get('reason','unknown')}")
            year_results.append(dict(year=yr, n_data=len(sub_y),
                                      tau=res.get('tau',np.nan),
                                      tau_lo=np.nan, tau_hi=np.nan,
                                      r2=res.get('r2',np.nan),
                                      le0=res.get('le0',np.nan),
                                      le_inf=res.get('le_inf',np.nan),
                                      failure_rate=res.get('failure_rate',np.nan),
                                      n_boot_success=0))
            continue
        print(f"    τ = {res['tau']:.2f} d  [{res['tau_lo']:.2f}, {res['tau_hi']:.2f}]  "
              f"R²={res['r2']:.2f}  失敗率={res['failure_rate']*100:.1f}%")
        year_results.append(dict(
            year=yr, n_data=res["n_data"],
            tau=res["tau"], tau_lo=res["tau_lo"], tau_hi=res["tau_hi"],
            le0=res["le0"], le_inf=res["le_inf"],
            r2=res["r2"],
            failure_rate=res["failure_rate"],
            n_boot_success=res["n_boot_success"]))
    year_df = pd.DataFrame(year_results)

    # ── CSV 保存 ──
    summary_row = dict(scope="全期間",
                       n_data=boot_res["n_data"],
                       tau=boot_res["tau"],
                       tau_lo=boot_res["tau_lo"], tau_hi=boot_res["tau_hi"],
                       le0=boot_res["le0"],
                       le0_lo=boot_res["le0_lo"], le0_hi=boot_res["le0_hi"],
                       le_inf=boot_res["le_inf"],
                       le_inf_lo=boot_res["le_inf_lo"], le_inf_hi=boot_res["le_inf_hi"],
                       r2=boot_res["r2"],
                       failure_rate=boot_res["failure_rate"],
                       n_boot_success=boot_res["n_boot_success"])
    if "tau_bca_lo" in boot_res:
        summary_row["tau_bca_lo"] = boot_res["tau_bca_lo"]
        summary_row["tau_bca_hi"] = boot_res["tau_bca_hi"]

    combined_df = pd.concat([
        pd.DataFrame([summary_row]),
        year_df.assign(scope=year_df["year"].astype(str))
    ], ignore_index=True)
    combined_df.to_csv(out_dir/"v21_recovery_bootstrap.csv", index=False)

    # ── 可視化 ──
    if not args.no_plots:
        print(f"\n--- 可視化 ---")
        plot_recovery_with_ci(boot_res, out_dir)
        plot_tau_bootstrap_dist(boot_res, out_dir)
        plot_per_year_tau_ci(year_df, boot_res["tau"], out_dir)
        plot_bin_samples(boot_res["grouped"], out_dir)

    # ============================================================
    # 論文報告用フォーマット
    # ============================================================
    print(f"\n{'='*60}\n★ 論文記述用フォーマット\n{'='*60}")
    print(f"\n  全期間 τ:")
    print(f"    'τ = {boot_res['tau']:.2f} d (95% CI: {boot_res['tau_lo']:.2f}–{boot_res['tau_hi']:.2f}),")
    print(f"     fit on {len(boot_res['grouped'])} bin medians from N={boot_res['n_data']} daily")
    print(f"     observations, R² = {boot_res['r2']:.3f},")
    print(f"     bootstrap based on raw-data resampling (B = {boot_res['n_boot_total']},")
    print(f"     {boot_res['failure_rate']*100:.1f}% iteration failure)'")

    trusted = year_df[year_df["r2"] >= R2_MIN_TRUST]
    if len(trusted) >= 2:
        print(f"\n  年間変動:")
        print(f"    'Per-year τ ranged {trusted['tau'].min():.1f}–{trusted['tau'].max():.1f} d")
        print(f"     across {len(trusted)} reliable years (R² ≥ {R2_MIN_TRUST})'")

    print(f"\n[done] outputs in {out_dir}/")


if __name__ == "__main__":
    main()
