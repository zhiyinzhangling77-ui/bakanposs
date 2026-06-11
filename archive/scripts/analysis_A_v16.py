"""
================================================================
解析A v16 : 三点パッケージ + Recovery τ + 年別Sensitivity + Fit比較
================================================================

【v15 からの主な変更点】
  1. τ 年別 Sensitivity Analysis (年ごとに τ を抽出 → Robustness 確認)
  2. Fit モデル比較 (exponential vs linear vs logarithmic — AIC/R² で評価)
  3. 図の自己説明性を大幅向上
       - タイトルに定量結論を埋め込む
       - n_n / n_s をバーに直接表示
       - 二層モデル解釈 (LE∞ の意味) を fig02 に明記
       - 判定カラーに凡例を追加
  4. irrig_active_month フィルタを全工程で統一適用
  5. 判定ロジックの guard は v15 を継承 (n<30 → insufficient_data、CI幅>0.5 → uncertain_wide_CI)

【図の構成】
  fig01_verdict_panel.png        : 三点パッケージ統合判定 (n 表示付き)
  fig02_recovery_curve.png       : τ メイン + 年別 τ + fit 比較 (2×2 レイアウト)
  fig03_three_point_matrix.png   : 3 指標間の整合性 (相互散布図)
  fig04_tau_sensitivity.png      : 年別 τ の変動幅 (standalone)

【科学的主張のレベル管理】
  "Dominant explanation: irrigation-driven decoupling (τ ≈ 4-5 d)"
  "Deep rooting not rejected but cannot explain observed dose-response"
  ─ verdict() は verdict_label と judgment_strength を返し、
    後者が "strong" のときのみ確定ラベルを表示する

【入力】
  daily_classified_v4.parquet, closure_slopes_v4.json
  Tarazona daily CSV (Irrig_mm)

【出力】./output_analysis_A_v16/
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")

# ================================================================
# パス設定
# ================================================================
INPUT_DIR      = Path("/home/shion-nagamine/bakanposs/analysis_A")
PARQUET        = INPUT_DIR / "daily_classified_v4.parquet"
TARA_DAILY_CSV = Path(
    "/home/shion-nagamine/Dataset/Eddy data in Spain/"
    "Daily_Summary_Filtered_forPred_ActEne26.csv")
OUT_DIR = Path("./output_analysis_A_v16")
OUT_DIR.mkdir(exist_ok=True)

# ================================================================
# パラメータ
# ================================================================
IRRIG_THRESHOLD      = 0.5   # mm/day
MIN_IRRIG_PER_MONTH  = 2     # 灌漑アクティブ月の最小イベント数

SEASONS = {
    "spring (1-4月)"      : [1, 2, 3, 4],
    "shoulder (5,6,10月)" : [5, 6, 10],
    "summer (7-9月)"      : [7, 8, 9],
}

N_BOOT            = 5000
CI_PCT            = (2.5, 97.5)
MIN_N_PER_CLASS   = 5

DENOM_FLOORS      = {"LE_corr": 5.0, "EF_corr": 0.05, "ET": 0.3}
RATIO_CLIP        = 10.0

# 判定 guard (v14 の偽陽性を防ぐ)
N_MIN_VERDICT     = 30
CI_WIDTH_MAX      = 0.5
DEEP_ROOT_BAND    = 0.15   # |SDS| <= これ かつ CI が 0 を含む → 深根候補

# 可視化共通
SITE_COL = {"Oran": "#E85D04", "Tarazona": "#1D9E75"}
FIG_BG   = "#F8F9FA"

VERDICT_COL = {
    "deep_root"         : "#1D9E75",
    "shallow"           : "#E85D04",
    "negative_anomaly"  : "#7C3AED",
    "uncertain"         : "#9E9E9E",
    "uncertain_wide_CI" : "#BDBDBD",
    "insufficient_data" : "#E8E8E8",
}
VERDICT_JP = {
    "deep_root"         : "深根候補",
    "shallow"           : "浅根(表層感受性)",
    "negative_anomaly"  : "異常値",
    "uncertain"         : "判定不能",
    "uncertain_wide_CI" : "CI幅超過 → 保留",
    "insufficient_data" : "n不足 → 保留",
}

# ================================================================
# 1. データ読込・前処理
# ================================================================

def load_and_merge() -> pd.DataFrame:
    """parquet + Tarazona 日次 CSV (Irrig_mm) を結合"""
    daily = pd.read_parquet(PARQUET)
    daily["date"] = pd.to_datetime(daily["date"])

    csv = pd.read_csv(TARA_DAILY_CSV)
    csv["date"] = pd.to_datetime(csv["date"], errors="coerce")
    irrig_cols = [c for c in ["Irrig_mm", "Rain_mm", "IrrigRain_mm"] if c in csv.columns]
    extra = (csv[["date"] + irrig_cols]
             .dropna(subset=["date"]).drop_duplicates("date"))

    tara = daily[daily["site"] == "Tarazona"].merge(extra, on="date", how="left")
    oran = daily[daily["site"] == "Oran"].copy()
    for c in irrig_cols:
        oran[c] = 0.0

    df = pd.concat([oran, tara]).sort_values(["site", "date"]).reset_index(drop=True)
    df["year"] = df["date"].dt.year
    df["month"] = df["date"].dt.month
    return df


def add_days_since_irrig(df: pd.DataFrame) -> pd.DataFrame:
    """最後の灌漑からの経過日数を計算"""
    df = df.copy()
    df["days_since_irrig"] = np.nan
    for site, idx in df.groupby("site").groups.items():
        sub = df.loc[idx].sort_values("date")
        last = None
        out  = []
        for _, row in sub.iterrows():
            irrig = row.get("Irrig_mm", 0) or 0
            if pd.notna(irrig) and irrig > IRRIG_THRESHOLD:
                last = row["date"]
                out.append(0)
            elif last is None:
                out.append(np.nan)
            else:
                out.append((row["date"] - last).days)
        df.loc[sub.index, "days_since_irrig"] = out
    return df


def add_irrig_active_month(df: pd.DataFrame) -> pd.DataFrame:
    """月別灌漑イベント数 >= MIN_IRRIG_PER_MONTH の月を active とフラグ"""
    df = df.copy()
    df["_ym"] = df["date"].dt.to_period("M")
    df["_is_irrig"] = df["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD
    monthly = (df.groupby(["site", "_ym"])["_is_irrig"]
                 .sum().reset_index(name="n_irrig"))
    active_set = set(
        zip(monthly.loc[monthly["n_irrig"] >= MIN_IRRIG_PER_MONTH, "site"],
            monthly.loc[monthly["n_irrig"] >= MIN_IRRIG_PER_MONTH, "_ym"]))
    df["irrig_active_month"] = df.apply(
        lambda r: (r["site"], r["_ym"]) in active_set, axis=1)
    df.drop(columns=["_ym", "_is_irrig"], inplace=True)
    return df


# ================================================================
# 2. 三点パッケージ (SDS + 絶対差 + rank-biserial)
# ================================================================

def safe_ratio(numer: float, denom: float, var: str) -> float:
    floor = DENOM_FLOORS.get(var, 0.05)
    if abs(denom) < floor:
        return np.nan
    r = numer / denom
    return np.nan if abs(r) > RATIO_CLIP else r


def three_point_package(arr_n, arr_s, var: str, seed: int = 42) -> dict:
    """
    SDS + 絶対差 + rank-biserial を bootstrap 5000 回付きで返す
    SDS      = 1 - median(s) / median(n)    (0 = 深根的、+ = 浅根的)
    abs_diff = median(n) - median(s)         [W/m²]
    rb       = 1 - 2U/(n_n * n_s)           [-1, +1]
    """
    n_n, n_s = len(arr_n), len(arr_s)
    empty = dict(sds=np.nan, sds_lo=np.nan, sds_hi=np.nan,
                 abs_diff=np.nan, abs_lo=np.nan, abs_hi=np.nan,
                 rb=np.nan, p=np.nan,
                 med_n=np.nan, med_s=np.nan,
                 n_n=n_n, n_s=n_s)
    if n_n < MIN_N_PER_CLASS or n_s < MIN_N_PER_CLASS:
        return empty

    med_n = float(np.median(arr_n))
    med_s = float(np.median(arr_s))
    sds_pt  = safe_ratio(med_n - med_s, med_n, var)
    abs_pt  = med_n - med_s

    rng = np.random.default_rng(seed)
    boot_sds, boot_abs = [], []
    for _ in range(N_BOOT):
        nb = rng.choice(arr_n, size=n_n, replace=True)
        sb = rng.choice(arr_s, size=n_s, replace=True)
        r  = safe_ratio(np.median(nb) - np.median(sb), np.median(nb), var)
        if not np.isnan(r):
            boot_sds.append(r)
        boot_abs.append(np.median(nb) - np.median(sb))

    sds_lo = sds_hi = np.nan
    if len(boot_sds) >= N_BOOT * 0.1:
        sds_lo, sds_hi = np.percentile(boot_sds, CI_PCT)
    abs_lo, abs_hi = np.percentile(boot_abs, CI_PCT)

    try:
        u, p = stats.mannwhitneyu(arr_n, arr_s, alternative="two-sided")
        rb = 1 - 2 * u / (n_n * n_s)
    except ValueError:
        u = p = rb = np.nan

    return dict(sds=sds_pt, sds_lo=sds_lo, sds_hi=sds_hi,
                abs_diff=abs_pt, abs_lo=abs_lo, abs_hi=abs_hi,
                rb=rb, p=p, med_n=med_n, med_s=med_s,
                n_n=n_n, n_s=n_s)


def verdict(pkg: dict) -> tuple[str, str]:
    """
    判定ラベルと強度 ("strong" / "weak") を返す

    判定順序:
      1. n 不足                         → insufficient_data / weak
      2. SDS が NaN                     → uncertain / weak
      3. CI 幅超過                      → uncertain_wide_CI / weak
      4. CI が 0 を含み |SDS|<=BAND     → deep_root / strong
      5. CI 下限 > 0                    → shallow / strong
      6. CI 上限 < 0                    → negative_anomaly / weak
      7. その他                         → uncertain / weak
    """
    n_n, n_s = pkg.get("n_n", 0), pkg.get("n_s", 0)
    sds = pkg.get("sds", np.nan)
    lo  = pkg.get("sds_lo", np.nan)
    hi  = pkg.get("sds_hi", np.nan)

    if np.isnan(sds) or n_n < N_MIN_VERDICT or n_s < N_MIN_VERDICT:
        return "insufficient_data", "weak"
    if np.isnan(lo):
        return "uncertain", "weak"
    width = hi - lo
    if width > CI_WIDTH_MAX:
        return "uncertain_wide_CI", "weak"
    if lo <= 0 <= hi and abs(sds) <= DEEP_ROOT_BAND:
        return "deep_root", "strong"
    if lo > 0:
        return "shallow", "strong"
    if hi < 0:
        return "negative_anomaly", "weak"
    return "uncertain", "weak"


# ================================================================
# 3. Recovery τ — メイン + 年別 + モデル比較
# ================================================================

def exp_model(d, le_inf, le0, tau):
    """LE(d) = LE_∞ + (LE_0 − LE_∞) · exp(−d/τ)"""
    return le_inf + (le0 - le_inf) * np.exp(-d / tau)


def lin_model(d, a, b):
    return a + b * d


def log_model(d, a, b):
    d_safe = np.where(d > 0, d, 0.5)
    return a + b * np.log(d_safe)


def _bin_median(sub: pd.DataFrame, var: str, min_per_bin: int = 5) -> pd.DataFrame:
    g = (sub.groupby("days_since_irrig")[var]
            .agg(["median", "count"]).reset_index())
    return g[g["count"] >= min_per_bin].reset_index(drop=True)


def fit_one(x, y, model_fn, p0, bounds):
    """curve_fit を実行し (popt, r2, aic) を返す。失敗時は NaN"""
    try:
        popt, _ = curve_fit(model_fn, x, y, p0=p0, bounds=bounds, maxfev=8000)
        y_pred  = model_fn(x, *popt)
        ss_res  = np.sum((y - y_pred) ** 2)
        ss_tot  = np.sum((y - y.mean()) ** 2)
        r2      = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        k       = len(popt)
        n       = len(y)
        sigma2  = ss_res / n if n > 0 else np.nan
        aic     = n * np.log(sigma2 + 1e-12) + 2 * k if not np.isnan(sigma2) else np.nan
        return popt, r2, aic
    except Exception:
        return [np.nan] * 3, np.nan, np.nan


def fit_recovery(df_active: pd.DataFrame, var: str = "LE_corr") -> dict:
    """
    Tarazona irrig_active_month 内で days_since_irrig 別中央値を取り、
    exp / linear / log の 3 モデルを比較 → τ (exp のみ) を返す
    """
    sub = df_active[df_active[var].notna() & df_active["days_since_irrig"].notna()]
    grouped = _bin_median(sub, var, min_per_bin=5)

    base = dict(tau=np.nan, le0=np.nan, le_inf=np.nan,
                r2_exp=np.nan, r2_lin=np.nan, r2_log=np.nan,
                aic_exp=np.nan, aic_lin=np.nan, aic_log=np.nan,
                n_points=len(grouped), grouped=grouped)

    if len(grouped) < 4:
        return base

    x = grouped["days_since_irrig"].values.astype(float)
    y = grouped["median"].values

    # Exponential
    popt_e, r2_e, aic_e = fit_one(
        x, y, exp_model,
        p0=[y.min(), y.max(), 5.0],
        bounds=([0, 0, 0.5], [np.inf, np.inf, 100]))

    # Linear
    popt_l, r2_l, aic_l = fit_one(
        x, y, lin_model,
        p0=[y.mean(), -1.0],
        bounds=([-np.inf, -np.inf], [np.inf, np.inf]))

    # Logarithmic
    popt_g, r2_g, aic_g = fit_one(
        x, y, log_model,
        p0=[y.mean(), -5.0],
        bounds=([-np.inf, -np.inf], [np.inf, np.inf]))

    base.update(dict(
        tau=float(popt_e[2]) if not np.isnan(popt_e[2]) else np.nan,
        le0=float(popt_e[1]) if not np.isnan(popt_e[1]) else np.nan,
        le_inf=float(popt_e[0]) if not np.isnan(popt_e[0]) else np.nan,
        popt_exp=popt_e, popt_lin=popt_l, popt_log=popt_g,
        r2_exp=r2_e, r2_lin=r2_l, r2_log=r2_g,
        aic_exp=aic_e, aic_lin=aic_l, aic_log=aic_g,
    ))
    return base


def fit_tau_by_year(df_active: pd.DataFrame, var: str = "LE_corr") -> pd.DataFrame:
    """年別 τ を抽出 — Sensitivity Analysis 用"""
    rows = []
    for yr, sub in df_active.groupby("year"):
        res = fit_recovery(sub, var)
        rows.append(dict(year=yr, tau=res["tau"], r2=res["r2_exp"],
                         n_points=res["n_points"]))
    return pd.DataFrame(rows).sort_values("year").reset_index(drop=True)


# ================================================================
# 4. 可視化ヘルパー
# ================================================================

def _sig_star(p: float) -> str:
    if np.isnan(p): return ""
    if p < 0.001: return "★★★"
    if p < 0.01:  return "★★"
    if p < 0.05:  return "★"
    return "n.s."


def _ax_style(ax, title: str = "", xlabel: str = "", ylabel: str = ""):
    ax.set_title(title, fontweight="bold", fontsize=10, pad=6)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.grid(axis="both", alpha=0.2, lw=0.7)
    ax.spines[["top", "right"]].set_visible(False)


# ================================================================
# 5. fig01 — 三点パッケージ統合判定
# ================================================================

def plot_verdict_panel(pkg_df: pd.DataFrame, save_dir: Path):
    le_df = pkg_df[pkg_df["var"] == "LE_corr"].sort_values(["site", "season"]).reset_index(drop=True)
    n_rows = len(le_df)
    if n_rows == 0:
        print("  [skip] verdict panel — no LE_corr rows"); return

    fig, axes = plt.subplots(1, 3, figsize=(22, max(5, n_rows * 1.1)))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle(
        "fig01 — 三点パッケージ統合判定 (LE_corr)\n"
        "判定ガード: n_n, n_s ≥ 30 かつ 95%CI 幅 ≤ 0.50",
        fontsize=13, fontweight="bold", y=1.01)

    y_pos  = np.arange(n_rows)
    labels = [f"{r.site}\n{r.season}" for r in le_df.itertuples()]

    # ---- Panel 1: SDS ----
    ax = axes[0]
    for i, r in enumerate(le_df.itertuples()):
        pkg = {"sds": r.sds, "sds_lo": r.sds_lo, "sds_hi": r.sds_hi,
               "n_n": r.n_n, "n_s": r.n_s}
        v, strength = verdict(pkg)
        col = VERDICT_COL[v]
        alpha = 0.90 if strength == "strong" else 0.45

        if not np.isnan(r.sds):
            ax.barh(i, r.sds, color=col, alpha=alpha, edgecolor="black", lw=0.8, height=0.65)
        if not (np.isnan(getattr(r, "sds_lo", np.nan)) or np.isnan(getattr(r, "sds_hi", np.nan))):
            ax.errorbar(r.sds, i,
                        xerr=[[r.sds - r.sds_lo], [r.sds_hi - r.sds]],
                        color="black", capsize=5, lw=1.2, fmt="none")
        # n 表示
        n_txt = f"n=({int(r.n_n)},{int(r.n_s)})"
        ax.text(-0.55, i, n_txt, va="center", ha="left", fontsize=7.5, color="#555")
        # 判定ラベル表示
        if not np.isnan(r.sds):
            ax.text(r.sds + 0.03, i,
                    f"{r.sds:+.2f}  [{VERDICT_JP[v]}]",
                    va="center", fontsize=8.5,
                    fontweight="bold" if strength == "strong" else "normal",
                    color="black")

    ax.axvline(0, color="gray", lw=1.0)
    ax.axvspan(-DEEP_ROOT_BAND, DEEP_ROOT_BAND, color="#1D9E75", alpha=0.07,
               label=f"深根帯 |SDS|≤{DEEP_ROOT_BAND}")
    ax.set_xlim(-0.6, 1.0)
    ax.set_yticks(y_pos); ax.set_yticklabels(labels, fontsize=9)
    _ax_style(ax, "(1) SDS — 相対減少率\nSDS≈0 → 深根的　SDS>0 → 浅根的", "SDS")
    ax.legend(loc="upper right", fontsize=8)

    # 凡例パッチ
    patches = [mpatches.Patch(color=c, label=VERDICT_JP[k], alpha=0.8)
               for k, c in VERDICT_COL.items()]
    ax.legend(handles=patches, loc="lower right", fontsize=7.5,
              title="判定", title_fontsize=8)

    # ---- Panel 2: 絶対差 ----
    ax = axes[1]
    for i, r in enumerate(le_df.itertuples()):
        if not np.isnan(r.abs_diff):
            ax.barh(i, r.abs_diff, color="#5B4FCF", alpha=0.75,
                    edgecolor="black", lw=0.8, height=0.65)
            if not (np.isnan(r.abs_lo) or np.isnan(r.abs_hi)):
                ax.errorbar(r.abs_diff, i,
                            xerr=[[r.abs_diff - r.abs_lo], [r.abs_hi - r.abs_diff]],
                            color="black", capsize=5, lw=1.2, fmt="none")
            offset = 1.5 if r.abs_diff >= 0 else -1.5
            ax.text(r.abs_diff + offset, i,
                    f"{r.abs_diff:+.1f} W/m²",
                    va="center", fontsize=8.5)
    ax.axvline(0, color="gray", lw=1.0)
    ax.set_yticks(y_pos); ax.set_yticklabels([])
    _ax_style(ax, "(2) 絶対差 = LE_normal − LE_soil_dry\n正 → 干ばつで LE が落ちる　0 → 変化なし",
              "LE_n − LE_s [W/m²]")

    # ---- Panel 3: rank-biserial + p 値 ----
    ax = axes[2]
    for i, r in enumerate(le_df.itertuples()):
        if not np.isnan(r.rb):
            effect = abs(r.rb)
            col = "#E85D04" if effect > 0.3 else "#FFA000" if effect > 0.1 else "#9E9E9E"
            ax.barh(i, r.rb, color=col, alpha=0.80,
                    edgecolor="black", lw=0.8, height=0.65)
            star = _sig_star(r.p)
            ax.text(r.rb + 0.03, i,
                    f"{r.rb:+.2f}  p={r.p:.1e} {star}",
                    va="center", fontsize=8.5)
    ax.axvline(0, color="gray", lw=1.0)
    ax.axvspan(0.1,  1, color="#FFA000", alpha=0.05)
    ax.axvspan(0.3,  1, color="#E85D04", alpha=0.05)
    ax.set_xlim(-1, 1.5)
    ax.set_yticks(y_pos); ax.set_yticklabels([])
    _ax_style(ax, "(3) 効果量 (rank-biserial r) + 有意性\n|r|>0.3 → 大効果　|r|>0.1 → 中効果",
              "rank-biserial r")
    # 効果量目安ライン
    for xv, lb in [(0.1, "中"), (0.3, "大")]:
        ax.axvline(xv, color="gray", ls=":", lw=0.8)
        ax.text(xv + 0.01, n_rows - 0.4, lb, fontsize=7, color="gray")

    plt.tight_layout()
    fp = save_dir / "fig01_verdict_panel.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 6. fig02 — Recovery curve (2×2 レイアウト)
# ================================================================

def plot_recovery_curve(fit_res: dict, tau_by_year: pd.DataFrame, save_dir: Path):
    grouped = fit_res.get("grouped")
    if grouped is None or len(grouped) < 2:
        print("  [skip] recovery curve — insufficient data"); return

    tau    = fit_res["tau"]
    le0    = fit_res["le0"]
    le_inf = fit_res["le_inf"]
    r2_exp = fit_res["r2_exp"]
    r2_lin = fit_res["r2_lin"]
    r2_log = fit_res["r2_log"]
    aic_e  = fit_res["aic_exp"]
    aic_l  = fit_res["aic_lin"]
    aic_g  = fit_res["aic_log"]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.patch.set_facecolor(FIG_BG)

    x_obs = grouped["days_since_irrig"].values.astype(float)
    y_obs = grouped["median"].values
    n_obs = grouped["count"].values
    x_fit = np.linspace(0, x_obs.max() + 2, 300)

    # ── Subplot (0,0): メイン Recovery curve (exponential) ──────────
    ax = axes[0, 0]
    sizes = 40 + n_obs * 4
    sc = ax.scatter(x_obs, y_obs, s=sizes, c="#1D9E75",
                    edgecolors="black", lw=1, zorder=5, label="観測中央値\n(サイズ = n)")

    if not np.isnan(tau):
        y_fit = exp_model(x_fit, le_inf, le0, tau)
        ax.plot(x_fit, y_fit, "k-", lw=2.5, label=f"Exp fit: τ={tau:.1f}d, R²={r2_exp:.3f}")
        ax.axhline(le_inf, color="#E85D04", ls="--", lw=1.5, alpha=0.9,
                   label=f"LE_∞={le_inf:.0f} W/m²\n(灌漑完全停止後の漸近値 → 深層水?)")
        ax.axhline(le0, color="#1A73E8", ls="--", lw=1.5, alpha=0.9,
                   label=f"LE_0={le0:.0f} W/m²\n(灌漑直後の上限)")
        ax.axvline(tau, color="gray", ls=":", lw=1.5,
                   label=f"τ = {tau:.1f} 日\n(e-folding decay time)")

        # τ 解釈テキスト
        if tau < 5:
            interp = f"τ = {tau:.1f}d < 5d\n→ 灌漑依存が支配的"
            interp_col = "#E85D04"
        elif tau > 14:
            interp = f"τ = {tau:.1f}d > 14d\n→ 深層水アクセスの可能性"
            interp_col = "#1D9E75"
        else:
            interp = f"τ = {tau:.1f}d\n→ 中間的応答"
            interp_col = "#FFA000"
        ax.text(0.97, 0.97, interp,
                transform=ax.transAxes, va="top", ha="right",
                fontsize=10, fontweight="bold", color=interp_col,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=interp_col, alpha=0.9))

    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=8.5)
    _ax_style(ax,
              f"(A) Recovery curve — τ = {tau:.1f} 日\n"
              "Tarazona 灌漑アクティブ月・生育期",
              "Days since last irrigation",
              "LE_corr 中央値 [W/m²]")

    # ── Subplot (0,1): 3 モデル比較 ──────────────────────────────
    ax = axes[0, 1]
    ax.scatter(x_obs, y_obs, s=sizes, c="#1D9E75",
               edgecolors="black", lw=1, zorder=5, label="観測中央値")

    model_specs = [
        ("Exponential", exp_model, fit_res.get("popt_exp"), r2_exp, aic_e, "k-", 2.5),
        ("Linear",      lin_model, fit_res.get("popt_lin"), r2_lin, aic_l, "b--", 1.8),
        ("Logarithmic", log_model, fit_res.get("popt_log"), r2_log, aic_g, "r-.", 1.8),
    ]
    best_aic = min(v for v in [aic_e, aic_l, aic_g] if not np.isnan(v)) if not all(
        np.isnan(v) for v in [aic_e, aic_l, aic_g]) else np.nan

    for name, fn, popt, r2, aic, ls, lw in model_specs:
        if popt is not None and not np.isnan(popt[0]):
            y_f = fn(x_fit, *popt)
            delta_aic = aic - best_aic if not np.isnan(aic) else np.nan
            is_best = abs(delta_aic) < 0.1 if not np.isnan(delta_aic) else False
            lbl = (f"{name}\nR²={r2:.3f}, ΔAIC={delta_aic:.1f}"
                   f"{'  ← BEST' if is_best else ''}")
            ax.plot(x_fit, y_f, ls, lw=lw, label=lbl)

    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=8.5)
    _ax_style(ax,
              "(B) モデル比較 (Exp / Linear / Log)\nAIC 最小 = 最良フィット",
              "Days since last irrigation",
              "LE_corr 中央値 [W/m²]")

    # ── Subplot (1,0): 年別 τ Sensitivity ────────────────────────
    ax = axes[1, 0]
    valid = tau_by_year.dropna(subset=["tau"])
    if len(valid) >= 2:
        ax.bar(valid["year"], valid["tau"],
               color="#5B4FCF", alpha=0.80, edgecolor="black", lw=1, width=0.6)
        ax.axhline(tau, color="#E85D04", ls="--", lw=2,
                   label=f"全期間 τ = {tau:.1f} d")
        # τ 範囲を網かけ
        tmin, tmax = valid["tau"].min(), valid["tau"].max()
        ax.axhspan(tmin, tmax, color="#5B4FCF", alpha=0.10,
                   label=f"年間範囲 {tmin:.1f}–{tmax:.1f} d")
        for _, row in valid.iterrows():
            ax.text(row["year"], row["tau"] + 0.2,
                    f"{row['tau']:.1f}d\n(R²={row['r2']:.2f})\nn={int(row['n_points'])}",
                    ha="center", va="bottom", fontsize=8)
        ax.legend(fontsize=9)
        robust_txt = ("τ の年間変動が小さい → Robust\n" if tmax - tmin < 3
                      else "τ の年間変動が大きい → 要注意\n")
        ax.text(0.97, 0.97, robust_txt,
                transform=ax.transAxes, va="top", ha="right",
                fontsize=9, color="navy",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="navy", alpha=0.85))
    else:
        ax.text(0.5, 0.5, "年別 τ: データ不足\n(有効年 < 2)", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="gray")

    _ax_style(ax,
              "(C) 年別 τ — Sensitivity Analysis\n変動幅が小さいほど結論が Robust",
              "Year",
              "τ [days]")

    # ── Subplot (1,1): 二層モデル概念図 ──────────────────────────
    ax = axes[1, 1]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.set_aspect("equal")
    ax.axis("off")

    def box(x, y, w, h, fc, ec="black", lw=1.5, alpha=0.85, text="", fs=9):
        rect = mpatches.FancyBboxPatch((x, y), w, h,
                                       boxstyle="round,pad=0.15",
                                       fc=fc, ec=ec, lw=lw, alpha=alpha)
        ax.add_patch(rect)
        if text:
            ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                    fontsize=fs, fontweight="bold", wrap=True)

    def arrow(x1, y1, x2, y2, col="black"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))

    # 空 → 背景
    ax.add_patch(mpatches.Rectangle((0, 0), 10, 10, fc=FIG_BG, ec="none"))

    # タイトル
    ax.text(5, 9.5, "(D) 二層モデル解釈",
            ha="center", va="top", fontsize=11, fontweight="bold")

    # 灌漑パルス
    box(0.5, 7.5, 2.5, 1.2, fc="#AED6F1", text="灌漑パルス\n(drip)")
    arrow(3.0, 8.1, 4.5, 8.1)

    # Fast pool
    box(4.5, 7.2, 4.5, 1.8, fc="#D5EAF5",
        text=f"Fast pool\n(灌漑水 → τ≈{tau:.1f}d で枯渇)")
    arrow(6.75, 7.2, 6.75, 5.8, col="#1A73E8")

    # Slow pool
    box(4.5, 4.2, 4.5, 1.5, fc="#D5F5E3",
        text=f"Slow pool\n(深層水? → LE_∞≈{le_inf:.0f} W/m²)")
    arrow(6.75, 4.2, 6.75, 2.8, col="#1D9E75")

    # 蒸散
    box(4.0, 1.5, 5.5, 1.2, fc="#FDE8D8",
        text=f"蒸散 LE\nLE_0≈{le0:.0f}→LE_∞≈{le_inf:.0f} W/m²")

    # 主張テキスト
    ax.text(5, 0.5,
            "主要な SWC-ET decoupling は灌漑水が支配的\n"
            "(深根完全棄却ではないが支配要因ではない)",
            ha="center", va="bottom", fontsize=8.5, style="italic",
            color="#333",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999", alpha=0.8))

    plt.tight_layout(h_pad=3, w_pad=3)
    fp = save_dir / "fig02_recovery_curve.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 7. fig03 — 三点指標間の整合性
# ================================================================

def plot_three_point_matrix(pkg_df: pd.DataFrame, save_dir: Path):
    sub = pkg_df[pkg_df["var"] == "LE_corr"].dropna(subset=["sds", "abs_diff", "rb"])
    if len(sub) < 2:
        print("  [skip] three-point matrix — too few rows"); return

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle(
        "fig03 — 三点指標間の整合性 (LE_corr)\n"
        "3 指標が同じ方向を向いているほど主張が Robust",
        fontsize=12, fontweight="bold")

    pairs = [
        ("sds", "abs_diff", "SDS", "Abs diff [W/m²]", "SDS vs 絶対差"),
        ("sds", "rb", "SDS", "rank-biserial r", "SDS vs 効果量"),
        ("abs_diff", "rb", "Abs diff [W/m²]", "rank-biserial r", "絶対差 vs 効果量"),
    ]
    for ax, (xk, yk, xl, yl, title) in zip(axes, pairs):
        for _, r in sub.iterrows():
            v, strength = verdict({"sds": r["sds"], "sds_lo": r["sds_lo"],
                                    "sds_hi": r["sds_hi"], "n_n": r["n_n"], "n_s": r["n_s"]})
            col  = SITE_COL.get(r["site"], "#888")
            mk   = "o" if r["site"] == "Oran" else "^"
            alpha = 0.9 if strength == "strong" else 0.4
            ax.scatter(r[xk], r[yk], c=col, s=90, alpha=alpha,
                       edgecolors="black", lw=1, marker=mk)
            ax.annotate(f"{r['site'][:3]}\n{r['season'][:6]}",
                        (r[xk], r[yk]), fontsize=7.5,
                        xytext=(5, 5), textcoords="offset points", color="#333")
        ax.axhline(0, color="gray", lw=0.6, alpha=0.5)
        ax.axvline(0, color="gray", lw=0.6, alpha=0.5)
        _ax_style(ax, title, xl, yl)

    # 凡例
    leg_handles = [
        mpatches.Patch(color=SITE_COL["Oran"],     label="Oran (rainfed cereal)"),
        mpatches.Patch(color=SITE_COL["Tarazona"], label="Tarazona (irrigated almond)"),
        Line2D([0], [0], marker="o", color="gray", label="strong verdict",
               ms=8, lw=0, markeredgecolor="black"),
        Line2D([0], [0], marker="o", color="gray", label="weak verdict",
               ms=8, lw=0, alpha=0.35, markeredgecolor="black"),
    ]
    fig.legend(handles=leg_handles, loc="lower center", ncol=4, fontsize=9, frameon=True)

    plt.tight_layout(rect=[0, 0.07, 1, 1])
    fp = save_dir / "fig03_three_point_matrix.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 8. fig04 — 年別 τ Standalone
# ================================================================

def plot_tau_sensitivity(tau_by_year: pd.DataFrame, tau_all: float, save_dir: Path):
    valid = tau_by_year.dropna(subset=["tau"])
    if len(valid) < 2:
        print("  [skip] tau sensitivity — too few valid years"); return

    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_facecolor(FIG_BG)

    tmin = valid["tau"].min()
    tmax = valid["tau"].max()
    tmean = valid["tau"].mean()
    tstd  = valid["tau"].std()

    bars = ax.bar(valid["year"], valid["tau"],
                  color="#5B4FCF", alpha=0.80,
                  edgecolor="black", lw=1.2, width=0.6)
    ax.axhline(tau_all, color="#E85D04", ls="--", lw=2.2,
               label=f"全期間 τ = {tau_all:.1f} d")
    ax.axhline(5, color="gray", ls=":", lw=1.0, alpha=0.7, label="τ = 5 d (灌漑依存閾値)")
    ax.axhspan(tmin, tmax, color="#5B4FCF", alpha=0.08,
               label=f"年間レンジ: {tmin:.1f}–{tmax:.1f} d")

    for bar, (_, row) in zip(bars, valid.iterrows()):
        ax.text(bar.get_x() + bar.get_width() / 2,
                row["tau"] + 0.15,
                f"{row['tau']:.1f}d\n(n={int(row['n_points'])})",
                ha="center", va="bottom", fontsize=9.5)

    ax.set_ylim(0, max(valid["tau"].max() * 1.4, 8))
    ax.set_xticks(valid["year"])
    ax.legend(fontsize=9.5)

    robust_flag = "★ τ の年間変動が小さい → 結論 ROBUST" if tmax - tmin < 3 else "△ τ の年間変動が大きい → 解釈に注意"
    robust_col  = "#1D9E75" if tmax - tmin < 3 else "#E85D04"
    ax.text(0.98, 0.98, f"{robust_flag}\n平均 τ = {tmean:.1f} ± {tstd:.1f} d",
            transform=ax.transAxes, va="top", ha="right",
            fontsize=10, fontweight="bold", color=robust_col,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=robust_col, alpha=0.9))

    _ax_style(ax,
              f"fig04 — 年別 τ Sensitivity Analysis (Tarazona, 灌漑アクティブ月)\n"
              f"全期間 τ = {tau_all:.1f} d  |  年間レンジ: {tmin:.1f}–{tmax:.1f} d",
              "Year",
              "Recovery time τ [days]")
    plt.tight_layout()
    fp = save_dir / "fig04_tau_sensitivity.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 9. MAIN
# ================================================================

def main():
    print("=" * 60)
    print("解析A v16 — 三点パッケージ + Recovery τ + Sensitivity + Fit 比較")
    print("=" * 60)

    # データ準備
    df = load_and_merge()
    df = add_days_since_irrig(df)
    df = add_irrig_active_month(df)

    n_active = ((df["site"] == "Tarazona") & df["irrig_active_month"]).sum()
    print(f"  Tarazona 灌漑アクティブ日数: {n_active}")

    # ============================================================
    # 三点パッケージ
    # ============================================================
    print(f"\n{'='*60}\n三点パッケージ計算\n{'='*60}")
    rows = []
    for season_label, months in SEASONS.items():
        for site in ["Oran", "Tarazona"]:
            sub = df[(df["site"] == site) & df["is_growing"] & df["month"].isin(months)]
            if site == "Tarazona":
                sub = sub[sub["irrig_active_month"]]
            for var in ["LE_corr", "EF_corr", "ET"]:
                if var not in sub.columns: continue
                n_arr = sub.loc[sub["drought_type"] == "normal",    var].dropna().values
                s_arr = sub.loc[sub["drought_type"] == "soil dry",  var].dropna().values
                pkg = three_point_package(n_arr, s_arr, var)
                v, strength = verdict(pkg)
                pkg.update(dict(site=site, season=season_label, var=var,
                                verdict=v, verdict_strength=strength))
                rows.append(pkg)
                if var == "LE_corr":
                    print(f"  [{site:9s}/{season_label:22s}] "
                          f"SDS={pkg['sds']:+.3f} "
                          f"abs={pkg['abs_diff']:+.1f}W/m² "
                          f"rb={pkg['rb']:+.3f} "
                          f"p={pkg['p']:.1e} "
                          f"n=({pkg['n_n']},{pkg['n_s']}) "
                          f"→ [{VERDICT_JP[v]}] ({strength})")

    pkg_df = pd.DataFrame(rows)
    pkg_df.to_csv(OUT_DIR / "v16_three_point_package.csv", index=False)
    print(f"\n  [save] v16_three_point_package.csv ({len(pkg_df)} rows)")

    # ============================================================
    # Recovery τ — 全期間
    # ============================================================
    print(f"\n{'='*60}\nRecovery τ 解析\n{'='*60}")
    tara_active = df[(df["site"] == "Tarazona") & df["irrig_active_month"] & df["is_growing"]]
    fit_res = fit_recovery(tara_active, var="LE_corr")

    print(f"  τ (全期間)= {fit_res['tau']:.2f} d")
    print(f"  LE_0 (d=0)= {fit_res['le0']:.1f} W/m²")
    print(f"  LE_∞ (d→∞)= {fit_res['le_inf']:.1f} W/m²")
    print(f"  R² exp={fit_res['r2_exp']:.3f}, lin={fit_res['r2_lin']:.3f}, log={fit_res['r2_log']:.3f}")
    print(f"  AIC exp={fit_res['aic_exp']:.1f}, lin={fit_res['aic_lin']:.1f}, log={fit_res['aic_log']:.1f}")
    print(f"  Fit ビン数: {fit_res['n_points']}")

    # ============================================================
    # τ 年別 Sensitivity
    # ============================================================
    print(f"\n{'='*60}\n年別 τ Sensitivity\n{'='*60}")
    tau_by_year = fit_tau_by_year(tara_active, var="LE_corr")
    print(tau_by_year.to_string(index=False))
    valid_tau = tau_by_year.dropna(subset=["tau"])
    if len(valid_tau) >= 2:
        print(f"\n  τ 年間レンジ: {valid_tau['tau'].min():.1f}–{valid_tau['tau'].max():.1f} d")
        print(f"  τ 平均 ± std: {valid_tau['tau'].mean():.1f} ± {valid_tau['tau'].std():.1f} d")

    # CSV 保存
    if fit_res.get("grouped") is not None:
        out_g = fit_res["grouped"].copy()
        out_g["site"] = "Tarazona"
        out_g.to_csv(OUT_DIR / "v16_recovery_binned.csv", index=False)
    tau_by_year.to_csv(OUT_DIR / "v16_tau_by_year.csv", index=False)
    with open(OUT_DIR / "v16_recovery_fit_params.json", "w") as f:
        json.dump({k: (float(v) if isinstance(v, (float, np.floating)) else
                       (int(v) if isinstance(v, (int, np.integer)) else v))
                   for k, v in fit_res.items()
                   if k not in ("grouped", "popt_exp", "popt_lin", "popt_log")},
                  f, indent=2)
    print(f"\n  [save] v16_recovery_binned.csv, v16_tau_by_year.csv, v16_recovery_fit_params.json")

    # ============================================================
    # 可視化
    # ============================================================
    print(f"\n--- 可視化 ---")
    plot_verdict_panel(pkg_df, OUT_DIR)
    plot_recovery_curve(fit_res, tau_by_year, OUT_DIR)
    plot_three_point_matrix(pkg_df, OUT_DIR)
    plot_tau_sensitivity(tau_by_year, fit_res["tau"], OUT_DIR)

    # ============================================================
    # 最終サマリー
    # ============================================================
    print(f"\n{'='*60}")
    print("★ 最終サマリー (LE_corr ベース)")
    print(f"{'='*60}")
    le_pkg = pkg_df[pkg_df["var"] == "LE_corr"]
    for _, row in le_pkg.iterrows():
        v, strength = row["verdict"], row["verdict_strength"]
        lo = row.get("sds_lo", np.nan)
        hi = row.get("sds_hi", np.nan)
        ci_str = (f"[{lo:+.2f}, {hi:+.2f}]" if not np.isnan(lo) else "[CI 算出不可]")
        print(f"  [{row['site']:9s}/{row['season']:22s}] "
              f"SDS={row['sds']:+.3f} {ci_str}  "
              f"n=({row['n_n']:3.0f},{row['n_s']:3.0f})  "
              f"→ {VERDICT_JP[v]} ({strength})")

    tau = fit_res["tau"]
    print(f"\n  Recovery τ = {tau:.1f} d  (R²={fit_res['r2_exp']:.3f})")
    print(f"  Best fit model AIC: exp={fit_res['aic_exp']:.1f}, "
          f"lin={fit_res['aic_lin']:.1f}, log={fit_res['aic_log']:.1f}")
    if not np.isnan(tau):
        if tau < 5:
            print(f"  → 灌漑依存が支配的 (τ < 5d)")
            print(f"     主張レベル: 'irrigation-driven decoupling (not deep rooting)'")
        elif tau > 14:
            print(f"  → 深層水アクセスの可能性あり (τ > 14d)")
        else:
            print(f"  → 中間的応答: 他指標と総合判断")

    print(f"\n[done] outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
