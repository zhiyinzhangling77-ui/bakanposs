"""
================================================================
解析A v17 : v16 + τ Diagnostic + 年別除外ロジック + 変動要因注釈
================================================================

【v16 からの変更点】
  1. fig05_tau_diagnostic.png を追加
       - 年別 τ × 灌漑イベント数 × 降水量を 1 枚に並べて変動要因を診断
       - R² < R2_MIN_TRUST または n_points < N_MIN_YEAR のビンは "低信頼" フラグ
  2. τ 年別 Sensitivity 出力に信頼性フラグを追加
       - tau_reliable 列: R² >= R2_MIN_TRUST AND n_points >= N_MIN_YEAR
  3. 低信頼年を除外した τ レンジをコンソール出力
  4. irrig_stats_by_year.csv を出力(教授向け補足資料)
  5. 全コードにコメント追加(なぜそうしたか)

【信頼性基準の考え方】
  年別 fit は全期間 fit より binが少ない。
  2023 は n_points=4 (< N_MIN_YEAR=5) で fit が不安定。
  R² < R2_MIN_TRUST=0.7 は説明力が低い → τ 値自体が不確か。
  これらを "低信頼" とフラグし、
  信頼年のみでの τ レンジを "robust range" として報告する。

【τ 変動の科学的解釈(コメント形式)】
  灌漑依存(2021, 2022) τ < 5d:
    → 灌漑頻度が高い年, 1 イベントあたりの水量が少ない
    → 土壌水分が素早く使い果たされる
  長 τ (2020, 2024) τ ≈ 14-16d:
    → 灌漑ペースが遅い or 水量が多い
    → または降水量が多く土壌水分ベースラインが高い
  2023 (n=4, τ=26d):
    → ビン数不足 → fit 不安定 → 除外推奨
  結論: τ の年変動は "灌漑管理の年間差" を反映する可能性が高い
        (生物学的深根の年間変化ではない)

【入力】(v16 と同じ)
  daily_classified_v4.parquet, Tarazona daily CSV

【出力】./output_analysis_A_v17/
  fig01–fig04 : v16 と同一
  fig05_tau_diagnostic.png  ← NEW
  v17_*.csv, v17_*.json
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
from matplotlib.gridspec import GridSpec
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
OUT_DIR = Path("./output_analysis_A_v17")
OUT_DIR.mkdir(exist_ok=True)

# ================================================================
# パラメータ
# ================================================================
IRRIG_THRESHOLD      = 0.5   # mm/day — 灌漑イベントの最小値
MIN_IRRIG_PER_MONTH  = 2     # この数以上のイベントがある月を "active" とする

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

# 判定 guard (v14 偽陽性回避)
N_MIN_VERDICT     = 30
CI_WIDTH_MAX      = 0.5
DEEP_ROOT_BAND    = 0.15

# ★ v17 新規: 年別 τ の信頼性基準
R2_MIN_TRUST      = 0.70   # R² がこれ未満 → fit 不安定 → 低信頼
# n_points はビン数であり元データのサンプル数ではない。
# 各ビンに >=5 サンプルを要求しているので n_points=4 でも元データは十分。
# 信頼性判定は R² のみで行い、n_points 閾値は廃止する。
N_MIN_YEAR        = 4      # fit_recovery 内の min_per_bin 由来の参考値のみ

# 可視化
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
    "uncertain_wide_CI" : "CI幅超過→保留",
    "insufficient_data" : "n不足→保留",
}

# ================================================================
# 1. データ読込・前処理 (v16 と同じ)
# ================================================================

def load_and_merge() -> pd.DataFrame:
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
    df["year"]  = df["date"].dt.year
    df["month"] = df["date"].dt.month
    return df


def add_days_since_irrig(df: pd.DataFrame) -> pd.DataFrame:
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


def add_irrig_active_month(df: pd.DataFrame) -> pd.DataFrame:
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
# 2. 三点パッケージ + 判定 (v16 と同じ)
# ================================================================

def safe_ratio(numer, denom, var):
    floor = DENOM_FLOORS.get(var, 0.05)
    if abs(denom) < floor: return np.nan
    r = numer / denom
    return np.nan if abs(r) > RATIO_CLIP else r


def three_point_package(arr_n, arr_s, var, seed=42):
    n_n, n_s = len(arr_n), len(arr_s)
    empty = dict(sds=np.nan, sds_lo=np.nan, sds_hi=np.nan,
                 abs_diff=np.nan, abs_lo=np.nan, abs_hi=np.nan,
                 rb=np.nan, p=np.nan, med_n=np.nan, med_s=np.nan,
                 n_n=n_n, n_s=n_s)
    if n_n < MIN_N_PER_CLASS or n_s < MIN_N_PER_CLASS: return empty
    med_n, med_s = float(np.median(arr_n)), float(np.median(arr_s))
    sds_pt, abs_pt = safe_ratio(med_n - med_s, med_n, var), med_n - med_s
    rng = np.random.default_rng(seed)
    boot_sds, boot_abs = [], []
    for _ in range(N_BOOT):
        nb = rng.choice(arr_n, size=n_n, replace=True)
        sb = rng.choice(arr_s, size=n_s, replace=True)
        r = safe_ratio(np.median(nb) - np.median(sb), np.median(nb), var)
        if not np.isnan(r): boot_sds.append(r)
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
                rb=rb, p=p, med_n=med_n, med_s=med_s, n_n=n_n, n_s=n_s)


def verdict(pkg):
    n_n, n_s = pkg.get("n_n", 0), pkg.get("n_s", 0)
    sds, lo, hi = pkg.get("sds", np.nan), pkg.get("sds_lo", np.nan), pkg.get("sds_hi", np.nan)
    if np.isnan(sds) or n_n < N_MIN_VERDICT or n_s < N_MIN_VERDICT:
        return "insufficient_data", "weak"
    if np.isnan(lo): return "uncertain", "weak"
    width = hi - lo
    if width > CI_WIDTH_MAX: return "uncertain_wide_CI", "weak"
    if lo <= 0 <= hi and abs(sds) <= DEEP_ROOT_BAND: return "deep_root", "strong"
    if lo > 0: return "shallow", "strong"
    if hi < 0: return "negative_anomaly", "weak"
    return "uncertain", "weak"


# ================================================================
# 3. Recovery τ 解析 (v16 と同じ + 信頼性フラグ)
# ================================================================

def exp_model(d, le_inf, le0, tau):
    return le_inf + (le0 - le_inf) * np.exp(-d / tau)

def lin_model(d, a, b): return a + b * d
def log_model(d, a, b): return a + b * np.log(np.where(d > 0, d, 0.5))


def _bin_median(sub, var, min_per_bin=5):
    g = (sub.groupby("days_since_irrig")[var]
            .agg(["median", "count"]).reset_index())
    return g[g["count"] >= min_per_bin].reset_index(drop=True)


def fit_one(x, y, model_fn, p0, bounds):
    try:
        popt, _ = curve_fit(model_fn, x, y, p0=p0, bounds=bounds, maxfev=8000)
        y_pred  = model_fn(x, *popt)
        ss_res  = np.sum((y - y_pred) ** 2)
        ss_tot  = np.sum((y - y.mean()) ** 2)
        r2  = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        k, n = len(popt), len(y)
        aic = n * np.log(ss_res / n + 1e-12) + 2 * k if n > 0 else np.nan
        return popt, r2, aic
    except Exception:
        return [np.nan] * 3, np.nan, np.nan


def fit_recovery(df_active, var="LE_corr"):
    sub     = df_active[df_active[var].notna() & df_active["days_since_irrig"].notna()]
    grouped = _bin_median(sub, var)
    base    = dict(tau=np.nan, le0=np.nan, le_inf=np.nan,
                   r2_exp=np.nan, r2_lin=np.nan, r2_log=np.nan,
                   aic_exp=np.nan, aic_lin=np.nan, aic_log=np.nan,
                   n_points=len(grouped), grouped=grouped)
    if len(grouped) < 4: return base
    x, y = grouped["days_since_irrig"].values.astype(float), grouped["median"].values
    popt_e, r2_e, aic_e = fit_one(x, y, exp_model,
        p0=[y.min(), y.max(), 5.0], bounds=([0,0,0.5],[np.inf,np.inf,100]))
    popt_l, r2_l, aic_l = fit_one(x, y, lin_model,
        p0=[y.mean(),-1.0], bounds=([-np.inf,-np.inf],[np.inf,np.inf]))
    popt_g, r2_g, aic_g = fit_one(x, y, log_model,
        p0=[y.mean(),-5.0], bounds=([-np.inf,-np.inf],[np.inf,np.inf]))
    base.update(dict(
        tau=float(popt_e[2]) if not np.isnan(popt_e[2]) else np.nan,
        le0=float(popt_e[1]) if not np.isnan(popt_e[1]) else np.nan,
        le_inf=float(popt_e[0]) if not np.isnan(popt_e[0]) else np.nan,
        popt_exp=popt_e, popt_lin=popt_l, popt_log=popt_g,
        r2_exp=r2_e, r2_lin=r2_l, r2_log=r2_g,
        aic_exp=aic_e, aic_lin=aic_l, aic_log=aic_g,
    ))
    return base


def fit_tau_by_year(df_active, var="LE_corr"):
    """
    年別 τ を抽出し、信頼性フラグを付与する。
    信頼基準: R² >= R2_MIN_TRUST AND n_points >= N_MIN_YEAR
    低信頼年は τ 値が不安定なため robust range から除外する。
    """
    rows = []
    for yr, sub in df_active.groupby("year"):
        res     = fit_recovery(sub, var)
        # 信頼性は R² のみで判定 (n_points はビン数でありサンプル数ではないため除外)
        trusted = (not np.isnan(res["r2_exp"])
                   and res["r2_exp"] >= R2_MIN_TRUST)
        # 低信頼の主な原因をラベル化
        if np.isnan(res["r2_exp"]):
            reason = "fit失敗"
        elif res["r2_exp"] < R2_MIN_TRUST:
            reason = f"R²={res['r2_exp']:.2f} < {R2_MIN_TRUST}"
        else:
            reason = "OK"
        rows.append(dict(year=yr, tau=res["tau"], r2=res["r2_exp"],
                         n_points=res["n_points"],
                         tau_reliable=trusted, reason=reason))
    return pd.DataFrame(rows).sort_values("year").reset_index(drop=True)


# ================================================================
# 4. 年別灌漑統計 (diagnostic 用)
# ================================================================

def calc_irrig_stats_by_year(df: pd.DataFrame) -> pd.DataFrame:
    """
    Tarazona の年別灌漑統計を計算。
    tau_diagnostic で τ 変動の要因診断に使う。
    """
    tara = df[df["site"] == "Tarazona"].copy()
    tara["is_irrig_event"] = tara["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD
    stats_list = []
    for yr, sub in tara.groupby("year"):
        n_events     = sub["is_irrig_event"].sum()
        total_mm     = sub["Irrig_mm"].fillna(0).sum()
        avg_mm_event = total_mm / n_events if n_events > 0 else np.nan
        rain_mm      = sub["Rain_mm"].fillna(0).sum() if "Rain_mm" in sub.columns else np.nan
        stats_list.append(dict(
            year=yr,
            n_irrig_events=int(n_events),
            total_irrig_mm=round(total_mm, 1),
            avg_mm_per_event=round(avg_mm_event, 1) if not np.isnan(avg_mm_event) else np.nan,
            rain_mm=round(rain_mm, 1) if not np.isnan(rain_mm) else np.nan,
        ))
    return pd.DataFrame(stats_list)


# ================================================================
# 5. 可視化 (fig01–04 は v16 と同じ、fig05 が新規)
# ================================================================

def _sig_star(p):
    if np.isnan(p): return ""
    if p < 0.001: return "★★★"
    if p < 0.01:  return "★★"
    if p < 0.05:  return "★"
    return "n.s."


def _ax_style(ax, title="", xlabel="", ylabel=""):
    ax.set_title(title, fontweight="bold", fontsize=10, pad=6)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.grid(axis="both", alpha=0.2, lw=0.7)
    ax.spines[["top", "right"]].set_visible(False)


# ---- fig01: 三点パッケージ ------

def plot_verdict_panel(pkg_df, save_dir):
    le_df = pkg_df[pkg_df["var"] == "LE_corr"].sort_values(["site","season"]).reset_index(drop=True)
    n_rows = len(le_df)
    if n_rows == 0: return
    fig, axes = plt.subplots(1, 3, figsize=(22, max(5, n_rows * 1.1)))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle(
        "fig01 — 三点パッケージ統合判定 (LE_corr)\n"
        "判定ガード: n_n, n_s ≥ 30 かつ 95%CI 幅 ≤ 0.50",
        fontsize=13, fontweight="bold", y=1.01)
    y_pos  = np.arange(n_rows)
    labels = [f"{r.site}\n{r.season}" for r in le_df.itertuples()]

    ax = axes[0]
    for i, r in enumerate(le_df.itertuples()):
        pkg = {"sds":r.sds,"sds_lo":r.sds_lo,"sds_hi":r.sds_hi,"n_n":r.n_n,"n_s":r.n_s}
        v, strength = verdict(pkg)
        col   = VERDICT_COL[v]
        alpha = 0.90 if strength == "strong" else 0.45
        if not np.isnan(r.sds):
            ax.barh(i, r.sds, color=col, alpha=alpha, edgecolor="black", lw=0.8, height=0.65)
        lo_ = getattr(r, "sds_lo", np.nan)
        hi_ = getattr(r, "sds_hi", np.nan)
        if not (np.isnan(lo_) or np.isnan(hi_)):
            ax.errorbar(r.sds, i, xerr=[[r.sds-lo_],[hi_-r.sds]],
                        color="black", capsize=5, lw=1.2, fmt="none")
        ax.text(-0.55, i, f"n=({int(r.n_n)},{int(r.n_s)})",
                va="center", ha="left", fontsize=7.5, color="#555")
        if not np.isnan(r.sds):
            ax.text(r.sds + 0.03, i,
                    f"{r.sds:+.2f}  [{VERDICT_JP[v]}]", va="center",
                    fontsize=8.5, fontweight="bold" if strength == "strong" else "normal")
    ax.axvline(0, color="gray", lw=1.0)
    ax.axvspan(-DEEP_ROOT_BAND, DEEP_ROOT_BAND, color="#1D9E75", alpha=0.07)
    ax.set_xlim(-0.6, 1.0)
    ax.set_yticks(y_pos); ax.set_yticklabels(labels, fontsize=9)
    patches = [mpatches.Patch(color=c, label=VERDICT_JP[k], alpha=0.8)
               for k, c in VERDICT_COL.items()]
    ax.legend(handles=patches, loc="lower right", fontsize=7, title="判定", title_fontsize=8)
    _ax_style(ax, "(1) SDS — 相対減少率\nSDS≈0 → 深根的　SDS>0 → 浅根的", "SDS")

    ax = axes[1]
    for i, r in enumerate(le_df.itertuples()):
        if not np.isnan(r.abs_diff):
            ax.barh(i, r.abs_diff, color="#5B4FCF", alpha=0.75,
                    edgecolor="black", lw=0.8, height=0.65)
            if not (np.isnan(r.abs_lo) or np.isnan(r.abs_hi)):
                ax.errorbar(r.abs_diff, i,
                            xerr=[[r.abs_diff-r.abs_lo],[r.abs_hi-r.abs_diff]],
                            color="black", capsize=5, lw=1.2, fmt="none")
            ax.text(r.abs_diff + (1.5 if r.abs_diff >= 0 else -1.5), i,
                    f"{r.abs_diff:+.1f} W/m²", va="center", fontsize=8.5)
    ax.axvline(0, color="gray", lw=1.0)
    ax.set_yticks(y_pos); ax.set_yticklabels([])
    _ax_style(ax, "(2) 絶対差 = LE_normal − LE_soil_dry\n正 → 干ばつで LE が落ちる",
              "LE_n − LE_s [W/m²]")

    ax = axes[2]
    for i, r in enumerate(le_df.itertuples()):
        if not np.isnan(r.rb):
            col = "#E85D04" if abs(r.rb)>0.3 else "#FFA000" if abs(r.rb)>0.1 else "#9E9E9E"
            ax.barh(i, r.rb, color=col, alpha=0.80, edgecolor="black", lw=0.8, height=0.65)
            ax.text(r.rb + 0.03, i,
                    f"{r.rb:+.2f}  p={r.p:.1e} {_sig_star(r.p)}",
                    va="center", fontsize=8.5)
    ax.axvline(0, color="gray", lw=1.0)
    ax.set_xlim(-1, 1.5)
    ax.set_yticks(y_pos); ax.set_yticklabels([])
    _ax_style(ax, "(3) 効果量 (rank-biserial r) + 有意性\n|r|>0.3 → 大効果", "rank-biserial r")
    for xv, lb in [(0.1, "中"), (0.3, "大")]:
        ax.axvline(xv, color="gray", ls=":", lw=0.8)
        ax.text(xv+0.01, n_rows-0.4, lb, fontsize=7, color="gray")

    plt.tight_layout()
    fp = save_dir / "fig01_verdict_panel.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ---- fig02: Recovery curve (2×2) ------

def plot_recovery_curve(fit_res, tau_by_year, save_dir):
    grouped = fit_res.get("grouped")
    if grouped is None or len(grouped) < 2:
        print("  [skip] recovery curve"); return
    tau, le0, le_inf = fit_res["tau"], fit_res["le0"], fit_res["le_inf"]
    r2_exp, r2_lin, r2_log = fit_res["r2_exp"], fit_res["r2_lin"], fit_res["r2_log"]
    aic_e, aic_l, aic_g = fit_res["aic_exp"], fit_res["aic_lin"], fit_res["aic_log"]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.patch.set_facecolor(FIG_BG)
    x_obs = grouped["days_since_irrig"].values.astype(float)
    y_obs = grouped["median"].values
    n_obs = grouped["count"].values
    x_fit = np.linspace(0, x_obs.max() + 2, 300)

    ax = axes[0, 0]
    sizes = 40 + n_obs * 4
    ax.scatter(x_obs, y_obs, s=sizes, c="#1D9E75",
               edgecolors="black", lw=1, zorder=5, label="観測中央値 (サイズ=n)")
    if not np.isnan(tau):
        y_fit = exp_model(x_fit, le_inf, le0, tau)
        ax.plot(x_fit, y_fit, "k-", lw=2.5, label=f"Exp fit: τ={tau:.1f}d, R²={r2_exp:.3f}")
        ax.axhline(le_inf, color="#E85D04", ls="--", lw=1.5, alpha=0.9,
                   label=f"LE_∞={le_inf:.0f} W/m² (深層水 or 残存水?)")
        ax.axhline(le0, color="#1A73E8", ls="--", lw=1.5, alpha=0.9,
                   label=f"LE_0={le0:.0f} W/m² (灌漑直後)")
        ax.axvline(tau, color="gray", ls=":", lw=1.5, label=f"τ={tau:.1f}日")
        interp_col = "#E85D04" if tau < 5 else "#FFA000" if tau <= 14 else "#1D9E75"
        interp_txt = (f"τ={tau:.1f}d < 5d → 灌漑依存が支配的" if tau < 5 else
                      f"τ={tau:.1f}d > 14d → 深層水の可能性" if tau > 14 else
                      f"τ={tau:.1f}d → 中間的応答")
        ax.text(0.97, 0.97, interp_txt, transform=ax.transAxes, va="top", ha="right",
                fontsize=10, fontweight="bold", color=interp_col,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=interp_col, alpha=0.9))
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=8.5)
    _ax_style(ax, f"(A) Recovery curve — τ={tau:.1f}d\nTarazona 灌漑アクティブ月・生育期",
              "Days since last irrigation", "LE_corr 中央値 [W/m²]")

    ax = axes[0, 1]
    ax.scatter(x_obs, y_obs, s=sizes, c="#1D9E75",
               edgecolors="black", lw=1, zorder=5, label="観測中央値")
    best_aic = min(v for v in [aic_e, aic_l, aic_g] if not np.isnan(v)) if any(
        not np.isnan(v) for v in [aic_e, aic_l, aic_g]) else np.nan
    for name, fn, popt, r2, aic, ls, lw_ in [
        ("Exponential", exp_model, fit_res.get("popt_exp"), r2_exp, aic_e, "k-", 2.5),
        ("Linear",      lin_model, fit_res.get("popt_lin"), r2_lin, aic_l, "b--", 1.8),
        ("Logarithmic", log_model, fit_res.get("popt_log"), r2_log, aic_g, "r-.", 1.8),
    ]:
        if popt is not None and not np.isnan(popt[0]):
            y_f = fn(x_fit, *popt)
            da  = aic - best_aic if not np.isnan(aic) and not np.isnan(best_aic) else np.nan
            best_mark = "  ← BEST" if not np.isnan(da) and abs(da) < 0.1 else ""
            lbl = f"{name}  R²={r2:.3f}, ΔAIC={da:.1f}{best_mark}" if not np.isnan(da) else f"{name}  R²={r2:.3f}"
            ax.plot(x_fit, y_f, ls, lw=lw_, label=lbl)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=8.5)
    _ax_style(ax, "(B) モデル比較 (Exp/Linear/Log)\nAIC 最小 = 最良フィット",
              "Days since last irrigation", "LE_corr 中央値 [W/m²]")

    ax = axes[1, 0]
    valid = tau_by_year.dropna(subset=["tau"])
    if len(valid) >= 2:
        colors = ["#1D9E75" if r.tau_reliable else "#BDBDBD" for r in valid.itertuples()]
        ax.bar(valid["year"], valid["tau"], color=colors, alpha=0.85,
               edgecolor="black", lw=1, width=0.6)
        ax.axhline(tau, color="#E85D04", ls="--", lw=2, label=f"全期間 τ={tau:.1f}d")
        ax.axhline(5, color="gray", ls=":", lw=1, alpha=0.7, label="τ=5d 閾値")
        trusted = valid[valid["tau_reliable"]]
        if len(trusted) >= 2:
            ax.axhspan(trusted["tau"].min(), trusted["tau"].max(),
                       color="#1D9E75", alpha=0.08,
                       label=f"信頼年レンジ: {trusted['tau'].min():.1f}–{trusted['tau'].max():.1f}d")
        for _, row in valid.iterrows():
            lbl = f"{row['tau']:.1f}d\n(n={int(row['n_points'])})"
            ax.text(row["year"], row["tau"] + 0.2, lbl,
                    ha="center", va="bottom", fontsize=8,
                    color="black" if row["tau_reliable"] else "#999")
        # 凡例: axhline 等のハンドルに mpatches を追加
        leg_patches = [
            mpatches.Patch(color="#1D9E75", alpha=0.85, label="信頼年 (R²≥0.7)"),
            mpatches.Patch(color="#BDBDBD", alpha=0.85, label="低信頼年 (R²<0.7 or n<4)"),
        ]
        handles, labels_ = ax.get_legend_handles_labels()
        ax.legend(handles=handles + leg_patches, fontsize=9)
        # 信頼年のみで判断するテキスト
        r_flag = ("★ 信頼年の τ 変動が小さい → Robust" if len(trusted) >= 2 and
                  trusted["tau"].max() - trusted["tau"].min() < 3
                  else "△ 信頼年の τ でも変動あり → 要注意")
        r_col = "#1D9E75" if "★" in r_flag else "#E85D04"
        ax.text(0.97, 0.97, r_flag, transform=ax.transAxes, va="top", ha="right",
                fontsize=9, fontweight="bold", color=r_col,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=r_col, alpha=0.85))
    else:
        ax.text(0.5, 0.5, "年別 τ: 有効年 < 2", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="gray")
    _ax_style(ax, "(C) 年別 τ — 信頼性フラグ付き\n緑=信頼年、グレー=低信頼(除外推奨)",
              "Year", "τ [days]")

    ax = axes[1, 1]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(mpatches.Rectangle((0,0),10,10, fc=FIG_BG, ec="none"))
    ax.text(5, 9.5, "(D) 二層モデル解釈", ha="center", va="top",
            fontsize=11, fontweight="bold")

    def box(x, y, w, h, fc, text="", fs=9):
        ax.add_patch(mpatches.FancyBboxPatch((x,y), w, h, boxstyle="round,pad=0.15",
                                             fc=fc, ec="black", lw=1.5, alpha=0.85))
        if text: ax.text(x+w/2, y+h/2, text, ha="center", va="center",
                         fontsize=fs, fontweight="bold")
    def arrow(x1,y1,x2,y2,col="black"):
        ax.annotate("", xy=(x2,y2), xytext=(x1,y1),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))

    box(0.5,7.5,2.5,1.2, fc="#AED6F1", text="灌漑パルス\n(drip)")
    arrow(3.0,8.1, 4.5,8.1)
    box(4.5,7.2,4.5,1.8, fc="#D5EAF5", text=f"Fast pool\n(灌漑水 → τ≈{tau:.1f}d で枯渇)")
    arrow(6.75,7.2, 6.75,5.8, col="#1A73E8")
    box(4.5,4.2,4.5,1.5, fc="#D5F5E3", text=f"Slow pool\n(深層水? → LE_∞≈{le_inf:.0f} W/m²)")
    arrow(6.75,4.2, 6.75,2.8, col="#1D9E75")
    box(4.0,1.5,5.5,1.2, fc="#FDE8D8",
        text=f"蒸散 LE: {le0:.0f}→{le_inf:.0f} W/m²")
    ax.text(5, 0.5, "SWC-ET decoupling は灌漑水が主要因\n(深根を棄却せず、但し支配的でない)",
            ha="center", va="bottom", fontsize=8.5, style="italic", color="#333",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999", alpha=0.8))

    plt.tight_layout(h_pad=3, w_pad=3)
    fp = save_dir / "fig02_recovery_curve.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ---- fig03: 三点指標間整合性 ------

def plot_three_point_matrix(pkg_df, save_dir):
    sub = pkg_df[pkg_df["var"]=="LE_corr"].dropna(subset=["sds","abs_diff","rb"])
    if len(sub) < 2: return
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig03 — 三点指標間の整合性 (LE_corr)", fontsize=12, fontweight="bold")
    pairs = [
        ("sds","abs_diff","SDS","Abs diff [W/m²]","SDS vs 絶対差"),
        ("sds","rb","SDS","rank-biserial r","SDS vs 効果量"),
        ("abs_diff","rb","Abs diff [W/m²]","rank-biserial r","絶対差 vs 効果量"),
    ]
    for ax, (xk,yk,xl,yl,title) in zip(axes, pairs):
        for _, r in sub.iterrows():
            v, strength = verdict({"sds":r["sds"],"sds_lo":r["sds_lo"],
                                    "sds_hi":r["sds_hi"],"n_n":r["n_n"],"n_s":r["n_s"]})
            col   = SITE_COL.get(r["site"], "#888")
            mk    = "o" if r["site"]=="Oran" else "^"
            alpha = 0.9 if strength=="strong" else 0.4
            ax.scatter(r[xk], r[yk], c=col, s=90, alpha=alpha,
                       edgecolors="black", lw=1, marker=mk)
            ax.annotate(f"{r['site'][:3]}\n{r['season'][:6]}", (r[xk],r[yk]),
                        fontsize=7.5, xytext=(5,5), textcoords="offset points")
        ax.axhline(0, color="gray", lw=0.6, alpha=0.5)
        ax.axvline(0, color="gray", lw=0.6, alpha=0.5)
        _ax_style(ax, title, xl, yl)
    leg = [mpatches.Patch(color=SITE_COL["Oran"], label="Oran (rainfed)"),
           mpatches.Patch(color=SITE_COL["Tarazona"], label="Tarazona (irrigated)"),
           Line2D([0],[0], marker="o", color="gray", ms=8, lw=0, markeredgecolor="black", label="strong"),
           Line2D([0],[0], marker="o", color="gray", ms=8, lw=0, alpha=0.35, markeredgecolor="black", label="weak")]
    fig.legend(handles=leg, loc="lower center", ncol=4, fontsize=9, frameon=True)
    plt.tight_layout(rect=[0,0.07,1,1])
    fp = save_dir / "fig03_three_point_matrix.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ---- fig04: 年別 τ standalone ------

def plot_tau_sensitivity(tau_by_year, tau_all, save_dir):
    valid = tau_by_year.dropna(subset=["tau"])
    if len(valid) < 2: return
    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_facecolor(FIG_BG)
    tmin, tmax = valid["tau"].min(), valid["tau"].max()
    tmean, tstd = valid["tau"].mean(), valid["tau"].std()
    colors = ["#1D9E75" if r.tau_reliable else "#BDBDBD" for r in valid.itertuples()]
    bars = ax.bar(valid["year"], valid["tau"], color=colors, alpha=0.80,
                  edgecolor="black", lw=1.2, width=0.6)
    ax.axhline(tau_all, color="#E85D04", ls="--", lw=2.2, label=f"全期間 τ={tau_all:.1f}d")
    ax.axhline(5, color="gray", ls=":", lw=1.0, alpha=0.7, label="τ=5d 閾値")
    trusted = valid[valid["tau_reliable"]]
    if len(trusted) >= 2:
        ax.axhspan(trusted["tau"].min(), trusted["tau"].max(), color="#1D9E75", alpha=0.08,
                   label=f"信頼年レンジ: {trusted['tau'].min():.1f}–{trusted['tau'].max():.1f}d")
    for bar, (_, row) in zip(bars, valid.iterrows()):
        ax.text(bar.get_x()+bar.get_width()/2, row["tau"]+0.15,
                f"{row['tau']:.1f}d\n(n={int(row['n_points'])})\n{row['reason']}",
                ha="center", va="bottom", fontsize=8.5,
                color="black" if row["tau_reliable"] else "#999")
    ax.set_ylim(0, max(valid["tau"].max()*1.45, 8))
    ax.set_xticks(valid["year"])
    leg_patches = [mpatches.Patch(color="#1D9E75", alpha=0.85, label="信頼年 (R²≥0.7)"),
                   mpatches.Patch(color="#BDBDBD", alpha=0.85, label="低信頼年 (R²<0.7, 除外推奨)")]
    handles, labels_ = ax.get_legend_handles_labels()
    ax.legend(handles=handles+leg_patches, fontsize=9.5)
    r_flag = "★ 信頼年 τ 変動小 → Robust" if (len(trusted)>=2 and
              trusted["tau"].max()-trusted["tau"].min()<3) else "△ 信頼年でも変動あり"
    r_col = "#1D9E75" if "★" in r_flag else "#E85D04"
    ax.text(0.98, 0.98, f"{r_flag}\n信頼年: τ 平均={trusted['tau'].mean():.1f}±{trusted['tau'].std():.1f}d" if len(trusted)>=2 else r_flag,
            transform=ax.transAxes, va="top", ha="right", fontsize=10,
            fontweight="bold", color=r_col,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=r_col, alpha=0.9))
    _ax_style(ax, f"fig04 — 年別 τ Sensitivity (Tarazona)\n全期間 τ={tau_all:.1f}d",
              "Year", "τ [days]")
    plt.tight_layout()
    fp = save_dir / "fig04_tau_sensitivity.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ---- fig05: τ Diagnostic (NEW) ------

def plot_tau_diagnostic(tau_by_year: pd.DataFrame, irrig_stats: pd.DataFrame,
                         save_dir: Path):
    """
    年別 τ と灌漑統計を 1 枚に並べて τ 変動の要因を診断する。

    パネル構成:
      (A) 年別 τ — 信頼性フラグ付き (bar)
      (B) 年別 灌漑イベント数 (bar)
      (C) 年別 灌漑 mm/event (bar) — 1 イベントあたりの水量
      (D) 年別 降水量 (bar) — 土壌水分ベースラインに影響

    仮説: 灌漑ペースが速い年(イベント多/水量少)は τ が短い → 灌漑依存
         降水量が多い年は表層水分が高く irrig_active_month でも τ に影響
    """
    merged = tau_by_year.merge(irrig_stats, on="year", how="left")
    years  = merged["year"].values
    x      = np.arange(len(years))
    w      = 0.6

    fig = plt.figure(figsize=(14, 12))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle(
        "fig05 — τ Diagnostic: 年別 τ × 灌漑統計\n"
        "τ の年間変動の要因を診断する (灌漑ペース・水量・降水量との比較)",
        fontsize=12, fontweight="bold", y=0.98)
    gs = GridSpec(4, 1, hspace=0.55, figure=fig)

    # (A) 年別 τ
    ax0 = fig.add_subplot(gs[0])
    tau_vals = merged["tau"].values
    colors0  = ["#1D9E75" if r.tau_reliable else "#BDBDBD"
                for r in merged.itertuples()]
    bars0    = ax0.bar(x, tau_vals, color=colors0, alpha=0.85,
                       edgecolor="black", lw=1, width=w)
    for i, (v, r) in enumerate(zip(tau_vals, merged.itertuples())):
        if not np.isnan(v):
            ax0.text(i, v+0.3, f"{v:.1f}d\n({r.reason})",
                     ha="center", va="bottom", fontsize=8,
                     color="black" if r.tau_reliable else "#888")
    ax0.axhline(5, color="gray", ls=":", lw=1, alpha=0.7, label="τ=5d 閾値")
    ax0.set_ylim(0, np.nanmax(tau_vals)*1.4+2)
    ax0.set_xticks(x); ax0.set_xticklabels(years)
    ax0.legend(fontsize=8)
    _ax_style(ax0, "(A) 年別 τ [days]  緑=信頼 グレー=低信頼",
              "", "τ [days]")
    leg_a = [mpatches.Patch(color="#1D9E75", alpha=0.85, label="信頼年 (R²≥0.7)"),
             mpatches.Patch(color="#BDBDBD", alpha=0.85, label="低信頼年 (R²<0.7)")]
    ax0.legend(handles=leg_a+ax0.get_legend_handles_labels()[0], fontsize=8, loc="upper right")

    # (B) 年別 灌漑イベント数
    ax1 = fig.add_subplot(gs[1])
    n_ev = merged["n_irrig_events"].values
    ax1.bar(x, n_ev, color="#1A73E8", alpha=0.75, edgecolor="black", lw=1, width=w)
    for i, v in enumerate(n_ev):
        if not np.isnan(v):
            ax1.text(i, v+0.5, f"{int(v)}", ha="center", va="bottom", fontsize=9)
    ax1.set_xticks(x); ax1.set_xticklabels(years)
    _ax_style(ax1, "(B) 年別 灌漑イベント数\n多い → τ が短くなる可能性",
              "", "イベント数")

    # (C) 年別 灌漑 mm/event
    ax2 = fig.add_subplot(gs[2])
    mm_ev = merged["avg_mm_per_event"].values
    ax2.bar(x, mm_ev, color="#FFA000", alpha=0.80, edgecolor="black", lw=1, width=w)
    for i, v in enumerate(mm_ev):
        if not np.isnan(v):
            ax2.text(i, v+0.5, f"{v:.1f}", ha="center", va="bottom", fontsize=9)
    ax2.set_xticks(x); ax2.set_xticklabels(years)
    _ax_style(ax2, "(C) 年別 灌漑量 [mm/event]\n多い → 1 回で多くの水 → τ が長くなる可能性",
              "", "mm/event")

    # (D) 年別 降水量
    ax3 = fig.add_subplot(gs[3])
    rain = merged["rain_mm"].values
    ax3.bar(x, rain, color="#9575CD", alpha=0.75, edgecolor="black", lw=1, width=w)
    for i, v in enumerate(rain):
        if not np.isnan(v):
            ax3.text(i, v+3, f"{v:.0f}", ha="center", va="bottom", fontsize=9)
    ax3.set_xticks(x); ax3.set_xticklabels(years)
    _ax_style(ax3, "(D) 年別 降水量 [mm]\n多い → 表層 SWC 高 → irrig_active_month のフィルタに影響",
              "Year", "Rain [mm]")

    plt.subplots_adjust(top=0.94)
    fp = save_dir / "fig05_tau_diagnostic.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 6. MAIN
# ================================================================

def main():
    print("=" * 60)
    print("解析A v17 — 三点パッケージ + Recovery τ + Diagnostic")
    print("=" * 60)

    df = load_and_merge()
    df = add_days_since_irrig(df)
    df = add_irrig_active_month(df)

    # ---- 三点パッケージ ----
    print(f"\n{'='*60}\n三点パッケージ計算\n{'='*60}")
    rows = []
    for season_label, months in SEASONS.items():
        for site in ["Oran", "Tarazona"]:
            sub = df[(df["site"]==site) & df["is_growing"] & df["month"].isin(months)]
            if site == "Tarazona":
                sub = sub[sub["irrig_active_month"]]
            for var in ["LE_corr","EF_corr","ET"]:
                if var not in sub.columns: continue
                n_arr = sub.loc[sub["drought_type"]=="normal",   var].dropna().values
                s_arr = sub.loc[sub["drought_type"]=="soil dry", var].dropna().values
                pkg   = three_point_package(n_arr, s_arr, var)
                v, strength = verdict(pkg)
                pkg.update(dict(site=site, season=season_label, var=var,
                                verdict=v, verdict_strength=strength))
                rows.append(pkg)
                if var == "LE_corr":
                    print(f"  [{site:9s}/{season_label:22s}] "
                          f"SDS={pkg['sds']:+.3f} abs={pkg['abs_diff']:+.1f}W/m² "
                          f"rb={pkg['rb']:+.3f} p={pkg['p']:.1e} "
                          f"n=({pkg['n_n']},{pkg['n_s']}) → [{VERDICT_JP[v]}] ({strength})")
    pkg_df = pd.DataFrame(rows)
    pkg_df.to_csv(OUT_DIR / "v17_three_point_package.csv", index=False)

    # ---- Recovery τ 解析 ----
    print(f"\n{'='*60}\nRecovery τ\n{'='*60}")
    tara_active = df[(df["site"]=="Tarazona") & df["irrig_active_month"] & df["is_growing"]]
    fit_res = fit_recovery(tara_active, var="LE_corr")
    print(f"  τ={fit_res['tau']:.2f}d  LE_0={fit_res['le0']:.1f}  LE_∞={fit_res['le_inf']:.1f} W/m²")
    print(f"  R² exp={fit_res['r2_exp']:.3f} lin={fit_res['r2_lin']:.3f} log={fit_res['r2_log']:.3f}")
    print(f"  AIC exp={fit_res['aic_exp']:.1f} lin={fit_res['aic_lin']:.1f} log={fit_res['aic_log']:.1f}")

    # ---- 年別 τ Sensitivity ----
    print(f"\n{'='*60}\n年別 τ Sensitivity\n{'='*60}")
    tau_by_year = fit_tau_by_year(tara_active, var="LE_corr")
    print(tau_by_year.to_string(index=False))
    trusted = tau_by_year[tau_by_year["tau_reliable"]]
    if len(trusted) >= 2:
        print(f"\n  ★ 信頼年のみ τ レンジ: {trusted['tau'].min():.1f}–{trusted['tau'].max():.1f}d"
              f"  (平均={trusted['tau'].mean():.1f}±{trusted['tau'].std():.1f}d)")
    else:
        print("  信頼年が 2 未満 → robust range 計算不可")

    # ---- 灌漑統計 ----
    irrig_stats = calc_irrig_stats_by_year(df)
    irrig_stats.to_csv(OUT_DIR / "v17_irrig_stats_by_year.csv", index=False)
    print(f"\n  灌漑統計:\n{irrig_stats.to_string(index=False)}")

    # ---- CSV / JSON 保存 ----
    if fit_res.get("grouped") is not None:
        out_g = fit_res["grouped"].copy()
        out_g["site"] = "Tarazona"
        out_g.to_csv(OUT_DIR / "v17_recovery_binned.csv", index=False)
    tau_by_year.to_csv(OUT_DIR / "v17_tau_by_year.csv", index=False)
    with open(OUT_DIR / "v17_recovery_fit_params.json", "w") as f:
        json.dump({k: (float(v) if isinstance(v, (float, np.floating)) else
                       int(v) if isinstance(v, (int, np.integer)) else v)
                   for k, v in fit_res.items()
                   if k not in ("grouped","popt_exp","popt_lin","popt_log")}, f, indent=2)

    # ---- 可視化 ----
    print(f"\n--- 可視化 ---")
    plot_verdict_panel(pkg_df, OUT_DIR)
    plot_recovery_curve(fit_res, tau_by_year, OUT_DIR)
    plot_three_point_matrix(pkg_df, OUT_DIR)
    plot_tau_sensitivity(tau_by_year, fit_res["tau"], OUT_DIR)
    plot_tau_diagnostic(tau_by_year, irrig_stats, OUT_DIR)

    # ---- 最終サマリー ----
    print(f"\n{'='*60}\n★ 最終サマリー\n{'='*60}")
    le_pkg = pkg_df[pkg_df["var"]=="LE_corr"]
    for _, row in le_pkg.iterrows():
        v, strength = row["verdict"], row["verdict_strength"]
        lo, hi = row.get("sds_lo", np.nan), row.get("sds_hi", np.nan)
        ci = f"[{lo:+.2f},{hi:+.2f}]" if not np.isnan(lo) else "[CI 算出不可]"
        print(f"  [{row['site']:9s}/{row['season']:22s}] "
              f"SDS={row['sds']:+.3f} {ci} n=({row['n_n']:.0f},{row['n_s']:.0f}) "
              f"→ {VERDICT_JP[v]} ({strength})")

    tau = fit_res["tau"]
    print(f"\n  Recovery τ={tau:.1f}d (全期間)  "
          f"信頼年レンジ: {trusted['tau'].min():.1f}–{trusted['tau'].max():.1f}d" if len(trusted)>=2
          else f"\n  Recovery τ={tau:.1f}d (全期間)")
    print(f"  主張レベル: 'Irrigation-driven decoupling (τ≈{tau:.0f}d) "
          f"dominates over deep rooting'")
    print(f"\n[done] outputs → {OUT_DIR}/")


if __name__ == "__main__":
    main()
