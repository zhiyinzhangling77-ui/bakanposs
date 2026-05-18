"""
================================================================
解析A v15 : SDS 三点パッケージ + Recovery Time τ 抽出
================================================================

【目的】
  v14 までの解析を統合し、SDS 議論で出た改善点を全て反映:
  1. 判定ロジックのバグ修正 (CI 幅 + n_s ガード)
  2. 三点報告パッケージ (SDS + 絶対差 + rank-biserial)
  3. Recovery time τ の exp fit (灌漑効果の減衰時定数)
  4. irrig_active_month フィルタ(灌漑非アクティブ月の混入を全工程で除外)

【三点パッケージとは】
  単一の SDS だけだと不十分。3 つを併記することで結論を頑健化:
    - SDS         : 相対減少率 (%表現) — 直感的
    - 絶対差      : LE_n − LE_s [W/m²] — 物理量
    - rank-biserial r : ノンパラ効果量 [-1, +1] — 有意性 + 効果量

【判定ロジック修正】
  v14 の `in_band` は "CI が 0 を含む → 深根" を機械的に判定し、
  n_s=15 で CI が広がりすぎた Tarazona spring を偽陽性にした。
  v15 では:
    - n_s < N_MIN_VERDICT (=30) → insufficient_data
    - CI 幅 > CI_WIDTH_MAX (=0.5) → uncertain_wide_CI
    - 上記をクリアしてから "CI が 0 を含む" を深根帯と判定

【Recovery time τ】
  灌漑から経過するほど LE が低下するという v14 fig03 のパターンを、
  指数モデル `LE(d) = LE_∞ + (LE_0 − LE_∞) * exp(-d/τ)` で fit。
  τ が小さい(数日)= 灌漑効果が短期間で減衰 = 灌漑依存
  τ が大きい(>14日)or fit 失敗 = 灌漑効果が持続 or 関係なし

【入力】
  daily_classified_v4.parquet, closure_slopes_v4.json
  Tarazona daily CSV (Irrig_mm)

【出力】./output_analysis_A_v15/
  fig01_verdict_panel.png       : 三点パッケージ統合判定
  fig02_recovery_curve.png      : LE_median vs days_since + exp fit
  fig03_three_point_matrix.png  : 3指標の相互比較
  v15_three_point_package.csv
  v15_recovery_time.csv
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")


# 修正前：
# INPUT_DIR       = Path("./")
# PARQUET         = INPUT_DIR / "daily_classified_v4.parquet"

# 修正後：
INPUT_DIR       = Path("/home/shion-nagamine/bakanposs/analysis_A")
PARQUET         = INPUT_DIR / "daily_classified_v4.parquet"
TARA_DAILY_CSV  = Path(
    "/home/shion-nagamine/Dataset/Eddy data in Spain/"
    "Daily_Summary_Filtered_forPred_ActEne26.csv")
OUT_DIR = Path("./output_analysis_A_v15")
OUT_DIR.mkdir(exist_ok=True)

# 灌漑判定 & active-month フィルタ
IRRIG_THRESHOLD     = 0.5   # mm/day, 灌漑イベント
MIN_IRRIG_PER_MONTH = 2     # 月別灌漑イベント数(これ以上を active 判定)

# 季節定義(参考、Tarazona の灌漑カレンダーに基づく)
SEASONS = {
    "spring (1-4月)"        : [1, 2, 3, 4],
    "shoulder (5,6,10月)"   : [5, 6, 10],
    "summer (7-9月)"        : [7, 8, 9],
}

# Bootstrap
N_BOOT          = 5000
CI_PCT          = (2.5, 97.5)
MIN_N_PER_CLASS = 5

# 三点パッケージ用の guard
DENOM_FLOORS  = {"LE_corr": 5.0, "EF_corr": 0.05, "ET": 0.3}
RATIO_CLIP    = 10.0

# 判定ロジックの guard(v14 の偽陽性回避)
N_MIN_VERDICT   = 30     # n_n, n_s 両方 >= これ未満なら判定保留
CI_WIDTH_MAX    = 0.5    # CI 幅がこれを超えたら判定保留
DEEP_ROOT_BAND  = 0.15   # |SDS| <= これ かつ CI が 0 を含む → 深根

# 可視化
SITE_COL = {"Oran":"#E85D04", "Tarazona":"#1D9E75"}


# ================================================================
# 1. 入力読込
# ================================================================

def load_and_merge():
    """v4 parquet + Tarazona daily CSV (Irrig_mm) を結合"""
    daily = pd.read_parquet(PARQUET)
    daily["date"] = pd.to_datetime(daily["date"])
    csv = pd.read_csv(TARA_DAILY_CSV)
    csv["date"] = pd.to_datetime(csv["date"], errors="coerce")
    irrig_cols = [c for c in ["Irrig_mm","Rain_mm","IrrigRain_mm"]
                  if c in csv.columns]
    extra = (csv[["date"] + irrig_cols]
             .dropna(subset=["date"]).drop_duplicates("date"))
    tara = daily[daily["site"]=="Tarazona"].merge(extra, on="date", how="left")
    oran = daily[daily["site"]=="Oran"].copy()
    for c in irrig_cols:
        oran[c] = 0.0
    return pd.concat([oran, tara]).sort_values(["site","date"]).reset_index(drop=True)


def add_days_since_irrig(df):
    df = df.sort_values(["site","date"]).reset_index(drop=True).copy()
    df["days_since_irrig"] = np.nan
    for site, idx in df.groupby("site").indices.items():
        sub = df.loc[idx].sort_values("date")
        last = None; out = []
        for _, row in sub.iterrows():
            if pd.notna(row.get("Irrig_mm")) and row["Irrig_mm"] > IRRIG_THRESHOLD:
                last = row["date"]; out.append(0)
            elif last is None:
                out.append(np.nan)
            else:
                out.append((row["date"] - last).days)
        df.loc[sub.index, "days_since_irrig"] = out
    return df


def add_irrig_active_month(df):
    """月別灌漑イベント数 >= MIN_IRRIG_PER_MONTH の月を active と判定"""
    df = df.copy()
    df["year_month"] = df["date"].dt.to_period("M")
    df["is_irrig_event"] = df["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD
    monthly_events = (df.groupby(["site","year_month"])["is_irrig_event"]
                        .sum().reset_index(name="n_irrig"))
    active_mask = monthly_events["n_irrig"] >= MIN_IRRIG_PER_MONTH
    active_pairs = set(zip(monthly_events.loc[active_mask, "site"],
                            monthly_events.loc[active_mask, "year_month"]))
    df["irrig_active_month"] = df.apply(
        lambda r: (r["site"], r["year_month"]) in active_pairs, axis=1)
    return df.drop(columns=["year_month", "is_irrig_event"])


# ================================================================
# 2. 三点パッケージ (SDS + 絶対差 + rank-biserial)
# ================================================================

def safe_ratio(numer, denom, var):
    floor = DENOM_FLOORS.get(var, 0.05)
    if abs(denom) < floor: return np.nan
    r = numer / denom
    return np.nan if abs(r) > RATIO_CLIP else r


def three_point_package(arr_n, arr_s, var, seed=42):
    """
    SDS + 絶対差 + rank-biserial を bootstrap 付きで返す
    SDS = 1 - median(LE_s)/median(LE_n)
    abs_diff = median(LE_n) - median(LE_s)
    rb = 1 - 2U/(n_n * n_s)  (ノンパラ効果量、Mann-Whitney 由来)
    """
    if len(arr_n) < MIN_N_PER_CLASS or len(arr_s) < MIN_N_PER_CLASS:
        return dict(sds=np.nan, sds_lo=np.nan, sds_hi=np.nan,
                    abs_diff=np.nan, abs_lo=np.nan, abs_hi=np.nan,
                    rb=np.nan, p=np.nan,
                    med_n=np.nan, med_s=np.nan,
                    n_n=len(arr_n), n_s=len(arr_s))

    med_n, med_s = float(np.median(arr_n)), float(np.median(arr_s))
    sds_point = safe_ratio(med_n - med_s, med_n, var)
    abs_point = med_n - med_s

    # Bootstrap SDS と abs_diff
    rng = np.random.default_rng(seed)
    boot_sds, boot_abs = [], []
    for _ in range(N_BOOT):
        nb = rng.choice(arr_n, size=len(arr_n), replace=True)
        sb = rng.choice(arr_s, size=len(arr_s), replace=True)
        m_n, m_s = np.median(nb), np.median(sb)
        r = safe_ratio(m_n - m_s, m_n, var)
        if not np.isnan(r):
            boot_sds.append(r)
        boot_abs.append(m_n - m_s)

    sds_lo, sds_hi = (np.nan, np.nan)
    if len(boot_sds) >= N_BOOT * 0.1:
        sds_lo, sds_hi = np.percentile(boot_sds, CI_PCT)
    abs_lo, abs_hi = np.percentile(boot_abs, CI_PCT)

    # Mann-Whitney + rank-biserial
    try:
        u, p = stats.mannwhitneyu(arr_n, arr_s, alternative="two-sided")
        rb = 1 - 2 * u / (len(arr_n) * len(arr_s))
    except ValueError:
        u, p, rb = np.nan, np.nan, np.nan

    return dict(sds=sds_point, sds_lo=sds_lo, sds_hi=sds_hi,
                abs_diff=abs_point, abs_lo=abs_lo, abs_hi=abs_hi,
                rb=rb, p=p, med_n=med_n, med_s=med_s,
                n_n=len(arr_n), n_s=len(arr_s))


def verdict(pkg):
    """
    判定ロジック(v14 偽陽性回避):
      1. n_n または n_s < N_MIN_VERDICT     → insufficient_data
      2. CI 幅 > CI_WIDTH_MAX                → uncertain_wide_CI
      3. CI が 0 を含み |SDS|<=DEEP_ROOT_BAND → deep_root
      4. CI 下限 > 0                         → shallow
      5. CI 上限 < 0                         → negative (異常)
      6. その他                              → uncertain
    """
    if (np.isnan(pkg["sds"]) or pkg["n_n"] < N_MIN_VERDICT
            or pkg["n_s"] < N_MIN_VERDICT):
        return "insufficient_data"
    if np.isnan(pkg["sds_lo"]):
        return "uncertain"
    width = pkg["sds_hi"] - pkg["sds_lo"]
    if width > CI_WIDTH_MAX:
        return "uncertain_wide_CI"
    if (pkg["sds_lo"] <= 0 <= pkg["sds_hi"]) and abs(pkg["sds"]) <= DEEP_ROOT_BAND:
        return "deep_root"
    if pkg["sds_lo"] > 0:
        return "shallow"
    if pkg["sds_hi"] < 0:
        return "negative_anomaly"
    return "uncertain"


# ================================================================
# 3. Recovery time τ の指数 fit
# ================================================================

def exp_model(d, le_inf, le0, tau):
    """LE(d) = LE_∞ + (LE_0 − LE_∞) * exp(-d/τ)"""
    return le_inf + (le0 - le_inf) * np.exp(-d / tau)


def fit_recovery(df_active, var="LE_corr"):
    """
    Tarazona irrig_active_month 内で、days_since_irrig 別に LE 中央値を取り、
    指数モデルを fit して τ を抽出
    """
    sub = df_active[df_active[var].notna() &
                    df_active["days_since_irrig"].notna()]
    if len(sub) < 30:
        return dict(tau=np.nan, le0=np.nan, le_inf=np.nan,
                    r2=np.nan, n_points=0)

    # days_since_irrig 別の中央値 (>=5 サンプルあるビンのみ)
    grouped = (sub.groupby("days_since_irrig")[var]
                  .agg(["median","count"]).reset_index())
    grouped = grouped[grouped["count"] >= 5]
    if len(grouped) < 4:
        return dict(tau=np.nan, le0=np.nan, le_inf=np.nan,
                    r2=np.nan, n_points=len(grouped),
                    grouped=grouped)

    x = grouped["days_since_irrig"].values.astype(float)
    y = grouped["median"].values

    # 初期値: LE_0=最大値、LE_∞=最小値、τ=平均値
    try:
        popt, _ = curve_fit(
            exp_model, x, y,
            p0=[y.min(), y.max(), 5.0],
            bounds=([0, 0, 0.5], [np.inf, np.inf, 100]),
            maxfev=5000)
        le_inf, le0, tau = popt
        y_pred = exp_model(x, *popt)
        ss_res = np.sum((y - y_pred)**2)
        ss_tot = np.sum((y - y.mean())**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    except Exception as e:
        print(f"  fit error: {e}")
        return dict(tau=np.nan, le0=np.nan, le_inf=np.nan,
                    r2=np.nan, n_points=len(grouped), grouped=grouped)

    return dict(tau=tau, le0=le0, le_inf=le_inf,
                r2=r2, n_points=len(grouped), grouped=grouped,
                popt=popt)


# ================================================================
# 4. 可視化
# ================================================================

VERDICT_COL = {
    "deep_root"          : "#1D9E75",
    "shallow"            : "#E85D04",
    "negative_anomaly"   : "#7C3AED",
    "uncertain"          : "#9E9E9E",
    "uncertain_wide_CI"  : "#BDBDBD",
    "insufficient_data"  : "#E0E0E0",
}

VERDICT_LABEL = {
    "deep_root"          : "★深根",
    "shallow"            : "浅根",
    "negative_anomaly"   : "異常",
    "uncertain"          : "不明",
    "uncertain_wide_CI"  : "CI広",
    "insufficient_data"  : "n不足",
}


def plot_verdict_panel(pkg_df, save_dir):
    """
    各 (site, season) について三点パッケージを縦に並べた包括図
    """
    rows = pkg_df.sort_values(["site","season"])
    n_rows = len(rows)
    if n_rows == 0: return

    fig, axes = plt.subplots(1, 3, figsize=(20, max(5, n_rows*0.6)))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle(
        "v15 三点パッケージ — SDS + 絶対差 + rank-biserial を統合判定\n"
        f"判定ガード: n>={N_MIN_VERDICT}, CI幅<={CI_WIDTH_MAX}",
        fontsize=12, fontweight="bold")

    y_pos = np.arange(n_rows)
    labels = [f"{r.site}\n{r.season}" for r in rows.itertuples()]

    # Panel 1: SDS
    ax = axes[0]
    for i, r in enumerate(rows.itertuples()):
        v = verdict({"sds":r.sds, "sds_lo":r.sds_lo, "sds_hi":r.sds_hi,
                     "n_n":r.n_n, "n_s":r.n_s})
        col = VERDICT_COL[v]
        if not np.isnan(r.sds):
            ax.barh(i, r.sds, color=col, alpha=0.85, edgecolor="black", lw=0.8)
            if not np.isnan(r.sds_lo):
                ax.errorbar(r.sds, i,
                            xerr=[[r.sds-r.sds_lo],[r.sds_hi-r.sds]],
                            color="black", capsize=6, lw=1.2)
            ax.text(r.sds + 0.02, i, f"{r.sds:+.2f} [{VERDICT_LABEL[v]}]",
                    va="center", fontsize=8, fontweight="bold")
    ax.axvline(0, color="gray", lw=0.8)
    ax.axvspan(-DEEP_ROOT_BAND, DEEP_ROOT_BAND, color="#1D9E75", alpha=0.08)
    ax.set_yticks(y_pos); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("SDS"); ax.set_title("(1) SDS (相対減少率)", fontweight="bold")
    ax.grid(axis="x", alpha=0.25)

    # Panel 2: Absolute diff
    ax = axes[1]
    for i, r in enumerate(rows.itertuples()):
        if not np.isnan(r.abs_diff):
            ax.barh(i, r.abs_diff, color="#7C3AED", alpha=0.7,
                    edgecolor="black", lw=0.8)
            if not np.isnan(r.abs_lo):
                ax.errorbar(r.abs_diff, i,
                            xerr=[[r.abs_diff-r.abs_lo],[r.abs_hi-r.abs_diff]],
                            color="black", capsize=6, lw=1.2)
            ax.text(r.abs_diff + (1 if r.abs_diff >= 0 else -1), i,
                    f"{r.abs_diff:+.1f} W/m²",
                    va="center", fontsize=8)
    ax.axvline(0, color="gray", lw=0.8)
    ax.set_yticks(y_pos); ax.set_yticklabels([])
    ax.set_xlabel("LE_n − LE_s [W/m²]")
    ax.set_title("(2) 絶対差 (物理量)", fontweight="bold")
    ax.grid(axis="x", alpha=0.25)

    # Panel 3: rank-biserial + p
    ax = axes[2]
    for i, r in enumerate(rows.itertuples()):
        if not np.isnan(r.rb):
            col = ("#E85D04" if abs(r.rb) > 0.3 else
                   "#FFA000" if abs(r.rb) > 0.1 else "#9E9E9E")
            ax.barh(i, r.rb, color=col, alpha=0.85,
                    edgecolor="black", lw=0.8)
            star = "★★★" if r.p < 0.001 else ("★★" if r.p < 0.01 else
                    ("★" if r.p < 0.05 else "n.s."))
            ax.text(r.rb + 0.02, i,
                    f"{r.rb:+.2f}  p={r.p:.1e} {star}",
                    va="center", fontsize=8)
    ax.axvline(0, color="gray", lw=0.8)
    ax.set_xlim(-1, 1)
    ax.set_yticks(y_pos); ax.set_yticklabels([])
    ax.set_xlabel("rank-biserial r")
    ax.set_title("(3) 効果量 + 有意性", fontweight="bold")
    ax.grid(axis="x", alpha=0.25)

    plt.tight_layout()
    fp = save_dir / "fig01_verdict_panel.png"
    plt.savefig(fp, dpi=140, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [save] {fp}")


def plot_recovery_curve(fit_res, save_dir):
    """LE 中央値 vs days_since_irrig + 指数 fit"""
    if "grouped" not in fit_res or fit_res["grouped"] is None:
        print("  [skip] recovery fit が実行されず"); return
    g = fit_res["grouped"]
    tau, le0, le_inf, r2 = fit_res["tau"], fit_res["le0"], fit_res["le_inf"], fit_res["r2"]

    fig, ax = plt.subplots(figsize=(11, 6.5))
    fig.patch.set_facecolor("#F8F8F8")

    # 観測中央値(点の大きさはサンプル数)
    sizes = 30 + g["count"] * 3
    ax.scatter(g["days_since_irrig"], g["median"], s=sizes, c="#1D9E75",
                edgecolors="black", lw=1, zorder=5,
                label="観測中央値 (サイズ = n)")

    # Fit curve
    if not np.isnan(tau):
        x_fit = np.linspace(0, g["days_since_irrig"].max(), 200)
        y_fit = exp_model(x_fit, le_inf, le0, tau)
        ax.plot(x_fit, y_fit, "k-", lw=2.2, alpha=0.85,
                 label=f"Exp fit: τ = {tau:.1f}d, R² = {r2:.3f}")
        ax.axhline(le_inf, color="red", ls=":", lw=1.0, alpha=0.7,
                    label=f"LE_∞ = {le_inf:.1f} W/m² (枯渇)")
        ax.axhline(le0, color="blue", ls=":", lw=1.0, alpha=0.7,
                    label=f"LE_0 = {le0:.1f} W/m² (灌漑直後)")
        ax.axvline(tau, color="gray", ls="--", lw=1.0, alpha=0.6,
                    label=f"τ = {tau:.1f} 日")

    ax.set_xlabel("Days since last irrigation")
    ax.set_ylabel("LE_corr median (W/m²)")
    interp_txt = "灌漑依存(τ短い)" if tau < 7 else ("ほぼ無関係(τ大)" if tau > 14 else "中間")
    ax.set_title(
        f"[Tarazona, 灌漑アクティブ月] Recovery curve\n"
        f"τ = {tau:.1f}日 → 解釈: {interp_txt}",
        fontweight="bold")
    ax.legend(loc="best", fontsize=9)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    fp = save_dir / "fig02_recovery_curve.png"
    plt.savefig(fp, dpi=140, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [save] {fp}")


def plot_three_point_matrix(pkg_df, save_dir):
    """3指標の相互比較(SDS × abs × rb)"""
    sub = pkg_df.dropna(subset=["sds","abs_diff","rb"])
    if len(sub) < 2: return
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle("v15 三点パッケージ — 指標間の整合性確認",
                 fontsize=12, fontweight="bold")

    pairs = [
        ("sds", "abs_diff", "SDS", "Abs diff [W/m²]"),
        ("sds", "rb", "SDS", "rank-biserial r"),
        ("abs_diff", "rb", "Abs diff [W/m²]", "rank-biserial r"),
    ]
    for ax, (x, y, xl, yl) in zip(axes, pairs):
        for _, r in sub.iterrows():
            col = SITE_COL.get(r["site"], "#888")
            ax.scatter(r[x], r[y], c=col, s=80, alpha=0.85,
                        edgecolors="black", lw=1)
            ax.annotate(f"{r['site'][:3]}\n{r['season'][:6]}",
                         (r[x], r[y]), fontsize=7, xytext=(5,5),
                         textcoords="offset points")
        ax.axhline(0, color="gray", lw=0.6, alpha=0.5)
        ax.axvline(0, color="gray", lw=0.6, alpha=0.5)
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        ax.grid(alpha=0.25)
    plt.tight_layout()
    fp = save_dir / "fig03_three_point_matrix.png"
    plt.savefig(fp, dpi=140, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 5. MAIN
# ================================================================

def main():
    print("=" * 60)
    print("解析A v15 — 三点パッケージ + Recovery τ")
    print("=" * 60)

    df = load_and_merge()
    df = add_days_since_irrig(df)
    df = add_irrig_active_month(df)
    df["month"] = df["date"].dt.month

    n_active_tara = ((df["site"]=="Tarazona") & df["irrig_active_month"]).sum()
    print(f"  Tarazona 灌漑アクティブ日数: {n_active_tara}")

    # ============================================================
    # 三点パッケージを (site, season) 別に計算
    # ============================================================
    print(f"\n{'='*60}\n三点パッケージ計算\n{'='*60}")
    rows = []
    for season_label, months in SEASONS.items():
        for site in ["Oran","Tarazona"]:
            sub = df[(df["site"]==site) & df["is_growing"] &
                      df["month"].isin(months)]
            # Tarazona の場合は灌漑アクティブ月のみに限定
            if site == "Tarazona":
                sub = sub[sub["irrig_active_month"]]
            for var in ["LE_corr","EF_corr","ET"]:
                if var not in sub.columns: continue
                n_arr = sub.loc[sub["drought_type"]=="normal", var].dropna().values
                s_arr = sub.loc[sub["drought_type"]=="soil dry", var].dropna().values
                pkg = three_point_package(n_arr, s_arr, var)
                v = verdict(pkg)
                pkg.update(dict(site=site, season=season_label, var=var, verdict=v))
                rows.append(pkg)
                if var == "LE_corr":
                    print(f"  [{site}/{season_label}] LE_corr "
                          f"SDS={pkg['sds']:+.3f}  abs={pkg['abs_diff']:+.1f}W/m²  "
                          f"rb={pkg['rb']:+.3f}  p={pkg['p']:.1e}  "
                          f"n=({pkg['n_n']},{pkg['n_s']}) → [{VERDICT_LABEL[v]}]")
    pkg_df = pd.DataFrame(rows)
    pkg_df.to_csv(OUT_DIR / "v15_three_point_package.csv", index=False)
    print(f"\n  [save] v15_three_point_package.csv")

    # ============================================================
    # Recovery time τ の指数 fit (Tarazona 灌漑アクティブ月)
    # ============================================================
    print(f"\n{'='*60}\nRecovery time τ 抽出\n{'='*60}")
    tara_active = df[(df["site"]=="Tarazona") & df["irrig_active_month"] &
                      df["is_growing"]]
    fit_res = fit_recovery(tara_active, var="LE_corr")
    print(f"  τ = {fit_res['tau']:.2f} 日")
    print(f"  LE_0 (d=0) = {fit_res['le0']:.1f} W/m²")
    print(f"  LE_∞ (d→∞) = {fit_res['le_inf']:.1f} W/m²")
    print(f"  R² = {fit_res['r2']:.3f}")
    print(f"  fit に使ったビン数: {fit_res['n_points']}")

    # CSV: 観測中央値とfitパラメータ
    if "grouped" in fit_res and fit_res["grouped"] is not None:
        out_csv = fit_res["grouped"].copy()
        out_csv["site"] = "Tarazona"
        out_csv["filter"] = "irrig_active_month==True"
        out_csv.to_csv(OUT_DIR / "v15_recovery_time.csv", index=False)
        with open(OUT_DIR / "v15_recovery_fit_params.json", "w") as f:
            json.dump(dict(tau=float(fit_res["tau"]),
                           le0=float(fit_res["le0"]),
                           le_inf=float(fit_res["le_inf"]),
                           r2=float(fit_res["r2"]),
                           n_points=int(fit_res["n_points"])), f, indent=2)
        print(f"  [save] v15_recovery_time.csv + v15_recovery_fit_params.json")

    # ============================================================
    # 可視化
    # ============================================================
    print(f"\n--- 可視化 ---")
    plot_verdict_panel(pkg_df, OUT_DIR)
    plot_recovery_curve(fit_res, OUT_DIR)
    plot_three_point_matrix(pkg_df, OUT_DIR)

    # ============================================================
    # 最終判定
    # ============================================================
    print(f"\n{'='*60}\n★ 最終判定(LE_corr ベース)\n{'='*60}")
    le_pkg = pkg_df[pkg_df["var"]=="LE_corr"]
    for _, row in le_pkg.iterrows():
        v = row["verdict"]
        marker = VERDICT_LABEL.get(v, v)
        sds_ci = (f"[{row['sds_lo']:+.2f}, {row['sds_hi']:+.2f}]"
                   if not np.isnan(row.get("sds_lo", np.nan)) else "[CI算出不可]")
        print(f"  [{row['site']:9s}/{row['season']:20s}] "
              f"SDS={row['sds']:+.3f} {sds_ci}  "
              f"n=({row['n_n']:3d},{row['n_s']:3d})  → {marker}")

    print(f"\n  Recovery τ = {fit_res['tau']:.1f} 日")
    if fit_res['tau'] < 5:
        print(f"  → 灌漑から短期間で蒸散低下 = **灌漑依存**(深根仮説却下)")
    elif fit_res['tau'] > 14:
        print(f"  → 灌漑停止しても長期間蒸散維持 = **深根仮説支持**(灌漑非依存)")
    else:
        print(f"  → 中間的応答 ({fit_res['tau']:.1f}日)、解釈は他指標と総合で")

    print(f"\n[done] outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
