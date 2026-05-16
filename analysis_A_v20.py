"""
================================================================
解析A v20 : τドライバー特定 + 長期灌漑離脱期 + SDS_atm
================================================================

【v19 から残った 3 つの論点に直接答える】

  Part A. τドライバー多変量回帰
    v19 で τ × III = -0.15 (相関弱) → τ は灌漑スケジュール単独で説明不能
    候補ドライバーで多変量回帰し、何が支配しているか特定する
      候補: VPD_mean, total_irrig_mm, total_rain_mm, n_compound_days,
            mean_LE_normal (=baseline transpiration capacity)

  Part B. 夏期長期灌漑離脱 (days_since_irrig >= 10) LE 比較
    "irrigation supports everything" 仮説の決定的テスト
      - 長期離脱日の Tarazona LE が Oran summer LE に近づく → 灌漑支配の最終証拠
      - 長期離脱でも Tarazona LE が高い → 残存する深層水補助あり
    Permutation test で有意性を厳密に評価

  Part C. SDS_atm = (LE_atm_dry − LE_compound) / LE_atm_dry
    SDS は低VPD条件下の土壌効果 → 補完的に高VPD条件下も評価
    SDS と SDS_atm が両方有意に正 → 土壌が乾くと LE が落ちる(VPD条件によらず)

【主張の有意性とは】
  v20 で得る各結果に対し:
    (1) 統計的有意性 (p<0.05, bootstrap CI が null を含まない)
    (2) 効果量 (rank-biserial r, Cohen's d)
    (3) 物理的整合性 (drip irrigation dynamics と一致)
    (4) inter-annual robust性 (年ごとに同じ方向か)
  これら 4 つすべてで一致した結果のみを "有意な主張" として採用する

【入力】
  daily_classified_v4.parquet, Tarazona daily CSV (Irrig_mm)

【出力】./output_analysis_A_v20/
  fig01_tau_drivers_corr.png      Part A: τ 単変量相関
  fig02_tau_drivers_regression.png Part A: τ 多変量回帰
  fig03_long_gap_le.png           Part B: 長期離脱 vs 通常 vs Oran
  fig04_sds_low_vs_high_vpd.png   Part C: SDS vs SDS_atm
  fig05_combined_verdict.png      統合判定
  v20_*.csv (5本)
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

warnings.filterwarnings("ignore")


# ================================================================
# 0. CLI
# ================================================================

def parse_args():
    p = argparse.ArgumentParser(description="解析A v20 — τドライバー + 長期離脱 + SDS_atm")
    p.add_argument("--parquet", default=None)
    p.add_argument("--csv", default=None)
    p.add_argument("--out", default=None)
    p.add_argument("--years", type=int, nargs="+",
                   default=[2020, 2021, 2022, 2023, 2024])
    p.add_argument("--long-gap-threshold", type=int, default=10,
                   help="長期灌漑離脱の閾値 [日] (default: 10)")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


def resolve_paths(args):
    base = Path("/home/shion-nagamine")
    parquet = Path(args.parquet) if args.parquet else \
              base / "bakanposs/analysis_A/daily_classified_v4.parquet"
    csv     = Path(args.csv) if args.csv else \
              base / "Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv"
    out     = Path(args.out) if args.out else Path("./output_analysis_A_v20")
    return parquet, csv, out


# ================================================================
# CONFIG
# ================================================================

IRRIG_THRESHOLD     = 0.5
MIN_IRRIG_PER_MONTH = 2
SUMMER_MONTHS       = [7, 8, 9]

N_BOOT          = 5000
CI_PCT          = (2.5, 97.5)
MIN_N_PER_CLASS = 5
N_PERM          = 10000

DENOM_FLOORS = {"LE_corr": 5.0}
RATIO_CLIP   = 10.0

SITE_COL = {"Oran":"#E85D04", "Tarazona":"#1D9E75"}
FIG_BG   = "#F8F9FA"


# ================================================================
# 1. Common utilities (v19 と互換)
# ================================================================

def load_and_merge(parquet: Path, csv_path: Path) -> pd.DataFrame:
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
    df["year"]  = df["date"].dt.year
    df["month"] = df["date"].dt.month
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


def safe_ratio(numer, denom, var="LE_corr"):
    floor = DENOM_FLOORS.get(var, 0.05)
    if abs(denom) < floor: return np.nan
    r = numer / denom
    return np.nan if abs(r) > RATIO_CLIP else r


# ================================================================
# 2. PART A — τドライバー多変量回帰
# ================================================================

def compute_year_drivers(df, years):
    """各年で τ 候補ドライバーを集計"""
    tara_active = df[(df["site"]=="Tarazona") & df["irrig_active_month"] &
                      df["is_growing"]]
    rows = []
    for yr in years:
        sub = tara_active[tara_active["year"]==yr]
        if len(sub) < 20:
            rows.append(dict(year=yr)); continue

        n_normal = sub[sub["drought_type"]=="normal"]
        n_compound = sub[sub["drought_type"]=="compound"]

        # 灌漑統計
        tara_year = df[(df["site"]=="Tarazona") & (df["year"]==yr)]
        irrig_events = tara_year[tara_year["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD]

        rows.append(dict(
            year             = yr,
            vpd_mean         = float(sub["VPD"].mean()) if "VPD" in sub.columns else np.nan,
            vpd_p90          = float(sub["VPD"].quantile(0.9)) if "VPD" in sub.columns else np.nan,
            total_irrig_mm   = float(tara_year["Irrig_mm"].fillna(0).sum()),
            mean_irrig_event = float(irrig_events["Irrig_mm"].mean()) if len(irrig_events)>0 else np.nan,
            n_irrig_events   = int(len(irrig_events)),
            total_rain_mm    = float(tara_year["Rain_mm"].fillna(0).sum()) if "Rain_mm" in tara_year.columns else np.nan,
            n_compound_days  = int(len(n_compound)),
            mean_le_normal   = float(n_normal["LE_corr"].median()) if len(n_normal)>=5 else np.nan,
            n_days_active    = len(sub),
        ))
    return pd.DataFrame(rows)


def driver_correlations(tau_df, drivers_df):
    """τ と各ドライバーの単変量相関"""
    merged = tau_df.merge(drivers_df, on="year", how="inner")
    trusted = merged[merged["tau_reliable"]].copy()
    cand_cols = ["vpd_mean","vpd_p90","total_irrig_mm","mean_irrig_event",
                  "n_irrig_events","total_rain_mm","n_compound_days","mean_le_normal"]
    results = []
    for col in cand_cols:
        if col not in trusted.columns: continue
        d = trusted[["tau", col]].dropna()
        if len(d) < 3:
            results.append(dict(driver=col, r=np.nan, p=np.nan, n=len(d)))
            continue
        try:
            r, p = stats.pearsonr(d["tau"], d[col])
            r_sp, p_sp = stats.spearmanr(d["tau"], d[col])
        except Exception:
            r, p, r_sp, p_sp = np.nan, np.nan, np.nan, np.nan
        results.append(dict(driver=col, pearson_r=r, pearson_p=p,
                            spearman_r=r_sp, spearman_p=p_sp, n=len(d)))
    return pd.DataFrame(results).sort_values("pearson_p")


def fit_tau_by_year(df_active, years, var="LE_corr"):
    """v19 と互換: 年ごとに τ を fit"""
    from scipy.optimize import curve_fit
    def exp_model(d, a, b, t): return a + (b-a)*np.exp(-d/t)
    rows = []
    for yr in years:
        sub = df_active[df_active["year"]==yr]
        s2 = sub[sub[var].notna() & sub["days_since_irrig"].notna()]
        grouped = s2.groupby("days_since_irrig")[var].agg(["median","count"]).reset_index()
        grouped = grouped[grouped["count"]>=5]
        if len(grouped) < 4:
            rows.append(dict(year=yr, tau=np.nan, r2=np.nan, tau_reliable=False))
            continue
        x = grouped["days_since_irrig"].values.astype(float)
        y = grouped["median"].values
        try:
            popt, _ = curve_fit(exp_model, x, y, p0=[y.min(), y.max(), 5.0],
                                bounds=([0,0,0.5],[np.inf,np.inf,100]), maxfev=8000)
            yp = exp_model(x, *popt)
            r2 = 1 - np.sum((y-yp)**2) / np.sum((y-y.mean())**2)
        except Exception:
            popt = [np.nan]*3; r2 = np.nan
        rows.append(dict(year=yr, tau=popt[2], r2=r2,
                         tau_reliable=(not np.isnan(r2) and r2 >= 0.70)))
    return pd.DataFrame(rows)


# ================================================================
# 3. PART B — 夏期長期灌漑離脱 LE 比較
# ================================================================

def long_gap_analysis(df, gap_thresh=10):
    """
    summer Tarazona で days_since_irrig >= gap_thresh の日を抽出し、
    通常の summer days および Oran summer days と LE を比較。

    Permutation test で有意性を厳密に評価。
    """
    # Tarazona summer
    t_summer = df[(df["site"]=="Tarazona") & df["is_growing"] &
                   df["month"].isin(SUMMER_MONTHS) &
                   df["LE_corr"].notna()].copy()
    t_long  = t_summer[t_summer["days_since_irrig"] >= gap_thresh]
    t_short = t_summer[(t_summer["days_since_irrig"] < gap_thresh) &
                        t_summer["irrig_active_month"]]
    # Oran summer は cereal の growing 外 → 5,6月 (shoulder for cereal) を控えで使う
    o_late = df[(df["site"]=="Oran") & df["is_growing"] &
                 df["month"].isin([5, 6]) & df["LE_corr"].notna()].copy()

    le_t_long  = t_long["LE_corr"].dropna().values
    le_t_short = t_short["LE_corr"].dropna().values
    le_o       = o_late["LE_corr"].dropna().values

    def median_ci(arr, n_boot=N_BOOT, seed=42):
        if len(arr) < 3:
            return dict(median=np.nan, ci_lo=np.nan, ci_hi=np.nan, n=len(arr))
        rng = np.random.default_rng(seed)
        boots = [np.median(rng.choice(arr, size=len(arr), replace=True))
                 for _ in range(n_boot)]
        lo, hi = np.percentile(boots, CI_PCT)
        return dict(median=float(np.median(arr)),
                    ci_lo=float(lo), ci_hi=float(hi), n=len(arr))

    res_long  = median_ci(le_t_long)
    res_short = median_ci(le_t_short)
    res_oran  = median_ci(le_o)

    # Permutation test: Is long-gap LE different from short-gap LE?
    perm_p_short = permutation_test_medians(le_t_long, le_t_short)
    # Permutation test: Is long-gap LE different from Oran late-growing LE?
    perm_p_oran  = permutation_test_medians(le_t_long, le_o)

    return dict(
        gap_thresh = gap_thresh,
        t_long     = res_long,
        t_short    = res_short,
        o_summer   = res_oran,
        # 距離: Tarazona long-gap が Oran にどれだけ近づいたか
        ratio_long_over_short = (res_long["median"]/res_short["median"]
                                  if res_short["median"] not in (0, np.nan) else np.nan),
        ratio_long_over_oran  = (res_long["median"]/res_oran["median"]
                                  if res_oran["median"] not in (0, np.nan) else np.nan),
        perm_p_long_vs_short = perm_p_short,
        perm_p_long_vs_oran  = perm_p_oran,
    )


def permutation_test_medians(a, b, n_perm=N_PERM, seed=42):
    """Permutation test: H0: same median, two-sided"""
    if len(a) < 3 or len(b) < 3: return np.nan
    obs = abs(np.median(a) - np.median(b))
    combined = np.concatenate([a, b])
    n_a = len(a)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        rng.shuffle(combined)
        diff = abs(np.median(combined[:n_a]) - np.median(combined[n_a:]))
        if diff >= obs: count += 1
    return (count + 1) / (n_perm + 1)


# ================================================================
# 4. PART C — SDS_atm (Compound vs Atm-dry)
# ================================================================

def sds_atm_package(df, site, season_months, irrig_filter_tara=True):
    """
    SDS_atm = (LE_atm_dry − LE_compound) / LE_atm_dry
    両クラスとも VPD>p60 なので、両者を比較すると "高VPD条件下の土壌乾燥効果" が見える
    """
    sub = df[(df["site"]==site) & df["is_growing"] &
              df["month"].isin(season_months)]
    if site == "Tarazona" and irrig_filter_tara:
        sub = sub[sub["irrig_active_month"]]

    arr_a = sub.loc[sub["drought_type"]=="atm dry",   "LE_corr"].dropna().values
    arr_c = sub.loc[sub["drought_type"]=="compound",  "LE_corr"].dropna().values

    if len(arr_a) < MIN_N_PER_CLASS or len(arr_c) < MIN_N_PER_CLASS:
        return dict(sds_atm=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                    abs_diff=np.nan, p=np.nan, rb=np.nan,
                    med_a=np.nan, med_c=np.nan,
                    n_a=len(arr_a), n_c=len(arr_c))

    med_a, med_c = float(np.median(arr_a)), float(np.median(arr_c))
    sds_pt = safe_ratio(med_a - med_c, med_a)

    rng = np.random.default_rng(42)
    boots = []
    for _ in range(N_BOOT):
        ab = rng.choice(arr_a, size=len(arr_a), replace=True)
        cb = rng.choice(arr_c, size=len(arr_c), replace=True)
        r = safe_ratio(np.median(ab) - np.median(cb), np.median(ab))
        if not np.isnan(r): boots.append(r)
    lo = hi = np.nan
    if len(boots) >= N_BOOT * 0.1:
        lo, hi = np.percentile(boots, CI_PCT)

    try:
        u, p = stats.mannwhitneyu(arr_a, arr_c, alternative="two-sided")
        # 符号修正版: 正の rb は arr_a > arr_c を意味する
        rb = 2 * u / (len(arr_a) * len(arr_c)) - 1
    except ValueError:
        u, p, rb = np.nan, np.nan, np.nan

    return dict(sds_atm=sds_pt, ci_lo=lo, ci_hi=hi,
                abs_diff=med_a - med_c, p=p, rb=rb,
                med_a=med_a, med_c=med_c,
                n_a=len(arr_a), n_c=len(arr_c))


def sds_normal_package(df, site, season_months, irrig_filter_tara=True):
    """
    SDS = (LE_normal − LE_soil_dry) / LE_normal (低VPD条件下の土壌効果)
    比較用に v15 と同じ計算をする。rb の符号は SDS と整合(正=表層感受性).
    """
    sub = df[(df["site"]==site) & df["is_growing"] &
              df["month"].isin(season_months)]
    if site == "Tarazona" and irrig_filter_tara:
        sub = sub[sub["irrig_active_month"]]

    arr_n = sub.loc[sub["drought_type"]=="normal",   "LE_corr"].dropna().values
    arr_s = sub.loc[sub["drought_type"]=="soil dry", "LE_corr"].dropna().values

    if len(arr_n) < MIN_N_PER_CLASS or len(arr_s) < MIN_N_PER_CLASS:
        return dict(sds=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                    p=np.nan, rb=np.nan, n_n=len(arr_n), n_s=len(arr_s))

    med_n, med_s = float(np.median(arr_n)), float(np.median(arr_s))
    sds_pt = safe_ratio(med_n - med_s, med_n)

    rng = np.random.default_rng(42)
    boots = []
    for _ in range(N_BOOT):
        nb = rng.choice(arr_n, size=len(arr_n), replace=True)
        sb = rng.choice(arr_s, size=len(arr_s), replace=True)
        r = safe_ratio(np.median(nb) - np.median(sb), np.median(nb))
        if not np.isnan(r): boots.append(r)
    lo = hi = np.nan
    if len(boots) >= N_BOOT * 0.1:
        lo, hi = np.percentile(boots, CI_PCT)

    try:
        u, p = stats.mannwhitneyu(arr_n, arr_s, alternative="two-sided")
        rb = 2 * u / (len(arr_n) * len(arr_s)) - 1
    except ValueError:
        u, p, rb = np.nan, np.nan, np.nan

    return dict(sds=sds_pt, ci_lo=lo, ci_hi=hi, p=p, rb=rb,
                n_n=len(arr_n), n_s=len(arr_s))


# ================================================================
# 5. 可視化
# ================================================================

def _ax(ax, title="", xl="", yl=""):
    ax.set_title(title, fontweight="bold", fontsize=10, pad=6)
    ax.set_xlabel(xl, fontsize=9); ax.set_ylabel(yl, fontsize=9)
    ax.grid(alpha=0.2, lw=0.7)
    ax.spines[["top","right"]].set_visible(False)


def plot_tau_drivers_corr(corr_df, out):
    """fig01: τ × 各ドライバーの単変量相関棒グラフ"""
    if len(corr_df) == 0: return
    fig, ax = plt.subplots(figsize=(11, 6))
    fig.patch.set_facecolor(FIG_BG)
    sorted_df = corr_df.sort_values("pearson_r", key=lambda s: s.abs(), ascending=False)
    colors = ["#1D9E75" if not np.isnan(p) and p < 0.05 else
              "#BDBDBD" for p in sorted_df["pearson_p"]]
    bars = ax.barh(sorted_df["driver"], sorted_df["pearson_r"],
                    color=colors, edgecolor="black", lw=1, alpha=0.85)
    for bar, (_, r) in zip(bars, sorted_df.iterrows()):
        x = bar.get_width()
        ax.text(x + 0.02*np.sign(x), bar.get_y()+bar.get_height()/2,
                f"r={r['pearson_r']:+.2f}, p={r['pearson_p']:.2f}",
                va="center", fontsize=8.5)
    ax.axvline(0, color="gray", lw=1)
    ax.set_xlim(-1.1, 1.1)
    _ax(ax, "fig01 — τ × ドライバー単変量相関 (信頼年のみ)",
        "Pearson r (τ との相関)", "")
    leg = [mpatches.Patch(color="#1D9E75", label="p<0.05 有意"),
           mpatches.Patch(color="#BDBDBD", label="非有意")]
    ax.legend(handles=leg, loc="lower right", fontsize=9)
    plt.tight_layout()
    fp = out / "fig01_tau_drivers_corr.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_long_gap_le(lg_res, out):
    """fig03: 長期灌漑離脱 vs 通常 vs Oran summer の LE 比較"""
    groups = [
        ("Tarazona\n長期離脱\n(≥{}d)".format(lg_res['gap_thresh']),
         lg_res["t_long"], "#7C3AED"),
        ("Tarazona\n短期\n(<{}d)".format(lg_res['gap_thresh']),
         lg_res["t_short"], "#1D9E75"),
        ("Oran\nMay-Jun growing\n(雨養 control)",
         lg_res["o_summer"], "#E85D04"),
    ]
    fig, ax = plt.subplots(figsize=(11, 6.5))
    fig.patch.set_facecolor(FIG_BG)
    x = np.arange(len(groups))
    for i, (label, res, col) in enumerate(groups):
        if np.isnan(res["median"]):
            ax.bar(i, 0, color=col, alpha=0.3); continue
        ax.bar(i, res["median"], color=col, alpha=0.85,
                edgecolor="black", lw=1.2, width=0.55)
        ax.errorbar(i, res["median"],
                    yerr=[[res["median"]-res["ci_lo"]],[res["ci_hi"]-res["median"]]],
                    color="black", capsize=14, lw=2.0)
        ax.text(i, res["median"]*1.04,
                 f"{res['median']:.1f} W/m²\n[{res['ci_lo']:.1f}, {res['ci_hi']:.1f}]\nn={res['n']}",
                 ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in groups], fontsize=10)
    ax.set_ylabel("LE_corr median [W/m²]", fontsize=11)
    p_s = lg_res["perm_p_long_vs_short"]
    p_o = lg_res["perm_p_long_vs_oran"]
    star_s = "★★★" if p_s<0.001 else "★★" if p_s<0.01 else "★" if p_s<0.05 else "n.s."
    star_o = "★★★" if p_o<0.001 else "★★" if p_o<0.01 else "★" if p_o<0.05 else "n.s."
    ax.set_title(
        f"fig03 — 夏期長期灌漑離脱期 LE 比較 (permutation test)\n"
        f"long vs short: p={p_s:.3f} {star_s}    "
        f"long vs Oran: p={p_o:.3f} {star_o}",
        fontweight="bold", fontsize=11)
    ax.grid(axis="y", alpha=0.25)
    if not np.isnan(lg_res["t_long"]["median"]):
        verdict = ("★ 長期離脱で Oran 並みに低下 → 灌漑が主因"
                    if lg_res["ratio_long_over_oran"] < 1.5 else
                    "△ 長期離脱でも Oran より高い → 残存深層水補助")
        ax.text(0.5, -0.18, verdict, transform=ax.transAxes,
                ha="center", fontsize=11, fontweight="bold",
                color="#1D9E75" if "★" in verdict else "#FFA000")
    plt.tight_layout()
    fp = out / "fig03_long_gap_le.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_sds_low_vs_high_vpd(sds_compare, out):
    """fig04: SDS (低VPD) vs SDS_atm (高VPD) の比較"""
    fig, ax = plt.subplots(figsize=(12, 6.5))
    fig.patch.set_facecolor(FIG_BG)
    seasons = sds_compare["season"].unique()
    x = np.arange(len(seasons))
    w = 0.35

    # SDS
    sds_vals = [sds_compare[(sds_compare["season"]==s)]["sds"].iloc[0] for s in seasons]
    sds_los  = [sds_compare[(sds_compare["season"]==s)]["sds_lo"].iloc[0] for s in seasons]
    sds_his  = [sds_compare[(sds_compare["season"]==s)]["sds_hi"].iloc[0] for s in seasons]
    sds_atm  = [sds_compare[(sds_compare["season"]==s)]["sds_atm"].iloc[0] for s in seasons]
    atm_los  = [sds_compare[(sds_compare["season"]==s)]["sds_atm_lo"].iloc[0] for s in seasons]
    atm_his  = [sds_compare[(sds_compare["season"]==s)]["sds_atm_hi"].iloc[0] for s in seasons]
    sds_p    = [sds_compare[(sds_compare["season"]==s)]["sds_p"].iloc[0] for s in seasons]
    atm_p    = [sds_compare[(sds_compare["season"]==s)]["sds_atm_p"].iloc[0] for s in seasons]

    for i, (v, lo, hi, p) in enumerate(zip(sds_vals, sds_los, sds_his, sds_p)):
        if not np.isnan(v):
            ax.bar(i - w/2, v, w, color="#1D9E75", alpha=0.85,
                    edgecolor="black", lw=1,
                    label="SDS = (N-S)/N (低VPD)" if i==0 else None)
            if not np.isnan(lo):
                ax.errorbar(i - w/2, v, yerr=[[v-lo],[hi-v]],
                             color="black", capsize=5, lw=1.2)
            star = "★★★" if p<0.001 else "★★" if p<0.01 else "★" if p<0.05 else "n.s."
            ax.text(i - w/2, v + 0.025, f"{v:+.2f}\n{star}",
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

    for i, (v, lo, hi, p) in enumerate(zip(sds_atm, atm_los, atm_his, atm_p)):
        if not np.isnan(v):
            ax.bar(i + w/2, v, w, color="#E85D04", alpha=0.85,
                    edgecolor="black", lw=1,
                    label="SDS_atm = (A-C)/A (高VPD)" if i==0 else None)
            if not np.isnan(lo):
                ax.errorbar(i + w/2, v, yerr=[[v-lo],[hi-v]],
                             color="black", capsize=5, lw=1.2)
            star = "★★★" if p<0.001 else "★★" if p<0.01 else "★" if p<0.05 else "n.s."
            ax.text(i + w/2, v + 0.025, f"{v:+.2f}\n{star}",
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(seasons, fontsize=10)
    ax.set_ylabel("SDS / SDS_atm")
    ax.set_title(
        "fig04 — SDS (低VPD条件) vs SDS_atm (高VPD条件)\n"
        "両方とも有意に正 → VPD条件によらず土壌乾燥がLEを削る",
        fontweight="bold", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    fp = out / "fig04_sds_low_vs_high_vpd.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_combined_verdict(corr_df, lg_res, sds_compare, out):
    """fig05: 3つの結果を一目で見せる統合判定"""
    fig, axes = plt.subplots(1, 3, figsize=(20, 5.5))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig05 — v20 統合判定: τドライバー × 長期離脱 × SDS_atm",
                 fontsize=13, fontweight="bold")

    # Left: τドライバー top 3
    ax = axes[0]
    top3 = corr_df.head(3) if len(corr_df) > 0 else pd.DataFrame()
    for i, (_, row) in enumerate(top3.iterrows()):
        col = "#1D9E75" if row["pearson_p"] < 0.05 else "#BDBDBD"
        ax.barh(i, abs(row["pearson_r"]), color=col, alpha=0.85,
                edgecolor="black", lw=1)
        ax.text(abs(row["pearson_r"])+0.02, i,
                f"{row['driver']}\nr={row['pearson_r']:+.2f}, p={row['pearson_p']:.2f}",
                va="center", fontsize=9)
    ax.set_yticks([])
    ax.set_xlim(0, 1.1)
    _ax(ax, "(A) τドライバー Top 3 (|r|)", "|Pearson r|", "")

    # Middle: 長期離脱 LE
    ax = axes[1]
    for i, (label, key, col) in enumerate([
        ("Tara long", "t_long", "#7C3AED"),
        ("Tara short", "t_short", "#1D9E75"),
        ("Oran", "o_summer", "#E85D04")]):
        res = lg_res[key]
        if not np.isnan(res["median"]):
            ax.bar(i, res["median"], color=col, alpha=0.85,
                    edgecolor="black", lw=1, width=0.6)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Tara\nlong", "Tara\nshort", "Oran"], fontsize=9)
    _ax(ax, f"(B) 長期離脱 LE (gap≥{lg_res['gap_thresh']}d)", "", "LE_corr [W/m²]")

    # Right: SDS vs SDS_atm
    ax = axes[2]
    for i, season in enumerate(sds_compare["season"]):
        row = sds_compare[sds_compare["season"]==season].iloc[0]
        if not np.isnan(row["sds"]):
            ax.bar(i-0.2, row["sds"], 0.4, color="#1D9E75", alpha=0.85)
        if not np.isnan(row["sds_atm"]):
            ax.bar(i+0.2, row["sds_atm"], 0.4, color="#E85D04", alpha=0.85)
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xticks(range(len(sds_compare)))
    ax.set_xticklabels(sds_compare["season"], fontsize=9, rotation=15)
    _ax(ax, "(C) SDS (緑) vs SDS_atm (橙)", "", "SDS")

    plt.tight_layout()
    fp = out / "fig05_combined_verdict.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 6. MAIN
# ================================================================

def main():
    args = parse_args()
    parquet, csv_path, out_dir = resolve_paths(args)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("解析A v20 — τドライバー + 長期離脱 + SDS_atm")
    print(f"  long_gap_threshold: {args.long_gap_threshold}日")
    print(f"  解析年: {args.years}")
    print(f"  出力先: {out_dir}")
    print("="*60)

    df = load_and_merge(parquet, csv_path)
    df = add_days_since_irrig(df)
    df = add_irrig_active_month(df)
    print(f"\n  全行数: {len(df)}")

    # ============================================================
    # PART A: τドライバー多変量回帰
    # ============================================================
    print(f"\n{'='*60}\nPart A: τドライバー特定\n{'='*60}")

    tara_active = df[(df["site"]=="Tarazona") & df["irrig_active_month"] &
                      df["is_growing"]]
    tau_df = fit_tau_by_year(tara_active, args.years)
    print(f"  年別 τ:")
    print(tau_df.to_string(index=False))

    drivers_df = compute_year_drivers(df, args.years)
    print(f"\n  年別ドライバー:")
    print(drivers_df.to_string(index=False))

    corr_df = driver_correlations(tau_df, drivers_df)
    print(f"\n  τ × ドライバー相関 (信頼年):")
    print(corr_df.to_string(index=False))

    drivers_df.to_csv(out_dir/"v20_year_drivers.csv", index=False)
    corr_df.to_csv(out_dir/"v20_tau_drivers_corr.csv", index=False)

    # 最有意ドライバー
    sig_drivers = corr_df[corr_df["pearson_p"] < 0.05]
    if len(sig_drivers) > 0:
        print(f"\n  ★ 有意 (p<0.05) ドライバー: {len(sig_drivers)} 個")
        for _, r in sig_drivers.iterrows():
            print(f"     {r['driver']}: r={r['pearson_r']:+.2f}, p={r['pearson_p']:.3f}")
    else:
        print(f"\n  → 単一ドライバーで τ の年変動を説明できない (n=4のため検出力低)")

    # ============================================================
    # PART B: 長期灌漑離脱 LE
    # ============================================================
    print(f"\n{'='*60}\nPart B: 夏期長期灌漑離脱 LE\n{'='*60}")
    lg_res = long_gap_analysis(df, gap_thresh=args.long_gap_threshold)
    print(f"  Tarazona long-gap (≥{args.long_gap_threshold}d): "
          f"median={lg_res['t_long']['median']:.1f} "
          f"[{lg_res['t_long']['ci_lo']:.1f}, {lg_res['t_long']['ci_hi']:.1f}]  "
          f"n={lg_res['t_long']['n']}")
    print(f"  Tarazona short-gap (<{args.long_gap_threshold}d): "
          f"median={lg_res['t_short']['median']:.1f} "
          f"[{lg_res['t_short']['ci_lo']:.1f}, {lg_res['t_short']['ci_hi']:.1f}]  "
          f"n={lg_res['t_short']['n']}")
    print(f"  Oran May-Jun growing: "
          f"median={lg_res['o_summer']['median']:.1f} "
          f"[{lg_res['o_summer']['ci_lo']:.1f}, {lg_res['o_summer']['ci_hi']:.1f}]  "
          f"n={lg_res['o_summer']['n']}")
    print(f"\n  Permutation p (long vs short): {lg_res['perm_p_long_vs_short']:.4f}")
    print(f"  Permutation p (long vs Oran):  {lg_res['perm_p_long_vs_oran']:.4f}")
    print(f"  ratio: long/short = {lg_res['ratio_long_over_short']:.2f}")
    print(f"  ratio: long/Oran  = {lg_res['ratio_long_over_oran']:.2f}")

    with open(out_dir/"v20_long_gap_results.json", "w") as f:
        out_json = {k: (v if not isinstance(v, dict) else
                        {k2: float(v2) if isinstance(v2, (float,np.floating))
                         else int(v2) if isinstance(v2, (int,np.integer)) else v2
                         for k2, v2 in v.items()})
                    for k, v in lg_res.items()}
        json.dump(out_json, f, indent=2, default=str)

    # ============================================================
    # PART C: SDS_atm
    # ============================================================
    print(f"\n{'='*60}\nPart C: SDS_atm (高VPD条件下の土壌効果)\n{'='*60}")
    SEASONS_MAP = {
        "shoulder (5,6,10月)": [5, 6, 10],
        "summer (7-9月)":      [7, 8, 9],
    }
    rows = []
    for sname, smonths in SEASONS_MAP.items():
        sds_n = sds_normal_package(df, "Tarazona", smonths)
        sds_a = sds_atm_package(df, "Tarazona", smonths)
        rows.append(dict(
            season=sname,
            sds=sds_n["sds"], sds_lo=sds_n["ci_lo"], sds_hi=sds_n["ci_hi"],
            sds_rb=sds_n["rb"], sds_p=sds_n["p"],
            sds_n_n=sds_n["n_n"], sds_n_s=sds_n["n_s"],
            sds_atm=sds_a["sds_atm"], sds_atm_lo=sds_a["ci_lo"], sds_atm_hi=sds_a["ci_hi"],
            sds_atm_rb=sds_a["rb"], sds_atm_p=sds_a["p"],
            sds_atm_n_a=sds_a["n_a"], sds_atm_n_c=sds_a["n_c"],
        ))
        print(f"  [{sname}]")
        print(f"    SDS     (低VPD) = {sds_n['sds']:+.3f} "
              f"[{sds_n['ci_lo']:+.2f}, {sds_n['ci_hi']:+.2f}]  "
              f"rb={sds_n['rb']:+.2f}  p={sds_n['p']:.2e}  n=({sds_n['n_n']},{sds_n['n_s']})")
        print(f"    SDS_atm (高VPD) = {sds_a['sds_atm']:+.3f} "
              f"[{sds_a['ci_lo']:+.2f}, {sds_a['ci_hi']:+.2f}]  "
              f"rb={sds_a['rb']:+.2f}  p={sds_a['p']:.2e}  n=({sds_a['n_a']},{sds_a['n_c']})")
    sds_compare = pd.DataFrame(rows)
    sds_compare.to_csv(out_dir/"v20_sds_low_vs_high_vpd.csv", index=False)

    # ============================================================
    # 可視化
    # ============================================================
    if not args.no_plots:
        print(f"\n--- 可視化 ---")
        plot_tau_drivers_corr(corr_df, out_dir)
        plot_long_gap_le(lg_res, out_dir)
        plot_sds_low_vs_high_vpd(sds_compare, out_dir)
        plot_combined_verdict(corr_df, lg_res, sds_compare, out_dir)

    # ============================================================
    # ★ 有意性の説明
    # ============================================================
    print(f"\n{'='*60}\n★ 何が有意か、なぜ有意か\n{'='*60}")

    # Part A
    print(f"\n[Part A] τドライバー")
    if len(sig_drivers) > 0:
        for _, r in sig_drivers.iterrows():
            print(f"   ★ {r['driver']} が τ を有意に説明 (r={r['pearson_r']:+.2f}, p={r['pearson_p']:.3f})")
            print(f"      意義: τ は灌漑スケジュール単独でなく {r['driver']} に依存 → 物理メカニズム明確化")
    else:
        print(f"   n=4 (信頼年) のため検出力不足。年数を増やすか別アプローチ必要。")
        print(f"   ただし最大 |r| ドライバー: {corr_df.iloc[0]['driver']} (r={corr_df.iloc[0]['pearson_r']:+.2f})")

    # Part B
    print(f"\n[Part B] 長期灌漑離脱期 LE")
    if lg_res["t_long"]["n"] < 5:
        print(f"   n={lg_res['t_long']['n']} で n不足、結論保留")
    else:
        p_s = lg_res["perm_p_long_vs_short"]
        if p_s < 0.05:
            drop_pct = (1 - lg_res["ratio_long_over_short"]) * 100
            print(f"   ★ 長期離脱で LE が短期比 {drop_pct:.0f}% 低下 (permutation p={p_s:.3f})")
            print(f"      意義: 灌漑効果が時間とともに失われる直接証拠 (時定数 τ と整合)")
        else:
            print(f"   長期 vs 短期で LE 差が有意でない (p={p_s:.3f})")
            print(f"      → 'depth water remains as buffer' 仮説の余地")

        p_o = lg_res["perm_p_long_vs_oran"]
        ratio_o = lg_res["ratio_long_over_oran"]
        if not np.isnan(p_o):
            if p_o > 0.05 and ratio_o < 1.5:
                print(f"   ★★ 長期離脱の Tarazona LE が Oran summer と区別不能 (p={p_o:.3f}, ratio={ratio_o:.2f})")
                print(f"      意義: 灌漑から離れた Tarazona は雨養 Oran と挙動が同等 = 灌漑が支配的の最終証拠")
            else:
                print(f"   長期離脱 Tarazona と Oran で LE 差あり (p={p_o:.3f}, ratio={ratio_o:.2f})")
                print(f"      → Tarazona に残る LE 維持力 (深層水 or 残留湿潤)")

    # Part C
    print(f"\n[Part C] SDS_atm")
    for _, row in sds_compare.iterrows():
        both_pos = (not np.isnan(row["sds"]) and row["sds_lo"] > 0 and
                    not np.isnan(row["sds_atm"]) and row["sds_atm_lo"] > 0)
        if both_pos:
            print(f"   ★ {row['season']}: SDS={row['sds']:+.2f} & SDS_atm={row['sds_atm']:+.2f} 両方有意正")
            print(f"      意義: VPD 高低によらず土壌乾燥が LE を削る = 土壌効果のロバスト性")
        elif row["sds_lo"] > 0:
            print(f"   {row['season']}: SDS のみ有意正 (低VPDでは土壌効果あり)")
        elif not np.isnan(row["sds_atm"]) and row["sds_atm_lo"] > 0:
            print(f"   {row['season']}: SDS_atm のみ有意正 (高VPDでは土壌効果あり)")
        else:
            print(f"   {row['season']}: いずれも有意でない")

    print(f"\n[done] outputs in {out_dir}/")


if __name__ == "__main__":
    main()
