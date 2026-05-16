"""
================================================================
解析A v26 : Reviewer counter-attack — τ universality の徹底検証
================================================================

【v25 で出た結果】
  Tarazona τ = 2.51 d, Oran τ = 2.51 d (完全一致)
  → 「τ は climate 普遍定数」と reframe する誘惑が出るが、
     reviewer 指摘:
     "Same fitted τ does NOT imply same mechanism"
     - identifiability 問題
     - phenology/season/event size の cofounding
     - amplitude と τ の非独立性

【v26 で行う4つの追加検証】
  Part 1. Oran 季節別 τ stratification
    rain events を冬/春/夏/秋に分割、各季節で τ を別々に fit
    - 全季節で τ ≈ 2.5d → 真の普遍性
    - 季節で τ が大きく違う → "平均の偶然" → 棄却

  Part 2. Phenology-matched comparison
    Oran spring (Mar-May, 植物 active) ↔ Tarazona summer (Jun-Sep, almond active)
    両者とも植物 active = 物理的に同等の蒸散状態
    そこで τ を比較

  Part 3. Pulse-size scaling
    各 event を pulse size (mm) で分類
    pulse size と τ の相関 → 0 なら τ は event size 非依存(真の system 特性)
    相関あり → τ は emergent な fit 値で universal constant ではない

  Part 4. τ × amplitude independence test
    per-event 個別 fit で τ と amplitude を取り、Pearson r を見る
    - r ≈ 0 → τ と amplitude が独立 → 両者の解釈が分離可能
    - r ≠ 0 → fitted exponential の identifiability 問題、解釈に注意

【結論の framing】
  3つの test 全部で τ が 2.5d 近傍 → 'universal climatological τ' 主張可能
  どれか1つでも逸脱 → 'effective relaxation timescale (regime-dependent)' に格下げ

【入力】 v4 parquet + Tarazona daily CSV + Oran daily MASTER
【出力】./output_analysis_A_v26/
  v26_oran_seasonal_tau.csv         Part 1
  v26_phenology_matched_tau.csv     Part 2
  v26_pulse_size_tau.csv            Part 3
  v26_per_event_tau_amplitude.csv   Part 4
  fig01-04 + final verdict figure
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
from scipy import stats
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")


# ================================================================
# CONFIG
# ================================================================

def parse_args():
    p = argparse.ArgumentParser(description="解析A v26 — Reviewer attack: stratified τ")
    p.add_argument("--parquet", default=None)
    p.add_argument("--tara-csv", default=None)
    p.add_argument("--oran-daily", default=None)
    p.add_argument("--out", default=None)
    p.add_argument("--rain-threshold", type=float, default=3.0)
    p.add_argument("--min-window", type=int, default=4)
    p.add_argument("--max-window", type=int, default=14)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


def resolve_paths(args):
    base = Path("/home/shion-nagamine")
    parquet = Path(args.parquet) if args.parquet else \
              base / "bakanposs/analysis_A/daily_classified_v4.parquet"
    tara_csv = Path(args.tara_csv) if args.tara_csv else \
              base / "Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv"
    oran_daily = Path(args.oran_daily) if args.oran_daily else \
              base / "Dataset/Eddy data in Spain/Oran_EddyDaily_MASTER_2018_2020_correct"
    out = Path(args.out) if args.out else Path("./output_analysis_A_v26")
    return parquet, tara_csv, oran_daily, out


IRRIG_THRESHOLD = 0.5
CI_PCT = (2.5, 97.5)
FIG_BG = "#F8F9FA"
TAU_FLOOR = 0.5
TAU_CEIL = 60.0


# ================================================================
# Step 1: 共通: モデル + fit
# ================================================================

def exp_model_2p(d, le0, tau, le_inf_fixed):
    return le_inf_fixed + (le0 - le_inf_fixed) * np.exp(-d/tau)


def fit_pooled_tau(arr_d, arr_y, le_inf_fixed, min_per_bin=3, verbose=False):
    """pooled events のビン中央値で 2-param exp fit"""
    df = pd.DataFrame({"d": arr_d, "y": arr_y})
    g = df.groupby("d")["y"].agg(["median", "count"]).reset_index()
    g = g[g["count"] >= min_per_bin]
    if len(g) < 3: return None
    x = g["d"].values.astype(float)
    y = g["median"].values
    le_0_floor = max(le_inf_fixed + 1, 10.0)
    le_0_ceil = max(y.max() * 1.5, 350)
    def model(d, le0, tau): return exp_model_2p(d, le0, tau, le_inf_fixed)
    for p0 in [[y.max(), 5.0], [y.max() * 1.1, 3.0],
                [(y.max() + le_inf_fixed) / 2, 7.0]]:
        try:
            popt, _ = curve_fit(model, x, y, p0=p0,
                                bounds=([le_0_floor, TAU_FLOOR], [le_0_ceil, TAU_CEIL]),
                                maxfev=8000)
            yp = model(x, *popt)
            sse = float(np.sum((y - yp)**2))
            r2 = float(1 - sse / np.sum((y - y.mean())**2)) if np.var(y) > 0 else np.nan
            return dict(le0=popt[0], tau=popt[1], r2=r2, sse=sse,
                        n_bins=len(g), amplitude=popt[0] - le_inf_fixed)
        except Exception:
            continue
    return None


def bootstrap_tau(arr_d, arr_y, le_inf, n_boot=5000, min_per_bin=3, seed=42):
    rng = np.random.default_rng(seed)
    n = len(arr_d)
    taus = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        res = fit_pooled_tau(arr_d[idx], arr_y[idx], le_inf, min_per_bin)
        if res is None: continue
        taus.append(res["tau"])
    return np.array(taus)


def fit_single_event_tau(event_data, le_inf):
    """単一イベント内の生データで fit (per-event 解析用)"""
    d = event_data["days_since_event"].values.astype(float)
    y = event_data["LE"].values
    valid = ~np.isnan(y)
    d = d[valid]; y = y[valid]
    if len(d) < 4 or np.var(y) < 1.0: return None
    le_0_floor = max(le_inf + 1, 10.0)
    le_0_ceil = max(y.max() * 1.5, 350)
    def model(dd, le0, tau): return exp_model_2p(dd, le0, tau, le_inf)
    try:
        popt, _ = curve_fit(model, d, y, p0=[y.max(), 5.0],
                            bounds=([le_0_floor, TAU_FLOOR], [le_0_ceil, TAU_CEIL]),
                            maxfev=5000)
        yp = model(d, *popt)
        sse = float(np.sum((y-yp)**2))
        r2 = float(1 - sse/np.sum((y-y.mean())**2)) if np.var(y) > 0 else np.nan
        return dict(le0=popt[0], tau=popt[1], r2=r2,
                    amplitude=popt[0] - le_inf, n_obs=len(d))
    except Exception:
        return None


# ================================================================
# Step 2: 共通: data loaders
# ================================================================

def load_tarazona_data(parquet, tara_csv):
    daily = pd.read_parquet(parquet)
    daily["date"] = pd.to_datetime(daily["date"])
    raw = pd.read_csv(tara_csv)
    raw["date"] = pd.to_datetime(raw["date"], errors="coerce")
    irrig_cols = [c for c in ["Irrig_mm","Rain_mm","IrrigRain_mm"] if c in raw.columns]
    extra = raw[["date"] + irrig_cols].dropna(subset=["date"]).drop_duplicates("date")
    return daily[daily["site"]=="Tarazona"].merge(extra, on="date", how="left")


def load_oran_daily_master(path_no_ext):
    candidates = [path_no_ext.with_suffix(s) for s in ["", ".csv", ".xlsx", ".xls"]]
    fp = next((p for p in candidates if p.exists()), None)
    if fp is None: sys.exit(f"Oran daily not found near {path_no_ext}")
    if fp.suffix == ".csv":
        df = pd.read_csv(fp)
    else:
        df = pd.read_excel(fp)
    # detect columns
    date_col = next((c for c in ["Date","date","DATE","TIMESTAMP"] if c in df.columns), None)
    et_col   = next((c for c in ["Eddy_ET_mm_d","ET_mm_d","ET"] if c in df.columns), None)
    rain_col = next((c for c in ["Rain_mm","P_mm","Precip_mm"] if c in df.columns), None)
    nobs_col = next((c for c in ["n_obs_le"] if c in df.columns), None)
    if not all([date_col, et_col, rain_col]):
        sys.exit(f"Required cols missing: date={date_col}, ET={et_col}, rain={rain_col}")
    out = pd.DataFrame()
    out["date"] = pd.to_datetime(df[date_col], errors="coerce")
    out["LE"] = pd.to_numeric(df[et_col], errors="coerce") * 28.36  # mm/d → W/m²
    out["Rain_mm"] = pd.to_numeric(df[rain_col], errors="coerce")
    if nobs_col:
        n_obs = pd.to_numeric(df[nobs_col], errors="coerce")
        out.loc[n_obs.fillna(0) < 24, "LE"] = np.nan
    return out.dropna(subset=["date"]).drop_duplicates("date").sort_values("date").reset_index(drop=True)


# ================================================================
# Step 3: Event detection (両サイト共通)
# ================================================================

def detect_events(df_daily, water_col, threshold, min_window=4, max_window=14,
                  exclude_other_water_col=None):
    """
    water events 抽出。
    exclude_other_water_col: 他の水入力(灌漑など)が周辺にある event を除外
    """
    df = df_daily.sort_values("date").reset_index(drop=True).copy()
    event_dates = df[df[water_col].fillna(0) > threshold]["date"].values
    events = []
    for i, ed in enumerate(event_dates):
        ed = pd.Timestamp(ed)
        if i + 1 < len(event_dates):
            next_ed = pd.Timestamp(event_dates[i+1])
            window_end = min(next_ed - pd.Timedelta(days=1),
                             ed + pd.Timedelta(days=max_window))
        else:
            window_end = ed + pd.Timedelta(days=max_window)
        if (window_end - ed).days < min_window: continue
        evd = df[(df["date"] >= ed) & (df["date"] <= window_end)].copy()
        evd["days_since_event"] = (evd["date"] - ed).dt.days
        evd = evd[evd["LE"].notna()]
        if len(evd) < min_window: continue
        # exclude_other check
        if exclude_other_water_col and exclude_other_water_col in df.columns:
            # 周辺 ±14日に other_water > threshold あれば除外
            past = df[(df["date"] >= ed - pd.Timedelta(days=14)) & (df["date"] < ed)]
            future = df[(df["date"] > ed + pd.Timedelta(days=max_window)) &
                         (df["date"] <= ed + pd.Timedelta(days=max_window + 14))]
            if (past[exclude_other_water_col].fillna(0).max() > threshold or
                future[exclude_other_water_col].fillna(0).max() > threshold):
                continue
        events.append(dict(event_start=ed, data=evd,
                           pulse_mm=float(df.loc[df["date"]==ed, water_col].iloc[0])))
    return events


def estimate_le_inf(events, tail_threshold=7):
    """events を pool して d >= tail_threshold の median を LE_inf とする"""
    pooled = pd.concat([ev["data"] for ev in events])
    mask = pooled["days_since_event"] >= tail_threshold
    if mask.sum() < 5:
        return float(np.percentile(pooled["LE"], 25))
    return float(np.median(pooled.loc[mask, "LE"]))


# ================================================================
# PART 1: Oran 季節別 τ
# ================================================================

SEASONS = {
    "winter (Nov-Feb)" : [11, 12, 1, 2],
    "spring (Mar-May)" : [3, 4, 5],
    "summer (Jun-Aug)" : [6, 7, 8],
    "autumn (Sep-Oct)" : [9, 10],
}


def stratify_events_by_season(events):
    out = {s: [] for s in SEASONS}
    for ev in events:
        m = ev["event_start"].month
        for s, ms in SEASONS.items():
            if m in ms:
                out[s].append(ev); break
    return out


def compute_seasonal_taus(events_by_season, n_boot=5000):
    rows = []
    for season, evs in events_by_season.items():
        if len(evs) < 2:
            rows.append(dict(season=season, n_events=len(evs),
                             tau=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                             le_inf=np.nan, r2=np.nan,
                             reason="n<2 events"))
            continue
        pooled = pd.concat([ev["data"] for ev in evs])
        arr_d = pooled["days_since_event"].values.astype(float)
        arr_y = pooled["LE"].values
        le_inf = estimate_le_inf(evs)
        res = fit_pooled_tau(arr_d, arr_y, le_inf)
        if res is None:
            rows.append(dict(season=season, n_events=len(evs), n_data=len(pooled),
                             tau=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                             le_inf=le_inf, r2=np.nan, reason="fit failed"))
            continue
        boot_taus = bootstrap_tau(arr_d, arr_y, le_inf, n_boot=n_boot)
        if len(boot_taus) > 100:
            lo, hi = np.percentile(boot_taus, CI_PCT)
        else:
            lo, hi = np.nan, np.nan
        rows.append(dict(season=season, n_events=len(evs), n_data=len(pooled),
                         tau=res["tau"], ci_lo=lo, ci_hi=hi,
                         le_inf=le_inf, le0=res["le0"],
                         amplitude=res["amplitude"], r2=res["r2"]))
    return pd.DataFrame(rows)


# ================================================================
# PART 2: Phenology-matched comparison
# ================================================================

def phenology_matched(oran_events, tara_events, n_boot=5000):
    """
    Oran spring (Mar-May, cereal vegetative/reproductive) ↔
    Tarazona summer (Jun-Sep, almond peak transpiration)

    両者とも植物が "active growing" 状態 → physiological matching
    """
    oran_active = [ev for ev in oran_events if ev["event_start"].month in [3, 4, 5]]
    tara_active = [ev for ev in tara_events if ev["event_start"].month in [6, 7, 8, 9]]

    rows = []
    for label, evs in [("Oran spring (Mar-May, active cereal)", oran_active),
                        ("Tarazona summer (Jun-Sep, active almond)", tara_active)]:
        if len(evs) < 2:
            rows.append(dict(group=label, n_events=len(evs),
                             tau=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                             reason="n<2 events"))
            continue
        pooled = pd.concat([ev["data"] for ev in evs])
        arr_d = pooled["days_since_event"].values.astype(float)
        arr_y = pooled["LE"].values
        le_inf = estimate_le_inf(evs)
        res = fit_pooled_tau(arr_d, arr_y, le_inf)
        if res is None:
            rows.append(dict(group=label, n_events=len(evs),
                             tau=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                             reason="fit failed"))
            continue
        boot_taus = bootstrap_tau(arr_d, arr_y, le_inf, n_boot=n_boot)
        lo, hi = (np.percentile(boot_taus, CI_PCT)
                   if len(boot_taus) > 100 else (np.nan, np.nan))
        rows.append(dict(group=label, n_events=len(evs), n_data=len(pooled),
                         tau=res["tau"], ci_lo=lo, ci_hi=hi,
                         le_inf=le_inf, le0=res["le0"],
                         amplitude=res["amplitude"], r2=res["r2"]))
    return pd.DataFrame(rows)


# ================================================================
# PART 3: Per-event individual τ + pulse-size scaling
# ================================================================

def per_event_fits(events, site_name, le_inf):
    rows = []
    for ev in events:
        res = fit_single_event_tau(ev["data"], le_inf)
        if res is None:
            rows.append(dict(site=site_name, event_start=ev["event_start"],
                             pulse_mm=ev["pulse_mm"], month=ev["event_start"].month,
                             tau=np.nan, le0=np.nan, amplitude=np.nan, r2=np.nan,
                             n_obs=np.nan))
            continue
        rows.append(dict(site=site_name, event_start=ev["event_start"],
                         pulse_mm=ev["pulse_mm"], month=ev["event_start"].month,
                         **res))
    return pd.DataFrame(rows)


def pulse_size_test(per_event_df, r2_min=0.3):
    """pulse_mm と τ の相関 + amplitude vs τ"""
    good = per_event_df.dropna(subset=["tau","pulse_mm","amplitude"])
    good = good[good["r2"] >= r2_min]
    if len(good) < 5:
        return None
    r_pulse, p_pulse = stats.pearsonr(good["pulse_mm"], good["tau"])
    r_amp, p_amp = stats.pearsonr(good["amplitude"], good["tau"])
    return dict(n=len(good),
                pulse_tau_r=r_pulse, pulse_tau_p=p_pulse,
                amplitude_tau_r=r_amp, amplitude_tau_p=p_amp)


# ================================================================
# 可視化
# ================================================================

def _ax(ax, title="", xl="", yl=""):
    ax.set_title(title, fontweight="bold", fontsize=11, pad=6)
    ax.set_xlabel(xl, fontsize=10); ax.set_ylabel(yl, fontsize=10)
    ax.grid(alpha=0.2, lw=0.7); ax.spines[["top","right"]].set_visible(False)


def plot_seasonal_taus(seasonal_df, pooled_tau, out, title_prefix="Oran"):
    fig, ax = plt.subplots(figsize=(12, 6.5))
    fig.patch.set_facecolor(FIG_BG)
    df = seasonal_df.copy()
    cols = ["#1976D2", "#43A047", "#FB8C00", "#7C3AED"]
    x = np.arange(len(df))
    for i, (_, row) in enumerate(df.iterrows()):
        if np.isnan(row["tau"]):
            ax.text(i, 1, f"× {row.get('reason','')}", ha="center", fontsize=9)
            continue
        ax.bar(i, row["tau"], color=cols[i], alpha=0.85, edgecolor="black", lw=1.2, width=0.55)
        if not np.isnan(row["ci_lo"]):
            ax.errorbar(i, row["tau"],
                         yerr=[[row["tau"]-row["ci_lo"]],[row["ci_hi"]-row["tau"]]],
                         color="black", capsize=10, lw=1.6)
        ax.text(i, row["tau"]+0.2,
                f"τ={row['tau']:.2f}d\n[{row['ci_lo']:.2f},{row['ci_hi']:.2f}]\n"
                f"n_evt={int(row['n_events'])}",
                ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.axhline(pooled_tau, color="red", ls="--", lw=2,
                label=f"Pooled τ ({title_prefix}) = {pooled_tau:.2f}d")
    ax.set_xticks(x); ax.set_xticklabels(df["season"], fontsize=10, rotation=15)
    valid_taus = df["tau"].dropna().values
    if len(valid_taus) >= 2:
        tau_range = valid_taus.max() - valid_taus.min()
        verdict = ("★ 季節間で τ 安定 → universal 候補" if tau_range < 1.0 else
                   f"⚠ 季節間で τ {tau_range:.1f}d 変動 → universality 棄却")
        vcol = "#1D9E75" if "★" in verdict else "#E85D04"
        ax.text(0.5, 0.97, verdict, transform=ax.transAxes, ha="center", va="top",
                fontsize=12, fontweight="bold", color=vcol,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=vcol, alpha=0.95))
    ax.legend(fontsize=10)
    _ax(ax, f"fig01 — {title_prefix} 季節別 τ\n"
            "季節間で安定なら 'climate universal' 主張可能", "", "τ [days]")
    plt.tight_layout()
    fp = out/f"fig01_seasonal_tau_{title_prefix.lower()}.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_phenology_match(pheno_df, out):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    fig.patch.set_facecolor(FIG_BG)
    df = pheno_df.copy()
    cols = ["#E85D04", "#1D9E75"]
    x = np.arange(len(df))
    for i, (_, row) in enumerate(df.iterrows()):
        if np.isnan(row["tau"]): continue
        ax.bar(i, row["tau"], color=cols[i], alpha=0.85,
                edgecolor="black", lw=1.2, width=0.5)
        ax.errorbar(i, row["tau"],
                     yerr=[[row["tau"]-row["ci_lo"]],[row["ci_hi"]-row["tau"]]],
                     color="black", capsize=12, lw=2)
        ax.text(i, row["tau"]+0.3,
                f"τ={row['tau']:.2f}d\n[{row['ci_lo']:.2f},{row['ci_hi']:.2f}]\n"
                f"n_evt={int(row['n_events'])}",
                ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(df["group"], fontsize=10, rotation=10)
    valid_taus = df["tau"].dropna().values
    valid_cis = df[["ci_lo","ci_hi"]].dropna().values
    if len(valid_taus) == 2 and len(valid_cis) == 2:
        diff = abs(valid_taus[0] - valid_taus[1])
        ci_overlap = not (valid_cis[0,1] < valid_cis[1,0] or valid_cis[1,1] < valid_cis[0,0])
        if diff < 0.5 and ci_overlap:
            verdict = "★ phenology 揃えても τ 一致 → universality 強く支持"; vc = "#1D9E75"
        elif diff > 1.0 and not ci_overlap:
            verdict = "⚠ phenology 揃えると τ 異なる → pooled 値は coincidence"; vc = "#E85D04"
        else:
            verdict = "△ 部分的一致"; vc = "#FFA000"
        ax.text(0.5, 0.97, verdict, transform=ax.transAxes, ha="center", va="top",
                fontsize=12, fontweight="bold", color=vc,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=vc, alpha=0.95))
    _ax(ax, "fig02 — Phenology-matched τ comparison\n"
            "両者とも 'active growing' 季節で比較 (vegetation state を揃える)",
        "", "τ [days]")
    plt.tight_layout()
    fp = out/"fig02_phenology_matched.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_pulse_size_scaling(per_event_df, pulse_res, out):
    good = per_event_df.dropna(subset=["tau","pulse_mm"])
    good = good[good["r2"] >= 0.3]
    if len(good) < 5: return
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.patch.set_facecolor(FIG_BG)

    # pulse_mm vs τ
    ax = axes[0]
    for site, color in [("Oran","#E85D04"), ("Tarazona","#1D9E75")]:
        sub = good[good["site"]==site]
        if len(sub) == 0: continue
        ax.scatter(sub["pulse_mm"], sub["tau"], c=color, s=80, alpha=0.75,
                    edgecolors="black", lw=1, label=f"{site} (n={len(sub)})")
    if pulse_res:
        ax.text(0.97, 0.97,
                 f"Pearson r = {pulse_res['pulse_tau_r']:+.3f}\n"
                 f"p = {pulse_res['pulse_tau_p']:.3f}\n"
                 f"n events = {pulse_res['n']}",
                 transform=ax.transAxes, va="top", ha="right",
                 fontsize=10, family="monospace",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    ax.legend(fontsize=9)
    _ax(ax, "fig03a — τ vs pulse size\n"
            "r ≈ 0 → τ は event size 非依存", "Pulse size [mm]", "τ [days]")

    # amplitude vs τ (identifiability check)
    ax = axes[1]
    good2 = good.dropna(subset=["amplitude"])
    for site, color in [("Oran","#E85D04"), ("Tarazona","#1D9E75")]:
        sub = good2[good2["site"]==site]
        if len(sub) == 0: continue
        ax.scatter(sub["amplitude"], sub["tau"], c=color, s=80, alpha=0.75,
                    edgecolors="black", lw=1, label=f"{site} (n={len(sub)})")
    if pulse_res:
        ax.text(0.97, 0.97,
                 f"Pearson r (amp,τ) = {pulse_res['amplitude_tau_r']:+.3f}\n"
                 f"p = {pulse_res['amplitude_tau_p']:.3f}",
                 transform=ax.transAxes, va="top", ha="right",
                 fontsize=10, family="monospace",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    ax.legend(fontsize=9)
    _ax(ax, "fig03b — τ vs amplitude (identifiability check)\n"
            "r ≈ 0 → τ と amplitude が独立 → 解釈が分離可能",
        "Amplitude (LE_0 − LE_∞) [W/m²]", "τ [days]")
    plt.tight_layout()
    fp = out/"fig03_pulse_size_amplitude.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_final_verdict(verdicts, out):
    fig, ax = plt.subplots(figsize=(12, 7))
    fig.patch.set_facecolor(FIG_BG)
    ax.axis("off")
    fig.suptitle("fig04 — Reviewer counter-attack final verdict\n"
                 "Universal-τ 主張は 3 つの test 全て pass 必要",
                 fontsize=13, fontweight="bold")
    tbl = []
    col_labels = ["Test", "判定", "結論"]
    for v in verdicts:
        tbl.append([v["name"], v["sym"], v["conclusion"]])
    tb = ax.table(cellText=tbl, colLabels=col_labels, loc="center",
                   cellLoc="left", colColours=["#1D9E75"]*3, colWidths=[0.32, 0.1, 0.55])
    tb.auto_set_font_size(False); tb.set_fontsize(11); tb.scale(1, 2.5)
    col_map = {"★": "#C8E6C9", "△": "#FFE0B2", "⚠": "#FFCDD2"}
    for i, v in enumerate(verdicts):
        for j in range(3):
            tb[(i+1, j)].set_facecolor(col_map.get(v["sym"][0], "white"))

    all_pass = all(v["sym"] == "★" for v in verdicts)
    if all_pass:
        msg = "★★★ 3 tests pass → 'effective universal τ ≈ 2.5d' 主張可"
    elif any(v["sym"] == "⚠" for v in verdicts):
        msg = "⚠ いずれか fail → universality 棄却、'effective relaxation timescale' に格下げ"
    else:
        msg = "△ 部分的支持 → 慎重な reframe が必要"
    ax.text(0.5, 0.03, msg, transform=ax.transAxes, ha="center", va="bottom",
             fontsize=13, fontweight="bold",
             color="#1D9E75" if "★★★" in msg else "#E85D04" if "⚠" in msg else "#FFA000",
             bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.95))
    plt.tight_layout()
    fp = out/"fig04_final_verdict.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# MAIN
# ================================================================

def main():
    args = parse_args()
    if args.quick: args.n_boot = 500
    else: args.n_boot = 5000

    parquet, tara_csv, oran_path, out_dir = resolve_paths(args)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("解析A v26 — Reviewer counter-attack")
    print(f"  出力先: {out_dir}")
    print("="*60)

    # ── データ読込 ──
    print(f"\n--- データ読込 ---")
    tara = load_tarazona_data(parquet, tara_csv)
    print(f"  Tarazona: {len(tara)} 日")
    oran = load_oran_daily_master(oran_path)
    print(f"  Oran: {len(oran)} 日, LE 有効={oran['LE'].notna().sum()}")

    # Tarazona format を統一: irrigation events
    tara_for_event = tara.rename(columns={"LE_corr": "LE"}).copy()
    tara_for_event["LE"] = pd.to_numeric(tara_for_event["LE"], errors="coerce")
    tara_for_event["Irrig_mm"] = tara_for_event["Irrig_mm"].fillna(0)

    # Event detection
    print(f"\n--- Event detection ---")
    oran_events = detect_events(oran, "Rain_mm", args.rain_threshold,
                                 min_window=args.min_window, max_window=args.max_window)
    print(f"  Oran rain events:        {len(oran_events)}")
    tara_events = detect_events(tara_for_event, "Irrig_mm", IRRIG_THRESHOLD,
                                 min_window=args.min_window, max_window=args.max_window)
    print(f"  Tarazona irrig events:   {len(tara_events)}")

    # 季節分布
    print(f"\n  Oran events 季節分布:")
    for s, ms in SEASONS.items():
        n = sum(1 for ev in oran_events if ev["event_start"].month in ms)
        print(f"    {s}: {n} events")

    # LE_inf 推定
    le_inf_oran = estimate_le_inf(oran_events) if oran_events else 20
    le_inf_tara = estimate_le_inf(tara_events) if tara_events else 120
    print(f"\n  LE_∞ (Oran)     = {le_inf_oran:.2f} W/m²")
    print(f"  LE_∞ (Tarazona) = {le_inf_tara:.2f} W/m²")

    # ============================================================
    # PART 1: Oran 季節別 τ
    # ============================================================
    print(f"\n{'='*60}\nPART 1: Oran 季節別 τ stratification\n{'='*60}")
    seasonal_oran = stratify_events_by_season(oran_events)
    oran_seasonal_df = compute_seasonal_taus(seasonal_oran, n_boot=args.n_boot)
    print(oran_seasonal_df.to_string(index=False))
    oran_seasonal_df.to_csv(out_dir/"v26_oran_seasonal_tau.csv", index=False)

    valid_taus = oran_seasonal_df["tau"].dropna().values
    if len(valid_taus) >= 2:
        tau_range = float(valid_taus.max() - valid_taus.min())
        s1_sym = "★" if tau_range < 1.0 else ("△" if tau_range < 2.0 else "⚠")
        s1_conc = f"τ range {tau_range:.2f}d across seasons"
    else:
        s1_sym = "—"; s1_conc = "n季節 < 2 で判定不可"
    print(f"\n  Part 1 判定: {s1_sym} {s1_conc}")

    # ============================================================
    # PART 2: Phenology-matched
    # ============================================================
    print(f"\n{'='*60}\nPART 2: Phenology-matched comparison\n{'='*60}")
    pheno_df = phenology_matched(oran_events, tara_events, n_boot=args.n_boot)
    print(pheno_df.to_string(index=False))
    pheno_df.to_csv(out_dir/"v26_phenology_matched_tau.csv", index=False)

    valid = pheno_df.dropna(subset=["tau"])
    if len(valid) == 2:
        taus = valid["tau"].values
        lis = valid["ci_lo"].values
        his = valid["ci_hi"].values
        diff = abs(taus[0] - taus[1])
        ci_overlap = not (his[0] < lis[1] or his[1] < lis[0])
        if diff < 0.5 and ci_overlap:
            s2_sym = "★"; s2_conc = f"両 τ={taus[0]:.2f}/{taus[1]:.2f}d 一致 + CI 重なる"
        elif diff > 1.0 and not ci_overlap:
            s2_sym = "⚠"; s2_conc = f"τ 差 {diff:.2f}d + CI 重ならず → pooled は coincidence"
        else:
            s2_sym = "△"; s2_conc = f"差 {diff:.2f}d, 部分的一致"
    else:
        s2_sym = "—"; s2_conc = "n_group < 2 で判定不可"
    print(f"\n  Part 2 判定: {s2_sym} {s2_conc}")

    # ============================================================
    # PART 3: Per-event + pulse-size scaling
    # ============================================================
    print(f"\n{'='*60}\nPART 3: Per-event + pulse-size scaling\n{'='*60}")
    oran_per_event = per_event_fits(oran_events, "Oran", le_inf_oran)
    tara_per_event = per_event_fits(tara_events, "Tarazona", le_inf_tara)
    per_event = pd.concat([oran_per_event, tara_per_event], ignore_index=True)
    per_event.to_csv(out_dir/"v26_per_event_tau.csv", index=False)

    good = per_event.dropna(subset=["tau"])
    good = good[good["r2"] >= 0.3]
    print(f"  per-event fit 成功 (R²≥0.3): {len(good)}/{len(per_event)}")
    if len(good) > 0:
        print(f"    Oran:     {(good['site']=='Oran').sum()}")
        print(f"    Tarazona: {(good['site']=='Tarazona').sum()}")

    pulse_res = pulse_size_test(per_event)
    if pulse_res:
        print(f"\n  pulse_mm × τ: r = {pulse_res['pulse_tau_r']:+.3f}, p = {pulse_res['pulse_tau_p']:.3f}")
        print(f"  amplitude × τ: r = {pulse_res['amplitude_tau_r']:+.3f}, p = {pulse_res['amplitude_tau_p']:.3f}")
        s3_sym = "★" if abs(pulse_res['pulse_tau_r']) < 0.3 else (
                 "△" if abs(pulse_res['pulse_tau_r']) < 0.5 else "⚠")
        s3_conc = f"pulse-τ |r|={abs(pulse_res['pulse_tau_r']):.2f}"
    else:
        s3_sym = "—"; s3_conc = "サンプル不足"
    print(f"\n  Part 3 判定: {s3_sym} {s3_conc}")

    # ============================================================
    # 統合 verdict
    # ============================================================
    print(f"\n{'='*60}\n★ Final Verdict\n{'='*60}")
    verdicts = [
        dict(name="P1 Oran seasonal τ stratification", sym=s1_sym, conclusion=s1_conc),
        dict(name="P2 Phenology-matched τ", sym=s2_sym, conclusion=s2_conc),
        dict(name="P3 Pulse-size scaling", sym=s3_sym, conclusion=s3_conc),
    ]
    for v in verdicts:
        print(f"  {v['name']:42s} {v['sym']} {v['conclusion']}")

    pd.DataFrame(verdicts).to_csv(out_dir/"v26_final_verdict.csv", index=False)

    all_star = all(v["sym"] == "★" for v in verdicts)
    any_warn = any(v["sym"] == "⚠" for v in verdicts)

    print(f"\n{'='*60}\n★ 論文用 framing 推奨\n{'='*60}")
    if all_star:
        print("""
'Across season, phenology, and pulse-size strata, the fitted exponential
 relaxation timescale τ remained close to 2.5 days at both irrigated
 (Tarazona) and rainfed (Oran) Mediterranean semi-arid sites. The
 invariance of τ across these regimes supports its interpretation as
 an effective ecosystem-atmosphere relaxation timescale characteristic
 of this climate. Irrigation amplifies the ET response magnitude (LE_0 − LE_∞)
 by ~4×, but does not modify the relaxation rate.'""")
    elif any_warn:
        print("""
'Both sites exhibited a fitted exponential relaxation time τ near 2.5 d in
 the pooled analysis. However, stratified comparisons (by season, phenology,
 and event pulse size) revealed that τ varies systematically with these
 factors, indicating that the pooled equality is an aggregate average rather
 than evidence of a single universal mechanism. τ should therefore be
 interpreted as an EFFECTIVE relaxation timescale conditional on the
 ecosystem state at the time of the water input. The amplitude of ET
 response (LE_0 − LE_∞) emerged as the more robust management-distinguishing
 signal, being ~4× larger at the irrigated site.'""")
    else:
        print("""
'Pooled fits gave similar τ ≈ 2.5 d at both sites. Stratified diagnostics
 gave partially consistent results. τ is best reported as an effective
 ecosystem relaxation timescale of order 2–3 d in this Mediterranean
 system, with amplitude (LE_0 − LE_∞) carrying the irrigation-specific
 signal (~4× larger at the drip-irrigated site).'""")

    # ── 可視化 ──
    if not args.no_plots:
        print(f"\n--- 可視化 ---")
        plot_seasonal_taus(oran_seasonal_df,
                            pooled_tau=2.51, out=out_dir, title_prefix="Oran")
        plot_phenology_match(pheno_df, out_dir)
        if pulse_res:
            plot_pulse_size_scaling(per_event, pulse_res, out_dir)
        plot_final_verdict(verdicts, out_dir)

    print(f"\n[done] outputs in {out_dir}/")


if __name__ == "__main__":
    main()
