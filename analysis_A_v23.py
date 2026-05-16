"""
================================================================
解析A v23 : Reviewer-Bulletproof — 5つの S級懸念を全て潰す
================================================================

【目的】
  v22 の τ = 3.33 d [2.62, 5.24] BCa CI を、最強の reviewer
  攻撃から守るための網羅的 robustness 検証。

【S級懸念とその対応】
  S1. Pseudo-replication
      → PART 1: Event-level τ 分布(各灌漑イベントで個別に fit)
      → 真の独立単位は "事件" なのか "日" なのか?

  S2. Independence violation (autocorrelation)
      → PART 2: Block bootstrap (3, 5, 7 日ブロック)
      → iid 仮定を捨てて時系列構造を保持

  S3. LE_∞ 固定の circular reasoning
      → PART 3: Asymptote (d ≥ X) sensitivity scan
      → X = 7, 8, 10, 12, 15 で τ がどれだけ動くか

  S4. Alternative model comparison がない
      → PART 4: exponential / linear / power / logistic を AIC 比較
      → τ という概念は指数関数でしか定義されない

  S5. Seasonality confounding
      → PART 5: Pseudo-event control (Tarazona winter or Oran)
      → 灌漑のない期間で "偽τ" を作って同じ値が出たら artifact

【補助 A級】
  A1. Identifiability: τ × LE_∞ scatter from bootstrap (PART 1 内)
  A2. Aggregation artifact: raw scatter overlay (PART 4 内)

【入力】 v4 出力 parquet + Tarazona daily CSV
【出力】./output_analysis_A_v23/
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

warnings.filterwarnings("ignore")


# ================================================================
# CONFIG
# ================================================================

def parse_args():
    p = argparse.ArgumentParser(description="解析A v23 — Reviewer-bulletproof robustness checks")
    p.add_argument("--parquet", default=None)
    p.add_argument("--csv", default=None)
    p.add_argument("--out", default=None)
    p.add_argument("--n-boot", type=int, default=5000)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


def resolve_paths(args):
    base = Path("/home/shion-nagamine")
    parquet = Path(args.parquet) if args.parquet else \
              base / "bakanposs/analysis_A/daily_classified_v4.parquet"
    csv     = Path(args.csv) if args.csv else \
              base / "Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv"
    out     = Path(args.out) if args.out else Path("./output_analysis_A_v23")
    return parquet, csv, out


IRRIG_THRESHOLD     = 0.5
MIN_IRRIG_PER_MONTH = 2
CI_PCT              = (2.5, 97.5)
FIG_BG              = "#F8F9FA"

# 固定 bounds (v22 と一致)
LE_INF_FLOOR  = 50.0
LE_INF_CEIL   = 200.0
LE_0_FLOOR    = 100.0
LE_0_CEIL     = 350.0
TAU_FLOOR     = 0.5
TAU_CEIL      = 60.0
MIN_PER_BIN   = 5


# ================================================================
# 共通: 読込・前処理(v22 と互換)
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
# モデル定義
# ================================================================

def exp_model_2p(d, le0, tau, le_inf_fixed):
    return le_inf_fixed + (le0 - le_inf_fixed) * np.exp(-d/tau)

def exp_model_3p(d, le_inf, le0, tau):
    return le_inf + (le0 - le_inf) * np.exp(-d/tau)

def linear_model(d, a, b):
    return a + b*d

def power_model(d, a, b, c):
    return a + b * ((d + 1.0) ** c)  # +1 で d=0 を許容

def logistic_model(d, le_inf, le0, d0, k):
    return le_inf + (le0 - le_inf) / (1 + np.exp((d - d0)/k))


# ================================================================
# 共通: bin + fit
# ================================================================

def _bin(arr_d, arr_y, min_per_bin=MIN_PER_BIN):
    df = pd.DataFrame({"d": arr_d, "y": arr_y})
    g = df.groupby("d")["y"].agg(["median","count"]).reset_index()
    return g[g["count"] >= min_per_bin]


def compute_le_inf_observed(arr_d, arr_y, gap_threshold):
    mask = arr_d >= gap_threshold
    n = int(mask.sum())
    if n < 5:
        return np.nan, n
    return float(np.median(arr_y[mask])), n


def fit_2param(arr_d, arr_y, le_inf, min_per_bin=MIN_PER_BIN):
    g = _bin(arr_d, arr_y, min_per_bin)
    if len(g) < 3: return None
    x = g["d"].values.astype(float)
    y = g["median"].values
    def model(d, le0, tau): return exp_model_2p(d, le0, tau, le_inf)
    try:
        popt, _ = curve_fit(model, x, y, p0=[y.max(), 5.0],
                            bounds=([LE_0_FLOOR, TAU_FLOOR], [LE_0_CEIL, TAU_CEIL]),
                            maxfev=8000)
        yp = model(x, *popt)
        sse = float(np.sum((y - yp)**2))
        r2  = float(1 - sse / np.sum((y - y.mean())**2))
        return dict(le0=popt[0], tau=popt[1], le_inf=le_inf,
                    r2=r2, sse=sse, n_bins=len(g),
                    aic=len(g)*np.log(sse/len(g)+1e-12) + 2*2)
    except Exception:
        return None


# ================================================================
# PART 1: Event-level τ (Pseudo-replication 対策)
# ================================================================

def detect_irrigation_events(df_tara, min_window=4, max_window=14):
    """
    各灌漑イベントを抽出し、その後 min_window 日以上 LE データがあるものを返す。
    window_end = 次の灌漑日 − 1 か、event + max_window のうち早い方。
    """
    df_tara = df_tara.sort_values("date").reset_index(drop=True)
    irrig_dates = df_tara[df_tara["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD]["date"].sort_values().values

    events = []
    for i, ed in enumerate(irrig_dates):
        ed = pd.Timestamp(ed)
        # 次の灌漑日
        if i + 1 < len(irrig_dates):
            next_ed = pd.Timestamp(irrig_dates[i+1])
            window_end = min(next_ed - pd.Timedelta(days=1),
                             ed + pd.Timedelta(days=max_window))
        else:
            window_end = ed + pd.Timedelta(days=max_window)

        window_days = (window_end - ed).days + 1
        if window_days < min_window:
            continue
        # LE データを抽出
        evd = df_tara[(df_tara["date"] >= ed) & (df_tara["date"] <= window_end)].copy()
        evd["days_since_event"] = (evd["date"] - ed).dt.days
        evd = evd[evd["LE_corr"].notna()]
        if len(evd) < min_window:
            continue
        events.append(dict(event_start=ed, n_days=window_days, n_le=len(evd), data=evd))
    return events


def fit_event_tau(event, le_inf_fixed):
    """単一イベントで 2-param fit"""
    d = event["data"]["days_since_event"].values.astype(float)
    y = event["data"]["LE_corr"].values
    if len(d) < 4 or np.var(y) < 1.0:
        return None
    def model(dd, le0, tau): return exp_model_2p(dd, le0, tau, le_inf_fixed)
    try:
        popt, _ = curve_fit(model, d, y, p0=[y.max(), 5.0],
                            bounds=([LE_0_FLOOR, TAU_FLOOR], [LE_0_CEIL, TAU_CEIL]),
                            maxfev=5000)
        yp = model(d, *popt)
        sse = float(np.sum((y - yp)**2))
        r2  = float(1 - sse / np.sum((y - y.mean())**2)) if np.var(y) > 0 else np.nan
        return dict(event_start=event["event_start"],
                    n_le=len(d), le0=popt[0], tau=popt[1], r2=r2)
    except Exception:
        return None


# ================================================================
# PART 2: Block bootstrap (Independence violation 対策)
# ================================================================

def block_bootstrap_tau(arr_d, arr_y, le_inf_fixed, block_size, n_boot, seed=42):
    """
    Circular block bootstrap: 連続 block_size 日を1つの単位として復元抽出。
    これで時系列の自己相関が CI に反映される。
    """
    rng = np.random.default_rng(seed)
    n = len(arr_d)
    n_blocks = int(np.ceil(n / block_size))
    taus = []
    n_fail = 0
    for _ in range(n_boot):
        # ブロック開始位置を n_blocks 個ランダム選択
        starts = rng.integers(0, n, size=n_blocks)
        idx = []
        for s in starts:
            idx.extend([(s + j) % n for j in range(block_size)])
        idx = idx[:n]
        b_d, b_y = arr_d[idx], arr_y[idx]
        res = fit_2param(b_d, b_y, le_inf_fixed)
        if res is None:
            n_fail += 1; continue
        taus.append(res["tau"])
    return np.array(taus), n_fail


# ================================================================
# PART 3: Asymptote (LE_∞) sensitivity (Circular reasoning 対策)
# ================================================================

def asymptote_sensitivity(arr_d, arr_y, thresholds):
    """LE_∞ を d ≥ X の中央値とする時、X を変えた時の τ 変動"""
    rows = []
    for thr in thresholds:
        le_inf, n_inf = compute_le_inf_observed(arr_d, arr_y, thr)
        if np.isnan(le_inf):
            rows.append(dict(threshold=thr, le_inf=np.nan, n_inf=n_inf,
                             tau=np.nan, r2=np.nan, success=False))
            continue
        res = fit_2param(arr_d, arr_y, le_inf)
        if res is None:
            rows.append(dict(threshold=thr, le_inf=le_inf, n_inf=n_inf,
                             tau=np.nan, r2=np.nan, success=False))
            continue
        rows.append(dict(threshold=thr, le_inf=le_inf, n_inf=n_inf,
                         tau=res["tau"], r2=res["r2"], success=True))
    return pd.DataFrame(rows)


# ================================================================
# PART 4: Alternative models comparison (τ exists only if exp justified)
# ================================================================

def fit_all_models(arr_d, arr_y, le_inf, min_per_bin=MIN_PER_BIN):
    """4 つのモデルを比較"""
    g = _bin(arr_d, arr_y, min_per_bin)
    if len(g) < 4: return None
    x = g["d"].values.astype(float)
    y = g["median"].values
    n = len(y)
    results = {}

    def _calc(model_name, popt, k):
        yp = globals()[f"_{model_name}_eval"](x, popt, le_inf)
        sse = float(np.sum((y - yp)**2))
        r2  = float(1 - sse / np.sum((y - y.mean())**2))
        aic = n*np.log(sse/n + 1e-12) + 2*k
        return dict(name=model_name, popt=popt, sse=sse, r2=r2, aic=aic, k=k)

    # exponential 2p
    def m_exp(d, le0, tau): return exp_model_2p(d, le0, tau, le_inf)
    try:
        popt, _ = curve_fit(m_exp, x, y, p0=[y.max(), 5.0],
                            bounds=([LE_0_FLOOR, TAU_FLOOR], [LE_0_CEIL, TAU_CEIL]),
                            maxfev=8000)
        yp = m_exp(x, *popt)
        sse = float(np.sum((y-yp)**2))
        r2 = float(1 - sse/np.sum((y-y.mean())**2))
        aic = n*np.log(sse/n+1e-12) + 2*2
        results["exponential"] = dict(name="exponential", popt=popt, sse=sse, r2=r2, aic=aic, k=2)
    except Exception:
        results["exponential"] = None

    # linear
    try:
        popt, _ = curve_fit(linear_model, x, y, p0=[y.mean(), -1.0], maxfev=5000)
        yp = linear_model(x, *popt)
        sse = float(np.sum((y-yp)**2))
        r2 = float(1 - sse/np.sum((y-y.mean())**2))
        aic = n*np.log(sse/n+1e-12) + 2*2
        results["linear"] = dict(name="linear", popt=popt, sse=sse, r2=r2, aic=aic, k=2)
    except Exception:
        results["linear"] = None

    # power (a + b*(d+1)^c)
    try:
        popt, _ = curve_fit(power_model, x, y, p0=[y.min(), y.max()-y.min(), -0.5],
                            bounds=([0, 0, -2], [500, 500, 0]), maxfev=8000)
        yp = power_model(x, *popt)
        sse = float(np.sum((y-yp)**2))
        r2 = float(1 - sse/np.sum((y-y.mean())**2))
        aic = n*np.log(sse/n+1e-12) + 2*3
        results["power"] = dict(name="power", popt=popt, sse=sse, r2=r2, aic=aic, k=3)
    except Exception:
        results["power"] = None

    # logistic
    try:
        popt, _ = curve_fit(logistic_model, x, y,
                            p0=[le_inf, y.max(), x.mean(), 2.0],
                            bounds=([0, 0, 0, 0.1], [500, 500, x.max(), 50]),
                            maxfev=8000)
        yp = logistic_model(x, *popt)
        sse = float(np.sum((y-yp)**2))
        r2 = float(1 - sse/np.sum((y-y.mean())**2))
        aic = n*np.log(sse/n+1e-12) + 2*4
        results["logistic"] = dict(name="logistic", popt=popt, sse=sse, r2=r2, aic=aic, k=4)
    except Exception:
        results["logistic"] = None

    return results, g


# ================================================================
# PART 5: Pseudo-event control (Seasonality confounding 対策)
# ================================================================

def pseudo_event_analysis(df, le_inf_fixed):
    """
    Oran(雨養 → 灌漑なし)で擬似イベントを作り、同じ "見かけ τ" が出るか検証。
    Oran の rain event を擬似灌漑として扱い、降水後の LE 減衰を fit。
    もし τ ≈ 3.3 d が出たら → 単に "イベント後の自然減衰" を見ていた疑い。
    異なる τ なら → Tarazona の τ は灌漑固有の現象。
    """
    oran = df[(df["site"]=="Oran") & df["is_growing"] & df["LE_corr"].notna()].copy()
    if "Rain_mm" not in oran.columns or oran["Rain_mm"].fillna(0).sum() == 0:
        # Rain がなくても、各日を擬似イベントの d=0 として扱う方法
        # ここでは単純に Oran growing 日を疑似 days_since として扱う
        return None
    oran = oran.sort_values("date").reset_index(drop=True)
    # 擬似イベント: Rain > 1mm の日
    rain_thresh = 1.0
    pseudo_starts = oran[oran["Rain_mm"].fillna(0) > rain_thresh].sort_values("date")
    pseudo_dates = pseudo_starts["date"].values

    events = []
    for i, ed in enumerate(pseudo_dates):
        ed = pd.Timestamp(ed)
        if i + 1 < len(pseudo_dates):
            next_ed = pd.Timestamp(pseudo_dates[i+1])
            window_end = min(next_ed - pd.Timedelta(days=1), ed + pd.Timedelta(days=14))
        else:
            window_end = ed + pd.Timedelta(days=14)
        window_days = (window_end - ed).days + 1
        if window_days < 4:
            continue
        evd = oran[(oran["date"] >= ed) & (oran["date"] <= window_end)].copy()
        evd["days_since_event"] = (evd["date"] - ed).dt.days
        evd = evd[evd["LE_corr"].notna()]
        if len(evd) < 4:
            continue
        events.append(evd)

    # 全擬似イベントを pool して fit
    if len(events) == 0:
        return dict(n_events=0, tau=np.nan, r2=np.nan)
    pooled = pd.concat(events)
    arr_d = pooled["days_since_event"].values.astype(float)
    arr_y = pooled["LE_corr"].values
    # Oran 用 LE_∞ (異なるはず — cereal は枯れる)
    oran_le_inf, _ = compute_le_inf_observed(arr_d, arr_y, 7)
    if np.isnan(oran_le_inf):
        oran_le_inf = float(np.percentile(arr_y, 25))
    res = fit_2param(arr_d, arr_y, oran_le_inf)
    return dict(n_events=len(events),
                n_data=len(pooled),
                le_inf=oran_le_inf,
                tau=res["tau"] if res else np.nan,
                r2=res["r2"] if res else np.nan,
                le0=res["le0"] if res else np.nan)


# ================================================================
# 可視化
# ================================================================

def _ax(ax, title="", xl="", yl=""):
    ax.set_title(title, fontweight="bold", fontsize=11, pad=6)
    ax.set_xlabel(xl, fontsize=10); ax.set_ylabel(yl, fontsize=10)
    ax.grid(alpha=0.2, lw=0.7); ax.spines[["top","right"]].set_visible(False)


def plot_event_tau_dist(event_results, pooled_tau, out):
    """fig01: Event-level τ 分布"""
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.patch.set_facecolor(FIG_BG)

    ax = axes[0]
    taus = [r["tau"] for r in event_results if r is not None]
    r2s  = [r["r2"]  for r in event_results if r is not None]
    if len(taus) == 0:
        ax.set_title("No event fits succeeded"); return
    good = [(t, r) for t, r in zip(taus, r2s) if not np.isnan(r) and r >= 0.3]
    good_taus = [t for t, _ in good]

    ax.hist(taus, bins=20, color="#1D9E75", alpha=0.7, edgecolor="black", lw=0.8,
            label=f"All events (n={len(taus)})")
    if good_taus:
        ax.hist(good_taus, bins=20, color="#2E7D32", alpha=0.9, edgecolor="black", lw=0.8,
                label=f"R²≥0.3 events (n={len(good_taus)})")
    ax.axvline(pooled_tau, color="#E85D04", ls="--", lw=2.5,
               label=f"Pooled (v22) τ = {pooled_tau:.2f}d")
    med_evt = np.median(good_taus) if good_taus else np.median(taus)
    ax.axvline(med_evt, color="#7C3AED", ls="-", lw=2,
               label=f"Event median τ = {med_evt:.2f}d")
    _ax(ax, f"fig01a — Event-level τ distribution\n"
            f"擬似反復懸念への直接回答", "τ [days]", "Frequency")
    ax.legend(fontsize=9)

    # サマリー
    ax = axes[1]
    if good_taus:
        stats_dict = dict(
            n_events=len(good_taus),
            mean=np.mean(good_taus),
            median=np.median(good_taus),
            std=np.std(good_taus),
            p25=np.percentile(good_taus, 25),
            p75=np.percentile(good_taus, 75),
        )
        ax.boxplot([good_taus], widths=0.5, patch_artist=True,
                    boxprops=dict(facecolor="#2E7D32", alpha=0.7),
                    medianprops=dict(color="white", lw=2.5))
        ax.axhline(pooled_tau, color="#E85D04", ls="--", lw=2.5,
                    label=f"Pooled (v22) τ = {pooled_tau:.2f}d")
        ax.scatter([1]*len(good_taus), good_taus, color="black", s=20, alpha=0.5, zorder=10)
        ax.set_xticks([1]); ax.set_xticklabels(["Per-event τ"])

        txt = (f"n events  = {stats_dict['n_events']}\n"
               f"mean τ    = {stats_dict['mean']:.2f} d\n"
               f"median τ  = {stats_dict['median']:.2f} d\n"
               f"std        = {stats_dict['std']:.2f} d\n"
               f"IQR       = [{stats_dict['p25']:.2f}, {stats_dict['p75']:.2f}]\n"
               f"vs pooled = {pooled_tau:.2f} d")
        ax.text(1.4, ax.get_ylim()[1]*0.95, txt, family="monospace", fontsize=10,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.9))
        ax.legend(fontsize=9, loc="lower right")
    _ax(ax, "fig01b — 統計サマリー", "", "τ [days]")

    plt.tight_layout()
    fp = out/"fig01_event_tau_distribution.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_block_bootstrap(block_results, iid_taus, pooled_tau, out):
    """fig02: Block bootstrap vs iid bootstrap"""
    fig, ax = plt.subplots(figsize=(11, 6))
    fig.patch.set_facecolor(FIG_BG)
    labels = ["iid\n(v22 baseline)"] + [f"Block\nL={bs}d" for bs in block_results]
    data = [iid_taus] + [block_results[bs]["taus"] for bs in block_results]
    bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False,
                    medianprops=dict(color="white", lw=2.5))
    colors = ["#1A73E8"] + ["#E85D04", "#FFA000", "#7C3AED"][:len(block_results)]
    for box, c in zip(bp["boxes"], colors):
        box.set_facecolor(c); box.set_alpha(0.75)

    ax.axhline(pooled_tau, color="black", ls="--", lw=1.8,
                label=f"Point τ = {pooled_tau:.2f}d")
    ax.set_xticks(range(1, len(labels)+1)); ax.set_xticklabels(labels, fontsize=9)
    _ax(ax, "fig02 — iid bootstrap vs Block bootstrap\n"
            "Block CI が広ければ iid は楽観バイアス", "", "τ [days]")
    ax.legend(fontsize=9)

    # CI 幅を表示
    info_lines = []
    iid_lo, iid_hi = np.percentile(iid_taus, CI_PCT)
    info_lines.append(f"iid : [{iid_lo:.2f}, {iid_hi:.2f}] width={iid_hi-iid_lo:.2f}d")
    for bs in block_results:
        t = block_results[bs]["taus"]
        lo, hi = np.percentile(t, CI_PCT)
        info_lines.append(f"L={bs}: [{lo:.2f}, {hi:.2f}] width={hi-lo:.2f}d")
    ax.text(0.02, 0.97, "\n".join(info_lines), transform=ax.transAxes,
             va="top", ha="left", fontsize=9, family="monospace",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    plt.tight_layout()
    fp = out/"fig02_block_bootstrap.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_asymptote_sensitivity(asym_df, pooled_tau, out):
    """fig03: LE_∞ 閾値感度"""
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    fig.patch.set_facecolor(FIG_BG)

    valid = asym_df[asym_df["success"]]
    ax = axes[0]
    ax.bar(valid["threshold"].astype(str), valid["le_inf"], color="#7C3AED",
            alpha=0.75, edgecolor="black", lw=1)
    for x_str, v, n in zip(valid["threshold"].astype(str), valid["le_inf"], valid["n_inf"]):
        ax.text(x_str, v+1.5, f"{v:.1f}\n(n={int(n)})", ha="center", va="bottom", fontsize=9)
    _ax(ax, "fig03a — LE_∞ vs threshold", "d ≥ X 日", "LE_∞ observed [W/m²]")

    ax = axes[1]
    ax.bar(valid["threshold"].astype(str), valid["tau"], color="#1D9E75",
            alpha=0.85, edgecolor="black", lw=1)
    ax.axhline(pooled_tau, color="#E85D04", ls="--", lw=2,
                label=f"v22 baseline (d≥10): τ={pooled_tau:.2f}d")
    for x_str, t in zip(valid["threshold"].astype(str), valid["tau"]):
        ax.text(x_str, t+0.1, f"{t:.2f}d", ha="center", va="bottom",
                 fontsize=10, fontweight="bold")
    tau_range = valid["tau"].max() - valid["tau"].min()
    verdict_col = "#1D9E75" if tau_range < 0.5 else "#E85D04" if tau_range < 1.0 else "#7C3AED"
    verdict = ("✓ ROBUST" if tau_range < 0.5 else
               "△ moderate" if tau_range < 1.0 else "✗ FRAGILE")
    ax.text(0.97, 0.97,
             f"τ range = {tau_range:.2f}d\n{verdict}",
             transform=ax.transAxes, va="top", ha="right",
             fontsize=11, fontweight="bold", color=verdict_col,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=verdict_col, alpha=0.9))
    ax.legend(fontsize=9)
    _ax(ax, "fig03b — τ vs threshold (LE_∞ definition)", "d ≥ X 日", "τ [days]")

    plt.tight_layout()
    fp = out/"fig03_asymptote_sensitivity.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_model_comparison(model_results, grouped, le_inf, out):
    """fig04: 4 モデル AIC + curve overlay + raw scatter"""
    if model_results is None: return
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.patch.set_facecolor(FIG_BG)

    # 左: curve overlay + raw scatter (集約 artifact 確認)
    ax = axes[0]
    x_obs = grouped["d"].values.astype(float)
    y_obs = grouped["median"].values
    n_per = grouped["count"].values
    x_fit = np.linspace(0, x_obs.max()+1, 200)
    colors = {"exponential":"#1D9E75","linear":"#1A73E8",
              "power":"#FFA000","logistic":"#7C3AED"}
    for name, res in model_results.items():
        if res is None: continue
        if name == "exponential":
            y_fit = exp_model_2p(x_fit, *res["popt"], le_inf)
        elif name == "linear":
            y_fit = linear_model(x_fit, *res["popt"])
        elif name == "power":
            y_fit = power_model(x_fit, *res["popt"])
        elif name == "logistic":
            y_fit = logistic_model(x_fit, *res["popt"])
        ax.plot(x_fit, y_fit, color=colors[name], lw=2.0,
                label=f"{name} (AIC={res['aic']:.1f}, R²={res['r2']:.3f})")
    ax.scatter(x_obs, y_obs, s=40+n_per*3, c="black", edgecolors="white", lw=1,
                zorder=10, label="bin medians")
    ax.legend(fontsize=9, loc="upper right")
    _ax(ax, "fig04a — Alternative models comparison",
        "Days since last irrigation", "LE_corr [W/m²]")

    # 右: AIC ranking
    ax = axes[1]
    valid = {k: v for k, v in model_results.items() if v is not None}
    if valid:
        sorted_models = sorted(valid.items(), key=lambda kv: kv[1]["aic"])
        names = [k for k, _ in sorted_models]
        aics  = [v["aic"] for _, v in sorted_models]
        bars = ax.barh(names, aics, color=[colors[n] for n in names], alpha=0.8,
                        edgecolor="black", lw=1)
        for bar, aic in zip(bars, aics):
            ax.text(aic+0.3, bar.get_y()+bar.get_height()/2, f"{aic:.2f}",
                     va="center", fontsize=10)
        best_aic = aics[0]
        for bar, aic, name in zip(bars, aics, names):
            delta = aic - best_aic
            if delta > 0:
                ax.text(0, bar.get_y()+bar.get_height()/2, f"ΔAIC={delta:.2f}",
                         va="center", ha="left", fontsize=9, color="white", fontweight="bold")
    _ax(ax, "fig04b — AIC ranking\n指数モデルが ΔAIC>2 で勝てば τ 概念は justified",
        "AIC (lower = better)", "")

    plt.tight_layout()
    fp = out/"fig04_model_comparison.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


def plot_pseudo_event(pseudo_res, pooled_tau, out):
    """fig05: Pseudo-event control (Oran で擬似 τ を計算)"""
    if pseudo_res is None or np.isnan(pseudo_res.get("tau", np.nan)):
        return
    fig, ax = plt.subplots(figsize=(11, 6))
    fig.patch.set_facecolor(FIG_BG)
    labels = ["Tarazona\n灌漑後 (真)", "Oran\n降水後 (偽)"]
    taus   = [pooled_tau, pseudo_res["tau"]]
    cols   = ["#1D9E75", "#9E9E9E"]
    bars = ax.bar(labels, taus, color=cols, alpha=0.85, edgecolor="black", lw=1.2, width=0.55)
    for bar, t in zip(bars, taus):
        ax.text(bar.get_x()+bar.get_width()/2, t+0.15, f"τ = {t:.2f}d",
                 ha="center", va="bottom", fontsize=12, fontweight="bold")
    diff = abs(pooled_tau - pseudo_res["tau"])
    verdict = "★ 異なる → 灌漑固有" if diff > 1.0 else "⚠ 似ている → seasonal artifact 疑い"
    vcol = "#1D9E75" if diff > 1.0 else "#E85D04"
    ax.text(0.5, 0.97, verdict, transform=ax.transAxes, ha="center", va="top",
             fontsize=13, fontweight="bold", color=vcol,
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=vcol, alpha=0.95))
    info = (f"Pseudo events: {pseudo_res['n_events']} 個\n"
            f"Pseudo data points: {pseudo_res['n_data']}\n"
            f"Oran LE_∞: {pseudo_res['le_inf']:.1f} W/m²\n"
            f"R²: {pseudo_res['r2']:.3f}")
    ax.text(0.02, 0.97, info, transform=ax.transAxes, va="top", ha="left",
             fontsize=9, family="monospace",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
    _ax(ax, "fig05 — Pseudo-event control\n"
            "Oran 降水後 τ ≠ Tarazona 灌漑後 τ なら seasonality は除外できる",
        "", "τ [days]")
    plt.tight_layout()
    fp = out/"fig05_pseudo_event_control.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# MAIN
# ================================================================

def main():
    args = parse_args()
    if args.quick:
        args.n_boot = 500
        print("[quick mode] n_boot = 500")
    parquet, csv_path, out_dir = resolve_paths(args)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("解析A v23 — Reviewer-bulletproof robustness checks")
    print(f"  出力先: {out_dir}")
    print("="*60)

    df = load_and_merge(parquet, csv_path)
    df = add_days_since_irrig(df)
    df = add_irrig_active_month(df)

    # ベースサブセット (v22 と同一)
    sub = df[(df["site"]=="Tarazona") & df["irrig_active_month"] & df["is_growing"] &
              df["LE_corr"].notna() & df["days_since_irrig"].notna()].copy()
    arr_d = sub["days_since_irrig"].values.astype(float)
    arr_y = sub["LE_corr"].values
    print(f"\n  解析対象: Tarazona 灌漑アクティブ生育日 = {len(sub)} 日")

    le_inf_obs, n_inf = compute_le_inf_observed(arr_d, arr_y, 10)
    print(f"  LE_∞_obs (d≥10) = {le_inf_obs:.2f} W/m² (n={n_inf})")

    # v22 baseline τ
    base_res = fit_2param(arr_d, arr_y, le_inf_obs)
    pooled_tau = base_res["tau"]
    print(f"  v22 baseline (pooled): τ = {pooled_tau:.3f} d, R²={base_res['r2']:.3f}")

    # ============================================================
    # PART 1: Event-level τ (Pseudo-replication 対策)
    # ============================================================
    print(f"\n{'='*60}\nPART 1: Event-level τ (S1: pseudo-replication)\n{'='*60}")
    tara_full = df[(df["site"]=="Tarazona") & df["is_growing"]].copy()
    events = detect_irrigation_events(tara_full, min_window=4, max_window=14)
    print(f"  検出された灌漑イベント: {len(events)} 個")

    event_results = []
    for ev in events:
        res = fit_event_tau(ev, le_inf_obs)
        event_results.append(res)
    successful = [r for r in event_results if r is not None]
    print(f"  Fit 成功: {len(successful)}/{len(events)}")
    if successful:
        ev_df = pd.DataFrame(successful)
        ev_df.to_csv(out_dir/"v23_event_taus.csv", index=False)
        good = ev_df[ev_df["r2"] >= 0.3]
        print(f"  R²≥0.3 のイベント: {len(good)}/{len(successful)}")
        if len(good) > 0:
            print(f"  Event-level τ 統計 (R²≥0.3):")
            print(f"    n      = {len(good)}")
            print(f"    median = {good['tau'].median():.2f} d")
            print(f"    mean   = {good['tau'].mean():.2f} d ± {good['tau'].std():.2f}")
            print(f"    IQR    = [{good['tau'].quantile(0.25):.2f}, {good['tau'].quantile(0.75):.2f}]")
            print(f"    vs pooled = {pooled_tau:.2f} d")
            diff = good["tau"].median() - pooled_tau
            if abs(diff) < 1.0:
                print(f"    ★ Event-level median ≈ pooled (差 {diff:+.2f}d) → "
                      f"pseudo-replication は深刻ではない")
            else:
                print(f"    ⚠ Event-level median ≠ pooled (差 {diff:+.2f}d) → "
                      f"pseudo-replication 影響あり、再検討必要")

    # ============================================================
    # PART 2: Block bootstrap (Independence violation 対策)
    # ============================================================
    print(f"\n{'='*60}\nPART 2: Block bootstrap (S2: independence)\n{'='*60}")
    # iid bootstrap (v22 同等)
    print("  iid bootstrap...")
    rng = np.random.default_rng(42)
    iid_taus = []
    n_fail_iid = 0
    for _ in range(args.n_boot):
        idx = rng.integers(0, len(arr_d), size=len(arr_d))
        r = fit_2param(arr_d[idx], arr_y[idx], le_inf_obs)
        if r is None: n_fail_iid += 1; continue
        iid_taus.append(r["tau"])
    iid_taus = np.array(iid_taus)
    iid_lo, iid_hi = np.percentile(iid_taus, CI_PCT)
    print(f"    iid CI: [{iid_lo:.2f}, {iid_hi:.2f}]  width={iid_hi-iid_lo:.2f}d")

    # block bootstrap
    block_sizes = [3, 5, 7]
    block_results = {}
    for bs in block_sizes:
        print(f"  block bootstrap L={bs}...")
        taus, n_fail = block_bootstrap_tau(arr_d, arr_y, le_inf_obs,
                                            block_size=bs, n_boot=args.n_boot)
        lo, hi = np.percentile(taus, CI_PCT)
        block_results[bs] = dict(taus=taus, n_fail=n_fail, lo=lo, hi=hi, width=hi-lo)
        print(f"    L={bs}: [{lo:.2f}, {hi:.2f}]  width={hi-lo:.2f}d  失敗 {n_fail}")

    # 比較サマリー
    max_block_width = max(block_results[bs]["width"] for bs in block_sizes)
    inflation = max_block_width / (iid_hi - iid_lo)
    print(f"\n  CI width inflation: {inflation:.2f}x (block max / iid)")
    if inflation < 1.3:
        print(f"  ★ Block CI が iid とほぼ同じ → 時系列構造は CI に大きな影響なし")
    elif inflation < 2.0:
        print(f"  △ Block CI が中程度に広い → iid bootstrap は楽観気味")
    else:
        print(f"  ⚠ Block CI が iid の {inflation:.1f}倍 → 時系列構造の影響が大、Block を採用")

    # CSV
    pd.DataFrame([
        dict(method="iid", lo=iid_lo, hi=iid_hi, width=iid_hi-iid_lo, n_fail=n_fail_iid)] +
        [dict(method=f"block_L{bs}", lo=r["lo"], hi=r["hi"], width=r["width"], n_fail=r["n_fail"])
         for bs, r in block_results.items()]).to_csv(out_dir/"v23_block_bootstrap.csv", index=False)

    # ============================================================
    # PART 3: Asymptote sensitivity (Circular reasoning 対策)
    # ============================================================
    print(f"\n{'='*60}\nPART 3: Asymptote sensitivity (S3: circular reasoning)\n{'='*60}")
    asym_df = asymptote_sensitivity(arr_d, arr_y, thresholds=[7, 8, 10, 12, 15])
    print(asym_df.to_string(index=False))
    asym_df.to_csv(out_dir/"v23_asymptote_sensitivity.csv", index=False)
    valid = asym_df[asym_df["success"]]
    tau_range = valid["tau"].max() - valid["tau"].min()
    print(f"\n  τ range across thresholds: {tau_range:.2f}d")
    if tau_range < 0.5:
        print(f"  ★ ROBUST: τ ≈ {valid['tau'].median():.2f} regardless of asymptote definition")
    elif tau_range < 1.0:
        print(f"  △ moderate sensitivity")
    else:
        print(f"  ⚠ FRAGILE: τ depends on threshold → "
              f"reviewer 攻撃あり、報告に幅を持たせるべき")

    # ============================================================
    # PART 4: Alternative models comparison
    # ============================================================
    print(f"\n{'='*60}\nPART 4: Alternative models (S4: τ exists only if exp)\n{'='*60}")
    model_res, grouped = fit_all_models(arr_d, arr_y, le_inf_obs)
    if model_res:
        valid_models = {k: v for k, v in model_res.items() if v is not None}
        sorted_m = sorted(valid_models.items(), key=lambda kv: kv[1]["aic"])
        print(f"  AIC ランキング:")
        best_aic = sorted_m[0][1]["aic"]
        for name, res in sorted_m:
            d = res["aic"] - best_aic
            print(f"    {name:12s}: AIC={res['aic']:7.2f}  ΔAIC={d:+.2f}  R²={res['r2']:.3f}  k={res['k']}")

        # CSV
        pd.DataFrame([
            dict(model=k, aic=v["aic"], r2=v["r2"], sse=v["sse"], k=v["k"])
            for k, v in valid_models.items()
        ]).to_csv(out_dir/"v23_model_comparison.csv", index=False)

        # 判定
        if sorted_m[0][0] == "exponential":
            second = sorted_m[1]
            second_d = second[1]["aic"] - best_aic
            if second_d > 2:
                print(f"  ★ Exponential が AIC で最良 (vs {second[0]} ΔAIC={second_d:+.2f})")
                print(f"    → τ という概念は justified")
            else:
                print(f"  △ Exponential 最良だが {second[0]} と区別不能 (ΔAIC={second_d:+.2f})")
        else:
            print(f"  ⚠ Exponential が最良ではない (best = {sorted_m[0][0]}) → τ の存在意義に疑問")

    # ============================================================
    # PART 5: Pseudo-event control
    # ============================================================
    print(f"\n{'='*60}\nPART 5: Pseudo-event control (S5: seasonality)\n{'='*60}")
    pseudo_res = pseudo_event_analysis(df, le_inf_obs)
    if pseudo_res:
        print(f"  Oran 擬似イベント: {pseudo_res['n_events']} 個 (合計 {pseudo_res['n_data']} 日)")
        print(f"  Oran 擬似 τ = {pseudo_res['tau']:.2f} d  R²={pseudo_res['r2']:.3f}")
        print(f"  Tarazona 真 τ = {pooled_tau:.2f} d")
        diff = abs(pooled_tau - pseudo_res["tau"]) if not np.isnan(pseudo_res["tau"]) else np.nan
        if not np.isnan(diff):
            if diff > 1.0:
                print(f"  ★ τ 異なる (差 {diff:.2f}d) → seasonal artifact ではない、灌漑固有")
            else:
                print(f"  ⚠ τ 似ている (差 {diff:.2f}d) → seasonal artifact の可能性、再検討必要")
        pd.DataFrame([pseudo_res]).to_csv(out_dir/"v23_pseudo_event_control.csv", index=False)

    # ============================================================
    # 可視化
    # ============================================================
    if not args.no_plots:
        print(f"\n--- 可視化 ---")
        plot_event_tau_dist(event_results, pooled_tau, out_dir)
        plot_block_bootstrap(block_results, iid_taus, pooled_tau, out_dir)
        plot_asymptote_sensitivity(asym_df, pooled_tau, out_dir)
        plot_model_comparison(model_res, grouped, le_inf_obs, out_dir)
        plot_pseudo_event(pseudo_res, pooled_tau, out_dir)

    # ============================================================
    # 統合判定
    # ============================================================
    print(f"\n{'='*60}\n★ レビュアー攻撃 5つに対する答え\n{'='*60}")

    # S1
    if len(successful) > 0:
        good = pd.DataFrame(successful)
        good = good[good["r2"] >= 0.3]
        if len(good) > 0:
            evt_med = good["tau"].median()
            diff = abs(evt_med - pooled_tau)
            sym = "★" if diff < 1.0 else "⚠"
            print(f"  S1 Pseudo-replication: {sym} Event-level median τ={evt_med:.2f}d "
                  f"vs pooled {pooled_tau:.2f}d (n={len(good)})")

    # S2
    sym = "★" if inflation < 1.3 else "△" if inflation < 2.0 else "⚠"
    print(f"  S2 Independence    : {sym} Block CI inflation {inflation:.2f}x")

    # S3
    sym = "★" if tau_range < 0.5 else "△" if tau_range < 1.0 else "⚠"
    print(f"  S3 Circular reasoning: {sym} τ range across asymptote thresholds = {tau_range:.2f}d")

    # S4
    if model_res:
        best = sorted_m[0][0]
        second = sorted_m[1] if len(sorted_m) > 1 else None
        if second:
            second_d = second[1]["aic"] - sorted_m[0][1]["aic"]
            sym = "★" if (best == "exponential" and second_d > 2) else "△"
            print(f"  S4 Alternative models: {sym} Best={best}, "
                  f"runner-up {second[0]} ΔAIC={second_d:+.2f}")

    # S5
    if pseudo_res and not np.isnan(pseudo_res.get("tau", np.nan)):
        diff = abs(pooled_tau - pseudo_res["tau"])
        sym = "★" if diff > 1.0 else "⚠"
        print(f"  S5 Seasonality       : {sym} Oran pseudo-τ={pseudo_res['tau']:.2f}d "
              f"vs Tarazona {pooled_tau:.2f}d (差 {diff:.2f}d)")

    print(f"\n[done] outputs in {out_dir}/")


if __name__ == "__main__":
    main()
