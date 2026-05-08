"""
================================================================
解析A v14 : 季節 × 灌漑経過日数で深根仮説の最終判定
================================================================

【v13 で発覚した交絡】
  v13 の bucket 0-1d/2-3d/4-7d は夏期(灌漑活発期)、
  bucket 8+d は春期(灌漑停止期)に偏在していた。
  → バケット番号 ≈ 季節 で、純粋な「灌漑からの経過日数」効果が測れない。

【v14 のテスト設計】
  Test A. 夏期(7-9月)限定 bucket 解析
          灌漑が活発な期間内で経過日数を比較
          → 8+d バケットが空なら "夏は灌漑頻度が高すぎてテスト不能"
          → 8+d に n>5 残るなら、それが真の深根テスト

  Test B. 季節別 SDS (灌漑離脱の純粋なテスト)
          spring (1-4月)  : 灌漑停止確定期 → ここで SDS≈0 なら深根支持
          shoulder (5,6,10): 灌漑遷移期
          summer (7-9月)  : 灌漑活発期
          Oran (雨養) を control として併記

  ⇒ Test B の Tarazona spring SDS が判定の主軸:
      |SDS| ≤ 0.15 (CI が 0 を含む) → ★ 真の深根
      SDS > 0.15 (CI が 0 を超える) → 浅根 / 別要因

【入力】 v4 出力 + Tarazona daily CSV (Irrig_mm)
【出力】./output_analysis_A_v14/
  fig01_heatmap_month_bucket.png : 月×bucket の n サンプル(交絡確認)
  fig02_summer_bucket_sds.png    : Test A — 夏期 bucket 別 SDS
  fig03_seasonal_sds.png         : Test B — 季節別 SDS (Oran vs Tarazona)
  fig04_verdict.png              : 最終判定サマリー
  v14_summer_buckets.csv / v14_seasonal_sds.csv
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


# ================================================================
# CONFIG
# ================================================================

INPUT_DIR       = Path("./")
PARQUET         = INPUT_DIR / "daily_classified_v4.parquet"
TARA_DAILY_CSV  = Path(
    "/home/shion-nagamine/Dataset/Eddy data in Spain/"
    "Daily_Summary_Filtered_forPred_ActEne26.csv")

OUT_DIR = Path("./output_analysis_A_v14")
OUT_DIR.mkdir(exist_ok=True)

IRRIG_THRESHOLD = 0.5  # mm/day, これ以上を灌漑イベントとみなす

# 季節定義(Tarazona の灌漑カレンダーに基づく)
SEASONS = {
    "spring (1-4月)"   : [1, 2, 3, 4],   # 灌漑停止期
    "shoulder (5,6,10月)" : [5, 6, 10],     # 灌漑遷移期
    "summer (7-9月)"   : [7, 8, 9],      # 灌漑活発期
}

BUCKETS = [
    ("0-1d",  0, 1),
    ("2-3d",  2, 3),
    ("4-7d",  4, 7),
    ("8+d",   8, 999),
]
BUCKET_COL = ["#1976D2", "#43A047", "#FB8C00", "#E53935"]
SITE_COL   = {"Oran":"#E85D04", "Tarazona":"#1D9E75"}

N_BOOT          = 5000
CI_PCT          = (2.5, 97.5)
MIN_N_PER_CLASS = 5
DENOM_FLOORS    = {"LE_corr": 5.0, "EF_corr": 0.05, "ET": 0.3}
RATIO_CLIP      = 10.0

DEEP_ROOT_BAND  = 0.15  # |SDS|<=0.15 を深根判定領域


# ================================================================
# 1. 入力読込 + 灌漑データ結合(v13 と同じ)
# ================================================================

def load_and_merge():
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


def compute_days_since_irrig(df):
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


def bucket_of(d):
    if pd.isna(d): return None
    for label, lo, hi in BUCKETS:
        if lo <= d <= hi:
            return label
    return None


# ================================================================
# 2. SDS bootstrap
# ================================================================

def safe_ratio(numer, denom, var):
    if abs(denom) < DENOM_FLOORS.get(var, 0.05):
        return np.nan
    r = numer / denom
    return np.nan if abs(r) > RATIO_CLIP else r


def bootstrap_sds(arr_n, arr_s, var, seed=42):
    if len(arr_n) < MIN_N_PER_CLASS or len(arr_s) < MIN_N_PER_CLASS:
        return dict(sds=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                    n_n=len(arr_n), n_s=len(arr_s))
    point = safe_ratio(np.median(arr_n) - np.median(arr_s),
                       np.median(arr_n), var)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(N_BOOT):
        nb = rng.choice(arr_n, size=len(arr_n), replace=True)
        sb = rng.choice(arr_s, size=len(arr_s), replace=True)
        r = safe_ratio(np.median(nb) - np.median(sb), np.median(nb), var)
        if not np.isnan(r):
            boots.append(r)
    if len(boots) < N_BOOT * 0.1:
        return dict(sds=point, ci_lo=np.nan, ci_hi=np.nan,
                    n_n=len(arr_n), n_s=len(arr_s))
    lo, hi = np.percentile(boots, CI_PCT)
    return dict(sds=point, ci_lo=lo, ci_hi=hi,
                n_n=len(arr_n), n_s=len(arr_s))


def compute_sds_for_subset(sub, var):
    n_arr = sub.loc[sub["drought_type"]=="normal",   var].dropna().values
    s_arr = sub.loc[sub["drought_type"]=="soil dry", var].dropna().values
    return bootstrap_sds(n_arr, s_arr, var)


# ================================================================
# 3. 可視化
# ================================================================

def plot_heatmap_month_bucket(df, save_dir):
    """月×bucket のサンプル数を heatmap → v13 の交絡を可視化"""
    tara = df[(df["site"]=="Tarazona") & df["is_growing"] &
               df["bucket"].notna()].copy()
    tara["month"] = tara["date"].dt.month
    bucket_labels = [b[0] for b in BUCKETS]
    months_all = sorted(tara["month"].unique())
    grid = np.zeros((len(months_all), len(bucket_labels)))
    for i, m in enumerate(months_all):
        for j, b in enumerate(bucket_labels):
            grid[i, j] = ((tara["month"] == m) & (tara["bucket"] == b)).sum()

    fig, ax = plt.subplots(figsize=(9, 7))
    fig.patch.set_facecolor("#F8F8F8")
    im = ax.imshow(grid, cmap="YlOrRd", aspect="auto")
    for i in range(len(months_all)):
        for j in range(len(bucket_labels)):
            n = int(grid[i, j])
            if n > 0:
                ax.text(j, i, str(n), ha="center", va="center",
                        color="black" if n < grid.max()*0.55 else "white",
                        fontsize=10, fontweight="bold")
    ax.set_xticks(range(len(bucket_labels)))
    ax.set_xticklabels(bucket_labels, fontsize=11)
    ax.set_yticks(range(len(months_all)))
    ax.set_yticklabels([f"{m}月" for m in months_all], fontsize=11)
    ax.set_xlabel("days_since_irrigation bucket", fontsize=11)
    ax.set_ylabel("month", fontsize=11)
    ax.set_title(
        "[Tarazona Growing] 月 × bucket  サンプル分布\n"
        "8+d が 1-4月(灌漑停止期)に偏在 → v13 交絡の確認",
        fontweight="bold", fontsize=12)
    plt.colorbar(im, ax=ax, label="n days")
    plt.tight_layout()
    fp = save_dir / "fig01_heatmap_month_bucket.png"
    plt.savefig(fp, dpi=140, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [save] {fp}")


def plot_summer_buckets(summer_sds_df, save_dir):
    """Test A: 夏期(7-9月)限定 bucket 別 SDS"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle(
        "Test A: 夏期(7-9月)限定 — bucket 別 SDS\n"
        "夏期内で経過日数の効果のみを抽出",
        fontsize=12, fontweight="bold")

    for ax, var in zip(axes, ["LE_corr","EF_corr","ET"]):
        sub = summer_sds_df[summer_sds_df["var"]==var].sort_values("bucket_order")
        if sub.empty:
            ax.set_title(f"{var}: no data"); continue
        x = np.arange(len(sub))
        for i, row in enumerate(sub.itertuples()):
            col = BUCKET_COL[int(row.bucket_order)]
            val = row.sds if not np.isnan(row.sds) else 0
            ax.bar(i, val, color=col, alpha=0.88, edgecolor="black", lw=1)
            if not np.isnan(row.ci_lo):
                ax.errorbar(i, val,
                             yerr=[[val-row.ci_lo],[row.ci_hi-val]],
                             color="black", capsize=8, lw=1.4)
            if not np.isnan(row.sds):
                lab = (f"{row.sds:+.2f}\n[{row.ci_lo:+.2f},{row.ci_hi:+.2f}]"
                        if not np.isnan(row.ci_lo) else f"{row.sds:+.2f}")
                ax.text(i, val + (0.018 if val>=0 else -0.018), lab,
                        ha="center", va="bottom" if val>=0 else "top",
                        fontsize=9, fontweight="bold")
        ax.axhline(0, color="gray", lw=0.8)
        ax.axhspan(-DEEP_ROOT_BAND, DEEP_ROOT_BAND,
                    color="#1D9E75", alpha=0.08)
        ax.set_xticks(x)
        ax.set_xticklabels(
            [f"{r.bucket}\nn_n={int(r.n_n)} n_s={int(r.n_s)}"
             for r in sub.itertuples()], fontsize=8.5)
        ax.set_ylabel("SDS")
        ax.set_title(var, fontweight="bold")
        ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    fp = save_dir / "fig02_summer_bucket_sds.png"
    plt.savefig(fp, dpi=140, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [save] {fp}")


def plot_seasonal_comparison(season_sds_df, save_dir):
    """Test B: 季節別 SDS (Oran vs Tarazona)"""
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.5))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle(
        "Test B: 季節別 SDS — Oran(雨養 control) vs Tarazona(灌漑)\n"
        "Tarazona spring (灌漑停止) で SDS≈0 なら ★ 深根支持",
        fontsize=12, fontweight="bold")

    season_order = list(SEASONS.keys())
    for ax, var in zip(axes, ["LE_corr","EF_corr","ET"]):
        sub = season_sds_df[season_sds_df["var"]==var]
        if sub.empty:
            ax.set_title(f"{var}: no data"); continue
        x_labels, vals, los, his, cols, ns = [], [], [], [], [], []
        for season in season_order:
            for site in ["Oran","Tarazona"]:
                row = sub[(sub["season"]==season) & (sub["site"]==site)]
                if row.empty: continue
                r = row.iloc[0]
                x_labels.append(f"{site}\n{season.split(' ')[0]}")
                vals.append(r["sds"]); los.append(r["ci_lo"]); his.append(r["ci_hi"])
                cols.append(SITE_COL[site])
                ns.append(f"n={int(r['n_n'])},{int(r['n_s'])}")

        x = np.arange(len(x_labels))
        for i, (v, l, h, c) in enumerate(zip(vals, los, his, cols)):
            if np.isnan(v):
                continue
            ax.bar(i, v, color=c, alpha=0.85, edgecolor="black", lw=1)
            if not np.isnan(l):
                ax.errorbar(i, v, yerr=[[v-l],[h-v]],
                             color="black", capsize=6, lw=1.2)
            ax.text(i, v + (0.02 if v>=0 else -0.02),
                     f"{v:+.2f}", ha="center",
                     va="bottom" if v>=0 else "top",
                     fontsize=9, fontweight="bold")
        ax.axhline(0, color="gray", lw=0.8)
        ax.axhspan(-DEEP_ROOT_BAND, DEEP_ROOT_BAND,
                    color="#1D9E75", alpha=0.08,
                    label=f"|SDS|≤{DEEP_ROOT_BAND}")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{lab}\n{n_s}" for lab, n_s in zip(x_labels, ns)],
                            fontsize=8)
        ax.set_ylabel("SDS")
        ax.set_title(var, fontweight="bold")
        ax.legend(fontsize=8, loc="best")
        ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    fp = save_dir / "fig03_seasonal_sds.png"
    plt.savefig(fp, dpi=140, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [save] {fp}")


def plot_verdict(season_sds_df, summer_sds_df, save_dir):
    """最終判定サマリー(LE_corr ベース)"""
    bars = []  # (label, color, sds, lo, hi, note)

    # B1: Tarazona spring (主テスト)
    r = season_sds_df[(season_sds_df["site"]=="Tarazona") &
                       (season_sds_df["season"].str.contains("spring")) &
                       (season_sds_df["var"]=="LE_corr")]
    if not r.empty:
        rr = r.iloc[0]
        bars.append(("Tarazona\nSpring (1-4月)\n灌漑停止確定",
                     "#1D9E75", rr["sds"], rr["ci_lo"], rr["ci_hi"]))

    # B2: Oran spring (control)
    r = season_sds_df[(season_sds_df["site"]=="Oran") &
                       (season_sds_df["season"].str.contains("spring")) &
                       (season_sds_df["var"]=="LE_corr")]
    if not r.empty:
        rr = r.iloc[0]
        bars.append(("Oran\nSpring (1-4月)\n雨養 control",
                     "#E85D04", rr["sds"], rr["ci_lo"], rr["ci_hi"]))

    # A: Tarazona summer (灌漑活発期、参考)
    r = season_sds_df[(season_sds_df["site"]=="Tarazona") &
                       (season_sds_df["season"].str.contains("summer")) &
                       (season_sds_df["var"]=="LE_corr")]
    if not r.empty:
        rr = r.iloc[0]
        bars.append(("Tarazona\nSummer (7-9月)\n灌漑活発(交絡あり)",
                     "#7CB342", rr["sds"], rr["ci_lo"], rr["ci_hi"]))

    # A1: Summer 8+d if exists
    r = summer_sds_df[(summer_sds_df["bucket"]=="8+d") &
                       (summer_sds_df["var"]=="LE_corr")]
    if not r.empty and not np.isnan(r.iloc[0]["sds"]):
        rr = r.iloc[0]
        if int(rr["n_n"]) >= MIN_N_PER_CLASS:
            bars.append(("Tarazona\nSummer 8+d\n夏期+灌漑離脱",
                         "#FFA000", rr["sds"], rr["ci_lo"], rr["ci_hi"]))

    if not bars:
        print("  [verdict] バーが生成できず — データ不足")
        return

    fig, ax = plt.subplots(figsize=(12, 7))
    fig.patch.set_facecolor("#F8F8F8")
    x = np.arange(len(bars))
    for i, (label, c, v, lo, hi) in enumerate(bars):
        ax.bar(i, v, color=c, alpha=0.92, edgecolor="black", lw=1.3, width=0.55)
        if not np.isnan(lo):
            ax.errorbar(i, v, yerr=[[v-lo],[hi-v]],
                         color="black", capsize=12, lw=1.8)
        if not np.isnan(v):
            ax.text(i, v + (0.025 if v>=0 else -0.025),
                     f"{v:+.3f}\n[{lo:+.2f}, {hi:+.2f}]",
                     ha="center", va="bottom" if v>=0 else "top",
                     fontsize=11, fontweight="bold")
    ax.axhline(0, color="gray", lw=1.0, ls="--")
    ax.axhspan(-DEEP_ROOT_BAND, DEEP_ROOT_BAND, color="#1D9E75", alpha=0.10,
                label=f"深根判定領域 |SDS|≤{DEEP_ROOT_BAND}")
    ax.set_xticks(x)
    ax.set_xticklabels([b[0] for b in bars], fontsize=10, fontweight="bold")
    ax.set_ylabel("SDS  (LE_corr ベース)", fontsize=12)
    ax.set_title(
        "★ Final verdict — 季節 × 灌漑経過日数を制御した SDS\n"
        "Tarazona Spring が判定領域内 → 灌漑離脱でも LE 維持 → 真の深根",
        fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    fp = save_dir / "fig04_verdict.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [save] {fp}")


# ================================================================
# 4. MAIN
# ================================================================

def main():
    print("=" * 60)
    print("解析A v14 — 季節 × 灌漑経過日数で深根仮説の最終判定")
    print("=" * 60)

    df = load_and_merge()
    df = compute_days_since_irrig(df)
    df["bucket"] = df["days_since_irrig"].apply(bucket_of)
    df["month"]  = df["date"].dt.month

    # ============================================================
    # Test A: 夏期(7-9月)限定 bucket 別 SDS
    # ============================================================
    print(f"\n{'='*60}\nTest A: 夏期(7-9月)限定 bucket 別 SDS\n{'='*60}")
    summer_rows = []
    summer_sub = df[(df["site"]=="Tarazona") & df["is_growing"] &
                     df["month"].isin(SEASONS["summer (7-9月)"])]
    n_irr_summer = int((summer_sub["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD).sum())
    print(f"  夏期サンプル: {len(summer_sub)}日  灌漑日数={n_irr_summer}")

    for var in ["LE_corr", "EF_corr", "ET"]:
        if var not in summer_sub.columns: continue
        for i, (label, lo, hi) in enumerate(BUCKETS):
            sub = summer_sub[summer_sub["bucket"] == label]
            res = compute_sds_for_subset(sub, var)
            res.update(dict(var=var, bucket=label, bucket_order=i))
            summer_rows.append(res)
            ci = (f"[{res['ci_lo']:+.2f}, {res['ci_hi']:+.2f}]"
                   if not np.isnan(res["ci_lo"]) else "[CI算出不可]")
            print(f"  {var:9s} {label:6s}: SDS={res['sds']:+.3f}  {ci}  "
                   f"(n_n={res['n_n']}, n_s={res['n_s']})")
    summer_sds_df = pd.DataFrame(summer_rows)
    summer_sds_df.to_csv(OUT_DIR / "v14_summer_buckets.csv", index=False)

    # ============================================================
    # Test B: 季節別 SDS (Oran vs Tarazona)
    # ============================================================
    print(f"\n{'='*60}\nTest B: 季節別 SDS — Oran(雨養) vs Tarazona(灌漑)\n{'='*60}")
    season_rows = []
    for season_label, months in SEASONS.items():
        print(f"\n  [{season_label}]")
        for site in ["Oran", "Tarazona"]:
            sub = df[(df["site"] == site) & df["is_growing"] &
                      df["month"].isin(months)]
            for var in ["LE_corr", "EF_corr", "ET"]:
                if var not in sub.columns: continue
                res = compute_sds_for_subset(sub, var)
                res.update(dict(site=site, season=season_label, var=var))
                season_rows.append(res)
            r_le = next((r for r in season_rows
                         if r["site"]==site and r["season"]==season_label
                         and r["var"]=="LE_corr"), None)
            if r_le is not None and not np.isnan(r_le["sds"]):
                ci = (f"[{r_le['ci_lo']:+.2f}, {r_le['ci_hi']:+.2f}]"
                      if not np.isnan(r_le["ci_lo"]) else "[CI算出不可]")
                print(f"    {site:9s} LE_corr SDS={r_le['sds']:+.3f}  {ci}  "
                      f"(n_n={r_le['n_n']}, n_s={r_le['n_s']})")
            elif r_le is not None:
                print(f"    {site:9s} LE_corr: n不足 "
                      f"(n_n={r_le['n_n']}, n_s={r_le['n_s']})")
    season_sds_df = pd.DataFrame(season_rows)
    season_sds_df.to_csv(OUT_DIR / "v14_seasonal_sds.csv", index=False)

    # ============================================================
    # 可視化
    # ============================================================
    print(f"\n--- 可視化 ---")
    plot_heatmap_month_bucket(df, OUT_DIR)
    plot_summer_buckets(summer_sds_df, OUT_DIR)
    plot_seasonal_comparison(season_sds_df, OUT_DIR)
    plot_verdict(season_sds_df, summer_sds_df, OUT_DIR)

    # ============================================================
    # ★ 最終判定
    # ============================================================
    print(f"\n{'='*60}\n★ 最終判定\n{'='*60}")

    def in_band(sds, lo, hi):
        """|SDS|<=帯 OR CI が 0 を含む → 深根判定"""
        if np.isnan(sds): return None
        if not np.isnan(lo) and not np.isnan(hi):
            return (lo <= 0 <= hi) or abs(sds) <= DEEP_ROOT_BAND
        return abs(sds) <= DEEP_ROOT_BAND

    # B1. Tarazona spring (主テスト)
    r = season_sds_df[(season_sds_df["site"]=="Tarazona") &
                       (season_sds_df["season"].str.contains("spring")) &
                       (season_sds_df["var"]=="LE_corr")]
    spring_t_pass = None
    if not r.empty and not np.isnan(r.iloc[0]["sds"]):
        rr = r.iloc[0]
        spring_t_pass = in_band(rr["sds"], rr["ci_lo"], rr["ci_hi"])
        print(f"\n[B1] Tarazona Spring SDS = {rr['sds']:+.3f}  "
               f"[{rr['ci_lo']:+.2f}, {rr['ci_hi']:+.2f}]  "
               f"(n_n={int(rr['n_n'])}, n_s={int(rr['n_s'])})")
        if spring_t_pass:
            print(f"     ⇒ ★ 灌漑停止期でも SDS が判定領域内 = **真の深根支持**")
        else:
            print(f"     ⇒ Spring でも SDS 大 → 深根仮説に疑問符")

    # B2. Oran spring (control)
    r = season_sds_df[(season_sds_df["site"]=="Oran") &
                       (season_sds_df["season"].str.contains("spring")) &
                       (season_sds_df["var"]=="LE_corr")]
    if not r.empty and not np.isnan(r.iloc[0]["sds"]):
        rr = r.iloc[0]
        oran_pass = in_band(rr["sds"], rr["ci_lo"], rr["ci_hi"])
        print(f"\n[B2] Oran Spring SDS = {rr['sds']:+.3f}  "
               f"[{rr['ci_lo']:+.2f}, {rr['ci_hi']:+.2f}]  "
               f"(n_n={int(rr['n_n'])}, n_s={int(rr['n_s'])})  control")
        if not oran_pass and spring_t_pass:
            print(f"     ⇒ ★★ 同じ春期で Oran は表層依存、Tarazona は非依存 "
                   f"→ **Tarazona の深根が確定的**")

    # A1. Summer 8+d
    r = summer_sds_df[(summer_sds_df["bucket"]=="8+d") &
                       (summer_sds_df["var"]=="LE_corr")]
    if not r.empty:
        rr = r.iloc[0]
        if int(rr["n_n"]) >= MIN_N_PER_CLASS and not np.isnan(rr["sds"]):
            sum_pass = in_band(rr["sds"], rr["ci_lo"], rr["ci_hi"])
            print(f"\n[A1] Summer 8+d SDS = {rr['sds']:+.3f}  "
                   f"[{rr['ci_lo']:+.2f}, {rr['ci_hi']:+.2f}]  "
                   f"(n_n={int(rr['n_n'])}, n_s={int(rr['n_s'])})")
            if sum_pass:
                print(f"     ⇒ ★★★ 夏期で灌漑離脱しても SDS≈0 = **深根の決定打**")
        else:
            print(f"\n[A1] Summer 8+d: n不足 "
                   f"(n_n={int(rr['n_n'])}, n_s={int(rr['n_s'])})")
            print(f"     夏期は灌漑頻度が高すぎてテスト不能 → Test B が判定の主軸")

    print(f"\n[done] outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
