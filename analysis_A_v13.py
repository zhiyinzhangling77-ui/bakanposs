"""
================================================================
解析A v13 : 灌漑データによる深根仮説の真偽判定
================================================================

【背景】
  v12 までで Tarazona は表層SWC低下時に蒸散を維持することを示したが、
  Tarazona は灌漑あり。表層センサー(5cm)より深い領域に水が届く点滴
  灌漑は、見かけ上「土壌SWC非感受 = 深根」と区別不能。

【コア論理】
  最後の灌漑からの経過日数(days_since_irrig)で日次データを階層化し、
  各階層内で SDS を再計算する:
    - SDS が経過日数で平坦 → 灌漑から離れても LE 維持 = 真の深根
    - SDS が経過日数で上昇 → 灌漑直後は維持、時間とともに低下 = 灌漑依存

【入力】
  daily_classified_v4.parquet           (v4 出力)
  Tarazona daily CSV (Irrig_mm 含む)    (← ここに灌漑記録)

【出力】./output_analysis_A_v13/
  v13_irrig_classified.csv      : days_since_irrig タグ付き全データ
  v13_sds_by_days_since.csv     : バケット別 SDS + 95%CI
  fig01_timeseries.png          : 時系列(Irrig+Rain / SWC / LE / days_since)
  fig02_le_vs_days_since.png    : LE vs 経過日数 scatter + 中央値線
  fig03_sds_by_days_since.png   : ★主結果 SDS vs 経過日数バケット
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

warnings.filterwarnings("ignore")


# ================================================================
# CONFIG
# ================================================================

INPUT_DIR       = Path("./")
PARQUET         = INPUT_DIR / "daily_classified_v4.parquet"
TARA_DAILY_CSV  = Path(
    "/home/shion-nagamine/Dataset/Eddy data in Spain/"
    "Daily_Summary_Filtered_forPred_ActEne26.csv")

OUT_DIR = Path("./output_analysis_A_v13")
OUT_DIR.mkdir(exist_ok=True)

# 灌漑判定: Irrig_mm > 0.5 を灌漑イベントとみなす(計測ノイズ除外)
IRRIG_THRESHOLD = 0.5  # mm/day

# 経過日数のバケット定義 (灌漑からどれだけ離れたか)
BUCKETS = [
    ("0-1d (直後)",   0, 1),
    ("2-3d (短期)",   2, 3),
    ("4-7d (中期)",   4, 7),
    ("8+d (長期)",    8, 999),
]
BUCKET_COL = ["#1976D2", "#43A047", "#FB8C00", "#E53935"]

# Bootstrap
N_BOOT          = 5000
CI_PCT          = (2.5, 97.5)
MIN_N_PER_CLASS = 5
DENOM_FLOORS    = {"LE_corr": 5.0, "EF_corr": 0.05, "ET": 0.3}
RATIO_CLIP      = 10.0

DROUGHT_COL = {"normal":"#2196F3","soil dry":"#FF9800",
               "atm dry":"#4CAF50","compound":"#F44336"}


# ================================================================
# 1. 入力読込 + 灌漑データ結合
# ================================================================

def load_and_merge():
    daily = pd.read_parquet(PARQUET)
    daily["date"] = pd.to_datetime(daily["date"])

    if not TARA_DAILY_CSV.exists():
        raise FileNotFoundError(f"Tarazona daily CSV not found: {TARA_DAILY_CSV}")

    csv = pd.read_csv(TARA_DAILY_CSV)
    csv["date"] = pd.to_datetime(csv["date"], errors="coerce")
    irrig_cols = [c for c in ["Irrig_mm","Rain_mm","IrrigRain_mm"] if c in csv.columns]
    if "Irrig_mm" not in csv.columns:
        raise ValueError("'Irrig_mm' 列が CSV にありません")
    print(f"[load] Tarazona CSV 灌漑関連列: {irrig_cols}")

    extra = (csv[["date"] + irrig_cols]
             .dropna(subset=["date"])
             .drop_duplicates("date")
             .sort_values("date"))

    # Tarazona は灌漑データを結合、Oran は雨養なので 0 を入れる
    tara = (daily[daily["site"]=="Tarazona"]
            .merge(extra, on="date", how="left"))
    oran = daily[daily["site"]=="Oran"].copy()
    for c in irrig_cols:
        oran[c] = 0.0

    df = pd.concat([oran, tara]).sort_values(["site","date"]).reset_index(drop=True)
    n_ev = (df.loc[df["site"]=="Tarazona","Irrig_mm"].fillna(0) > IRRIG_THRESHOLD).sum()
    print(f"[load] Tarazona 灌漑イベント検出: {n_ev} 日 "
          f"(threshold={IRRIG_THRESHOLD}mm)")
    return df


# ================================================================
# 2. days_since_irrig 計算 + バケット割当
# ================================================================

def compute_days_since_irrig(df):
    """サイト別に最後の灌漑からの経過日数を付与(灌漑日=0)"""
    df = df.sort_values(["site","date"]).reset_index(drop=True).copy()
    df["days_since_irrig"] = np.nan
    for site, idx in df.groupby("site").indices.items():
        sub = df.loc[idx].sort_values("date")
        last = None
        out = []
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
# 3. SDS bootstrap
# ================================================================

def safe_ratio(numer, denom, var):
    if abs(denom) < DENOM_FLOORS.get(var, 0.05):
        return np.nan
    r = numer / denom
    return np.nan if abs(r) > RATIO_CLIP else r


def bootstrap_sds(arr_n, arr_s, var, n_boot=N_BOOT, seed=42):
    """SDS = (median(N) - median(S)) / median(N)、bootstrap 95%CI"""
    if len(arr_n) < MIN_N_PER_CLASS or len(arr_s) < MIN_N_PER_CLASS:
        return dict(sds=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                    n_n=len(arr_n), n_s=len(arr_s))
    point = safe_ratio(np.median(arr_n) - np.median(arr_s),
                       np.median(arr_n), var)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        nb = rng.choice(arr_n, size=len(arr_n), replace=True)
        sb = rng.choice(arr_s, size=len(arr_s), replace=True)
        r = safe_ratio(np.median(nb) - np.median(sb), np.median(nb), var)
        if not np.isnan(r):
            boots.append(r)
    if len(boots) < n_boot * 0.1:
        return dict(sds=point, ci_lo=np.nan, ci_hi=np.nan,
                    n_n=len(arr_n), n_s=len(arr_s))
    lo, hi = np.percentile(boots, CI_PCT)
    return dict(sds=point, ci_lo=lo, ci_hi=hi,
                n_n=len(arr_n), n_s=len(arr_s))


# ================================================================
# 4. 可視化
# ================================================================

def plot_timeseries(df, save_dir):
    """Tarazona の時系列: 水入力 / SWC / LE / days_since_irrig"""
    tara = df[df["site"]=="Tarazona"].sort_values("date").copy()
    fig, axes = plt.subplots(4, 1, figsize=(16, 12), sharex=True)
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle("[Tarazona] 灌漑+降水 / SWC / LE / 経過日数",
                 fontsize=13, fontweight="bold")

    # (1) 水入力
    if "Rain_mm" in tara.columns:
        axes[0].bar(tara["date"], tara["Rain_mm"].fillna(0),
                    color="#3B8BD4", alpha=0.7, label="Rain", width=1.0)
    axes[0].bar(tara["date"], tara["Irrig_mm"].fillna(0),
                color="#FF6F00", alpha=0.85, label="Irrigation", width=1.0)
    axes[0].set_ylabel("mm/day"); axes[0].legend(loc="upper right")
    axes[0].grid(alpha=0.25); axes[0].set_title("Water inputs")

    # (2) SWC
    axes[1].plot(tara["date"], tara["SWC_mean"], color="#E85D04", lw=0.9)
    axes[1].set_ylabel("SWC (%)"); axes[1].grid(alpha=0.25)
    axes[1].set_title("Surface SWC (5cm)")

    # (3) LE
    axes[2].plot(tara["date"], tara["LE_corr"], color="#1D9E75", lw=0.9)
    axes[2].set_ylabel("LE_corr (W/m²)"); axes[2].grid(alpha=0.25)
    axes[2].set_title("LE (closure-corrected)")

    # (4) days_since_irrig
    axes[3].plot(tara["date"], tara["days_since_irrig"],
                  color="#7C3AED", lw=0.9, drawstyle="steps-post")
    for _, _, lo in [(b, c, c) for b, c in zip(BUCKETS[1:], BUCKET_COL[1:])]:
        pass  # placeholder
    axes[3].axhline(7, color="gray", lw=0.8, ls="--", alpha=0.6,
                     label="7d (典型的灌漑間隔)")
    axes[3].set_ylabel("Days since irrigation")
    axes[3].set_xlabel("Date")
    axes[3].grid(alpha=0.25); axes[3].legend(loc="upper right")
    axes[3].set_title("Days since last irrigation event")

    for ax in axes:
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    plt.setp(axes[-1].xaxis.get_majorticklabels(), rotation=30, ha="right")

    plt.tight_layout()
    fp = save_dir / "fig01_timeseries.png"
    plt.savefig(fp, dpi=130, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [save] {fp}")


def plot_le_vs_days(df, save_dir):
    """LE_corr vs days_since_irrig (Tarazona growing) — 中央値の経時変化を見る"""
    tara = df[(df["site"]=="Tarazona") & df["is_growing"] &
              df["days_since_irrig"].notna() & df["LE_corr"].notna()]
    fig, ax = plt.subplots(figsize=(11, 6.5))
    fig.patch.set_facecolor("#F8F8F8")

    # drought_type 別に色分けして散布
    for dt, col in DROUGHT_COL.items():
        sub = tara[tara["drought_type"]==dt]
        ax.scatter(sub["days_since_irrig"], sub["LE_corr"],
                    c=col, alpha=0.55, s=24, label=f"{dt} (n={len(sub)})",
                    edgecolors="none")

    # 1日ごとの中央値線(>=3 サンプルあるビンのみ)
    bins = np.arange(0, int(tara["days_since_irrig"].max())+2)
    centers, meds = [], []
    for b in range(len(bins)-1):
        mask = ((tara["days_since_irrig"] >= bins[b]) &
                (tara["days_since_irrig"] <  bins[b+1]))
        if mask.sum() >= 3:
            centers.append((bins[b]+bins[b+1])/2)
            meds.append(tara.loc[mask,"LE_corr"].median())
    if centers:
        ax.plot(centers, meds, "k-", lw=2.5, alpha=0.85,
                 label="Median per 1d bin", zorder=5)

    ax.axvline(7, color="gray", lw=1, ls="--", alpha=0.6,
                label="7d (典型的灌漑間隔)")
    ax.set_xlabel("Days since last irrigation")
    ax.set_ylabel("LE_corr (W/m²)")
    ax.set_title("[Tarazona Growing] LE vs days since irrigation\n"
                 "中央値線が右下がり=灌漑依存 / 平坦=深根",
                 fontweight="bold")
    ax.legend(fontsize=9, loc="best")
    ax.grid(alpha=0.25)

    plt.tight_layout()
    fp = save_dir / "fig02_le_vs_days_since.png"
    plt.savefig(fp, dpi=130, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [save] {fp}")


def plot_sds_by_bucket(sds_df, save_dir):
    """★ 主結果: SDS を経過日数バケット別に表示"""
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.5))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle("★ SDS by days_since_irrigation — 灌漑依存 vs 深根の判定\n"
                 "右に行くほどバーが伸びる=灌漑依存 / 平坦=深根",
                 fontsize=13, fontweight="bold")

    for ax, var in zip(axes, ["LE_corr","EF_corr","ET"]):
        sub = sds_df[sds_df["var"]==var].sort_values("bucket_order")
        if sub.empty:
            ax.set_title(f"{var}: no data"); continue
        x = np.arange(len(sub))
        for i, row in enumerate(sub.itertuples()):
            col = BUCKET_COL[int(row.bucket_order)]
            val = row.sds if not np.isnan(row.sds) else 0
            ax.bar(i, val, color=col, alpha=0.88, edgecolor="black", lw=1.0)
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
        ax.set_xticks(x)
        ax.set_xticklabels(
            [f"{r.bucket}\nn_n={int(r.n_n)} n_s={int(r.n_s)}"
             for r in sub.itertuples()], fontsize=8.5)
        ax.set_ylabel("SDS")
        ax.set_title(var, fontweight="bold")
        ax.grid(axis="y", alpha=0.25)

    plt.tight_layout()
    fp = save_dir / "fig03_sds_by_days_since.png"
    plt.savefig(fp, dpi=140, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [save] {fp}")


# ================================================================
# 5. MAIN
# ================================================================

def main():
    print("="*60)
    print("解析A v13 — 灌漑データによる深根仮説の真偽判定")
    print("="*60)

    df = load_and_merge()
    df = compute_days_since_irrig(df)
    df["bucket"] = df["days_since_irrig"].apply(bucket_of)
    df.to_csv(OUT_DIR / "v13_irrig_classified.csv", index=False)

    # 灌漑統計の確認
    tara = df[df["site"]=="Tarazona"]
    if "Irrig_mm" in tara.columns:
        irrig_days = tara[tara["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD]
        if len(irrig_days):
            print(f"\n[灌漑統計]")
            print(f"  灌漑日数 : {len(irrig_days)}")
            print(f"  灌漑期間 : {irrig_days['date'].min().date()} 〜 "
                  f"{irrig_days['date'].max().date()}")
            print(f"  灌漑量(mm/event): "
                  f"med={irrig_days['Irrig_mm'].median():.1f}  "
                  f"max={irrig_days['Irrig_mm'].max():.1f}")
            # 月別灌漑頻度
            irrig_days["month"] = irrig_days["date"].dt.month
            mc = irrig_days["month"].value_counts().sort_index()
            print(f"  月別頻度 : " +
                  "  ".join(f"{m}月={n}" for m, n in mc.items()))

    # バケット別 SDS 計算 (Tarazona growing only)
    print(f"\n[Tarazona Growing] SDS by days_since_irrig bucket")
    print(f"  → SDS が経過日数で平坦なら深根、上昇なら灌漑依存")
    rows = []
    tg = df[(df["site"]=="Tarazona") & df["is_growing"]]
    for var in ["LE_corr","EF_corr","ET"]:
        if var not in tg.columns: continue
        for i, (label, lo, hi) in enumerate(BUCKETS):
            sub = tg[tg["bucket"]==label]
            n_arr = sub.loc[sub["drought_type"]=="normal", var].dropna().values
            s_arr = sub.loc[sub["drought_type"]=="soil dry", var].dropna().values
            res = bootstrap_sds(n_arr, s_arr, var)
            res.update(dict(var=var, bucket=label, bucket_order=i))
            rows.append(res)
            ci = (f"[{res['ci_lo']:+.2f}, {res['ci_hi']:+.2f}]"
                  if not np.isnan(res["ci_lo"]) else "[CI算出不可]")
            print(f"  {var:9s} {label:14s}: SDS={res['sds']:+.3f}  {ci}  "
                  f"(n_n={res['n_n']}, n_s={res['n_s']})")
    sds_df = pd.DataFrame(rows)
    sds_df.to_csv(OUT_DIR / "v13_sds_by_days_since.csv", index=False)

    # 可視化
    print(f"\n--- 可視化 ---")
    plot_timeseries(df, OUT_DIR)
    plot_le_vs_days(df, OUT_DIR)
    plot_sds_by_bucket(sds_df, OUT_DIR)

    # ★ 最終判定
    print("\n"+"="*60)
    print("★ 判定 (LE_corr ベース)")
    print("="*60)
    le = sds_df[sds_df["var"]=="LE_corr"].sort_values("bucket_order")
    valid = le["sds"].notna()
    if valid.sum() >= 2:
        x = le.loc[valid,"bucket_order"].values.astype(float)
        y = le.loc[valid,"sds"].values
        slope, intercept = np.polyfit(x, y, 1)
        r = np.corrcoef(x, y)[0,1]
        print(f"  SDS の経過日数バケット番号に対する傾き = {slope:+.4f}")
        print(f"  相関係数 r = {r:+.3f}")
        print(f"  各バケット SDS: " +
              "  ".join(f"{b[:7]}={v:+.2f}"
                        for b, v in zip(le.loc[valid,"bucket"], y)))
        if slope > 0.07 and r > 0.5:
            print("  ⇒ ★ SDS が経過日数で明確に上昇 → **灌漑依存**(深根仮説却下)")
        elif slope < -0.05:
            print("  ⇒ SDS が逆相関(不自然)→ 別の交絡要因疑い")
        else:
            print("  ⇒ ★ SDS が経過日数とほぼ無関係 → **深根仮説を支持**"
                  "(灌漑離脱後も LE 維持)")
    else:
        print("  バケット間で比較するには SDS のサンプルが不足")

    print(f"\n[done] outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
