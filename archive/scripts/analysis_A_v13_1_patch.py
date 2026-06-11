"""
v13.1 patch: LE vs days_since_irrigation の交絡除去版
=========================================================
v13 fig02 の問題: is_growing フラグ内に「灌漑オフシーズン」が混入する可能性
                  → days_since_irrig が大きい領域は単に非灌漑シーズンの低LE
                  を見ているだけかもしれない

修正: 月別灌漑頻度 >= 2 の月のみに限定して再プロット
      → 灌漑アクティブ期間内での純粋な経過日数効果を抽出

出力:
  output_analysis_A_v13/fig02_le_vs_days_since_filtered.png
    左パネル: 修正前(v13 オリジナル、is_growing のみ)
    右パネル: 修正後(irrig_active_month==True のみ)
"""

import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

INPUT_DIR      = Path("./")
PARQUET        = INPUT_DIR / "daily_classified_v4.parquet"
TARA_DAILY_CSV = Path(
    "/home/shion-nagamine/Dataset/Eddy data in Spain/"
    "Daily_Summary_Filtered_forPred_ActEne26.csv")
OUT_DIR        = Path("./output_analysis_A_v13")
OUT_DIR.mkdir(exist_ok=True)

IRRIG_THRESHOLD       = 0.5  # mm/day, 灌漑イベント判定
MIN_IRRIG_PER_MONTH   = 2    # 月別灌漑イベント数(これ以上の月を active 判定)

DROUGHT_COL = {"normal":"#2196F3","soil dry":"#FF9800",
               "atm dry":"#4CAF50","compound":"#F44336"}


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
    return tara


def compute_days_since_irrig(df):
    df = df.sort_values("date").reset_index(drop=True).copy()
    last = None; out = []
    for _, row in df.iterrows():
        if pd.notna(row.get("Irrig_mm")) and row["Irrig_mm"] > IRRIG_THRESHOLD:
            last = row["date"]; out.append(0)
        elif last is None:
            out.append(np.nan)
        else:
            out.append((row["date"] - last).days)
    df["days_since_irrig"] = out
    return df


def tag_irrig_active_month(df, min_per_month=MIN_IRRIG_PER_MONTH):
    """月別の灌漑イベント数を計算し、>= min_per_month を active と判定"""
    df = df.copy()
    df["year_month"] = df["date"].dt.to_period("M")
    monthly = (df.assign(is_irrig=df["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD)
                 .groupby("year_month")["is_irrig"].sum())
    active_months = monthly[monthly >= min_per_month].index
    df["irrig_active_month"] = df["year_month"].isin(active_months)
    df = df.drop(columns="year_month")
    return df, active_months


def plot_one_panel(ax, sub, title, show_legend=False):
    """1枚の散布図 + 中央値線を描く共通関数"""
    if len(sub) == 0:
        ax.set_title(f"{title}\n(no data)")
        return
    for dt, col in DROUGHT_COL.items():
        s = sub[sub["drought_type"]==dt]
        ax.scatter(s["days_since_irrig"], s["LE_corr"],
                   c=col, alpha=0.55, s=24, label=f"{dt} (n={len(s)})",
                   edgecolors="none")
    # 中央値線(1日刻み、>=3 サンプル)
    if sub["days_since_irrig"].notna().any():
        max_d = int(sub["days_since_irrig"].max())
        bins = np.arange(0, max_d + 2)
        centers, meds = [], []
        for b in range(len(bins)-1):
            mask = ((sub["days_since_irrig"] >= bins[b]) &
                    (sub["days_since_irrig"] <  bins[b+1]))
            if mask.sum() >= 3:
                centers.append((bins[b]+bins[b+1])/2)
                meds.append(sub.loc[mask, "LE_corr"].median())
        if centers:
            ax.plot(centers, meds, "k-", lw=2.5, alpha=0.85,
                    label="Median (1d bin)", zorder=5)
    ax.axvline(7, color="gray", lw=1, ls="--", alpha=0.6, label="7d 目安")
    ax.set_xlabel("Days since last irrigation")
    ax.set_ylabel("LE_corr (W/m²)")
    ax.set_title(title, fontweight="bold")
    if show_legend:
        ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.25)


def main():
    print("="*60)
    print("v13.1 patch — fig02 灌漑シーズン交絡除去")
    print("="*60)

    tara = load_and_merge()
    tara = compute_days_since_irrig(tara)
    tara, active_months = tag_irrig_active_month(tara)

    print(f"  全データ: {len(tara)} 日")
    print(f"  灌漑アクティブ月数: {len(active_months)} 月")
    print(f"  active な月: {[str(m) for m in active_months[:10]]}...")

    # フィルタ適用
    before = tara[tara["is_growing"] & tara["days_since_irrig"].notna() &
                   tara["LE_corr"].notna()]
    after  = before[before["irrig_active_month"]]

    print(f"\n  修正前(is_growing のみ)         : n = {len(before)}")
    print(f"  修正後(irrig_active_month=True): n = {len(after)}")
    print(f"  除外された日数: {len(before) - len(after)}")
    if len(before) > 0:
        max_days_before = int(before["days_since_irrig"].max())
        max_days_after  = int(after["days_since_irrig"].max()) if len(after) else 0
        print(f"  days_since_irrig max: before={max_days_before}d → after={max_days_after}d")

    # サイドバイサイドプロット
    fig, axes = plt.subplots(1, 2, figsize=(20, 7), sharey=True)
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle(
        "[Tarazona Growing] LE vs days_since_irrigation\n"
        f"灌漑非アクティブ月(灌漑イベント<{MIN_IRRIG_PER_MONTH}/月)の混入を除去",
        fontsize=12, fontweight="bold")

    plot_one_panel(
        axes[0], before,
        title=f"修正前: is_growing のみ  (n={len(before)})\n"
              f"max days_since={int(before['days_since_irrig'].max()) if len(before) else 0}d",
        show_legend=True)
    plot_one_panel(
        axes[1], after,
        title=f"修正後: irrig_active_month=True 限定  (n={len(after)})\n"
              f"max days_since={int(after['days_since_irrig'].max()) if len(after) else 0}d",
        show_legend=True)

    plt.tight_layout()
    fp = OUT_DIR / "fig02_le_vs_days_since_filtered.png"
    plt.savefig(fp, dpi=140, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"\n  [save] {fp}")

    # 数値サマリー: 灌漑非アクティブ月の混入度
    print("\n--- 解釈ヒント ---")
    if len(before) - len(after) > 100:
        print(f"  ⚠ {len(before)-len(after)} 日が「灌漑非アクティブ月の長間隔」として")
        print(f"    修正前プロットに混入していた。これらは days_since_irrig が大きく、")
        print(f"    かつ LE が phenology 由来で低い → 中央値線を右下がりに見せる artifact。")
    else:
        print(f"  混入は少量 ({len(before)-len(after)}日)、修正前後で大きな違いは出ない可能性")

    # 中央値線の比較
    print("\n--- 中央値線の比較(1日刻み、>=3 サンプルのみ) ---")
    for label, sub in [("修正前", before), ("修正後", after)]:
        if len(sub) == 0: continue
        max_d = int(sub["days_since_irrig"].max())
        for d_lo, d_hi, name in [(0,1,"0-1d"), (2,3,"2-3d"),
                                   (4,7,"4-7d"), (8,999,"8+d")]:
            mask = ((sub["days_since_irrig"] >= d_lo) &
                    (sub["days_since_irrig"] <= d_hi))
            if mask.sum() >= 5:
                m = sub.loc[mask, "LE_corr"].median()
                print(f"  {label} {name:6s}: LE_corr median = {m:6.1f} W/m²  (n={mask.sum()})")

    print(f"\n[done]")


if __name__ == "__main__":
    main()
