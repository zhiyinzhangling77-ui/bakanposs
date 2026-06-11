"""
================================================================
解析A v12 : 4象限分類からの干ばつ応答指標(深根/浅根仮説の検定)
================================================================

【目的】
  v4 の 4分類済日次データを入力に、サイトごとの「土壌乾燥に対する
  蒸散維持力」を 4 指標で定量化し、深根仮説を統計的に検定する。

【入力】(v4 が出力)
  daily_classified_v4.parquet : 日次×分類済
  closure_slopes_v4.json      : エネルギー閉合スロープ

【4指標の定義と意味】
  SDS = (LE_normal − LE_soil_dry) / LE_normal     表層SWC感受性
        ≈ 0  → 土壌乾燥でLE減らない (深根)
        > 0.3→ 強く減る (浅根)
  VAS = (LE_atm_dry − LE_normal) / LE_normal      VPD活性化
  DSO = (LE_compound − LE_soil_dry) / LE_soil_dry 土壌乾燥下のVPD駆動能
  CompoundDrop = (LE_normal − LE_compound) / LE_normal  複合干ばつ感受性
                 (負なら compound > normal = VPD効果が支配的)

【出力】./output_analysis_A_v12/
  fig01_headline.png            ★ 主結果(SDS, growing, LE)
  fig02_indices_LE_corr.png     LE 4指標 × 4区分
  fig03_indices_EF_corr.png     EF 同上
  fig04_indices_ET.png          ET 同上
  fig05_indices_Bowen_corr.png  Bowen 同上
  fig06_distributions_LE_corr.png   LE boxplot (site × season)
  fig07_distributions_EF_corr.png   EF 同上
  fig08_distributions_ET.png        ET 同上
  + CSV (indices / pairwise tests / class counts)
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

warnings.filterwarnings("ignore")


# ================================================================
#                   設定 & 物理的に意味のある閾値
# ================================================================

INPUT_DIR   = Path("./")
PARQUET     = INPUT_DIR / "daily_classified_v4.parquet"
SLOPES_JSON = INPUT_DIR / "closure_slopes_v4.json"
OUT_DIR     = Path("./output_analysis_A_v12")
OUT_DIR.mkdir(exist_ok=True)

# --- bootstrap / 検定のパラメータ ---
N_BOOT          = 5000     # 95%CI 安定化のため
CI_PCT          = (2.5, 97.5)
MIN_N_PER_CLASS = 5        # 中央値計算の最低n

# --- 比の暴走を防ぐためのガード ---
# 変数ごとに物理スケールが違う(LE: W/m², EF: 0-1, Bowen: ~1, ET: mm/day)ので
# 一律閾値ではなく変数別に分母最小値を設定する
DENOM_FLOORS = {
    "LE_corr"    : 5.0,    # W/m²
    "EF_corr"    : 0.05,   # 無次元 0-1
    "Bowen_corr" : 0.1,    # 無次元 ~1
    "ET"         : 0.3,    # mm/day
}
DENOM_FLOOR_DEFAULT = 0.05
RATIO_CLIP = 10.0          # 計算後の比の絶対値クリップ
BOWEN_CLIP = 20.0          # |Bowen|>20 は LE near-zero 由来の極値として除外

# --- 視覚化(配色は v4 と整合) ---
CLASSES   = ["normal", "soil dry", "atm dry", "compound"]
CLASS_COL = {"normal":"#2196F3", "soil dry":"#FF9800",
             "atm dry":"#4CAF50", "compound":"#F44336"}
SITE_COL  = {"Oran":"#E85D04", "Tarazona":"#1D9E75"}
VARS      = ["LE_corr", "EF_corr", "Bowen_corr", "ET"]

# 指標の説明(figure 上に表示する短い解釈)
INDEX_DESC = {
    "SDS"          : "(N − S) / N   ≈0 = 表層SWC非依存(深根)",
    "VAS"          : "(A − N) / N   高=VPDで蒸散ブースト",
    "DSO"          : "(C − S) / S   高=土壌乾燥下でもVPD駆動",
    "CompoundDrop" : "(N − C) / N   負=VPD効果が支配的",
}


# ================================================================
#                       入出力ヘルパ
# ================================================================

def load_inputs():
    """v4 が吐いた中間ファイルを読む。無ければ親切なエラーメッセージを返す。"""
    if not PARQUET.exists() or not SLOPES_JSON.exists():
        raise FileNotFoundError(
            f"中間ファイルがありません。v4 末尾に下記2行を追加してください:\n"
            f"  daily.to_parquet('{PARQUET}')\n"
            f"  json.dump(closure_slopes, open('{SLOPES_JSON}','w'))\n"
        )
    daily  = pd.read_parquet(PARQUET)
    slopes = json.loads(SLOPES_JSON.read_text())
    return daily, slopes


# ================================================================
#                       指標計算ロジック
# ================================================================

def safe_ratio(numer, denom, var):
    """
    比 numer/denom を計算するが、分母が物理的に意味のある最小値を
    下回るとき・結果の絶対値が大きすぎるときは NaN を返す。
    """
    floor = DENOM_FLOORS.get(var, DENOM_FLOOR_DEFAULT)
    if denom is None or np.isnan(denom) or abs(denom) < floor:
        return np.nan
    r = numer / denom
    return np.nan if abs(r) > RATIO_CLIP else r


def indices_from_medians(meds, var):
    """4 クラスの中央値辞書から 4 指標を一括計算"""
    n, s, a, c = (meds.get(k) for k in ["normal","soil dry","atm dry","compound"])
    has = lambda *xs: not any(np.isnan(x) for x in xs)
    return {
        "SDS"          : safe_ratio(n - s, n, var) if has(n, s) else np.nan,
        "VAS"          : safe_ratio(a - n, n, var) if has(a, n) else np.nan,
        "DSO"          : safe_ratio(c - s, s, var) if has(c, s) else np.nan,
        "CompoundDrop" : safe_ratio(n - c, n, var) if has(n, c) else np.nan,
    }


def class_arrays(sub_df, var):
    """サブセットからクラス別配列を取り出す。Bowen のみ極値を除外。"""
    out = {c: sub_df.loc[sub_df["drought_type"]==c, var].dropna().values
           for c in CLASSES}
    if var == "Bowen_corr":
        out = {c: v[np.abs(v) < BOWEN_CLIP] for c, v in out.items()}
    return out


def bootstrap_indices(data, var, n_boot=N_BOOT, seed=42):
    """各クラス内でリサンプリングして指標分布を生成し、点推定 + 95%CI を返す"""
    rng = np.random.default_rng(seed)

    # 点推定: 観測中央値からそのまま計算
    point_meds = {c: (np.median(v) if len(v) >= MIN_N_PER_CLASS else np.nan)
                   for c, v in data.items()}
    point = indices_from_medians(point_meds, var)

    # bootstrap 反復: 各クラスを独立に復元抽出 → 中央値 → 指標
    boots = {k: [] for k in point.keys()}
    for _ in range(n_boot):
        b_meds = {}
        for c, v in data.items():
            if len(v) < MIN_N_PER_CLASS:
                b_meds[c] = np.nan
            else:
                b_meds[c] = np.median(rng.choice(v, size=len(v), replace=True))
        for k, val in indices_from_medians(b_meds, var).items():
            if not np.isnan(val):
                boots[k].append(val)

    # CI: bootstrap 分布の 2.5/97.5 パーセンタイル
    out = {}
    for k, vals in boots.items():
        if len(vals) < n_boot * 0.1:
            out[k] = (point[k], np.nan, np.nan, len(vals))
        else:
            lo, hi = np.percentile(vals, CI_PCT)
            out[k] = (point[k], lo, hi, len(vals))
    return out, point_meds


def pairwise_tests(data, var):
    """全6ペアで Mann-Whitney(両側)+ rank-biserial 効果量"""
    pairs = [("normal","soil dry"),("normal","atm dry"),("normal","compound"),
             ("soil dry","compound"),("soil dry","atm dry"),("atm dry","compound")]
    rows = []
    for a, b in pairs:
        da, db = data[a], data[b]
        if len(da) < 3 or len(db) < 3:
            continue
        try:
            u, p = stats.mannwhitneyu(da, db, alternative="two-sided")
            r_rb = 1 - 2 * u / (len(da) * len(db))
        except ValueError:
            u, p, r_rb = np.nan, np.nan, np.nan
        rows.append(dict(class_a=a, class_b=b, var=var,
                         n_a=len(da), n_b=len(db),
                         med_a=float(np.median(da)), med_b=float(np.median(db)),
                         U=u, p=p, rank_biserial=r_rb))
    return rows


# ================================================================
#                        可視化(8枚の図)
# ================================================================

def _sig_marker(p):
    if np.isnan(p): return ""
    return "★★★" if p < 0.001 else ("★★" if p < 0.01 else
            ("★" if p < 0.05 else "n.s."))


# ---- fig01: ヘッドライン -----------------------------------------

def plot_headline(summary, save_dir):
    """
    SDS だけを Oran vs Tarazona で対比する1枚。
    研究の主結論をこの1枚で読み取れることを目指す。
    """
    sub = summary[(summary["var"]=="LE_corr") &
                   (summary["index"]=="SDS") &
                   (summary["season"]=="growing")].sort_values("site")
    if sub.empty:
        return

    fig, ax = plt.subplots(figsize=(9, 6.5))
    fig.patch.set_facecolor("#F8F8F8")

    x = np.arange(len(sub))
    for i, row in enumerate(sub.itertuples()):
        col = SITE_COL[row.site]
        ax.bar(i, row.point, color=col, alpha=0.9, edgecolor="black", lw=1.4,
               width=0.55)
        if not np.isnan(row.ci_lo):
            ax.errorbar(i, row.point,
                        yerr=[[row.point - row.ci_lo], [row.ci_hi - row.point]],
                        color="black", capsize=14, lw=2.0)
        # 値ラベル
        ax.text(i, row.point + 0.018,
                f"{row.point:+.3f}\n[{row.ci_lo:+.2f}, {row.ci_hi:+.2f}]",
                ha="center", va="bottom", fontsize=11, fontweight="bold")
        # 解釈ラベル
        verdict = ("浅根 (表層依存)" if row.point > 0.15 else
                    "深根 (非依存)"   if abs(row.point) <= 0.15 else "?")
        ax.text(i, -0.06, verdict, ha="center", va="top",
                fontsize=12, fontweight="bold", color=col,
                transform=ax.get_xaxis_transform())

    ax.axhline(0, color="gray", lw=1.0, ls="--", label="SDS = 0  (深根の理論値)")
    ax.axhspan(-0.15, 0.15, color="#1D9E75", alpha=0.08,
                label="深根判定領域 (|SDS| ≤ 0.15)")

    ax.set_xticks(x)
    ax.set_xticklabels([row.site for row in sub.itertuples()],
                       fontsize=13, fontweight="bold")
    ax.set_ylabel("SDS  =  (LE_normal − LE_soil_dry) / LE_normal", fontsize=12)
    ax.set_title(
        "★ Soil Drought Sensitivity (SDS) — Growing season\n"
        "Tarazona の 95%CI が 0 を含む = 表層SWC低下に対し LE 不変 = 深根の証拠",
        fontsize=12, fontweight="bold", pad=12)
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(axis="y", alpha=0.25)
    ax.set_ylim(-0.2, max(0.5, sub["ci_hi"].max() * 1.25))

    plt.tight_layout()
    fp = save_dir / "fig01_headline.png"
    plt.savefig(fp, dpi=160, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [save] {fp}")


# ---- fig02-05: 変数別の指標 bar chart ----------------------------

def plot_indices(summary, var, save_dir, fig_idx):
    """4指標 × (Oran/Tarazona × growing/non_growing) を 1×4 で並べる"""
    sub = summary[summary["var"] == var]
    if sub.empty:
        return
    fig, axes = plt.subplots(1, 4, figsize=(20, 5.8))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle(f"{var} — drought response indices  "
                 f"(bars = 中央値, error bars = 95% bootstrap CI)",
                 fontsize=13, fontweight="bold")

    indices = ["SDS", "VAS", "DSO", "CompoundDrop"]
    for ax, idx in zip(axes, indices):
        d = sub[sub["index"] == idx].copy()
        if d.empty:
            ax.set_title(f"{idx}\n(no data)"); continue

        # キーで並び替え: site→season(growingを左に)
        d["sort_key"] = d["site"] + d["season"].map({"growing":"0","non_growing":"1"})
        d = d.sort_values("sort_key")
        x = np.arange(len(d))

        for i, row in enumerate(d.itertuples()):
            col = SITE_COL[row.site]
            alpha = 0.92 if row.season == "growing" else 0.40
            ax.bar(i, row.point if not np.isnan(row.point) else 0,
                   color=col, alpha=alpha, edgecolor="black", lw=1.0)
            if not np.isnan(row.ci_lo):
                ax.errorbar(i, row.point,
                             yerr=[[row.point-row.ci_lo],[row.ci_hi-row.point]],
                             color="black", capsize=6, lw=1.2)
            if not np.isnan(row.point):
                ax.text(i, row.point + (0.015 if row.point>=0 else -0.015),
                        f"{row.point:+.2f}",
                        ha="center",
                        va="bottom" if row.point>=0 else "top",
                        fontsize=9, fontweight="bold")

        ax.axhline(0, color="gray", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r.site}\n{r.season}" for r in d.itertuples()],
                           fontsize=8.5)
        ax.set_ylabel(idx, fontsize=11)
        ax.set_title(f"{idx}\n{INDEX_DESC[idx]}", fontsize=10)
        ax.grid(axis="y", alpha=0.22)

    plt.tight_layout()
    fp = save_dir / f"fig{fig_idx:02d}_indices_{var}.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [save] {fp}")


# ---- fig06-09: 変数別の boxplot ----------------------------------

def plot_distributions(daily, var, save_dir, fig_idx, tests_df):
    """
    site × season の 2×2 boxplot。
    各 box の中央値を直接ラベルし、normal vs soil_dry の p 値を annotate。
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle(f"{var} distributions by drought class\n"
                 f"(中央値ラベル / N→S は normal vs soil_dry 比較)",
                 fontsize=13, fontweight="bold")

    for r, site in enumerate(["Oran","Tarazona"]):
        for c, is_g in enumerate([True, False]):
            ax = axes[r, c]
            sub = daily[(daily["site"]==site) & (daily["is_growing"]==is_g)]
            arrs = [class_arrays(sub, var)[cls] for cls in CLASSES]
            ns   = [len(a) for a in arrs]
            if all(n == 0 for n in ns):
                ax.set_title(f"[{site}] no data"); continue

            bp = ax.boxplot(arrs, patch_artist=True, widths=0.6, showfliers=True,
                             medianprops=dict(color="white", lw=2.4),
                             flierprops=dict(marker="o", markersize=3,
                                              alpha=0.4))
            for p, cls in zip(bp["boxes"], CLASSES):
                p.set_facecolor(CLASS_COL[cls]); p.set_alpha(0.78)

            # 中央値ラベル(boxの上)
            ymax = max((np.percentile(a, 95) if len(a) else 0) for a in arrs)
            for i, a in enumerate(arrs, start=1):
                if len(a) >= MIN_N_PER_CLASS:
                    ax.text(i, np.median(a), f"{np.median(a):.2f}",
                            ha="center", va="center", fontsize=8.5,
                            fontweight="bold", color="black",
                            bbox=dict(boxstyle="round,pad=0.15",
                                      facecolor="white", alpha=0.85,
                                      edgecolor="none"))

            ax.set_xticklabels([f"{c}\nn={n}" for c, n in zip(CLASSES, ns)],
                                fontsize=9)
            season = "Growing" if is_g else "Non-growing"
            ax.set_title(f"[{site}]  {season}", fontweight="bold", fontsize=11)
            ax.set_ylabel(var); ax.grid(axis="y", alpha=0.22)

            # N→S の p 値を annotate
            t = tests_df[(tests_df["site"]==site) &
                          (tests_df["season"]==("growing" if is_g else "non_growing")) &
                          (tests_df["var"]==var) &
                          (tests_df["class_a"]=="normal") &
                          (tests_df["class_b"]=="soil dry")]
            if not t.empty:
                p = t.iloc[0]["p"]; sig = _sig_marker(p)
                ax.annotate(f"N→S: p={p:.2e}  {sig}",
                            xy=(0.5, 0.97), xycoords="axes fraction",
                            ha="center", va="top",
                            fontsize=10, fontweight="bold",
                            color="#7C3AED",
                            bbox=dict(boxstyle="round,pad=0.3",
                                      facecolor="white", alpha=0.85,
                                      edgecolor="#7C3AED"))

    plt.tight_layout()
    fp = save_dir / f"fig{fig_idx:02d}_distributions_{var}.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
#                            MAIN
# ================================================================

def main():
    print("=" * 60)
    print("解析A v12 — 4象限分類からの干ばつ応答指標")
    print("=" * 60)

    daily, slopes = load_inputs()
    print(f"[input] daily rows = {len(daily)}")
    print(f"[input] closure slopes = {slopes}")

    # --------------------------------------------------------
    # Step 1. サイト × 季節 × 変数 で指標 + bootstrap CI を計算
    # --------------------------------------------------------
    summary_rows, test_rows, count_rows = [], [], []

    for site in sorted(daily["site"].unique()):
        for is_g in [True, False]:
            season = "growing" if is_g else "non_growing"
            sub = daily[(daily["site"]==site) & (daily["is_growing"]==is_g)]
            counts = {c: int((sub["drought_type"]==c).sum()) for c in CLASSES}
            count_rows.append(dict(site=site, season=season, **counts))
            print(f"\n[{site} / {season}] counts: {counts}")

            for var in VARS:
                if var not in sub.columns:
                    continue
                data = class_arrays(sub, var)
                idx_results, meds = bootstrap_indices(data, var)

                med_str = "  ".join(
                    f"{c}={meds[c]:.2f}" if not np.isnan(meds[c]) else f"{c}=NaN"
                    for c in CLASSES)
                print(f"  {var}: medians {med_str}")
                # 休眠期警告: LE_normal が低すぎると指標は phenology 由来になる
                if var == "LE_corr" and not np.isnan(meds["normal"]) \
                        and meds["normal"] < 30:
                    print(f"    ⚠ LE_normal={meds['normal']:.1f} < 30 W/m² "
                          f"→ 休眠/裸地期、指標解釈に注意")

                for idx_name, (pt, lo, hi, n_b) in idx_results.items():
                    summary_rows.append(dict(
                        site=site, season=season, var=var, index=idx_name,
                        point=pt, ci_lo=lo, ci_hi=hi, n_boot=n_b,
                        med_normal=meds["normal"], med_soil_dry=meds["soil dry"],
                        med_atm_dry=meds["atm dry"], med_compound=meds["compound"],
                    ))
                    if not np.isnan(pt):
                        ci = (f"[{lo:+.3f}, {hi:+.3f}]"
                               if not np.isnan(lo) else "[CI算出不可]")
                        print(f"    {idx_name:13s}: {pt:+.3f}  {ci}")

                for t in pairwise_tests(data, var):
                    t.update(dict(site=site, season=season))
                    test_rows.append(t)

    summary_df = pd.DataFrame(summary_rows)
    tests_df   = pd.DataFrame(test_rows)
    counts_df  = pd.DataFrame(count_rows)

    summary_df.to_csv(OUT_DIR / "v12_indices.csv", index=False)
    tests_df.to_csv(OUT_DIR / "v12_pairwise_tests.csv", index=False)
    counts_df.to_csv(OUT_DIR / "v12_class_counts.csv", index=False)
    print(f"\n[save] CSV 3 files in {OUT_DIR}/")

    # --------------------------------------------------------
    # Step 2. 8枚の図を生成
    # --------------------------------------------------------
    print(f"\n--- 可視化 (8枚) ---")
    plot_headline(summary_df, OUT_DIR)
    for i, var in enumerate(VARS, start=2):
        plot_indices(summary_df, var, OUT_DIR, fig_idx=i)
    for i, var in enumerate(VARS, start=6):
        plot_distributions(daily, var, OUT_DIR, fig_idx=i, tests_df=tests_df)

    # --------------------------------------------------------
    # Step 3. 主結果サマリー
    # --------------------------------------------------------
    print("\n" + "=" * 60)
    print("★ 主結果(Growing season, LE_corr ベース)")
    print("=" * 60)
    g = summary_df[(summary_df["var"]=="LE_corr") &
                    (summary_df["season"]=="growing")]
    for idx_name in ["SDS","VAS","DSO","CompoundDrop"]:
        d = g[g["index"]==idx_name]
        if d.empty: continue
        print(f"\n[{idx_name}]  {INDEX_DESC[idx_name]}")
        for _, row in d.iterrows():
            ci = (f"[{row['ci_lo']:+.3f}, {row['ci_hi']:+.3f}]"
                   if not np.isnan(row["ci_lo"]) else "[CI算出不可]")
            print(f"  {row['site']:9s}: {row['point']:+.3f}  {ci}")

    # 仮説検証(SDS で 0-overlap を判定)
    sds = g[g["index"]=="SDS"].set_index("site")
    if {"Oran","Tarazona"}.issubset(sds.index):
        oran = sds.loc["Oran"]; tara = sds.loc["Tarazona"]
        print(f"\n→ 仮説検証(SDS, growing season):")
        print(f"   Oran     : {oran['point']:+.3f}  "
              f"[{oran['ci_lo']:+.2f}, {oran['ci_hi']:+.2f}]")
        print(f"   Tarazona : {tara['point']:+.3f}  "
              f"[{tara['ci_lo']:+.2f}, {tara['ci_hi']:+.2f}]")
        if (tara["ci_hi"] < oran["ci_lo"]):
            print(f"   ⇒ ★ CI 重ならず: Tarazona の SDS が有意に低い = 深根仮説支持")
        elif tara["ci_lo"] <= 0 <= tara["ci_hi"] and oran["ci_lo"] > 0:
            print(f"   ⇒ ★ Tarazona の CI は 0 を含む(土壌乾燥に非感受) = 深根仮説支持")
        else:
            print(f"   ⇒ 統計的に決着せず — n や閾値の見直しを検討")

    print(f"\n[done] outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
