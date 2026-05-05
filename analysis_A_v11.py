"""
解析A v11: SDS / VAS / DSO 指標 — 4象限分類からの直接定量化
============================================================

【目的】
  ユーザー提示の日変化図(v4スクリプト)で既に視覚的に明らかな結論
  「Tarazona は soil_dry でも LE が落ちない、Oran は落ちる」
  を、bootstrap 信頼区間付きの数値指標として定量化する。

【入力】(v4 が出力する中間ファイル)
  - daily_classified_v4.parquet : 日次分類済み(site, is_growing, drought_type, LE_corr, ...)
  - closure_slopes_v4.json      : エネルギー閉合スロープ {"Oran": 0.74, ...}

【指標定義】(各サイト×生育期/非生育期 別)
  SDS = (LE_normal - LE_soil_dry) / LE_normal     表層土壌乾燥感受性
        ≈0  → 表層SWCに非依存(深根)
        >0.3 → 表層依存(浅根)
  VAS = (LE_atm_dry - LE_normal) / LE_normal      VPD活性化感受性
        高い → 大気乾燥でむしろ蒸散を増やす(strong stomatal-VPD coupling)
  DSO = (LE_compound - LE_soil_dry) / LE_soil_dry 土壌乾燥下のVPD上書き能力
        高い → 土壌乾燥でも大気が乾けばさらに蒸散できる = 深層水アクセスの強い証拠
  CompoundDrop = (LE_normal - LE_compound) / LE_normal  複合干ばつ下のLE減少率

【出力】./output_analysis_A_v11/
  - v11_indices.csv         : 全指標の点推定 + 95%CI
  - v11_pairwise_tests.csv  : Mann-Whitney(全6ペア) + rank-biserial
  - v11_class_counts.csv    : サイト×季節×classのn数
  - v11_indices_LE.png      : 主指標(LE_corr)の bar chart
  - v11_indices_EF.png      : 副指標(EF_corr)の bar chart
  - v11_distributions_*.png : class別 boxplot
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
# CONFIG
# ================================================================
INPUT_DIR   = Path("./")
PARQUET     = INPUT_DIR / "daily_classified_v4.parquet"
SLOPES_JSON = INPUT_DIR / "closure_slopes_v4.json"

OUT_DIR = Path("./output_analysis_A_v11")
OUT_DIR.mkdir(exist_ok=True)

N_BOOT          = 5000     # bootstrap 反復数
CI_PCT          = (2.5, 97.5)
MIN_N_PER_CLASS = 5        # クラス内最小サンプル(これ未満はNaN)
DENOM_FLOOR     = 5.0      # 比の分母最小値(W/m²等)。これ未満ならその iteration を破棄
RATIO_CLIP      = 10.0     # 比の絶対値クリップ(暴走防止)

CLASSES   = ["normal", "soil dry", "atm dry", "compound"]
COLORS    = {"normal":"#2196F3","soil dry":"#FF9800",
             "atm dry":"#4CAF50","compound":"#F44336"}
SITE_COL  = {"Oran":"#E85D04", "Tarazona":"#1D9E75"}
VARS      = ["LE_corr", "EF_corr", "Bowen_corr", "ET"]

INDEX_INTERP = {
    "SDS"          : "(LE_normal − LE_soil_dry) / LE_normal\n→ 0 = 表層SWC非依存(深根)",
    "VAS"          : "(LE_atm_dry − LE_normal) / LE_normal\n→ 高い = VPDで蒸散ブースト",
    "DSO"          : "(LE_compound − LE_soil_dry) / LE_soil_dry\n→ 高い = 土壌乾燥下でもVPD駆動可(深根)",
    "CompoundDrop" : "(LE_normal − LE_compound) / LE_normal\n→ 低い = 複合干ばつ耐性",
}


# ================================================================
# 入出力
# ================================================================

def load_inputs():
    if not PARQUET.exists():
        raise FileNotFoundError(
            f"{PARQUET} がありません。\n"
            f"v4 の末尾に下記2行を追加して再実行してください:\n"
            f"  daily.to_parquet('{PARQUET}')\n"
            f"  json.dump(closure_slopes, open('{SLOPES_JSON}','w'))\n"
        )
    if not SLOPES_JSON.exists():
        raise FileNotFoundError(f"{SLOPES_JSON} がありません(同上)")

    daily = pd.read_parquet(PARQUET)
    slopes = json.loads(SLOPES_JSON.read_text())
    print(f"[input] daily rows = {len(daily)}")
    print(f"[input] closure slopes = {slopes}")

    required = {"site","date","is_growing","drought_type","LE_corr"}
    missing = required - set(daily.columns)
    if missing:
        raise ValueError(f"daily に列が足りません: {missing}")

    return daily, slopes


# ================================================================
# 指標計算
# ================================================================

def safe_ratio(numer, denom):
    if denom is None or np.isnan(denom) or abs(denom) < DENOM_FLOOR:
        return np.nan
    r = numer / denom
    if abs(r) > RATIO_CLIP:
        return np.nan
    return r


def compute_indices_from_medians(meds):
    n = meds.get("normal");   s = meds.get("soil dry")
    a = meds.get("atm dry");  c = meds.get("compound")
    return {
        "SDS"          : safe_ratio(n - s, n) if not (np.isnan(n) or np.isnan(s)) else np.nan,
        "VAS"          : safe_ratio(a - n, n) if not (np.isnan(a) or np.isnan(n)) else np.nan,
        "DSO"          : safe_ratio(c - s, s) if not (np.isnan(c) or np.isnan(s)) else np.nan,
        "CompoundDrop" : safe_ratio(n - c, n) if not (np.isnan(n) or np.isnan(c)) else np.nan,
    }


def class_arrays(sub_df, var):
    return {c: sub_df.loc[sub_df["drought_type"]==c, var].dropna().values
            for c in CLASSES}


def bootstrap_indices(data, n_boot=N_BOOT, seed=42):
    rng = np.random.default_rng(seed)

    point_meds = {c: (np.median(v) if len(v) >= MIN_N_PER_CLASS else np.nan)
                   for c, v in data.items()}
    point = compute_indices_from_medians(point_meds)

    boots = {k: [] for k in point.keys()}
    for _ in range(n_boot):
        b_meds = {}
        valid = True
        for c, v in data.items():
            if len(v) < MIN_N_PER_CLASS:
                b_meds[c] = np.nan
                continue
            sample = rng.choice(v, size=len(v), replace=True)
            b_meds[c] = np.median(sample)
        idx = compute_indices_from_medians(b_meds)
        for k, val in idx.items():
            if not np.isnan(val):
                boots[k].append(val)

    out = {}
    for k, vals in boots.items():
        if len(vals) < N_BOOT * 0.1:
            out[k] = (point[k], np.nan, np.nan, len(vals))
        else:
            lo, hi = np.percentile(vals, CI_PCT)
            out[k] = (point[k], lo, hi, len(vals))
    return out, point_meds


# ================================================================
# pairwise 検定
# ================================================================

def pairwise_tests(data, var):
    pairs = [("normal","soil dry"),("normal","atm dry"),("normal","compound"),
             ("soil dry","compound"),("soil dry","atm dry"),("atm dry","compound")]
    rows = []
    for a, b in pairs:
        da, db = data[a], data[b]
        if len(da) < 3 or len(db) < 3:
            rows.append(dict(class_a=a, class_b=b, var=var,
                             n_a=len(da), n_b=len(db),
                             med_a=np.nan if not len(da) else float(np.median(da)),
                             med_b=np.nan if not len(db) else float(np.median(db)),
                             U=np.nan, p=np.nan, rank_biserial=np.nan))
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
# 可視化
# ================================================================

def plot_indices(summary, var, save_dir):
    """指標のバーチャート(95%CI付き)"""
    sub = summary[summary["var"] == var].copy()
    if sub.empty:
        return
    indices = ["SDS","VAS","DSO","CompoundDrop"]
    fig, axes = plt.subplots(1, 4, figsize=(22, 6))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle(f"v11 — Drought response indices on {var}",
                 fontsize=13, fontweight="bold")

    for ax, idx_name in zip(axes, indices):
        d = sub[sub["index"] == idx_name].copy()
        if d.empty:
            ax.set_title(f"{idx_name}\n(no data)"); continue
        d["key"] = d["site"] + "\n" + d["season"]
        d = d.sort_values(["site","season"], ascending=[True, False])
        x = np.arange(len(d))
        colors = [SITE_COL[s] for s in d["site"]]
        alphas = [0.95 if s == "growing" else 0.45 for s in d["season"]]
        for i, (val, lo, hi, col, alp) in enumerate(zip(
                d["point"], d["ci_lo"], d["ci_hi"], colors, alphas)):
            ax.bar(i, val, color=col, alpha=alp, edgecolor="black", lw=1.0)
            if not np.isnan(lo) and not np.isnan(hi):
                ax.errorbar(i, val, yerr=[[val-lo],[hi-val]],
                             color="black", capsize=8, lw=1.4)
            if not np.isnan(val):
                ax.text(i, val + (0.02 if val >= 0 else -0.04),
                        f"{val:+.2f}",
                        ha="center", va="bottom" if val >= 0 else "top",
                        fontsize=9, fontweight="bold")
        ax.axhline(0, color="gray", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(d["key"], fontsize=9)
        ax.set_ylabel(idx_name)
        ax.set_title(f"{idx_name}\n{INDEX_INTERP[idx_name]}", fontsize=10)
        ax.grid(axis="y", alpha=0.25)

    plt.tight_layout()
    fp = save_dir / f"v11_indices_{var}.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [save] {fp}")


def plot_distributions(daily, var, save_dir):
    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle(f"v11 — {var} distributions by drought class",
                 fontsize=13, fontweight="bold")

    for r, site in enumerate(["Oran","Tarazona"]):
        for c, is_g in enumerate([True, False]):
            ax = axes[r, c]
            sub = daily[(daily["site"] == site) & (daily["is_growing"] == is_g)]
            data = [sub.loc[sub["drought_type"]==cls, var].dropna().values
                    for cls in CLASSES]
            ns = [len(d) for d in data]
            if all(n == 0 for n in ns):
                ax.set_title(f"[{site}] {var} — no data"); continue
            bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                             medianprops=dict(color="white", lw=2.5))
            for p, cls in zip(bp["boxes"], CLASSES):
                p.set_facecolor(COLORS[cls]); p.set_alpha(0.78)
            ax.set_xticklabels([f"{c}\n(n={n})" for c, n in zip(CLASSES, ns)],
                               fontsize=9)
            season = "Growing" if is_g else "Non-growing"
            ax.set_title(f"[{site}] {var} — {season}", fontweight="bold")
            ax.set_ylabel(var); ax.grid(axis="y", alpha=0.25)

    plt.tight_layout()
    fp = save_dir / f"v11_distributions_{var}.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [save] {fp}")


def plot_summary_panel(summary, save_dir):
    """SDS/DSO を Oran vs Tarazona で対比する1枚サマリ"""
    sub = summary[(summary["var"] == "LE_corr") &
                   (summary["index"].isin(["SDS","DSO"])) &
                   (summary["season"] == "growing")].copy()
    if sub.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle("★ v11 Summary — Growing season, LE_corr based",
                 fontsize=13, fontweight="bold")

    for ax, idx_name in zip(axes, ["SDS","DSO"]):
        d = sub[sub["index"] == idx_name].sort_values("site")
        if d.empty:
            ax.set_title(f"{idx_name}: no data"); continue
        x = np.arange(len(d))
        colors = [SITE_COL[s] for s in d["site"]]
        for i, (val, lo, hi, col) in enumerate(zip(
                d["point"], d["ci_lo"], d["ci_hi"], colors)):
            ax.bar(i, val, color=col, alpha=0.85, edgecolor="black", lw=1.2)
            if not np.isnan(lo):
                ax.errorbar(i, val, yerr=[[val-lo],[hi-val]],
                             color="black", capsize=10, lw=1.5)
            ax.text(i, val + (0.02 if val >= 0 else -0.04),
                    f"{val:+.3f}\n[{lo:+.2f}, {hi:+.2f}]",
                    ha="center", va="bottom" if val >= 0 else "top",
                    fontsize=10, fontweight="bold")
        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.set_xticks(x); ax.set_xticklabels(d["site"], fontsize=11, fontweight="bold")
        ax.set_ylabel(idx_name, fontsize=11)
        ax.set_title(f"{idx_name}\n{INDEX_INTERP[idx_name]}", fontsize=10)
        ax.grid(axis="y", alpha=0.25)

    plt.tight_layout()
    fp = save_dir / "v11_summary_LE.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# MAIN
# ================================================================

def main():
    print("=" * 60)
    print("解析A v11: SDS/VAS/DSO 指標 — 4象限分類定量化")
    print("=" * 60)

    daily, slopes = load_inputs()

    summary_rows, test_rows, count_rows = [], [], []

    for site in sorted(daily["site"].unique()):
        for is_g in [True, False]:
            season = "growing" if is_g else "non_growing"
            sub = daily[(daily["site"] == site) & (daily["is_growing"] == is_g)]
            counts = {c: int((sub["drought_type"] == c).sum()) for c in CLASSES}
            count_rows.append(dict(site=site, season=season, **counts))

            print(f"\n[{site} / {season}]  counts: {counts}")

            for var in VARS:
                if var not in sub.columns:
                    continue
                data = class_arrays(sub, var)
                idx_results, meds = bootstrap_indices(data)
                med_str = "  ".join(
                    f"{c}={meds[c]:.2f}" if not np.isnan(meds[c]) else f"{c}=NaN"
                    for c in CLASSES)
                print(f"  {var}: medians {med_str}")
                for idx_name, (pt, lo, hi, n_b) in idx_results.items():
                    summary_rows.append(dict(
                        site=site, season=season, var=var, index=idx_name,
                        point=pt, ci_lo=lo, ci_hi=hi, n_boot=n_b,
                        med_normal=meds["normal"], med_soil_dry=meds["soil dry"],
                        med_atm_dry=meds["atm dry"], med_compound=meds["compound"],
                    ))
                    if not np.isnan(pt):
                        ci_str = (f"[{lo:+.3f}, {hi:+.3f}]"
                                   if not np.isnan(lo) else "[CI算出不可]")
                        print(f"    {idx_name:13s}: {pt:+.3f}  {ci_str}  (n_boot={n_b})")

                for t in pairwise_tests(data, var):
                    t.update(dict(site=site, season=season))
                    test_rows.append(t)

    summary_df = pd.DataFrame(summary_rows)
    tests_df   = pd.DataFrame(test_rows)
    counts_df  = pd.DataFrame(count_rows)

    summary_df.to_csv(OUT_DIR / "v11_indices.csv", index=False)
    tests_df.to_csv(OUT_DIR / "v11_pairwise_tests.csv", index=False)
    counts_df.to_csv(OUT_DIR / "v11_class_counts.csv", index=False)
    print(f"\n[save] {OUT_DIR}/v11_indices.csv")
    print(f"[save] {OUT_DIR}/v11_pairwise_tests.csv")
    print(f"[save] {OUT_DIR}/v11_class_counts.csv")

    print(f"\n--- 可視化 ---")
    for var in VARS:
        if var in daily.columns:
            plot_indices(summary_df, var, OUT_DIR)
            plot_distributions(daily, var, OUT_DIR)
    plot_summary_panel(summary_df, OUT_DIR)

    # 最終解釈
    print("\n" + "=" * 60)
    print("★ 最終解釈 (Growing season, LE_corr ベース)")
    print("=" * 60)
    g = summary_df[(summary_df["var"] == "LE_corr") &
                    (summary_df["season"] == "growing")]
    for idx_name in ["SDS","VAS","DSO","CompoundDrop"]:
        d = g[g["index"] == idx_name]
        if d.empty: continue
        print(f"\n[{idx_name}]  {INDEX_INTERP[idx_name].splitlines()[0]}")
        for _, row in d.iterrows():
            ci = (f"[{row['ci_lo']:+.3f}, {row['ci_hi']:+.3f}]"
                   if not np.isnan(row["ci_lo"]) else "[CI算出不可]")
            print(f"  {row['site']:9s}: {row['point']:+.3f}  {ci}")

    sds_o = g[(g["site"]=="Oran") & (g["index"]=="SDS")]["point"]
    sds_t = g[(g["site"]=="Tarazona") & (g["index"]=="SDS")]["point"]
    if len(sds_o) and len(sds_t) and not (np.isnan(sds_o.iloc[0]) or np.isnan(sds_t.iloc[0])):
        print(f"\n→ 仮説検証(SDS):")
        print(f"   Oran SDS = {sds_o.iloc[0]:+.3f}  (高いほど表層SWC依存=浅根)")
        print(f"   Tarazona SDS = {sds_t.iloc[0]:+.3f}")
        if sds_t.iloc[0] < sds_o.iloc[0] - 0.1:
            print(f"   ⇒ Tarazona の SDS が顕著に低い = 表層SWC低下に対しLE維持力が強い")
            print(f"     = 深根仮説を支持")
        else:
            print(f"   ⇒ 差が小さい/逆転 → 仮説を再検討")

    print(f"\n[done] outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
