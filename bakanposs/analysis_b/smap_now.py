"""
analysis_B_smap_now.py
======================
H-SAF アカウント承認待ちの間に、
master_full_v2.csv に既存の SMAP SM で driver attribution を先行実行する。

SMAP は H-SAF H141/H142/H26 と同じ役割（衛星 SM）を担う。
H-SAF が揃ったら analysis_B_v6_driver_attribution.py に置き換える。

Usage
-----
python analysis_B_smap_now.py \
    --master data/master_full_v2.csv \
    --output figures/driver_attribution_smap.pdf

列名 (master_full_v2.csv から確認済み):
  ET_mm          : EC フラックスタワー ET
  ET_metv3_mm    : ETv3 衛星 ET
  VPD_kPa_mean   : 日平均 VPD
  Rn_Wm2         : 正味放射
  SWC            : in-situ 土壌水分
  smap_surface   : SMAP 表層 SM
  smap_rootzone  : SMAP 根圏 SM
"""

import argparse
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler

# -----------------------------------------------------------------------
# 設定 (master_full_v2.csv の実際の列名)
# -----------------------------------------------------------------------
SITE        = "TzM"
JJA_MONTHS  = [6, 7, 8]

EC_ET_COL   = "ET_mm"
SAT_ET_COL  = "ET_metv3_mm"
VPD_COL     = "VPD_kPa_mean"
RN_COL      = "Rn_Wm2"
SWC_COL     = "SWC"
SMAP_S_COL  = "smap_surface"
SMAP_R_COL  = "smap_rootzone"

X_COLS = [VPD_COL, RN_COL, SWC_COL, SMAP_S_COL, SMAP_R_COL]

VAR_LABELS = {
    VPD_COL:    "VPD",
    RN_COL:     "Rn",
    SWC_COL:    "Tower SWC",
    SMAP_S_COL: "SMAP Surface",
    SMAP_R_COL: "SMAP Rootzone",
}

COLORS = {
    "EC ET":  "#1B5E20",
    "Sat ET": "#0D47A1",
}


# -----------------------------------------------------------------------
# ユーティリティ
# -----------------------------------------------------------------------
def weekly_agg(df: pd.DataFrame, cols: list) -> pd.DataFrame:
    df = df.copy()
    df["week"] = df["date"].dt.to_period("W").apply(lambda p: p.start_time)
    return df.groupby("week")[cols].mean().reset_index().rename(columns={"week": "date"})


def ols_stats(X: pd.DataFrame, y: pd.Series, label: str) -> dict:
    X_c = sm.add_constant(X, has_constant="add")
    mod = sm.OLS(y, X_c, missing="drop").fit()
    coefs = mod.params.drop("const", errors="ignore")
    pvals = mod.pvalues.drop("const", errors="ignore")
    cis   = mod.conf_int().drop("const", errors="ignore")
    return {"label": label, "r2": mod.rsquared, "r2_adj": mod.rsquared_adj,
            "n": int(mod.nobs), "coefs": coefs, "pvals": pvals,
            "ci_low": cis[0], "ci_high": cis[1], "model": mod}


def compute_vif(X: pd.DataFrame) -> pd.DataFrame:
    rows = [{"variable": c,
             "VIF": variance_inflation_factor(X.values, i)}
            for i, c in enumerate(X.columns)]
    return pd.DataFrame(rows).set_index("variable")


# -----------------------------------------------------------------------
# メイン
# -----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--master",  required=True)
    parser.add_argument("--output",  default="figures/driver_attribution_smap.pdf")
    parser.add_argument("--site",    default=SITE)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    base = args.output.replace(".pdf", "")

    # --- 読み込み ---
    print(f"Loading: {args.master}")
    df = pd.read_csv(args.master, parse_dates=["date"])
    df_site = df[df["site"] == args.site].copy()
    print(f"  {args.site}: {len(df_site)} rows total")

    # 列存在チェック
    missing = [c for c in [EC_ET_COL, SAT_ET_COL] + X_COLS if c not in df_site.columns]
    if missing:
        print(f"  [ERROR] 以下の列が見つかりません: {missing}")
        return

    # 夏季フィルタ
    df_jja = df_site[df_site["date"].dt.month.isin(JJA_MONTHS)].copy()
    print(f"  JJA rows: {len(df_jja)}")

    # 週次集計
    all_cols = [EC_ET_COL, SAT_ET_COL] + X_COLS
    df_w = weekly_agg(df_jja, all_cols)
    print(f"  Weekly rows: {len(df_w)}")

    # dropna
    valid = df_w[X_COLS + [EC_ET_COL, SAT_ET_COL]].dropna()
    print(f"  Valid rows (all vars non-NaN): {len(valid)}")
    if len(valid) < 10:
        print("  [WARN] 有効行が少なすぎます。SMAP データの欠損を確認してください。")
        print(df_w[SMAP_S_COL].describe())
        return

    # 標準化
    sc_X   = StandardScaler()
    sc_ec  = StandardScaler()
    sc_sat = StandardScaler()
    X_z      = pd.DataFrame(sc_X.fit_transform(valid[X_COLS]),
                            columns=X_COLS, index=valid.index)
    y_ec_z   = pd.Series(sc_ec.fit_transform(valid[[EC_ET_COL]]).ravel(),  index=valid.index)
    y_sat_z  = pd.Series(sc_sat.fit_transform(valid[[SAT_ET_COL]]).ravel(), index=valid.index)

    # VIF
    vif_df = compute_vif(X_z)
    print("\n--- VIF ---")
    print(vif_df.round(2).to_string())
    vif_df.to_csv(f"{base}_vif.csv")

    # OLS
    res_ec  = ols_stats(X_z, y_ec_z,  "EC ET")
    res_sat = ols_stats(X_z, y_sat_z, "Sat ET (ETv3)")

    print(f"\n--- EC ET  R²={res_ec['r2']:.3f}  adj-R²={res_ec['r2_adj']:.3f}  n={res_ec['n']} ---")
    print(res_ec["model"].summary().tables[1])
    print(f"\n--- Sat ET R²={res_sat['r2']:.3f}  adj-R²={res_sat['r2_adj']:.3f}  n={res_sat['n']} ---")
    print(res_sat["model"].summary().tables[1])

    # 係数テーブル保存
    rows = []
    for res in [res_ec, res_sat]:
        for var in X_COLS:
            if var in res["coefs"].index:
                rows.append({"target": res["label"], "variable": var,
                             "coef":    res["coefs"][var],
                             "pval":    res["pvals"][var],
                             "ci_low":  res["ci_low"][var],
                             "ci_high": res["ci_high"][var]})
    coef_df = pd.DataFrame(rows)
    coef_df.to_csv(f"{base}_coef.csv", index=False)
    print(f"\n[Saved] {base}_coef.csv")

    # -----------------------------------------------------------------------
    # Figure 1: 標準化係数バープロット (EC vs Sat 横並び)
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    x = np.arange(len(X_COLS))

    for ax, res, color in zip(axes,
                               [res_ec, res_sat],
                               [COLORS["EC ET"], COLORS["Sat ET"]]):
        coefs   = [res["coefs"].get(v, np.nan) for v in X_COLS]
        pvals   = [res["pvals"].get(v, np.nan) for v in X_COLS]
        ci_low  = [res["ci_low"].get(v, np.nan)  for v in X_COLS]
        ci_high = [res["ci_high"].get(v, np.nan) for v in X_COLS]
        yerr_lo = np.array(coefs) - np.array(ci_low)
        yerr_hi = np.array(ci_high) - np.array(coefs)

        ax.bar(x, coefs, color=color, alpha=0.75, width=0.55,
               yerr=[np.nan_to_num(yerr_lo), np.nan_to_num(yerr_hi)],
               error_kw={"capsize": 4, "lw": 1.5})

        for i, (c, p) in enumerate(zip(coefs, pvals)):
            if np.isnan(c) or np.isnan(p):
                continue
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            offset = 0.04 + abs(np.nan_to_num(yerr_hi[i]))
            ax.text(i, c + (offset if c >= 0 else -offset - 0.08),
                    sig, ha="center", fontsize=11, color=color)

        # SMAP グループをハイライト
        smap_idxs = [i for i, v in enumerate(X_COLS) if "smap" in v]
        if smap_idxs:
            ax.axvspan(min(smap_idxs) - 0.4, max(smap_idxs) + 0.4,
                       color="#FF9800", alpha=0.08, zorder=0, label="Satellite SM (SMAP)")

        ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels([VAR_LABELS[v] for v in X_COLS],
                           rotation=20, ha="right", fontsize=10)
        ax.set_ylabel("Standardized coefficient (β)", fontsize=10)
        ax.set_title(
            f"{res['label']}\nR²={res['r2']:.3f}  adj-R²={res['r2_adj']:.3f}  n={res['n']}",
            fontsize=11, fontweight="bold", color=color)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        f"Driver Attribution — {args.site} JJA weekly\n"
        f"(SMAP SM 使用: H-SAF H141/H142/H26 の暫定代替)",
        fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(args.output, bbox_inches="tight", dpi=200)
    print(f"[Saved] {args.output}")

    # -----------------------------------------------------------------------
    # Figure 2: EC vs Sat の係数を直接比較するミラープロット
    # -----------------------------------------------------------------------
    fig2, ax2 = plt.subplots(figsize=(9, 5))
    width = 0.35

    coefs_ec  = [res_ec["coefs"].get(v, np.nan)  for v in X_COLS]
    coefs_sat = [res_sat["coefs"].get(v, np.nan) for v in X_COLS]
    pvals_ec  = [res_ec["pvals"].get(v, np.nan)  for v in X_COLS]
    pvals_sat = [res_sat["pvals"].get(v, np.nan) for v in X_COLS]

    bars_ec  = ax2.bar(x - width/2, coefs_ec,  width, color=COLORS["EC ET"],
                       alpha=0.75, label=f"EC ET (R²={res_ec['r2']:.3f})")
    bars_sat = ax2.bar(x + width/2, coefs_sat, width, color=COLORS["Sat ET"],
                       alpha=0.75, label=f"Sat ET/ETv3 (R²={res_sat['r2']:.3f})")

    for i, (ce, cs, pe, ps) in enumerate(zip(coefs_ec, coefs_sat, pvals_ec, pvals_sat)):
        for c, p, offset in [(ce, pe, -width/2), (cs, ps, +width/2)]:
            if np.isnan(c) or np.isnan(p):
                continue
            sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
            if sig:
                ax2.text(i + offset, c + (0.03 if c >= 0 else -0.08),
                         sig, ha="center", fontsize=9)

    # SMAP グループ背景
    smap_idxs = [i for i, v in enumerate(X_COLS) if "smap" in v]
    if smap_idxs:
        ax2.axvspan(min(smap_idxs)-0.5, max(smap_idxs)+0.5,
                    color="#FF9800", alpha=0.07, zorder=0)
        ax2.text(np.mean(smap_idxs), ax2.get_ylim()[0] if ax2.get_ylim()[0] < 0 else -0.05,
                 "← Satellite SM →", ha="center", fontsize=8, color="#FF9800", style="italic")

    ax2.axhline(0, color="black", lw=0.8, ls="--", alpha=0.6)
    ax2.set_xticks(x)
    ax2.set_xticklabels([VAR_LABELS[v] for v in X_COLS], fontsize=11)
    ax2.set_ylabel("Standardized β", fontsize=11)
    ax2.set_title(
        f"EC ET vs ETv3 — Driver Attribution\n"
        f"{args.site} JJA weekly  (SMAP SM, 暫定)",
        fontsize=12, fontweight="bold")
    ax2.legend(fontsize=10, framealpha=0.85)
    ax2.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()

    mirror_out = args.output.replace(".pdf", "_mirror.pdf")
    fig2.savefig(mirror_out, bbox_inches="tight", dpi=200)
    print(f"[Saved] {mirror_out}")

    # -----------------------------------------------------------------------
    # 期待パターンとの一致確認
    # -----------------------------------------------------------------------
    print("\n--- 仮説との一致チェック ---")
    for var, label in VAR_LABELS.items():
        ce = res_ec["coefs"].get(var, np.nan)
        cs = res_sat["coefs"].get(var, np.nan)
        pe = res_ec["pvals"].get(var, np.nan)
        ps = res_sat["pvals"].get(var, np.nan)
        if "smap" in var:
            match_ec  = "✓" if (not np.isnan(ce) and abs(ce) < 0.1) else "✗"
            match_sat = "✓" if (not np.isnan(cs) and abs(cs) > 0.2) else "✗"
            print(f"  {label:20s}: EC β={ce:+.3f}(p={pe:.3f}) {match_ec}≈0  "
                  f"Sat β={cs:+.3f}(p={ps:.3f}) {match_sat}>0.2")
        elif var == SWC_COL:
            match_ec  = "✓" if (not np.isnan(ce) and ce > 0.15) else "✗"
            match_sat = "✓" if (not np.isnan(cs) and abs(cs) < 0.1) else "✗"
            print(f"  {label:20s}: EC β={ce:+.3f}(p={pe:.3f}) {match_ec}>0.15  "
                  f"Sat β={cs:+.3f}(p={ps:.3f}) {match_sat}≈0")


if __name__ == "__main__":
    main()
