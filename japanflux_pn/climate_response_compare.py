"""多サイト横断：情報結合の「状態依存」が生態系をまたいで不変か（＝転移するか）を見る。

climate_response が 1 サイトで出す「結合強度(アノマリ MI) vs その年の状態(生の VPD/土壌水分)」の
年々 Spearman 相関を、複数の生態系で並べる。同じ符号で有意なら**状態依存の法則が普遍**（転移する）、
サイトで符号が飛ぶなら**その生態系固有**。EXTRAPOLATION_SYNTHESIS の「強さは状態依存で転移しない」を、
1 サイト（JP-Tak: I(Rg;GEP) vs VPD r=-0.68）から**多生態系**へ広げて検証する薄いラッパ。

目玉の 2 本（乾湿の気候軸で効くはず）：
  - I(Rg;GEP) vs VPD  … 乾いた年ほど放射↔光合成が脱結合（負を期待）
  - I(th;GER) vs VPD  … 乾いた年ほど土壌水分↔呼吸が律速に（正を期待）

    # 既定の生態系多様サイト（湿潤日本↔乾燥モンゴル/青海の気候勾配）
    python -m japanflux_pn.climate_response_compare

    # サイトを明示 / 図を保存
    python -m japanflux_pn.climate_response_compare \
        --sites JP-Tak JP-Ta2 JP-BBY CN-HaM MN-Kbu --heatmap state_dep_cross.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .climate_response import scan, _spearman, _spearman_p, COUPLINGS


# 既定：水田を含まない生態系多様サイト（湿潤↔乾燥の VPD 勾配）
DEFAULT_SITES = ["JP-Tak", "JP-Ta2", "JP-BBY", "CN-HaM", "MN-Kbu"]

# 状態依存として見る (結合, 状態) の組。第 3 要素は期待する符号（物理仮説）。
PROBES = [
    ("I_Rg_GEP", "VPD_mean", "-"),   # 乾→放射↔光合成 脱結合
    ("I_Rg_GEP", "SWC_mean", "+"),   # 湿→放射↔光合成 結合
    ("I_th_GER", "VPD_mean", "+"),   # 乾→土壌水分↔呼吸 が律速
    ("I_Rg_gLE", "VPD_mean", "-"),   # 乾→放射↔蒸散 脱結合（気孔閉鎖）
]
STATE_LABEL = {"VPD_mean": "VPD", "SWC_mean": "土壌水分", "Ta_mean": "気温"}
COUP_LABEL = {f"I_{a}_{b}": f"{a}→{b}" for a, b in COUPLINGS}


def probe_label(coup: str, state: str) -> str:
    return f"{COUP_LABEL.get(coup, coup)} vs {STATE_LABEL.get(state, state)}"


def collect(sites: list[str], config: AnalysisConfig | None = None
            ) -> dict[str, pd.DataFrame]:
    """各サイトの年々メトリクスを集める（欠測/年数不足のサイトは飛ばす）。"""
    config = config or AnalysisConfig()
    per_site: dict[str, pd.DataFrame] = {}
    for s in sites:
        print(f"\n----- {s} 収集中 -----", flush=True)
        try:
            df = scan(s, config)
        except Exception as e:  # データ欠如・登録漏れは飛ばして続行
            print(f"  {s}: SKIP {type(e).__name__}: {e}", flush=True)
            continue
        if len(df) < 4:
            print(f"  {s}: 有効年 {len(df)} < 4 → 相関評価に不足、除外", flush=True)
            continue
        per_site[s] = df
    return per_site


def build_r_table(per_site: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """行=状態依存プローブ, 列=サイト の Spearman r 表。"""
    rows: dict[str, dict[str, float]] = {}
    for coup, state, _ in PROBES:
        label = probe_label(coup, state)
        rows[label] = {}
        for s, df in per_site.items():
            if coup in df and state in df:
                rows[label][s] = _spearman(df[coup].to_numpy(), df[state].to_numpy())
            else:
                rows[label][s] = np.nan
    return pd.DataFrame(rows).T  # 行=プローブ, 列=サイト


def draw_heatmap(rtab: pd.DataFrame, path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jp_path = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jp_path) if Path(jp_path).exists() else None
    m = rtab
    fig, ax = plt.subplots(figsize=(1.6 + 1.0 * m.shape[1], 0.7 * m.shape[0] + 1.6))
    im = ax.imshow(m.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=20, ha="right")
    ax.set_yticks(range(m.shape[0]))
    ax.set_yticklabels(m.index, fontproperties=jp)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.iat[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=9,
                        color="white" if abs(v) >= 0.55 else "#222")
    ax.set_title("状態依存の生態系間比較（結合強度 vs 気候状態の年々相関 r）",
                 fontproperties=jp)
    cb = fig.colorbar(im, ax=ax, fraction=0.03)
    cb.set_label("Spearman r（＋赤/−青）", fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=130); plt.close(fig)


def report(sites: list[str], config: AnalysisConfig | None = None,
           heatmap=None, out_csv=None, n_min: int = 10) -> pd.DataFrame:
    per_site = collect(sites, config)
    if len(per_site) < 2:
        print(f"\n有効サイト {len(per_site)} < 2。比較できません。")
        return pd.DataFrame()

    used = list(per_site.keys())
    rtab = build_r_table(per_site)

    print(f"\n===== 状態依存の生態系間比較（{len(used)} サイト: {', '.join(used)}）=====")
    print("  各セル = 結合強度(アノマリ MI) vs その年の状態(生平均) の年々 Spearman r\n")
    print("  " + " " * 22 + "".join(f"{s:>10}" for s in used))
    for coup, state, sign in PROBES:
        label = probe_label(coup, state)
        cells = "".join(f"{rtab.loc[label, s]:>+10.2f}"
                        if np.isfinite(rtab.loc[label, s]) else f"{'--':>10}"
                        for s in used)
        print(f"  期待{sign} {label:<16}{cells}")

    # 目玉 3 本を p 値つきで判定。★重要：符号の一貫性は「検出力のある年数を持つサイト」
    # だけで採る。n が小さいサイト（n<n_min）は |r| が大きくても偶然で符号が飛ぶので、
    # そこを一貫性の証拠に使うと誤った結論（見かけの生態系固有 or 見かけの普遍）を生む。
    print(f"\n=== 状態依存は転移するか（符号一貫性＋順列検定 p）===")
    print(f"  ※判定は検出力のあるサイト（有効年 n≥{n_min}）だけで採る。"
          f"n<{n_min} は年数不足で符号が偶然に飛ぶため判定保留。\n")
    print(f"  参考: p<0.05 に要る |r| ≈ n=5→0.88, n=10→0.63, n=11→0.60, n=21→0.43")
    for coup, state, sign in PROBES[:3]:
        label = probe_label(coup, state)
        print(f"\n  ■ {label}（期待符号 {sign}）")
        powered_signs = []          # (n≥n_min のサイトの) 符号一致 1/0
        for s in used:
            df = per_site[s]
            if coup not in df or state not in df:
                print(f"    {s:>8}: 列なし"); continue
            n = len(df)
            r, p = _spearman_p(df[coup].to_numpy(), df[state].to_numpy())
            agree = (sign == "-" and r < 0) or (sign == "+" and r > 0)
            hit = "○一致" if agree else "×逆"
            star = "有意" if (np.isfinite(p) and p < 0.05) else "n.s."
            if n >= n_min:
                powered_signs.append(1 if agree else 0)
                tag = "検出力OK"
            else:
                tag = "★年数不足=判定保留"
            print(f"    {s:>8}: r={r:+.2f}  p={p:.3f}  {hit} {star}  (n={n}, {tag})")
        k = len(powered_signs)
        if k == 0:
            print(f"    → 検出力のあるサイト無し（全て n<{n_min}）＝この軸は判定不能。")
        else:
            agree = sum(powered_signs)
            if agree == k and k >= 3:
                verd = "✅ 検出力サイトで符号一致＝状態依存の法則が転移する候補（要・有意性の積み増し）"
            elif agree == k:
                verd = "○ 符号一致だが検出力サイトが少ない（サイト/年を増やして確認）"
            elif agree >= max(1, k - 1):
                verd = "△ 概ね一致（1 例外）＝多くの生態系で転移の可能性、例外は要精査"
            else:
                verd = "× 検出力サイトでも符号が飛ぶ＝状態依存は生態系固有の疑い"
            print(f"    → 検出力サイト {agree}/{k} 符号一致: {verd}")

    print("\n  読み方: 判定は n≥%d のサイトだけで採る。符号が(検出力のあるサイトで)"
          "一致＝『乾いた年ほど脱結合』等が生態系を越える転移候補。" % n_min)
    print("         n が小さいサイトの符号の飛びは偶然＝『生態系固有』の証拠にしてはいけない。")
    print("         ＝EXTRAPOLATION_SYNTHESIS「強さは状態依存で転移しない」の多生態系での検証。")
    print("         注: 多くのサイトは健全な夏の年数が少ない（本質的な検出力不足）。")
    print("             真の転移判定には長期記録／多サイト（衛星プロキシで年数を稼ぐ）が要る。")

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        rtab.to_csv(out_csv)
        print(f"\n[saved] {out_csv}")
    if heatmap:
        draw_heatmap(rtab, heatmap)
        print(f"[saved] {heatmap}")
    return rtab


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="状態依存の生態系間比較（転移するか）")
    p.add_argument("--sites", nargs="+", default=DEFAULT_SITES,
                   help=f"比較するサイト（既定: {' '.join(DEFAULT_SITES)}）")
    p.add_argument("--heatmap", default=None)
    p.add_argument("--out-csv", default=None)
    p.add_argument("--n-min", type=int, default=10,
                   help="符号一貫性を採るのに要る最小の有効年数（既定10。これ未満は判定保留）")
    a = p.parse_args(argv)
    report(a.sites, heatmap=a.heatmap, out_csv=a.out_csv, n_min=a.n_min)


if __name__ == "__main__":
    main()
