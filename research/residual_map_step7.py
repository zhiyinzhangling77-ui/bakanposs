"""旗7：残差の構造を"場所ごと"に地図化し、未観測原因のホットスポットを特定する。

問2a の手順骨格の次段：
  既知変数で説明 → 残差を場所・季節ごとに地図化 → 偏りのある場所を特定
  → その場所に共通する未観測要因を仮説化 → プロキシで検証。

ここでは複数の"場所(region)"を作り、隠れ駆動 Z の強さを場所ごとに変える（中央にホットスポット）。
既知の X で Y をモデル化した残差の構造（自己相関）を場所ごとに測り、
「どこに未観測原因の足跡が強いか」を地図（バー）で当てられるかを確認する。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from residual_footprint_step5 import _linfit_resid, _autocorr1
from unobserved_common_cause_step4 import mutual_info


def make_region(n, hidden_strength, seed):
    """1 つの場所のデータ。X は観測できる駆動、Z は未観測の隠れ駆動（強さ可変）。"""
    rng = np.random.default_rng(seed)
    z = np.zeros(n); x = np.zeros(n); y = np.zeros(n)
    ez, ex, ey = (rng.normal(0, 1, n) for _ in range(3))
    for t in range(1, n):
        z[t] = 0.7 * z[t - 1] + ez[t]                    # 自己相関のある隠れ駆動
        x[t] = 0.7 * x[t - 1] + ex[t]                    # 観測できる駆動
        y[t] = 0.8 * x[t - 1] + hidden_strength * z[t - 1] + ey[t]
    return x, y, z


def draw_map(strengths, ac_map, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jp_path = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jp_path) if Path(jp_path).exists() else None
    r = np.arange(len(ac_map))
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.bar(r, ac_map, color="#2e8b57")
    ax.plot(r, strengths * max(ac_map), "o--", color="#c0392b",
            label="仕込んだ隠れ駆動の強さ（正解・スケール調整）")
    ax.set_xlabel("場所（region）", fontproperties=jp)
    ax.set_ylabel("残差の構造（自己相関）＝未観測原因の足跡", fontproperties=jp)
    ax.set_title("残差の構造マップ：ホットスポット＝未観測原因が強い場所",
                 fontproperties=jp)
    ax.legend(prop=jp, frameon=False)
    fig.savefig(path, bbox_inches="tight", dpi=130); plt.close(fig)


if __name__ == "__main__":
    # 中央にホットスポットを仕込む（両端は隠れ駆動なし）
    strengths = np.array([0.0, 0.4, 1.0, 1.0, 0.4, 0.0])
    n = 4000

    print("=== 残差の構造マップ（未観測原因のホットスポット特定）===\n")
    print(f"  {'場所':>4}  {'隠れ駆動(正解)':>12}  {'残差の自己相関':>12}  {'I(残差;Z)':>10}")
    ac_map = []
    for r, h in enumerate(strengths):
        x, y, z = make_region(n, h, seed=r)
        resid = _linfit_resid(y[1:], x[:-1])          # 既知 X で Y をモデル化→残差
        ac = _autocorr1(resid)
        mi = mutual_info(resid[1:], z[1:-1], 8)        # 種明かし：残差 と 隠れ Z
        ac_map.append(ac)
        print(f"  {r:>4}  {h:>12.1f}  {ac:>12.3f}  {mi:>10.4f}")

    ac_map = np.array(ac_map)
    hotspot = int(np.argmax(ac_map))
    true_hot = int(np.argmax(strengths))
    print(f"\n  残差の構造が最大の場所 = {hotspot}（仕込んだ最大 = {true_hot} 付近）")
    print("  判定: " + ("✅ ホットスポットを当てられた（残差マップ＝未観測原因の在り処）"
                        if abs(hotspot - 2.5) <= 1.5 and ac_map[hotspot] > ac_map[0] + 0.05
                        else "⚠ うまく出ない（強さ/ノイズ/N を調整）"))

    out = Path(__file__).resolve().parent / "residual_map_step7.png"
    draw_map(strengths, ac_map, out)
    print(f"  [図] {out}")
    print("\n  → 実データでは region=サイト/グリッド, season=月。残差の構造が偏る場所を特定し、")
    print("     そこに共通する未観測要因を生態学の知識で仮説化→プロキシ（衛星等）で検証（問2a）。")
