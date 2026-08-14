"""旗11：不変性を実データで。年（or サイト）を"環境"にし、各予測子→目的変数の関係が
環境をまたいで不変か（＝転移する因果か）を係数の安定性でスコア化する。

指針の外挿の道：多環境で"不変な関係"ほど因果で、新しい条件に転移しやすい（Peters 2016）。
compare_sites（リンクの出現頻度）を「関係の強さの不変性(CV)」の観点で精緻化する。

各候補 V について、環境ごとに Y_t ~ V_t + Y_{t-1..}（自己ラグで Y の記憶を除く）を回帰し、
V の係数 b を得る。環境をまたいだ b のばらつき CV=std/|mean| が小さいほど不変＝転移候補。

    # 年を環境に（JP-Tak の各年で、GEP への各予測子の不変性）
    python research/invariance_real_step11.py --site JP-Tak \
        --years 2003 2004 2005 2006 2007 2008 --month 7 8 --target GEP --self-lag 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def slope_of(Y: np.ndarray, V: np.ndarray, self_lag: int) -> float:
    """Y_t ~ V_t + Y_{t-1..t-L} を回帰し、V の係数（contemporaneous な効果）を返す。"""
    L = self_lag
    Yt = Y[L:]
    cols = [V[L:]]
    for k in range(1, L + 1):
        cols.append(Y[L - k:len(Y) - k])
    A = np.column_stack(cols + [np.ones(len(Yt))])
    coef, *_ = np.linalg.lstsq(A, Yt, rcond=None)
    return float(coef[0])


def main() -> None:
    from japanflux_pn.config import RK_VARS
    from japanflux_pn.preprocess import load_corevars_hh

    p = argparse.ArgumentParser(description="実データの不変性スコア（年を環境に）")
    p.add_argument("--site", required=True)
    p.add_argument("--years", type=int, nargs="+", required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--target", required=True)
    p.add_argument("--self-lag", type=int, default=1)
    p.add_argument("--include-carbon", action="store_true",
                   help="GEP/GER/NEE も候補に含める（既定は除外：分割の代数的相関を避ける）")
    a = p.parse_args()

    CARBON = {"GEP", "GER", "NEE"}
    drop = set() if a.include_carbon else (CARBON - {a.target})
    candidates = [v for v in RK_VARS if v != a.target and v not in drop]

    # 環境（年）ごとに、各候補の係数を集める
    slopes: dict[str, list[float]] = {v: [] for v in candidates}
    used_years = []
    for y in a.years:
        try:
            pre = load_corevars_hh(a.site, y, a.month, None)
            vf = pre.valid_frame
            Y = vf[a.target].to_numpy(dtype=float)
            for v in candidates:
                slopes[v].append(slope_of(Y, vf[v].to_numpy(dtype=float), a.self_lag))
            used_years.append(y)
        except Exception as e:
            print(f"  {y}: SKIP {type(e).__name__}: {e}")

    print(f"\n=== 不変性スコア: {a.target} への各予測子 / {a.site}"
          f"（環境={len(used_years)}年）===")
    print("  b=関係の強さ（環境平均）, CV=環境をまたいだ b のばらつき（小=不変=転移候補）\n")
    print(f"  {'予測子':>6}  {'平均b':>8}  {'|b|':>7}  {'CV(不変性)':>10}  判定")
    rows = []
    for v in candidates:
        s = np.array(slopes[v])
        mb = float(s.mean()); ab = abs(mb)
        cv = float(s.std() / ab) if ab > 1e-9 else np.inf
        rows.append((v, mb, ab, cv))
    # 効果が弱すぎ(|b|小)ると CV は無意味なので、|b| で足切りしてから不変性で並べる
    strong = [r for r in rows if r[2] >= 0.05]
    weak = [r for r in rows if r[2] < 0.05]
    for v, mb, ab, cv in sorted(strong, key=lambda r: r[3]):
        mark = "✅不変=転移候補" if cv < 0.5 else ("△やや変動" if cv < 1.0 else "×環境で変わる")
        print(f"  {v:>6}  {mb:>8.3f}  {ab:>7.3f}  {cv:>10.3f}  {mark}")
    for v, mb, ab, cv in sorted(weak, key=lambda r: -r[2]):
        print(f"  {v:>6}  {mb:>8.3f}  {ab:>7.3f}  {'  --':>10}  （効果が弱く不変性は評価外）")

    print("\n  → 不変な（CV小）強い関係ほど、時間（年）をまたいで転移する因果の候補。")
    print("     ※これは1サイトの『時間』不変性。真の転移にはサイト間の環境も要る（--site を変えて比較）。")
    print("     さらに機構（物理・生物法則）と整合すれば主張が最強（不変性×機構×比較設計）。")


if __name__ == "__main__":
    main()
