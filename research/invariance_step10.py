"""旗10：不変性による因果予測（Invariant Causal Prediction の最小版）。

外挿（将来予測）を"データ側"から正当化する道＝**多環境で関係が不変な予測子は因果で、
新しい条件にも転移しやすい**（Peters, Bühlmann & Meinshausen 2016）。逆に、環境ごとに
関係が変わる予測子は見かけ（交絡や結果）で、転移しない。

最小実験：
  - 真の因果 X1→Y（係数は全環境で不変）。X2 は Y の結果（子）。
  - 環境ごとに X1 の分散を変える"介入"をする（気候レジームが違うイメージ）。
  - 予測子ごとに「残差が環境をまたいで不変か」を測る。
    X1 で説明した残差は環境不変（＝因果・転移する）／X2 では環境で変わる（＝転移しない）。
"""

from __future__ import annotations

import numpy as np


def make_environments(sigmas=(0.5, 2.0, 5.0), n=3000, seed=0):
    """環境ごとに X1 の分散を変える（介入）。Y の生成則 Y=X1+NY は全環境で不変。"""
    rng = np.random.default_rng(seed)
    envs = []
    for e, s in enumerate(sigmas):
        x1 = rng.normal(0, s, n)                 # ★環境で分散が変わる（介入）
        ny = rng.normal(0, 1.0, n)               # Y のノイズは不変
        y = 1.0 * x1 + ny                        # 真の因果 X1→Y（係数1.0で不変）
        x2 = y + rng.normal(0, 1.0, n)           # X2 は Y の"結果"（子）
        envs.append({"X1": x1, "X2": x2, "Y": y})
    return envs


def invariance_score(envs, predictors) -> dict:
    """予測子集合で Y を回帰し、残差が環境をまたいで不変かを測る。
    不変性スコア = 環境ごとの残差分散のばらつき（変動係数）。小さいほど不変＝因果候補。
    """
    # 全環境プールで係数を1つ推定（同じ関係を全環境に当てる）
    X = np.concatenate([np.column_stack([env[p] for p in predictors]) for env in envs])
    Y = np.concatenate([env["Y"] for env in envs])
    A = np.column_stack([X, np.ones(len(Y))])
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    # 環境ごとの残差分散
    resvars = []
    for env in envs:
        Xe = np.column_stack([env[p] for p in predictors] + [np.ones(len(env["Y"]))])
        r = env["Y"] - Xe @ coef
        resvars.append(float(r.var()))
    resvars = np.array(resvars)
    cv = float(resvars.std() / resvars.mean()) if resvars.mean() > 0 else np.inf
    return {"resvars": resvars, "cv": cv}


if __name__ == "__main__":
    envs = make_environments()
    print("=== 不変性による因果予測（ICP 最小版）===")
    print("  真の因果 X1→Y（不変）、X2 は Y の結果（子）。環境ごとに X1 の分散を変えて介入。\n")
    cands = [["X1"], ["X2"], ["X1", "X2"]]
    res = {tuple(s): invariance_score(envs, s) for s in cands}
    scores = {k: v["cv"] for k, v in res.items()}
    best = min(scores, key=scores.get)
    best_cv = scores[best]

    print(f"  {'予測子':>10}  {'環境ごとの残差分散':>26}  {'不変性CV':>8}  判定")
    for s in cands:
        r = res[tuple(s)]
        rv = "  ".join(f"{v:5.2f}" for v in r["resvars"])
        # 相対判定: 最小CVの ~3倍以内なら不変、それ以上に大きければ環境で変わる
        mark = ("✅不変=因果候補" if r["cv"] <= 3 * best_cv
                else "×環境で変わる=非因果（子/交絡）")
        print(f"  {'+'.join(s):>10}  {rv:>26}  {r['cv']:>8.3f}  {mark}")

    print(f"\n  最も不変な予測子集合 = {'+'.join(best)}")
    print("  判定: " + ("✅ X1 が選ばれた＝多環境で不変な因果＝将来に転移する最良の候補"
                        if best == ("X1",) else "⚠ 期待どおりでない（介入の強さ/ノイズを調整）"))
    print("\n  → 意味: 『将来予測の保証』は当てはめ関数からは出ない。だが")
    print("     『多くの環境で不変な因果関係』は、新しい条件にも転移しやすい＝データ側の外挿の道。")
    print("     実データでは環境=サイト/年/気候レジーム。不変な因果骨格ほど信頼して外挿でき、")
    print("     さらに機構（物理・生物法則）と整合すれば主張が最強になる（不変性×機構×比較設計）。")
