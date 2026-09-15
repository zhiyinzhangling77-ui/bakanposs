"""**旗157 の後知恵（★数を見てから書いた。事前登録には無い）。**

**この 2 つは判定ではない。** 事前登録 step157 の判定は `year_share_none_step157.py` が出し切っており、
**D-1 の合否も H29〜H31 の勝敗も、本スクリプトでは一切動かさない。**
**ここでやるのは「D-1 の逆転を、どう読んではいけないか」を詰めることだけである。**

## A：**D-1 の逆転は、大きさで選んだことの機械的な帰結ではないか**

**事前登録 2 節は「取り分 `s` は `|Δ̂|` の大きさと交絡しない」と書いた。それは正しくない。**
**有意で選ぶことは `|Δ̂|` の大きい側で選ぶことであり、`Δ̂ = Δ_b + Δ_w` の二成分のうち
分散の大きいほう（本周は `Δ_w`・77%）が `|Δ̂|` を大きくしやすい。**
**＝帰無が完全に較正されていても、有意群の `s` は低いほうに偏る。**
**だから「`s` が有意群で小さい」だけでは「年内成分が偽陽性を運んでいる」と言えない。**

**切り分け**：**`|Δ̂|` を揃えてなお `s` が効くか**を見る。
  ・A-1 `|Δ̂|` の上位 `n_sig` 本（**`p` を使わず大きさだけで選ぶ**）の `s` と、有意群の `s` を比べる
  ・A-2 `sig ~ log|Δ̂| + s` のロジスティック回帰で、`|Δ̂|` を入れたあと `s` の係数が残るか

## B：**帰無が狭いぶんだけで、旗156 の 0.073／0.097 が出るか**

**D-5 は成分ごとの SD しか出していない**（`Δ_b` 比 0.63・`Δ_w` 比 1.0）。
**`Δ*` そのものの SD は記録していない**ので、ここで測る。
そのうえで、**「帰無が `c` 倍だけ狭い正規分布」という粗い近似で偽陽性率を計算し**、
旗156 の実測（0.073・0.097）と並べる。**近似なので一致を期待しない。桁が合うかだけを見る。**

    .venv/bin/python research/posthoc_step157.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import phase_asymmetry_step139 as M                      # noqa: E402
from year_share_none_step157 import _perm_replay, SEED0, YEARS, LOGS  # noqa: E402

PARTS = LOGS / "step157_parts_none.csv"
N_B, NPERM_B = 200, 500        # B：後知恵なので本数は「数分で終わる」を基準に選んだ


def a_checks(t: pd.DataFrame) -> None:
    print("\n===== A：D-1 の逆転は大きさで選んだことの機械的な帰結か（**後知恵**）=====")
    for k in ("h", "le"):
        s = t[f"s_{k}"].to_numpy(float)
        ad = np.abs(t[f"delta_{k}"].to_numpy(float))
        g = t[f"csv_sig_{k}"].to_numpy(float) > 0.5
        n_sig = int(g.sum())
        top = np.zeros(len(t), bool)
        top[np.argsort(-ad)[:n_sig]] = True
        print(f"\n  θ→γ{k.upper()}（有意 {n_sig} 本）")
        print(f"    A-1 有意群の s 中央 {np.median(s[g]):.3f}"
              f" ／ **|Δ̂| 上位 {n_sig} 本**の s 中央 {np.median(s[top]):.3f}"
              f" ／ 残り {len(t) - n_sig} 本 {np.median(s[~top]):.3f}")
        print(f"        二群の重なり：有意かつ上位 {int((g & top).sum())} 本"
              f"（有意の {(g & top).sum() / max(n_sig, 1):.0%}）")
        u = stats.mannwhitneyu(s[top], s[~top], alternative="two-sided")
        print(f"        **`p` を使わず大きさだけで選んでも、s は同じ向きに下がるか**："
              f"上位 対 残り の両側 p = {u.pvalue:.2e}"
              f"／rank-biserial {2 * u.statistic / (top.sum() * (~top).sum()) - 1:+.3f}")

        # A-2：|Δ̂| を入れたあと s の係数が残るか（ロジスティック・Newton 法・切片つき）
        X = np.column_stack([np.ones(len(t)), np.log(ad + 1e-12), s])
        y = g.astype(float)
        b = np.zeros(3)
        for _ in range(60):
            p = 1.0 / (1.0 + np.exp(-(X @ b)))
            W = p * (1 - p) + 1e-12
            H = X.T @ (X * W[:, None])
            b = b + np.linalg.solve(H, X.T @ (y - p))
        p = 1.0 / (1.0 + np.exp(-(X @ b)))
        se = np.sqrt(np.diag(np.linalg.inv(X.T @ (X * (p * (1 - p) + 1e-12)[:, None]))))
        z = b / se
        pv = 2 * (1 - stats.norm.cdf(np.abs(z)))
        for nm, i in (("切片", 0), ("log|Δ̂|", 1), ("**s**", 2)):
            print(f"        A-2 {nm:>8}：係数 {b[i]:+.3f}（SE {se[i]:.3f}・z {z[i]:+.2f}"
                  f"・p {pv[i]:.3g}）")
        print(f"        → **`|Δ̂|` を入れたあと `s` は "
              f"{'まだ効いている' if pv[2] < 0.05 else '**効かなくなる**'}**")


def b_check(t: pd.DataFrame) -> None:
    print(f"\n===== B：`Δ*` そのものの SD と、帰無の狭さで出る偽陽性率（**後知恵**）=====")
    print(f"  {N_B} 本 × {NPERM_B} 置換（**D-5 と同じ replicate の並び・種も同じ**）")
    sd_null = {"h": [], "le": []}
    for i in range(N_B):
        d = M.synth("none", years=YEARS, seed=SEED0 + i, thin=False)
        cols, blocks = M._pool_arrays(d)
        rng = np.random.default_rng(i)
        acc = {"h": [], "le": []}
        for _ in range(NPERM_B):
            sp_parts, au_parts = [], []
            for idx, kk in blocks:
                perm = rng.permutation(idx)
                sp_parts.append(perm[:kk])
                au_parts.append(perm[kk:])
            star = M._delta_from_idx(cols, np.sort(np.concatenate(sp_parts)),
                                     np.sort(np.concatenate(au_parts)))
            for k in ("h", "le"):
                if np.isfinite(star[k][2]):
                    acc[k].append(star[k][2])
        for k in ("h", "le"):
            sd_null[k].append(float(np.std(acc[k])))
        if (i + 1) % 50 == 0:
            print(f"    … {i + 1}/{N_B} 本")
    obs156 = {"h": 0.073, "le": 0.097}
    for k in ("h", "le"):
        nn = float(np.median(sd_null[k]))
        oo = float(np.std(t[f"delta_{k}"].to_numpy(float)))
        c = nn / oo
        # **粗い近似**：帰無が c 倍狭いなら、両側 0.05 の境は観測の分布で 1.96c に当たる
        rate = float(2 * (1 - stats.norm.cdf(1.959964 * c)))
        print(f"\n  θ→γ{k.upper()}：帰無 SD(Δ*) 中央 {nn:.4f} ／ 600 本の SD(Δ̂) {oo:.4f}"
              f" → **比 {c:.3f}**")
        print(f"    **正規近似で出る偽陽性率 {rate:.3f}** ／ 旗156 の実測 {obs156[k]:.3f}")
        print("    **近似なので一致は期待しない。桁が合うかだけを見る。**")


def main() -> int:
    if not PARTS.exists():
        print(f"**{PARTS} が無い**——先に `year_share_none_step157.py --all` を走らせる")
        return 1
    t = pd.read_csv(PARTS)
    print("**旗157 の後知恵**（★数を見てから書いた・事前登録には無い・判定を動かさない）")
    print(f"**標本：{PARTS.name}（{len(t)} 行）**")
    a_checks(t)
    b_check(t)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
