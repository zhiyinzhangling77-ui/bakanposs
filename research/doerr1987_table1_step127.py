"""**旗127**：Dörr & Münnich 1987 (Tellus 39B, 114–121) の **Table 1 を一次から書き写し**、
**「水分（夏季降水量）が見かけ Q10 を動かす」という先行の主張が、
彼ら自身の表の中で夏季平均気温と交絡していないか**を数で見る。

## なぜこれをやるのか

**この論文は、A-2（水分依存の見かけ Q10）の最も古い経験的な先行候補である。**
**旗124 は Davidson 1998 の孫引きで「Q10 1.4〜3.1」とだけ記録していた。**
**本周は Tellus B の OA スキャン PDF から本文 8 ページに一次で到達し、Table 1 を読んだ。**

**表には `Nd`（4〜9 月降水量の平年比 %）と `T`（夏季平均気温 °C）が両方載っている。**
**＝先行が「降水のせい」と読んだ Q10 の年々変動が、同じ表の中で気温とも動いていないか
を、我々の A-2 の差分（温度エイリアシングの分離）そのものの形で確かめられる。**

## 何を計算するか（すべて出版された数値のみ。新しい測定は無い）

  1. `Nd` と `T` の相関——**交絡の有無そのもの。**
  2. `Q10` と `Nd`、`Q10` と `T` の相関——**どちらの説明が強いか。**
  3. 両方を入れた重回帰の係数——**n が極小なので参考値としてのみ。**
  4. **並べ替え検定（全組合せの厳密 p）**——**n≤9 で漸近 p を使わないため。**

## やらないこと（★ここを守る）

- **n は最大 9 サイト年・独立な年は 7 つしかない。有意・非有意を主張しない。**
- **我々の 36 データセットの符号と直接比べない**（**量の組も層別の単位も違う**——
  先行は「年をまたぐ夏季降水の平年比」、旗44 は「1 サイトの通年プールを θ でビン分け」）。
- **この計算で先行の結論を否定しない。** 出せるのは「**先行の表からは、降水と気温が
  分離されていない**」という一点だけである。
"""
from __future__ import annotations

import itertools

import numpy as np

from runlog import tee_stdout

# --- Dörr & Münnich 1987, Table 1 を一次（Tellus B の OA スキャン PDF・p.118）から書き写した ------
# 列は論文の定義どおり:
#   T0   : CO2 生成がゼロになる温度 (°C)
#   m    : 回帰直線の傾き (mmol m^-2 h^-1 °C^-1)
#   R2   : 相関係数（論文の表記は R。Bevington 1969 に依拠）
#   Nd   : 夏季降水量 (4–9 月) の長期平均に対する比 (%)
#   Q10  : 5–15 °C 区間で回帰直線から換算した Q10
#   T    : 夏季平均気温 (°C)
#   j0   : 夏季平均 CO2 生成量 (mmol m^-2 h^-1)
TABLE1 = [
    # (site, year,  T0,    m,    R2,   Nd,  Q10,   T,    j0)
    ("NU", 1979, -0.7, 0.66, 0.95, 70, 2.7, 15.4, 10.8),
    ("NU", 1980, -7.3, 0.27, 0.83, 122, 1.8, 15.0, 6.2),
    ("NU", 1981, -1.0, 0.59, 0.81, 89, 2.7, 16.0, 10.1),
    ("NU", 1982, -1.9, 0.59, 0.93, 79, 2.4, 16.8, 11.4),
    ("NU", 1983, -7.1, 0.38, 0.78, 161, 1.8, 17.1, 8.4),
    ("NU", 1984, -1.8, 0.35, 0.96, 119, 2.5, 15.5, 5.5),
    ("SA", 1983, -27.0, 0.15, 0.46, 161, 1.4, 17.1, 5.7),
    ("SA", 1984, 0.2, 0.56, 0.81, 119, 3.1, 15.5, 5.8),
    ("SA", 1985, -6.3, 0.30, 0.60, 102, 1.9, 16.5, 7.1),
]

# NU=Rhine 谷のローム質・非耕作・草地 / SA=砂質土の混交ブナ-トウヒ林。**どちらも温帯**である。
# **`Nd` と `T` は年ごとの気象なので、同じ年の NU と SA は同一値**（表でもそうなっている）。
# ＝**独立な年は 7 つ（1979–1985）しかない。**


def _pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float(np.corrcoef(x, y)[0, 1])


def _spearman(x, y):
    def rank(v):
        v = np.asarray(v, float)
        order = v.argsort()
        r = np.empty(len(v), float)
        r[order] = np.arange(len(v), dtype=float)
        # 同順位は平均順位に潰す
        for val in np.unique(v):
            m = v == val
            r[m] = r[m].mean()
        return r

    return _pearson(rank(x), rank(y))


def _exact_perm_p(x, y, stat=_pearson, n_max=9):
    """**全順列で厳密 p（両側）**。n が小さいので漸近近似を使わない。
    n が `n_max` を超えたら `None` を返す（今回は超えない）。"""
    x, y = list(x), list(y)
    if len(x) > n_max:
        return None
    obs = abs(stat(x, y))
    hits = total = 0
    for perm in itertools.permutations(range(len(y))):
        total += 1
        if abs(stat(x, [y[i] for i in perm])) >= obs - 1e-12:
            hits += 1
    return hits / total


def _ols(y, X, names):
    """定数項つき最小二乗。**n が極小なので係数は参考値**。"""
    X = np.column_stack([np.ones(len(y)), np.asarray(X, float)])
    y = np.asarray(y, float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    out = {"dof": dof, "coef": dict(zip(["const"] + list(names), beta.tolist()))}
    if dof > 0:
        s2 = float(resid @ resid) / dof
        cov = s2 * np.linalg.pinv(X.T @ X)
        out["se"] = dict(zip(["const"] + list(names), np.sqrt(np.diag(cov)).tolist()))
        # 分散拡大係数（回帰子が 2 本なので回帰子間相関から直接）
        r = _pearson(X[:, 1], X[:, 2])
        out["vif"] = 1.0 / (1.0 - r * r)
        out["corr_regressors"] = r
    return out


def _report(label, rows):
    nd = [r[5] for r in rows]
    q10 = [r[6] for r in rows]
    t = [r[7] for r in rows]
    print(f"\n===== {label}（n={len(rows)}）=====")
    for r in rows:
        print(f"  {r[0]} {r[1]}  Nd={r[5]:5.0f}%  T={r[7]:5.1f}C  Q10={r[6]:4.1f}  m={r[3]:.2f}")
    print("\n  -- 交絡の有無（先行が分離していない 2 本の回帰子どうし）--")
    print(f"  corr(Nd, T)  Pearson={_pearson(nd, t):+.3f}  "
          f"Spearman={_spearman(nd, t):+.3f}  厳密p={_exact_perm_p(nd, t)}")
    print("\n  -- Q10 をどちらが説明するか（単独）--")
    print(f"  corr(Q10, Nd) Pearson={_pearson(q10, nd):+.3f}  "
          f"Spearman={_spearman(q10, nd):+.3f}  厳密p={_exact_perm_p(q10, nd)}")
    print(f"  corr(Q10, T)  Pearson={_pearson(q10, t):+.3f}  "
          f"Spearman={_spearman(q10, t):+.3f}  厳密p={_exact_perm_p(q10, t)}")
    if len(rows) >= 5:
        print("\n  -- 両方を入れた重回帰（★n が極小。係数は参考値であって検定ではない）--")
        fit = _ols(q10, np.column_stack([nd, t]), ["Nd", "T"])
        print(f"  dof={fit['dof']}  corr(回帰子)={fit.get('corr_regressors'):+.3f}  "
              f"VIF={fit.get('vif'):.2f}")
        for k in ("const", "Nd", "T"):
            se = fit.get("se", {}).get(k)
            se_s = f" ± {se:.4f}" if se is not None else ""
            print(f"    {k:6s} = {fit['coef'][k]:+.4f}{se_s}")


def main():
    tee_stdout("step127_doerr1987")
    print(__doc__)
    print("\n★ 出典：Dörr, H. & Münnich, K.O. (1987), Tellus 39B, 114–121, Table 1（p.118）。")
    print("★ 本文 8 ページに一次で到達して書き写した数値のみを使う。推定・補完は一切していない。")

    _report("全 9 サイト年（NU 6 年 + SA 3 年）", TABLE1)
    _report("NU のみ（ローム質・草地・1979–1984）", [r for r in TABLE1 if r[0] == "NU"])

    # 気象（Nd・T）は年の量なので、独立な単位は「年」である。
    years = {}
    for r in TABLE1:
        years.setdefault(r[1], []).append(r)
    print("\n\n===== 独立な年は 7 つしかない（1979–1985）=====")
    print("  年ごとの (Nd, T) と、その年に得られている Q10:")
    nd_y, t_y = [], []
    for y in sorted(years):
        rs = years[y]
        nd_y.append(rs[0][5])
        t_y.append(rs[0][7])
        q = ", ".join(f"{r[0]}={r[6]:.1f}" for r in rs)
        print(f"  {y}  Nd={rs[0][5]:5.0f}%  T={rs[0][7]:5.1f}C   Q10: {q}")
    print(f"\n  corr(Nd, T) 年単位 n={len(nd_y)}: Pearson={_pearson(nd_y, t_y):+.3f}  "
          f"Spearman={_spearman(nd_y, t_y):+.3f}  厳密p={_exact_perm_p(nd_y, t_y)}")

    print("\n\n===== 単調性の点検（先行自身が山型を示唆している）=====")
    print("  論文 p.118 本文：多雨年で Q10 は『平年条件の値のほぼ半分まで下がる』(NU)。")
    print("  同 p.118：SA の月別データは『夏季降水が平年より少なくても Q10 は下がる』ことを示唆。")
    print("  ＝ 先行の主張は『乾ほど高 Q10』の単調関係ではなく、平年付近を頂点とする山型である。")
    for site in ("NU", "SA"):
        rows = sorted([r for r in TABLE1 if r[0] == site], key=lambda r: r[5])
        seq = " -> ".join(f"({r[5]:.0f}%, {r[6]:.1f})" for r in rows)
        print(f"  {site}: Nd 昇順の (Nd, Q10) = {seq}")

    print("\n\n===== この計算から言えること／言えないこと =====")
    print("  言える : 先行の表の中で Nd と T は分離されていない（上の corr を見よ）。")
    print("  言える : 先行は Nd と T を同時に入れた検定を行っていない（本文・表とも単変量）。")
    print("  言えない: どちらが真の駆動かは、この n では決まらない。")
    print("  言えない: 我々の 36 データセットの符号との整合／不整合（層別の単位が違う）。")


if __name__ == "__main__":
    main()
