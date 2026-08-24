"""旗48（前提監査⑥）：out-of-sample 検証を初めて入れる＝「後付けの説明」を予測で試す。

前提の穴：本研究の結果はすべて**同一データ上の記述統計**で、保留サンプル予測が一つも無い。
とくに旗44 の「**水分依存Q10 の符号は生態系依存**（温帯・冷温帯林で正、熱帯排水泥炭・凍土で逆転）」は
**データを見てから気づいた説明**＝典型的な後付けであり、記述統計のままでは反証も確認もできない。

テスト：COSORE の各サイトについて、旗44 の曲率制御後の交互作用係数 d の**符号**をラベルとし、
**サイト属性（緯度・気候・IGBP）だけから、そのサイトを見ずに符号を当てられるか**を
leave-one-**cluster**-out で測る（後述）。比較対象を2つ置く：
  ・**多数決ベースライン**（常に多数派の符号と答える）… これを超えなければ「属性で説明できている」とは言えない
  ・**ラベル並べ替えヌル**（符号をシャッフルして同じ手続き）… 偶然当たる確率の分布
判定：LOO 正答率が多数決を超え、並べ替えヌルに対して p<0.05 なら＝**符号の生態系依存は転移可能な本物**。
超えなければ＝**後付けの物語**であり、そう記す。

分類器は k-NN（標準化属性空間, k=3）＝当てはめの自由度が小さく、少数サイトでも解釈が容易。

**擬似反復の罠（第1回で実際に踏んだ）**：使える属性が緯度・経度・標高＝純粋な地理だったため、
素朴な leave-one-site-out は「**近いサイト同士で当て合う**」だけで高い正答率を出す。実際、負符号
7 サイトのうち 3 つは Hirano のインドネシア泥炭（ほぼ同一地点）だった。そこで **50km 以内を1クラスタと
みなし、予測時にそのクラスタを丸ごと除外する leave-one-cluster-out** を主判定にする。
(1) サイト単位と (2) クラスタ単位の両方を出し、その差＝近接で当て合っていた分を可視化する。

    python research/out_of_sample_step48.py                                    # 合成で検証
    python research/out_of_sample_step48.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

# description.csv から探す属性の候補（実際に在るものだけ使う）
ATTR_CANDIDATES = ["CSR_LATITUDE", "CSR_LONGITUDE", "CSR_ELEVATION",
                   "CSR_MAT", "CSR_MAP", "CSR_ANNUAL_PRECIP", "CSR_ANNUAL_TEMP"]


def _haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def geo_clusters(lat, lon, km=50.0):
    """近接サイトを単連結でまとめる＝擬似反復の単位（同一地点の反復を1つと数える）。"""
    n = len(lat); lab = -np.ones(n, int); c = 0
    for i in range(n):
        if lab[i] >= 0:
            continue
        stack, lab[i] = [i], c
        while stack:
            j = stack.pop()
            for m in range(n):
                if lab[m] < 0 and _haversine(lat[j], lon[j], lat[m], lon[m]) <= km:
                    lab[m] = c; stack.append(m)
        c += 1
    return lab


def loo_knn(X, y, k=3, groups=None):
    """leave-one-out の k-NN 正答率。groups を渡すと**同じクラスタを丸ごと除外**して予測する
    （近接サイト同士で当て合う＝擬似反復による水増しを防ぐ）。"""
    n = len(y)
    if n < k + 2:
        return np.nan
    hit = 0
    for i in range(n):
        d = np.linalg.norm(X - X[i], axis=1)
        if groups is None:
            d[i] = np.inf
        else:
            d[groups == groups[i]] = np.inf     # 自分のクラスタは全部使わない
        if not np.isfinite(d).any():
            continue
        kk = min(k, int(np.isfinite(d).sum()))
        nb = np.argsort(d)[:kk]
        hit += int(np.sign(y[nb].sum() or 1) == y[i])
    return hit / n


def evaluate(X, y, k=3, nperm=2000, seed=0, groups=None):
    acc = loo_knn(X, y, k, groups)
    base = max((y > 0).mean(), (y < 0).mean())          # 多数決ベースライン
    rng = np.random.default_rng(seed)
    null = np.array([loo_knn(X, rng.permutation(y), k, groups) for _ in range(nperm)])
    p = float((null >= acc).mean()) if np.isfinite(acc) else np.nan
    return acc, base, p, null


def report(acc, base, p, y, used, n_used):
    print(f"\n  === 結果（n={len(y)} サイト・属性 {used}）===")
    print(f"  LOO 正答率        {acc:.3f}")
    print(f"  多数決ベースライン  {base:.3f}   （常に多数派と答えた場合）")
    print(f"  並べ替えヌル p      {p:.3f}   （符号をシャッフルして同じ手続き, 2000回）")
    if not np.isfinite(acc):
        print("  → 判定不能（サイト数不足）"); return
    if acc > base and p < 0.05:
        print("  → ★**符号の生態系依存は転移可能**：見ていないサイトの符号を属性から当てられる。")
    elif p < 0.05:
        print("  → ○ヌルは超えるが多数決を超えない＝属性の寄与は弱い。「生態系依存」とまでは言えない。")
    else:
        print("  → ▲**後付けの物語**：見ていないサイトの符号は属性から当てられない。")
        print("     ＝旗44 の『排水泥炭・凍土で逆転』は観察された事実ではあるが、"
              "**一般化できる規則としては未確立**と記すべき。")


# ---------- 合成：転移可能な場合と後付けの場合を作り分ける -------------------------
def _synth(kind, n=36, seed=0):
    rng = np.random.default_rng(seed)
    lat = rng.uniform(-10, 65, n); mat = 25 - 0.35 * lat + rng.normal(0, 2, n)
    X = np.column_stack([lat, mat])
    if kind == "transferable":
        y = np.where(mat > 12, -1.0, 1.0)               # 暖かいほど逆符号＝属性で決まる
        flip = rng.random(n) < 0.10                     # 1割はノイズ
        y[flip] *= -1
    else:
        y = rng.choice([-1.0, 1.0], n)                  # 属性と無関係＝後付け
    return X, y


def run_synth():
    print("=== 旗48 合成検証：転移可能な規則と後付けを見分けられるか ===")
    for kind, lab in [("transferable", "符号が属性で決まる（転移可能）"),
                      ("post_hoc", "符号が属性と無関係（後付け）")]:
        X, y = _synth(kind)
        Xs = (X - X.mean(0)) / X.std(0)
        acc, base, p, _ = evaluate(Xs, y, nperm=1000)
        mark = "★" if (acc > base and p < 0.05) else ("▲" if p >= 0.05 else "○")
        print(f"  {lab:<28} LOO={acc:.3f} 多数決={base:.3f} p={p:.3f}  {mark}")
    print("\n  → 上が★（転移可能）、下が▲（後付け）と出れば検出器は妥当。")


# ---------- 実データ ---------------------------------------------------------------
def run_real(cosore_dir, igbp, months):
    import pandas as pd
    from q10_confound_step44 import analyze
    from cosore_memory_step40 import load_cosore

    root = Path(cosore_dir)
    desc = pd.read_csv(root / "description.csv")
    attrs = [c for c in ATTR_CANDIDATES if c in desc.columns]
    print(f"=== 旗48 実データ：旗44 の『符号は生態系依存』を out-of-sample で試す ===")
    print(f"  description.csv で使える属性：{attrs or '（数値属性なし）'}")
    if not attrs:
        print("  数値属性が無いので予測できない。列名を確認のこと。"); return

    rows = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        if igbp and igbp.lower() not in ig.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f, months)
            if "Tsoil" not in df or "SM" not in df:
                continue
            r = analyze(df["Tsoil"].to_numpy(), df["SM"].to_numpy(), df["Rs"].to_numpy())
        except Exception:
            continue
        if "note" in r or r.get("ci") is None:
            continue
        lo, hi = r["ci"]
        if lo <= 0 <= hi:              # 符号が確定しないサイトはラベルを付けない
            continue
        a = [pd.to_numeric(d.get(c), errors="coerce") for c in attrs]
        if not all(np.isfinite(v) for v in a):
            continue
        rows.append({"ds": ds, "igbp": ig, "y": 1.0 if lo > 0 else -1.0,
                     **{c: float(v) for c, v in zip(attrs, a)}})
        print(f"  {ds:<32}{ig[:12]:<13} d={r['d_quad']:+.4f} 符号={'+' if lo>0 else '−'}",
              flush=True)

    if len(rows) < 8:
        print(f"\n  ラベル付けできたサイトが {len(rows)} 件＝少なすぎて予測を評価できない。")
        return
    df = pd.DataFrame(rows)
    X = df[attrs].to_numpy(float)
    keep = [j for j in range(X.shape[1]) if X[:, j].std() > 0]
    X = X[:, keep]; used = [attrs[j] for j in keep]
    Xs = (X - X.mean(0)) / X.std(0)
    y = df["y"].to_numpy()

    print("\n  --- (1) leave-one-SITE-out（近接サイトを使ってよい＝擬似反復を許す）---")
    acc1, base, p1, _ = evaluate(Xs, y)
    print(f"  LOO 正答率 {acc1:.3f} ／ 多数決 {base:.3f} ／ 並べ替え p {p1:.3f}")

    if "CSR_LATITUDE" in attrs and "CSR_LONGITUDE" in attrs:
        g = geo_clusters(df["CSR_LATITUDE"].to_numpy(float),
                         df["CSR_LONGITUDE"].to_numpy(float), km=50.0)
        sizes = np.bincount(g)
        print(f"\n  --- (2) leave-one-CLUSTER-out（50km以内を1クラスタとして丸ごと除外）---")
        print(f"  {len(sizes)} クラスタ（最大 {sizes.max()} サイト・"
              f"2以上のクラスタ {int((sizes>1).sum())} 個）＝これが実効的な独立標本数")
        acc2, _, p2, _ = evaluate(Xs, y, groups=g)
        print(f"  LOO 正答率 {acc2:.3f} ／ 多数決 {base:.3f} ／ 並べ替え p {p2:.3f}")
        report(acc2, base, p2, y, used, len(y))
        print(f"  ※(1)と(2)の差 {acc1-acc2:+.3f} ＝近接サイト同士で当て合っていた分。")
    else:
        report(acc1, base, p1, y, used, len(y))
    print(f"  内訳：正 {(y>0).sum()} サイト／負 {(y<0).sum()} サイト")
    print("  留保：属性は description.csv にある数値のみ（泥炭排水・凍土という"
          "『我々が後から気づいた区分』そのものは属性に入っていない）。")


def main():
    p = argparse.ArgumentParser(description="後付け説明をout-of-sample予測で試す")
    p.add_argument("--cosore-dir"); p.add_argument("--igbp", default=None)
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()
    if a.cosore_dir:
        run_real(a.cosore_dir, a.igbp, a.month)
    else:
        run_synth()


if __name__ == "__main__":
    main()
