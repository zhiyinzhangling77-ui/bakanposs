"""旗41：主要な順列検定にFDR(Benjamini-Hochberg)補正＝穴⑤(多重比較)を叩く。

旗13/15 の事前指定・層別順列検定の p 値を1つの族として BH-FDR で補正し、
q=0.05 / 0.10 で生き残る主張を出す。効果量+ブートCIで判定した発見(θ→γLE 水マスター・
Bowen反転・チャンバーメモリ・SIF棄却)は別枠(独立サイト再現+CI)なのでここには入れない。

    python research/fdr_correction_step41.py
"""
from __future__ import annotations

# (ラベル, p値[両側], 事前指定か) — セッションの実データ順列検定出力より
TESTS = [
    ("旗13 Rg→GEP vs VPD (乾で光合成脱結合, 層別併合)", 0.009, True),
    ("旗15 炭素コア冗長 vs VPD (乾で冗長弱まる, 事前指定)", 0.033, True),
    ("旗13 Rg→GEP vs 土壌水分 (VPDの裏返し)",           0.11,  False),
    ("旗13 θ→GER vs VPD (呼吸律速の状態依存)",          0.24,  True),
    ("旗15 呼吸相乗 vs VPD (状態非依存=ヌル, 事前指定)",   0.27,  True),
    ("旗15 光合成サブ系 vs VPD (撤回済)",                0.71,  True),
]


def bh(pvals, q):
    """Benjamini-Hochberg：p(i) ≤ (i/m)q を満たす最大 i までを有意とする。"""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    kmax = -1
    for rank, idx in enumerate(order, 1):
        if pvals[idx] <= rank / m * q:
            kmax = rank
    survive = set(order[:kmax]) if kmax > 0 else set()
    return survive, kmax


def main():
    ps = [t[1] for t in TESTS]
    m = len(ps)
    print(f"=== 旗41 FDR補正（Benjamini-Hochberg, 族サイズ m={m}）===\n")
    print(f"  {'検定':<44}{'p':>7}  素の0.05  BH-crit(q=.05/.10)")
    order = sorted(range(m), key=lambda i: ps[i])
    rank = {idx: r for r, idx in enumerate(order, 1)}
    for i, (lab, p, pre) in enumerate(TESTS):
        r = rank[i]
        c05, c10 = r / m * 0.05, r / m * 0.10
        raw = "有意" if p < 0.05 else "・"
        print(f"  {lab:<44}{p:>7.3f}  {raw:>4}     {c05:.3f}/{c10:.3f}")
    for q in (0.05, 0.10):
        surv, k = bh(ps, q)
        names = [TESTS[i][0].split(' (')[0] for i in sorted(surv, key=lambda i: ps[i])]
        print(f"\n  BH-FDR q={q}: 生存 {len(surv)}/{m} → {names or '（なし）'}")
    print("\n  結論：素で有意だった2件(0.009, 0.033)は、6検定の族としてFDR補正すると")
    print("    q=0.10 では両方生存(乾→光合成脱結合・乾→炭素コア冗長弱化)、q=0.05 では非生存（周辺的）。")
    print("    ＝状態依存の主張は『中程度の証拠・要年数積み増し』が正直。ヌル(呼吸相乗 状態非依存)は補正後も頑健。")
    print("  留保：効果量+ブートCIの発見(θ→γLE水マスター・Bowen反転・チャンバー~4日メモリ・SIF棄却)は")
    print("    独立サイト再現+CIの別枠＝この族に含めない。それらは単一p値でなく多サイト一貫性で支える。")


if __name__ == "__main__":
    main()
