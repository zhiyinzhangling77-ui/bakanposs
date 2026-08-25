"""旗47（前提監査②）：因果十分性の自己矛盾を、潜在交絡を許す LPCMCI で正す。

前提の穴：旗25 が「**未観測の遅い駆動がある**」と結論した瞬間、PCMCI の前提である
**因果十分性（未観測交絡なし）は壊れている**。つまり旗35 の「岩盤リンク」は、厳密には
*同定された因果*ではなく条件付き独立のスキーマにすぎない。自分の発見が自分の推定の前提を壊している。

試みた対処：**LPCMCI（Gerhardus & Runge 2020）は潜在交絡を許す**因果探索で、リンクに型を付ける
（`-->` 因果／`<->` 潜在交絡／`o->`,`o-o` 判断保留）。岩盤リンクが `-->` のまま残るかを見る計画だった。

**結果（重要）：この対処は効かない。** 潜在駆動を仕込んだ合成で検証したところ（ParCorr・既定設定・
n=3000・tau_max=4）：
  ・**`<->` は 1 本も出ない**（4条件すべて）。潜在交絡は `o-o`/`o->` という判断保留として現れる。
  ・**真の因果すら `o->` に格下げ**される＝「原因か、未観測の共通原因か、決められない」。
  ・**潜在駆動が遅いほど、偽リンクが `-->`＝因果と確定されてしまう**（φ=0.85 で顕著）。
    ＝旗25 が示した「遅い未観測駆動」は、この手法群にとって**最悪ケース**。
＝**穴②はアルゴリズムでは塞げない**。正しい塞ぎ方は認識論的：因果骨格を「同定された因果」と呼ぶのを
やめ、**条件付き独立のスキーマ**と言い直す。物理的妥当性（放射→顕熱・潜熱）は先験知識が担保しており、
アルゴリズムの手柄ではない。本ツールはその根拠を実データでなく**制御された合成で示す**装置である。

**計算量の断り**：LPCMCI は PCMCI より遥かに重く、本研究の既定 tau_max=36 × 11 変数は現実的でない。
本ツールは **tau_max を小さく（既定 4）して PCMCI と LPCMCI を同じ tau で走らせ**、公平に比較する。
＝「長いラグのリンク」は評価できない、という限定付きの検証である。

    python research/latent_confounder_step47.py                      # 合成で検証
    python research/latent_confounder_step47.py --site JP-Tak --year 2010
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# 旗35 で「岩盤」とした独立測定間リンク（+ 対照として恒等式リンク GEP→NEE）
BEDROCK = [("Rg", "gH"), ("Rg", "gLE"), ("Rg", "Ta"), ("Ta", "Ts"), ("Rg", "Ts")]
CONTROL = [("GEP", "NEE")]


def _classify(graph, names, src, dst):
    """src→dst のリンク型。**PAG の記号を正しく読む**（旗47 第1回はここを誤読していた）：
      `-->` 祖先関係を主張（因果）／`o->` 終点に矢だけ確定（部分的）／
      `o-o` 隣接するが**向きは決まらない**（＝潜在交絡を排除できない・因果を主張しない）／
      `<->` 潜在交絡が確定。
    第1回は `<->` だけを潜在交絡の証拠とみなし、`o-o` を「向き未定」と軽く扱っていたが、
    **`o-o` こそが「因果と言えない」の主要な出力**である（合成で確認）。"""
    i, j = names.index(src), names.index(dst)
    marks = [str(graph[i, j, t]) for t in range(graph.shape[2])
             if graph[i, j, t] not in ("", None)]
    if not marks:
        return "—"
    for pref in ("-->", "<->", "o->", "o-o"):      # 強い主張から順に拾う
        if pref in marks:
            return pref
    return marks[0]


def _verdict(pc, lp):
    """問いを「どのリンクが交絡か」から**「PCMCI の因果主張が潜在交絡を許しても残るか」**へ。
    LPCMCI は交絡を積極的に同定できるとは限らず（合成では `o-o` を返す）、
    **向きを付けない**ことで「因果と言えない」を表明する。そこを判定に使う。"""
    if pc == "—" and lp == "—":
        return "―どちらも無し"
    if pc != "-->":
        return f"―PCMCIが因果を主張していない（{pc}）"
    if lp == "-->":
        return "★因果主張が残る（潜在交絡を許しても祖先関係）"
    if lp == "<->":
        return "▲潜在交絡が確定＝因果ではない"
    if lp == "—":
        return "▲LPCMCIでは辺自体が消える"
    return f"▲**向きを失う（{lp}）＝潜在交絡を排除できない＝因果と言えない**"


def run_pair(data, names, tau_max, pc_alpha):
    """同じデータ・同じ tau で PCMCI と LPCMCI を走らせ、graph を返す。"""
    from tigramite import data_processing as tp
    from tigramite.pcmci import PCMCI
    from tigramite.lpcmci import LPCMCI
    from tigramite.independence_tests.parcorr import ParCorr

    df = tp.DataFrame(data, var_names=list(names))
    pc = PCMCI(dataframe=df, cond_ind_test=ParCorr(), verbosity=0)
    r_pc = pc.run_pcmciplus(tau_min=0, tau_max=tau_max, pc_alpha=pc_alpha)
    lp = LPCMCI(dataframe=df, cond_ind_test=ParCorr(), verbosity=0)
    r_lp = lp.run_lpcmci(tau_max=tau_max, pc_alpha=pc_alpha)
    return r_pc["graph"], r_lp["graph"]


# ---------- 合成：潜在駆動の「遅さ」を掃引し、両手法が何と言うかを実証 -------------
def _synth(phi=0.85, contemporaneous=False, n=3000, seed=0):
    """H(未観測駆動, AR係数φ)→X,Y ／ A→B(真の因果)。観測は X,Y,A,B のみ。
    contemporaneous=False は H の1ステップ遅れが X,Y を駆動＝旗25 が見つけた状況の型。"""
    rng = np.random.default_rng(seed)
    H = np.zeros(n)
    for t in range(1, n):
        H[t] = phi * H[t - 1] + rng.normal(0, 1)
    X = np.zeros(n); Y = np.zeros(n); A = rng.normal(0, 1, n); B = np.zeros(n)
    for t in range(1, n):
        h = H[t] if contemporaneous else H[t - 1]
        X[t] = 0.6 * h + rng.normal(0, 0.6)
        Y[t] = 0.6 * h + rng.normal(0, 0.6)      # X,Y に直接の因果は無い（偽リンクの素）
        B[t] = 0.7 * A[t - 1] + rng.normal(0, 0.6)   # A→B は本物
    return np.column_stack([X, Y, A, B]), ["X", "Y", "A", "B"]


def _marks(g, names, src, dst):
    i, j = names.index(src), names.index(dst)
    return [str(g[i, j, t]) for t in range(g.shape[2]) if g[i, j, t] not in ("", None)]


def run_synth(tau_max, pc_alpha):
    print("=== 旗47 合成：未観測駆動があるとき、因果探索は何と言うか ===")
    print("  仕込み：H(未観測)→X,Y（X-Y に直接の因果は無い＝偽リンクの素）／A→B は本物。")
    print("  期待：LPCMCI は X-Y を <->（潜在交絡）と印し、A→B は --> のまま残すはず。\n")
    print(f"  {'潜在駆動の型':<26}{'X→Y(偽)':>22}{'A→B(真)':>12}{'<->':>6}")
    for phi, contemp, lab in [(0.0, True, "白色・同時刻 φ=0"),
                              (0.5, True, "中庸・同時刻 φ=0.5"),
                              (0.85, True, "遅い・同時刻 φ=0.85"),
                              (0.85, False, "遅い・ラグ φ=0.85（旗25型）")]:
        data, names = _synth(phi, contemp)
        _, g_lp = run_pair(data, names, tau_max, pc_alpha)
        xy = ",".join(_marks(g_lp, names, "X", "Y")) or "—"
        ab = ",".join(_marks(g_lp, names, "A", "B")) or "—"
        nbi = sum(1 for i in range(len(names)) for j in range(len(names))
                  for t in range(g_lp.shape[2]) if g_lp[i, j, t] == "<->")
        print(f"  {lab:<26}{xy:>22}{ab:>12}{nbi:>6}")
    print("\n  → 結果（ParCorr・既定設定・n=3000・tau_max=4）：")
    print("     ・**<-> は 1 本も出ない**。潜在交絡は <-> でなく o-o / o-> という『判断保留』として現れる。")
    print("     ・**真の因果 A→B すら o-> に格下げ**される＝『A が原因、または未観測の共通原因がある』。")
    print("     ・**潜在駆動が遅いほど、偽リンク X→Y が --> として"
          "『因果と確定』されてしまう**（φ=0.85 で顕著）。")
    print("     ＝我々の状況（旗25＝遅い未観測駆動）は、この手法群にとって最悪ケースにあたる。")
    print("     ＝**穴②はアルゴリズムでは塞げない**。塞ぎ方は認識論的なもの（下記の記録を参照）。")


# ---------- 実データ ---------------------------------------------------------------
def run_real(site, year, months, tau_max, pc_alpha):
    from japanflux_pn.config import AnalysisConfig, RK_VARS
    from japanflux_pn.preprocess import load_corevars_hh

    cfg = AnalysisConfig()
    pre = load_corevars_hh(site, year, months[0], cfg) if len(months) == 1 else None
    if pre is None:
        from japanflux_pn.sites import get_site
        from japanflux_pn.preprocess import (load_raw_all, slice_span_and_anomaly,
                                             PreprocessResult)
        raw = load_raw_all(get_site(site), cfg)
        anom, valid = slice_span_and_anomaly(raw, year, months, cfg)
        pre = PreprocessResult(anomaly=anom, valid=valid, site=site, year=year,
                               month=months[0], config=cfg, months=months)
    if not bool(pre.valid.all()):
        print(f"  {site} {year}: 欠測 {int((~pre.valid).sum())} 点＝完全被覆でない。"
              f"別の年を選ぶこと（PCMCI/LPCMCI とも完全被覆が前提）。")
        return
    data = pre.anomaly[RK_VARS].to_numpy(dtype=float)
    names = list(RK_VARS)
    print(f"=== 旗47 実データ：{site} {year} 月{months} tau_max={tau_max} α={pc_alpha} "
          f"(N={len(data)}) ===")
    print("  LPCMCI は重いので数分〜かかる。同じ tau で PCMCI と並べる。\n", flush=True)
    g_pc, g_lp = run_pair(data, names, tau_max, pc_alpha)

    print(f"  {'リンク':<12} {'PCMCI':>8} {'LPCMCI':>8}  判定")
    print("  --- 旗35 で『岩盤』とした独立測定間リンク ---")
    tally = {}
    for src, dst in BEDROCK:
        pc, lp = _classify(g_pc, names, src, dst), _classify(g_lp, names, src, dst)
        v = _verdict(pc, lp); tally[v[:12]] = tally.get(v[:12], 0) + 1
        print(f"  {src+'→'+dst:<12} {pc:>8} {lp:>8}  {v}")
    print("  --- 対照：恒等式リンク（因果ではない） ---")
    for src, dst in CONTROL:
        pc, lp = _classify(g_pc, names, src, dst), _classify(g_lp, names, src, dst)
        print(f"  {src+'→'+dst:<12} {pc:>8} {lp:>8}  {_verdict(pc, lp)}")

    n_bi = sum(1 for i in range(len(names)) for j in range(len(names))
               for t in range(g_lp.shape[2]) if g_lp[i, j, t] == "<->")
    n_dir = sum(1 for i in range(len(names)) for j in range(len(names))
                for t in range(g_lp.shape[2]) if g_lp[i, j, t] == "-->")
    print(f"\n  === 全体 ===  LPCMCI: 因果 --> {n_dir} 本／潜在交絡 <-> {n_bi} 本")
    for k, v in sorted(tally.items(), key=lambda x: -x[1]):
        print(f"    岩盤リンクの内訳 {k:<14} {v}")
    print("  読み（重要な限定）：合成検証の通り、この手法群は遅い未観測駆動の下で"
          "**偽リンクを --> と確定してしまう**。")
    print("        よって --> が残っても『因果が確認された』とは読めない。逆に o-> / o-o への"
          "格下げは『因果か潜在交絡か決められない』という情報として読める。")
    print(f"  限定：tau_max={tau_max}（計算量の制約）＝これより長いラグのリンクは評価していない。")


def main():
    p = argparse.ArgumentParser(description="潜在交絡を許すLPCMCIで因果骨格を測り直す")
    p.add_argument("--site"); p.add_argument("--year", type=int)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--tau-max", type=int, default=4)
    p.add_argument("--pc-alpha", type=float, default=0.01)
    a = p.parse_args()
    if a.site and a.year:
        run_real(a.site, a.year, a.month, a.tau_max, a.pc_alpha)
    else:
        run_synth(a.tau_max, a.pc_alpha)


if __name__ == "__main__":
    main()
