"""旗34：炭素コア O-info の「冗長」は NEE=GER−GEP の定義で膨らんでいないか（C の検証）。

既存の主要結論：炭素コア {Rg,GEP,GER,NEE} は冗長支配（z>0）＝放射共通駆動の背骨。だが旗32 で
**GEP は NEE に係数 −1 で入る＝定義的リンク**が確定した。so この冗長は "本物の共通駆動" でなく
"NEE=GER−GEP の代数" で膨らんでいる疑いがある。切り分ける：

  1. 合成：GEP・GER を独立乱数、NEE=GER−GEP（＋微ノイズ）、Rg 独立。**共通原因が無いのに**
     {Rg,GEP,GER,NEE} が冗長 z>0 を出すなら＝代数だけで冗長が捏造される（懸念は本物）。
  2. 実データ：各サイトで O-info z を 3 サブ系で比較
     ・炭素三角 {Rg,GEP,GER,NEE}（疑わしい）
     ・エネルギー {Rg,Ta,γH,γLE}（炭素変数ゼロ＝綺麗な背骨）
     ・単一炭素 {Rg,Ta,θ,GER}（炭素1つ・三角なし）
     綺麗なサブ系でも冗長が出れば＝背骨は本物（代数のせいでない）。三角だけなら＝膨らみ。

    python research/carbon_algebra_oinfo_step34.py                       # 合成で検証
    python research/carbon_algebra_oinfo_step34.py --sites JP-Tak CN-HaM MN-Kbu --qc-max 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from o_information_step14 import _digitize, o_information, surrogate_z

M = 6           # 粗ビン（既存 oinfo_analysis と同じ）
NSURR = 200

SUBSYS = {
    "炭素三角{Rg,GEP,GER,NEE}": ["Rg", "GEP", "GER", "NEE"],
    "エネルギー{Rg,Ta,γH,γLE}": ["Rg", "Ta", "gH", "gLE"],
    "単一炭素{Rg,Ta,θ,GER}": ["Rg", "Ta", "th", "GER"],
}


def oinfo_z(df, cols, m=M, nsurr=NSURR):
    """df の指定列を等幅ビン化し O-info と シャッフル z を返す。"""
    digs = [_digitize(df[c].to_numpy(float), m) for c in cols]
    obs, mu, z = surrogate_z(digs, m, nsurr, correct=True, seed=0)
    return obs, z


def _judge(z):
    if not np.isfinite(z):
        return "—"
    if z > 2.36:
        return "冗長支配(z>0)"
    if z < -2.36:
        return "★相乗支配(z<0)"
    return "・非有意"


def make_synth(kind, n=12000, seed=0):
    rng = np.random.default_rng(seed)
    if kind == "algebra":        # 共通原因なし・NEE=GER−GEP の代数だけ
        Rg = rng.normal(0, 1, n)
        GEP = rng.normal(0, 1, n); GER = rng.normal(0, 1, n)
        NEE = GER - GEP + rng.normal(0, 1e-3, n)
        return {"Rg": Rg, "GEP": GEP, "GER": GER, "NEE": NEE}
    if kind == "common":         # 放射共通駆動（本物の冗長）
        Rg = rng.normal(0, 1, n)
        return {"Rg": Rg, "Ta": Rg + rng.normal(0, .3, n),
                "gH": Rg + rng.normal(0, .3, n), "gLE": Rg + rng.normal(0, .3, n)}
    # independent
    return {k: rng.normal(0, 1, n) for k in ("Rg", "Ta", "gH", "gLE")}


def _synth_df(d):
    import pandas as pd
    return pd.DataFrame(d)


def load_pooled(site, months, qc_max):
    """健全年の夏アノマリ（valid_frame）をプールした DataFrame。"""
    import pandas as pd
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.preprocess import load_corevars_hh
    from japanflux_pn.run_robustness import get_site_years
    years, mo = get_site_years(site)
    ms = sorted(months or mo)
    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    frames = []
    for y in years:
        try:
            pr = load_corevars_hh(site, y, ms, cfg)
            vf = pr.valid_frame
            if len(vf) > 200:
                frames.append(vf)
        except Exception:
            continue
    if not frames:
        return None
    return pd.concat(frames)


def main():
    p = argparse.ArgumentParser(description="炭素コアO-infoの冗長は定義で膨らむか")
    p.add_argument("--sites", nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    if not a.sites:
        print("=== 旗34 合成検証：代数だけで冗長が捏造されるか ===")
        for kind, lab, cols in [
            ("algebra", "代数のみ NEE=GER−GEP（共通原因なし）", ["Rg", "GEP", "GER", "NEE"]),
            ("common", "放射共通駆動（本物の冗長）", ["Rg", "Ta", "gH", "gLE"]),
            ("independent", "独立4本（帰無）", ["Rg", "Ta", "gH", "gLE"])]:
            df = _synth_df(make_synth(kind))
            obs, z = oinfo_z(df, cols)
            print(f"  {lab:<32} Ω={obs:+.3f} z={z:+7.1f}  {_judge(z)}")
        print("\n  → 代数のみでも z>0 なら『炭素三角の冗長は定義で膨らむ』が実証される。")
        print("    共通駆動は本物の冗長 z>0、独立は z≈0 が期待。")
        return

    print(f"=== 旗34 実データ 炭素コアの冗長は定義由来か（QC≤{a.qc_max}, 月={a.month}）===")
    print("  各サブ系の O-info z（>0冗長/<0相乗）。炭素三角(疑)とエネルギー・単一炭素(綺麗)を比較。\n")
    hdr = "".join(f"{k:>26}" for k in SUBSYS)
    print(f"  {'サイト':<8}{hdr}")
    rows = []
    for s in a.sites:
        try:
            df = load_pooled(s, a.month, a.qc_max)
        except Exception as e:
            print(f"  {s:<8} SKIP {type(e).__name__}: {e}"); continue
        if df is None or len(df) < 1000:
            print(f"  {s:<8} データ不足"); continue
        cells = []
        rec = {"site": s}
        for name, cols in SUBSYS.items():
            if all(c in df for c in cols):
                obs, z = oinfo_z(df, cols)
                rec[name] = z
                cells.append(f"z={z:+6.1f} {_judge(z)[:6]:>8}")
            else:
                cells.append("—")
        rows.append(rec)
        print(f"  {s:<8}" + "".join(f"{c:>26}" for c in cells))

    # まとめ：綺麗なサブ系でも冗長が出るか
    print("\n  === まとめ ===")
    for name in SUBSYS:
        zs = [r[name] for r in rows if name in r and np.isfinite(r[name])]
        if zs:
            nred = sum(1 for z in zs if z > 2.36)
            print(f"  {name:<26} 冗長 {nred}/{len(zs)}  平均z={np.mean(zs):+.1f}")
    print("\n  読み方：エネルギー・単一炭素（炭素三角なし）でも冗長が出れば＝背骨は本物（代数でない）。")
    print("    炭素三角だけ突出して冗長なら＝その分は NEE=GER−GEP の定義で膨らんでいる＝割引が必要。")


if __name__ == "__main__":
    main()
