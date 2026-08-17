"""旗32：各変数は独立観測か、他変数から算出された派生量か（各サイトで監査）。

問題意識：ある変数がもとから他変数の計算結果なら、その2変数の「関係」が出るのは当たり前
（＝定義による自動相関）。どの関係が"発見"でどれが"定義"かを、各フラックスサイトで切り分ける。

コードから確定している派生関係を、実データで定量化する：
  1. **炭素の恒等式**：NEE = GER − GEP（分割の定義）。NEE を GER,GEP で回帰し R²≈1・係数≈(+1,−1)なら
     ＝GER/GEP/NEE は独立3変数でなく実質1自由度（NEE 測定＋分割で2つに割っただけ）。
  2. **VPD は気温の関数**：VPD = es(Ta)·(1−RH/100)。VPD を es(Ta) で回帰した R² が高いほど
     VPD は Ta で決まる（残りが RH）。＝VPD は独立変数として綺麗でない。
  3. **総合の冗長度**：各変数を「残り全部」で線形予測した R²。1 に近い変数は独立情報が薄い
     （派生 or 強共線）＝関係を見ても"当たり前"になりやすい。

これで「θ→炭素はセーフ／背骨の炭素同士リンクは定義／VPD は Ta 混じり」を各サイトで裏取りする。

    python research/derivation_audit_step32.py                     # 合成で検証
    python research/derivation_audit_step32.py --sites JP-Tak CN-HaM MN-Kbu JP-Mse --qc-max 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

RK = ["Rg", "Ta", "VPD", "Ts", "P", "th", "gH", "gLE", "GER", "NEE", "GEP"]


def _r2(y, X):
    """y を X（1D 配列のリスト, 切片込み）で最小二乗回帰した R² と係数。"""
    y = np.asarray(y, float)
    A = np.column_stack([np.asarray(x, float) for x in X] + [np.ones(len(y))])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    ss = np.sum((y - y.mean()) ** 2)
    r2 = 1 - np.sum(resid ** 2) / ss if ss > 0 else np.nan
    return r2, coef


def _es(Ta):
    """飽和水蒸気圧 es(Ta) [hPa]（Magnus-Tetens）。"""
    return 6.1094 * np.exp(17.625 * Ta / (Ta + 243.04))


def audit(df):
    """派生関係の定量化。df は RK 列を持つ生データ（欠測除去済み）。"""
    out = {}
    # 1. 炭素の恒等式 NEE ~ GER + GEP
    if all(c in df for c in ("NEE", "GER", "GEP")):
        r2, coef = _r2(df["NEE"].to_numpy(),
                       [df["GER"].to_numpy(), df["GEP"].to_numpy()])
        out["carbon"] = {"r2": r2, "b_GER": coef[0], "b_GEP": coef[1]}
    # 2. VPD ~ es(Ta)
    if all(c in df for c in ("VPD", "Ta")):
        es = _es(df["Ta"].to_numpy())
        r2, coef = _r2(df["VPD"].to_numpy(), [es])
        out["vpd"] = {"r2": r2, "slope": coef[0]}
    # 3. 各変数の「残り全部」に対する R²（冗長度）
    red = {}
    cols = [c for c in RK if c in df]
    M = {c: df[c].to_numpy(float) for c in cols}
    # 標準化（スケール差を除く）
    Z = {c: (M[c] - M[c].mean()) / (M[c].std() or 1) for c in cols}
    for v in cols:
        others = [Z[c] for c in cols if c != v]
        r2, _ = _r2(Z[v], others)
        red[v] = r2
    out["redundancy"] = red
    return out


def make_synth(n=8000, seed=0):
    """NEE=GER−GEP・VPD=es(Ta)(1−RH/100)・独立θ を仕込んだ合成。"""
    import pandas as pd
    rng = np.random.default_rng(seed)
    Ta = 20 + rng.normal(0, 5, n)
    RH = np.clip(70 + rng.normal(0, 15, n), 10, 100)
    VPD = _es(Ta) * (1 - RH / 100)                     # 定義通り Ta,RH から
    Rg = np.clip(300 + rng.normal(0, 80, n), 0, None)
    th = np.clip(0.3 + rng.normal(0, 0.05, n), 0.05, 0.6)   # 独立
    Ts = Ta - 1 + rng.normal(0, 0.5, n)
    GEP = np.clip(0.01 * Rg + rng.normal(0, 0.5, n), 0, None)
    GER = np.exp(0.06 * (Ta - 20)) + rng.normal(0, 0.2, n)
    NEE = GER - GEP + rng.normal(0, 1e-3, n)            # 分割の定義（ほぼ厳密）
    return pd.DataFrame({"Rg": Rg, "Ta": Ta, "VPD": VPD, "Ts": Ts,
                         "P": rng.normal(0, 0.5, n), "th": th,
                         "gH": Rg * 0.2, "gLE": Rg * 0.3,
                         "GER": GER, "NEE": NEE, "GEP": GEP})


def _report(res, tag):
    print(f"\n  === {tag} ===")
    c = res.get("carbon")
    if c:
        flag = "★定義(独立でない)" if c["r2"] > 0.98 else "分割に誤差あり"
        print(f"  炭素 NEE~GER+GEP: R²={c['r2']:.3f} 係数 GER={c['b_GER']:+.2f} "
              f"GEP={c['b_GEP']:+.2f}  {flag}")
    v = res.get("vpd")
    if v:
        print(f"  VPD~es(Ta):      R²={v['r2']:.3f}（残り{(1-v['r2'])*100:.0f}%が RH）"
              f"  {'★Ta依存が強い' if v['r2'] > 0.6 else 'RH寄与大'}")
    red = res.get("redundancy", {})
    if red:
        ordered = sorted(red.items(), key=lambda kv: -kv[1])
        print("  冗長度 R²(各変数|残り全部)  高い＝独立情報が薄い:")
        print("   " + "  ".join(f"{k}={val:.2f}" for k, val in ordered))


def main():
    p = argparse.ArgumentParser(description="各変数は独立観測か派生量か（サイト監査）")
    p.add_argument("--sites", nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    if not a.sites:
        print("=== 旗32 合成検証：派生関係を定量化できるか ===")
        _report(audit(make_synth()), "合成（NEE=GER−GEP, VPD=es(Ta)(1−RH/100), θ独立）")
        print("\n  → 炭素 R²≈1・係数≈(+1,−1)＝定義／VPD R²は中〜高（Ta依存）／θ の冗長度は低い（独立）が期待。")
        return

    import pandas as pd
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    from japanflux_pn.run_robustness import get_site_years
    qtag = f"QC≤{a.qc_max}" if a.qc_max is not None else "gap-fill込み"
    print(f"=== 旗32 実データ 各変数の素性（独立観測か派生量か, {qtag}, 月={a.month}）===")
    for s in a.sites:
        try:
            cfg = AnalysisConfig(qc_max=a.qc_max) if a.qc_max is not None else AnalysisConfig()
            years, mo = get_site_years(s)
            ms = sorted(a.month or mo)
            raw = load_raw_all(get_site(s), cfg)
            raw = raw[raw.index.month.isin(ms)]
            df = raw[[c for c in RK if c in raw.columns]].dropna()
        except Exception as e:
            print(f"\n  {s}: SKIP {type(e).__name__}: {e}"); continue
        if len(df) < 500:
            print(f"\n  {s}: データ不足({len(df)})"); continue
        _report(audit(df), f"{s}（N={len(df)}）")
    print("\n  読み方：炭素 R²≈1＝GER/GEP/NEE は実質1自由度（炭素同士の『関係』は定義, 発見でない）。")
    print("    VPD の R² が高いほど VPD は Ta で決まる＝VPD を独立駆動として使うと Ta 混じり。")
    print("    冗長度が高い変数は独立情報が薄い＝その変数との関係は『当たり前』になりやすい。")
    print("  →アトラスの軸は、両変数が独立測定 or 派生入力に相手を含まないペアに限るべき"
          "（θ→炭素はセーフ, 炭素同士・VPD絡みは要注意）。")


if __name__ == "__main__":
    main()
