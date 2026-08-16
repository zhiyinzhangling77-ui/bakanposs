"""旗23：呼吸(GER)の"真に未観測な原因"を、非加法込みで観測駆動を使い尽くしてから探す。

北極星：豊かな因果構造→貧しい観測に写らない部分＝未観測。旗9 は GEP で線形残差 R²≈0.84 まで来た。
だが旗20-22 で呼吸は非加法(Ta×θ)と分かった＝**線形残差だと非加法分が"未観測"に化ける**。
so 呼吸では観測駆動を「線形＋主要交互作用(Ta×θ, Ta², θ² 等)」で使い尽くし、
  ・R²(線形) vs R²(＋交互作用)：非加法がどれだけ効くか（旗20 と整合するか）
  ・残差の天井 1−R²：観測(非加法込み)で説明できない割合＝真に未観測 or さらに高次
  ・残差 vs 他の観測変数：入れ忘れの観測駆動か（炭素分割は除外）
  ・残差の自己相関(self-lag 制御)：未観測の"遅い駆動"(基質・フェノロジー・深層土壌水)の足跡か
を見る。旗9(線形・GEP)を非加法対応で呼吸へ。

    python research/respiration_unobserved_step23.py                     # 合成で検証
    python research/respiration_unobserved_step23.py --site JP-Tak --years 1999 ... --month 7 8
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _z(x):
    x = np.asarray(x, float); s = x.std()
    return (x - x.mean()) / s if s > 0 else x - x.mean()


def _r2(y, X):
    """最小二乗 y~X(+切片) の R²。"""
    A = np.column_stack([X, np.ones(len(y))])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    ss = np.sum((y - y.mean()) ** 2)
    return 1 - np.sum(resid ** 2) / ss if ss > 0 else np.nan, resid


def _autocorr1(x):
    x = np.asarray(x, float) - np.mean(x)
    if len(x) < 3 or x.std() == 0:
        return np.nan
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def analyze(Y, drivers: dict, self_lag=1):
    """drivers: {名前: 配列}(観測駆動, 炭素分割は呼び出し側で除外済み)。"""
    names = list(drivers)
    Xlin = np.column_stack([_z(drivers[v]) for v in names])
    r2_lin, _ = _r2(Y, Xlin)
    # 主要交互作用＋二次（Ta×θ 等）。名前にある駆動のみ追加。
    zc = {v: _z(drivers[v]) for v in names}
    extra = []
    def add(a, b, lab):
        if a in zc and b in zc:
            extra.append((lab, zc[a] * zc[b]))
    add("Ta", "th", "Ta*th"); add("Rg", "VPD", "Rg*VPD"); add("Ta", "VPD", "Ta*VPD")
    for v in ("Ta", "th", "VPD"):
        if v in zc:
            extra.append((v + "^2", zc[v] ** 2))
    Xnl = np.column_stack([Xlin] + [e[1] for e in extra]) if extra else Xlin
    r2_nl, resid = _r2(Y, Xnl)
    # 残差 vs 他の観測変数（入れ忘れの観測駆動）
    leftover = []
    for v in names:
        leftover.append((v, abs(float(np.corrcoef(resid, _z(drivers[v]))[0, 1]))))
    leftover.sort(key=lambda t: -t[1])
    # 残差の自己相関（self-lag で Y 自身の記憶を除いてから）
    L = self_lag
    Ycols = [Xnl[L:]]
    for k in range(1, L + 1):
        Ycols.append(Y[L - k:len(Y) - k][:, None])
    A = np.column_stack(Ycols + [np.ones(len(Y) - L)])
    coef, *_ = np.linalg.lstsq(A, Y[L:], rcond=None)
    resid_sl = Y[L:] - A @ coef
    return {"r2_lin": r2_lin, "r2_nl": r2_nl, "resid": resid,
            "leftover": leftover, "resid_autocorr": _autocorr1(resid_sl),
            "n_extra": len(extra)}


def make_synth(hidden=True, n=8000, seed=0):
    rng = np.random.default_rng(seed)
    Rg = rng.normal(0, 1, n); Ta = rng.normal(0, 1, n); VPD = rng.normal(0, 1, n); th = rng.normal(0, 1, n)
    # 遅い未観測駆動（自己相関つき）＝基質/フェノロジー的
    h = np.zeros(n)
    for t in range(1, n):
        h[t] = 0.9 * h[t - 1] + rng.normal(0, 0.4)
    GER = 0.6 * Ta + 0.3 * th + 0.8 * (Ta * th) + (1.0 * h if hidden else 0) + rng.normal(0, 0.3, n)
    return GER, {"Rg": Rg, "Ta": Ta, "VPD": VPD, "th": th}, h


def _report(res, tag):
    print(f"\n  === {tag} ===")
    print(f"  R²(線形)={res['r2_lin']:.3f}  →  R²(＋交互作用{res['n_extra']}項)={res['r2_nl']:.3f}"
          f"  （回復={res['r2_nl']-res['r2_lin']:+.3f}）")
    print(f"  残差の天井 1−R² = {1-res['r2_nl']:.3f}（観測を非加法込みで使っても残る割合）")
    print(f"  残差 vs 観測変数（入れ忘れ候補・上位）: " +
          ", ".join(f"{v}:{c:.2f}" for v, c in res['leftover'][:3]))
    print(f"  残差の自己相関(self-lag後)={res['resid_autocorr']:+.2f}"
          f"（大＝未観測の『遅い』駆動の足跡 / ≈0＝ノイズ）")


def main():
    p = argparse.ArgumentParser(description="呼吸の真に未観測な原因を非加法込みで探す")
    p.add_argument("--site")
    p.add_argument("--years", type=int, nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--self-lag", type=int, default=1)
    a = p.parse_args()

    if not a.site:
        print("=== 旗23 合成検証：呼吸の未観測原因（非加法込み残差）===")
        g1, d1, h1 = make_synth(hidden=True)
        _report(analyze(g1, d1, a.self_lag), "隠れ駆動あり（自己相関つき未観測）")
        g0, d0, _ = make_synth(hidden=False)
        _report(analyze(g0, d0, a.self_lag), "隠れ駆動なし（観測＋交互作用で説明可）")
        print("\n  → 期待：交互作用でR²回復（Ta×θ）。隠れ駆動ありは残差天井が高く残差自己相関が大、")
        print("     なしは残差天井が低く自己相関≈0。＝真に未観測 vs 説明可 を分離。")
        return

    from japanflux_pn.config import RK_VARS
    from japanflux_pn.preprocess import load_corevars_hh
    CARBON = {"GEP", "GER", "NEE"}
    cols = {v: [] for v in RK_VARS}
    used = []
    for y in a.years or []:
        try:
            vf = load_corevars_hh(a.site, y, a.month, None).valid_frame
            for v in RK_VARS:
                cols[v].append(vf[v].to_numpy(float))
            used.append(y)
        except Exception as e:
            print(f"  {y}: SKIP {type(e).__name__}: {e}")
    if not used:
        print("有効年なし"); return
    data = {v: np.concatenate(cols[v]) for v in RK_VARS}
    Y = data["GER"]
    drivers = {v: data[v] for v in RK_VARS if v not in CARBON}  # 炭素分割は除外
    print(f"=== 旗23 実データ {a.site}（GER の未観測原因, プール {len(used)}年, N={len(Y)}）===")
    print(f"  観測駆動（炭素分割 GEP/GER/NEE を除外）: {', '.join(drivers)}")
    _report(analyze(Y, drivers, a.self_lag), f"{a.site} 呼吸 GER")
    print("\n  読み方：R²が交互作用で回復＝非加法(旗20)が効く。残差天井が高く残差自己相関が大なら、")
    print("     観測(非加法込み)でも説明しきれない『遅い未観測駆動』(基質/フェノロジー/深層土壌水)の候補。")
    print("     →そこが領域知識で仮説化し衛星プロキシ(SIF/SMAP等)で検証する入口(北極星)。")


if __name__ == "__main__":
    main()
