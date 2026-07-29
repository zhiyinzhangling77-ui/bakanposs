"""共通駆動の条件付け診断 (妥当性チェック)。

ペアワイズ MI/TE は共通駆動 (放射 Rg が 11 変数を同時駆動) を分離できないため、
「同期支配」という結論が共通駆動のアーティファクトである恐れがある。ここでは各ペア
(X,Y) について:

    I(X;Y)         … 素の相互情報 (共通駆動込み)
    I(X;Y | Rg)    … 駆動変数 Rg を与えた上での直接依存

を比較し、「Rg を除いても残る共有情報は何割か」を測る。条件付け後も有意なペアだけが
Rg で説明できない直接結合の候補。

    python -m japanflux_pn.condition_driver --site JP-Tak --year 2003 --month 7 8

驚くほど多くのペアが Rg 条件付けで消えるなら、同期支配の相当部分は共通駆動の影で
あり、次は PCMCI/多変量 TE で本格的に条件付けするべき、という判断材料になる。
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS
from .preprocess import load_corevars_hh, PreprocessResult
from . import information_theory as it


@dataclass
class ConditioningResult:
    driver: str
    mi: pd.DataFrame          # I'(X;Y) [%]、対称
    cmi: pd.DataFrame         # I'(X;Y|driver) [%]、対称
    mi_sig: pd.DataFrame      # I(X;Y) が有意 (bool)
    cmi_sig: pd.DataFrame     # I(X;Y|driver) が有意 (bool)
    config: AnalysisConfig

    def save(self, outdir) -> None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.mi.to_csv(outdir / f"MI_pct.csv")
        self.cmi.to_csv(outdir / f"CMI_given_{self.driver}_pct.csv")
        self.cmi_sig.to_csv(outdir / f"CMI_given_{self.driver}_significant.csv")


def condition_on_driver(
    pre: PreprocessResult, driver: str = "Rg"
) -> ConditioningResult:
    """全ペアの I(X;Y) と I(X;Y|driver) を計算し、条件独立ヌルで有意性判定。"""
    cfg = pre.config
    m = cfg.n_bins
    vf = pre.valid_frame
    idx = {v: it.digitize_series(vf[v].to_numpy(dtype=float), m) for v in RK_VARS}
    z_cols = [idx[driver]]
    rng = np.random.default_rng(cfg.seed)
    logm = cfg.log_m

    nan = np.nan
    mi = pd.DataFrame(nan, index=RK_VARS, columns=RK_VARS, dtype=float)
    cmi = pd.DataFrame(nan, index=RK_VARS, columns=RK_VARS, dtype=float)
    mi_sig = pd.DataFrame(False, index=RK_VARS, columns=RK_VARS)
    cmi_sig = pd.DataFrame(False, index=RK_VARS, columns=RK_VARS)

    # Miller-Madow 補正で 2D(MI) と 3D(CMI) の推定バイアスを揃え、条件付けで
    # 見かけ上 MI が増える (負の drop) 次元アーティファクトを抑える。
    MM = True
    for a, b in combinations(RK_VARS, 2):
        i_ab = it.mutual_information_indices(idx[a], idx[b], m, MM)
        mi_stats = it.surrogate_mi_stats(idx[a], idx[b], m, cfg.n_surrogates,
                                          cfg.sig_c, rng, MM)
        mi.loc[a, b] = mi.loc[b, a] = 100.0 * i_ab / logm
        sig_i = i_ab > mi_stats["threshold"]
        mi_sig.loc[a, b] = mi_sig.loc[b, a] = bool(sig_i)

        if a == driver or b == driver:
            continue  # 駆動変数自身を条件付けるのは自明
        c_ab = it.conditional_mutual_information_indices(idx[a], idx[b], z_cols, m, MM)
        c_stats = it.surrogate_cmi_stats(idx[a], idx[b], z_cols, m,
                                         cfg.n_surrogates, cfg.sig_c, rng, MM)
        cmi.loc[a, b] = cmi.loc[b, a] = 100.0 * c_ab / logm
        sig_c = c_ab > c_stats["threshold"]
        cmi_sig.loc[a, b] = cmi_sig.loc[b, a] = bool(sig_c)

    return ConditioningResult(driver=driver, mi=mi, cmi=cmi,
                              mi_sig=mi_sig, cmi_sig=cmi_sig, config=cfg)


def report(site: str, year: int, months: list[int], driver: str = "Rg",
           config: AnalysisConfig | None = None,
           outroot: str | Path | None = None) -> ConditioningResult:
    config = config or AnalysisConfig()
    pre = load_corevars_hh(site, year, months, config)
    print(f"[preprocess] {site} {year}-{pre.month_label}: n_points={pre.n_points}")
    res = condition_on_driver(pre, driver)

    pairs = list(combinations([v for v in RK_VARS if v != driver], 2))
    sig_pairs = [(a, b) for a, b in pairs if bool(res.mi_sig.loc[a, b])]
    # 減少率 (drop) を主指標に。>= DROP_CUT なら「大半が共通駆動で説明」。
    DROP_CUT = 60.0

    def drop_pct(a, b) -> float:
        i, c = res.mi.loc[a, b], res.cmi.loc[a, b]
        return 100.0 * (1 - c / i) if i > 0 else 0.0

    n_common = sum(1 for a, b in sig_pairs if drop_pct(a, b) >= DROP_CUT)
    print(f"\n=== 共通駆動 {driver} の条件付け診断 "
          f"(driver 以外の {len(pairs)} ペア) ===")
    print(f"  MI 有意ペア: {len(sig_pairs)}")
    print(f"  うち {driver} 条件付けで MI が {DROP_CUT:.0f}%以上 減少 "
          f"(=大半が共通駆動の影): {n_common}")
    print(f"  残り (直接依存が相当残る候補): {len(sig_pairs) - n_common}")
    print(f"  ※ ビン条件付けは連続 {driver} の残差を取りきれないため、減少率を主指標"
          f"とする。厳密には KNN-CMI (Tigramite) で追試。\n")

    rows = sorted(((drop_pct(a, b), res.mi.loc[a, b], res.cmi.loc[a, b], a, b)
                   for a, b in sig_pairs), reverse=True)
    print(f"  {'pair':<12} {'I%':>6} {'I|'+driver+'%':>8} {'drop':>6}  判定")
    for dp, i_ab, c_ab, a, b in rows:
        tag = f"共通駆動 ({driver}で説明)" if dp >= DROP_CUT else "DIRECT 候補 (残存大)"
        print(f"  {a+'-'+b:<12} {i_ab:6.1f} {c_ab:8.1f} {dp:5.0f}%  {tag}")

    if outroot is not None:
        outdir = Path(outroot) / f"{site}_{year}{pre.month_label}_cond{driver}"
        res.save(outdir)
        print(f"\n[output] {outdir}")
    return res


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="共通駆動の条件付け診断")
    p.add_argument("--site", required=True)
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7])
    p.add_argument("--driver", default="Rg", help="条件付ける共通駆動変数 (既定 Rg)")
    p.add_argument("--bins", type=int, default=None)
    p.add_argument("--surrogates", type=int, default=None)
    p.add_argument("--qc-max", type=int, default=None,
                   help="QC≤この値のみ残す (0=実測, 1=実測+良質補完)。既定 None=gap-fill込み")
    p.add_argument("--outroot", default=None)
    args = p.parse_args(argv)

    kw = {}
    if args.bins is not None:
        kw["n_bins"] = args.bins
    if args.surrogates is not None:
        kw["n_surrogates"] = args.surrogates
    if args.qc_max is not None:
        kw["qc_max"] = args.qc_max
    config = AnalysisConfig(**kw)
    report(args.site, args.year, args.month, args.driver, config, args.outroot)


if __name__ == "__main__":
    main()
