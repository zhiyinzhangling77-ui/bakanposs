"""旗29：呼吸の結論は「昼分割(DT) vs 夜分割(NT)」でひっくり返らないか。

DATA_QUALITY §6 の核心的懸念：GER(生態系呼吸)は NEE を分割アルゴリズムで割った
派生量で、既定は昼分割(DT, RECO_DT)。夜分割(NT, RECO_NT)は別の仮定で NEE を割る。
旗26 の「水分依存 Q10(湿るほど温度感度が上がる)」が **分割法に依存する見かけ** なのか、
**分割法をまたいで残る本物** なのかを直接叩く。

DT と NT で Q10 vs θ の Spearman r を並べて出す。
  ・両方 r>0 → 分割に依らず水分依存 Q10 は本物（結論は分割アーティファクトでない）
  ・符号が割れる → 分割アルゴリズム由来の見かけ（旗26 の結論は要留保）
参考に DT と NT の GER 自体の相関 corr(RECO_DT, RECO_NT) も出す（分割間の一致度）。

RK_VARS には NT 列が無いので生 CSV から RECO_DT/RECO_NT を直接読む。japanflux 形式のみ。

    python research/dt_nt_partition_step29.py                       # 合成で検証
    python research/dt_nt_partition_step29.py --site JP-Tak --qc-max 1
    python research/dt_nt_partition_step29.py --sites JP-Tak JP-Tef JP-Ta2 --qc-max 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from respiration_q10_moisture_step26 import q10_by_moisture, _boot_trend


def _nt_of(dt_col):
    """"RECO_DT_vUT" → "RECO_NT_vUT"。DT を含まなければ None。"""
    return dt_col.replace("_DT_", "_NT_").replace("_DT", "_NT") if "DT" in dt_col else None


def load_dt_nt(site, months, qc_max):
    """生 CSV から RECO_DT/RECO_NT・TA・SWC を直接読む(japanflux 形式)。"""
    import pandas as pd
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.preprocess import (
        _read_table_header, _read_table_columns, find_corevars_files)

    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    vmap = site.var_map()
    dt_col = vmap["GER"]; nt_col = _nt_of(dt_col)
    ta_col = vmap["Ta"]; th_col = vmap["th"]
    files = find_corevars_files(site)
    header0 = _read_table_header(files[0])
    hset = set(header0)
    if nt_col is None or nt_col not in hset:
        return None, dt_col, nt_col, "NT列なし"
    for c in (dt_col, ta_col, th_col):
        if c not in hset:
            return None, dt_col, nt_col, f"{c} なし"

    want = {"TIMESTAMP_START", dt_col, nt_col, ta_col, th_col}
    parts = []
    for f in files:
        df = _read_table_columns(f, want)
        ts = pd.to_datetime(
            pd.to_numeric(df["TIMESTAMP_START"]).astype("int64").astype(str),
            format="%Y%m%d%H%M")
        df = df.drop(columns=["TIMESTAMP_START"]); df.index = ts
        parts.append(df)
    raw = pd.concat(parts)
    raw = raw[~raw.index.duplicated(keep="first")].sort_index()
    raw = raw.replace(cfg.na_sentinel, np.nan)
    raw = raw.rename(columns={dt_col: "GER_DT", nt_col: "GER_NT",
                              ta_col: "Ta", th_col: "th"})
    if months:
        raw = raw[raw.index.month.isin(months)]
    return raw.dropna(), dt_col, nt_col, None


def make_synth(n=60000, seed=0):
    """水分で温度感度が上がる真の呼吸に、DT/NT で別々のノイズを載せる。"""
    rng = np.random.default_rng(seed)
    Ta = rng.uniform(8, 30, n); th = rng.uniform(0.1, 0.5, n)
    thn = (th - 0.3) / 0.2
    b = 0.05 + 0.04 * thn
    true = np.exp(b * (Ta - 20) + 0.5 * thn)
    import pandas as pd
    idx = pd.date_range("2020-07-01", periods=n, freq="30min")
    return pd.DataFrame({
        "GER_DT": np.clip(true * (1 + rng.normal(0, 0.05, n)), 1e-3, None),
        "GER_NT": np.clip(true * (1 + rng.normal(0, 0.08, n)), 1e-3, None),
        "Ta": Ta, "th": th}, index=idx)


def q10_trend(df, col, nbin):
    Ta = df["Ta"].to_numpy(); th = df["th"].to_numpy(); GER = df[col].to_numpy()
    r, ci, _ = _boot_trend(Ta, th, GER, nbin, 200)
    return r, ci


def _mark(r, ci):
    if not isinstance(ci, tuple):
        return f"{r:+.2f} [—]"
    tag = "水分依存Q10(正)" if ci[0] > 0 else "逆" if ci[1] < 0 else "CI0跨ぎ"
    return f"{r:+.2f} [{ci[0]:+.2f},{ci[1]:+.2f}] {tag}"


def main():
    p = argparse.ArgumentParser(description="呼吸の結論は DT/NT 分割でひっくり返るか")
    p.add_argument("--site")
    p.add_argument("--sites", nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--nbin", type=int, default=5)
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    if not a.site and not a.sites:
        print("=== 旗29 合成検証：DT/NT で水分依存 Q10 は一致するか ===")
        df = make_synth()
        for col, lab in [("GER_DT", "昼分割 DT"), ("GER_NT", "夜分割 NT")]:
            r, ci = q10_trend(df, col, a.nbin)
            print(f"  {lab:<10} Q10 vs θ: {_mark(r, ci)}")
        cc = float(np.corrcoef(df['GER_DT'], df['GER_NT'])[0, 1])
        print(f"  corr(DT,NT)={cc:.2f}")
        print("\n  → 真は同じ呼吸なので DT/NT とも r>0(水分依存)で一致するのが期待。")
        return

    from japanflux_pn.sites import get_site
    sites = a.sites or [a.site]
    qtag = f"QC≤{a.qc_max}" if a.qc_max is not None else "gap-fill込み"
    print(f"=== 旗29 実データ 呼吸 Q10 の DT/NT 分割感度（{qtag}, 月={a.month}）===")
    print("  水分依存 Q10(旗26)が分割法をまたいで残るか。両方 r>0＝本物 / 割れる＝分割由来\n")
    print(f"  {'サイト':<8} {'DT: Q10 vs θ':<28} {'NT: Q10 vs θ':<28} {'corr':>5}  一致")
    n_agree = n_tot = 0
    for s in sites:
        try:
            df, dt_col, nt_col, err = load_dt_nt(get_site(s), a.month, a.qc_max)
        except Exception as e:
            print(f"  {s:<8} SKIP {type(e).__name__}: {e}"); continue
        if err:
            print(f"  {s:<8} {err}（DT={dt_col} NT={nt_col}）"); continue
        if len(df) < 3000:
            print(f"  {s:<8} データ不足(N={len(df)})"); continue
        rd, cd = q10_trend(df, "GER_DT", a.nbin)
        rn, cn = q10_trend(df, "GER_NT", a.nbin)
        cc = float(np.corrcoef(df["GER_DT"], df["GER_NT"])[0, 1])
        pos_d = isinstance(cd, tuple) and cd[0] > 0
        pos_n = isinstance(cn, tuple) and cn[0] > 0
        both_ci = isinstance(cd, tuple) and isinstance(cn, tuple)
        agree = "✅一致" if (pos_d and pos_n) else \
                ("△符号は同じ" if both_ci and np.sign(rd) == np.sign(rn) else "⚠割れる")
        if both_ci:
            n_tot += 1; n_agree += int(pos_d and pos_n)
        print(f"  {s:<8} {_mark(rd, cd):<28} {_mark(rn, cn):<28} {cc:>5.2f}  {agree}")
    if n_tot:
        print(f"\n  DT/NT 両方で水分依存 Q10(r>0,CI>0)：{n_agree}/{n_tot}")
        if n_agree == n_tot:
            print("    → ✅ 分割法に依らず結論一致＝旗26 は分割アーティファクトでない")
        elif n_agree >= max(1, n_tot - 1):
            print("    → ○ 大半で一致（例外は分割感度あり）")
        else:
            print("    → ⚠ 分割で割れる＝旗26 の水分依存 Q10 は分割法に依存する疑い")
    print("  留保：DT/NT は同じ NEE を別仮定で割った派生量。相関が高くても独立検証では")
    print("    ない（真の夜間チャンバー実測ではない）。あくまで分割法依存の下限チェック。")


if __name__ == "__main__":
    main()
