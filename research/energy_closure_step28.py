"""旗28：エネルギー収支の閉合(EBR)を site-year の品質指標として使う。

DATA_QUALITY §6 で挙げた「未対応の品質軸」の筆頭。渦相関(EC)は乱流で運ばれる
熱・水・炭素を測るが、乱流が足りないと H+LE が過小評価され、正味放射 Rn−G を
下回る。この不閉合率 EBR = Σ(H+LE) / Σ(Rn−G) は、その site-year がどれだけ
ちゃんと乱流を捉えられていたかの物理的なサイン。EBR が低い(<0.7)＝乱流不足＝
炭素フラックスも取りこぼしている疑い＝除外候補、という**物理の裏付けのある**除外基準。

SSITC(定常性テスト)や u* フィルタは点ごとの品質だが、EBR は site-year 丸ごとの
「測れ具合」を 1 数字で出せる。gap-fill 混じりだと 1 に寄る(補完で埋まる)ため、
--qc-max 1 で実測寄りにして EBR を見るのが本筋。

RK_VARS(11変数)に Rn/G は無いので、生 CSV から直接 H/LE/Rn/G 列を読む。
japanflux 形式のみ対象(BASE/ChinaFlux は列規約が違うため別途)。

    python research/energy_closure_step28.py                       # 合成で検証
    python research/energy_closure_step28.py --site JP-Tak --qc-max 1
    python research/energy_closure_step28.py --sites JP-Tak JP-Tef JP-Ta2 --qc-max 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# H/LE/Rn/G の実カラム候補(japanflux2024 / FLUXNET2015 の表記揺れを吸収)。
# gap-fill 済み(_F_MDS)を優先し、無ければ生列(_1_1_1)へフォールバック。
CAND = {
    "H":  ["H_F_MDS", "H_F", "H_1_1_1", "H"],
    "LE": ["LE_F_MDS", "LE_F", "LE_1_1_1", "LE"],
    "Rn": ["NETRAD_F_MDS", "NETRAD_F", "NETRAD_1_1_1", "NETRAD", "RN", "NET_RAD"],
    "G":  ["G_F_MDS", "G_F_MDS_1", "G_F", "G_1_1_1", "G"],
}


def _resolve(header, cands):
    hset = set(header)
    for c in cands:
        if c in hset:
            return c
    return None


def load_energy(site, months, qc_max):
    """生 CSV から H/LE/Rn/G を直接読み、対象月に絞った 30 分値 DataFrame。"""
    import pandas as pd
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.preprocess import (
        _read_table_header, _read_table_columns, find_corevars_files)

    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    files = find_corevars_files(site)
    cols = {}
    header0 = _read_table_header(files[0])
    for k, cands in CAND.items():
        cols[k] = _resolve(header0, cands)
    # 必須は H/LE/Rn。G(地中熱)は無いサイトがある(例 JP-Ta2)ので任意扱い＝G=0近似。
    missing = [k for k in ("H", "LE", "Rn") if cols[k] is None]
    if missing:
        return None, cols, missing, False
    g_approx = cols["G"] is None
    present = {k: c for k, c in cols.items() if c is not None}

    # QC 列 (存在するもののみ)。--qc-max 指定時に低品質 gap-fill を落とす。
    qc_of = {}
    if qc_max is not None:
        for k, c in present.items():
            qc = c + "_QC"
            if qc in set(header0):
                qc_of[k] = qc
    want = {"TIMESTAMP_START", *present.values(), *qc_of.values()}
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
    # QC>閾値 (低品質補完) を NaN 化してから改名。
    for k, qc in qc_of.items():
        raw[present[k]] = raw[present[k]].where(raw[qc] <= qc_max)
    ren = {v: k for k, v in present.items()}
    raw = raw.rename(columns=ren)
    if g_approx:                       # G 列が無いサイトは G=0 で近似
        raw["G"] = 0.0
    raw = raw[["H", "LE", "Rn", "G"]]
    if months:
        raw = raw[raw.index.month.isin(months)]
    return raw, cols, [], g_approx


def closure(df):
    """EBR=Σ(H+LE)/Σ(Rn−G) と、原点通し回帰の傾き・R²・n を返す。"""
    ok = df[["H", "LE", "Rn", "G"]].notna().all(axis=1)
    d = df.loc[ok]
    turb = (d["H"] + d["LE"]).to_numpy()
    avail = (d["Rn"] - d["G"]).to_numpy()
    n = len(d)
    if n < 200:
        return None
    ebr = float(turb.sum() / avail.sum()) if avail.sum() != 0 else np.nan
    # 原点通し回帰 turb = s·avail
    s = float((avail @ turb) / (avail @ avail)) if (avail @ avail) != 0 else np.nan
    ss_res = float(((turb - s * avail) ** 2).sum())
    ss_tot = float(((turb - turb.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return {"ebr": ebr, "slope": s, "r2": r2, "n": n}


def by_year(df):
    out = {}
    for y, g in df.groupby(df.index.year):
        c = closure(g)
        if c:
            out[int(y)] = c
    return out


def make_synth(kind, n=40000, seed=0):
    rng = np.random.default_rng(seed)
    avail = np.abs(rng.normal(300, 120, n))            # Rn−G [W/m^2]
    G = rng.normal(20, 10, n)
    frac = 1.0 if kind == "closed" else 0.75           # 良閉合 vs 25%不閉合
    turb = frac * avail + rng.normal(0, 25, n)
    import pandas as pd
    idx = pd.date_range("2020-07-01", periods=n, freq="30min")
    return pd.DataFrame({"H": turb * 0.4, "LE": turb * 0.6,
                         "Rn": avail + G, "G": G}, index=idx)


def verdict(ebr):
    if not np.isfinite(ebr):
        return "—"
    if ebr < 0.7:
        return "⚠ 不閉合(乱流不足の疑い/除外候補)"
    if ebr <= 1.05:
        return "✅ 良好(0.7–1.05)"
    return "△ 過閉合(Rn/G測器 or 貯留無視)"


def main():
    p = argparse.ArgumentParser(description="エネルギー収支閉合を品質指標に")
    p.add_argument("--site")
    p.add_argument("--sites", nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    if not a.site and not a.sites:
        print("=== 旗28 合成検証：エネルギー収支閉合 EBR ===")
        for kind, lab in [("closed", "良閉合(真値1.00)"), ("open", "25%不閉合(真値0.75)")]:
            c = closure(make_synth(kind))
            print(f"  {lab:<22} EBR={c['ebr']:.2f} 傾き={c['slope']:.2f} "
                  f"R²={c['r2']:.2f}  {verdict(c['ebr'])}")
        print("\n  → 良閉合 EBR≈1.0、不閉合 EBR≈0.75 が復元できれば OK。")
        return

    from japanflux_pn.sites import get_site
    sites = a.sites or [a.site]
    qtag = f"QC≤{a.qc_max}" if a.qc_max is not None else "gap-fill込み"
    print(f"=== 旗28 実データ エネルギー収支閉合（{qtag}, 月={a.month}）===")
    print("  EBR=Σ(H+LE)/Σ(Rn−G)。<0.7=乱流不足の疑い(除外候補), 0.7–1.05=良好\n")
    print(f"  {'サイト':<8} {'EBR':>5} {'傾き':>5} {'R²':>5} {'N':>7}  判定")
    for s in sites:
        try:
            df, cols, missing, g_approx = load_energy(get_site(s), a.month, a.qc_max)
        except Exception as e:
            print(f"  {s:<8} SKIP {type(e).__name__}: {e}"); continue
        if missing:
            print(f"  {s:<8} 必須列欠如 {missing}（H/LE/Rn が無く EBR 算出不可）"); continue
        c = closure(df)
        if c is None:
            print(f"  {s:<8} データ不足"); continue
        gtag = "  ⚠G無し(G=0近似)" if g_approx else ""
        print(f"  {s:<8} {c['ebr']:>5.2f} {c['slope']:>5.2f} {c['r2']:>5.2f} "
              f"{c['n']:>7}  {verdict(c['ebr'])}{gtag}")
        yr = by_year(df)
        bad = [y for y, cc in yr.items() if np.isfinite(cc['ebr']) and cc['ebr'] < 0.7]
        if yr:
            evals = [cc['ebr'] for cc in yr.values() if np.isfinite(cc['ebr'])]
            print(f"           年別 EBR: {min(evals):.2f}–{max(evals):.2f}"
                  f"（{len(yr)}年）" + (f"  ⚠ <0.7 の年: {bad}" if bad else ""))
        print(f"           使用列: H={cols['H']} LE={cols['LE']} "
              f"Rn={cols['Rn']} G={cols['G'] if cols['G'] else '（無し→0近似）'}")
    print("\n  留保：貯留補正(canopy/soil storage)を引いていないので夏の EBR は")
    print("    やや低めに出うる。gap-fill 込みは補完で 1 に寄るため --qc-max 1 が本筋。")
    print("    G無しサイトは G=0 近似（分母 Rn−G の G を落とす＝分母が大きく＝EBR はやや"
          "低めに出る。真値はもう少し高い）。")


if __name__ == "__main__":
    main()
