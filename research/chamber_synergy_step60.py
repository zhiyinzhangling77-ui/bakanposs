"""旗60（策B の最後の項目）：呼吸の「相乗」を、チャンバーの多深度データ＝測定量だけで測る。

`MEASURED_ONLY_SPINE.md` の区分 B（計算量に依存したまま）に残っている主張のうち、
**呼吸の相乗（旗14：{Rg,Ta,θ,GER} が O-info で相乗支配、mean z=−7.1・19/21）**は
GER＝分割派生量に依存している。旗56 でタワー側メモリを測定量へ移そうとして失敗したので、
「相乗も無理」と書いたが、根拠は**「チャンバーは土壌温度と水分の2変数しかない」**だった。

**それは正しくない**：COSORE には**土壌温度を複数深度で測るサイト**がある
（例：SAVAGE_hf006-03 は 8 深度、CARBONE_SC_EMBUDO は T 4深度＋SM 3深度）。
＝**{Rs, T浅, T深, SM}** の4変数がそろい、**分割も穴埋めも通さずに O-information が計算できる**。

  ・チャンバーでも相乗支配なら → 旗14 の主張が**測定量だけで裏取りされる**（B→A へ移せる）。
  ・そうでなければ → 旗14 の相乗は**分割派生量に固有**の可能性＝そう記す。

**判定は必ず z（シャッフルヌルからの距離）で行う。Ω の生の符号で判定してはならない**——
4変数・8ビン・N数百では有限標本バイアスでヌルの平均自体が負になり、Ω<0 は相乗を意味しない。
旗14 も z で報告している（mean z=−7.1）。旗60 第1回はここを誤り、結論を逆に出した。

前処理はタワー側と揃える：**5日アノマリ**（値から5日中心移動平均を引く）＝季節共変動を除く。
O-information と シャッフルヌルからの z は本研究の実装（`japanflux_pn.information_theory`）をそのまま使う。

    python research/chamber_synergy_step60.py                                    # 合成で検証
    python research/chamber_synergy_step60.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore

BINS = 8
NSUR = 200
MIN_N = 400


def _anomaly(s, win=5):
    """5日アノマリ＝値 − 5日中心移動平均（タワー側の前処理と揃える）。"""
    return s - s.rolling(win, center=True, min_periods=max(2, win // 2)).mean()


def omega_z(cols, bins=BINS, nsur=NSUR, seed=0, block_len=None):
    """O-information Ω とヌルからの z。Ω<0＝相乗支配, Ω>0＝冗長支配。

    ``block_len`` を渡すと**自己相関を保つブロック並べ替え**のヌルを使う。
    素の並べ替え（``block_len=None``）は**自己相関のある系列では正しくない**——
    旗72 の監査で **偽陽性率 5%→27%・z が「冗長」側へ偏る**ことが確定している
    （AR(1) φ=0.8・4変数・8ビン・N=500・反復60 で z の平均 +1.21）。
    """
    from japanflux_pn import information_theory as it
    idx = [it.digitize_series(np.asarray(c, float), bins) for c in cols]
    om = it.o_information_indices(idx, bins, correct=True)
    rng = np.random.default_rng(seed)
    if block_len:
        st = it.surrogate_o_information_stats_block(
            idx, bins, nsur, 0.0, rng, correct=True, block_len=block_len)
    else:
        st = it.surrogate_o_information_stats(idx, bins, nsur, 0.0, rng, correct=True)
    z = (om - st["mu"]) / st["sigma"] if st["sigma"] > 0 else np.nan
    return float(om), float(z)


def depth_cols(df_raw_cols, prefix):
    out = []
    for c in df_raw_cols:
        m = re.fullmatch(rf"CSR_{prefix}(\d+\.?\d*)", c)
        if m:
            out.append((float(m.group(1)), c))
    return sorted(out)


def load_multi(path, months=None):
    """{Rs, T浅, T深, SM}（SM が多深度なら {Rs, T浅, SM浅, SM深}）の日次アノマリを返す。"""
    raw = pd.read_csv(path)
    cols = list(raw.columns)
    if "CSR_FLUX_CO2" not in cols:
        return None, None
    tc = "CSR_TIMESTAMP_BEGIN" if "CSR_TIMESTAMP_BEGIN" in cols else "CSR_TIMESTAMP_END"
    ts = pd.to_datetime(raw[tc], errors="coerce")
    T, S = depth_cols(cols, "T"), depth_cols(cols, "SM")
    if len(T) < 2 or len(S) < 1:
        return None, {"T": [d for d, _ in T], "SM": [d for d, _ in S]}
    sh = min(T, key=lambda x: abs(x[0] - 5))
    dp = T[-1]
    if dp[1] == sh[1]:
        return None, {"T": [d for d, _ in T], "SM": [d for d, _ in S]}
    smsh = min(S, key=lambda x: abs(x[0] - 5))
    d = pd.DataFrame({"Rs": pd.to_numeric(raw["CSR_FLUX_CO2"], errors="coerce").to_numpy(),
                      "T_sh": pd.to_numeric(raw[sh[1]], errors="coerce").to_numpy(),
                      "T_dp": pd.to_numeric(raw[dp[1]], errors="coerce").to_numpy(),
                      "SM": pd.to_numeric(raw[smsh[1]], errors="coerce").to_numpy()},
                     index=ts)
    d = d[d.index.notna()]
    if months:
        d = d[d.index.month.isin(months)]
    daily = d.groupby(d.index.normalize()).mean()
    grid = pd.date_range(daily.index.min(), daily.index.max(), freq="D")
    daily = daily.reindex(grid)
    anom = daily.apply(_anomaly).dropna()
    meta = {"T": [d0 for d0, _ in T], "SM": [d0 for d0, _ in S],
            "used": (sh[0], dp[0], smsh[0])}
    return (anom if len(anom) >= MIN_N else None), meta


# ---------- 合成検証 ---------------------------------------------------------------
def _synth(kind, n=1500, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.normal(0, 1, n); b = rng.normal(0, 1, n); c = rng.normal(0, 1, n)
    if kind == "redundant":                       # 共通駆動＝冗長支配 (Ω>0)
        h = rng.normal(0, 1, n)
        x1, x2, x3 = h + 0.4 * a, h + 0.4 * b, h + 0.4 * c
        y = h + 0.4 * rng.normal(0, 1, n)
    else:                                          # XOR 的な相互作用＝相乗支配 (Ω<0)
        x1, x2, x3 = a, b, c
        y = np.sign(a) * np.sign(b) * np.abs(c) + 0.3 * rng.normal(0, 1, n)
    return [y, x1, x2, x3]


def run_synth():
    print("=== 旗60 合成検証：Ω の符号で相乗／冗長を見分けられるか ===")
    for kind, lab in [("redundant", "共通駆動（冗長支配 Ω>0 のはず）"),
                      ("synergy", "掛け算的相互作用（相乗支配 Ω<0 のはず）")]:
        om, z = omega_z(_synth(kind))
        print(f"  {lab:<34} Ω={om:+.4f}  z={z:+.1f}  "
              f"→ {'冗長' if om > 0 else '相乗'}")
    print("\n  → 上が Ω>0、下が Ω<0 なら推定は妥当（タワー側と同じ実装を使っている）。")


# ---------- 実データ ---------------------------------------------------------------
def run_real(cosore_dir, igbp, months, block_lens=(0, 5, 10, 20)):
    root = Path(cosore_dir); desc = pd.read_csv(root / "description.csv")
    print("=== 旗60 実データ：チャンバー多深度で呼吸の相乗を測る（測定量だけ）===")
    print("  変数 {Rs, T浅, T深, SM}（5日アノマリ）。Ω<0＝相乗支配、Ω>0＝冗長支配。")
    print(f"  比較：タワー側の呼吸系 {{Rg,Ta,θ,GER}} は **相乗 19/21・mean z=−7.1**（旗14, GER依存）")
    print("  **旗72 の修正**：素の並べ替えヌルは自己相関を壊し、偽陽性率 5%→27%・z が冗長側へ偏る。")
    print("  → **ブロック長を変えた複数のヌルを併記**し、**結論がブロック長に依らないか**を見る。")
    print("     ブロック長 0＝従来の素の並べ替え（**参考値・信用しない**）。\n")
    zc = "".join(f"{('z(素)' if L == 0 else f'z(塊{L})'):>9}" for L in block_lens)
    print(f"  {'dataset':<30}{'深度':>14}{'N':>6}{'Ω':>9}{zc}  判定(塊10)")
    tally = {L: [0, 0, 0] for L in block_lens}       # [相乗, 冗長, 有意でない]
    n_syn = n_red = n_ns = 0
    skipped = 0
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        if igbp and igbp.lower() not in ig.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            anom, meta = load_multi(f, months)
        except Exception:
            continue
        if anom is None:
            skipped += 1
            continue
        cols = [anom[c].to_numpy() for c in ("Rs", "T_sh", "T_dp", "SM")]
        try:
            zs = {}
            for L in block_lens:
                om, zz = omega_z(cols, block_len=(L or None))
                zs[L] = zz
                if np.isfinite(zz):
                    tally[L][0 if zz < -2 else (1 if zz > 2 else 2)] += 1
                else:
                    tally[L][2] += 1
        except Exception:
            continue
        z = zs.get(10, np.nan)          # **判定はブロック長10（呼吸の4日より長い）で行う**
        # **判定は z で行う**（Ω の生の符号ではない）。4変数・8ビン・N数百では有限標本バイアスで
        # シャッフルヌルの平均自体が負になるため、Ω を 0 と比べるのは無効。旗14 も z で報告している。
        # （旗60 第1回はここを誤り、Ω の符号で「相乗5/6」と出していた＝結論が逆転する誤り）
        syn = np.isfinite(z) and z < -2
        red = np.isfinite(z) and z > 2
        n_syn += int(syn); n_red += int(red); n_ns += int(not (syn or red))
        u = meta["used"]
        zstr = "".join(f"{zs[L]:>9.1f}" if np.isfinite(zs[L]) else f"{'—':>9}"
                       for L in block_lens)
        print(f"  {ds:<30}{f'{u[0]:.0f},{u[1]:.0f},{u[2]:.0f}':>14}{len(anom):>6}"
              f"{om:>9.4f}{zstr}  "
              f"{'★相乗' if syn else ('·冗長' if red else '△有意でない')}", flush=True)
    print(f"\n  === まとめ ===")
    print(f"  多深度がそろって計算できたサイト：{n_syn + n_red + n_ns}（不足で除外 {skipped}）")
    if n_syn + n_red + n_ns == 0:
        print("  → **測定量だけでは計算できない**＝旗14 は B（計算量依存）に留まる。")
        return
    print(f"  **ヌルの作り方ごとの内訳（相乗／冗長／有意でない）**")
    for L in block_lens:
        t = tally[L]
        tag = "素の並べ替え（**信用しない**）" if L == 0 else f"ブロック長 {L}"
        print(f"    {tag:<28} {t[0]:>3} ／ {t[1]:>3} ／ {t[2]:>3}")
    print(f"  判定に採るのは**ブロック長10**：★相乗 {n_syn}／·冗長 {n_red}／△有意でない {n_ns}")
    print("  **ブロック長を変えて結論が変わるなら、その結論はヌルの作り方に依存している**＝そう記す。")
    print("\n  読み：相乗が多数＝**旗14 の相乗支配が測定量だけで裏取りされた**＝B から A へ移せる。")
    print("        冗長が多数／割れる＝相乗は**分割派生量に固有**の可能性＝旗14 は B に留めそう記す。")
    print("  留保：チャンバーは土壌呼吸のみ・点測定で、タワーの生態系呼吸とは対象が違う。")
    print("        変数も {Rs,T浅,T深,SM} で {Rg,Ta,θ,GER} と一対一ではない＝**同型の検定であって同一ではない**。")
    print("        **旗72 の帰結は両刃**：素のヌルは z を冗長側へ偏らせるので、")
    print("        **冗長は偽陽性が出やすく、相乗は逆に隠れる**。旗60 第2回の『相乗0／冗長3』は")
    print("        **どちらの向きにも安全でなかった**＝本実行がその作り直しである。")


def main():
    p = argparse.ArgumentParser(description="チャンバー多深度で呼吸の相乗を測る")
    p.add_argument("--cosore-dir"); p.add_argument("--igbp", default=None)
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()
    if a.cosore_dir:
        run_real(a.cosore_dir, a.igbp, a.month)
    else:
        run_synth()


if __name__ == "__main__":
    main()
