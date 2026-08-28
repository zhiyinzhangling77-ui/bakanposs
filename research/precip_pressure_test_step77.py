"""旗77：**降水と気圧でメモリは説明できるか**——六つ目と七つ目の対抗仮説。

旗76 で確かめたこと：**同一地点4組すべてで、チャンバー観測期間の全日に**
タワーの**降水**が重なり（年 997〜1742 mm ＝日本の森林として妥当な値域）、
**気圧 `PA_F`** も在る。＝**新規データなしで、二つの対抗仮説を検定できる**。

## ⑥ 降水（＝旗45/75 の最大の限定を外す）

旗45 は **Birch 効果**（乾いた土に雨が降ると呼吸が跳ねる）を検定したが、
**COSORE に降水が無いため θ の増加を代理**にせざるを得なかった。旗75 も同じ限定を書いた。
ここでは**本物の降水**で、しかも**Birch の形（乾燥期間 × 降雨量）**を直接作って当てる。

## ⑦ 気圧（＝物理側の対抗仮説）

**気圧変動はチャンバー測定に既知のアーティファクトを生む**（pressure pumping：
気圧が下がると土壌孔隙から CO₂ が吸い出される）。これは**生物ではなく物理**の機構である。
気圧とその変化率を当てて残差メモリが落ちるなら、**メモリの一部は物理**と言える。

## 足す列（段階的）

  ・**＋降水**：当日 P・遅れ 1–3 日・積算 7/30 日
  ・**＋Birch**：**最後の降雨からの日数（乾燥期間）**・**乾燥期間 × 当日降雨**（＝Birch の形）
  ・**＋気圧**：当日 PA・**日変化 ΔPA**・**|ΔPA|**・遅れ 1–2 日
  ・**＋全部**

基準線は旗74 の最も豊かな段階（T×W テンソルビン・同時刻のみ）。
**外挿残差**（時間ブロック交差検証・ブロック中心化）で測り、
**プラセボ**（タワー側の列だけを位相ずらしして作り直す・4 通りのうち最も下げた回）を併走させる。

## 合成での判別試験

  ・`hidden`：気象と無関係な隠れ過程のメモリ → **どれを足しても落ちないはず**
  ・`birch`：**降雨イベント後に減衰するパルス**でメモリを作る → **Birch 列で落ちるはず**
  ・`pump`：**気圧変化に比例した見かけの変動** → **気圧列で落ちるはず**

    python research/precip_pressure_test_step77.py --synth
    python research/precip_pressure_test_step77.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore, _acf_gap, _efold_gap
from model_richness_step74 import design, residuals
from same_site_arc_step66 import PAIRS

SETS = ("同時刻のみ", "＋降水", "＋Birch", "＋気圧", "＋全部")
SHIFTS = (100, 180, 260, 500)
RAIN_MM = 1.0          # 「雨が降った日」の閾値


def dryspell(P, thr=RAIN_MM):
    """**最後の降雨からの日数**（Birch 効果の中心量）。"""
    out = np.full(len(P), np.nan)
    d = np.nan
    for i, p in enumerate(P):
        if not np.isfinite(p):
            out[i] = d
            continue
        d = 0.0 if p > thr else (d + 1.0 if np.isfinite(d) else 0.0)
        out[i] = d
    return out


def _lags(x, lags):
    out = []
    for L in lags:
        v = np.full(len(x), np.nan)
        v[L:] = x[:-L]
        out.append(v)
    return out


def _cum(x, wins):
    s = pd.Series(x)
    return [s.shift(1).rolling(w, min_periods=max(3, w // 3)).sum().to_numpy() for w in wins]


def tower_cols(P, PA, kind):
    """段階ごとに足す列（**タワー由来のものだけ**）。"""
    cols = []
    if kind in ("＋降水", "＋全部") and P is not None:
        cols += [P] + _lags(P, (1, 2, 3)) + _cum(P, (7, 30))
    if kind in ("＋Birch", "＋全部") and P is not None:
        ds = dryspell(P)
        # **雨が降る「直前」の乾燥期間を使う**（旗77 第1版の誤り＝自分の道具の欠陥15件目）。
        # `dryspell` は**雨の日に 0 へリセット**するので、同じ日の添字で `ds × P` を作ると
        # **降雨のある日は必ず 0** になり、列が構造上ほぼゼロで意味を持たなかった。
        # 合成（Birch を植てた場合）で「説明せず」と出たことで気づいた——
        # **実データでも黙って偽の否定を返していた**ところだった。
        ds_prev = np.concatenate([[np.nan], ds[:-1]])
        pulse = ds_prev * np.nan_to_num(P, nan=0.0)   # Birch＝(雨前の乾燥期間) × 降雨
        cols += [ds, ds_prev, pulse, np.log1p(np.clip(ds, 0, None))]
        # **蓄積の形も入れる**（合成で捕まえた設計の穴）：
        # Birch 効果は「降雨の瞬間」ではなく**数日かけて減衰するパルス**として現れる。
        # 当日の量だけでは表現できないので、**減衰時定数を数通り**用意して当てる。
        # 真の時定数は分からないので τ=2/5/10 日を並べ、どれが効くかはデータに決めさせる。
        for tau in (2.0, 5.0, 10.0):
            dec = np.exp(-1.0 / tau)
            acc = np.zeros(len(P))
            pv = np.nan_to_num(pulse, nan=0.0)
            for i in range(1, len(P)):
                acc[i] = dec * acc[i - 1] + pv[i]
            cols.append(acc)
    if kind in ("＋気圧", "＋全部") and PA is not None:
        d1 = np.concatenate([[np.nan], np.diff(PA)])
        cols += [PA, d1, np.abs(d1)] + _lags(PA, (1, 2)) + _lags(d1, (1,))
    return cols


def measure(y, T, W, P, PA, kind, shift=None):
    base = design("テンソルビン", T, W)
    if base is None:
        return None
    if kind != "同時刻のみ":
        Ps, PAs = P, PA
        if shift is not None:
            # **タワー側の生の列をずらしてから派生量を作る**（派生量をずらすと構造が壊れる）
            Ps = np.roll(P, shift) if P is not None else None
            PAs = np.roll(PA, shift) if PA is not None else None
        ex = tower_cols(Ps, PAs, kind)
        if ex:
            base = np.hstack([base, np.column_stack(ex)])
    with np.errstate(divide="ignore", invalid="ignore"):
        ly = np.log(np.where(y > 0, y, np.nan))
    res = residuals(base, ly, True)
    if res is None or np.isfinite(res).sum() < 40:
        return None
    ss = np.nansum((ly - np.nanmean(ly)) ** 2)
    r2 = float(1 - np.nansum(res ** 2) / ss) if ss > 0 else np.nan
    a1 = _acf_gap(res, 1)
    v = float(np.nanvar(res))
    # **記憶量＝ラグ1 の自己共分散**（＝ACF1 × 残差分散）。
    # **ACF1 だけでは対象概念を測れない**（旗61 と同じ形の誤り）：
    # 説明変数が効いて残差が小さくなると、**残りかすの相対的な自己相関はむしろ上がる**。
    # 実際、合成で Birch 列は R² を 0.14→0.79 に上げたのに ACF1 は 0.49→0.63 に上がった。
    # 「自己相関している分散が**どれだけ減ったか**」を見るには**絶対量**が要る。
    mem = a1 * v if (np.isfinite(a1) and np.isfinite(v)) else np.nan
    return {"r2": r2, "acf1": a1, "efold": _efold_gap(res), "var": v, "mem": mem,
            "n": int(np.isfinite(res).sum())}


def placebo_best(y, T, W, P, PA, kind):
    """**記憶量を最も減らした**シフトを採る（プラセボに有利な、厳しい比較）。"""
    best = None
    for sh in SHIFTS:
        if sh >= len(T):
            continue
        m = measure(y, T, W, P, PA, kind, shift=sh)
        if m is None or not np.isfinite(m.get("mem", np.nan)):
            continue
        if best is None or m["mem"] < best["mem"]:
            best = m
    return best


def analyze(y, T, W, P, PA):
    base = measure(y, T, W, P, PA, "同時刻のみ")
    if base is None or not np.isfinite(base.get("mem", np.nan)):
        return None
    rows = {"同時刻のみ": (base, None)}
    for kind in SETS[1:]:
        if kind == "＋気圧" and PA is None:
            rows[kind] = (None, None); continue
        if kind in ("＋降水", "＋Birch") and P is None:
            rows[kind] = (None, None); continue
        rows[kind] = (measure(y, T, W, P, PA, kind),
                      placebo_best(y, T, W, P, PA, kind))
    return rows


def _cut(base, m):
    """**記憶量の削減率**（1.0＝完全に消えた・0＝変わらない・負＝増えた）。"""
    if m is None or not np.isfinite(m.get("mem", np.nan)) or not np.isfinite(base["mem"]):
        return np.nan
    if base["mem"] <= 0:
        return np.nan
    return 1.0 - m["mem"] / base["mem"]


BEAT = 0.20        # プラセボの削減率をこれだけ上回れば「説明する」


def show(rows):
    base = rows["同時刻のみ"][0]
    print(f"    {'足した列':<12}{'ACF1':>7}{'記憶量':>9}{'**削減率**':>11}"
          f"{'プラセボ':>10}{'R²':>7}  判定")
    print(f"    {'同時刻のみ':<12}{base['acf1']:>7.2f}{base['mem']:>9.4f}"
          f"{'—':>11}{'—':>10}{base['r2']:>7.2f}")
    for kind in SETS[1:]:
        m, pl = rows.get(kind, (None, None))
        if m is None or not np.isfinite(m.get("mem", np.nan)):
            print(f"    {kind:<12}{'—':>7}{'—':>9}{'—':>11}{'—':>10}{'—':>7}  測れず"); continue
        c, cp = _cut(base, m), _cut(base, pl)
        v = ("**説明する**" if np.isfinite(cp) and c > cp + BEAT else
             "プラセボと同程度＝説明せず" if np.isfinite(cp) else "プラセボ無し＝判定保留")
        print(f"    {kind:<12}{m['acf1']:>7.2f}{m['mem']:>9.4f}{c:>+10.0%}"
              f"{(f'{cp:+.0%}' if np.isfinite(cp) else '—'):>10}{m['r2']:>7.2f}  {v}")


# ---------- 合成 -----------------------------------------------------------------
def make_synth(kind, n=800, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(n)
    seas = np.sin(2 * np.pi * t / 365.25)
    wT = np.zeros(n); wW = np.zeros(n)
    for i in range(1, n):
        wT[i] = 0.85 * wT[i - 1] + rng.normal(0, 1)
        wW[i] = 0.90 * wW[i - 1] + rng.normal(0, 1)
    T = 12 + 10 * seas + 1.2 * wT
    W = np.clip(0.25 + 0.05 * seas + 0.02 * wW, 0.03, 0.6)
    # 降水：3割の日に降る（指数分布）
    P = np.where(rng.random(n) < 0.22, rng.exponential(9.0, n), 0.0)   # 乾燥期間を長めに
    PA = 1000 + 8 * np.sin(2 * np.pi * t / 11.0) + 4 * rng.standard_normal(n)
    ly = 0.05 * T + 1.5 * W
    if kind == "hidden":
        phi = np.exp(-1.0 / 4.0)
        h = np.zeros(n)
        for i in range(1, n):
            h[i] = phi * h[i - 1] + rng.normal(0, 1)
        ly = ly + 0.30 * h / max(np.std(h), 1e-9)
    elif kind == "birch":
        ds = dryspell(P)
        pulse = np.zeros(n)
        for i in range(n):
            if P[i] > RAIN_MM:
                # **陽性対照として意味を持つ強さにする**。弱いと基準線の残差 ACF1 が
                # 0.2 台にしかならず、「取り除くべきメモリがほとんど無い」状態になって
                # 検定の合否を判定できない（第2版で気づいた）。
                pulse[i] += 0.30 * min(ds[i - 1] if i else 0.0, 20.0)
        # 降雨後に指数減衰するパルス（e-fold 3日）
        dec = np.exp(-1.0 / 3.0)
        acc = np.zeros(n)
        for i in range(1, n):
            acc[i] = dec * acc[i - 1] + pulse[i]
        ly = ly + acc
    else:                                   # pump（気圧由来の見かけの変動）
        d1 = np.concatenate([[0.0], np.diff(PA)])
        s = pd.Series(d1).rolling(4, min_periods=1).mean().to_numpy()
        ly = ly - 0.06 * s
    ly = ly + rng.normal(0, 0.10, n)
    return np.exp(ly), T, W, P, PA


def run_synth():
    print("  ── 合成：隠れ過程・Birch・気圧ポンピングを見分けられるか ──")
    for kind, lab in [("hidden", "**隠れ過程**（気象と無関係・e-fold 4日）"),
                      ("birch", "**Birch**（乾燥期間に応じた降雨後パルス・e-fold 3日）"),
                      ("pump", "**気圧ポンピング**（ΔPA に比例した見かけの変動）")]:
        y, T, W, P, PA = make_synth(kind)
        r = analyze(y, T, W, P, PA)
        print(f"  ━ {lab} ━")
        if r is None:
            print("    測れず\n"); continue
        show(r); print()
    print("  → 期待：隠れ過程は**どれを足しても落ちない**。")
    print("     Birch は**＋Birch で落ちる**。気圧ポンピングは**＋気圧で落ちる**。\n")


# ---------- 実データ -------------------------------------------------------------
def _pressure_from_csv(code, data_dir):
    """気圧を**CSV から直接読む**（旗77 第1版の欠陥＝自分の道具の欠陥16件目）。

    第1版は `load_raw_all` の返り値に `PA_F` を探した。だが `load_raw_all` は
    **変数マップに写像された列（RK_VARS）だけ**を返すので、**マップに無い `PA_F` は
    そこで捨てられている**。「生列から直接読む」とコメントしながら、
    **実際には捨てられた後のフレームを見ていた**＝全サイトで「気圧あり 0 日」になっていた。
    旗76 の下調べでは `PA_F` が確かに在ると確認できていたので、**食い違いで気づけた**。
    """
    root = Path(data_dir)
    csvs = sorted([p for p in root.rglob(f"*{code}*")
                   if p.is_file() and p.suffix.lower() == ".csv"
                   and "ALLVARS_HH" in p.name and "__MACOSX" not in p.parts],
                  key=lambda p: p.stat().st_size, reverse=True)
    for f in csvs[:2]:
        try:
            head = pd.read_csv(f, nrows=1)
        except Exception:
            continue
        col = next((c for c in ("PA_F", "PA_1_1_1") if c in head.columns), None)
        tcol = next((c for c in ("TIMESTAMP_START", "TIMESTAMP_END") if c in head.columns), None)
        if col is None or tcol is None:
            continue
        try:
            d = pd.read_csv(f, usecols=[tcol, col])
        except Exception:
            continue
        v = pd.to_numeric(d[col], errors="coerce")
        v[v <= -9000] = np.nan                      # JapanFLUX の欠測コード
        idx = pd.to_datetime(d[tcol].astype("Int64").astype(str),
                             format="%Y%m%d%H%M", errors="coerce")
        ser = pd.Series(v.to_numpy(), index=idx).dropna()
        if ser.empty:
            continue
        return ser.groupby(ser.index.normalize()).mean(), col
    return None, None


def tower_series(code, data_dir):
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    raw = load_raw_all(get_site(code), AnalysisConfig())
    out = {}
    if "P" in raw.columns:
        out["P"] = raw["P"].groupby(raw.index.normalize()).sum()
    pa, col = _pressure_from_csv(code, data_dir)
    if pa is not None:
        out["PA"], out["PA_col"] = pa, col
    return out


def run_real(cosore_dir, data_dir):
    root = Path(cosore_dir)
    print("  ── 実データ（同一地点4組）──")
    print("     基準線は旗74 の最も豊かな段階（T×W テンソルビン）。そこにタワー由来の列を足す。")
    print("     プラセボは**タワー側の生列を位相ずらししてから派生量を作り直す**"
          "（4 通りのうち最も下げた回）。\n")
    for code, ds, km in PAIRS:
        f = root / "datasets" / f"data_{ds}.csv"
        print(f"  ━ {code} ↔ {ds}（{km:.2f} km）━")
        if not f.exists():
            print("    チャンバーのデータが無い\n"); continue
        try:
            df, st, sm = load_cosore(f, None)
        except Exception as e:
            print(f"    チャンバー読み込み失敗 {type(e).__name__}: {str(e)[:90]}\n"); continue
        if df is None or "Tsoil" not in df or "Rs" not in df:
            print("    チャンバー側に Rs / Tsoil が無い\n"); continue
        cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
        daily = df[cols].groupby(df.index.normalize()).mean()
        daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
        try:
            tw = tower_series(code, data_dir)
        except Exception as e:
            print(f"    タワー読み込み失敗 {type(e).__name__}: {str(e)[:90]}\n"); continue
        P = tw["P"].reindex(daily.index).to_numpy() if "P" in tw else None
        PA = tw["PA"].reindex(daily.index).to_numpy() if "PA" in tw else None
        nP = int(np.isfinite(P).sum()) if P is not None else 0
        nPA = int(np.isfinite(PA).sum()) if PA is not None else 0
        print(f"    日数 {len(daily)}／降水あり {nP} 日／気圧あり {nPA} 日"
              f"（列 {tw.get('PA_col', '—')}）")
        if P is not None and nP:
            wet = float(np.nanmean(P > RAIN_MM))
            print(f"    >1mm の日 {wet:.0%}・最長乾燥期間 {np.nanmax(dryspell(P)):.0f} 日")
        y = daily["Rs"].to_numpy(); T = daily["Tsoil"].to_numpy()
        W = daily["SM"].to_numpy() if "SM" in daily else None
        r = analyze(y, T, W, P if nP > 60 else None, PA if nPA > 60 else None)
        if r is None:
            print("    基準線が測れない\n"); continue
        show(r)
        print()
    print("  === 読み方 ===")
    print("  **＋降水・＋Birch がプラセボと同程度なら、旗45/75 の最大の限定が外れる**")
    print("  ——『降水が無いから言えない』ではなく『**本物の降水でも説明しない**』になる。")
    print("  **＋気圧で落ちるなら、メモリの一部は物理（pressure pumping）**であり、")
    print("  旗59 が未決着にした『生物か物理か』に、**新規観測なしで一つ答えが出る**。")
    print("  留保：")
    print("   ・タワーの降水・気圧は**タワー位置**の観測で、チャンバー直上ではない（0.00〜0.69 km）。")
    print("   ・`PA_F` は**穴埋め済み**（`PA_ERA5` による再解析補完を含みうる）。")
    print("     気圧は再解析でよく拘束される量だが、**実測のみではない**ことは記しておく。")
    print("   ・pressure pumping は**チャンバーの設計（開放系か閉鎖系か）**に強く依存する。")
    print("     ここで落ちなくても『物理機構が無い』ではなく『**この形では検出されない**』である。")
    print("   ・4 組しかない＝**集計ではなく事例**として読むこと。")


def main():
    p = argparse.ArgumentParser(description="降水と気圧でメモリは説明できるか")
    p.add_argument("--cosore-dir")
    p.add_argument("--data-dir", default="/mnt/hdd/JAPANFLUX")
    p.add_argument("--synth", action="store_true")
    a = p.parse_args()
    print("=== 旗77：降水と気圧でメモリは説明できるか（六つ目と七つ目の対抗仮説）===")
    print("  旗76 で**4組すべて・全期間**で降水が重なり、**気圧 PA_F も在る**と確認済み。\n")
    if a.synth or not a.cosore_dir:
        run_synth()
        if not a.cosore_dir:
            return
    run_real(a.cosore_dir, a.data_dir)


if __name__ == "__main__":
    main()
