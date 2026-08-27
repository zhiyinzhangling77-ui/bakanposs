"""旗75：**過去の気象は、どんな形でもメモリを説明できないか**——題名そのものの検定。

旗74 で A-1（チャンバー呼吸の多日メモリ）は**雑音でも誤特定でも作れない**と分かった。
だが旗74 の留保に書いたとおり、**ビンは「同時刻の T と W」しか使っていない**。
＝「メモリ」が**過去の気象の関数**である可能性は、そこでは潰していない。

旗45/54 が同じ問いを扱ったが、モデル族は**線形＋Birch 非線形**に限られていた。
いまは**非パラメトリックな基底・外挿残差・複数プラセボ**が揃っている。**やり直す価値がある**。

## 何を足すか（段階的に）

  1. **同時刻のみ**：T×W のテンソルビン（旗74 の最も豊かな段階＝これが基準線）
  2. **＋遅れ 1–3 日**：T・W の 1〜3 日前（一次と二次）
  3. **＋積算 7/30 日**：過去 7 日・30 日の移動平均（一次と二次）＝「先行湿潤」「熱履歴」の一般形
  4. **＋両方**

## プラセボが要る理由

**説明変数を足せば、それだけで残差の自己相関は多少下がる**（自由度が増えるため）。
だから**同じ次元・同じ作り方で、時間の対応だけを壊した変数**を併走させる。
本ツールは**複数の位相ずらし**（100・180・260・500 日）を試し、**最も ACF1 を下げた回**を基準線に採る
＝**プラセボに有利な、厳しい比較**にする。

  ・**実データの低下がプラセボの低下を明確に上回る** → 過去の気象が説明している。
  ・**同程度** → **どんな形の過去の気象でも説明できない**＝メモリは観測の外側にある。

## 合成での判別試験

  ・`hidden`：メモリを**気象と無関係な隠れ過程**から作る → 遅れを足しても落ちないはず。
  ・`lagged`：メモリを**過去7日の水分の関数**として作る → 遅れ・積算を足すと落ちるはず。

    python research/lagged_drivers_step75.py --synth
    python research/lagged_drivers_step75.py --cosore-dir /mnt/hdd/cosore-0.7.0 --igbp forest
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
from model_richness_step74 import design, residuals, star

SETS = ("同時刻のみ", "＋遅れ1-3日", "＋積算7/30日", "＋両方")
SHIFTS = (100, 180, 260, 500)          # プラセボの位相ずらし（日）
LAGS = (1, 2, 3)
WINS = (7, 30)


def _lagmat(x, lags):
    """遅れ変数（一次と二次）。先頭は NaN になる。"""
    out = []
    for L in lags:
        v = np.full(len(x), np.nan)
        v[L:] = x[:-L]
        out += [v, v ** 2]
    return out


def _rollmat(x, wins):
    """過去 w 日の移動平均（一次と二次）。**当日を含めない**（過去だけ）。"""
    s = pd.Series(x)
    out = []
    for w in wins:
        v = s.shift(1).rolling(w, min_periods=max(3, w // 3)).mean().to_numpy()
        out += [v, v ** 2]
    return out


def extra_cols(T, W, kind):
    """段階ごとに足す列。``kind`` は SETS のいずれか。"""
    cols = []
    src = [T] + ([W] if W is not None else [])
    if kind in ("＋遅れ1-3日", "＋両方"):
        for x in src:
            cols += _lagmat(x, LAGS)
    if kind in ("＋積算7/30日", "＋両方"):
        for x in src:
            cols += _rollmat(x, WINS)
    return cols


def measure(y, T, W, kind, shift=None):
    """テンソルビン（同時刻）に段階ごとの列を足して、**外挿残差**の ACF1 を測る。

    ``shift`` を渡すと、**足す列だけ**を位相ずらしした T/W から作る＝プラセボ。
    同時刻のビンは動かさないので、**増えた自由度は同じで、時間の対応だけが壊れる**。
    """
    base = design("テンソルビン", T, W)
    if base is None:
        return None
    if kind != "同時刻のみ":
        if shift is None:
            Ts, Ws = T, (W if W is not None else None)
        else:
            Ts = np.roll(T, shift)
            Ws = np.roll(W, shift) if W is not None else None
        ex = extra_cols(Ts, Ws, kind)
        if ex:
            base = np.hstack([base, np.column_stack(ex)])
    with np.errstate(divide="ignore", invalid="ignore"):
        ly = np.log(np.where(y > 0, y, np.nan))
    res = residuals(base, ly, True)
    if res is None or np.isfinite(res).sum() < 40:
        return None
    ss = np.nansum((ly - np.nanmean(ly)) ** 2)
    r2 = float(1 - np.nansum(res ** 2) / ss) if ss > 0 else np.nan
    return {"r2": r2, "acf1": _acf_gap(res, 1), "efold": _efold_gap(res),
            "n": int(np.isfinite(res).sum())}


def placebo_best(y, T, W, kind):
    """複数の位相ずらしを試し、**最も ACF1 を下げた回**（プラセボに有利な側）を返す。"""
    best = None
    for sh in SHIFTS:
        if sh >= len(T):
            continue
        m = measure(y, T, W, kind, shift=sh)
        if m is None or not np.isfinite(m["acf1"]):
            continue
        if best is None or m["acf1"] < best["acf1"]:
            best = m
    return best


def analyze(y, T, W):
    base = measure(y, T, W, "同時刻のみ")
    if base is None or not np.isfinite(base["acf1"]):
        return None
    rows = {"同時刻のみ": (base, None)}
    for kind in SETS[1:]:
        rows[kind] = (measure(y, T, W, kind), placebo_best(y, T, W, kind))
    return rows


def _d(base, m):
    if m is None or not np.isfinite(m["acf1"]):
        return np.nan
    return base["acf1"] - m["acf1"]


# ---------- 合成 -----------------------------------------------------------------
def make_synth(kind, n=900, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(n)
    seas = np.sin(2 * np.pi * t / 365.25)
    wT = np.zeros(n); wW = np.zeros(n)
    for i in range(1, n):
        wT[i] = 0.85 * wT[i - 1] + rng.normal(0, 1)
        wW[i] = 0.90 * wW[i - 1] + rng.normal(0, 1)
    T = 12 + 10 * seas + 1.2 * wT
    W = np.clip(0.25 + 0.05 * seas + 0.02 * wW, 0.03, 0.6)
    ly = 0.05 * T + 1.5 * W
    if kind == "hidden":
        # **気象と無関係な隠れ過程**（e-fold 4日）
        phi = np.exp(-1.0 / 4.0)
        h = np.zeros(n)
        for i in range(1, n):
            h[i] = phi * h[i - 1] + rng.normal(0, 1)
        ly = ly + 0.30 * h / max(np.std(h), 1e-9)
    else:
        # **過去7日の水分の関数**（＝観測から作れるメモリ）
        past = pd.Series(W).shift(1).rolling(7, min_periods=3).mean().to_numpy()
        past = np.where(np.isfinite(past), past, np.nanmean(W))
        ly = ly + 6.0 * (past - np.nanmean(past))
    ly = ly + rng.normal(0, 0.10, n)
    return np.exp(ly), T, W


def show(rows):
    base = rows["同時刻のみ"][0]
    print(f"    {'足した列':<16}{'ACF1':>8}{'低下':>8}{'プラセボ':>10}{'低下':>8}{'R²':>8}  判定")
    print(f"    {'同時刻のみ（基準）':<16}{base['acf1']:>8.2f}{'—':>8}{'—':>10}{'—':>8}"
          f"{base['r2']:>8.2f}")
    for kind in SETS[1:]:
        m, pl = rows[kind]
        if m is None or not np.isfinite(m["acf1"]):
            print(f"    {kind:<16}{'—':>8}{'—':>8}{'—':>10}{'—':>8}{'—':>8}  測れず"); continue
        d = _d(base, m); dp = _d(base, pl)
        v = ("**過去の気象が説明**" if np.isfinite(dp) and d > dp + 0.10 else
             "プラセボと同程度＝説明せず" if np.isfinite(dp) else "プラセボ無し＝判定保留")
        print(f"    {kind:<16}{m['acf1']:>8.2f}{d:>+8.2f}"
              f"{(pl['acf1'] if pl else np.nan):>10.2f}{dp:>+8.2f}{m['r2']:>8.2f}  {v}")


def run_synth():
    print("  ── 合成：隠れ過程のメモリと、過去の気象から作れるメモリを見分けられるか ──")
    for kind, lab in [("hidden", "**隠れ過程**（気象と無関係・e-fold 4日）"),
                      ("lagged", "**過去7日の水分の関数**（＝観測から作れる）")]:
        y, T, W = make_synth(kind)
        r = analyze(y, T, W)
        print(f"  ━ {lab} ━")
        if r is None:
            print("    測れず\n"); continue
        show(r); print()
    print("  → 期待：隠れ過程は**どれを足しても落ちない**（プラセボと同程度）。")
    print("     過去7日の水分から作った場合は**積算を足すと大きく落ち、プラセボを明確に上回る**。\n")


# ---------- 実データ -------------------------------------------------------------
def run_real(cosore_dir, igbp, months):
    root = Path(cosore_dir); desc = pd.read_csv(root / "description.csv")
    print(f"  ── 実データ（{igbp or '全'}）──")
    print("     基準線は旗74 の最も豊かな段階（テンソルビン・同時刻のみ）。そこに遅れ・積算を足す。")
    print("     プラセボは**同じ次元で時間の対応だけ壊した列**（4 通りの位相ずらしのうち"
          "**最も下げた回**を採る）。\n")
    hdr = "".join(f"{k:>12}" for k in SETS)
    print(f"  {'dataset':<30}{hdr}   プラセボ最良")
    tally = {k: [] for k in SETS[1:]}
    npl = {k: [] for k in SETS[1:]}
    beat = {k: 0 for k in SETS[1:]}
    nsite = 0
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        if igbp and igbp.lower() not in ig.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f, months)
        except Exception:
            continue
        if df is None or "Tsoil" not in df or "Rs" not in df:
            continue
        cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
        daily = df[cols].groupby(df.index.normalize()).mean()
        daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
        y = daily["Rs"].to_numpy(); T = daily["Tsoil"].to_numpy()
        W = daily["SM"].to_numpy() if "SM" in daily else None
        r = analyze(y, T, W)
        if r is None:
            continue
        nsite += 1
        base = r["同時刻のみ"][0]
        cells = f"{base['acf1']:>12.2f}"
        pls = []
        for kind in SETS[1:]:
            m, pl = r[kind]
            cells += (f"{m['acf1']:>12.2f}" if m and np.isfinite(m["acf1"]) else f"{'—':>12}")
            dd, dp = _d(base, m), _d(base, pl)
            if np.isfinite(dd):
                tally[kind].append(dd)
            if np.isfinite(dp):
                npl[kind].append(dp)
            if np.isfinite(dd) and np.isfinite(dp) and dd > dp + 0.10:
                beat[kind] += 1
            pls.append(f"{dp:+.2f}" if np.isfinite(dp) else "—")
        print(f"  {ds:<30}{cells}   {' '.join(pls)}")
    print(f"\n  === まとめ（{nsite} サイト）===")
    print(f"  {'足した列':<16}{'実データの低下':>16}{'プラセボの低下':>16}"
          f"{'プラセボを 0.10 超えた数':>24}")
    for kind in SETS[1:]:
        a = np.asarray(tally[kind], float); b = np.asarray(npl[kind], float)
        if len(a) == 0:
            print(f"  {kind:<16}{'—':>16}{'—':>16}{'—':>24}"); continue
        print(f"  {kind:<16}{np.median(a):>+15.3f}"
              f"{(np.median(b) if len(b) else np.nan):>+15.3f}"
              f"{beat[kind]:>20}/{nsite}")
    print("\n  読み方：")
    print("   ・**実データの低下がプラセボと同程度なら、どんな形の過去の気象でも説明できない**")
    print("     ＝メモリは**観測の外側**にある＝この研究の題名そのものが支持される。")
    print("   ・**プラセボを明確に上回るサイトが多いなら、そのぶんは過去の気象で説明できる**")
    print("     ＝『観測の外側』の範囲を狭め、どの遅れ・どの窓が効いたかを報告し直す。")
    print("  留保：")
    print("   ・遅れは 1–3 日、窓は 7・30 日に限った。**もっと長い記憶は試していない**。")
    print("   ・**降水は COSORE に無い**（旗45 と同じ限定）。ここでの『過去の気象』は")
    print("     **土壌温度と土壌水分の過去**であって、降水イベントそのものではない。")
    print("   ・列を足すと外挿 R² が落ちるサイトがある（自由度の代償）。")
    print("     R² が大きく落ちた場合、ACF1 の比較も不安定になる。")


def main():
    p = argparse.ArgumentParser(description="過去の気象はメモリを説明するか")
    p.add_argument("--cosore-dir")
    p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    p.add_argument("--synth", action="store_true")
    a = p.parse_args()
    print("=== 旗75：過去の気象は、どんな形でもメモリを説明できないか ===")
    print("  旗74 は**同時刻の T・W だけ**を使った。ここでは**遅れと積算**を足して、")
    print("  **外挿残差**の ACF1 が**プラセボを超えて**落ちるかを見る。\n")
    if a.synth or not a.cosore_dir:
        run_synth()
        if not a.cosore_dir:
            return
    run_real(a.cosore_dir, a.igbp, a.month)


if __name__ == "__main__":
    main()
