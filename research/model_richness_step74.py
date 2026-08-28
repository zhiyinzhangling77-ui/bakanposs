"""旗74：**駆動モデルを豊かにしてもメモリは残るか**——中核主張 A-1 の決定的な検定。

旗73 で分かったこと：
  ・**雑音では A-1 を説明できない**（実データ由来の帰無で偽陽性率 0.0%）。
  ・**だがその帰無は当てはめモデル自身から作るので、誤特定に対して盲目**である（循環している）。
  ・合成では**誤特定があれば帰無 ACF1 が +0.59 まで上がる**＝**閾値 0.64 のすぐ下**。

＝**残る唯一の重大な対抗仮説は「駆動モデルの誤特定」**である。
「メモリ」に見えているものが、実は**地温・水分の効き方を取り違えた分の残りかす**かもしれない。

## やり方：モデルを段階的に豊かにして、ACF1 がどう動くかを見る

  1. **線形**：[T, W]（旗40 の形。旗52 でこれが自己相関残差を作ると判明済み）
  2. **非線形基底**：Lloyd-Taylor 項・二乗・交互作用・log W（**旗53/54 の現行**）
  3. **加法ビン**：T を10分位・W を6分位に切った**指示変数**（形を仮定しない・加法的）
  4. **テンソルビン**：T×W を 8×4 の**セル平均**（形も交互作用も仮定しない・完全に非パラメトリック）

**ACF1 が階段状に落ちるなら誤特定**（豊かにするほど説明され、残差から消える）。
**+0.8 付近で頭打ちなら本物**（どんな形の駆動関数でも説明できない成分が残る）。

## **外挿残差を使う**（ここが要点）

豊かなモデルは**メモリそのものを吸ってしまう**。日ごとにパラメータを増やせば残差はゼロにできる。
＝**内挿残差で見ると、豊かにするほど ACF1 が落ちる**のは**当たり前**であって、誤特定の証拠にならない。

そこで**時間ブロックの交差検証**（連続した5ブロックに分け、他の4つで当てはめて残りを予測）を使い、
**外挿残差**で ACF1 を測る。過学習は**予測を悪くする**が、**自己相関は作らない**。
＝**外挿で ACF1 が落ちないなら、それは吸えなかった本物の成分**である。

内挿・外挿の**両方を出す**ので、読者はその差そのものを見られる。

## 合成での判別試験

  ・`memory`：駆動は基底で表せる形＋**本物のメモリ**を植える
    → **どの豊かさでも ACF1 は落ちないはず**。
  ・`misspec`：**メモリは無い**が、駆動の効き方を**基底で表せない形**にする
    → **線形・非線形基底では高く、ビンに移すと落ちるはず**。

これが分かれれば、実データの結果を読む資格がある。

    python research/model_richness_step74.py --synth
    python research/model_richness_step74.py --cosore-dir /mnt/hdd/cosore-0.7.0 --igbp forest
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
from memory_attribution_flex_step54 import flex_basis, ACF_THR, EFOLD_MAX

LEVELS = ("線形", "非線形基底", "加法ビン", "テンソルビン")
NFOLD = 5


def _bins(x, nb):
    """分位でビン分けした**指示変数**（形を仮定しない）。"""
    ok = np.isfinite(x)
    if ok.sum() < nb * 5:
        return None
    edges = np.nanpercentile(x[ok], np.linspace(0, 100, nb + 1)[1:-1])
    lab = np.digitize(x, edges)
    M = np.zeros((len(x), nb))
    M[np.arange(len(x)), np.clip(lab, 0, nb - 1)] = 1.0
    M[~ok] = np.nan
    return M


def design(level, T, W):
    """豊かさの段階ごとの計画行列（定数項は各関数の末尾で足す）。"""
    n = len(T)
    one = np.ones((n, 1))
    if level == "線形":
        cols = [T.reshape(-1, 1)] + ([W.reshape(-1, 1)] if W is not None else [])
        return np.hstack(cols + [one])
    if level == "非線形基底":
        return np.column_stack(flex_basis(T, W) + [np.ones(n)])
    if level == "加法ビン":
        mt = _bins(T, 10)
        if mt is None:
            return None
        parts = [mt]
        if W is not None:
            mw = _bins(W, 6)
            if mw is not None:
                parts.append(mw)
        return np.hstack(parts + [one])
    if level == "テンソルビン":
        mt = _bins(T, 8)
        if mt is None:
            return None
        if W is None:
            return np.hstack([mt, one])
        mw = _bins(W, 4)
        if mw is None:
            return np.hstack([mt, one])
        # 8×4 のセル指示（交互作用も形も仮定しない）
        cells = np.einsum("ij,ik->ijk", mt, mw).reshape(len(T), -1)
        return np.hstack([cells, one])
    raise ValueError(level)


def _fit_predict(X, y, rows_fit, rows_pred):
    ok = rows_fit & np.isfinite(y) & np.isfinite(X).all(axis=1)
    if ok.sum() < X.shape[1] + 10:
        return None
    coef = np.linalg.lstsq(X[ok], y[ok], rcond=None)[0]
    out = np.full(len(y), np.nan)
    pr = rows_pred & np.isfinite(X).all(axis=1)
    pred = X[pr] @ coef
    # **外挿の暴走を止める**（自分の道具の欠陥17件目）。
    # 列を増やすと、学習範囲の外で予測が桁違いの値を返すことがある。
    # 残差分散が爆発し、**記憶量（ACF1×分散）の比が発散**して
    # 削減率が 10^15 % のような無意味な数値になっていた（旗75 の実データで発覚）。
    # 学習した y の範囲を**その幅の半分だけ広げた区間**に予測を丸める＝
    # 「この当てはめが語れる範囲」を超えた外挿を、**無かったことにせず、頭打ちにする**。
    lo, hi = np.nanmin(y[ok]), np.nanmax(y[ok])
    pad = 0.5 * (hi - lo)
    out[pr] = np.clip(pred, lo - pad, hi + pad)
    return out


def residuals(X, y, oos):
    """``oos`` なら**時間ブロック交差検証の外挿残差**、でなければ内挿残差。"""
    n = len(y)
    if not oos:
        pred = _fit_predict(X, y, np.ones(n, bool), np.ones(n, bool))
        return None if pred is None else y - pred
    edges = np.linspace(0, n, NFOLD + 1).astype(int)
    res = np.full(n, np.nan)
    for k in range(NFOLD):
        test = np.zeros(n, bool); test[edges[k]:edges[k + 1]] = True
        pred = _fit_predict(X, y, ~test, test)
        if pred is None:
            continue
        r = y[test] - pred[test]
        if not np.isfinite(r).any():        # 有効な残差が無いブロック（警告を出さない）
            continue
        # **ブロックごとに中心化する**（旗74 第1版の欠陥＝14件目）。
        # 交差検証では各ブロックを**別々のデータで当てはめる**ため、
        # ブロック間に**系統的な水準差（段差）**が生じる。段差は長い時間尺度の成分であり、
        # **e-fold を大きく見せて「短メモリ」の判定を潰す**（★が不当に減る）。
        # ブロック長は全体の1/5（数百日）なので、**4〜7日のメモリには影響しない**。
        m = np.nanmean(r)
        res[test] = r - (m if np.isfinite(m) else 0.0)
    return res


def measure(y, T, W, level, oos):
    X = design(level, T, W)
    if X is None:
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        ly = np.log(np.where(y > 0, y, np.nan))
    res = residuals(X, ly, oos)
    if res is None or np.isfinite(res).sum() < 40:
        return None
    ss = np.nansum((ly - np.nanmean(ly)) ** 2)
    r2 = float(1 - np.nansum(res ** 2) / ss) if ss > 0 else np.nan
    return {"r2": r2, "acf1": _acf_gap(res, 1), "efold": _efold_gap(res),
            "n": int(np.isfinite(res).sum())}


def star(m):
    if m is None or not np.isfinite(m.get("r2", np.nan)) or m["r2"] < 0.3:
        return None
    if not (np.isfinite(m["acf1"]) and np.isfinite(m["efold"])):
        return None
    return bool(m["acf1"] >= ACF_THR and m["efold"] <= EFOLD_MAX)


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
    if kind == "memory":
        # 基底で表せる駆動 ＋ **本物のメモリ**
        ly = 0.05 * T + 1.5 * W
        phi = np.exp(-1.0 / 4.0)
        h = np.zeros(n)
        for i in range(1, n):
            h[i] = phi * h[i - 1] + rng.normal(0, 1)
        ly = ly + 0.30 * h / max(np.std(h), 1e-9)
    else:
        # **メモリは無い**が、駆動の効き方が基底で表せない（階段＋鋭い折れ）
        step = 0.5 * (T > 14) + 0.5 * (T > 20)
        kink = np.where(W < 0.22, 6.0 * (0.22 - W), 0.0)
        ly = 0.02 * T + step - kink
    ly = ly + rng.normal(0, 0.10, n)
    return np.exp(ly), T, W


def show_row(lab, ms):
    cells = "".join(
        (f"{m['acf1']:>7.2f}" if m and np.isfinite(m['acf1']) else f"{'—':>7}") for m in ms)
    r2s = "".join(
        (f"{m['r2']:>7.2f}" if m and np.isfinite(m['r2']) else f"{'—':>7}") for m in ms)
    print(f"    {lab:<22}ACF1{cells}      R²{r2s}")


def run_synth():
    print("  ── 合成：本物のメモリと、誤特定だけの場合を見分けられるか ──")
    hdr = "".join(f"{l:>7}" for l in LEVELS)
    for kind, lab in [("memory", "**本物のメモリ**（e-fold 4日・駆動は基底で表せる）"),
                      ("misspec", "**誤特定だけ**（メモリ無し・駆動が階段＋折れ）")]:
        y, T, W = make_synth(kind)
        print(f"  ━ {lab} ━")
        print(f"    {'':<22}    {hdr}")
        show_row("内挿残差", [measure(y, T, W, L, False) for L in LEVELS])
        show_row("**外挿残差（CV）**", [measure(y, T, W, L, True) for L in LEVELS])
        print()
    print("  → 期待：**本物のメモリはどの豊かさでも ACF1 が落ちない**。")
    print("     **誤特定だけの場合は、ビンに移すと外挿でも ACF1 が落ちる**。")
    print("     内挿だけを見ると両方落ちうる（豊かなモデルがメモリを吸うため）＝**外挿が要**。\n")


# ---------- 実データ -------------------------------------------------------------
def run_real(cosore_dir, igbp, months):
    root = Path(cosore_dir); desc = pd.read_csv(root / "description.csv")
    print(f"  ── 実データ（{igbp or '全'}）：外挿残差（{NFOLD}ブロック交差検証）の ACF1 ──")
    print("     豊かにしても落ちないなら、それは**どんな形の駆動関数でも説明できない成分**である。\n")
    hdr = "".join(f"{l:>9}" for l in LEVELS)
    print(f"  {'dataset':<30}{'N':>5}{hdr}   ★の推移")
    counts = {L: 0 for L in LEVELS}
    judged = 0
    drops = []
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
        ms = [measure(y, T, W, L, True) for L in LEVELS]
        sts = [star(m) for m in ms]
        if all(s is None for s in sts):
            continue
        judged += 1
        for L, s in zip(LEVELS, sts):
            counts[L] += int(bool(s))
        cells = "".join((f"{m['acf1']:>9.2f}" if m and np.isfinite(m['acf1'])
                         else f"{'—':>9}") for m in ms)
        trail = "".join("★" if s else ("·" if s is False else "?") for s in sts)
        print(f"  {ds:<30}{int(np.isfinite(y).sum()):>5}{cells}   {trail}")
        a0 = ms[1]["acf1"] if ms[1] and np.isfinite(ms[1]["acf1"]) else np.nan
        a3 = ms[3]["acf1"] if ms[3] and np.isfinite(ms[3]["acf1"]) else np.nan
        if np.isfinite(a0) and np.isfinite(a3):
            drops.append(a0 - a3)
    print(f"\n  === まとめ ===")
    print(f"  判定できたサイト {judged}")
    print("  **★の件数（外挿残差・規則は旗53/54 と同じ）**")
    for L in LEVELS:
        print(f"    {L:<12} {counts[L]:>3} 件")
    if drops:
        dr = np.asarray(drops)
        print(f"  **非線形基底 → テンソルビン の ACF1 低下**："
              f"中央 {np.median(dr):+.3f}（四分位 {np.percentile(dr,25):+.3f}〜{np.percentile(dr,75):+.3f}）")
    print("\n  読み方：")
    print("   ・**★の件数が豊かさによらずほぼ一定なら、A-1 は誤特定では説明できない**")
    print("     ＝中核主張は最も強い形で残る。")
    print("   ・**ビンに移して件数が大きく減るなら、その分は誤特定だった**")
    print("     ＝A-1 の件数を減らし、**どの形の駆動関数が効いていたか**を報告し直す。")
    print("   ・低下の中央値が 0 付近なら、非線形基底は**すでに十分豊か**だったことになる。")
    print("  留保：")
    print("   ・ビンは**同時刻の T と W だけ**を使う。**遅れ**（過去の T/W）は入れていない。")
    print("     ＝『**過去の気象で説明できるか**』は別の問い（旗45/54 で扱った）。")
    print("   ・テンソルビンはセルが 32 個あり、**日数の少ないサイトでは外挿が不安定**になる。")
    print("     R² が大きく落ちているサイトはそう読むこと。")
    print("   ・時間ブロック交差検証は**ブロック境界で残差がつながる**。5分割なので")
    print("     境界は4点であり ACF1 への影響は小さいが、ゼロではない。")


def main():
    p = argparse.ArgumentParser(description="駆動モデルを豊かにしてもメモリは残るか")
    p.add_argument("--cosore-dir")
    p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    p.add_argument("--synth", action="store_true")
    a = p.parse_args()
    print("=== 旗74：駆動モデルを豊かにしてもメモリは残るか（A-1 の決定的検定）===")
    print("  旗73 で『雑音では説明できない』が確定し、**残る対抗仮説は誤特定**だけになった。")
    print("  線形→非線形基底→加法ビン→テンソルビンと豊かにし、**外挿残差**の ACF1 を見る。\n")
    if a.synth or not a.cosore_dir:
        run_synth()
        if not a.cosore_dir:
            return
    run_real(a.cosore_dir, a.igbp, a.month)


if __name__ == "__main__":
    main()
