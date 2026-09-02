"""旗102：**★という判定は、駆動センサの正体に依存していないか**（事前登録 step102）。

**新しい仮説ではない。既存の主張（A-1：約 1/3 のチャンバーに多日記憶）の頑健性検査である。**

旗101 で二つ分かった：
  ・**駆動に使われた層は 1 cm〜30 cm に散らばっている**（「5 cm に最も近い層」を選ぶだけ）
  ・**土壌温度列が無い地点では `CSR_TAIR`（気温）に落ちる**——**一度も報告していない**

**当初の見立て**：土壌は気温を減衰・遅延させる → **気温で引いた地点ほど★になりやすい**。
**合成でこれは否定された**（追補）——**気温は「鈍った土壌温度」ではなく、
白色雑音を余分に持つ**。**その雑音が残差を薄め、★を減らす。**
**どちらの向きも起こりうるので、判定規則は両側のまま走らせる。**

**事前登録 step102 で固定済み**：
  ・★の判定は**旗74 と同一**（テンソルビン＋外挿・ACF1 ≥0.64・e-fold ≤7 日）／**窓は全期間**
  ・**軸A：土壌温度 vs 気温**／**軸B：浅い（≤5cm） vs 深い（>5cm）**
  ・**気温群の★率が土壌群の 2 倍以上 → A-1 は水増しされている**
  ・**逆に低ければ「予想と逆」とそのまま書く**（**合成では気温群が★を失った**——追補参照）
  ・**どちらかの群が 3 本未満 → その軸は判定しない**
  ・**本数と地点数を必ず併記する**（SZUTU 5 本・KAYE 8 本・SIHI 3 本は**それぞれ 1 地点**）

    python research/driver_sensor_step102.py                      # 合成で検証（既定）
    python research/driver_sensor_step102.py --real --cosore-dir /mnt/hdd/cosore-0.7.0
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
from model_richness_step74 import measure, star
from runlog import tee_stdout      # **出力を最初からファイルに残す**（旗110 の反省）

MIN_DAYS, MIN_YEARS = 60, 3
SHALLOW_CM = 5.0          # 軸B の境（**≤5 cm を浅いとする**・事前登録で固定）
RATIO_HI = 2.0            # 「水増し」の閾値（事前登録で固定）


def classify(col):
    """駆動列の名前から **(種別, 深さ)** を返す。`CSR_T5`→('土壌', 5.0)／`CSR_TAIR`→('気温', nan)。"""
    if not col:
        return "無し", np.nan
    m = re.fullmatch(r"CSR_T(\d+\.?\d*)", col)
    if m:
        return "土壌", float(m.group(1))
    if col.startswith("CSR_TAIR"):
        return "気温", np.nan
    m = re.fullmatch(r"CSR_SM(\d+\.?\d*)", col)
    if m:
        return "水分", float(m.group(1))
    return "その他", np.nan


def verdict(daily):
    """★かどうか（**旗74 と同一の物差し・窓は全期間**）。"""
    # **道具の欠陥 #34（旗102 の実行で判明）**：第1版はこの二つを同じ文言で報告していた。
    # **「駆動列が無い」と「日数が足りない」は別の理由である。** 分けて書く。
    if "Tsoil" not in daily:
        return None, "土壌温度の列が無い"
    if len(daily) < MIN_DAYS:
        return None, f"{len(daily)} 日（下限未満）"
    if daily.index.year.nunique() < MIN_YEARS:
        return None, f"{daily.index.year.nunique()} 暦年（下限未満）"
    y = daily["Rs"].to_numpy()
    T = daily["Tsoil"].to_numpy()
    W = daily["SM"].to_numpy() if "SM" in daily else None
    if np.isfinite(T).sum() < MIN_DAYS:
        return None, "土壌温度の有限値が足りない"
    m = measure(y, T, W, "テンソルビン", True)
    if not m:
        return None, "残差を作れない"
    return (star(m), m), None


def rate(rows, key, val):
    sub = [r for r in rows if r[key] == val]
    n = len(sub); k = sum(r["star"] for r in sub)
    return n, k, (k / n if n else np.nan), sub


def report_axis(rows, name, key, ga, gb):
    print(f"\n  === {name} ===")
    out = {}
    for g in (ga, gb):
        n, k, p, sub = rate(rows, key, g)
        sites = len({r["site"] for r in sub})
        out[g] = (n, k, p)
        print(f"    {g:<10}{n:>3} 本（**{sites} 地点**）／★ {k:>3} 本／"
              f"**★率 {p:.0%}**" if n else f"    {g:<10}  0 本")
        for r in sorted(sub, key=lambda r: -r["star"]):
            print(f"        {r['ds'][:34]:<36}{r['col']:<12}"
                  f"{'**★**' if r['star'] else '  －'}  "
                  f"ACF1 {r['acf1']:+.2f}／e-fold {r['efold']:.0f}日／{r['days']:,}日/{r['yrs']}年")
    na, ka, pa = out[ga]; nb, kb, pb = out[gb]
    if na < 3 or nb < 3:
        print(f"    → **判定しない**（{ga} {na} 本・{gb} {nb} 本——"
              f"**どちらかが 3 本未満**）")
        return
    if not (np.isfinite(pa) and np.isfinite(pb)) or pa == 0:
        print("    → **判定しない**（比を出せない）"); return
    ratio = pb / pa
    print(f"    **比（{gb} ÷ {ga}）＝ {ratio:.2f}**")
    if ratio >= RATIO_HI:
        print(f"    → **★{gb} の方が★率が {RATIO_HI:.0f} 倍以上高い**"
              f"——**A-1 は水増しされている**")
    elif ratio <= 1 / RATIO_HI:
        print(f"    → **○予想と逆**（{gb} の方が★率が低い）。**そのまま書く**")
    else:
        print(f"    → **○この軸では★率は変わらない**（比 {1/RATIO_HI:.1f}〜{RATIO_HI:.1f} の間）")


def synth(kind, seed=0, years=5, amp=0.30):
    """**真の駆動は土壌温度＋隠れた ~4 日駆動**。`airdriven` は**気温しか観測できない**。

    **第1版の予測（気温で引くと★が増える）は間違いだった**——追補に記録。
    **気温は土壌温度の「鈍った版」ではなく、白色雑音を余分に持つ。**
    **その雑音が残差に混ざり、ACF1 を薄めて★を減らす。**
    """
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2010-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Tair = 12 + 14 * np.sin(2 * np.pi * (doy - 100) / 365) + rng.normal(0, 4.0, len(idx))
    # **土壌＝気温の減衰つき遅延**（指数平滑＝熱慣性）。**これが真の駆動。**
    Tsoil = np.zeros(len(idx)); Tsoil[0] = Tair[0]
    for i in range(1, len(idx)):
        Tsoil[i] = 0.82 * Tsoil[i - 1] + 0.18 * Tair[i]
    W = np.clip(.25 + .05 * np.sin(2 * np.pi * (doy - 200) / 365)
                + rng.normal(0, .02, len(idx)), .02, .6)
    # **本物の ~4 日メモリ**（未観測駆動）——**これが★の源**
    hid = np.convolve(rng.normal(0, 1, len(idx)), np.ones(6) / 6, "same")
    hid /= hid.std()
    lrs = -1.0 + 0.07 * Tsoil + 1.5 * W + amp * hid + rng.normal(0, .03, len(idx))
    obsT = Tair if kind == "airdriven" else Tsoil
    return pd.DataFrame({"Rs": np.exp(lrs), "Tsoil": obsT, "SM": W}, index=idx)


def main():
    ap = argparse.ArgumentParser(description="旗102：★は駆動センサの正体に依存するか")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--cosore-dir", default="/mnt/hdd/cosore-0.7.0")
    a = ap.parse_args()

    tee_stdout("step102")
    print("=== 旗102：★という判定は、駆動センサの正体に依存していないか ===")
    print("  **新しい仮説ではない。A-1（約 1/3 に多日記憶）の頑健性検査である。**")
    print("  **当初の見立て（気温で引くと★が増える）は合成で否定された**（追補）——")
    print("  **合成では気温で引くと★を失った**。**判定規則は両側のまま走らせる。**")
    print(f"  ★は**旗74 と同一**（ACF1 ≥0.64・e-fold ≤7 日）／**窓は全期間**／"
          f"下限 {MIN_DAYS} 日・{MIN_YEARS} 暦年。")

    if not a.real:
        print("\n  【合成データで道具を検証する】")
        print("  **真の駆動は土壌温度＋隠れた ~4 日駆動**（＝★の源を仕込む）。")
        print("  **`soildriven`（土壌温度を観測できる）→ ★になるべき**")
        print("  **`airdriven`（気温しか観測できない）→ ★を失うべき**")
        print("  **追補どおり、第1版の予測（気温で★が増える）は逆だった。**")
        print("  **気温は土壌温度の「鈍った版」ではなく、白色雑音を余分に持つ**——")
        print("  **その雑音が残差に混ざり、ACF1 を薄める。**")
        for kind, want in (("soildriven", "**★になるべき**"),
                           ("airdriven", "**★を失うべき**")):
            d = synth(kind)
            v, why = verdict(d)
            if v is None:
                print(f"\n    `{kind}`：**判定できない**（{why}）"); continue
            st, m = v
            print(f"\n    `{kind}`（期待：{want}）")
            print(f"      R²={m['r2']:.2f}／ACF1={m['acf1']:+.2f}／e-fold={m['efold']:.0f}日"
                  f" → {'**★**' if st else '★でない'}")
        print("\n  → **soildriven→★・airdriven→★でない** なら、道具は向きを捉えている。")
        return

    root = Path(a.cosore_dir)
    desc = pd.read_csv(root / "description.csv")
    site_of = {}
    for _, r in desc.iterrows():
        try:
            site_of[str(r["CSR_DATASET"])] = (round(float(r["CSR_LATITUDE"]), 3),
                                              round(float(r["CSR_LONGITUDE"]), 3))
        except (TypeError, ValueError, KeyError):
            pass

    rows, skipped = [], []
    for f in sorted((root / "datasets").glob("data_*.csv")):
        ds = f.stem[5:]
        try:
            d, st, sm = load_cosore(f, None)
        except Exception as e:
            skipped.append((ds, f"{type(e).__name__}")); continue
        daily = d.groupby(d.index.normalize()).mean() if len(d) else d
        kind, depth = classify(st)
        v, why = verdict(daily)
        if v is None:
            skipped.append((ds, f"{why}／駆動列 {st or '無し'}")); continue
        star_, m = v
        rows.append({"ds": ds, "site": site_of.get(ds, ds), "col": st or "無し",
                     "kind": kind, "depth": depth, "star": bool(star_),
                     "acf1": m["acf1"], "efold": m["efold"], "r2": m["r2"],
                     "days": len(daily), "yrs": daily.index.year.nunique(),
                     "shallow": ("浅い(≤5cm)" if np.isfinite(depth) and depth <= SHALLOW_CM
                                 else ("深い(>5cm)" if np.isfinite(depth) else "—"))})

    sites = len({r["site"] for r in rows})
    n_star = sum(r["star"] for r in rows)
    print(f"\n  判定できた：**{len(rows)} 本／{sites} 地点**"
          f"（★ {n_star} 本＝**{n_star/len(rows):.0%}**）")
    print(f"  判定しなかった：{len(skipped)} 本")
    for ds, w in skipped:
        print(f"    {ds[:34]:<36}{w}")

    report_axis(rows, "軸A：土壌温度で駆動 vs 気温で駆動", "kind", "土壌", "気温")
    report_axis(rows, "軸B：浅い（≤5cm） vs 深い（>5cm）", "shallow", "浅い(≤5cm)", "深い(>5cm)")

    print("\n  === 参考（**主判定には使わない**）===")
    print("  **JASSAL は水分 3 層・温度 3 層を持つ唯一のデータセット**（旗101）。")
    print("  **同じチャンバーで駆動の深さだけ変えると★判定が変わるか**——**n=1 なので参考**。")
    j = root / "datasets" / "data_d20200108_JASSAL.csv"
    if j.exists():
        raw = pd.read_csv(j, low_memory=False)
        ts = pd.to_datetime(raw.get("CSR_TIMESTAMP_BEGIN", raw.get("CSR_TIMESTAMP_END")),
                            errors="coerce")
        rs = pd.to_numeric(raw["CSR_FLUX_CO2"], errors="coerce")
        for tc in [c for c in raw.columns if re.fullmatch(r"CSR_T\d+\.?\d*", c)]:
            for wc in [c for c in raw.columns if re.fullmatch(r"CSR_SM\d+\.?\d*", c)]:
                dd = pd.DataFrame({"Rs": rs.to_numpy(),
                                   "Tsoil": pd.to_numeric(raw[tc], errors="coerce").to_numpy(),
                                   "SM": pd.to_numeric(raw[wc], errors="coerce").to_numpy()},
                                  index=ts).dropna(subset=["Rs"])
                dd = dd[dd.index.notna()]
                dd = dd.groupby(dd.index.normalize()).mean()
                v, why = verdict(dd)
                if v is None:
                    print(f"    {tc:<10}×{wc:<10}判定できない（{why}）"); continue
                s2, m2 = v
                print(f"    {tc:<10}×{wc:<10}ACF1 {m2['acf1']:+.2f}／"
                      f"e-fold {m2['efold']:.0f}日／R² {m2['r2']:.2f} → "
                      f"{'**★**' if s2 else '★でない'}")
    else:
        print("    JASSAL が見つからない。")

    print("\n  留保（事前登録どおり）：")
    print("   ・**群分けは観測的である**——**気温しか無い地点は生態系も測定体制も違いうる。**")
    print("     **★率の差が「センサの正体」によるとは、この設計では言えない。**")
    print("     **言えるのは「★率が群で違うかどうか」だけ。交絡は解けない。**")
    print("   ・**擬似反復**：**SZUTU 5 本・KAYE 8 本・SIHI 3 本はそれぞれ 1 地点**。")
    print("     **本数と地点数を上に併記した。**")
    print("   ・**深さの表記は列名依存**——**実際の設置深度は確かめられない**（旗97/101 と同じ）。")
    print("   ・**★率が変わらなくても、深さが揃っていない事実は残る。**")
    print("     **それは結果に関わらず、限界として本文に書く。**")


if __name__ == "__main__":
    main()
