"""旗106：**手C の事前登録が前提にする三つの事実**を、分布で確かめる——**検定はしない**。

**旗105 で分かったこと**（`FLAGS_LOG.md:4361`）：手C の前提（春と秋で雨の入り方が違う）は
**三度目にして初めて成り立った**が、**θ を揃える形（重なり帯）は 6 通り中 5 通りで下限未満**だった。
＝**θ の帯は捨て、θ は共変量に入れる形で登録する**——**これは旗105 の時点で決まっている。**

**その形にすると、旗105 が測っていない事実が三つ、前提になる。**
**事前登録を書く前に、それを分布で確かめる**（**旗97 で足した作法**。
**旗96/97 はどちらも一行で確かめられる分布を見ずに設計に入って外した**）。

  ① **門②の前提**：セル内で `dryspell` と θ が**同じ量になっていないか**（Spearman）。
     **「雨の直後は θ が高い」のは当たり前**なので、**言い換えなら層別する意味が無い。**
     **閾値は事前登録で固定する**（この道具は数を出すだけで、判定しない）。
  ② **層別軸の選択の前提**：`直前 7 日積算 P` は**春で中央 0.0 mm**（旗105）＝**同順位が多い**。
     **タイの割合**を数える。**多ければ層別軸としては `dryspell` の方が素直**である。
  ③ **設計の核の前提**：**秋の `遠い` 腕は、春のセルと θ 水準が揃っているか**。
     旗105 の出力では θ 中央が **秋遠い 10.583/4.220/4.856 対 春セル 9.781/4.521/5.030**（Wkg/Whs/SRM）
     と**近い**——**近いなら「同じ θ 水準で雨履歴だけが違う対」が季節をまたいで作れる。**
     **ここでは中央値と 10/90 を併記するだけで、帯は作らない**（**旗97/105 の自己矛盾を繰り返さない**）。

  さらに、**腕の間で Rg・Ta が偏っていないか**も併記する——**偏っていれば共変量に入れる必要がある。**

**この道具は θ→γLE・θ→γH を一行も計算しない。★も▲も出さない。**

    .venv/bin/python research/premise_probe_step106.py --check
    .venv/bin/python research/premise_probe_step106.py --real --sites US-Wkg US-Whs US-SRM
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaporation_regime_step36 import daily_energy
from stratified_bowen_step89 import cell_of, MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import SPRING, AUTUMN
from vpd_match_step96 import spearman
# **旗103 の道具をそのまま呼ぶ**——**しきい値も数え方も作り直さない**（旗105 の出力と突き合わせるため）
from rain_history_probe_step103 import (rain_history, daily_precip, synth,
                                        PRIMARY_THR, RECENT_MAX, REMOTE_MIN)


def _d(x, tag, unit="", fmt="{:.3f}"):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return f"{tag} なし"
    f = fmt.format
    return (f"{tag} 中央 {f(np.median(x))}{unit} "
            f"[{f(np.percentile(x, 10))}–{f(np.percentile(x, 90))}]")


def arms(u):
    """`直後`／`遠い` に割る。**旗103 と同一の境（≤3 日／≥7 日）**。"""
    return u[u["dry"] <= RECENT_MAX], u[u["dry"] >= REMOTE_MIN]


def ok_(g):
    return len(g) >= MIN_DAYS and (g.index.year.nunique() if len(g) else 0) >= MIN_YEARS


def season_premise(name, u):
    """1 季節。**①②と腕ごとの θ・Rg・Ta** を出す。**検定はしない。**"""
    if u.empty:
        print(f"      {name}：使える日が無い")
        return None
    r_dry = spearman(u["dry"], u["th"])
    r_c7 = spearman(u["cum7"], u["th"])
    tie = float((u["cum7"] <= 0).mean())
    print(f"      {name}：使える日 {len(u)}／{u.index.year.nunique()} 年")
    print(f"        ① セル内 Spearman r(dryspell, θ) = {r_dry:+.2f}"
          f"／r(直前7日積算P, θ) = {r_c7:+.2f}")
    print(f"        ② 直前 7 日積算 P が 0 の日＝{tie:.0%}"
          f"{'  ← **同順位が多い**' if tie >= 0.30 else ''}")
    rec, rem = arms(u)
    out = {"r_dry_th": r_dry, "r_cum7_th": r_c7, "tie0": tie,
           "rec": rec, "rem": rem, "all": u}
    for tag, g in ((f"直後(≤{RECENT_MAX})", rec), (f"遠い(≥{REMOTE_MIN})", rem)):
        if not len(g):
            print(f"        {tag}：0 日")
            continue
        print(f"        {tag}：{len(g)} 日／{g.index.year.nunique()} 年"
              f"  {'**下限を満たす**' if ok_(g) else '**下限未満**'}")
        print(f"          {_d(g['th'], 'θ')}／{_d(g['Rg'], 'Rg', fmt='{:.0f}')}"
              f"／{_d(g['Ta'], 'Ta', fmt='{:.1f}')}" if "Ta" in g else
              f"          {_d(g['th'], 'θ')}／{_d(g['Rg'], 'Rg', fmt='{:.0f}')}")
    return out


def run_site(tag, d, P):
    print(f"\n  ━━ {tag} ━━")
    if P is None or P.dropna().empty:
        print("    **降水 P が無い**")
        return None
    lab, tmed, rmed = cell_of(d)
    hh = d[lab == "θ高×Rg高"].copy()
    print(f"    しきい値（旗89 と同一・作り直さない）：θ={tmed:.3f}／Rg={rmed:.1f}"
          f"／θ高×Rg高 {len(hh)} 日")
    j = hh.join(rain_history(P, PRIMARY_THR), how="left")
    j["usable"] = j["usable"].fillna(False).astype(bool)
    out = {}
    for nm, mon in (("春", SPRING), ("秋", AUTUMN)):
        u = j[np.array([m in mon for m in j.index.month]) & j["usable"].to_numpy()]
        out[nm] = season_premise(nm, u)
    sp, au = out.get("春"), out.get("秋")
    print("    ── ③ **秋の `遠い` 腕は、春のセルと θ 水準が揃っているか** ──")
    if sp is None or au is None or not len(au["rem"]):
        print("      片方が空＝比べられない")
    else:
        a, b = au["rem"]["th"].to_numpy(), sp["all"]["th"].to_numpy()
        print(f"      {_d(a, '秋 遠い θ')}")
        print(f"      {_d(b, '春 セル θ')}")
        print(f"      中央の差 = {np.median(a) - np.median(b):+.3f} %vol"
              f"／秋遠いの中央は春の {np.mean(b <= np.median(a)):.0%} 分位")
        print(f"      {_d(au['rec']['th'], '（参考）秋 直後 θ')}")
    return out


def check():
    """**自己試験**——合成 `separable` で、①②③の枝すべてに到達するか。"""
    print("\n  【自己試験】合成 `separable` で**印字する枝に到達するか**を見る")
    print("  （**旗104 の欠陥 #36＝陽性の対照が判定枝に一度も到達していなかった**の再発防止）")
    d, P = synth("separable")
    out = run_site("合成-separable", d, P)
    hit = []
    for nm in ("春", "秋"):
        s = out.get(nm) if out else None
        if s is None:
            hit.append(f"{nm}:なし"); continue
        hit.append(f"{nm}:r_dry={s['r_dry_th']:+.2f} 直後{len(s['rec'])}/遠い{len(s['rem'])}")
    print(f"  → 到達した枝：{'／'.join(hit)}")
    good = out is not None and all(
        out.get(nm) is not None and len(out[nm]["rec"]) and len(out[nm]["rem"])
        and np.isfinite(out[nm]["r_dry_th"]) for nm in ("春", "秋"))
    print(f"  → 自己試験：{'**通過**' if good else '**未到達の枝がある**'}")
    return good


def main():
    ap = argparse.ArgumentParser(description="旗106：手C の事前登録が前提にする事実の下調べ")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--sites", nargs="+", default=["US-Wkg", "US-Whs", "US-SRM"])
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()

    print("=== 旗106：**手C の事前登録が前提にする三つの事実**（**検定はしない**）===")
    print("  ① dryspell と θ が同じ量でないか／② 積算 P の同順位／③ 秋の `遠い` と春の θ 水準。")
    print(f"  下限は旗58 以来と同じ（{MIN_DAYS} 日・{MIN_YEARS} 年）。"
          f"境は旗103 と同一（直後 ≤{RECENT_MAX} 日／遠い ≥{REMOTE_MIN} 日・イベント {PRIMARY_THR:.0f} mm）。")

    if a.check or not a.real:
        check()
        if not a.real:
            return

    for s in a.sites:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
            P = daily_precip(s, a.qc_max)
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読み込み失敗 {type(e).__name__}: {str(e)[:110]}")
            continue
        run_site(s, d, P)

    print("\n  === この下調べが事前登録に決めさせること（**結論ではない**）===")
    print("  ① r(dryspell, θ) が **1 に近ければ、雨履歴は θ の言い換え**＝手C は成立しない。")
    print("     **閾値は事前登録で固定する**（この道具は判定しない）。")
    print("  ② 同順位が多い軸は**主軸にしない**（感度解析に回す）。")
    print("  ③ **秋の `遠い` と春のセルの θ 水準が近ければ**、")
    print("     **「同じ θ 水準・違う雨履歴」の対が季節をまたいで作れる**＝設計の核が立つ。")
    print("\n  留保：**`P` はギャップフィル済みの可能性がある**（旗46 と同じ問題・区別できない）。")
    print("   **独立クラスタは 2 つ**（Walnut Gulch・Santa Rita）＝3 サイト≠3 反復。")
    print("   **この道具は何も検定していない。★も▲も出さない。**")


if __name__ == "__main__":
    main()
