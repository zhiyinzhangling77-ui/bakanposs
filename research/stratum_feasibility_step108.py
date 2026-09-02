"""旗108：**手C の層別が、他のサイトでもできるか**（下調べ・**検定はしない**）。

## なぜこれを見るのか

**旗106・旗107 で、二つのクラスタが一貫して違う振る舞いを見せた**——
Walnut Gulch は「雨で説明できる」側、Santa Rita は「季節差が残る」側。
**そして Santa Rita は、事前に用意した三つの世界のどれとも合わなかった**
（`遠い` で差が残り、`直後` では差が無い＝**「雨の直後 **または** 秋」という論理和の型**）。

**だが独立クラスタは 2 しかない。**
**同じ 2 クラスタを何度叩いても、この型が本物かは決まらない。**
**足りないのはクラスタ数である。**

## **この道具は実行可能性だけを出す。相関も Δ も計算しない。**

**それは検定の答えそのものであり、事前登録の前に見てはいけない**
（旗94/98/101/103 と同じ作法）。**出すのは日数・年数・θ の中央値だけである。**

数えるのは **4 群**：**春×`直後`／春×`遠い`／秋×`直後`／秋×`遠い`**。
**旗107 で分かったとおり、`直後` 層が `rain_only` と `both` を区別する唯一の軸**なので、
**4 群すべてが下限（60 日・3 暦年）を満たすサイト**を探す。

    python research/stratum_feasibility_step108.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stratified_bowen_step89 import cell_of, MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import SPRING, AUTUMN
from rain_history_probe_step103 import (rain_history, daily_precip,
                                        PRIMARY_THR, RECENT_MAX, REMOTE_MIN)
from evaporation_regime_step36 import daily_energy

# **旗79/80 で登録した全サイト**（`sites.py` から機械的に拾う。手で並べない）
DONE = {"US-Wkg", "US-Whs", "US-SRM"}      # 旗106/107 で既に測った＝3 クラスタ目の候補から除く
# **実サイトではない登録**。`sites.py` の `A3-base` は雛形で、
# `data_dir` が "← 実パスに変更" のまま置かれている（**読めないのではなく、存在しない**）。
# **「読めない」に混ぜると原因を取り違えるので、名指しで除く。**
TEMPLATES = {"A3-base"}


def all_sites():
    """**`get_site` が知っているサイトを全部**拾う（**手で並べない**）。

    **道具の欠陥 #39（旗109 の後に判明）**：第1版は**静的な `SITES` 辞書だけ**を読んでいた。
    **`sites.py` には登録簿が四つある**——
    `SITES`（手登録）・`discover_japanflux_sites`（`/mnt/hdd/JAPANFLUX`）・
    `discover_chinaflux_sites`（`/mnt/hdd/ChinaFlux`）・`discover_koflux_sites`。
    **`get_site` はこの四つを順に見る**（見つからないときの例外に
    `sorted(set(SITES) | set(disc) | set(cn) | set(kr))` と書いてある）。

    **＝第1版は「A-3 の元のクラスタであるモンゴルのサイト」も
    「ChinaFlux の乾燥草原」も数えていなかった。**
    **旗108 の結論「手元の乾燥地クラスタは 2 が上限」は、この不完全な数え上げの上に立っていた。**
    **本版で数え直すまで、その結論は保留する。**

    **`get_site` 自身の定義をそのまま使う**——**別の並べ方を自分で書かない。**
    """
    from japanflux_pn import sites as S
    names = set(S.SITES)
    for fn in ("discover_japanflux_sites", "discover_chinaflux_sites",
               "discover_koflux_sites"):
        f = getattr(S, fn, None)
        if f is None:
            continue
        try:
            got = f()
        except Exception as e:            # ルートが無い環境では空 dict が返る設計
            print(f"  （{fn} が失敗：{type(e).__name__}）")
            continue
        names |= set(got)
    return sorted(names)


def four_groups(d, P):
    """**4 群の日数・年数・θ 中央**を返す。**相関は計算しない。**"""
    lab, tmed, rmed = cell_of(d)
    hh = d[lab == "θ高×Rg高"]
    j = hh.join(rain_history(P, PRIMARY_THR), how="left")
    j["usable"] = j["usable"].fillna(False).astype(bool)
    u = j[j["usable"]]
    out = {}
    for sn, mon in (("春", SPRING), ("秋", AUTUMN)):
        s = u[[m in mon for m in u.index.month]]
        for ln, sel in (("直後", s[s["dry"] <= RECENT_MAX]),
                        ("遠い", s[s["dry"] >= REMOTE_MIN])):
            out[f"{sn}×{ln}"] = (len(sel),
                                 sel.index.year.nunique() if len(sel) else 0,
                                 float(np.nanmedian(sel["th"])) if len(sel) else np.nan)
    return out, tmed, rmed, len(hh)


def main():
    ap = argparse.ArgumentParser(description="旗108：層別が他のサイトでもできるか")
    ap.add_argument("--qc-max", type=int, default=None)
    ap.add_argument("--cosore-dir", default=None, help="使わない（引数の互換のため）")
    a = ap.parse_args()

    print("=== 旗108：手C の層別が、他のサイトでもできるか（下調べ・検定はしない）===")
    print("  **相関も Δ も計算しない**——**それは検定の答えであり、事前登録の前に見てはいけない。**")
    print(f"  数えるのは **4 群**（春×`直後`／春×`遠い`／秋×`直後`／秋×`遠い`）。")
    print(f"  下限は **{MIN_DAYS} 日・{MIN_YEARS} 暦年**／`直後` ≤{RECENT_MAX} 日・"
          f"`遠い` ≥{REMOTE_MIN} 日／イベント {PRIMARY_THR:.0f} mm（旗103–107 と同一）。")
    print("  **旗107 で `直後` 層が唯一の識別軸**と分かったので、**4 群すべて**を要求する。\n")

    sites = [s for s in all_sites() if s not in TEMPLATES]
    print(f"  登録サイト：{len(sites)} 件"
          f"（**手登録＋JapanFlux／ChinaFlux／KoFlux の自動発見をすべて合わせた**）"
          f"／雛形を除いた：{sorted(TEMPLATES)}")
    print("  **第1版は手登録の 31 件しか数えていなかった**（道具の欠陥 #39）。\n")
    full, partial, no_p, no_th, unread = [], [], [], [], []
    for s in sites:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
        except Exception as e:
            unread.append((s, f"{type(e).__name__}: {str(e)[:50]}")); continue
        try:
            P = daily_precip(s, a.qc_max)
        except Exception as e:
            unread.append((s, f"P: {type(e).__name__}")); continue
        if P is None or P.dropna().empty:
            no_p.append(s); continue
        # **道具の欠陥 #40（旗108 の再実行で判明）**：
        # **`KeyError: 'dry'` を「読めない」に混ぜていた**。
        # 実際の原因は**土壌水分 θ が全 NaN で `θ高×Rg高` セルが空**になり、
        # 空の枠に降雨履歴を join しても列が付かないことである
        # （警告に「変数 th の写像先 'SWC_F_MDS' が … に無い＝**この変数は全 NaN**」と出ている）。
        # **「読めない」と「θ が無い」は別の事実**なので分けて数える（欠陥 #34 と同じ形・2 度目）。
        if d.empty or "th" not in d.columns or not np.isfinite(d["th"]).any():
            no_th.append(s); continue
        try:
            g, tmed, rmed, ncell = four_groups(d, P)
        except KeyError as e:
            if str(e).strip("'\"") == "dry":
                no_th.append(s); continue
            unread.append((s, f"層別: KeyError: {str(e)[:40]}")); continue
        except Exception as e:
            unread.append((s, f"層別: {type(e).__name__}: {str(e)[:40]}")); continue
        ok = {k: (n >= MIN_DAYS and y >= MIN_YEARS) for k, (n, y, _) in g.items()}
        n_ok = sum(ok.values())
        mark = ("**4 群すべて**" if n_ok == 4 else
                f"{n_ok}/4 群" if n_ok else "**0 群**")
        tag = "  ← **既に測った**" if s in DONE else ""
        print(f"  ━━ {s} ━━ {mark}{tag}")
        print(f"    θ高×Rg高 {ncell} 日／しきい値 θ={tmed:.3f}・Rg={rmed:.1f}"
              f"／P {len(P)} 日・{P.index.year.nunique()} 年")
        for k in ("春×直後", "春×遠い", "秋×直後", "秋×遠い"):
            n, y, th = g[k]
            print(f"      {k:<8}{n:>5} 日／{y:>2} 年／θ 中央 "
                  f"{th:>7.2f}  {'**使える**' if ok[k] else '下限未満'}")
        if n_ok == 4:
            full.append(s)
        elif n_ok:
            partial.append((s, n_ok))
        print()

    print("  === まとめ ===")
    print(f"  **4 群すべて使える：{len(full)} 件** → {full}")
    print(f"  一部だけ：{len(partial)} 件 → {partial}")
    print(f"  降水 P が無い：{len(no_p)} 件 → {no_p}")
    print(f"  **土壌水分 θ が無い（＝セルを作れない）：{len(no_th)} 件** → {no_th}")
    print("    （**第1版はこれを『読めない』に混ぜていた**＝道具の欠陥 #40）")
    if unread:
        print(f"  読めない：{len(unread)} 件")
        for s, w in unread:
            print(f"    {s:<10}{w}")
    new = [s for s in full if s not in DONE]
    print(f"\n  **既に測った 3 サイトを除いた「新しく使えるサイト」：{len(new)} 件** → {new}")

    print("\n  === 次の判断（**事前登録の前に決める**）===")
    print("  ・**新しく使えるサイトが 1 件以上**（**かつ既存 2 クラスタと 50 km 以上離れている**）")
    print("    → **旗109 を事前登録**して、**「雨の直後 または 秋」という型を検定する**")
    print("      （**Santa Rita で見えた型が、独立なクラスタで再現するか**）")
    print("  ・**新しく使えるサイトが 0 件** → **手元では 2 クラスタが限界**と記し、")
    print("    **「雨または秋」は事後のパターン合わせのまま残す**と明記して新規観測へ渡す")
    print("  ・**一部の群だけ使えるサイトがある** → **どの群が落ちるかを書く。**")
    print("    **下限は緩めない**（旗93 で緩めなかったのと同じ）。")
    print("\n  留保：")
    print("   ・**相関も Δ も一度も計算していない**（事前登録の前に答えを見ない）。")
    print("   ・**4 群が揃うことと、そのサイトで Bowen 反転が起きることは別**である。")
    print("     **反転しないサイトでは、層別しても何も言えない**（旗106 の門①-a）。")
    print("     **それは旗109 の門で受ける**——**ここでは反転を見ない。**")
    print("   ・**クラスタは座標で決める**（旗82/107 と同じ 50 km 単連結）。")
    print("     **同じ流域の 2 サイトは 1 クラスタである**（Walnut Gulch の教訓）。")


if __name__ == "__main__":
    main()
