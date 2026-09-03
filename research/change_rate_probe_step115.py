"""旗115：**手B の下調べ**——θ の「変化の向き」で層別できるか（feasibility のみ・検定しない）。

## なぜ手B か

**季節依存の説明を、θ・Rg・Ts・GCC・VPD・深層θ・雨からの日数の 7 変数で試して外した**
（旗89–113）。**すべて「ある量の“水準”を揃える」形だった。**
**手B は「水準」ではなく「傾き（変化の向き）」**——**春は昇温・乾燥へ向かう途中、
秋は降温・湿潤の直後で、同じ θ でも軌道が違う**（`OPEN_QUESTIONS_OPTIONS.md` 手B）。

**旗107 で Santa Rita は雨からの日数を揃えても季節差が残った（`or` 型）。**
**手B が問うのは：その残差は「軌道（Δθ の向き）」で説明されるか。**

## この道具がすること（**下調べのみ**）

**各サイトで、10 日間の θ の変化量 Δθ を作り、その符号で層別したとき、
（春／秋）×（Δθ>0／Δθ<0）の 4 群が下限（60 日・3 暦年）を満たすかだけを数える。**
**相関も Δ（結果側）も一度も計算しない**（旗96/97 の轍を踏まないため・事前登録の前に答えを見ない）。

**Δθ は層別変数であって結果ではない**——`rain_history` の `dryspell`（旗103）と同じ位置づけ。

    python research/change_rate_probe_step115.py            # 対象クラスタ（既定）
    python research/change_rate_probe_step115.py --only US-SRM
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
from stratified_bowen_step89 import MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import SPRING, AUTUMN
from runlog import tee_stdout

WINDOW = 10  # 日（OPEN_QUESTIONS 手B：「10 日間の変化量」）

# 反転が起きた乾燥クラスタ（旗106/107/111）。**手B は反転する所でしか意味がない。**
SITES = ("US-Wkg", "US-Whs", "US-SRM", "CN-Du2")
CLUSTER = {"US-Wkg": "Walnut Gulch", "US-Whs": "Walnut Gulch",
           "US-SRM": "Santa Rita", "CN-Du2": "Duolun"}


def delta_theta(d: pd.DataFrame, window: int = WINDOW) -> pd.Series:
    """日次 θ を暦日グリッドに並べ直し、window 日前との差 Δθ を返す（欠測は NaN）。

    ``daily_energy`` は dropna 済みで日付が飛ぶので、そのまま shift すると
    「10 行前」＝「10 暦日前」にならない。暦日で reindex してから差をとる。
    """
    th = d["th"].copy()
    full = pd.date_range(th.index.min(), th.index.max(), freq="D")
    thf = th.reindex(full)
    dth = thf - thf.shift(window)     # window 日前との差（間に欠測があれば NaN）
    return dth.reindex(th.index)      # 解析日（dropna 済み）だけに戻す


def counts(dth_season: pd.Series) -> dict:
    """ある季節の Δθ について、正/負それぞれの日数と暦年数を返す。"""
    out = {}
    for name, mask in (("上昇 Δθ>0", dth_season > 0), ("下降 Δθ<0", dth_season < 0)):
        s = dth_season[mask]
        yrs = sorted(set(s.index.year))
        out[name] = (len(s), len(yrs))
    return out


def run_site(tag: str, d: pd.DataFrame):
    print(f"\n  ━━ {tag}（{CLUSTER.get(tag, '?')}）━━")
    dth = delta_theta(d)
    j = pd.DataFrame({"dth": dth, "month": d.index.month}, index=d.index).dropna()
    ok = True
    for season, months in (("春", SPRING), ("秋", AUTUMN)):
        sub = j[[m in months for m in j["month"]]]["dth"]
        c = counts(sub)
        print(f"    {season}：Δθ 有効 {len(sub)} 日")
        for name, (n, ny) in c.items():
            flag = "" if (n >= MIN_DAYS and ny >= MIN_YEARS) else "  **下限未満**"
            if not (n >= MIN_DAYS and ny >= MIN_YEARS):
                ok = False
            print(f"      {season}×{name}：{n:4d} 日／{ny:2d} 年{flag}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="カンマ区切りでサイトを限定")
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()
    tee_stdout("step115")

    sites = SITES
    if a.only:
        want = {x.strip() for x in a.only.split(",") if x.strip()}
        sites = tuple(s for s in SITES if s in want) or tuple(sorted(want))

    print("=== 旗115：手B の下調べ（θ の変化の向きで層別できるか）===")
    print(f"  Δθ = θ(t) − θ(t−{WINDOW}日)。符号で層別。**検定はしない。**")
    print(f"  下限 {MIN_DAYS} 日・{MIN_YEARS} 暦年。季節 春{tuple(SPRING)}・秋{tuple(AUTUMN)}。")
    print(f"  対象＝反転が起きた乾燥クラスタ（旗106/107/111）：{', '.join(sites)}")

    result = {}
    for s in sites:
        try:
            # ★全12か月を明示する（months=None はサイト既定＝夏のみで春秋が空になる・
            #   ゲート道具 step109/113 と同じ呼び方にそろえる）
            d, _nyr = daily_energy(s, list(range(1, 13)), a.qc_max)
        except Exception as e:                # noqa: BLE001
            print(f"\n  ━━ {s} ━━  読み込み失敗：{type(e).__name__}: {e}")
            result[s] = None
            continue
        if "th" not in d or len(d) < MIN_DAYS:
            print(f"\n  ━━ {s} ━━  θ 無し/日数不足（{len(d)}）")
            result[s] = None
            continue
        result[s] = run_site(s, d)

    print("\n  === まとめ：4 群（春×秋 × Δθ 上昇×下降）が下限を満たすクラスタ ===")
    by_cluster: dict[str, list] = {}
    for s, ok in result.items():
        by_cluster.setdefault(CLUSTER.get(s, s), []).append((s, ok))
    viable = 0
    for cl, items in by_cluster.items():
        any_ok = any(ok for _s, ok in items if ok is not None)
        tag = "**下限を満たすサイトあり**" if any_ok else "**満たすサイトなし**"
        if any_ok:
            viable += 1
        detail = "・".join(f"{s}={'○' if ok else '×' if ok is not None else '—'}"
                           for s, ok in items)
        print(f"    {cl}：{tag}（{detail}）")
    print(f"\n  → **下限を満たす独立クラスタ数：{viable}**")
    print("    ・2 以上 → 手B を事前登録できる（旗116）。**ただし何が変わるかを先に書く。**")
    print("    ・1 以下 → 手B も 2 クラスタの壁に当たる。**軌道の軸は手元では検定できないと記す。**")


if __name__ == "__main__":
    main()
