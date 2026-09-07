"""旗133：**しきい値感度**——旗106/107 の型（雨からの日数の層別）は `直後≤3日 / 遠い≥7日` に依るか。

旗107 は **`遠い` 層（雨から ≥7 日）の中の季節差 Δ_H** を主判定にし、
**Walnut Gulch は季節差なし（＝雨で説明・`rain_only`）、Santa Rita は季節差が残る（`or`）**と結論した。
**7 日（と 3 日）は旗103 で置いた候補値**で、事前登録 step107 の限界に「しきい値感度は後で」と書いた。**本旗はそれを果たす。**

**判定規則は一切変えない**——**旗107 の実際の道具 `season_within_stratum_step107.run_site` をそのまま使い、
しきい値だけ振る。** `dry`（雨からの日数）の計算＝イベント 5mm は不変で、**動くのは層の境界だけ**。
**良い値を選んで拾わない**——**全しきい値の結果を並べる。**

  返り値：True=季節差が残る（or 型）／False=季節差なし（rain_only 型）／None=判定不能（下限未満）。
  主軸は `遠い` 層なので、効くのは主に REMOTE_MIN（`遠い` 境界）。

    python research/threshold_sensitivity_step133.py
    python research/threshold_sensitivity_step133.py --only US-SRM --pairs 2/5,3/7,4/10
"""
from __future__ import annotations

import argparse
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import season_within_stratum_step107 as _s7  # noqa: E402  (RECENT_MAX/REMOTE_MIN を実行時に差し替える)

SITES = ("US-Wkg", "US-Whs", "US-SRM", "CN-Du2")   # 反転が起きた乾燥クラスタ（旗106/107/111）
CLUSTER = {"US-Wkg": "Walnut Gulch", "US-Whs": "Walnut Gulch",
           "US-SRM": "Santa Rita", "CN-Du2": "Duolun"}
PAIRS = ((2, 5), (3, 7), (4, 10))   # (RECENT_MAX, REMOTE_MIN)・(3,7) が旗106/107 の基準


def _label(v):
    if v is True:
        return "季節差が残る（or 型）"
    if v is False:
        return "季節差なし（rain_only 型）"
    return "判定不能（遠い層が下限未満）"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--pairs", default=None, help="例 2/5,3/7,4/10")
    ap.add_argument("--qc-max", type=int, default=None)
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    _s7.tee_stdout("step133")

    pairs = PAIRS
    if a.pairs:
        pairs = tuple(tuple(int(x) for x in p.split("/")) for p in a.pairs.split(","))
    sites = SITES
    if a.only:
        want = {x.strip() for x in a.only.split(",") if x.strip()}
        sites = tuple(s for s in SITES if s in want) or tuple(sorted(want))

    print("=== 旗133：しきい値感度（旗107 の型は 直後≤3/遠い≥7 に依るか）===")
    print(f"  しきい値対（直後≤R / 遠い≥M）＝{pairs}。**(3,7) が旗106/107 の基準。**")
    print(f"  旗107 の run_site をそのまま・しきい値だけ振る（主軸＝遠い層 Δ_H）。対象＝{', '.join(sites)}。")

    loaded = {}
    for s in sites:
        try:
            d, _ = _s7.daily_energy(s, list(range(1, 13)), a.qc_max)
            P = _s7.daily_precip(s, a.qc_max)
            loaded[s] = (d, P)
        except Exception as e:                       # noqa: BLE001
            print(f"\n  ━━ {s} ━━ 読み込み失敗：{type(e).__name__}: {e}")
            loaded[s] = None

    table: dict[str, dict] = {}
    for s in sites:
        table[s] = {}
        if loaded[s] is None:
            for pr in pairs:
                table[s][pr] = "読み込み失敗"
            continue
        d, P = loaded[s]
        for (R, M) in pairs:
            _s7.RECENT_MAX, _s7.REMOTE_MIN = R, M     # pick() が実行時に参照する名前を差し替え
            buf = io.StringIO()
            with redirect_stdout(buf):
                v = _s7.run_site(s, d, P)
            table[s][(R, M)] = _label(v)
            if a.verbose:
                print(buf.getvalue())

    print("\n  === しきい値感度の表（クラスタごと）===")
    by_cl: dict[str, list] = {}
    for s in sites:
        by_cl.setdefault(CLUSTER.get(s, s), []).append(s)
    for cl, ss in by_cl.items():
        print(f"\n  【{cl}】")
        for s in ss:
            print(f"    {s}:")
            for pr in pairs:
                star = " ←基準" if pr == (3, 7) else ""
                print(f"       直後≤{pr[0]}/遠い≥{pr[1]} → {table[s][pr]}{star}")

    print("\n  === 読み ===")
    print("    ・各クラスタで型がしきい値を跨いで安定なら、旗106/107 の結論はしきい値にロバスト。")
    print("    ・型が動くなら、3日・7日という選択が結論を作っていた（旗107 の限界に追記）。")


if __name__ == "__main__":
    main()
