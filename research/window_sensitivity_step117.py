"""旗117：**手B の窓長感度**——旗116 の ▲ は Δθ の窓（10 日）に依るか。

旗116 は Δθ = θ(t) − θ(t−10日) の符号で層別し、**季節差は軌道を揃えても残る（▲）**と結論した。
**10 日は `OPEN_QUESTIONS_OPTIONS.md` の記述に合わせた恣意的な選択**で、事前登録 step116 の限界に
「窓長感度は主検定の後に別途」と書いた。**本旗はそれを果たす。**

**主検定の判定規則は一切変えない**（`trajectory_stratified_step116` をそのまま import・窓だけ振る）。
**窓を変えて結論（クラスタごとの分類）が動くかだけを見る。** **良い窓を選んで拾わない**——**全部並べる。**

    python research/window_sensitivity_step117.py                 # 既定クラスタ・窓 5/10/15/20
    python research/window_sensitivity_step117.py --only US-SRM --windows 5,10,15,20,30
"""
from __future__ import annotations

import argparse
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from trajectory_stratified_step116 import SITES, run_site, classify, daily_energy
from change_rate_probe_step115 import CLUSTER
from runlog import tee_stdout

WINDOWS = (5, 10, 15, 20)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--windows", default=None, help="カンマ区切り（例 5,10,15,20）")
    ap.add_argument("--qc-max", type=int, default=None)
    ap.add_argument("--verbose", action="store_true", help="各窓の全出力も表示")
    a = ap.parse_args()
    tee_stdout("step117")

    windows = tuple(int(x) for x in a.windows.split(",")) if a.windows else WINDOWS
    sites = SITES
    if a.only:
        want = {x.strip() for x in a.only.split(",") if x.strip()}
        sites = tuple(s for s in SITES if s in want) or tuple(sorted(want))

    print("=== 旗117：手B の窓長感度（旗116 の ▲ は窓に依るか）===")
    print(f"  窓＝{windows} 日。**判定規則は旗116 と同一・窓だけ振る。**")
    print(f"  対象＝{', '.join(sites)}。**結論が窓で動くかだけを見る（良い窓を拾わない）。**")

    # 日次データは 1 回だけ読む（読み込みが重い）。窓ごとに run_site を回す。
    table: dict[str, dict[int, str]] = {}
    for s in sites:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
        except Exception as e:                # noqa: BLE001
            print(f"\n  ━━ {s} ━━  読み込み失敗：{type(e).__name__}: {e}")
            table[s] = {w: "読み込み失敗" for w in windows}
            continue
        table[s] = {}
        for w in windows:
            buf = io.StringIO()
            with redirect_stdout(buf):
                diffs = run_site(s, d, window=w)
            verdict = "門①-a を通らない/測れない" if diffs is None else classify(diffs)
            table[s][w] = verdict
            if a.verbose:
                print(buf.getvalue())

    print("\n  === 窓長感度の表（クラスタごと）===")
    by_cl: dict[str, list] = {}
    for s in sites:
        by_cl.setdefault(CLUSTER.get(s, s), []).append(s)
    for cl, ss in by_cl.items():
        print(f"\n  【{cl}】")
        for s in ss:
            print(f"    {s}:")
            for w in windows:
                print(f"       窓 {w:>2} 日 → {table[s][w]}")

    print("\n  === 安定性の判定 ===")
    print("    ・各サイトで、窓を変えても分類（○/▲説明されない/▲部分的/門①落ち）の"
          "**主要な向きが変わらなければ、旗116 の結論は窓にロバスト**。")
    print("    ・変われば、10 日という選択が結論を作っていたことになる（その旨を旗116 の限界に追記）。")


if __name__ == "__main__":
    main()
