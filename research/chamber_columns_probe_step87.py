"""旗87：**判定できなかった対のチャンバーに、何の列が在るのか**を確かめる（下調べ・検定はしない）。

旗86 の結論は **★で選ばない群 n=8** に載っている。判定できなかったのは 12 組で、
**最大の理由は「チャンバー側に土壌温度が無い」7 組**（カリフォルニアの 6 組＋MATHES）。
**救えれば n が倍近くになる**——**旗86 の 1/8 という数字が偶然かどうかを決められる。**

## **まず自分の誤りを訂正する**

私は前回「`CSR_TAIR` など別の温度列を使えば救えるかもしれない」と書いたが、
**`load_cosore._pick_soil_temp` は既に `CSR_TAIR` / `CSR_TAIR_AMB` へ落ちる**。
＝**その手はもう入っている**。だから残る可能性は次の三つで、**どれかは中を見ないと分からない**：

  1. **本当に温度列が無い**（＝救えない。そう確定して記録する）
  2. **列名が正規表現に当たらない**——`_pick_soil_temp` は
     ``re.fullmatch(r"CSR_T(\\d+\\.?\\d*)")`` という**厳しい形**しか拾わない。
     `CSR_TSOIL`・`CSR_T_5`・`CSR_T5CM`・`CSR_T5_1` はどれも**落ちる**。
     ＝**旗64/82/86 と同じ「名前で当てにいく」型**であり、**実際に起きていてもおかしくない**。
  3. **別ファイル（ancillary 等）に在る**＝データセット csv の外を見る必要がある。

## やること（**三つだけ。判定はしない**）

  1. 対象データセットの **全列名をそのまま出す**（**推測せず、在るものを見る**）
  2. **現在の抽出規則が何を拾い、何を捨てたか**を並べる
     ＝**温度・水分らしいのに拾われていない列**を名指しする
  3. COSORE ディレクトリに**データセット csv 以外のファイル**が在るかを見る

**救える列が見つかっても、ここでは使わない。** 拾い方を直すかどうかは**人が決める**——
`CSR_TAIR`（気温）と `CSR_T5`（地温 5cm）は**別の量**であり、
**気温で代用した結果を土壌温度と同じ物差しに並べてよいかは、自明ではない**。

    python research/chamber_columns_probe_step87.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import _pick_soil_temp, _pick_sm

# 旗86 で判定できなかった 12 組のチャンバー。**なぜ落ちたか**で分けて見る。
NO_TSOIL = [                      # ①「チャンバー側に土壌温度が無い」＝本ツールの主対象
    "d20200214_SZUTU_TONZI", "d20191017_BALDOCCHI", "d20200214_SZUTU_VAIRA",
    "d20200214_SZUTU_BOULDINA", "d20200214_SZUTU_BOULDINC",
    "d20200214_SZUTU_TWITCHELL", "d20200224_MATHES",
]
OTHER_FAIL = [                    # ②温度は在ったが判定に届かなかった組（比較のため並べる）
    "d20190527_GOULDEN",          # R²=0.28（駆動弱）
    "d20190610_SIHI_H2wetland",   # 判定不能
    "d20200221_MATHES",           # 判定不能
    "d20190520_RUEHR",            # タワー側の期間が重ならない
    "d20190610_SIHI_H1",          # タワー側の期間が重ならない
]
# **判定できた対**も 3 件だけ並べる＝**「在るとはどう見えるか」の対照**（旗52 の作法）
CONTROL = ["d20190526_PENNINGTON", "d20200419_PEREZ-QUEZADA", "d20200220_GAVAZZI"]

# 温度・水分「らしい」列を、**現在の規則より広く**拾う（**採用はしない。名指しするだけ**）
TEMP_LIKE = re.compile(r"(TEMP|_T[_0-9]|TSOIL|TS_|TAIR|_TS$|_T$)", re.IGNORECASE)
SM_LIKE = re.compile(r"(SM|SWC|MOIST|VWC|WATER)", re.IGNORECASE)


def probe_one(path: Path):
    """1 データセットの列を、**現在の規則が拾ったか**で仕分けて返す。"""
    try:
        head = pd.read_csv(path, nrows=200, low_memory=False)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {str(e)[:80]}"}
    cols = list(head.columns)
    st, sm = _pick_soil_temp(cols), _pick_sm(cols)
    # **拾われなかったが温度・水分らしい列**＝ここが本ツールの眼目
    miss_t = [c for c in cols if c != st and TEMP_LIKE.search(str(c))]
    miss_s = [c for c in cols if c != sm and SM_LIKE.search(str(c))]
    # 中身が在るかも見る（**列名だけでは足りない**——旗80 で「列は在るが空」を踏んだ）
    def n_ok(c):
        if c is None or c not in head:
            return None
        return int(pd.to_numeric(head[c], errors="coerce").notna().sum())
    return {"cols": cols, "picked_t": st, "picked_sm": sm,
            "picked_t_n": n_ok(st), "picked_sm_n": n_ok(sm),
            "miss_t": [(c, n_ok(c)) for c in miss_t],
            "miss_s": [(c, n_ok(c)) for c in miss_s]}


def show(label, names, root, full_cols):
    print(f"\n  ── {label} ──")
    for ds in names:
        f = root / "datasets" / f"data_{ds}.csv"
        print(f"  ━━ {ds} ━━")
        if not f.exists():
            print("    **ファイルが無い**"); continue
        r = probe_one(f)
        if "error" in r:
            print(f"    **読めない**：{r['error']}"); continue
        t, s = r["picked_t"], r["picked_sm"]
        t_txt = f"{t}（有効 {r['picked_t_n']}）" if t else "**無し**"
        s_txt = f"{s}（有効 {r['picked_sm_n']}）" if s else "**無し**"
        print(f"    現在の規則が拾った：温度={t_txt}／水分={s_txt}")
        if r["miss_t"]:
            print(f"    **温度らしいのに拾われていない列**："
                  f"{', '.join(f'{c}(有効{n})' for c, n in r['miss_t'])}")
        elif not t:
            print(f"    → **温度らしい列は本当に一つも無い**（救えない）")
        if r["miss_s"]:
            print(f"    水分らしいのに拾われていない列："
                  f"{', '.join(f'{c}(有効{n})' for c, n in r['miss_s'])}")
        if full_cols:
            print(f"    全列（{len(r['cols'])}）：{r['cols']}")


def main():
    p = argparse.ArgumentParser(description="チャンバー側の列を確かめる（下調べ）")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--full-cols", action="store_true", help="全列名を出す")
    a = p.parse_args()
    root = Path(a.cosore_dir)

    print("=== 旗87：判定できなかった対のチャンバーに何の列が在るか（下調べ・検定はしない）===")
    print("  **訂正**：前回『CSR_TAIR を使えば救えるかも』と書いたが、")
    print("  **`_pick_soil_temp` は既に `CSR_TAIR`/`CSR_TAIR_AMB` へ落ちる**＝その手は入っている。")
    print("  残る可能性は **(1) 本当に無い (2) 名前が正規表現に当たらない (3) 別ファイルに在る**。")
    print("  **現在の規則**：温度 `CSR_T<数字>`（5cm に最も近い層）→ 無ければ `CSR_TAIR`。")
    print("             ＝`CSR_TSOIL`・`CSR_T_5`・`CSR_T5CM` は**どれも落ちる**。\n")

    show("① 土壌温度が無いとされた 7 組（**主対象**）", NO_TSOIL, root, a.full_cols)
    show("② 温度は在ったが判定に届かなかった 5 組（比較）", OTHER_FAIL, root, a.full_cols)
    show("③ 判定できた 3 組（**対照＝在るとはどう見えるか**）", CONTROL, root, a.full_cols)

    # ── データセット csv 以外に何が在るか ──
    print(f"\n  ── COSORE ディレクトリの中身（**csv 以外に手がかりが無いか**）──")
    tops = sorted(p_ for p_ in root.iterdir() if p_.name != "datasets")
    for t in tops[:20]:
        kind = "ディレクトリ" if t.is_dir() else f"{t.stat().st_size/1e6:.1f} MB"
        print(f"    {t.name:<40}{kind}")
    if len(tops) > 20:
        print(f"    …ほか {len(tops)-20} 件")

    # ── 第2段：ancillary.csv に土壌温度が在るか ──
    # **第1段で「7 組とも本当に温度列が無い」と出た**＝名前の問題ではなかった。
    # 残るのは可能性(3)＝**別ファイル**。`ancillary.csv` は 87MB あり、**一番大きい**。
    anc = root / "ancillary.csv"
    print(f"\n  ── 第2段：`ancillary.csv` に土壌温度が在るか ──")
    if not anc.exists():
        print("    **ancillary.csv が無い**＝可能性(3) も潰れた＝**救えない**")
    else:
        # 小さい説明ファイルを先に見る（**87MB を読む前に、何が入っているかを知る**）
        for meta in ("columns.csv", "CSR_COLUMN_UNITS.csv"):
            f = root / meta
            if not f.exists():
                continue
            try:
                m = pd.read_csv(f)
            except Exception as e:
                print(f"    {meta} を読めない（{type(e).__name__}）"); continue
            txt = m.astype(str).apply(lambda c: c.str.upper())
            hit = m[txt.apply(lambda c: c.str.contains("TEMP|TSOIL|SOIL T", na=False)).any(axis=1)]
            print(f"    {meta}：{len(m)} 行／**温度に触れる行 {len(hit)}**")
            for _, r in hit.head(8).iterrows():
                print(f"       {dict(list(r.items())[:4])}")
        try:
            head = pd.read_csv(anc, nrows=5, low_memory=False)
        except Exception as e:
            print(f"    ancillary.csv を読めない（{type(e).__name__}: {str(e)[:80]}）")
            head = None
        if head is not None:
            print(f"    ancillary.csv の列（{len(head.columns)}）：{list(head.columns)}")
            key = next((c for c in head.columns
                        if str(c).upper() in ("CSR_DATASET", "DATASET")), None)
            tcols = [c for c in head.columns if TEMP_LIKE.search(str(c))]
            print(f"    データセット列：{key or '**無し**'}"
                  f"／温度らしい列：{tcols if tcols else '**無し**'}")
            if key is None:
                print("    → **データセットを特定する列が無い**＝**対応付けられない**")
            elif not tcols:
                print("    → **温度らしい列が無い**＝**ここにも土壌温度は無い**")
            else:
                # **7 組の行が在るか、在るなら中身が在るか**を分割読みで数える
                want = set(NO_TSOIL)
                use = [key] + tcols
                n_row = {d: 0 for d in want}
                n_ok = {d: {c: 0 for c in tcols} for d in want}
                try:
                    for ch in pd.read_csv(anc, usecols=lambda c: c in set(use),
                                          chunksize=200_000, low_memory=False):
                        sub = ch[ch[key].astype(str).isin(want)]
                        if sub.empty:
                            continue
                        for d, g in sub.groupby(sub[key].astype(str)):
                            n_row[d] += len(g)
                            for c in tcols:
                                n_ok[d][c] += int(pd.to_numeric(g[c], errors="coerce")
                                                  .notna().sum())
                except Exception as e:
                    print(f"    分割読みに失敗（{type(e).__name__}: {str(e)[:80]}）")
                else:
                    print(f"    {'chamber':<32}{'行数':>9}  温度列の有効数")
                    for d in NO_TSOIL:
                        got = ", ".join(f"{c}={n_ok[d][c]:,}" for c in tcols
                                        if n_ok[d][c] > 0) or "**すべて 0**"
                        print(f"    {d:<32}{n_row[d]:>9,}  {got}")
                    print("    → **行が在り温度に中身が在る組だけ**が救える見込み。")
                    print("      **時刻列が無ければ時系列にできない**ので、列名も上で確かめること。")

    print("\n  === 次の判断 ===")
    print("  ・**温度らしい列が在るのに拾われていなかった** → **抽出規則の欠陥**。")
    print("    直せば **★で選ばない群の n が 8 → 最大 15** になる＝**旗86 の 1/8 が偶然か決まる**。")
    print("  ・**本当に無い** → **救えないと確定して記録する**。n=8 のまま幅を明示して述べる。")
    print("  留保：")
    print("   ・**`CSR_TAIR`（気温）と `CSR_T5`（地温5cm）は別の量**である。")
    print("     気温で代用した対を、地温で測った対と**同じ表に並べてよいかは自明でない**")
    print("     （旗33 の『深度不統一』と同じ問題が、**深度どころか媒体の違い**として出る）。")
    print("     **代用するなら、代用した対には印を付けて別に数える。**")
    print("   ・**列が在ることと中身が在ることは別**（旗80 で踏んだ）＝有効数も併記してある。")
    print("     ただし**先頭 200 行しか見ていない**＝**全期間の有効数ではない**。")


if __name__ == "__main__":
    main()
