"""旗68：未登録サイトを登録するための下調べ——実ファイルを見て SiteSpec の中身を決める。

旗67 で同一地点の対が 2→7 に増えたが、**インドネシア泥炭3組（ID-PaB / ID-PaD / ID-Pag）は
タワー側が `KeyError` で読めなかった**。原因は `sites.py` の自動発見が
`**/*ALLVARS_HH_*.csv` しか探さないため——これらのサイトは**別の命名**らしい。

この3サイトは旗42/44 で**水分依存Q10 が逆符号**だった場所であり、同一地点のタワーで
確かめられれば価値が高い。だが**登録には実ファイルの形を知る必要がある**（推測で書かない）。

本ツールは各コードについて：
  1. データ配下で**そのコードを含むファイル**を種類別に数える
  2. HH（30分値）らしき csv を選び、**先頭行（列名）を出す**
  3. 既存の2つの変数マップ（`japanflux` / `base`）の**どちらがどれだけ一致するか**を採点
  4. その結果から **SiteSpec の登録案**を提示する（**適用は人が判断**）

    python research/register_site_probe_step68.py --codes ID-PaB ID-PaD ID-Pag \
        --data-dir /mnt/hdd/JAPANFLUX
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def probe(code, data_dir, max_show=12):
    root = Path(data_dir)
    hits = [p for p in root.rglob(f"*{code}*") if p.is_file()
            and "__MACOSX" not in p.parts and not p.name.startswith(("._", "~$"))]
    print(f"  ━━ {code} ━━")
    if not hits:
        print("    このコードを含むファイルが無い\n"); return
    ext = Counter(p.suffix.lower() for p in hits)
    print(f"    ファイル {len(hits)} 件  拡張子 {dict(ext)}")
    # 代表的なパス（相対）を数件
    for p in hits[:4]:
        print(f"      {p.relative_to(root)}")
    if len(hits) > 4:
        print(f"      …他 {len(hits)-4}")

    # HH らしき csv を選ぶ（サイズが大きい csv を優先）
    csvs = sorted([p for p in hits if p.suffix.lower() == ".csv"],
                  key=lambda p: p.stat().st_size, reverse=True)
    if not csvs:
        print("    csv が無い＝30分値の本体が別形式かもしれない\n"); return
    big = csvs[0]
    print(f"    最大の csv：{big.relative_to(root)}（{big.stat().st_size/1e6:.1f} MB）")
    try:
        import pandas as pd
        head = pd.read_csv(big, nrows=3)
    except Exception as e:
        print(f"    読み込み失敗 {type(e).__name__}: {e}\n"); return
    cols = list(head.columns)
    print(f"    列数 {len(cols)}／先頭 {cols[:max_show]}")

    # どちらの変数マップに合うか採点
    from japanflux_pn.sites import DEFAULT_VAR_MAP, DEFAULT_VAR_MAP_BASE
    best = None
    for name, vm in (("japanflux", DEFAULT_VAR_MAP), ("base", DEFAULT_VAR_MAP_BASE)):
        need = list(vm.values())
        hit = [v for v in need if v in cols]
        miss = [k for k, v in vm.items() if v not in cols]
        print(f"    変数マップ '{name}'：{len(hit)}/{len(need)} 一致"
              f"{'／欠け ' + ','.join(miss) if miss else ''}")
        if best is None or len(hit) > best[1]:
            best = (name, len(hit), vm, miss)

    # **欠けている変数の候補列を探す**（旗68 第1版はここを出しておらず、
    # 「列が無い」で終わって var_overrides を書けなかった）
    TOKENS = {"Ts": ("TS", "SOIL_T", "TSOIL"), "th": ("SWC", "SM_", "SOIL_W", "VWC"),
              "Rg": ("SW_IN", "RG"), "Ta": ("TA_", "TAIR"), "VPD": ("VPD",),
              "P": ("P_", "PREC"), "gH": ("H_",), "gLE": ("LE_",),
              "GER": ("RECO",), "NEE": ("NEE",), "GEP": ("GPP",)}
    if best and best[3]:
        print(f"    → 採用候補 '{best[0]}'。**欠けている変数の候補列**：")
        for k in best[3]:
            toks = TOKENS.get(k, (k.upper(),))
            cand = [c for c in cols if any(t in c.upper() for t in toks)]
            print(f"       {k:<4}: {cand[:10] if cand else '**候補なし**'}"
                  f"{' …他 %d' % (len(cand)-10) if len(cand) > 10 else ''}")
    # 登録案
    rel = big.relative_to(root)
    site_root = rel.parts[0] if len(rel.parts) > 1 else "."
    print(f"    **登録案**：SiteSpec(code={code!r}, data_dir={str(root / site_root)!r},")
    print(f"                      hh_glob={'**/' + big.name.replace(code, '*')!r}, fmt=<上の採点で決める>)")
    print()


def main():
    p = argparse.ArgumentParser(description="未登録サイトの下調べ")
    p.add_argument("--codes", nargs="+", required=True)
    p.add_argument("--data-dir", default="/mnt/hdd/JAPANFLUX")
    a = p.parse_args()
    print("=== 旗68：未登録サイトを登録するための下調べ ===")
    print("  推測で SiteSpec を書かず、**実ファイルの列名を見てから**決めるための下調べ。\n")
    for c in a.codes:
        probe(c, a.data_dir)
    print("  === 次にやること ===")
    print("  変数マップの一致が高い方（japanflux / base）を fmt に採り、hh_glob を実ファイル名に合わせて")
    print("  `japanflux_pn/sites.py` の SITES に登録する。**欠けている変数があれば var_overrides で補う**。")
    print("  登録後、旗66 を再実行すればインドネシア泥炭3組のタワー側が読める。")
    print("  留保：列名が合っても**単位・符号規約・u* フィルタの有無**が違う可能性がある。")
    print("        登録したら、まず既知サイトと**同じ季節の値域**を見て妥当性を確かめること。")


if __name__ == "__main__":
    main()
