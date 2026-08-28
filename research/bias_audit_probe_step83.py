"""旗83：**選定バイアスを測るための材料が手元にあるか**を確かめる（下調べ・検定はしない）。

前提の監査で四つの穴を挙げた。うち**三つは手元で測れる可能性がある**——だが
**推測で進めない**（旗68/69/76/79 と同じ作法）。**まず材料の有無を確かめる。**

## ① 私が作ったバイアス：**★の順で取得サイトを選んだ**

旗78 は**外挿基準で★が出るチャンバー**を上位に並べ、私はその上位を取得するよう勧めた。
＝旗66 拡張版の「両側★ 9 組」は **「チャンバーが★だった場所で、タワーも★か」**を測っており、
**「一般に一致するか」ではない**。

**これを消す方法がある**：**AmeriFlux 全サイトの座標**があれば、
**★で選ばずに COSORE 全チャンバーと総当たり**できる。
＝**取得しなかった対がいくつあり、そこで一致率が違うか**が分かる。
取得ページで **multi-site BADM（`AA-Flx_BIF`）** を含める選択肢があった——**あるかどうかを見る**。

## ② 研究者・機材の独立性

旗43 のクラスタ縮約は**地理だけ**で、**同じ研究者・同じ機材・同じ処理選択**を独立と数えている。
北米乾燥3サイトは**同一PI**らしい＝**地理クラスタ3でも研究群は2**かもしれない。
COSORE の `description.csv` に**寄与者の列**があれば、**研究者単位で縮約し直せる**。

## ③ ふるいの脱落

`R²≥0.3` は「駆動でよく説明できる」サイトを選ぶふるいで、**メモリの出方と相関しうる**。
**各段階で何件落ちたか**を数える材料は、COSORE のファイルがあれば揃っている。

本ツールは**この三つの材料の有無だけ**を報告する。**無ければ「測れない」と確定して記録する。**

    python research/bias_audit_probe_step83.py --cosore-dir /mnt/hdd/cosore-0.7.0 \
        --amf-dir /mnt/hdd/AmeriFlux_FLUXNET
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

# 寄与者・研究者を示しうる列名（COSORE の規約と、実データで見かけうる別名）
WHO_COLS = ("CSR_CONTRIBUTOR", "CSR_PI", "CSR_INVESTIGATOR", "CSR_AUTHOR",
            "CSR_DATASET_CONTRIBUTOR", "CSR_PRIMARY_PUB", "CSR_CITATION",
            "CSR_SITE_NAME", "CSR_NETWORK")


def main():
    p = argparse.ArgumentParser(description="選定バイアスを測る材料があるか")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--amf-dir", default="/mnt/hdd/AmeriFlux_FLUXNET")
    a = p.parse_args()

    print("=== 旗83：選定バイアスを測る材料が手元にあるか（下調べ・検定はしない）===")
    print("  **推測で進めない**。材料が無ければ「測れない」と確定して記録する。\n")

    # ── ① AmeriFlux 全サイトの座標（★で選ばない総当たりに要る）──
    print("  ── ① AmeriFlux **全サイト**の座標（★選択バイアスを消すのに要る）──")
    amf = Path(a.amf_dir)
    multi = [p_ for p_ in amf.rglob("*")
             if p_.is_file() and "AA-Flx" in p_.name and "BIF" in p_.name.upper()]
    per_site = [p_ for p_ in amf.rglob("*") if p_.is_file()
                and "BIF" in p_.name.upper() and "AA-Flx" not in p_.name]
    print(f"    multi-site BADM（AA-Flx_BIF）：{len(multi)} 件"
          f"{'  ' + multi[0].name if multi else ''}")
    print(f"    サイト個別の BIF：{len(per_site)} 件")
    src = multi[0] if multi else (per_site[0] if per_site else None)
    if src is None:
        print("    → **どちらも無い＝全サイト総当たりはできない**")
        print("       （取得ページの『Include multi-site BADM file』を入れて取り直すと可能）")
    else:
        try:
            # **BIF は utf-8 でないことがある**（旗83 で UnicodeDecodeError）。
            # 旗79 も同じ弱さを持ち、**例外を握って黙って飛ばしていた**。
            df = pd.read_csv(src, low_memory=False,
                             encoding="latin-1", encoding_errors="replace")
            cols = {c.upper(): c for c in df.columns}
            sid = next((cols[k] for k in ("SITE_ID", "SITEID") if k in cols), None)
            var = next((cols[k] for k in ("VARIABLE", "VARIABLE_NAME") if k in cols), None)
            if sid and var:
                n_site = df[sid].nunique()
                has_lat = df[var].astype(str).str.upper().eq("LOCATION_LAT").sum()
                print(f"    {src.name}：**{n_site} サイト**／LOCATION_LAT の行 {has_lat} 件")
                if n_site > 50:
                    print(f"    → **全サイト総当たりができる**（{n_site} サイト分の座標がある）")
                else:
                    print(f"    → **取得した分だけ**（{n_site} サイト）＝総当たりには足りない")
            else:
                print(f"    {src.name}：SITE_ID/VARIABLE 列が無く読めない")
        except Exception as e:
            print(f"    読み込み失敗 {type(e).__name__}: {str(e)[:100]}")

    # ── ② COSORE の寄与者情報（研究者単位の縮約に要る）──
    print("\n  ── ② COSORE の**寄与者**情報（研究者単位の縮約に要る）──")
    desc = pd.read_csv(Path(a.cosore_dir) / "description.csv")
    found = [c for c in WHO_COLS if c in desc.columns]
    print(f"    description.csv の列数 {len(desc.columns)}")
    print(f"    寄与者らしき列：{found if found else '**無し**'}")
    if not found:
        print(f"    参考・全列名：{list(desc.columns)}")
    for c in found:
        u = desc[c].astype(str).nunique()
        print(f"      {c}：{u} 通り／{len(desc)} データセット")
        if 1 < u <= 40:
            vc = desc[c].astype(str).value_counts()
            top = ", ".join(f"{k}×{v}" for k, v in vc.head(8).items())
            print(f"        上位：{top}")
    if found:
        print("    → **研究者単位のクラスタ縮約ができる**")
    else:
        print("    → **できない**。代案：`CSR_DATASET` の接頭辞（提供者名らしき部分）で近似する")
        pref = desc["CSR_DATASET"].astype(str).str.extract(r"^d\d+_([A-Za-z\-]+)")[0]
        print(f"       接頭辞で分けると **{pref.nunique()} 群**／{len(desc)} データセット")
        vc = pref.value_counts()
        print(f"       上位：{', '.join(f'{k}×{v}' for k, v in vc.head(10).items())}")
        print("       ＝**データセット名は提供者名を含む**ので、**近似としては使える**")
        print("         （**同姓・共同研究は区別できない**＝そう明記して使うこと）")

    # ── ③ ふるいの脱落（材料は COSORE のファイルだけで足りる）──
    print("\n  ── ③ ふるいの脱落を数える材料 ──")
    ds_dir = Path(a.cosore_dir) / "datasets"
    n_files = len(list(ds_dir.glob("data_*.csv"))) if ds_dir.exists() else 0
    print(f"    description の行数 {len(desc)}／datasets のファイル数 {n_files}")
    igbp = desc.get("CSR_IGBP")
    if igbp is not None:
        vc = igbp.astype(str).value_counts()
        print(f"    IGBP の内訳：{', '.join(f'{k}×{v}' for k, v in vc.head(12).items())}")
        print("    → **森林／非森林の脱落を別々に数えられる**")
    print("    → 材料は揃っている（各段階の脱落は旗84 で数える）")

    print("\n  === 次の判断 ===")
    print("  ①が可能なら：**★で選ばない総当たり**で、取得しなかった対の数と一致率を見積もる")
    print("    ＝**私が作ったバイアスの大きさを、数字で出せる**。")
    print("  ②が可能なら：**研究者単位で縮約**して、A-1 の 15/44 と A-3 の『クラスタ3』が残るか見る。")
    print("  ③は材料が揃っているので、**各段階の脱落と、落ちたサイトの性質**を数える。")
    print("  留保：")
    print("   ・**総当たりで見つかる対も、結局は『タワーがある場所』に限られる**＝")
    print("     **タワーの設置バイアス（平坦・均質・フェッチ良好）は消えない**。")
    print("   ・研究者単位の縮約は**共同研究や機材の共有を捉えない**＝**近似**である。")


if __name__ == "__main__":
    main()
