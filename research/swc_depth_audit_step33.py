"""旗33：各サイトの土壌水分(SWC)・地温(TS)の"深度列"を棚卸しする（比較の交絡確認）。

Perplexity/文献（Pastorello 2020, Wang 2025 ほか）で確認：SWC の測定深度はサイトで異なり
標準深度は無い。表層(0–10cm)=微生物分解／根圏(10–30cm)=根呼吸と、捕まえるプロセスが違うので、
深度が揃わないと生態系間の「θ→呼吸」比較が歪む。我々は `SWC_F_MDS`（層インデックス無し＝単一深度）
を使っているが、その深度がサイトで揃う保証はない。

このツールは各サイトの生 CSV ヘッダを読み、SWC*/TS* の列（層インデックス _1,_2… の有無）を
棚卸しする。複数深度を持つサイトは「1つしか使っていない」ことが分かり、単一 `SWC_F_MDS` の
サイトは深度がヘッダからは不明（BADM メタデータが要る）と明示する。値は読まない（ヘッダのみ＝高速）。

    python research/swc_depth_audit_step33.py --sites JP-Tak CN-HaM MN-Kbu JP-BBY JP-Mse KR-CRK
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SWC_RE = re.compile(r"^SWC", re.IGNORECASE)
TS_RE = re.compile(r"^TS", re.IGNORECASE)
LAYER_RE = re.compile(r"_(\d+)$")     # 末尾の層インデックス


def _layers(cols):
    """層インデックス（_1,_2…）の集合。無ければ空。"""
    ls = set()
    for c in cols:
        m = LAYER_RE.search(c)
        if m:
            ls.add(int(m.group(1)))
    return sorted(ls)


def audit_site(site):
    from japanflux_pn.preprocess import _read_table_header, find_corevars_files
    files = find_corevars_files(site)
    header = _read_table_header(files[0])
    swc = [c for c in header if SWC_RE.match(c) and "QC" not in c]
    ts = [c for c in header if TS_RE.match(c) and "QC" not in c and "TSTAMP" not in c.upper()]
    used = "SWC_F_MDS" if "SWC_F_MDS" in header else \
           ("SWC_F_MDS_1" if "SWC_F_MDS_1" in header else (swc[0] if swc else "—"))
    return {"swc": swc, "ts": ts, "swc_layers": _layers(swc),
            "ts_layers": _layers(ts), "used": used}


def main():
    from japanflux_pn.sites import get_site
    p = argparse.ArgumentParser(description="SWC/TS の深度列を棚卸し")
    p.add_argument("--sites", nargs="+", required=True)
    a = p.parse_args()

    print("=== 旗33 SWC/TS の深度列 棚卸し（θ→呼吸の深度交絡チェック）===")
    print("  used=解析で使う SWC 列／層=末尾 _1,_2… の有無（複数=別深度あり, 単一 SWC_F_MDS=深度不明）\n")
    print(f"  {'サイト':<8} {'使用SWC列':<14} {'SWC層':>8} {'TS層':>8}  SWC列一覧")
    multi = single = 0
    for s in a.sites:
        try:
            r = audit_site(get_site(s))
        except Exception as e:
            print(f"  {s:<8} SKIP {type(e).__name__}: {e}"); continue
        lay = r["swc_layers"]
        tag = f"{len(lay)}層{lay}" if lay else "単一(層なし)"
        if lay and len(lay) > 1:
            multi += 1
        else:
            single += 1
        print(f"  {s:<8} {r['used']:<14} {tag:>8} "
              f"{(str(r['ts_layers']) if r['ts_layers'] else '単一'):>8}  {r['swc'][:4]}")
    print(f"\n  複数深度あり={multi} / 単一深度={single}")
    print("  読み方：単一 SWC_F_MDS のサイトは深度がヘッダからは不明＝BADM(Height/Depth)が要る。")
    print("    複数層のサイトは 1 深度しか使っておらず、サイト間で深度が揃う保証はない。")
    print("  ＝旗31 の θ→呼吸 比較には除去できない深度交絡（表層=微生物 vs 根圏=根呼吸）がある。")
    print("    符号の質的差（草原＋ vs 湿地−）は深度で反転しにくいが、順位・絶対値は深度依存＝留保。")


if __name__ == "__main__":
    main()
