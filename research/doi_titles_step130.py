#!/usr/bin/env python3
"""旗130 段0 — ⚪ 5 件と対照の書誌を DOI から確定する（欠陥 #54 の作法）。

旗129 で「⚪ の題名を DOI で確かめずに前の周のクエリ文字列を引き継ぎ、別論文を探していた」
（欠陥 #54）ため、索引を引く前に必ずここを通す。Crossref API のみ（生 HTTP で通る＝旗129 段5）。

出力: research/logs/step130_titles_<timestamp>.txt
"""

import json
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runlog import tee_stdout  # noqa: E402

# ⚪（本文未到達）と、対照（本文まで一次到達した実績のあるもの）を**同じ手で**引く。
DOIS = [
    ("WHITE Xu-Qi2001", "10.1029/2000GB001365"),
    ("WHITE Gaumont-Guay2006", "10.1016/j.agrformet.2006.08.003"),
    ("WHITE Tucker-Reed2016", "10.1007/s10533-016-0200-1"),
    ("WHITE Suseela2012", "10.1111/j.1365-2486.2011.02516.x"),
    ("WHITE Bunnell1977", "10.1016/0038-0717(77)90058-X"),
    ("CTRL Davidson1998", "10.1046/j.1365-2486.1998.00128.x"),
    ("CTRL Wang2014", "10.5194/bg-11-6205-2014"),
    ("CTRL Zhang2018", "10.1016/j.agrformet.2018.05.005"),
]

UA = "japanflux-pn-loop/1.0 (mailto:zhiyinzhangling77@gmail.com)"


def one(doi):
    req = urllib.request.Request(
        "https://api.crossref.org/works/" + doi, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)["message"]


def main():
    tee_stdout("step130_titles")
    print("=== 旗130 段0 — DOI から書誌を確定（Crossref）===")
    for tag, doi in DOIS:
        print(f"\n### {tag}  doi:{doi}")
        try:
            m = one(doi)
        except Exception as e:
            print(f"  !! {type(e).__name__}: {e}")
            continue
        title = (m.get("title") or ["?"])[0]
        cont = (m.get("container-title") or ["?"])[0]
        year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
        authors = "; ".join(
            f"{a.get('family', '?')} {a.get('given', '')}".strip()
            for a in m.get("author", []))
        print(f"  題名: {title}")
        print(f"  誌名: {cont} ({year})  vol={m.get('volume')} page={m.get('page')}")
        print(f"  著者: {authors}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
