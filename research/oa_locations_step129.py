#!/usr/bin/env python3
"""旗129 段5 — **まだ当てていない索引 2 つ**を、対照つきで全 ⚪ に当てる。

旗128 が見たのは Unpaywall の `is_oa` と OpenAlex の `best_oa_location` だけだった。
本段が足すのは：

  (1) **OpenAlex の `locations[]` 全部** — `best_oa_location` が null でも、
      `locations` にリポジトリ版が並ぶことがある（`best_` は「出版社版に最も近い 1 つ」）。
  (2) **Semantic Scholar の `openAccessPdf`** — Unpaywall/OpenAlex とは別に集めた索引。
      旗121〜127 は S2 を**抄録の取得にしか使っていなかった**。

**旗128 の作法どおり、既知の陽性（本文まで一次到達した 5 件）に先に当てる。**
感度が 0 なら、⚪ 側の空欄は到達不能の根拠にならない。

出力: research/logs/step129d_<timestamp>.txt
"""

import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runlog import tee_stdout  # noqa: E402

MAIL = "zhiyinzhangling77@gmail.com"
UA = f"bakanposs-research/1.0 (mailto:{MAIL})"

CONTROLS = [
    ("Zhang2018 AFM 259 (旗121 到達)", "10.1016/j.agrformet.2018.05.005"),
    ("Wang2014 BG 11 (旗122 到達)", "10.5194/bg-11-259-2014"),
    ("Wen2006 AFM 137 (旗123 到達)", "10.1016/j.agrformet.2006.02.005"),
    ("Davidson1998 GCB 4 (旗124 到達)", "10.1046/j.1365-2486.1998.00128.x"),
    ("Doerr&Muennich1987 Tellus 39B (旗127 到達)", "10.3402/tellusb.v39i1-2.15331"),
]
OPEN = [
    ("Xu&Qi2001 GBC 15 687-696", "10.1029/2000GB001365"),
    ("Gaumont-Guay2006 AFM 140 220-235", "10.1016/j.agrformet.2006.08.003"),
    ("Tucker&Reed2016 Biogeochem 128 155-169", "10.1007/s10533-016-0200-1"),
    ("Suseela2012 GCB 18 336-348", "10.1111/j.1365-2486.2011.02516.x"),
    ("Bunnell1977 SBB 9 33-40", "10.1016/0038-0717(77)90058-X"),
]


def get_json(url, timeout=45):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode("utf-8", "replace")), None
    except urllib.error.HTTPError as e:
        return None, f"HTTP {e.code} {e.reason}"
    except Exception as e:
        return None, f"{type(e).__name__}: {e}"


def openalex_locations(doi):
    d, err = get_json(f"https://api.openalex.org/works/doi:{urllib.parse.quote(doi)}"
                      f"?mailto={MAIL}")
    if err:
        return None, err
    out = []
    for loc in (d.get("locations") or []):
        src = (loc.get("source") or {}).get("display_name")
        out.append({"is_oa": loc.get("is_oa"), "type": loc.get("version"),
                    "source": src, "pdf": loc.get("pdf_url"),
                    "landing": loc.get("landing_page_url")})
    return out, None


def s2_oa(doi):
    d, err = get_json("https://api.semanticscholar.org/graph/v1/paper/DOI:"
                      f"{urllib.parse.quote(doi)}?fields=title,isOpenAccess,openAccessPdf,externalIds")
    if err:
        return None, err
    return {"isOpenAccess": d.get("isOpenAccess"), "openAccessPdf": d.get("openAccessPdf"),
            "pmid": (d.get("externalIds") or {}).get("PubMed")}, None


def run(items, header):
    print(f"\n--- {header} ---")
    verd = {}
    for label, doi in items:
        print(f"\n  {label}\n    doi={doi}")
        locs, err = openalex_locations(doi)
        pdfs = []
        if err:
            print(f"    [OpenAlex] ✗ {err}")
        else:
            print(f"    [OpenAlex] locations={len(locs)}")
            for L in locs:
                mark = "PDF" if L["pdf"] else "   "
                print(f"       {mark} is_oa={L['is_oa']!s:5s} ver={L['type']!s:12s} "
                      f"src={L['source']}")
                if L["pdf"]:
                    print(f"           pdf: {L['pdf']}")
                    pdfs.append(L["pdf"])
        time.sleep(0.4)
        s2, err2 = s2_oa(doi)
        if err2:
            print(f"    [S2]       ✗ {err2}")
        else:
            oap = s2.get("openAccessPdf") or {}
            print(f"    [S2]       isOpenAccess={s2['isOpenAccess']} "
                  f"pdf={oap.get('url')} pmid={s2.get('pmid')}")
            if oap.get("url"):
                pdfs.append(oap["url"])
        time.sleep(0.4)
        verd[label] = ("PDF_URL" if pdfs else "NO_PDF_URL")
        if pdfs:
            print(f"    ★ 取りにいける PDF URL {len(pdfs)} 本")
    return verd


def main():
    tee_stdout("step129d")
    print("=== 旗129 段5：OpenAlex `locations[]` 全部 ＋ S2 `openAccessPdf` ===")
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    print("**対照を先に通す。感度 0 なら ⚪ の空欄は情報を持たない**（旗128 の読み方）。")
    c = run(CONTROLS, "対照（本文まで一次到達した実績あり）")
    o = run(OPEN, "残る ⚪ 5 件")
    hit = sum(1 for v in c.values() if v == "PDF_URL")
    print(f"\n=== まとめ ===\n  感度: {hit}/{len(CONTROLS)}")
    for k, v in c.items():
        print(f"    [対照] {v:12s} {k}")
    for k, v in o.items():
        print(f"    [⚪]   {v:12s} {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
