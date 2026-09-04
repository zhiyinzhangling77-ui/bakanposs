#!/usr/bin/env python3
"""旗132 段4: 到達性スクリーンに **OpenAIRE** を足せるか（旗129 の三索引の和の続き）。

旗129 の和 = Unpaywall ∪ OpenAlex.locations[] ∪ S2.openAccessPdf（感度 3/5）。
OpenAIRE は欧州のリポジトリ収集に強く、上の 3 つが見ていない機関リポジトリを持つことがある。
**門①**: 本文到達済みの陽性（Suseela 2012・Gaumont-Guay 2006）と、
到達不能で確定した陰性（Wildung 1975）に同じ問い合わせを当てる。
陽性で何も返らないなら、この索引は ⚪ の判断に使えない。
"""
import re
import sys
import urllib.parse
import urllib.request
import urllib.error

UA = "Mozilla/5.0 (X11; Linux x86_64) research-primary-check/1.0"

CASES = [
    ("TARGET_tucker2016", "10.1007/s10533-016-0200-1"),
    ("CTRL_POS_suseela2012", "10.1111/j.1365-2486.2011.02516.x"),
    ("CTRL_POS_gaumontguay2006", "10.1016/j.agrformet.2006.08.003"),
    ("CTRL_NEG_fakedoi", "10.9999/this-doi-does-not-exist-9999"),
]


def http(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, (e.read()[:800] if e.fp else b"")
    except Exception as e:
        return None, repr(e).encode()


for tag, doi in CASES:
    url = "https://api.openaire.eu/search/publications?doi=%s&format=json" % urllib.parse.quote(doi)
    st, body = http(url)
    print("=" * 74)
    print("[%s] doi=%s" % (tag, doi))
    print("  http:", st, "bytes:", len(body))
    if st != 200:
        print("  body head:", body[:300])
        continue
    txt = body.decode("utf-8", "replace")
    # 生 JSON の構造に依存せず、fulltext/PDF らしい URL を全部拾う（要約モデルを介さない）
    urls = sorted(set(re.findall(r'https?://[^"\s\\<>]+', txt)))
    pdfish = [u for u in urls if ".pdf" in u.lower() or "/download" in u.lower()
              or "fulltext" in u.lower() or "bitstream" in u.lower()]
    print("  総 URL:", len(urls), "/ PDF らしい:", len(pdfish))
    for u in pdfish[:15]:
        print("    *", u[:150])
    if not pdfish:
        for u in urls[:12]:
            print("    -", u[:150])
sys.stdout.flush()
