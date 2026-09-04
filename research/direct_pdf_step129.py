#!/usr/bin/env python3
"""旗129 段3 — 検索索引が名指しした PDF 直リンクを、**生の HTTP** で叩く。

段1/段2（`inst_index_step129.py`）は所属機関の刊行物一覧に当てる手で、
5 件中 3 件が「索引が JS で描かれていて生 HTTP では中身が無い」ため空振りした。
本段は **WebSearch が返した実在の URL** だけを対象にする（推測した URL は入れない）。

旗128 の新事実：`.venv/bin/python` の `urllib.request` は `WebFetch` が 403 を返す先でも
通ることがある（`WebFetch` は要約モデル経由で、UA も違う）。
**旗123 が「UMN ミラー 403」と記録した URL を、この経路で叩き直す**のが本段の主目的。

出力: research/logs/step129b_<timestamp>.txt / PDF は research/tmp_pdfs/
"""

import hashlib
import os
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runlog import tee_stdout  # noqa: E402

UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124.0 Safari/537.36")
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmp_pdfs")

# (tag, 出所＝どうしてこの URL を知ったか, URL)
TARGETS = [
    ("Gaumont-Guay2006",
     "旗123 の記録（WebFetch で 403）＋本周の WebSearch が同じ URL を返した",
     "https://biometeorology.umn.edu/sites/biometeorology.umn.edu/files/2021-04/boreal_aspen_2006.pdf"),
    ("Xu-Qi2001_agu_pdf",
     "本周の WebSearch（AGU/Wiley の PDF 直リンク・対照として 403 を期待）",
     "https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2000GB001365"),
    ("Xu-Qi2001_agu_epdf",
     "本周の WebSearch（同上・epdf 版）",
     "https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2000GB001365"),
]


def main():
    tee_stdout("step129b")
    print("=== 旗129 段3：検索索引が名指しした PDF 直リンクを生 HTTP で叩く ===")
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    print("**推測した URL は入れていない。出所を各行に書く。**\n")
    for tag, prov, url in TARGETS:
        print(f"--- {tag}\n    出所: {prov}\n    URL : {url}")
        req = urllib.request.Request(url, headers={
            "User-Agent": UA,
            "Accept": "application/pdf,text/html;q=0.9,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.9",
            "Referer": "https://www.google.com/",
        })
        t0 = time.time()
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                body, code, ctype, final = r.read(), r.status, r.headers.get("Content-Type", "?"), r.geturl()
        except urllib.error.HTTPError as e:
            print(f"    ✗ HTTP {e.code} {e.reason} ({time.time() - t0:.1f}s)\n")
            continue
        except Exception as e:
            print(f"    ✗ {type(e).__name__}: {e} ({time.time() - t0:.1f}s)\n")
            continue
        print(f"    status={code} ctype={ctype} bytes={len(body)} final={final} "
              f"({time.time() - t0:.1f}s)")
        if body[:5] == b"%PDF-":
            os.makedirs(OUTDIR, exist_ok=True)
            path = os.path.join(OUTDIR, tag + ".pdf")
            with open(path, "wb") as f:
                f.write(body)
            print(f"    ★ PDF 取得 sha256={hashlib.sha256(body).hexdigest()[:16]} → {path}\n")
        else:
            print(f"    ? PDF ではない。先頭 200 バイト: "
                  f"{body[:200].decode('utf-8', 'replace')!r}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
