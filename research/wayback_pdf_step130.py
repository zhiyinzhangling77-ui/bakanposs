#!/usr/bin/env python3
"""旗130 段2 — scholar.archive.org が出した全文リンクを実際に取りに行く。

段1（WebFetch で scholar.archive.org の検索ページを読む）で、⚪ 5 件のうち 2 件に
全文リンクが付いていた。段2 はその本体を取り、**PDF であることと中身の語**まで確かめる。

**旗129 の申し送り**：生 HTTP は API には通るがファイルサーバ（Wiley・UMN の Drupal）には 403。
**本周の当ては「wayback の捕獲物なら出版社/研究室のサーバを経由しない」である。**
当たっても外れても、**落ちた段を分けて記録する**（欠陥 #40）。

対照（門①）：本文まで一次到達した実績のある 2 件を、**同じ経路で同じように**取る。
  取れなければ「経路が悪い」、取れれば「その ⚪ に本文が無い」と分けられる。

出力: research/logs/step130_wayback_<timestamp>.txt / PDF は research/tmp_pdfs/（版管理しない）
"""

import hashlib
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runlog import tee_stdout  # noqa: E402

UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124.0 Safari/537.36")
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmp_pdfs")

SA = "https://scholar.archive.org"

# (tag, 種別, 期待される本文の語, [候補 URL を上から順に])
# 種別 WHITE = 本文未到達の ⚪ / CTRL = 本文まで一次到達済みの対照
# **`id_` を付けると wayback は中継 HTML ではなく捕獲した原本を返す**（1 走目の欠陥 #55）。
CANDIDATES = [
    ("Gaumont-Guay2006", "WHITE", ["aspen", "soil respiration", "water content"], [
        "https://web.archive.org/web/20080708193216id_/"
        "http://www.biometeorology.umn.edu/pdf/boreal_aspen_2006.pdf",
        SA + "/work/janteh56h5aqfdjjr2dqqb2rbq/access/wayback/"
             "http://www.biometeorology.umn.edu/pdf/boreal_aspen_2006.pdf",
    ]),
    ("Suseela2012", "WHITE", ["heterotrophic", "old-field", "temperature sensitivity"], [
        "https://web.archive.org/web/20171211163614id_/"
        "https://dge.carnegiescience.edu/DGE/Dukes/SuseelaEtAl2012.pdf",
        "https://dge.carnegiescience.edu/DGE/Dukes/SuseelaEtAl2012.pdf",
    ]),
    # --- 門①（対照）: 同じ経路で、本文到達済みの 2 件を取る ---
    ("CTRL_Davidson1998", "CTRL", ["confounded", "soil respiration"], [
        "https://web.archive.org/web/2019id_/http://pdfs.semanticscholar.org/f72a/"
        "8e9f706df71f0b7ce7750157f43792a2e66b.pdf",
        "http://pdfs.semanticscholar.org/f72a/"
        "8e9f706df71f0b7ce7750157f43792a2e66b.pdf",
    ]),
    ("CTRL_Wen2006", "CTRL", ["soil moisture", "ecosystem respiration", "pinus"], [
        # 段1 の対照検索が「Archived PDF (459 kB, captured 2017)」を出した項目。
        # 直リンクは段1 の要約に無いので、scholar.archive.org 経由の 1 本だけを試す。
        SA + "/search?q=%22Soil+moisture+effect+on+the+temperature+dependence+of+"
             "ecosystem+respiration%22",
    ]),
]


def fetch(url, timeout=60):
    req = urllib.request.Request(url, headers={
        "User-Agent": UA,
        "Accept": "application/pdf,text/html,*/*",
        "Accept-Language": "en-US,en;q=0.9",
    })
    t0 = time.time()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return {"ok": True, "status": r.status, "body": r.read(),
                    "ctype": r.headers.get("Content-Type", "?"),
                    "final": r.geturl(), "dt": round(time.time() - t0, 1)}
    except urllib.error.HTTPError as e:
        return {"ok": False, "status": e.code, "why": f"HTTPError {e.reason}",
                "dt": round(time.time() - t0, 1)}
    except Exception as e:
        return {"ok": False, "status": None, "why": f"{type(e).__name__}: {e}",
                "dt": round(time.time() - t0, 1)}


def pdf_text(path):
    """pdftotext があれば本文の語を確かめる。無ければ None（**「確かめた」と書かない**）。"""
    try:
        out = subprocess.run(["pdftotext", "-l", "4", path, "-"],
                             capture_output=True, timeout=120)
        if out.returncode == 0:
            return out.stdout.decode("utf-8", "replace")
    except Exception as e:
        print(f"    (pdftotext 不可: {type(e).__name__}: {e})")
    return None


def run(tag, kind, words, urls):
    print(f"\n{'=' * 78}\n### {tag} [{kind}]")
    for url in urls:
        r = fetch(url)
        if not r["ok"]:
            print(f"  ✗ {url}\n    status={r['status']} {r['why']} ({r['dt']}s)")
            continue
        body = r["body"]
        print(f"  ✓ {url}\n    status={r['status']} ctype={r['ctype']} "
              f"bytes={len(body)} final={r['final'][:120]} ({r['dt']}s)")
        if body[:5] != b"%PDF-":
            head = body[:200].decode("utf-8", "replace").replace("\n", " ")
            print(f"    → PDF ではない。先頭: {head}")
            continue
        os.makedirs(OUTDIR, exist_ok=True)
        path = os.path.join(OUTDIR, tag + ".pdf")
        with open(path, "wb") as f:
            f.write(body)
        sha = hashlib.sha256(body).hexdigest()[:16]
        print(f"    ★ PDF 取得 bytes={len(body)} sha256={sha} → {path}")
        txt = pdf_text(path)
        if txt is None:
            print("    判定: GOT_PDF_UNVERIFIED（PDF は取れたが中身の語を確かめていない）")
            return "GOT_PDF_UNVERIFIED"
        low = txt.lower()
        hits = [w for w in words if w.lower() in low]
        print(f"    本文の語 {hits} / {words}  （抽出 {len(txt)} 文字）")
        print("    先頭 400 字: " + re.sub(r"\s+", " ", txt[:400]))
        return "GOT_PDF" if len(hits) >= 2 else "GOT_PDF_WRONG_CONTENT"
    return "NO_ROUTE"


def main():
    tee_stdout("step130_wayback")
    print("=== 旗130 段2 — scholar.archive.org が出した全文リンクを取りに行く ===")
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    v = {}
    for tag, kind, words, urls in CANDIDATES:
        try:
            v[tag] = run(tag, kind, words, urls)
        except Exception as e:
            print(f"  !! 例外 {type(e).__name__}: {e}")
            v[tag] = f"ERROR({type(e).__name__})"
    print(f"\n{'=' * 78}\n### まとめ")
    for k, s in v.items():
        print(f"  {k:22s} {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
