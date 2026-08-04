#!/usr/bin/env python3
"""Unpaywall 経由で「合法的にオープンアクセスな」PDFだけを取得する。

- 有料(ペイウォール)論文は取得しない → "要・機関アクセス" として一覧に出す。
- ペイウォール回避は一切しない(合法な OA コピーのみ)。有料のものは Zotero＋
  大学プロキシ等の正規ルートで取得してください。

使い方:
  python fetch_oa_pdfs.py dois.txt your_email@example.com [出力フォルダ]

dois.txt: 1行1 DOI(# で始まる行はコメント)。DOI は事前に doi.org で検証すること。
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

UA = "Mozilla/5.0 (research-oa-fetcher; mailto:%s)"


def _get(url: str, email: str, binary: bool = False):
    req = urllib.request.Request(url, headers={"User-Agent": UA % email})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read() if binary else json.load(r)


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    dois_file, email = sys.argv[1], sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else "oa_pdfs"
    os.makedirs(outdir, exist_ok=True)

    dois = [
        l.strip() for l in open(dois_file, encoding="utf-8")
        if l.strip() and not l.lstrip().startswith("#")
    ]
    got: list[str] = []
    paywalled: list[str] = []
    failed: list[tuple[str, str]] = []

    for doi in dois:
        try:
            api = (f"https://api.unpaywall.org/v2/{urllib.parse.quote(doi)}"
                   f"?email={urllib.parse.quote(email)}")
            data = _get(api, email)
            loc = data.get("best_oa_location") or {}
            pdf_url = loc.get("url_for_pdf")
            if pdf_url:
                safe = re.sub(r"[^\w.-]", "_", doi) + ".pdf"
                dest = os.path.join(outdir, safe)
                content = _get(pdf_url, email, binary=True)
                with open(dest, "wb") as f:
                    f.write(content)
                got.append(doi)
                print(f"OA      {doi}  -> {dest}")
            else:
                paywalled.append(doi)
                print(f"PAYWALL {doi}  (要・機関アクセス / Zotero で)")
        except Exception as e:  # noqa: BLE001
            failed.append((doi, str(e)[:80]))
            print(f"ERROR   {doi}  {e}")
        time.sleep(1)  # Unpaywall へのマナー

    print(f"\n=== 合計: OA取得 {len(got)} / 有料 {len(paywalled)} / エラー {len(failed)} ===")
    if paywalled:
        print("\n[要・機関アクセス(Zotero＋大学プロキシで)]")
        for d in paywalled:
            print(f"  https://doi.org/{d}")
    if failed:
        print("\n[エラー(DOI検証を推奨)]")
        for d, e in failed:
            print(f"  {d}  {e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
