#!/usr/bin/env python3
"""旗130 補助 — 取得した PDF の本文テキストを書き出す（読むのは人／エージェント）。

使い方: .venv/bin/python research/pdf_read_step130.py <tag>
  research/tmp_pdfs/<tag>.pdf → research/tmp_pdfs/<tag>.txt
本文そのものは生成物なので版管理しない。**要約に頼らず、この .txt を読んで引用する。**
"""

import os
import subprocess
import sys

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmp_pdfs")


def main():
    tag = sys.argv[1]
    src, dst = os.path.join(D, tag + ".pdf"), os.path.join(D, tag + ".txt")
    out = subprocess.run(["pdftotext", src, dst], capture_output=True, timeout=180)
    if out.returncode != 0:
        print("pdftotext 失敗:", out.stderr.decode("utf-8", "replace"))
        return 1
    n = os.path.getsize(dst)
    with open(dst, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    print(f"{dst}  {n} bytes  {len(lines)} 行")
    return 0


if __name__ == "__main__":
    sys.exit(main())
