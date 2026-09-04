#!/usr/bin/env python3
"""旗130 段2b — Gaumont-Guay 2006 の wayback 捕獲物を、途中で切れた読みから取り直す。

1 走目は `IncompleteRead(130825 bytes read, 1138759 more expected)`。
**これは「本文が無い」ではなく「転送が途中で切れた」である**（欠陥 #40 の作法で分ける）。
Range ヘッダで続きを継ぎ、5 回まで試す。**それでも埋まらなければ `TRUNCATED` と書き、
「取れた」とは書かない。**
"""

import hashlib
import os
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runlog import tee_stdout  # noqa: E402

URL = ("https://web.archive.org/web/20080708193216id_/"
       "http://www.biometeorology.umn.edu/pdf/boreal_aspen_2006.pdf")
UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124.0 Safari/537.36")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "tmp_pdfs", "Gaumont-Guay2006.pdf")


def get(url, start=None):
    h = {"User-Agent": UA, "Accept": "application/pdf,*/*"}
    if start:
        h["Range"] = f"bytes={start}-"
    req = urllib.request.Request(url, headers=h)
    with urllib.request.urlopen(req, timeout=120) as r:
        try:
            return r.read(), r.status, None
        except Exception as e:      # IncompleteRead は部分データを持っている
            part = getattr(e, "partial", b"")
            return part, r.status, f"{type(e).__name__}: {e}"


def main():
    tee_stdout("step130_retry")
    print("=== 旗130 段2b — Gaumont-Guay 2006 を Range で継ぐ ===")
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    buf = b""
    for attempt in range(1, 6):
        try:
            chunk, status, why = get(URL, start=len(buf) if buf else None)
        except urllib.error.HTTPError as e:
            print(f"  試行{attempt}: HTTPError {e.code} {e.reason}")
            time.sleep(2)
            continue
        except Exception as e:
            print(f"  試行{attempt}: {type(e).__name__}: {e}")
            time.sleep(2)
            continue
        buf += chunk
        print(f"  試行{attempt}: status={status} +{len(chunk)}B 累計={len(buf)}B "
              f"{'切断: ' + why if why else '完了'}")
        if why is None and chunk:
            break
        if not chunk:
            print("    → これ以上増えない")
            break
        time.sleep(1)
    if not buf:
        print("  判定: NO_ROUTE（1 バイトも取れない）")
        return 0
    if buf[:5] != b"%PDF-":
        print(f"  判定: NOT_PDF（先頭 {buf[:20]!r}）")
        return 0
    tail_ok = b"%%EOF" in buf[-2048:]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "wb") as f:
        f.write(buf)
    print(f"  bytes={len(buf)} sha256={hashlib.sha256(buf).hexdigest()[:16]} "
          f"末尾に %%EOF: {tail_ok} → {OUT}")
    print("  判定: " + ("GOT_PDF" if tail_ok else "TRUNCATED（読める範囲だけ・完本ではない）"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
