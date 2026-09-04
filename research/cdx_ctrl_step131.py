#!/usr/bin/env python3
"""旗131 段2 — **門①の陽性対照だけを引き直す。**

1 走目は CDX API が 503 を返し、**陽性対照が「取れない」ではなく「情報が無い」で落ちた**。
**この状態で白（Gaumont-Guay 2006）の成功を主張すると、門①を省いたことになる。**
CDX は混み合うと 503 を返すので、間を置いて数回試す。**それでも 503 のままなら、
`CDX_FAIL` と書き、門①は満たされなかったと正直に記す。**
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cdx_fulltext_step131 as M  # noqa: E402
from runlog import tee_stdout  # noqa: E402


def main():
    tee_stdout("step131_ctrl")
    print("=== 旗131 段2 — 門①の陽性対照（Davidson 1998）を引き直す ===")
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    name, kind, orig, words = M.TARGETS[1]
    print(f"  元 URL: {orig}")
    caps = None
    for attempt in range(1, 6):
        caps, note = M.cdx_list(orig)
        print(f"  CDX 試行{attempt}: {note if caps is None else str(len(caps)) + ' 件'}")
        if caps:
            break
        time.sleep(10)
    if not caps:
        print("  判定: CDX_FAIL のまま —— **門①は満たされていない**")
        return 0
    for c in caps[:8]:
        print(f"    {c['ts']}  {c['mime']:<20} len={c['len']}")
    for i, c in enumerate(caps, 1):
        url, buf, dec, err = M.fetch_capture(c["ts"], c["orig"])
        v = M.verdict(buf, dec) if not err else f"ERR({err})"
        print(f"  [{i}/{len(caps)}] {c['ts']} → {len(buf)}B 宣言={dec} 判定={v}")
        if v != "GOT_PDF":
            time.sleep(2)
            continue
        path = os.path.join(M.OUTDIR, f"{name}.pdf")
        with open(path, "wb") as f:
            f.write(buf)
        txt, hit = M.pdf_words(path, words)
        print(f"    ★ 完本 {len(buf)}B / 本文の語 {hit} / {words}")
        print("    先頭 200 字: " + " ".join((txt or "")[:200].split()))
        print("  判定: GOT_PDF —— **門①（陽性）は通った**")
        return 0
    print("  判定: ALL_CAPTURES_FAILED —— **門①は通らなかった**")
    return 0


if __name__ == "__main__":
    sys.exit(main())
