"""Markdown を日本語対応 PDF に変換する (markdown → HTML → headless Chromium → PDF)。

pandoc/LaTeX 無し・IPAGothic フォント前提の環境向け。表・コードブロック対応。
Chromium が日本語フォントをネイティブに扱うので CJK が確実に出る。

    python scripts/md_to_pdf.py japanflux_pn/UNDERSTANDING.md [...]

各 <name>.md を <name>.pdf に変換して同じディレクトリに置く。
"""

from __future__ import annotations

import glob
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import markdown

CSS = """
@page { size: A4; margin: 18mm 15mm; }
html { -webkit-print-color-adjust: exact; }
* { font-family: "IPAGothic","IPAPGothic","Noto Sans CJK JP",sans-serif; }
body { font-size: 10.5pt; line-height: 1.65; color: #111; }
h1 { font-size: 20pt; border-bottom: 2px solid #333; padding-bottom: 5px; }
h2 { font-size: 15pt; border-bottom: 1px solid #999; padding-bottom: 3px; margin-top: 1.2em; }
h3 { font-size: 12.5pt; margin-top: 1em; }
code { background: #f0f0f0; padding: 1px 3px; border-radius: 3px; font-size: 9.5pt; }
pre { background: #f5f5f5; border: 1px solid #ddd; padding: 8px; border-radius: 4px;
      font-size: 9pt; white-space: pre-wrap; word-wrap: break-word; }
pre code { background: none; padding: 0; }
table { border-collapse: collapse; width: 100%; margin: 8px 0; font-size: 9.3pt;
        page-break-inside: avoid; }
th, td { border: 1px solid #bbb; padding: 4px 7px; text-align: left; vertical-align: top; }
th { background: #ececec; }
blockquote { border-left: 3px solid #bbb; margin: 8px 0; padding: 2px 12px; color: #444; }
h1, h2, h3 { page-break-after: avoid; }
"""

HTML_TMPL = ('<!doctype html><html lang="ja"><head><meta charset="utf-8">'
             '<style>{css}</style></head><body>{body}</body></html>')


def _find_chrome() -> str:
    for pat in ("/opt/pw-browsers/chromium-*/chrome-linux/chrome",
                "/opt/pw-browsers/chromium_headless_shell-*/chrome-linux/headless_shell"):
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[-1]
    for name in ("chromium", "chromium-browser", "google-chrome"):
        from shutil import which
        p = which(name)
        if p:
            return p
    raise RuntimeError("Chromium が見つかりません")


def convert(md_path: Path, chrome: str) -> Path:
    body = markdown.markdown(
        md_path.read_text(encoding="utf-8"),
        extensions=["tables", "fenced_code", "sane_lists"],
    )
    html = HTML_TMPL.format(css=CSS, body=body)
    pdf_path = md_path.with_suffix(".pdf")
    with tempfile.TemporaryDirectory() as td:
        html_path = Path(td) / (md_path.stem + ".html")
        html_path.write_text(html, encoding="utf-8")
        subprocess.run(
            [chrome, "--headless", "--no-sandbox", "--disable-gpu",
             "--no-pdf-header-footer", f"--print-to-pdf={pdf_path}",
             html_path.as_uri()],
            check=True, capture_output=True, timeout=180,
            env={**os.environ, "HOME": td},
        )
    if not pdf_path.exists():
        raise RuntimeError(f"PDF 生成失敗: {md_path}")
    return pdf_path


def main(argv: list[str]) -> None:
    if not argv:
        print("usage: md_to_pdf.py FILE.md [FILE2.md ...]")
        raise SystemExit(1)
    chrome = _find_chrome()
    for a in argv:
        out = convert(Path(a), chrome)
        print(f"[pdf] {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main(sys.argv[1:])
