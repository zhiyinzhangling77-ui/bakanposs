"""PDF の目次(見出し構造)を抽出する。

第一に PDF に埋め込まれたアウトライン(しおり)を使う。
無ければ Markdown 変換後の見出しから推定する(run 側で利用)。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class TocEntry:
    level: int
    title: str
    page: int  # 1 始まり


def extract_pdf_outline(pdf: Path) -> list[TocEntry]:
    """PyMuPDF で PDF 埋め込みアウトラインを取得。無ければ空リスト。"""
    try:
        import fitz  # PyMuPDF
    except Exception:
        return []
    try:
        doc = fitz.open(str(pdf))
    except Exception:
        return []
    entries: list[TocEntry] = []
    for level, title, page in doc.get_toc(simple=True):
        entries.append(TocEntry(level=level, title=title.strip(), page=page))
    doc.close()
    return entries


def top_level_chapters(entries: list[TocEntry]) -> list[TocEntry]:
    """アウトラインから最上位(章)レベルだけを抜き出す。"""
    if not entries:
        return []
    min_level = min(e.level for e in entries)
    return [e for e in entries if e.level == min_level]


def format_toc(entries: list[TocEntry], max_lines: int = 200) -> str:
    """目次を人間が読める形に整形。章には [番号] を付けて --chapters 指定に使える。"""
    if not entries:
        return "  (埋め込み目次なし — 変換後の見出しで章を判定します)"
    lines: list[str] = []
    min_level = min(e.level for e in entries)
    chap_no = 0
    for e in entries[:max_lines]:
        indent = "  " * (e.level - min_level)
        if e.level == min_level:
            chap_no += 1
            marker = f"[{chap_no:>2}] "
        else:
            marker = "     "
        lines.append(f"  {marker}{indent}{e.title}  (p.{e.page})")
    if len(entries) > max_lines:
        lines.append(f"  … 他 {len(entries) - max_lines} 項目")
    return "\n".join(lines)
