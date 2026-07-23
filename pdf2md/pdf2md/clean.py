"""Markdown のノイズ除去(トークン削減用)。

除去するもの:
  - ページ番号だけの行 / ヘッダ・フッタの繰り返し行
  - 索引(Index / 索引)セクション
  - 参考文献の羅列(References / Bibliography / 参考文献)セクション

必ず残すもの:
  - 数式(LaTeX: $...$, $$...$$, \\(...\\), \\[...\\])
  - 表(Markdown table: | ... |)
  - コードブロック(``` ... ```)
"""

from __future__ import annotations

import re
from dataclasses import dataclass

HEADING_RE = re.compile(r"^(#{1,6})\s+(.*\S)\s*$")
# 見出し行の Markdown 強調記号(** * __ `)を除去するための正規表現
_EMPHASIS_RE = re.compile(r"(\*\*|__|\*|`)")

# セクション見出しが「ノイズ節」かどうか(索引・参考文献)
NOISE_SECTION_RE = re.compile(
    r"^(index|indices|bibliography|references?|works\s+cited|"
    r"索\s*引|参\s*考\s*文\s*献|引\s*用\s*文\s*献|文\s*献)\s*$",
    re.IGNORECASE,
)

# ページ番号だけの行(「12」「- 12 -」「p. 12」「12 / 340」など)
PAGE_NUM_RE = re.compile(
    r"^\s*(?:[-–—]\s*)?(?:page\s+|p\.?\s*)?\d{1,4}(?:\s*[-–—/]\s*\d{1,4})?\s*(?:[-–—])?\s*$",
    re.IGNORECASE,
)

MULTI_BLANK_RE = re.compile(r"\n{3,}")


@dataclass
class CleanStats:
    removed_page_numbers: int = 0
    removed_repeated: int = 0
    dropped_sections: int = 0


def _is_protected(line: str) -> bool:
    """表・数式・箇条書きなど、消してはいけない行か。"""
    s = line.strip()
    if not s:
        return False
    if s.startswith("|") or s.startswith("$") or s.startswith("\\["):
        return True
    if s.startswith(("#", ">", "-", "*", "+")):
        return True
    if "$" in s or "|" in s:  # インライン数式・表セル
        return True
    return False


def _drop_noise_sections(text: str, stats: CleanStats) -> str:
    """索引・参考文献の見出しから、同レベル以上の次の見出しまでを丸ごと落とす。"""
    lines = text.splitlines()
    out: list[str] = []
    drop_until_level: int | None = None
    in_code = False

    for line in lines:
        if line.lstrip().startswith("```"):
            in_code = not in_code
            out.append(line)
            continue
        if in_code:
            out.append(line)
            continue

        m = HEADING_RE.match(line)
        if m:
            level = len(m.group(1))
            title = m.group(2).strip()
            if drop_until_level is not None:
                # ドロップ中: 同レベル以上の見出しが来たらドロップ終了
                if level <= drop_until_level:
                    drop_until_level = None
                else:
                    continue  # ノイズ節の下位見出しも捨てる
            if drop_until_level is None and NOISE_SECTION_RE.match(title):
                drop_until_level = level
                stats.dropped_sections += 1
                continue
        if drop_until_level is not None:
            continue
        out.append(line)
    return "\n".join(out)


def _remove_page_numbers(text: str, stats: CleanStats) -> str:
    out: list[str] = []
    in_code = False
    for line in text.splitlines():
        if line.lstrip().startswith("```"):
            in_code = not in_code
            out.append(line)
            continue
        s = line.strip()
        # 表・数式・見出しは対象外。それ以外でページ番号パターンなら除去
        # (PAGE_NUM_RE は数字と区切りのみに一致するので "- 340 -" も拾える)
        is_tableormath = ("|" in s) or ("$" in s) or s.startswith("#")
        if not in_code and not is_tableormath and PAGE_NUM_RE.match(line):
            stats.removed_page_numbers += 1
            continue
        out.append(line)
    return "\n".join(out)


def _remove_repeated_headers(
    text: str, stats: CleanStats, min_count: int = 6, max_len: int = 60
) -> str:
    """全体で何度も現れる短い行(=ヘッダ/フッタ)を除去。"""
    lines = text.splitlines()
    freq: dict[str, int] = {}
    for line in lines:
        s = line.strip()
        if s and len(s) <= max_len and not _is_protected(line):
            freq[s] = freq.get(s, 0) + 1
    noisy = {s for s, c in freq.items() if c >= min_count}
    if not noisy:
        return text
    out: list[str] = []
    in_code = False
    for line in lines:
        if line.lstrip().startswith("```"):
            in_code = not in_code
            out.append(line)
            continue
        if not in_code and line.strip() in noisy and not _is_protected(line):
            stats.removed_repeated += 1
            continue
        out.append(line)
    return "\n".join(out)


def _normalize_headings(text: str) -> str:
    """見出し行(#〜######)から強調記号を除去。`## **X**` → `## X`。"""
    out: list[str] = []
    in_code = False
    for line in text.splitlines():
        if line.lstrip().startswith("```"):
            in_code = not in_code
            out.append(line)
            continue
        m = HEADING_RE.match(line) if not in_code else None
        if m:
            title = _EMPHASIS_RE.sub("", m.group(2)).strip()
            out.append(f"{m.group(1)} {title}")
        else:
            out.append(line)
    return "\n".join(out)


def clean_markdown(
    text: str,
    drop_index_refs: bool = True,
    remove_page_numbers: bool = True,
    remove_repeated: bool = True,
    normalize_headings: bool = True,
) -> tuple[str, CleanStats]:
    stats = CleanStats()
    if normalize_headings:
        text = _normalize_headings(text)
    if drop_index_refs:
        text = _drop_noise_sections(text, stats)
    if remove_page_numbers:
        text = _remove_page_numbers(text, stats)
    if remove_repeated:
        text = _remove_repeated_headers(text, stats)
    text = MULTI_BLANK_RE.sub("\n\n", text).strip() + "\n"
    return text, stats
