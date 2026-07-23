"""sample / toc / run の処理本体。"""

from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass, field
from pathlib import Path

from . import clean, convert, split, summarize, toc


def find_pdfs(input_path: Path) -> list[Path]:
    """フォルダなら中の *.pdf を、単一の .pdf ファイルならそれ1つを返す。"""
    if input_path.is_file() and input_path.suffix.lower() == ".pdf":
        return [input_path]
    return sorted(p for p in input_path.glob("*.pdf") if p.is_file())


# ---------------------------------------------------------------- sample

def run_sample(
    input_dir: Path,
    pdf: str | None,
    pages: int,
    prefer: str,
    save_to: Path | None,
    device: str = "auto",
) -> str:
    pdfs = find_pdfs(input_dir)
    if not pdfs:
        raise FileNotFoundError(f"{input_dir} に PDF がありません。")
    # --input がフォルダで --pdf 指定があればその1冊、無ければ先頭。単一ファイル入力ならそれ。
    if pdf and input_dir.is_dir():
        target = input_dir / pdf
    else:
        target = pdfs[0]
    if not target.exists():
        raise FileNotFoundError(target)

    eff_device = convert.detect_device() if device == "auto" else device
    page_range = f"0-{max(0, pages - 1)}"
    res = convert.convert(target, page_range=page_range, device=device, prefer=prefer)
    cleaned, _ = clean.clean_markdown(res.markdown)

    header = (
        f"# 品質サンプル: {target.name}\n\n"
        f"- バックエンド: {res.backend}\n"
        f"- デバイス: {eff_device}\n"
        f"- ページ: 先頭 {pages} ページ\n\n---\n\n"
    )
    sample = header + cleaned
    if save_to:
        save_to.parent.mkdir(parents=True, exist_ok=True)
        save_to.write_text(sample, encoding="utf-8")
    return sample


# ---------------------------------------------------------------- toc

def run_toc(input_dir: Path) -> str:
    pdfs = find_pdfs(input_dir)
    if not pdfs:
        raise FileNotFoundError(f"{input_dir} に PDF がありません。")
    blocks: list[str] = []
    for p in pdfs:
        entries = toc.extract_pdf_outline(p)
        chapters = toc.top_level_chapters(entries)
        n = len(chapters)
        blocks.append(
            f"■ {p.name}  (章候補: {n} 章)\n" + toc.format_toc(entries)
        )
    return "\n\n".join(blocks)


# ---------------------------------------------------------------- run

@dataclass
class FileOutcome:
    pdf: str
    ok: bool
    chapters_written: int = 0
    backend: str = ""
    summaries_ok: int = 0
    summaries_failed: int = 0
    note: str = ""


@dataclass
class RunReport:
    outcomes: list[FileOutcome] = field(default_factory=list)

    def as_markdown(self) -> str:
        ts = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        lines = [f"# 変換ログ ({ts})\n", "| PDF | 結果 | 章数 | backend | 要約OK/NG | 備考 |",
                 "|---|---|---|---|---|---|"]
        for o in self.outcomes:
            status = "✅成功" if o.ok else "❌失敗"
            lines.append(
                f"| {o.pdf} | {status} | {o.chapters_written} | {o.backend} "
                f"| {o.summaries_ok}/{o.summaries_failed} | {o.note} |"
            )
        ok = sum(1 for o in self.outcomes if o.ok)
        ng = len(self.outcomes) - ok
        lines.append(f"\n**合計: 成功 {ok} / 失敗 {ng}**\n")
        return "\n".join(lines)


def _targets_would_overwrite(book_dir: Path) -> list[Path]:
    if not book_dir.exists():
        return []
    return [p for p in book_dir.glob("*.md")]


def process_one_pdf(
    pdf: Path,
    output_root: Path,
    chapters_spec: str | None,
    prefer: str,
    model: str,
    make_summary: bool,
    device: str,
    page_range: str | None = None,
    split_mode: str = "h1",
    chapter_pattern: str | None = None,
) -> FileOutcome:
    outcome = FileOutcome(pdf=pdf.name, ok=False)
    try:
        res = convert.convert(
            pdf, page_range=page_range, device=device, prefer=prefer
        )
        outcome.backend = res.backend
        cleaned, _ = clean.clean_markdown(res.markdown)
        frontmatter, chapters = split.split_chapters(
            cleaned, mode=split_mode, pattern=chapter_pattern, fallback_title=pdf.stem
        )

        if not chapters:
            # 見出しが取れない場合は 1 ファイルとして扱う
            chapters = [split.Chapter(index=1, title=pdf.stem, body=cleaned)]

        selected = split.parse_selection(chapters_spec, len(chapters))
        book_dir = output_root / split._sanitize(pdf.stem, maxlen=80)
        book_dir.mkdir(parents=True, exist_ok=True)

        for ch in chapters:
            if ch.index not in selected:
                continue
            if make_summary:
                summ = summarize.summarize_chapter(ch.body, model=model)
            else:
                summ = summarize.SummaryResult(
                    text=summarize._fallback_header("要約なしで実行"), ok=False
                )
            outcome.summaries_ok += int(summ.ok)
            outcome.summaries_failed += int(not summ.ok)
            header = summarize.build_chapter_header(pdf.name, ch.title, ch.index, summ)
            (book_dir / split.chapter_filename(ch)).write_text(
                header + ch.body, encoding="utf-8"
            )
            outcome.chapters_written += 1

        outcome.ok = outcome.chapters_written > 0
        if not outcome.ok:
            outcome.note = "書き出した章がありません"
    except Exception as e:  # noqa: BLE001 — 失敗ファイルはスキップして続行
        outcome.note = str(e)[:200]
    return outcome


def run_pipeline(
    input_dir: Path,
    output_dir: Path,
    chapters_spec: str | None,
    only: str | None,
    prefer: str,
    model: str,
    make_summary: bool,
    device: str = "auto",
    page_range: str | None = None,
    split_mode: str = "h1",
    chapter_pattern: str | None = None,
) -> RunReport:
    pdfs = find_pdfs(input_dir)
    if only:
        pdfs = [p for p in pdfs if p.name == only]
    if not pdfs:
        raise FileNotFoundError("対象 PDF がありません。")

    # device は各 convert 呼び出しへ渡す("auto" は convert 内で判定・OOM時CPU再試行)
    output_dir.mkdir(parents=True, exist_ok=True)
    report = RunReport()
    for p in pdfs:
        report.outcomes.append(
            process_one_pdf(
                p, output_dir, chapters_spec, prefer, model,
                make_summary, device, page_range, split_mode, chapter_pattern,
            )
        )  # device="auto" 可。convert 内で判定し、CUDA OOM 時は自動で CPU 再試行
    log_path = output_dir / "_conversion_log.md"
    log_path.write_text(report.as_markdown(), encoding="utf-8")
    return report
