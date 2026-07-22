"""コマンドラインインターフェース。

  python -m pdf2md sample --input <PDFフォルダ> [--pages 5]
  python -m pdf2md toc    --input <PDFフォルダ>
  python -m pdf2md run    --input <PDFフォルダ> --output <Obsidian保存先> [オプション]

方針: まず sample → 目視で確認 → run(本処理)の順。
      既存ファイルを上書きする前には必ず確認を取る(--yes で自動承認)。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import pipeline
from .split import _sanitize


def _confirm(question: str) -> bool:
    try:
        ans = input(f"{question} [y/N]: ").strip().lower()
    except EOFError:
        return False
    return ans in ("y", "yes")


def cmd_sample(args: argparse.Namespace) -> int:
    save_to = Path(args.save) if args.save else None
    sample = pipeline.run_sample(
        Path(args.input), args.pdf, args.pages, args.prefer, save_to
    )
    print(sample)
    if save_to:
        print(f"\n[サンプルを保存しました: {save_to}]", file=sys.stderr)
    print(
        "\n──────────────────────────────────────────\n"
        "この品質でよければ次に進んでください:\n"
        "  1) 目次確認:  python -m pdf2md toc --input <フォルダ>\n"
        "  2) 本処理:    python -m pdf2md run --input <フォルダ> --output <保存先> --reviewed-sample\n"
        "──────────────────────────────────────────",
        file=sys.stderr,
    )
    return 0


def cmd_toc(args: argparse.Namespace) -> int:
    print(pipeline.run_toc(Path(args.input)))
    print(
        "\n特定の章だけ変換するには run に --chapters '1,3,5-7' を渡してください"
        "(指定なしなら全章)。",
        file=sys.stderr,
    )
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    input_dir = Path(args.input)
    output_dir = Path(args.output)

    # サンプル→確認→本処理 のガード
    if not args.reviewed_sample and not args.yes:
        print(
            "本処理の前に、まず `sample` で品質を確認してください。\n"
            "確認済みなら --reviewed-sample を付けて再実行してください。",
            file=sys.stderr,
        )
        return 2

    # 既存ファイルの上書き確認
    pdfs = pipeline.find_pdfs(input_dir)
    if args.only:
        pdfs = [p for p in pdfs if p.name == args.only]
    existing: list[Path] = []
    for p in pdfs:
        book_dir = output_dir / _sanitize(p.stem, maxlen=80)
        existing += pipeline._targets_would_overwrite(book_dir)
    if existing and not args.yes:
        print(
            f"⚠ 保存先に既存の .md が {len(existing)} 件あります(上書きされる可能性があります):",
            file=sys.stderr,
        )
        for e in existing[:10]:
            print(f"    {e}", file=sys.stderr)
        if len(existing) > 10:
            print(f"    … 他 {len(existing) - 10} 件", file=sys.stderr)
        if not _confirm("続行して上書きしますか?"):
            print("中止しました。", file=sys.stderr)
            return 1

    report = pipeline.run_pipeline(
        input_dir=input_dir,
        output_dir=output_dir,
        chapters_spec=args.chapters,
        only=args.only,
        prefer=args.prefer,
        model=args.model,
        make_summary=not args.no_summary,
    )
    print(report.as_markdown())
    print(f"\n[ログ: {output_dir / '_conversion_log.md'}]", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pdf2md",
        description="PDF(教科書・論文) → AIで勉強しやすい Markdown 変換パイプライン",
    )
    sub = p.add_subparsers(dest="command", required=True)

    ps = sub.add_parser("sample", help="1冊の先頭数ページだけ変換して品質を確認")
    ps.add_argument("--input", required=True, help="PDF が入ったフォルダ")
    ps.add_argument("--pdf", help="対象PDFのファイル名(省略時は最初の1冊)")
    ps.add_argument("--pages", type=int, default=5, help="変換するページ数(既定 5)")
    ps.add_argument("--prefer", choices=["marker", "mineru"], default="marker")
    ps.add_argument("--save", help="サンプルの保存先 .md パス(任意)")
    ps.set_defaults(func=cmd_sample)

    pt = sub.add_parser("toc", help="各PDFの目次(見出し構造)を一覧表示")
    pt.add_argument("--input", required=True, help="PDF が入ったフォルダ")
    pt.set_defaults(func=cmd_toc)

    pr = sub.add_parser("run", help="本処理: 変換→ノイズ除去→章分割→要約付与→保存")
    pr.add_argument("--input", required=True, help="PDF が入ったフォルダ")
    pr.add_argument("--output", required=True, help="Obsidian の保存先フォルダ")
    pr.add_argument("--chapters", help="変換する章の指定 例:'1,3,5-7'(省略時は全章)")
    pr.add_argument("--only", help="このファイル名の1冊だけ処理")
    pr.add_argument("--prefer", choices=["marker", "mineru"], default="marker")
    pr.add_argument("--model", default="claude-opus-4-8", help="要約に使う Claude モデル")
    pr.add_argument("--no-summary", action="store_true", help="Claude 要約を付けない")
    pr.add_argument("--reviewed-sample", action="store_true",
                    help="sample を確認済み(このフラグが無いと本処理しない)")
    pr.add_argument("--yes", action="store_true", help="上書き確認等を自動承認")
    pr.set_defaults(func=cmd_run)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except KeyboardInterrupt:
        return 130
    except Exception as e:  # noqa: BLE001
        print(f"エラー: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
