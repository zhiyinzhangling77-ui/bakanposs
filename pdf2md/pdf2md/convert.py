"""PDF → Markdown 変換バックエンド。

marker-pdf を主バックエンド、MinerU をフォールバックとして使う。
どちらも CLI 経由で呼ぶ(Python API はバージョン間で変わりやすいため、CLI の方が安定)。
GPU があれば使い、無ければ CPU で動く。
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ConversionResult:
    markdown: str
    backend: str          # "marker" / "mineru"
    pages: str | None     # 変換したページ範囲(例: "0-4")、全ページなら None


class ConversionError(RuntimeError):
    pass


def detect_device() -> str:
    """利用可能な計算デバイスを返す(cuda / mps / cpu)。"""
    try:
        import torch  # marker/mineru が入っていれば torch も入る

        if torch.cuda.is_available():
            return "cuda"
        mps = getattr(torch.backends, "mps", None)
        if mps is not None and mps.is_available():
            return "mps"
    except Exception:
        pass
    return "cpu"


def backend_available(name: str) -> bool:
    exe = {"marker": "marker_single", "mineru": "mineru"}.get(name)
    return exe is not None and shutil.which(exe) is not None


def _find_markdown(root: Path) -> str | None:
    """出力ディレクトリ以下から最初に見つかった .md を読む。"""
    mds = sorted(root.rglob("*.md"))
    if not mds:
        return None
    # 一番大きい .md を本文とみなす(付随する小さいメタ .md を避ける)
    md = max(mds, key=lambda p: p.stat().st_size)
    return md.read_text(encoding="utf-8", errors="replace")


def _run(cmd: list[str], env: dict, timeout: int) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd, env=env, capture_output=True, text=True, timeout=timeout
    )


def convert_with_marker(
    pdf: Path, page_range: str | None, device: str, timeout: int = 3600
) -> str:
    """marker_single CLI で変換して Markdown 文字列を返す。"""
    env = dict(os.environ)
    env["TORCH_DEVICE"] = device
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        cmd = [
            "marker_single",
            str(pdf),
            "--output_dir",
            str(out),
            "--output_format",
            "markdown",
        ]
        if page_range:
            # marker のページ指定は 0 始まり(例: "0-4" で先頭5ページ)
            cmd += ["--page_range", page_range]
        proc = _run(cmd, env, timeout)
        md = _find_markdown(out)
        if md is None:
            raise ConversionError(
                f"marker が Markdown を生成しませんでした。\n"
                f"stdout: {proc.stdout[-500:]}\nstderr: {proc.stderr[-500:]}"
            )
        return md


def convert_with_mineru(
    pdf: Path, page_range: str | None, device: str, timeout: int = 3600
) -> str:
    """mineru CLI で変換して Markdown 文字列を返す(日本語・中国語に強い)。"""
    env = dict(os.environ)
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        cmd = ["mineru", "-p", str(pdf), "-o", str(out)]
        # デバイス指定(mineru 2.x: -d cuda|cpu)。未対応バージョンでも -d は概ね受理される。
        cmd += ["-d", "cuda" if device == "cuda" else "cpu"]
        if page_range:
            # mineru は開始/終了ページ(0 始まり)を個別指定する
            start, _, end = page_range.partition("-")
            if start:
                cmd += ["-s", start]
            if end:
                cmd += ["-e", end]
        proc = _run(cmd, env, timeout)
        md = _find_markdown(out)
        if md is None:
            raise ConversionError(
                f"mineru が Markdown を生成しませんでした。\n"
                f"stdout: {proc.stdout[-500:]}\nstderr: {proc.stderr[-500:]}"
            )
        return md


def convert(
    pdf: Path,
    page_range: str | None = None,
    device: str | None = None,
    prefer: str = "marker",
    timeout: int = 3600,
) -> ConversionResult:
    """1冊を変換する。prefer を先に試し、失敗したら他方にフォールバック。"""
    device = device or detect_device()
    order = ["marker", "mineru"] if prefer == "marker" else ["mineru", "marker"]
    order = [b for b in order if backend_available(b)]
    if not order:
        raise ConversionError(
            "marker_single も mineru も見つかりません。setup_venv.sh で環境を作ってください。"
        )

    errors: list[str] = []
    for backend in order:
        try:
            fn = convert_with_marker if backend == "marker" else convert_with_mineru
            md = fn(pdf, page_range, device, timeout)
            if md and md.strip():
                return ConversionResult(markdown=md, backend=backend, pages=page_range)
            errors.append(f"{backend}: 空の出力")
        except subprocess.TimeoutExpired:
            errors.append(f"{backend}: タイムアウト({timeout}s)")
        except Exception as e:  # noqa: BLE001 — どのバックエンド失敗もスキップ対象
            errors.append(f"{backend}: {e}")

    raise ConversionError(" / ".join(errors))
