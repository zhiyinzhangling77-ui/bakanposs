"""PDF → Markdown 変換バックエンド。

marker-pdf を主バックエンド、MinerU をフォールバックとして使う。
どちらも CLI 経由で呼ぶ(Python API はバージョン間で変わりやすいため、CLI の方が安定)。
GPU があれば使い、無ければ CPU で動く。
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

# CUDA メモリ不足を示す stderr のキーワード(自動 CPU 再試行の判定に使う)
_OOM_MARKERS = ("out of memory", "outofmemory", "cuda error", "cublas")


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


def resolve_backend_exe(name: str) -> str | None:
    """バックエンドの実行ファイルを探す。

    venv を activate せず `.venv/bin/python -m pdf2md` のように呼ぶと PATH に
    marker_single が無いことがある。その場合に備え、実行中 Python の隣
    (通常は venv の bin/Scripts)も探す。
    """
    exe = {"marker": "marker_single", "mineru": "mineru"}.get(name)
    if exe is None:
        return None
    found = shutil.which(exe)
    if found:
        return found
    bindir = Path(sys.executable).parent  # venv/conda の bin(Windows は Scripts)
    for cand in (bindir / exe, bindir / f"{exe}.exe"):
        if cand.exists():
            return str(cand)
    return None


def backend_available(name: str) -> bool:
    return resolve_backend_exe(name) is not None


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


def _device_env(device: str) -> dict:
    env = dict(os.environ)
    env["TORCH_DEVICE"] = device
    if device == "cpu":
        # 低VRAM GPU での CUDA 使用を確実に無効化(CUDA OOM 回避)
        env["CUDA_VISIBLE_DEVICES"] = ""
    return env


def convert_with_marker(
    pdf: Path, page_range: str | None, device: str, timeout: int = 3600
) -> str:
    """marker_single CLI で変換して Markdown 文字列を返す。"""
    exe = resolve_backend_exe("marker") or "marker_single"
    env = _device_env(device)
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        cmd = [
            exe,
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
    exe = resolve_backend_exe("mineru") or "mineru"
    env = _device_env(device)
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        cmd = [exe, "-p", str(pdf), "-o", str(out)]
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


def _looks_like_oom(msg: str) -> bool:
    m = msg.lower()
    return any(k in m for k in _OOM_MARKERS)


def convert(
    pdf: Path,
    page_range: str | None = None,
    device: str | None = None,
    prefer: str = "marker",
    timeout: int = 3600,
) -> ConversionResult:
    """1冊を変換する。prefer を先に試し、失敗したら他方にフォールバック。

    device: "auto"/None なら自動判定。GPU で CUDA メモリ不足になった場合は
    自動的に CPU で 1 回だけ再試行する(4GB 級の GPU 対策)。
    """
    if device in (None, "auto"):
        device = detect_device()
    order = ["marker", "mineru"] if prefer == "marker" else ["mineru", "marker"]
    order = [b for b in order if backend_available(b)]
    if not order:
        raise ConversionError(
            "marker_single も mineru も見つかりません。setup_venv.sh で環境を作ってください。"
        )

    errors: list[str] = []
    for backend in order:
        fn = convert_with_marker if backend == "marker" else convert_with_mineru
        for attempt_device in ([device, "cpu"] if device != "cpu" else [device]):
            try:
                md = fn(pdf, page_range, attempt_device, timeout)
                if md and md.strip():
                    return ConversionResult(
                        markdown=md, backend=backend, pages=page_range
                    )
                errors.append(f"{backend}[{attempt_device}]: 空の出力")
                break
            except subprocess.TimeoutExpired:
                errors.append(f"{backend}[{attempt_device}]: タイムアウト({timeout}s)")
                break
            except Exception as e:  # noqa: BLE001 — 失敗はスキップ対象
                errors.append(f"{backend}[{attempt_device}]: {e}")
                # CUDA OOM っぽければ CPU で再試行(次ループ)、それ以外は次バックエンドへ
                if attempt_device != "cpu" and _looks_like_oom(str(e)):
                    continue
                break

    raise ConversionError(" / ".join(errors))
