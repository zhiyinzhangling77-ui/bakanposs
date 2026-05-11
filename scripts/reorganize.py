#!/usr/bin/env python3
"""
bakanposs ディレクトリ整理スクリプト

実行前に git status でコミット漏れがないか確認してください。
dry_run=True で変更内容を確認してから、False に切り替えて実行。
"""

import shutil
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
DRY_RUN = True  # False に変えると実際に移動する


def mv(src, dst, dry=DRY_RUN):
    src, dst = Path(src), Path(dst)
    if not src.exists():
        print(f"  SKIP  (not found) {src.relative_to(ROOT)}")
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        print(f"  SKIP  (dest exists) {src.relative_to(ROOT)} → {dst.relative_to(ROOT)}")
        return
    action = "DRY" if dry else "MOVE"
    print(f"  {action}  {src.relative_to(ROOT)}  →  {dst.relative_to(ROOT)}")
    if not dry:
        shutil.move(str(src), str(dst))


def main():
    print(f"ROOT: {ROOT}")
    print(f"DRY_RUN: {DRY_RUN}")
    print("=" * 60)

    # ----------------------------------------------------------------
    # 1. analysis_A_v10〜v14.py をルートから analysis_A/ へ
    # ----------------------------------------------------------------
    print("\n[1] analysis_A_v10〜v14.py → analysis_A/")
    for v in range(10, 15):
        mv(ROOT / f"analysis_A_v{v}.py",
           ROOT / "analysis_A" / f"analysis_A_v{v}.py")

    # ----------------------------------------------------------------
    # 2. output_analysis_C_v1/ を analysis_A/ → analysis_C/ へ
    #    (解析Cの出力が誤って analysis_A/ 以下に入っている)
    # ----------------------------------------------------------------
    print("\n[2] analysis_A/output_analysis_C_v1/ → analysis_C/output_analysis_C_v1/")
    mv(ROOT / "analysis_A" / "output_analysis_C_v1",
       ROOT / "analysis_C" / "output_analysis_C_v1")

    # ----------------------------------------------------------------
    # 3. ルートの run_C_v9.log → analysis_C/
    # ----------------------------------------------------------------
    print("\n[3] run_C_v9.log → analysis_C/")
    mv(ROOT / "run_C_v9.log",
       ROOT / "analysis_C" / "run_C_v9.log")

    # ----------------------------------------------------------------
    # 残す場所の確認 (移動しない)
    # ----------------------------------------------------------------
    print("\n[確認] ルートに残すファイル:")
    keep = [
        "analysis_A_v9.py",    # analysis_C_v1.py が from analysis_A_v9 import している
        "analysis_C_v1.py",    # メインスクリプト (data_loaders と同じ場所が必要)
        "data_loaders.py",     # 共通 loader (A/B/C から import)
        "paper_methods_results.md",
        "paper_outline.md",
        "requirements.txt",
    ]
    for f in keep:
        p = ROOT / f
        status = "OK" if p.exists() else "MISSING"
        print(f"  [{status}] {f}")

    print("\n完了。DRY_RUN=False で再実行すると実際に移動します。")


if __name__ == "__main__":
    main()
