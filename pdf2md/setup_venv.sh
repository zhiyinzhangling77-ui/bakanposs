#!/usr/bin/env bash
# 専用 venv を作り、marker-pdf(主) + MinerU(フォールバック)などを入れる。
# GPU があれば自動で使い、無ければ CPU で動く。
#
# 使い方:
#   bash setup_venv.sh            # CPU / 既定
#   PDF2MD_CUDA=1 bash setup_venv.sh   # CUDA 対応 torch を先に入れてから依存を入れる
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
VENV="$HERE/.venv"

echo "==> venv を作成: $VENV"
python3 -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"

python -m pip install --upgrade pip wheel

if [ "${PDF2MD_CUDA:-0}" = "1" ]; then
  echo "==> CUDA 対応 torch を先にインストール(CUDA 12.1 例)"
  pip install torch --index-url https://download.pytorch.org/whl/cu121
fi

echo "==> 依存パッケージをインストール"
pip install -r "$HERE/requirements.txt"

echo
echo "==> インストール確認"
python - <<'PY'
import shutil
print(f"  marker_single: {'OK' if shutil.which('marker_single') else '見つかりません'}")
print("  mineru: 既定では入れません(marker と Pillow が衝突するため。別venv=requirements-mineru.txt)")
try:
    import torch
    dev = "cuda" if torch.cuda.is_available() else (
        "mps" if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available() else "cpu")
    print(f"  torch device: {dev}")
except Exception as e:
    print(f"  torch: 未確認 ({e})")
PY

cat <<'EOF'

セットアップ完了。使い方:
  source .venv/bin/activate
  export ANTHROPIC_API_KEY=sk-ant-...          # 要約に必要(未設定なら要約はスキップ)

  # 1) まず品質サンプル(先頭5ページ)。--input は単一 .pdf でも可
  python -m pdf2md sample --input <PDFまたはフォルダ> --pages 5
  #   低VRAMのGPUで CUDA OOM になる場合は --device cpu を付ける(自動再試行もあり)

  # 2) 目次を確認
  python -m pdf2md toc --input <PDFフォルダ>

  # 3) 本処理(サンプル確認後)
  python -m pdf2md run --input <PDFフォルダ> --output <Obsidian保存先> --reviewed-sample
EOF
