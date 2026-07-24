#!/usr/bin/env bash
# 複数の本を順番に(--by-toc で)Markdown 化するキュー。
#
# 使い方:
#   1) 下の「本の一覧」を編集(1行1冊: convert_book "PDFパス" "出力フォルダ")
#   2) nohup で起動:  nohup bash ~/convert_books.sh >/dev/null 2>&1 &
#   すでに走っている変換ジョブ(run_dl_convert.sh)が終わるまで待ってから開始する。
#
# 前提: 各PDFに「埋め込み目次」があること(--by-toc はそれで章分割する)。
#       目次が無い本は失敗としてログに残るので、その本だけ後で --page-range で個別変換。
set -uo pipefail

REPO="$HOME/bakanposs/pdf2md"      # pdf2md パッケージのある場所(環境に合わせて)
PY="$REPO/.venv/bin/python"        # marker 入りの venv
LOG="$HOME/md_queue.log"

convert_book() {
  local pdf="$1" out="$2"
  echo ""                                              >> "$LOG"
  echo ">>> $(date '+%F %T')  START: $pdf"             >> "$LOG"
  if [ ! -f "$pdf" ]; then
    echo "    SKIP: PDF が見つかりません: $pdf"          >> "$LOG"; return
  fi
  ( cd "$REPO" && CUDA_VISIBLE_DEVICES="" "$PY" -m pdf2md run \
      --input "$pdf" --output "$out" \
      --device cpu --reviewed-sample --no-summary --by-toc --yes \
      >> "$LOG" 2>&1 )
  echo ">>> $(date '+%F %T')  END:   $pdf"             >> "$LOG"
}

echo "=== md queue 開始: $(date) ===" >> "$LOG"

# 先に走っている変換ジョブが終わるまで待機(無ければ即開始)
while pgrep -f 'run_dl_convert.sh' >/dev/null; do
  sleep 60
done

# ==================== 本の一覧(ここを編集) ====================
# convert_book "PDFのフルパス" "Obsidianの出力フォルダ"
convert_book "/home/shion-nagamine/Desktop/BOOK1.pdf" "/home/shion-nagamine/ObsidianVault/BOOK1"
convert_book "/home/shion-nagamine/Desktop/BOOK2.pdf" "/home/shion-nagamine/ObsidianVault/BOOK2"
# =============================================================

echo "=== md queue 完了: $(date) ===" >> "$LOG"
touch "$HOME/md_queue.DONE"
command -v notify-send >/dev/null && \
  notify-send "本のMD化キュー完了" "$(grep -c 'END:' "$LOG") 冊処理" 2>/dev/null || true
