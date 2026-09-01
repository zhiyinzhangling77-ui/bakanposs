#!/usr/bin/env bash
# run_loop.sh — 研究の自走ループ・ドライバ
#
# 1 周ごとに claude を別プロセスで起動する。各周はまっさらな文脈で
# research/SESSION_STATE.md から起動するため、**文脈の枯渇が構造的に起きない**。
# 状態はすべてディスク（SESSION_STATE / FLAGS_LOG / HUMAN_GATES）にある。
#
# ★ 実データ（/mnt/hdd）のあるローカルで回すと、HUMAN_GATES の D 分類
#   （GATE-01〜04, 06）はゲートでなくなり、Claude が自分で実行できる。
#   コンテナで回す場合は合成検証・道具作り・文章だけが進む。
#
#   : > loop.log
#   nohup research/run_loop.sh -n 30 >> loop.log 2>&1 &
#   tail -f loop.log
#
# 止め方: Ctrl-C か、`touch research/.loop_stop`（次の周の前に止まる）

set -uo pipefail

BRANCH="claude/new-branch-creation-heq0i1"
ANY_BRANCH=0
MAX_FLAGS=20
TIMEOUT_SEC=3600
MODEL="opus"
PERM_MODE="acceptEdits"
DO_PUSH=1
DRY_RUN=0
STALL_LIMIT=2

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STOP_FILE="$REPO/research/.loop_stop"
PROMPT_FILE="$REPO/research/ITERATION_PROMPT.md"
LOGDIR="$REPO/research/loop_runs"

usage() {
  sed -n '2,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
  cat <<'USAGE'

Options:
  -n, --max-flags N    最大何周回すか (default 20)
  -t, --timeout SEC    1 周あたりの上限秒 (default 3600)
  -m, --model NAME     opus / sonnet / fable など (default opus)
  -b, --branch NAME    回してよいブランチ (default claude/new-branch-creation-heq0i1)
      --any-branch     今いるブランチが何であれ回す（push もそのブランチへ）
      --yolo           全権限を素通し (--dangerously-skip-permissions)。
                       無人放置向けだが、実データのあるマシンでは中身を理解してから使うこと
      --no-push        push しない（ローカルで試すとき）
      --dry-run        プロンプトと起動コマンドを表示して終了
  -h, --help           これ
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--max-flags) MAX_FLAGS="$2"; shift 2 ;;
    -t|--timeout)   TIMEOUT_SEC="$2"; shift 2 ;;
    -m|--model)     MODEL="$2"; shift 2 ;;
    -b|--branch)    BRANCH="$2"; shift 2 ;;
    --any-branch)   ANY_BRANCH=1; shift ;;
    --yolo)         PERM_MODE="__yolo__"; shift ;;
    --no-push)      DO_PUSH=0; shift ;;
    --dry-run)      DRY_RUN=1; shift ;;
    -h|--help)      usage; exit 0 ;;
    *) echo "不明な引数: $1" >&2; usage; exit 2 ;;
  esac
done

log() { printf '%s | %s\n' "$(date '+%F %T')" "$*"; }
die() { log "✗ $*"; exit 1; }

# ---------- 事前チェック（一度だけ。ここで弾けば無人放置しても事故らない） ----------

cd "$REPO" || die "リポジトリに移動できない: $REPO"

command -v claude >/dev/null || die "claude が PATH にない"
command -v git    >/dev/null || die "git が PATH にない"
[[ -f "$PROMPT_FILE" ]] || die "$PROMPT_FILE が無い"

for f in research/SESSION_STATE.md research/HUMAN_GATES.md research/LOOP_PROTOCOL.md; do
  [[ -f "$REPO/$f" ]] || die "$f が無い。COLD_START の前提が崩れている"
done

CLAUDE_ARGS=( -p --model "$MODEL" --output-format stream-json --verbose )
if [[ "$PERM_MODE" == "__yolo__" ]]; then
  CLAUDE_ARGS+=( --dangerously-skip-permissions )
else
  CLAUDE_ARGS+=( --permission-mode "$PERM_MODE" )
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
  log "起動コマンド:"
  printf '  claude'; printf ' %q' "${CLAUDE_ARGS[@]}"; printf ' "$(cat %s)"\n' "$PROMPT_FILE"
  echo; echo "───── プロンプト ─────"; cat "$PROMPT_FILE"
  exit 0
fi

CUR_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
[[ "$CUR_BRANCH" == "main" || "$CUR_BRANCH" == "master" ]] && \
  die "main では回さない（LOOP_PROTOCOL の禁止事項）。git checkout $BRANCH"
if [[ "$CUR_BRANCH" != "$BRANCH" && $ANY_BRANCH -eq 0 ]]; then
  die "想定ブランチと違う。
     いる:   $CUR_BRANCH
     想定:   $BRANCH
     Claude は今いるブランチへ push する。意図しないブランチに研究の記録が積まれるので、
     ここでは進めない。どれかを選ぶこと:
       git checkout $BRANCH          # 想定ブランチへ移る
       research/run_loop.sh -b $CUR_BRANCH   # このブランチで回してよいと明示する
       research/run_loop.sh --any-branch     # ブランチを問わない"
fi
log "ブランチ: $CUR_BRANCH"

if [[ -f "$STOP_FILE" ]]; then
  echo "──────────────────────────────────────────"
  cat "$STOP_FILE"
  echo "──────────────────────────────────────────"
  die "前回のループが人間待ちで止まっている（上記）。
     対応してから rm $STOP_FILE で再開すること。"
fi

# 見たいのは「中断された編集の痕跡」＝**追跡下のファイルの変更**。
# 解析の生成物（outputs_*/・*.log・図）は untracked で常にそこにあるので、門にしない。
if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
  echo "── 追跡下の未コミット変更 ──"
  git status --short --untracked-files=no
  die "追跡下に未コミットの変更がある。前の周が中断された痕跡かもしれない。
     git diff で中身を確認し、完結させるか捨ててから回すこと。"
fi

untracked_n="$(git ls-files --others --exclude-standard | wc -l | tr -d ' ')"
if [[ "$untracked_n" -gt 0 ]]; then
  log "未追跡ファイル $untracked_n 件（解析の生成物とみなして無視する）"
fi

# timeout（macOS は coreutils の gtimeout）
TIMEOUT_BIN=""
if command -v timeout  >/dev/null; then TIMEOUT_BIN="timeout"
elif command -v gtimeout >/dev/null; then TIMEOUT_BIN="gtimeout"
else log "⚠ timeout が無いので 1 周あたりの上限をかけられない"; fi

# 実データの有無 = このマシンで D 分類のゲートが外れるか
if [[ -d /mnt/hdd ]]; then
  log "環境: LOCAL（/mnt/hdd あり）→ D 分類のゲートは外れる。実データ解析まで自走できる"
else
  log "環境: CONTAINER（/mnt/hdd なし）→ 合成検証・道具作り・文章のみ進む"
fi

mkdir -p "$LOGDIR"

# ---------- ループ ----------

PROMPT="$(cat "$PROMPT_FILE")"
START_HEAD="$(git rev-parse HEAD)"
stall=0
done_flags=0
stop_reason="最大周回数に到達"

log "開始: 最大 $MAX_FLAGS 周 / 1 周上限 ${TIMEOUT_SEC}s / model=$MODEL / branch=$CUR_BRANCH"
log "起点: $START_HEAD"

for (( i=1; i<=MAX_FLAGS; i++ )); do
  echo
  log "━━━━━━━━━━ 周 $i / $MAX_FLAGS ━━━━━━━━━━"

  before="$(git rev-parse HEAD)"
  raw="$LOGDIR/$(date '+%Y%m%d-%H%M%S')-iter$(printf '%03d' "$i").jsonl"

  if [[ -n "$TIMEOUT_BIN" ]]; then
    "$TIMEOUT_BIN" "$TIMEOUT_SEC" claude "${CLAUDE_ARGS[@]}" "$PROMPT" > "$raw" 2>&1
  else
    claude "${CLAUDE_ARGS[@]}" "$PROMPT" > "$raw" 2>&1
  fi
  rc=$?

  # 生 JSONL から人が読める形を取り出す（python3 が無ければ素通し）
  if command -v python3 >/dev/null; then
    python3 "$REPO/research/loop_tail.py" "$raw" || tail -20 "$raw"
  else
    tail -20 "$raw"
  fi

  if [[ $rc -eq 124 ]]; then
    log "⚠ 周 $i は ${TIMEOUT_SEC}s で打ち切られた"
  elif [[ $rc -ne 0 ]]; then
    log "⚠ claude が rc=$rc で終了（ログ: $raw）"
  fi

  after="$(git rev-parse HEAD)"

  # 打ち切り・異常終了で中途半端な変更が残っていたら、そこで止める（勝手に捨てない）
  if [[ -n "$(git status --porcelain)" ]]; then
    stop_reason="未コミットの変更が残った（周 $i・rc=$rc）。git diff を確認すること"
    log "✗ $stop_reason"
    break
  fi

  if [[ -f "$STOP_FILE" ]]; then
    stop_reason="人間待ち（.loop_stop）"
    log "■ Claude が人間待ちで停止を要求した"
    break
  fi

  if [[ "$after" == "$before" ]]; then
    stall=$(( stall + 1 ))
    log "⚠ 周 $i はコミットを生まなかった（$stall/$STALL_LIMIT）"
    if [[ $stall -ge $STALL_LIMIT ]]; then
      stop_reason="$STALL_LIMIT 周続けて進捗が無い。プロンプトか状態ファイルを見直すこと"
      log "✗ $stop_reason"
      break
    fi
  else
    stall=0
    done_flags=$(( done_flags + 1 ))
    log "✓ 周 $i 完了: $(git log --oneline -1)"
    if [[ $DO_PUSH -eq 1 ]]; then
      for attempt in 1 2 3 4; do
        if git push -u origin "$CUR_BRANCH" >/dev/null 2>&1; then
          log "  push 済"; break
        fi
        log "  push 失敗（$attempt 回目）。$(( 2 ** attempt ))s 待って再試行"
        sleep $(( 2 ** attempt ))
      done
    fi
  fi
done

# ---------- まとめ ----------

echo
log "━━━━━━━━━━ 終了 ━━━━━━━━━━"
log "停止理由: $stop_reason"
log "進んだ周: $done_flags"
if [[ "$(git rev-parse HEAD)" != "$START_HEAD" ]]; then
  echo
  echo "このループで積んだコミット:"
  git log --oneline "$START_HEAD..HEAD"
fi
if [[ -f "$STOP_FILE" ]]; then
  echo
  echo "══════════ あなたの手が必要です ══════════"
  cat "$STOP_FILE"
  echo "════════════════════════════════════════"
  echo "対応したら: rm $STOP_FILE && research/run_loop.sh"
fi
echo
echo "開いている GATE:"
grep -n '^### GATE-' research/HUMAN_GATES.md | sed 's/^/  /' || true
