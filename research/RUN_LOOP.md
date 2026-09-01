# 自走ループの回し方

## いちばん短い説明

```bash
# 手元（実データのあるマシン）で
cd ~/bakanposs
git pull origin claude/new-branch-creation-heq0i1
: > loop.log                                      # 先に作る（tail の競争負けを防ぐ）
nohup research/run_loop.sh -n 30 >> loop.log 2>&1 &
tail -f loop.log
```

**人手が要る作業だけを残して止まります。** 止まったら `loop.log` の末尾か
`research/.loop_stop` に「あなたに何を頼みたいか」がバッチにまとまって出ています。
対応したら `rm research/.loop_stop` して同じコマンドで再開。

---

## なぜプロセスを 1 周ごとに分けるのか

**これが文脈枯渇への答えです。**

1 つの長いセッションで回すと、いずれコンテキスト上限に当たって作業が途切れます。
このドライバは **1 周ごとに `claude` を別プロセスで起動**します。各周はまっさらな文脈で
`research/SESSION_STATE.md` から状態を復元して始まるので、**そもそも上限に届きません**。

そのため、**状態がディスクに正しく書かれていることが全て**です。
だから `ITERATION_PROMPT.md` は「`SESSION_STATE.md` を更新せずに周を終えるな」を最上位の規則にしています。

```
周1: claude 起動 → SESSION_STATE 読む → 1旗 進める → 記録 → commit+push → 終了
周2: claude 起動 → SESSION_STATE 読む → 1旗 進める → 記録 → commit+push → 終了
...
周N: claude 起動 → 全部 BLOCKED と判断 → .loop_stop を書く → 終了 → ドライバが止まる
```

## ★ 手元で回すとゲートが減る

`HUMAN_GATES.md` の **D 分類（GATE-01〜04, 06）は「コンテナから `/mnt/hdd` に届かない」ことが理由**です。
**実データのあるマシンで回せば、この理由は消えます。** ドライバは起動時に `/mnt/hdd` の有無を見て
どちらの環境かを表示し、`ITERATION_PROMPT.md` も Claude 自身に確認させています。

| 回す場所 | 進むもの | 残るゲート |
|---|---|---|
| **手元（実データあり）** | 手 C・手 B の実検定、`RUN_HEAVY.md` ①〜④、PhenoCam 再解析まで**自走** | A（取得）・P（外界）・J（判断）・S（仕様）だけ |
| コンテナ | 事前登録・合成での帰無較正・道具作り・文章・記録の整理 | D も残る |

## 止まる条件（4 つ）

| 条件 | ドライバの挙動 |
|---|---|
| Claude が `research/.loop_stop` を書いた | **正常終了。** 内容（あなたへの依頼）を表示 |
| 2 周続けてコミットが生まれない | 停止。プロンプトか状態ファイルを見直す合図 |
| 未コミットの変更が残った（打ち切り・異常終了） | **停止して、勝手に捨てない。** `git diff` で確認 |
| `-n` の周回数に到達 | 停止 |

`touch research/.loop_stop` すれば、次の周の前に自分で止められます。

## 事前チェック（無人放置で事故らないための門）

起動時に一度だけ確認し、1 つでも欠ければ**回さずに落ちます**。

- `claude` / `git` が PATH にあるか
- `SESSION_STATE.md` / `HUMAN_GATES.md` / `LOOP_PROTOCOL.md` が揃っているか
- **`main` にいないか**（`LOOP_PROTOCOL.md` の禁止事項）
- **想定ブランチにいるか**（`-b` / `--any-branch` で変えられる）。
  Claude は**今いるブランチ**へ push するので、ここを緩めると記録が散らばる
- **作業ツリーがクリーンか**（汚れていれば前の周の中断の痕跡なので、人が見るまで進まない）
- **`.loop_stop` が残っていないか**（人間待ちのまま回すのを防ぐ）

## オプション

```
-n, --max-flags N    最大何周 (default 20)
-t, --timeout SEC    1 周の上限秒 (default 3600)。超えたら打ち切り
-m, --model NAME     opus / sonnet / fable (default opus)
-b, --branch NAME    回してよいブランチ (default claude/new-branch-creation-heq0i1)
    --any-branch     今いるブランチが何であれ回す
    --yolo           --dangerously-skip-permissions。無人放置向け
    --no-push        push しない
    --dry-run        プロンプトと起動コマンドを見るだけ
```

**`--yolo` について正直に**: 無人で放置するには権限の確認を素通しする必要がありますが、
これは**実データのあるあなたのマシンで、Claude が任意のコマンドを確認なしに実行する**ということです。
`ITERATION_PROMPT.md` の禁止事項（過去の記録を書き換えない・main に push しない・PR を作らない）は
守らせていますが、**まず `--yolo` なしで数周見てから**判断してください。
既定は `--permission-mode acceptEdits`（ファイル編集は通す、コマンドは確認）です。

## うまく動かないとき

### `tail: cannot open 'loop.log'` と出る
`[1] <PID>` が出ていれば**ジョブは起動している**。親シェルが PID を表示した時点では、
子プロセスがまだ `loop.log` を開いていないことがある（競争負け）。
上の起動コマンドのように **`: > loop.log` で先にファイルを作ってから** `nohup` すれば起きない。
既に起きてしまったら `ls -la loop.log && tail -f loop.log` でよい。

### すぐ終わってしまう
事前チェックの門に引っかかっている。**理由は `loop.log` の先頭 1〜2 行に出ている。**

| `loop.log` の中身 | 意味 | 直し方 |
|---|---|---|
| `claude が PATH にない` | venv や nvm の関係で見えていない | `which claude` で確認。無ければ PATH を通す |
| `main では回さない` | main にいる | `git checkout claude/new-branch-creation-heq0i1` |
| `想定ブランチと違う` | 別ブランチにいる。**Claude は今いるブランチへ push する**ので、意図しない所に研究の記録が積まれる | `git checkout` で移るか、`-b <今のブランチ>` / `--any-branch` で明示的に許可 |
| `未コミットの変更がある` | 前の周の中断の痕跡かもしれない | `git status` / `git diff` で確認し、完結させるか捨てる |
| `前回のループが人間待ちで止まっている` | `.loop_stop` が残っている | 内容に対応して `rm research/.loop_stop` |
| `<file> が無い` | 状態ファイルが揃っていない | `git pull` |
| `Permission denied` | 実行ビットが落ちている | `chmod +x research/run_loop.sh` |

### 動いているか確かめる
```bash
jobs                      # このシェルから起動した場合
pgrep -af run_loop.sh     # シェルを閉じた後でも
tail -f loop.log
```

### 止める
```bash
touch research/.loop_stop   # 次の周の前に、きれいに止まる（推奨）
pkill -f run_loop.sh        # 今すぐ止める（周の途中なら未コミットの変更が残りうる）
```

## ログ

- `loop.log` — ドライバの進行（周ごとの一行と最後のまとめ）
- `research/loop_runs/*.jsonl` — 各周の生ログ（gitignore 済）
- `research/loop_tail.py` — 生ログを読める形に潰す。ドライバが自動で呼ぶ

## ファイル

| ファイル | 役割 |
|---|---|
| `research/run_loop.sh` | ドライバ本体 |
| `research/ITERATION_PROMPT.md` | 1 周ぶんの指示（これを claude に渡す） |
| `research/LOOP_PROTOCOL.md` | ループ規約・禁止事項・上限プロトコル |
| `research/HUMAN_GATES.md` | 人手が要る箇所の台帳 |
| `research/COLD_START.md` | 手で再開するとき用（対話セッションに貼る） |
| `research/GATE_AUDIT_PROMPT.md` | 台帳を洗い直すとき |
| `research/SESSION_STATE.md` | **現在地。ループの生命線** |
