# 無人・自動実行のためのパターン(時間/作業の効率化)

方針: **重い作業を Claude から切り離し、非対話スクリプトにして、スケジューラで無人完走**
させる。こうするとトークン不要・レート制限に無関係・プロンプトで止まらない。

---

## 原則

1. **Claude を処理ループに入れない** — 変換・集計・ダウンロード等の重い処理は
   「普通のコマンド/スクリプト」にする(Claude は“作り方を考える”ときだけ使う)。
2. **完全に非対話にする** — 確認プロンプトが出ない形にする。pdf2md なら
   `--yes`(上書き確認の自動承認)＋`--reviewed-sample`。入力待ちが1つでもあると無人実行は止まる。
3. **ログに落とす** — `> job.log 2>&1`。後で `tail` で結果確認。
4. **スケジューラに載せる** — `at`(一回)/`cron`(定期)/`systemd`。

このパターンなら「全部 Yes で進める」を安全に実現できる。スクリプトは**決まった1つの
ことしかしない**ので、全承認でも想定外の破壊は起きない。

---

## テンプレ: ワンショット・ジョブを作って時刻予約

```bash
# 1) スクリプト化(パス等は直書き。環境変数は予約時に引き継がれないため)
cat > ~/jobs/myjob.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail
cd ~/project || exit 1
# …重い処理を非対話で。例:
# ./.venv/bin/python -m pdf2md run --input "..." --output "..." \
#   --device cpu --reviewed-sample --no-summary --by-toc --yes
SH
chmod +x ~/jobs/myjob.sh

# 2) 一度だけ動作確認(短い範囲で)
bash ~/jobs/myjob.sh

# 3) 時刻予約(例: 深夜2:51)
echo "bash ~/jobs/myjob.sh > ~/jobs/myjob.log 2>&1" | at 02:51
atq        # 予約確認

# デーモンが無ければ: sudo apt install -y at && sudo systemctl enable --now atd
```

`at` が使えない環境の代替(デーモン不要・ログアウトに強い):
```bash
nohup bash -c '
  t=$(date -d "02:51 today" +%s); n=$(date +%s)
  [ "$t" -le "$n" ] && t=$(date -d "02:51 tomorrow" +%s)
  sleep $((t-n)); bash ~/jobs/myjob.sh > ~/jobs/myjob.log 2>&1
' >/dev/null 2>&1 &
```

定期実行なら `crontab -e` に例:
```
51 2 * * *  bash ~/jobs/myjob.sh > ~/jobs/myjob.log 2>&1
```

> ⚠️ **PC がスリープ/電源オフだと動かない**。予約時刻に起きている必要がある
> (自動サスペンド無効化 or 電源つけっぱなし)。

---

## Claude Code 側のプロンプトを減らす(対話作業を速く)

無人ジョブ(上)とは別に、**普段の対話**での確認を減らしたいときは、プロジェクトの
`.claude/settings.json` に**よく使うコマンドの許可リスト**を置く。破壊的操作は
ゲートを残せるので安全寄り。

```jsonc
// .claude/settings.json (例。実際に許可したいものだけ)
{
  "permissions": {
    "allow": [
      "Bash(git status)",
      "Bash(git add:*)",
      "Bash(git commit:*)",
      "Bash(git pull:*)",
      "Bash(./.venv/bin/python -m pdf2md:*)",
      "Bash(ls:*)", "Bash(cat:*)", "Bash(grep:*)", "Bash(tail:*)"
    ]
  }
}
```

- **許可リスト方式(推奨)**: 上記のように「安全でよく使うもの」だけ無確認にする。
- **完全自動承認(bypass)**: すべての確認をスキップするモードも存在するが、
  `rm`・上書き・`push`・外部送信も無確認になるため、**使い捨て/サンドボックス等の
  信頼できる文脈に限る**。常用は非推奨。

---

## まとめ(あなたの効率化の型)

- 重い/長い/Claude不要な作業 → **①非対話スクリプト＋スケジューラ**で無人完走(これが主役)
- 対話作業の確認削減 → **②許可リスト**
- **③完全自動承認は最終手段**。安全に「全部Yes」したいなら、①の“1機能スクリプト”に閉じ込めるのが正解
