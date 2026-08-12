# 時間のかかる解析の実行シート（あなたのマシンで夜間・放置実行）

このリポジトリの重い解析を、**データのある自分のマシン**（`/mnt/hdd/...`）で、
`nohup ... &` で起動して放置するためのコマンド集。すべて**ログ保存・再開/隔離あり**。
※このリモート環境にはデータが無いので実行はローカルで。

前提：リポジトリ直下で `python -m japanflux_pn.<module>`。出力は各 `--outroot` に CSV。
進捗は `tail -f <log>` で確認。

---

## 優先度つき（上ほど価値・重さのバランスが良い）

### ① fig2b のビン数感度スイープ（Phase 1・新規）
「m=11 を都合よく選んだ」を封じる。全ビンで flux 群が崩れ thermal 群が残れば結論不変。
```bash
nohup python -m japanflux_pn.run_sensitivity --site JP-Tak \
  --years 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 \
          2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 \
  --bins 7 9 11 13 15 --month 7 8 \
  --outroot ~/bakanposs/japanflux_pn/outputs_sensitivity \
  > sens_Tak.log 2>&1 &
tail -f sens_Tak.log
```
→ 出力 `sensitivity_bins_JP-Tak.csv`（列：n_bins, flux平均drop%, thermal平均drop%, separation）。
最後に「✅ 全ビンで結論不変」か「⚠ 崩れる」を自動判定。

### ② 非線形（CMIknn）での因果骨格の再現（Phase 1・最重量）
線形 ParCorr で出た骨格が非線形でも残るか。**一番時間がかかる**（サロゲート×近傍探索）。
```bash
# サイトごとに（森林から）。夜通し想定。
nohup python -m japanflux_pn.run_robustness --site JP-Tak --test cmiknn \
  --tau-max 6 --sig-samples 200 \
  --outroot ~/bakanposs/japanflux_pn/outputs_robust_cmiknn \
  > robust_cmiknn_Tak.log 2>&1 &
tail -f robust_cmiknn_Tak.log
```
（`--test parcorr` 版を先に回して基準を作り、cmiknn 版と骨格を突き合わせる。）

### ③ 完全版 PCMCI+（tau_max=36・sig=500・CMIknn、チェックポイント付き）
最も厳密な単発ネットワーク。落ちても済んだサイトは残る。
```bash
nohup python -m japanflux_pn.run_causal_all --test cmiknn \
  --outroot ~/bakanposs/japanflux_pn/outputs_pcmci_full \
  > pcmci_overnight.log 2>&1 &
tail -f pcmci_overnight.log
```

### ④ O-information バッチ（全サイト・高次の相乗/冗長）
```bash
nohup python -m japanflux_pn.batch_oinfo --all \
  --csv ~/bakanposs/japanflux_pn/outputs_oinfo/batch_oinfo.csv \
  > oinfo_all.log 2>&1 &
tail -f oinfo_all.log
```

### ⑤ ラグ（tau_max）感度：因果骨格が窓長に依存しないか
```bash
for TM in 4 6 12; do
  nohup python -m japanflux_pn.run_robustness --site JP-Tak --test parcorr \
    --tau-max $TM --outroot ~/bakanposs/japanflux_pn/outputs_tau_$TM \
    > robust_tau_$TM.log 2>&1 &
done
```

---

## 運用のコツ
- **1つずつ確かめてから並列に**：まず①を回して形式・出力を確認 → 重い②③を夜に。
- **CPU 数に注意**：②③は重いので同時起動は 1〜2 本まで。`nproc` を確認。
- **落ちても平気**：③④はサイト単位でチェックポイント。②は年単位で自動スキップ。
- **結果の見方**：
  - ① `separation > 0` が全 m で成り立てば「ビン非依存」。
  - ② parcorr と cmiknn で**同じコアリンク**が残れば「線形に依らず頑健」。
  - ⑤ tau を変えてもコアが不変なら「窓長に依らず頑健」。
- これらが揃うと Phase 1（結果を"堅く"する）の主要な頑健性チェックが完成。

## 出力の回収
各 `--outroot` の CSV をこのリポジトリにコピー → `scripts/make_slides_figs.py` の
数値を実測に差し替えれば、図が「埋め込み値」から「実データ」に更新できる。
（どの CSV をどの図に使うかは相談してくれれば対応します。）
