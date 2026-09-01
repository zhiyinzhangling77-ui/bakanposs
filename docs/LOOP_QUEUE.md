# LOOP_QUEUE — 作業キュー

**status**: `TODO` / `DOING` / `DONE` / `BLOCKED(GATE-xx)`
**種別**: `AUTO`（クラウドの Claude だけで完結） / `GATE`（人間待ち）
**出典**: `SESSION_SUMMARY.md` §6 (U1–U15), §7 (優先度)

> 毎周の最後に必ず更新してコミットすること。**このファイルが進捗の唯一の記録**であり、会話履歴は記録として信頼しない。

---

## AUTO — 自走で進められる（これを先に食い潰す）

| ID | 内容 | 由来 | 優先 | status | 備考 |
|---|---|---|---|---|---|
| Q01 | `analysis_A_v15.py` 執筆：`in_band()` を `(lo≤0≤hi) AND (hi−lo<0.5) AND n_s≥30` に修正 | U1 | ★★★ | TODO | 実走は GATE-02 |
| Q02 | A v11–v14 を `data_loaders` に統合（v10 の commit `2ade97e` と同パターン） | U3 | ★★ | TODO | grep で旧ローダ呼び出しを検出してから |
| Q03 | 解析B のファイル索引 `docs/analysis_B_index.md` 作成 | U4 | ★★ | TODO | リネームは GATE-15 待ち。索引だけなら AUTO |
| Q04 | Abstract 草稿（250 語想定）を `paper_methods_results.md` から起こす | U7 | ★★★ | TODO | 主張4 は GATE-13 の結論待ちなので**主張1–3 のみで一旦書く** |
| Q05 | `fig_H1/H4/H6` のラベル・凡例・n 表示を整備するスクリプト | U8 | ★★★ | TODO | dpi/サイズは GATE-18 待ち。文言整備は仕様非依存 |
| Q06 | 引用の BibTeX 化（Mu 2011, Zhang 2019, Trigo 2018, Reichle 2018, Pettorelli 2005, Burnham & Anderson 2002） | U9 | ★★ | TODO | 出力スタイルは GATE-19 待ちだが .bib 自体は非依存 |
| Q07 | `reports/analysis_C_report.md §8.6 / §8.7` の表枠・文章構造を作る（**値は空欄**） | U5,U6 | ★★ | TODO | 値の充填は GATE-04 待ち。**空欄を推測で埋めない** |
| Q08 | ポスター原稿の構成案（`SESSION_SUMMARY §10.9` の 6 パート） | — | ★★★ | TODO | サイズ・言語は GATE-12 待ち |
| Q09 | H2（灌漑タイプ別 τ）のサイト選定基準 + τ-fit スクリプト | U10 | ★★ | TODO | データは GATE-07。`§10.3-8` bin 中央値 fit 禁止、τ∈(0,60] |
| Q10 | H7（SDS 広域マッピング）の抽出設計 | U12 | ★ | TODO | `§10.3-9` SMAP は 9 km grid で |
| Q11 | `scripts/gee_extract.js` 改修（SMAP 9 km 対応、non-ASCII 除去） | — | ★★ | TODO | 実行は GATE-06。`§10.3-10, -11` 参照 |
| Q12 | U14 スクリプト：TzM の lag>14d のみで EF 集計 | U14 | ★ | TODO | 実走は GATE-01 |

## BLOCKED — 人間待ち（`docs/HUMAN_GATES.md` 参照）

| ID | 内容 | status |
|---|---|---|
| Q13 | A v15 の実データ再走と結果反映 | BLOCKED(GATE-02) |
| Q14 | Oran 新ローダでの 3 値再走（SDS +0.43 / n=203 r=+0.80 等） | BLOCKED(GATE-03) |
| Q15 | C 分析8/9 の実行と §8.6/8.7 充填 | BLOCKED(GATE-04) |
| Q16 | 中間ファイルの所在確認 | BLOCKED(GATE-05) |
| Q17 | GEE 抽出の実行 | BLOCKED(GATE-06) |
| Q18 | FLUXNET/AmeriFlux データ取得 | BLOCKED(GATE-07) |
| Q19 | 主張4 を維持するかの判断 | BLOCKED(GATE-13) |
| Q20 | 投稿先確定 → 語数・図仕様・引用スタイル | BLOCKED(GATE-12) |
| Q21 | 多深度 SWC センサー照会 | BLOCKED(GATE-10) |

---

## 完了ログ

| 周 | 日付 | ID | コミット | 一行要約 |
|---|---|---|---|---|
| — | — | — | — | （まだなし） |
