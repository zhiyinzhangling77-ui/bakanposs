# Repository Guide for Claude

> ⚠️ **このセッションでまず最初にやること**:
> [`SESSION_SUMMARY.md`](./SESSION_SUMMARY.md) を Read tool で**全文読んでから**作業を開始してください。
> 過去セッションの全判断履歴・失敗・未解決課題・制約が記載されています。
> これを読まずに作業を進めると、既に失敗確認された手法を再試行してしまう恐れがあります。

## クイックリンク
- 詳細ハンドオフ: `SESSION_SUMMARY.md`
- 解析A 最終版: `analysis_A_v14.py`
- 共通ローダ: `data_loaders.py`(v9-v11 で使うバグなし版)
- 並行解析C: `analysis_C_v1.py`

## 絶対にしてはいけないこと(再試行禁止)
1. Oran と Tarazona の **絶対 LE/ET の cross-site 比較**(種・LAI・季節すべて違う、`SESSION_SUMMARY.md` §5 F2)
2. **NDWI 絶対閾値** での `deep_access` 抽出(F1)
3. **季節間プーリング**での SDS 計算(F3, artifact を生む)
4. **n < 30 のサンプルで "CI が 0 を含む = 深根支持" 判定**(F6, 偽陽性)

## 現在のスタンス
- 当初の "Tarazona 深根" 仮説は **棄却寄り**
- 真相は **灌漑による 3-4 日の SWC-ET decoupling**(v14 で判明)
- 論文化方向は "Drip irrigation decouples surface SWC from canopy ET"

## 作業前に確認すべきこと
- ユーザーに優先タスクを尋ねる:
  - v15(verdict bug fix + recovery time analysis)?
  - 論文骨子の draft?
  - 衛星補強(MODIS LST / ECOSTRESS / Sentinel-1)?
  - 解析A と解析C の統合?
