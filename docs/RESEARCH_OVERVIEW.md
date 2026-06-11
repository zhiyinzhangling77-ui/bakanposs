# Research Overview — Mediterranean ET-SWC Decoupling Project

> 📍 **どの session でもまずここから読む**
> 全体像 → 各解析の目的 → 現在のステータス → 次のステップ

---

## 1. 研究全体の問い

**「Mediterranean 半乾燥地で、灌漑農地 (Tarazona almond) と雨養農地 (Oran cereal) の蒸発散応答はどう違うか?」**

具体的には:
- 水入力イベント(雨/灌漑)後の蒸発散 ET 回復ダイナミクス
- 表層土壌水分 SWC と ET の関係
- 衛星 ET 推定モデル(MOD16, ECOSTRESS, SMAP)の灌漑バイアスへの含意

---

## 2. 3 つの解析の関係図

```
┌──────────────────────────────────────────────────────────────┐
│  解析A (Eddy Covariance ベース) ★ COMPLETED                  │
│    ↓ τ ≈ 3-4 d は 2 サイトで普遍                              │
│    ↓ 振幅 (LE_0 - LE_∞) は灌漑で 4× scale                     │
│    ↓ Publication-ready (validated, power-checked)            │
└──────────────────────────────────────────────────────────────┘
                ↓ 検証/補強
┌────────────────────────┐  ┌────────────────────────┐
│ 解析B (衛星 ET 検証)    │  │ 解析C (NDVI phenology)  │
│   PLANNED              │  │   IN PROGRESS (v1)      │
│                        │  │                        │
│ • MOD16/ECOSTRESS で τ │  │ • NDVI で生育期客観定義 │
│   を再現できるか?       │  │ • NDVI×Rn でLE予測精度  │
│ • 灌漑バイアスの時間構造 │  │ • 仮説のフェノロジー側面 │
│ • 空間スケール拡張       │  │                        │
└────────────────────────┘  └────────────────────────┘
```

**解析A は完成**(v27 で reviewer-bulletproof)、これを **独立データで検証**するのが解析B、**別の物理量で補強**するのが解析C。

---

## 3. 各解析の最終目的

### 解析A: 地上 Eddy Covariance ベースの仮説検証 ★ 完了
**問い**: 灌漑 (drip) と 雨水 (rain) で ET 回復ダイナミクスは違うか?
**答え**:
- **τ ≈ 3.0-3.8 d は両サイト・両季節で統計的に区別不能**
- **振幅 (LE_0 − LE_∞) は灌漑が約 4× 大きい**
- 結論: 「時間スケールは普遍、振幅は management signal」
→ 詳細は `ANALYSIS_A_FINAL.md`

### 解析B: 衛星 ET プロダクトでの検証(これからやる)
**問い**: 衛星から推定される ET は同じ τ ≈ 3-4 d を再現するか?
**含意**: できれば → 普遍性が独立に確認、できなければ → モデル系統誤差の発見
→ 詳細は `ANALYSIS_B_PLAN.md`

### 解析C: NDVI による生理学的補強(v1 進行中)
**問い**: NDVI で見た植被活性が 解析A の結論を支持するか?
**含意**: NDVI dynamics で τ ≈ 3-4 d の生理メカニズムを補強
→ 詳細は `ANALYSIS_C_PLAN.md`

---

## 4. 主要データソース(共通参照)

### Eddy Covariance タワー
| サイト | 種別 | 期間 | 灌漑 | パス |
|---|---|---|---|---|
| **Oran** | rainfed winter cereal | 2018-2020 | なし | `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_EddyDaily_MASTER_2018_2020_correct.csv` |
| **Tarazona** | drip-irrigated almond | 2020-2024 | あり (~138 events) | `/home/shion-nagamine/Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv` |

座標:
- Oran: lat 38.82, lon −1.86
- Tarazona: lat 39.266, lon −1.9397

### 衛星データ(解析B 候補)
| プロダクト | 解像度 | 期間 | 用途 |
|---|---|---|---|
| MODIS MOD16A2 | 500m / 8-day | 2000-現在 | アルゴリズム ET |
| ECOSTRESS | 70m / 1-5 day | 2018-現在 | 熱赤外 ET (高解像度) |
| PML-V2 | 500m / 8-day | 2000-現在 | プロセスベース ET |
| MODIS MOD13Q1 NDVI | 250m / 16-day | 2000-現在 | 植被(解析C) |
| Sentinel-2 NDVI/NDWI | 10m / 5-day | 2017-現在 | 高解像度植被 |
| SMAP root-zone SM | 9km / 2-3 day | 2015-現在 | 土壌水分(独立) |

### 補助
- ERA5 (気象再解析): `/mnt/hdd/Dataset/ERA5_*`
- CHIRPS (降水): `/home/shion-nagamine/Dataset/Chirps v2 for persipitation`
- TerraClimate: `/home/shion-nagamine/Dataset/TerraClimate_ppt`

---

## 5. 現在のステータス(2026-05 時点)

| 解析 | 状況 | 主結果 |
|---|---|---|
| **A** | ✅ COMPLETED | τ=3.0-3.8d (validated, MDE-checked) / 振幅 4× scaling |
| **B** | ⬜ PLANNED | これから着手(satellite ET 検証) |
| **C** | 🔶 v1 in progress | NDVI 抽出 + フラックス相関(統合待ち) |

---

## 6. 重要な制約と注意(全解析共通)

### 比較禁止事項
- ❌ Oran と Tarazona の **絶対 LE 値** cross-site 比較(種・LAI・季節差で交絡)
- ❌ NDWI/NDVI の **絶対閾値**(サイト依存性強い)
- ❌ 季節間プーリングしたままの分析(artifact 生む)

### 推奨手法
- ✅ Within-site での相対指標(SDS, DRR, τ)
- ✅ Season-stratified 解析
- ✅ Bootstrap CI + MDE での power 確認
- ✅ Event-level fitting(pseudo-replication 回避)
- ✅ Validated fits のみで verdict

### Tarazona に関する重要事実
- アーモンド orchard、樹間 4.5m × 6.5m、密度 342 trees/ha
- **2 本の drip lines を樹列ごとに配置**(commercial drip irrigation)
- 灌漑量 12-15 mm/event、5-10月に集中(月 7-47 events)
- → 解析A の "灌漑後 τ" は drip wetted-bulb dynamics を反映

---

## 7. 次セッションでまずやること

1. 本ファイル `RESEARCH_OVERVIEW.md` を読む
2. 個別解析を進める前に対象解析の `ANALYSIS_X_PLAN.md` を読む
3. 解析A の最新結果が必要なら `ANALYSIS_A_FINAL.md`
4. 過去の判断履歴・失敗が必要なら `SESSION_SUMMARY.md`
5. 環境構築は `requirements.txt`

### 優先順
- 解析B が**有効性検証として最強**(独立データソース) → 次に着手推奨
- 解析C はフェノロジー軸での補強 → 解析B と並行可

---

## 8. ファイル構成

```
bakanposs/                            (Claude が作業する repo)
├── CLAUDE.md                          自動読込されるガイド
├── pyproject.toml                     パッケージメタデータ
├── docs/
│   ├── RESEARCH_OVERVIEW.md           ★ このファイル
│   ├── ANALYSIS_A_FINAL.md            解析A 最終結果
│   ├── ANALYSIS_B_PLAN.md             解析B 設計図
│   ├── ANALYSIS_C_PLAN.md             解析C 設計図
│   └── RUN_ANALYSIS_C.md              解析C 実行ガイド
├── reports/
│   ├── SESSION_SUMMARY.md             過去セッション履歴
│   └── analysis_C_report.md           解析C v1 報告
├── bakanposs/                         Python パッケージ
│   ├── loaders.py                     共通ローダ
│   ├── analysis_a.py                  解析A (v9 が最新公開版)
│   └── analysis_c/
│       ├── v1_legacy.py               解析C v1 (多目的)
│       └── v2_phenology.py            解析C v2 (集約版)
├── analyses/
│   ├── run_analysis_A.py              A エントリポイント
│   └── run_analysis_C.py              C エントリポイント (--version v1|v2)
└── output_analysis_*_v*/               実行ごとの出力 (figures + CSV, gitignore)
```

---

## 9. Publication プラン

**解析A の論文化フェーズ**:
1. 解析B と解析C で解析A を補強 →
2. 統合論文 draft →
3. 投稿

**想定 Title**:
> *"Universal post-water-input ET relaxation timescale (~3.3 d) in Mediterranean drylands, with management-scaled amplitude: implications for satellite ET retrieval"*

**想定 Journal candidates**: Agricultural and Forest Meteorology, Remote Sensing of Environment, Water Resources Research
