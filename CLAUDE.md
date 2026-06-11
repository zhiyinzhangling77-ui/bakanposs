# Repository Guide for Claude

> ⚠️ **このセッションで最初にやること**:
>
> 1. **[`docs/RESEARCH_OVERVIEW.md`](./docs/RESEARCH_OVERVIEW.md)** — 研究全体図、3 つの解析の関係(必読)
> 2. **[`reports/SESSION_SUMMARY.md`](./reports/SESSION_SUMMARY.md)** — 過去セッションの全判断履歴・失敗・制約
> 3. 進める解析に応じて以下のいずれか:
>    - 解析A 結果参照: [`docs/ANALYSIS_A_FINAL.md`](./docs/ANALYSIS_A_FINAL.md)
>    - 解析B(衛星 ET): [`docs/ANALYSIS_B_PLAN.md`](./docs/ANALYSIS_B_PLAN.md)
>    - 解析C(NDVI): [`docs/ANALYSIS_C_PLAN.md`](./docs/ANALYSIS_C_PLAN.md)

---

## クイックステータス(2026-05 時点)

| 解析 | 状態 | 主結果 |
|---|---|---|
| **A** | ✅ COMPLETED (v27) | τ ≈ 3.0-3.8d universal、振幅 4× scaling |
| **B** | ⬜ PLANNED 次着手 | 衛星 ET(MOD16/ECOSTRESS)で τ 検証 |
| **C** | 🔶 v1 in progress | NDVI で active period 客観定義 |

→ **次セッションの優先タスク**: 解析B(衛星 ET 検証)着手

---

## 主要ファイル

```
bakanposs/                                  ← repo root
├── CLAUDE.md                                ★ このファイル(自動読込)
├── pyproject.toml                           パッケージメタデータ
├── docs/
│   ├── RESEARCH_OVERVIEW.md                 ★ 全体図
│   ├── ANALYSIS_A_FINAL.md                  解析A 最終結果
│   ├── ANALYSIS_B_PLAN.md                   解析B 設計図
│   ├── ANALYSIS_C_PLAN.md                   解析C 設計図
│   └── RUN_ANALYSIS_C.md                    解析C 実行ガイド
├── reports/
│   ├── SESSION_SUMMARY.md                   過去履歴
│   └── analysis_C_report.md                 解析C v1 報告
├── bakanposs/                               Python パッケージ
│   ├── loaders.py                           共通ローダ
│   ├── analysis_a.py                        解析A (v9 最新)
│   └── analysis_c/
│       ├── v1_legacy.py                     解析C v1 (多目的)
│       └── v2_phenology.py                  解析C v2 (集約版)
└── analyses/
    ├── run_analysis_A.py                    A エントリポイント
    └── run_analysis_C.py                    C エントリポイント (--version v1|v2)
```

---

## 絶対に再試行禁止(過去失敗確認済)

1. ❌ Oran と Tarazona の **絶対 LE/ET の cross-site 比較**(種・LAI・季節差で交絡、F2)
2. ❌ **NDWI 絶対閾値**での deep_access 抽出(F1)
3. ❌ **季節間プーリング**での SDS 計算(F3, artifact)
4. ❌ **n < 30** のサンプルで CI が 0 を含む = 深根支持 判定(F6, 偽陽性)
5. ❌ Fit τ ≥ 50d を有効値として扱う(v26 で発覚、v27 で validation 追加)
6. ❌ `LE_0_FLOOR` 一律 50(Oran 雨養 cereal の LE 27 W/m² で fit 不可)
   → adaptive bound 必須

---

## 確立された結論(解析A、変更厳禁)

### Tarazona は灌漑あり(ユーザ情報、v13 で組込)
- アーモンド orchard、drip irrigation
- 灌漑期: 5-10月(月 7-47 events、12-15 mm/event)
- 非灌漑期: 11-4月

### τ の最終値(v27 validated)
- Oran active (Nov-Jun, n=10): τ=2.82d [1.82, 5.34]
- Tarazona active (Jun-Sep, n=41): τ=3.36d [2.46, 4.90]
- Pairwise diff < MDE = 統計的に区別不能
- → "τ ≈ 3-4d effective universal in Mediterranean drylands"

### Amplitude (LE_0 − LE_∞)
- Oran ~23 W/m² vs Tarazona ~94 W/m²
- **~4× scaling は management signal**

---

## 作業前に確認すべきこと(ユーザに尋ねる)

- 解析B 着手時:
  - MOD16A2 / ECOSTRESS データは取得済みか?
  - GEE アカウントの有無
  - 既存のフェッチコードの所在
- 解析C 着手時:
  - MOD13Q1 NDVI は取得済みか?
  - `bakanposs/analysis_c/v1_legacy.py` を baseline にするか、`v2_phenology.py` を拡張するか
- 論文化フェーズ:
  - Target journal
  - 共著者の有無

---

## コーディング規約(過去セッションから)

1. **コメントは "なぜ" を書く**(what は code を読めばわかる)
2. **自己説明的 figure**(タイトルに結論、軸ラベル明確、サンプル数表示)
3. **段階的バージョン**(v_X+1 は v_X の問題を1つ解決)
4. **CSV + figure を必ずペア出力**
5. **失敗 fit を validation で除外**(τ at boundary, R²<0.3, n<3)
6. **bootstrap CI を必ず併記**(point estimate だけは禁止)
7. **MDE で power 明示**(差が「ない」結果は MDE 計算してこそ意味)

---

## Tone

- 自分の結果を疑う(反論を先取りして潰す)
- ユーザに「これでいいか?」を確認してから先に進む(特に方針転換時)
- 失敗を隠さない(間違ったときは正直に認める)
- 過剰な楽観・過剰な悲観を避ける
