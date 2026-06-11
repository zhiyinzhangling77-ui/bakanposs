# Repository Guide for Claude

> ⚠️ **このセッションで最初にやること(順番厳守)**:
>
> 1. **[`RESEARCH_OVERVIEW.md`](./docs/RESEARCH_OVERVIEW.md)** — 研究全体図、3 つの解析の関係
> 2. **[`SESSION_SUMMARY.md`](./docs/sessions/SESSION_SUMMARY.md)** — 旧セッション履歴(v9-v27)
> 3. **[`SESSION_UPDATE_v28_v32.md`](./docs/sessions/SESSION_UPDATE_v28_v32.md)** — v28-v32 poster figure + 解析B 橋渡し
> 4. **★[`SESSION_UPDATE_v32_blindspot.md`](./docs/sessions/SESSION_UPDATE_v32_blindspot.md)** — **最新議論(絶対参照)**: τ_bias の位置づけ訂正、ETv3 mechanism 訂正、SMAP/H26 評価、narrative 転換、Tarazona サイズ訂正(1 ha)
> 5. **[`SATELLITE_ET_NOTES.md`](./docs/SATELLITE_ET_NOTES.md)** — 衛星 ET caveats(§4 SM 入力ありに訂正済)
> 6. 進める解析に応じて:
>    - 解析A 結果参照: [`ANALYSIS_A_FINAL.md`](./docs/ANALYSIS_A_FINAL.md)
>    - 解析B(衛星 ET): [`ANALYSIS_B_PLAN.md`](./docs/ANALYSIS_B_PLAN.md)
>    - 解析C(NDVI): [`ANALYSIS_C_PLAN.md`](./docs/ANALYSIS_C_PLAN.md)

---

## クイックステータス(2026-05 時点、v32 blindspot session 後)

| 解析 | 状態 | 主結果 |
|---|---|---|
| **A** | ✅ COMPLETED (v31 poster figure) | τ ≈ 3.0-3.8d universal、振幅 4.5× scaling |
| **B** | 🟡 v32 + bias_stats 完了、次は narrative 反映 (v33) | Tarazona: r=0.07, bias −70%。amp_EC=95 vs amp_Sat≈0 |
| **C** | 🔶 v1 in progress | NDVI で active period 客観定義 |
| **Poster** | 🟡 A0 縦テンプレ生成済、v2 narrative 反映待ち | `poster/build_poster_template.py` |

→ **次セッションの最初**: `SESSION_UPDATE_v32_blindspot.md` §10 の決定肢 (A)–(D) のどれで進めるかをユーザーに確認

---

## 主要ファイル(Phase 1+2+3 整理後、2026-06)

```
/home/user/bakanposs/
├── CLAUDE.md                            ★ このファイル(自動読込)
├── requirements.txt                     Python 環境
├── data_loaders.py                      共通ローダ
├── run_analysis_C.py                    C runner
│
├── docs/                                ★ 全マークダウン文書
│   ├── RESEARCH_OVERVIEW.md             ★ 研究全体図
│   ├── ANALYSIS_A_FINAL.md              解析A 最終結果
│   ├── ANALYSIS_A_FAQ.md                解析A FAQ
│   ├── ANALYSIS_B_PLAN.md               解析B 設計図
│   ├── ANALYSIS_C_PLAN.md               解析C 設計図
│   ├── SATELLITE_ET_NOTES.md            衛星 ET caveats
│   ├── REPORT_to_site_collaborators.md  site PI 向け報告
│   ├── paper_methods_results.md         論文 draft
│   ├── paper_outline.md                 論文 outline
│   ├── RUN_ANALYSIS_C.md                C 実行手順
│   └── sessions/                        セッション履歴
│       ├── SESSION_SUMMARY.md
│       ├── SESSION_UPDATE_v28_v32.md
│       └── SESSION_UPDATE_v32_blindspot.md
│
├── analysis_A_v31.py                    ★ 現役 (poster Fig 4 main)
├── analysis_A_v32.py                    ★ 現役 (Tarazona blind-spot)
├── analysis_B_v3_bias_tau.py            ★ 現役 (bias recovery pool)
├── analysis_B_v6_driver_attribution.py  ★ 現役 (driver bars)
├── analysis_C_v2_ndvi_phenology.py      ★ 現役 (NDVI phenology)
│
├── archive/                             旧 script の保管(履歴保持)
│   └── scripts/
│       ├── analysis_A_v10–v30.py        21 版
│       ├── analysis_A_v13_1_patch.py
│       ├── analysis_B_v1, v2, v4.py
│       └── analysis_C_v1.py
│
├── data/                                入力データ
├── output/                              ★ 全出力を統合
│   ├── analysis_A/v15–v30/              旧 output
│   ├── analysis_A/v31/, v32/            ★ 現役出力
│   ├── analysis_B/v1, v2, v4/           旧 output
│   ├── analysis_B/v3/, v6/              ★ 現役出力
│   ├── analysis_C/v1/, v2/, last/
│   └── bias_stats/                       EC vs ETv3 scatter
│
├── scripts/                             helper scripts
├── poster/                              poster template
└── reports/                             reports
```

★ パス更新: 全現役 script は `./output/analysis_*/v*/` を default に。
   `analysis_A_v32.py` は `./output/analysis_B/v3/v3_bias_*.csv` から bias pool を読む。

---

## 絶対に再試行禁止(過去失敗確認済)

1. ❌ Oran と Tarazona の **絶対 LE/ET の cross-site 比較**(種・LAI・季節差で交絡、F2)
2. ❌ **NDWI 絶対閾値**での deep_access 抽出(F1)
3. ❌ **季節間プーリング**での SDS 計算(F3, artifact)
4. ❌ **n < 30** のサンプルで CI が 0 を含む = 深根支持 判定(F6, 偽陽性)
5. ❌ Fit τ ≥ 50d を有効値として扱う(v26 で発覚、v27 で validation 追加)
6. ❌ `LE_0_FLOOR` 一律 50(Oran 雨養 cereal の LE 27 W/m² で fit 不可)
   → adaptive bound 必須
7. ❌ **τ_bias を独立な物理測定として扱う**(v32 blindspot 議論)
   → bias = EC − Sat ≈ EC − const なので τ_bias = τ_EC は数学的構造、
     新規情報ではない。consistency check のみ。
8. ❌ **SMAP 9 km で空間希釈論を実証**(v32 blindspot 議論)
   → orchard 1 ha との比較で循環論証。in-situ tower SWC を ground truth
     に据えること。
9. ❌ **「ETv3 は降水だけで ET 計算」と書く**(v32 blindspot 議論)
   → H-SAF H141/H142/H26 の衛星 SM を入力に持つ。drip blind spot は
     SM 入力チェーンの 3 重の解像度限界による(空間/物理/SVAT 構造)。
10. ❌ **Tarazona orchard を「0.1 km²」と書く**
    → 正しくは **1 ha = 0.01 km²**。過去メモの 0.1 km² は誤り。
11. ❌ **「Meteosat ET は不正確」と総括**(v32 blindspot 議論)
    → rainfed Oran では r=0.82, bias +1%, KGE=0.61 で良好。drip 特異の
      blind spot として書くこと。

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
  - `analysis_C_v1.py` を baseline にするか書き直すか
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
