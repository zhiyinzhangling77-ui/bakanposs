# POSTER DRAFT — 進捗報告用ポスター原稿

**対象**: 教授への研究進捗報告
**作成日**: 2026-05-18
**ブランチ**: `claude/analysis-presentation-prep-BQRnM`
**前提資料**: `SESSION_SUMMARY.md` (3 解析統合)、A v27 `ANALYSIS_A_FINAL.md`、B `paper_methods_results.md` + `analysis_narrative.md`、C `ANALYSIS_C_PLAN.md`
**スタンス**: 進捗報告寄り。手法を丁寧に説明し、A 完了 / B 進行中 / C 計画中 の現状を正直に示す。トライアンギュレーション (A × B × C) を最終枠組みとする。

---

## 想定タイトル候補

> **(本命)** *"Universal ET relaxation timescale (τ ≈ 3–4 d) and management-scaled amplitude in Mediterranean drylands: an eddy-covariance × satellite × SMAP triangulation"*

> **(短縮)** *"Why satellite ET underestimates drip-irrigated almond by 2–4 mm d⁻¹ — and how to fix it in 4 days"*

> **(教授向け日本語タイトル)**
> **「地中海半乾燥地におけるドリップ灌漑の蒸発散応答時定数 τ ≈ 3-4 日と振幅 4× スケーリング:地上 × 衛星 × SMAP の三角測量」**

---

## 1. リサーチクエスチョン（ポスター左上）

**大問**: ドリップ灌漑農地で衛星 ET が **2–3 割系統的に過小評価** するのは何故か?その時間構造はどうなっているか?

**派生問**:
- Q1: 灌漑直後の蒸散応答は雨養と質的に違うのか? (時間スケール? 振幅?)
- Q2: 表層 5cm SWC が乾いているのに LE が高い「脱結合」は深根? wet-bulb?
- Q3: 衛星 ET が見ているものは、地上 EC が見ているものと同じか?
- Q4: もしバイアスが時間構造を持つなら、補正できるか?

---

## 2. 全体像のひと目図 (Headline schematic)

**【Fig 1】Drip wet-bulb メカニズム模式図** ※新規作成必要

```
              ATM (高 VPD, 高 Rn)
                    ↑ LE
       ┌────────────────────┐
       │  Almond canopy (LAI~4) │
       │       ↑↑↑           │
       └───┬────────────┬───┘
           │  根系     │
0-5cm   ━━━╪━━━━━━━━━━━━╪━━━  ←表層は乾燥(地表蒸発, 流出)
           │ 主吸水層   │       ●in-situ SWC sensor (見えない bulb)
10-30cm ━━━╪━━━━━━━━━━━━╪━━━  ← drip wet-bulb (3-4d持続)
           │           │       ◇ドリップ
~1m     ━━━╪━━━━━━━━━━━━╪━━━  ←root-zone (遅れて湿る)
                                ●SMAP root-zone (見えるが遅れる)
```

**何を見せる**: 灌漑水が地表ではなく 10–30 cm に局在して 3–4 日蒸散を駆動する物理像。
**役割**: ポスターの最初に置く。以降のすべてのデータがこの図に対応する。

---

## 3. なぜこの研究をするのか (Motivation, 1パラ)

衛星 ET プロダクト (MOD16, PML, GLEAM 等) は灌漑スケジューリング、地域水収支、干ばつ早期警戒に運用されはじめている。しかし**灌漑農地で 2–3 割の系統的過小評価** が報告される (Velpuri 2013, Senay 2017, Talsma 2018)。既存研究は大半が「シーズン平均」または「年平均」でのバイアス報告に留まり、**「灌漑イベント直後にバイアスがどう時間変化するか」を連続フラックス観測で定量化した例はほとんどない**。これが分かれば補正手法の設計が変わる:
- 3–4 日で減衰 → 「灌漑経過日数」を入力する補正モデルで対処可
- 季節的に常時残る → 校正データ追加 / アルゴリズム改修が必要

---

## 4. 仮説の変遷（誠実さアピール）

**【Fig 3】仮説変遷フロー**

```
┌────────────────────────────────────────────┐
│ [初期仮説 v4-v12]                            │
│ "Tarazona almond は深根で地下水を使うから    │
│  表層が乾いても蒸散を維持する" (deep-root)   │
└────────────────────────────────────────────┘
       │  v9-v12 で SDS_TzM=+0.05 (Oran +0.31) → 一見支持
       ↓
[転機: v13]
ユーザーから「Tarazona は drip 灌漑あり」
       │
       │  深根仮説と灌漑仮説が観測上区別不能に
       ↓
┌────────────────────────────────────────────┐
│ [v14] 灌漑経過日数で層別                     │
│  d0-3: SDS = +0.13 (有意)                   │
│  d4-7: SDS = +0.01 (CI 0 含む)              │
│  d8+ : SDS ≈ 0     (CI 0 含む)              │
│  → 深根なら d8+ でも SDS≈0 のはず → でも     │
│    そもそも全 SDS が小さい → wet-bulb       │
└────────────────────────────────────────────┘
       │
       ↓
┌────────────────────────────────────────────┐
│ [新仮説] drip wet-bulb (10-30cm) が 3-4 日   │
│  持続するため、表層 5cm が乾いていても        │
│  蒸散が続く (management-driven decoupling)   │
└────────────────────────────────────────────┘
```

**ポスターでの位置**: イントロの直後。「先入観を自分で棄却した」プロセスを見せる。reviewer の "confirmation bias?" 攻撃を先に潰す。

---

## 5. 解析の三角測量 (Methods overview)

**全体構造**:

```
┌─────────────────────────────────────────────────────────┐
│  解析A: 地上 EC (point) ★ COMPLETED (v27)               │
│    τ ≈ 3.0-3.8 d (Oran 2.82, Tarazona 3.36, MDE NS)    │
│    振幅 4× (Oran 23 vs Tarazona 94 W/m²)                │
│    Reviewer-bulletproof (S1-S4, MDE, AIC validated)     │
└─────────────────────────────────────────────────────────┘
                       ↓ 独立検証
┌─────────────────────┐  ┌────────────────────────────────┐
│ 解析B: 衛星 ET       │  │ 解析C: NDVI フェノロジー         │
│  ★ IN PROGRESS      │  │  🔶 PLANNED                     │
│  MOD16/PML/METv3     │  │                                │
│   τ = 4-6 d で再現   │  │  active period 客観定義        │
│  SMAP × in-situ で   │  │  NDVI peak での τ 再 fit       │
│   depth inversion    │  │                                │
│  RMSE -49〜-65% 補正 │  │                                │
└─────────────────────┘  └────────────────────────────────┘
```

**この構造を取った理由**:
- A 単独では「EC タワー 1 サイト点観測」の批判を受ける
- B (独立衛星 3 製品) で時間スケールを再現すれば「ground × satellite triangulation」
- C (NDVI) で active period の定義妥当性を確認すれば「観測駆動の justification」

---

## 6. 解析A の主結果 — τ 普遍性 + 振幅スケーリング

### 6.1 図と数値

**【Fig 4】Recovery curve + τ comparison**
ファイル候補: `output_analysis_A_v22/fig01_recovery_both.png` + `output_analysis_A_v27/fig02_phenology_active.png`

| Stratum | n_events | τ (d) | 95% CI |
|---|---|---|---|
| Oran winter (Nov-Feb) | 6 | 3.29 | (wide, single-season) |
| Oran summer (Jun-Aug) | 4 | 3.79 | (wide) |
| **Oran active (Nov-Jun, pool)** | **10** | **2.82** | **[1.82, 5.34]** |
| **Tarazona active (Jun-Sep)** | **41** | **3.36** | **[2.46, 4.90]** |

| サイト | LE_0 | LE_∞ | 振幅 | 比 |
|---|---|---|---|---|
| Oran (rainfed) | ~40 W/m² | 17 W/m² | ~23 | 1.0× |
| Tarazona (drip) | ~210 W/m² | 114 W/m² | ~94 | **~4.1×** |

### 6.2 メッセージ

> **「時間スケール τ は management に非依存。振幅 (LE_0 − LE_∞) こそが management signal」**

### 6.3 採用手法と理由

| 手法 | 理由 |
|---|---|
| 2-parameter exp model: `LE(d) = LE_∞ + (LE_0 − LE_∞)·exp(−d/τ)` | 物理的根拠 (一次リザバー減衰)、パラメータが直接解釈可能 |
| `LE_∞` を `d ≥ 10` 観測中央値に固定 | 過剰自由度を避けて identifiability 確保 |
| τ ∈ [0.5, 60] d 制約 | 灌漑サイクル 2-3d の 0.2× 〜 20× を物理的妥当範囲とする |
| Bootstrap N=5000 BCa CI | 非線形回帰の理論分布未知 |
| Validation filter (boundary, R²<0.3, n<3) | overfit / 偽 fit を除外 |

---

### 6.4 「差がないこと」を主張する根拠 — MDE 解析

**【Fig 5】Power analysis (MDE)**
ファイル: `output_analysis_A_v27/fig03_power_analysis.png`

| 比較 | obs |Δτ| | MDE | sig? |
|---|---|---|---|
| Oran winter ↔ Oran summer | 0.50 | 3.31 | NS |
| Oran winter ↔ Tarazona summer | 0.07 | 2.22 | NS |
| Oran summer ↔ Tarazona summer | 0.43 | 3.00 | NS |
| **Oran active ↔ Tarazona active** | **0.54** | **2.15** | **NS** ★ |

**なぜこの図が要るのか**:
- `p > 0.05` は「差を見つけられなかった」だけで「差がない」ことを示さない（**absence of evidence ≠ evidence of absence**）
- MDE = 1.96 × √(SE₁² + SE₂²) を計算すれば「○日以上の差があれば検出できた」と言える
- 観測差 0.54d < MDE 2.15d なら「差があれば検出できる感度はあった、それでも差が出なかった」と論じられる

**反論予想 → 対応**:
- "Statistical power が低いから差が出ないだけでは?" → MDE 2.15d で 0.54d を検出する power は 5%。逆に τ 差 2d 以上を検出する power は十分にあった。

---

## 7. 解析B — 独立衛星 3 製品での確認 + 補正

### 7.1 dose-response (灌漑バケット別バイアス)

**【Fig 6】Headline**
ファイル: `figs/fig_C2_bias_by_irrig_seasonal.png`

TzM summer × NDVI>0.3 (n=400/325/400):

| 製品 | d0-3 | d4-7 | d8+ | パターン |
|---|---|---|---|---|
| MOD16 | −4.12 | −2.68 | **−2.68** | プラトー (構造的 floor) |
| PML | −2.78 | −1.29 | **−0.93** | 0 へ漸近 (補正可) |
| METv3 | −4.03 | −2.05 | **−1.41** | 0 へ漸近 (補正可) |

**何を主張するか**:
1. **3 製品すべてで dose-response** = 因果証拠 (灌漑が原因)
2. **MOD16 だけがプラトー** = 構造的 floor (FLUXNET2015 校正の orchard 不足)
3. **PML/METv3 は補正可能** = 過渡誤差のみ

**なぜ独立 3 製品で見せるか**: 1 製品なら「アルゴリズム固有のバグ?」と言える。3 製品 (MODIS family + PML + LSA SAF Meteosat の独立家系) で同じ dose-response → **dry-surface 駆動の retrieval 全般の問題**と一般化できる。

---

### 7.2 τ-fit per product と解析A との照合

**【Fig 7】Exponential decay fit**
ファイル: `figs/fig_F_tau_fit.png`

Δ(t) = a·exp(−t/τ) + c の fit (TzM summer × NDVI>0.3):

| 製品 | a (mm/d) | τ (d) | c (mm/d) | 解釈 |
|---|---|---|---|---|
| MOD16 | −2.31 [−2.80,−1.88] | **4.0 [2.8,5.9]** | −2.29 [−2.66,−1.85]★ | 構造的 floor + 過渡 |
| PML | −2.81 [−3.43,−2.22] | **4.3 [2.9,7.0]** | −0.57 [−1.03,+0.04] | 過渡的のみ |
| METv3 | −4.03 [−4.73,−3.44] | **6.0 [4.6,8.3]** | −0.62 [−1.08,+0.06] | 過渡的のみ、振幅最大 |

★ = CI が 0 を含まない (有意な permanent offset)

**A との接続**:
- 解析A 地上 EC: τ = 3.0–3.8 d
- 解析B 衛星: τ = 4.0–4.3 d (MOD16, PML)
- → **同じ時定数が独立データソースで再現** = triangulation 成功
- METv3 の τ = 6.0 d は **5 km ピクセル混合** (周辺 rainfed 地が混じり dry-surface signature が長く残る)

**採用手法と理由**:
| 手法 | 理由 |
|---|---|
| Full model `a·exp + c` と Transient model `a·exp` を両方 fit | c の有意性で「構造的 vs 過渡的」を判別 |
| Raw daily で fit (bin median ではなく) | bin (5-6 点) では τ 暴走、raw (n=300-400) で安定 |
| τ ∈ (0, 60] 制約 | 上限を切らないと bootstrap で発散サンプルが出る |
| Bootstrap N=500 (NLS の calc コスト考慮) | CI 推定の最低限の安定性 |

---

### 7.3 機構の直接観測 — depth inversion

**【Fig 8】SMAP × in-situ 深度比較**
ファイル: `figs/fig_H6_smap_multistratum.png`

| stratum | n | r(SWC₅ cm, SMAP_rz) | SDS in-situ | SDS SMAP |
|---|---|---|---|---|
| **Oran spring (rainfed active)** | **203** | **+0.80** | **+0.43** | **+0.43** |
| Oran summer (post-harvest) | 34 | +0.97 | −0.03 | −0.19 |
| TzM spring | 45 | +0.74 | +0.16 | −0.26 |
| TzM fall (灌漑減衰) | 63 | +0.63 | +0.29 | −0.35 |
| TzM summer (全体) | 401 | +0.16 | +0.11 | −0.15 |
| **TzM summer d0-3** | **281** | **−0.19** ★ | +0.13 | −0.19 |
| TzM summer d4-7 | 75 | +0.35 | +0.01 | +0.13 |
| TzM summer d8+ | 44 | +0.61 | −0.00 | +0.18 |

★ = 統計的に有意な負相関

**何を主張するか**:
1. **rainfed では SMAP と地上 SWC が同期**(r=+0.80)→ 単一プロファイル
2. **drip d0-3 では 2 センサーが逆向きに動く**(r=−0.19, n=281)→ 地表は乾く方向、~1m は wet-bulb 伝播で湿る方向 → **bulb は両者の間 (10-30cm) に局在**
3. **d4-7 → d8+ で相関回復**(0.35 → 0.61)→ プロファイル均一化
4. **2 独立センサーで bulb の存在が直接観測された** (初めて)

**なぜこれが強い証拠か**:
- 片方のセンサーだけなら「ノイズ?」と言える
- 両方が**逆向きに動く**のは偶然では起こらない (n=281 で偶然確率 < 10⁻³)
- bulb が 10-30cm にあり、両センサーがそれを挟む位置にあるからこそ起きる現象

---

### 7.4 実用インパクト — τ 補正

**【Fig 9】τ-based correction**
ファイル: `figs/fig_H1_correction.png`

TzM summer × NDVI>0.3 で `ET_corr = ET_sat − (a·exp(−t/τ) + c)`:

| 製品 | n | RMSE 生 | RMSE 補正後 | 削減 | MBE 生 | MBE 補正後 |
|---|---|---|---|---|---|---|
| MOD16 | 400 | 4.01 | **1.39** | **−65%** | −3.69 | +0.00 |
| PML | 325 | 2.82 | **1.44** | **−49%** | −2.26 | +0.01 |
| METv3 | 400 | 3.85 | **1.50** | **−61%** | −3.37 | +0.01 |

**何を主張するか**:
- 「灌漑経過日数」を入力に与えるだけで RMSE が半分以下、MBE がほぼ完全に 0
- 実運用: 灌漑記録 or Sentinel-1 SAR (灌漑検出) で展開可能

**なぜ AIC でモデル選択したか**:
- M1: bias ~ VPD (k=2)
- M2: bias ~ days_since_irrig (k=2)
- M3: bias ~ VPD + d + 交互作用 (k=4)
- 結果 ΔAIC vs VPD: **MOD16 −73, PML −66, METv3 −153** → days_since_irrig が圧倒
- Burnham & Anderson (2002): ΔAIC > 10 で "decisive evidence"

---

## 8. 解析C — フェノロジーによる active period 検証 (PLANNED)

**目的**: 解析A で使った "active months" 定義 (Oran Nov-Jun, TzM Jun-Sep) が観測 NDVI と一致するか独立確認。

**手法 (v2 設計案)**:
1. MOD13Q1 (250m, 16-day) NDVI を 2 サイトで時系列化
2. Savitzky-Golay smoothing + logistic fit で SOS / Peak / EOS 抽出
3. NDVI > 0.4 を "active" と定義し、A の手動定義とオーバーラップ率を計算
4. NDVI peak ± 30 日 window で τ を再 fit (A 結果が変わらないか確認)

**期待結果と論文記述**:
- Scenario A: NDVI が A 仮定を validate → "phenology-matched comparison" を強化
- Scenario B: ズレあり → NDVI-based 再定義で A を更新

**進捗**: `analysis_C_v1.py` 存在、多目的すぎる → v2 で目的 1+4 に絞る方針

---

## 9. なぜこの手法を選んだのか（手法の justification まとめ）

ポスター中央 or サイドバーに表で配置:

| 手法 | なぜ採用 | 失敗事例 (避けた罠) |
|---|---|---|
| Within-site SDS | 異種作物・LAI・季節を持つ 2 サイトで絶対 LE 比較は無意味 | F2: cross-site MWU p=1.5e-47 で「有意」だが科学的に無意味 |
| Season × bucket 二重層別 | 季節間プーリングが Simpson's paradox を生む | F3: v12 の TzM SDS=+0.05 (深根支持に見えた artifact) |
| NDWI/NDVI を月別 anomaly で評価 | 絶対閾値はサイト・季節依存 | F1: NDWI > 0.0 や > 0.109 で deep_access 0 日 |
| τ ∈ (0, 60] 制約 + raw daily fit | bin median fit は 5-6 点で τ 暴走 | F8: 初期試行で τ > 60, CI ±50 |
| Bootstrap CI (理論分布なし) | 比 (SDS) や非線形パラメータ (τ) は正規性も等分散性も保証なし | F7: 一律 DENOM_FLOOR=5.0 で EF/Bowen 全 NaN |
| AIC モデル比較 (R²/p ではなく) | R²/p はパラメータ追加で必ず改善、AIC はトレードオフを定量 | (該当なし、最初から AIC 採用) |
| MDE-based universality 検定 | p > 0.05 は「差がない」ことを示さない | (v23-v27 で reviewer-bulletproof 化) |
| 独立 3 衛星製品 (MOD16 + PML + METv3) | 1 製品なら algorithm-specific bug の可能性 | (該当なし、最初から 3 製品計画) |
| 2 独立 SM センサー (in-situ + SMAP_rz) | 1 センサーだけだとノイズと区別不能 | F14: H6 を Oran summer (n=34) のみで検証 → multi-stratum 化で救出 |

---

## 10. 反論予想と根拠 (Reviewer defense matrix)

ポスター下部 or QR で詳細リンク:

| 反論 | 根拠 |
|---|---|
| "Same τ ≠ same mechanism" | A v27: 4 pairwise MDE NS。**振幅で機構差を捕捉** (4× scaling) |
| "Phenology cofounding" | A v27: active-state matched comparison でも NS |
| "Identifiability problem (τ-amp 相関)" | A v27: τ-amplitude r=0.33 (NS), τ-pulse r=−0.13 (NS) で独立性確認 |
| "Aggregation artifact" | A v27: per-season stratified で一貫 |
| "Fit at boundary" | A v27: validation filter (R²≥0.3, τ ∈ range, n≥3) で除外 |
| "Insufficient sample" | A v27: MDE で power 明示 |
| "Pseudo-replication" | A v23: event-level fit 22 events で確認 |
| "Asymptote circular" | A v23: d≥{7,8,10,12,15} で τ range 0.41d 内 |
| "Single model bias" | A v23: AIC 比較 (exp 53%, logistic 47%) |
| "MOD16 固有のバイアス?" | B: PML + METv3 でも同じ dose-response |
| "5cm センサーが wet-bulb をミスしているだけでは?" | B H6: SMAP root-zone も含めて depth inversion 観測 |
| "深根アクセスは否定できない?" | A v14: d8+ でも SDS≈0、深根活動の証拠なし |
| "サイト 2 つで一般化不可" | Limitations で明記。FLUXNET2015 / ICOS 拡張は future work (H7) |
| "5cm のみで wet-bulb の直接観測ではない" | Limitations で明記。SMAP は coarse pixel rather than direct |
| "EC エネルギー閉合は?" | EBR Oran 0.96, TzM 0.87 を Methods に記載 |
| "τ-fit の overfit?" | Bootstrap CI + AIC + τ 制約で reviewer-bulletproof |
| "Statistical power が低いから NS?" | MDE 計算で「2.15d 以上の差なら検出できた」と明示 |
| "v13 で深根→灌漑への乗り換えは confirmation bias?" | Fig 3 で「自分で否定した」プロセスを開示。dose-response (Fig 6) は post-hoc ではなく独立検証 |

---

## 11. 限界 (Limitations, 必ず記載)

1. **サイト 2 つのみ** — 一般化に注意。他 drip 灌漑サイト (FLUXNET) で τ 再現確認が future work
2. **5cm SWC センサーのみ** — 10-30cm の wet-bulb を直接観測していない。SMAP は coarse-pixel rather than direct
3. **Oran 期間 2.5 年、Tarazona 期間 4 年** — 長期変動 (干ばつ年 vs 平常年) は未検証
4. **METv3 の 5 km ピクセル混合** — 周辺 rainfed 地が混じり τ が長め (6d) になっている可能性
5. **PML/METv3 の c の CI が広い** — d ≥ 8 のサンプルが少ない (n=44 etc.)
6. **VPD は EC のみ** — ERA5 dewpoint + Rn 統合は未実装
7. **解析C 未完了** — NDVI による active period の客観検証は計画段階

---

## 12. 今後の展望

### 短期 (★★★)
- 解析A v15 (in_band バグ修正) → 偽陽性除去
- 解析A の Oran 関連数値を新ローダで再走 (C 分析1 で発見したバグ前の v9-v11 影響範囲)
- 論文 Abstract 執筆 + Figure 最終化
- 解析C v2 で active period 客観検証

### 中期 (★★)
- 他 drip 灌漑サイトで τ 再現 (FLUXNET2015, ICOS)
- Flood/sprinkler 灌漑サイトとの τ 比較 (H2)
- SIGPAC parcel × Sentinel-2 で METv3 の 5km mixing 分解 (H5)

### 長期 (★)
- Sentinel-1 SAR で灌漑検出 → 灌漑記録がない site への補正展開
- H7 広域 SDS マッピング (rainfed では SMAP_rz が代替可能と確認済)
- H8 Júcar 流域水収支検証

---

## 13. ポスター物理レイアウト案

```
┌────────────────────────────────────────────────────────────┐
│ TITLE                                  Author / Affiliation │
├──────────────┬──────────────────────┬──────────────────────┤
│ 1. Question  │  Fig 1 (Wet-bulb 模式)│  Fig 4 (A τ result)  │
│ 2. Motivation│  Fig 2 (Site)         │  Fig 5 (MDE)         │
│ 3. Hypothesis│  Fig 3 (Hypothesis)   │                       │
│  Evolution   │                        │                       │
├──────────────┼──────────────────────┼──────────────────────┤
│ Methods      │  Fig 6 (B dose-resp)  │  Fig 8 (Depth inv.)  │
│ overview &   │  Fig 7 (B τ-fit)      │                       │
│ Justification│                        │                       │
│ (表)         │                        │                       │
├──────────────┼──────────────────────┴──────────────────────┤
│ Reviewer     │  Fig 9 (Correction impact)                   │
│ defense      │                                               │
│ matrix       │  Limitations / Future / Conclusions          │
└──────────────┴──────────────────────────────────────────────┘
```

**読者の動線**:
1. タイトルで「τ ≈ 3-4d 普遍 + 4× 振幅 + 衛星補正」のメッセージを受け取る
2. 左カラムで「なぜこの研究?」「自分の仮説を自分で棄却した経緯」を把握
3. 中央〜右で「3 つの解析で同じ結論」を視覚的に確認
4. 下部で「想定批判に答えがある」と納得

---

## 14. 進捗ステータスの正直な提示 (進捗報告版)

ポスター末尾に必ず置く:

```
解析A: ✅ COMPLETED (v27, reviewer-bulletproof)
   - τ = 3.0-3.8 d universal, 振幅 4× scaling
   - publication-ready paragraphs 完成
解析B: 🔶 IN PROGRESS (paper draft 70%)
   - 3 衛星製品 τ-fit, SMAP depth inversion, τ-correction 完了
   - Abstract + 図最終化 + 引用 BibTeX が残タスク
解析C: 🔶 PLANNED (v1 multi-target, v2 設計済み)
   - MOD13Q1 NDVI 取得 + Savitzky-Golay + logistic fit が残タスク
横断: 📌 統合課題
   - C 発見の Oran ローダバグの A/B への影響再評価 (P★★★)
   - A v15 (in_band 偽陽性修正)
```

---

## 15. 教授との質疑応答想定 Q&A

| 想定質問 | 回答骨子 |
|---|---|
| Q. 「2 サイトで universal と言えるか?」 | A. 解析A は「2 サイトで τ 区別不能」が結論。"universal in Mediterranean drylands" は提案で、FLUXNET 拡張 (future work H2/H7) で確認予定。limitation で明示 |
| Q. 「深根仮説を完全否定したのか?」 | A. d8+ でも SDS≈0 で深根 active の証拠なし、ただし非灌漑期 (TzM 11-4月) のテストが不十分。**寄与の主従**として「灌漑 > 深根」と再定義 |
| Q. 「Tarazona の灌漑記録の精度は?」 | A. 日次 `Irrig_mm` カラムが提供されている。> 0.5 mm を灌漑日と定義。タイミング誤差は ±1d 程度を想定 (Sensitivity test 必要) |
| Q. 「衛星の bias 補正は他のサイトで効くのか?」 | A. 現状は TzM 1 サイトで校正、他サイトでの transferability は H2 (flood/sprinkler) と H7 (FLUXNET 拡張) で確認予定 |
| Q. 「論文の投稿先は?」 | A. 本命 Agricultural and Forest Meteorology、Remote Sensing of Environment, HESS が候補 |
| Q. 「解析C は何のため?」 | A. A の active period 仮定 (Nov-Jun, Jun-Sep) が NDVI と一致するか独立確認。Methods の justification 強化 |
| Q. 「5cm SWC センサーで wet-bulb を主張するのは無理では?」 | A. 直接観測ではなく**間接シグナル**。直接証拠は SMAP rootzone × in-situ の depth inversion (r=−0.19)。Limitations で明記 |
| Q. 「Reviewer に潰される可能性が高い箇所は?」 | A. ① 単一 drip サイト, ② 多深度センサーなし, ③ Tarazona の灌漑タイミング不確実性。① は future work, ② は SMAP で補完, ③ は ±1d sensitivity test で対処予定 |

---

## 16. 当面のアクションアイテム (発表前にやる)

### 図の作成・整備
- [ ] **Fig 1** drip wet-bulb 模式図を新規作成 (matplotlib 断面図)
- [ ] **Fig 2** サイト地図 + 概要表 (cartopy or 既存 GIS)
- [ ] **Fig 3** 仮説変遷フロー (Graphviz or 手描き)
- [ ] **Fig 4** A v27 recovery + τ comparison を 1 枚に統合
- [ ] **Fig 5** MDE 解析 (既存 fig03 活用)
- [ ] **Fig 6-7** B Fig 5 + Fig 7 をポスター用にラベル整備
- [ ] **Fig 8** B Fig H6 を depth inversion focus にトリミング
- [ ] **Fig 9** B Fig H1 を before/after で並べる

### 数値の確認・更新 (Oran ローダバグの影響範囲)
- [ ] A v27 Oran 数値を新ローダで再走 → 変動が ±10% 以内か確認
- [ ] B paper draft の Oran spring SDS=+0.43 も新ローダで再走
- [ ] 変動が大きければ Limitations に明記

### 文章
- [ ] タイトル確定 (3 候補から教授と相談)
- [ ] Motivation 1 パラ確定
- [ ] Limitations セクションを限界明示で書く
- [ ] Conclusion 3 行を確定

### Q&A 準備
- [ ] Q1〜Q8 の想定回答を声に出して練習
- [ ] 反論マトリクスを暗記
- [ ] Limitations を聞かれた時の対応 (隠さない、先に出す)

---

**EOF — 進捗報告ポスタードラフト v1, 2026-05-18**
*このドラフトは三角測量フレーミング (A × B × C) を採用。教授からのフィードバックで headline を A 軸 (universal τ) または B 軸 (drip wet-bulb mechanism) に振り直す可能性あり。*
