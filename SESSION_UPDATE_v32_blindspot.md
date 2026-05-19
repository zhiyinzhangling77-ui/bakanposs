# Session Update — v32 satellite blind spot + bias decomposition

> 続編: `SESSION_UPDATE_v28_v32.md` の後、Tarazona 衛星 ET bias の
> mechanism を深掘りした議論ログ。**次セッションで必読**。
>
> Last updated: 2026-05-19
> Branch: `claude/quantify-water-divergence-LxPGr`

---

## 0. このセッションで起きた重要な転換 (TL;DR)

1. **τ_bias は独立な測定ではない** — 数学的に EC の exp fit と等価。consistency check として扱う。main argument は `amp_EC = 95 vs amp_Sat ≈ 0`。
2. **ETv3 は衛星 SM (H-SAF) を入力に持つ** — 「降水だけで ET 計算」は誤り。drip blind spot は SM 入力チェーンの 3 重の解像度限界による。
3. **SMAP は 9 km なので灌漑検証には不適** — orchard 1 ha との比較で循環論証になる。
4. **Tarazona orchard サイズ訂正** — 1 ha (0.01 km²)、過去メモの「0.1 km²」は誤り。
5. **narrative を "blame satellite" から "attribute bias to irrigation" に転換** — 仮説検証型の方が強い。

---

## 1. v32 evolution(figure 改訂の履歴)

### v32.0 (initial)
3-panel: time series (a) + bias recovery (b) + τ comparison (c)。

### v32.1 (bug fixes)
- METv3 fit が R²≈0 で τ 床に張り付く → panel (c) で hatched + "FIT INVALID" 表示
- errorbar の負 yerr clip
- panel top alignment を `constrained_layout=True` に切替

### v32.2 (final, 現状)
**panel (a) を削除** → 2-panel side-by-side(bias recovery + τ comparison)。理由: parquet engine 未インストール + tara CSV fallback も EC 列発見失敗で空白になっていた。

ファイル: `analysis_A_v32.py`、`output_analysis_A_v32/fig04b_tarazona_blindspot_v32.{png,pdf}`

### v32.3 (✅ 実装済)
**τ_bias を main argument から外し、management amplitude 比較に置換**:
- panel (a) bias recovery: 変更なし(bias 内に pulse が保存されることを可視化)
- panel (b) `draw_panel_c` (τ bars) → `draw_panel_amp` (amplitude bars)
  - 3 bars: amp_EC = 95 / amp_Sat ≈ 0 / amp_bias = 95 W/m²
  - "Sat captures 0%", "Bias inherits 100%" の annotation
- suptitle・caption も amplitude-first に書き換え
- τ_bias は caption の (b) 節に **"consistency check, not independent evidence"** と注記
- CSV 出力名 `v32_panel_b_tau_bars.csv` → `v32_panel_b_amp_bars.csv`

---

## 2. EC vs Meteosat ETv3 bias 統計 (確定数値)

`scripts/bias_stats_satellite_et.py` で `data/master_full_v2.csv` から算出。
qflag フィルタは無し(daily mean qflag は bitmask 平均で 1240–1983、§5 参照)。
Freeze-day (Ta_min ≤ 0 °C) 除外、paired daily only。

| サイト | n paired | mean EC | mean Sat | bias | bias_pct | RMSE | r | KGE |
|---|---|---|---|---|---|---|---|---|
| **Oran (rainfed)** | 605 | 1.17 | 1.18 | +0.01 | **+1 %** | 0.58 | **0.82** | **+0.61** |
| **Tarazona (drip)** | 699 | 3.37 | 1.03 | −2.34 | **−70 %** | 3.08 | **0.07** | **−0.21** |

→ **rainfed では excellent agreement / drip では r も bias も KGE も全壊**。位相すら合っていない。

---

## 3. τ_bias 解析の位置づけ(重要な再整理)

### Math

実測で `LE_Sat(d) ≈ const`(Sat 灌漑後ほぼフラット)。よって:

```
bias(d) = LE_EC(d) − LE_Sat(d)
       ≈ LE_∞ + (LE_0 − LE_∞)·exp(−d/τ_EC) − const
       = (const′) + EC_amplitude·exp(−d/τ_EC)
```

→ **bias の exp fit は EC fit と数学的に等価**(定数項だけ違う)。
   τ_bias = 4.57 d ≈ τ_EC = 3.36 d、amp_bias = 95 ≈ EC amp = 95 W/m² は **construct のためそうなる**。

### 結論

- **❌**: 「衛星 bias も解析A の τ で減衰する」を独立証拠として扱う
- **✅**: 「Sat は本当に flat、解析A amplitude が丸ごと bias に流れ込んでいる」の **consistency check** として使う

### Main claim を切り替える

| 旧 main claim | 新 main claim |
|---|---|
| τ_bias = 4.57 d → 解析A timescale を再現 | **amp_EC = 95 W/m², amp_Sat ≈ 0 → 衛星は管理 amplitude をゼロ検出** |
| τ_bias 図がメイン | **amp_EC vs amp_Sat の bar が直接の主張** |
| τ_bias は appendix で「math 整合性」だけ示す |

→ v33 (or v32.3) で実装予定。

---

## 4. 6-panel overview figure(`scripts/poster_overview_figure.py`)— narrative の核

3 行 × 2 列(NDVI / GPP / ET × Oran / Tarazona)。

| 量 | Oran | Tarazona |
|---|---|---|
| NDVI (Sentinel-2 vs EC tower) | 良く一致 | 良く一致 |
| GPP (Meteosat MGPP vs EC tower) | 概ね一致 | 概ね一致(MGPP やや低めだが phase OK) |
| **ET (Meteosat ETv3 vs EC tower)** | **一致** | **乖離: amplitude −70%、phase ズレ** |

### Tarazona ET の phase ズレ(重要)

- **Sat ETv3**: 春 (Apr–Jun) に小ピーク、夏は flat
  → NDVI 駆動の potential ET。canopy cycle を追っている
- **EC tower**: 夏 (Jul–Sep) に急峻ピーク(4–7 mm/d)
  → 実際の water use。irrigation cycle を追っている

→ **drip 灌漑は "greenness と water use の decoupling" を生む。衛星 SVAT はこれを分離できない**。

### Narrative bridge(ポスター Fig 1 直下に置く文、英)

```
At drip-irrigated Tarazona, Meteosat ETv3
  •  under-estimates the tower by ~70 % in magnitude, AND
  •  peaks in spring (NDVI cycle) instead of summer (irrigation cycle).
The satellite is tracking greenness; the tower is tracking water.

Question — is this residual a random mismatch, or the irrigation
pulse itself? If the bias carries the Analysis-A τ ≈ 3 d signature,
it IS the missing irrigation response.  → Section 4c
```

---

## 5. ETv3 mechanism — 正しい説明

### 前回の誤り

「ETv3 は降水と土壌水分のバランスでしか ET を計算しない、灌漑を入力に持たない」と書いたが、これは**間違い**。

### 正しい説明

ETv3 の入力には **衛星土壌水分** が含まれる:
- **H141**: ASCAT 表層 SSM (~12.5 km)
- **H142**: ASCAT 根圏 SWI (~12.5 km, derived)
- **H26**: SCATSAR-SWI (1 km、ASCAT + Sentinel-1 融合)

SVAT で `ET = ET_pot × f_canopy(LAI/FVC) × f_water(SM/SM_capacity)` を計算。

### でも drip が見えない 3 つの理由

| # | 限界 | Tier | 根拠 |
|---|---|---|---|
| 1 | **空間スケール**: H141/142 = 12.5 km、H26 = 1 km vs orchard 1 ha (0.01 km²) → 希釈 100×–156,000× | **A** (documented) | Wagner et al. 2013; Bauer-Marschallinger et al. 2018 |
| 2 | **マイクロ波物理**: 表層 ~5 cm 感度、Mediterranean 夏で数時間で乾燥、成熟 canopy の VOD 減衰 | **A** (textbook) | Ulaby et al.; Vreugdenhil et al. 2020 |
| 3 | **SVAT `f(SM)` 構造**: generic beta function、orchard root system の bi-modal SM に非適合 | **C** (inference) | Ghilain et al. 2011; Balsamo et al. 2009 |

加えて **灌漑量・タイミングそのものは入力に含まれない**(management data なし)。

`SATELLITE_ET_NOTES.md` §4 もこの内容に更新済み。

---

## 6. 「衛星 SM をテストすべきか?」議論の結論

### 候補と評価

| product | 解像度 | Tarazona (1 ha) 比 | 在庫 | 評価 |
|---|---|---|---|---|
| SMAP L4 | 9 km | **90,000×** | ✅ `smap_daily.csv` | **使うと循環論証** — coarse pixel で見えなかったから coarse pixel がダメ、と言うことになる |
| SMAP L3 | 36 km | 360,000× | ❌ | 論外 |
| ASCAT H141 | 12.5 km | 156,000× | ❌ | 同上、coarse すぎ |
| ASCAT H142 | 12.5 km | 同上 | ❌ | 同上 |
| **H26 (SCATSAR-SWI)** | **1 km** | **100×** | ❌ | **唯一可能性あり**、ただし以下条件付き |
| **Tower SWC (in-situ)** | 点 | ground truth | ✅ `master_full_v2.csv` の `SWC` 列 | **これが必須の reference** |

### H26 の信号検出可能性試算

灌漑前後の orchard 内 ΔSWC ≈ 0.20 m³/m³。1 H26 pixel に orchard が 1% 占めるとき:
```
pixel-mean ΔSM = 0.01 × 0.20 = 0.002 m³/m³
```
これは典型 SM noise floor (~0.04 m³/m³) より 1 桁低い → **単独 orchard では検出不可**。

ただし pixel 内に複数の灌漑農地がある場合は累積で立ち上がる可能性:
- 20% 灌漑 → 0.04 m³/m³(noise と同程度、辛うじて検出)
- 50% 灌漑 → 0.10 m³/m³(検出可能)

→ **H26 fetch する価値は Tarazona 周辺(±500 m〜±2 km)の土地利用次第**。

### 判断フロー(✅ 実施済、結果は candidate C)

1. **GEE Sentinel-2 + ESA WorldCover 解析** `scripts/check_tarazona_pixel_landcover.js`、実施済
2. 結果 (summer 2020–2024 mean):

   | buffer | mean NDVI | fraction NDVI > 0.5 |
   |---|---|---|
   | 200 m (tower) | 0.31 | 4.5 % |
   | **1 km (~H26)** | 0.31 | **4.5 %** |
   | 5 km (~ETv3) | 0.22 | 0.7 % |
   | 12.5 km (~ASCAT) | 0.22 | 1.6 % |

3. **判定**: 1 km buffer の active fraction = 4.5 % は判定閾値 5 % をわずかに下回る。
   pixel-mean ΔSM = 0.013 × 0.2 ≈ **0.003 m³/m³** で noise floor の 1/15。
   → **H26 fetch しても検出不能**(事前に経験的に確定)
   → **candidate C 採用**: 衛星 SM 解析は外し、`amp_EC vs amp_Sat` 中心の
     argument + tower SWC で論を閉じる。

### Bonus: 空間希釈論の証拠 tier 格上げ

これまで「空間希釈で見えない」は **Tier A**(general 文献からの推論)だったが、
本研究 Tarazona pixel で **実測 4.5 % active fraction** が示せたことで、
**Tier A + 本研究での Tarazona 直接実測** に格上げ。
reviewer に「あなたのサイトで本当に空間希釈が効いているのか?」と聞かれたら即答可能。

---

## 7. Poster template の現状と次の改訂

### 現状(`poster/build_poster_template.py`、`poster_template_A0.pptx`)
- A0 縦、2 列レイアウト
- title strip → Intro+Methods (左) + Results A+B (右) → Take-home strip
- すべて英語、本文 30 pt 以上

### 未反映の改訂(v33 で入れる予定)
1. **Title**: narrative 反映版「Quantifying the irrigation contribution to satellite–tower ET bias in Mediterranean drylands — a 3-day-clock decomposition」
2. **Fig 1 (6-panel overview)** placeholder を Result B 冒頭に追加
3. **Fig 1 直下の narrative bridge**(§4 参照)
4. **"Why exp fit?" methods box**(consistency check として位置付け)
5. **"Why drip is invisible to ETv3" mechanism box**(3 つの解像度限界、§5 参照)
6. **Result B の Key numbers** を amp 比較中心に書き換え

---

## 8. やってはいけないこと(失敗回避)

1. ❌ **τ_bias を独立な物理測定として扱う** — math 上 EC fit と等価、consistency check のみ
2. ❌ **SMAP 9 km で空間希釈論を実証しようとする** — 循環論証
3. ❌ **「ETv3 は降水だけで ET を計算」と書く** — SM は入力に入っている
4. ❌ **Tarazona orchard を「0.1 km²」と書く** — 正しくは 1 ha = 0.01 km²
5. ❌ **理由 3 (SVAT 構造) を観測実証扱いする** — SVAT を独自に走らせていないので mechanistic background に留める
6. ❌ **「Meteosat ET は不正確」と総括する** — rainfed では r=0.82 で良好、製品自体は機能する。drip 特異の blind spot

---

## 9. 確定数値(コピペ用、updated)

```
─────────────────────────────────────────
解析A v31 main result(変更なし)
─────────────────────────────────────────
Tarazona active (drip-irrigated almond, Jun-Sep):
  n events = 41
  τ = 3.36 d (95% CI: 2.44–4.90), SE = 0.62
  LE_0 ≈ 210 W/m², LE_∞ = 114 W/m²
  Amplitude = 95 W/m² (≈ 3.4 mm/day)

Oran active (rainfed cereal, Nov-Jun):
  n events = 10
  τ = 2.82 d (95% CI: 1.83–5.48), SE = 0.90
  LE_0 ≈ 36 W/m², LE_∞ = 17 W/m²
  Amplitude = 21 W/m² (≈ 0.7 mm/day)

|Δτ| = 0.54 d < MDE = 2.15 d (α=0.05) → NS
Amplitude ratio = 4.5×

─────────────────────────────────────────
解析B v3 bias-recovery fit(consistency check として位置づけ)
─────────────────────────────────────────
Tarazona bias-pool (n events = 48):
  τ_bias = 4.57 d (95% CI: 3.03–13.10), SE = 2.41
  amp_bias = 95 W/m² (= EC amplitude)
  R² = 0.86
  → Sat fit は τ 床で R² ≈ 0 (FIT INVALID)
  → bias が EC の管理 amplitude をそのまま継承

─────────────────────────────────────────
解析B v5 EC vs ETv3 bulk stats(qflag 制約なし、freeze day 除外)
─────────────────────────────────────────
Oran (n=605):       r = 0.82,  bias = +1 %,   KGE = +0.61
Tarazona (n=699):   r = 0.07,  bias = −70 %,  KGE = −0.21

─────────────────────────────────────────
サイト固定情報(訂正版)
─────────────────────────────────────────
Tarazona orchard size = 1 ha = 0.01 km² (※過去メモの 0.1 km² は誤り)
H-SAF pixel sizes:
  H141/H142 = 12.5 km;  H26 = 1 km
SMAP L4 = 9 km
─────────────────────────────────────────
```

---

## 10. 次セッションの最初に決めること(更新版)

✅ (A) と (C) は完了:
- (A) GEE で Tarazona pixel 周辺確認 → **active fraction 4.5 %、candidate C 確定**
- (C) v32 を amp comparison 中心に refactor 済(panel (b) τ bars → amplitude bars)

残作業:

| 選択肢 | 内容 | 工数 |
|---|---|---|
| (B) **Tower SWC event-recovery 解析** — `master_full_v2.csv` の `SWC` を v32 と同枠組みで | candidate C の証拠補強(in-situ が irrigation を確認) | ~2 時間 |
| (D) **poster template v2** — v32 refactor + Fig 1 (6-panel) + GEE 周辺地図 + narrative | ポスター ready | ~30 分 |
| (E) **論文 Methods 段落原稿** — SATELLITE_ET_NOTES.md と本ファイルから組み立て | 論文化準備 | ~30 分 |

推奨順: **(B) → (D) → (E)**。

---

## 11. 重要ファイル一覧(本セッション分)

### コード
- `analysis_A_v32.py` — 2-panel poster figure (bias recovery + τ comparison)
- `scripts/bias_stats_satellite_et.py` — EC vs ETv3 bulk stats reproducer
- `scripts/poster_overview_figure.py` — 6-panel NDVI/GPP/ET overview(既存、本セッションで再評価)
- `poster/build_poster_template.py` — A0 縦 pptx テンプレ生成器

### ドキュメント
- `SATELLITE_ET_NOTES.md` — 衛星 ET caveats memo(§4 を SM 入力ありに訂正済)
- `SESSION_UPDATE_v28_v32.md` — 前セッション(v28-v32 initial)
- **`SESSION_UPDATE_v32_blindspot.md` — 本ファイル**

### Outputs(再生成可能)
- `output_analysis_A_v32/fig04b_tarazona_blindspot_v32.{png,pdf}`
- `output_bias_stats/bias_stats_summary.csv`
- `output_bias_stats/fig_scatter_ec_vs_metv3.png`
- `poster/poster_template_A0.pptx`

---

## 12. 次セッション開始時の必読

1. `CLAUDE.md`(navigation hub)
2. `SESSION_UPDATE_v28_v32.md`(v28-v32 までの履歴)
3. **本ファイル `SESSION_UPDATE_v32_blindspot.md`** ← 最新議論
4. `SATELLITE_ET_NOTES.md`(衛星 ET caveats、§4 訂正済)

これら 4 つを読めば本セッションの全 context が再現される。
