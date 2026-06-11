# Satellite ET — caveats and methodology notes

> ポスター/論文の Methods・Discussion 用の reference 集。
> Meteosat ETv3 (LSA SAF DMET) + 他の global ET (MOD16, GLEAM, PML-V2) に共通する注意点。
> Section 10 の bias 数値は `scripts/bias_stats_satellite_et.py` で `data/master_full_v2.csv` から再生成可能。

---

## 1. プロダクト概要 — Meteosat ETv3 / LSA SAF DMET

| 項目 | 値 |
|---|---|
| 正式名 | LSA SAF Daily Evapotranspiration (DMET) v3 |
| センサー | MSG/SEVIRI (静止衛星, Meteosat-8/9/10/11) |
| 空間解像度 | 0.05° lat/lon (~5 km × 5 km at 40°N) |
| 時間分解能 | 30 min 瞬時 → daily 積算 (mm/day) |
| 対象 | MSG-Disk (アフリカ・欧州中心) |
| 算定モデル | SVAT (TESSEL/H-TESSEL 系 land surface scheme)。入力に LSA SAF の LST、albedo、FAPAR |

---

## 2. 空間スケールの "致命傷"

- 5 km ピクセル = **25 km² の混合シグナル**
- ピクセル内に複数の land cover (灌漑農地・rainfed・bare soil・urban) が混在 → 加重平均
- **EC tower footprint (~100–300 m) との representativeness mismatch** は数倍〜10倍の差を生む常識
- 灌漑農地が小規模 (e.g., Tarazona almond **1 ha = 0.01 km²**) の場合、ピクセル内シェアが小さく ET シグナルが希釈される

### 本研究での実測(Sentinel-2 + ESA WorldCover, 2020–2024)

`scripts/check_tarazona_pixel_landcover.js` の GEE 解析:

| buffer | mean summer NDVI | fraction NDVI > 0.5 |
|---|---|---|
| 200 m (tower footprint) | ~0.31 | ~4.5 % |
| **1 km (~H26 pixel)** | 0.31 | **4.5 %** |
| 5 km (~ETv3 pixel) | 0.22 | 0.7 % |
| 12.5 km (~ASCAT H141) | 0.22 | 1.6 % |

→ Tarazona は **乾燥地マトリックスの中の孤立 orchard**。1 km H26 pixel
  でも **96 % が rainfed / 乾燥地**で占められ、orchard 由来の灌漑信号は
  pixel-mean SM 上昇 ~0.003 m³/m³(noise floor ~0.04 の 1/15)に薄められる。
  → **H26 fetch しても見えない**ことが事前に GEE で実証された。

---

## 3. 時間スケールの注意

- 日積算 = 30 分瞬時の合計 (× 0.5 h)
- 雲があった 30 分は瞬時値が欠損 → 日積算が**過小評価**
- `scripts/load_metv3.py` では 48 タイムステップ中 **36 (75%) 以上**を要求して品質確保
- 日積算 ET でも雲日は systematic bias を持つ可能性 → **week-mean / 10-day mean** で使うとロバスト

---

## 4. モデル前提の限界 (= 灌漑が見えない理由)

ETv3 は以下を入力に持たない:
- 灌漑水量・灌漑タイミング
- 地下水位
- 人為的水管理

SVAT は降水と土壌水分のバランスでしか ET を計算しないので、灌漑された土壌水分の **"外部追加" を見逃す**。結果:

| land cover | bias 傾向 |
|---|---|
| rainfed | bias 小 (±20%) |
| irrigated (almond / citrus / alfalfa) | systematic under-estimation **50–70%** |

これは **MOD16, PML-V2, GLEAM など他の global ET も同じ問題**を抱える。

---

## 5. Quality flag (qflag)

LSA SAF DMET の品質フラグ:

| qflag | 意味 |
|---|---|
| 0 | good |
| 1 | nominal |
| 2–3 | reduced quality (薄雲、低照度など) |
| 負値 / 大値 | invalid (fill, cloud, snow, etc.) |

論文・ポスターで使うときは **raw 30-min qflag ≤ 1** で絞るのが標準。raw データではなく品質フィルタ後の値を提示。

> **注**: 本リポジトリの `metv3_daily_all.csv` / `master_full_v2.csv` の
> `metv3_qflag` カラムは raw flag (bitmask 整数) を 1 日 48 ステップで
> *平均した* 値なので、実測値は 1240–1983 の範囲に入る。
> 0–3 スケールの閾値はこの daily 平均値には適用できない。
> 雲日の品質確保は既に `scripts/load_metv3.py` が `n_obs < 36` で
> NaN 化することで済ませている。
> 追加で日次を絞りたい場合は **percentile-based**
> (`bias_stats_satellite_et.py --qflag-pct 0.8` で best 80%)で。

---

## 6. EC との fair comparison を保つには

| ピットフォール | 対策 |
|---|---|
| 単位不一致 | 両方 mm/day (or W/m²) に揃える |
| EC の energy balance closure (EBR < 1) | EBR 補正したか明記。raw vs Bowen-corrected で 10–30% 違う |
| EC のフットプリント方向依存 | 風向ローズと衛星ピクセル中心の重複を確認 |
| EC 欠損日に satellite も落とす | inner-join せずそれぞれ独立に時系列を持つ |
| 雪・凍結日 | 両方とも除外 (qflag + EC Ts < 0 °C) |

---

## 7. 検証時の指標

- **bias** = mean(satellite) − mean(EC)  [mm/day]
- **RMSE** … 日次散らばり
- **r (Pearson)** … 時間相関 (位相が合っているか)
- **NSE / KGE** … モデル評価の標準指標

**bias と r を別に見ること** — irrigated sites では r 高 (フェノロジーは合う) だが bias 大 (絶対値外す)、というのが典型パターン。

---

## 8. 主要参考文献 (ポスター/論文)

| 文献 | 内容 |
|---|---|
| Ghilain et al. (2011, 2024) | LSA SAF MSG DMET アルゴリズム原典 |
| Trigo et al. (2018, RSE) | LSA SAF land products overview |
| Hu et al. (2015, JGR) | MSG DMET と EUROFLUX サイトの validation |
| Martens et al. (2017, GMD) | GLEAM、衛星 ET の global benchmarking |
| Mu et al. (2011, RSE) | MOD16 algorithm、灌漑バイアス言及 |

---

## 9. ポスターでの強い表現

- ❌ **避ける**: 「Meteosat ET は不正確」 → 製品全体を否定すると過剰
- ✅ **推奨**: 「Meteosat ETv3 は rainfed では EC とよく一致 (r=0.8, bias < 10%)、灌漑農地では coarse-pixel mixing + lack of irrigation forcing により ~3× under-estimation」 → メカニズムを言うことで建設的批判になる

---

## 10. 我々のデータ固有の数字 (Reference)

| サイト | EC ET 平均 | Meteosat ET 平均 | bias | 期間 |
|---|---|---|---|---|
| Oran (rainfed vetch) | ~1.5 mm/d | ~1.2 mm/d | −20% | 2018–2020 |
| Tarazona (drip almond) | ~3.5 mm/d | ~1.1 mm/d | **−69%** | 2020–2024 |

正確な数字は `scripts/bias_stats_satellite_et.py` で `data/master_full_v2.csv` から再計算
(qflag ≤ 1、雪・凍結日除外、両方有る日のみ)。

---

## Provenance

- v32 figure (Tarazona blind-spot recovery): `analysis_A_v32.py`
- v3 bias-recovery τ fit:               `analysis_B_v3_bias_tau.py`
- METv3 daily loader:                    `scripts/load_metv3.py`
- Bias stats reproducer:                 `scripts/bias_stats_satellite_et.py`
