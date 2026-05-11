# Session Summary — Deep Root vs Irrigation in Mediterranean Almond/Cereal EC Sites

Last updated: 2026-05-11
Branch: `claude/quantify-water-divergence-LxPGr`
Repo: `zhiyinzhangling77-ui/bakanposs`

---

## 1. 研究背景・目的

### 当初の研究目的
Mediterranean の 2 つの Eddy Covariance(EC)サイトで観測された **「表層 SWC と蒸散(LE)の解離(decoupling)」**現象を定量化し、その原因を特定する。

- **Oran**(スペイン Albacete 近郊、雨養 winter cereal): 表層 SWC が乾くと LE が顕著に落ちる
- **Tarazona**(同地域、灌漑 almond): 表層 SWC が乾いても LE を維持する

### 当初仮説
「Tarazona アーモンドは **生物学的深根** を持ち、表層センサー(5cm 深)では捉えられない地下水・深層浸透水を利用して蒸散を維持している」

### なぜこの分析を行ったか
- FLUX 時系列の日変化(diurnal cycle)解析で、Tarazona が VPD 高条件で LE を増やす(典型的木本性反応)現象が確認されていた
- 衛星 SWC product(SMAP 等)から ET を推定するモデルが、灌漑農地で系統誤差を出す問題への科学的貢献の可能性
- 「深根」を表層 EC 観測のみで定量化できるかは方法論的に興味深い

### 前提条件
- Oran: 雨養 cereal、生育期 11–6月、SWC 単位は %
- Tarazona: 灌漑 almond、生育期 1–10月、SWC 単位は m³/m³(→ ×100 で% に変換)
- ERA5 で VPD を統一(EC タワー観測の VPD はサイト間で校正不一致)
- エネルギー閉合誤差を slope 補正で対応(Oran 0.74, Tarazona 0.71)

### 仮説の最終的な変遷(重要)
```
深根仮説(初期)
    ↓ v9-v12 で SDS=+0.05 が観測 → 一見支持
    ↓
Tarazona に灌漑あり(ユーザーから情報)
    ↓ 仮説と区別不能化
    ↓
v13: days_since_irrig 階層化 → 8+d で SDS 負の奇妙な値
    ↓ Claude が「8+d は実質"春期"」と気づく
    ↓
v14: 季節分離 + 夏期内 bucket → dose-response 判明
    ↓
最終結論: 灌漑依存(時定数 3-4日)、生物学的深根仮説は棄却
```

---

## 2. 使用データ

### 2.1 EC データ
| サイト | ファイル | 形式 | 期間 | 用途 |
|---|---|---|---|---|
| Oran | `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.xlsx` | xlsx | 2018-01 – 2020-12 | v4 が読込 |
| Oran | `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv` | csv | 同上 | v9/v10 が読込 |
| Tarazona | `/home/shion-nagamine/Dataset/Eddy data in Spain/EddyAlmond_Raw5years_withG.xlsx` | xlsx | 2020-06 – 2024-10 | v4 が読込(半時間データ) |
| Tarazona | `/home/shion-nagamine/Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv` | csv | 同上 | v9/v10/v13/v14 が読込、**`Irrig_mm`, `Rain_mm`, `IrrigRain_mm` を含む** |

### 2.2 衛星・補助データ
| 種別 | パス | 用途 |
|---|---|---|
| Sentinel-2 NDWI Oran | `/mnt/hdd/Dataset/Sentinel2_NDWI/Oran_NDWI_Export.csv` | v9/v10 で NDWI マッチ |
| Sentinel-2 NDWI Tarazona | `/mnt/hdd/Dataset/Sentinel2_NDWI/TzM_NDWI_Export.csv` | 同上 |
| ERA5 2m T/Td | `/mnt/hdd/Dataset/ERA5_2m_Temperature/{year}.nc` (2018-2024) | VPD 統一計算 |
| GRACE-FO | `/mnt/hdd/Dataset/GRACE-FO_TWL/GRCTellus.JPL.200204_202602.GLO.RL06.3M.MSCNv04CRI.nc` | v9/v10 で地下水貯留量参照(空間解像度 ~300km、参考程度) |

### 2.3 サイト座標
- Oran: lat 38.82, lon −1.86
- Tarazona: lat 39.266, lon −1.9397
- GRACE は経度を 0–360° 系に変換(−1.86 → 358.14)

### 2.4 使用カラム
**v4(半時間データ)**
- Oran: `TIMESTAMP`, `SWC_1_1_1`, `SWC_2_1_1`, `SWC_3_1_1`, `H`, `LE`, `G`, `ET`, `VPD`, `FC_mass`(GPP proxy)
- Tarazona: `TIMESTAMP`, `H`, `LE`, `G`, `ET`, `VPD`, `VPD_kPa`, `co2_flux`, `SWC_avg` または `SWC_1_1_1`, `SWC_2_1_1`

**v13/v14(Tarazona 日次)**
- `date`, `Irrig_mm`, `Rain_mm`, `IrrigRain_mm`

### 2.5 前処理(v4 が実施)
1. **物理的 hard-clip**: LE ∈ [-100, 900], H ∈ [-300, 1000], G ∈ [-300, 600] (W/m²)
2. **昼夜分離 MAD ×7**: 昼夜混合では昼間の正常高値が誤除去される → 9-17h と それ以外で分離して MAD 計算、閾値 7
3. **時間帯別 z-score 4σ**: 24時間別の正規分布フィルタ
4. **エネルギー閉合補正**: OLS slope を計算し、`LE_corr = LE/slope`, `H_corr = H/slope`
5. **日次集約**: 30 min → daily(平均または合計)
6. **EF 計算**: `EF = LE_corr / (LE_corr + H_corr)` (closure-corrected)
7. **4象限分類**: SWC<p40 & VPD>p60 で normal/soil_dry/atm_dry/compound に分類
   - 閾値はサイト×is_growing 別に計算

---

## 3. 実施した分析(時系列順)

### 分析 v4(土台、ユーザー作成)

**目的**
半時間データから日次・分類済データセットを生成し、4 象限干ばつ分類で日変化(diurnal cycle)を可視化する。

**実施内容**
- xlsx 読込 → クレンジング → 外れ値除去 → 閉合補正 → 日次集約 → 4 象限分類 → diurnal 図生成
- パッチ追加(Claude 提案):
  ```python
  daily.to_parquet(f"{OUT_DIR}daily_classified_v4.parquet")
  json.dump(closure_slopes, open(f"{OUT_DIR}closure_slopes_v4.json","w"))
  ```

**結果**
- Oran: closure slope=0.7376, R²~0.92
- Tarazona: closure slope=0.7098
- daily 1939 行(Oran 914, Tarazona 689 だが一部期間外で集約後 1939)
- 4 象限分類で diurnal が綺麗に分離

**解釈**
- Oran growing diurnal: soil_dry < normal < compound < atm_dry(典型的浅根反応)
- Tarazona growing diurnal: soil_dry ≈ normal、atm_dry ≈ compound(典型的木本性、ただし灌漑効果含む)

**問題点**
- 直接 hypothesis test の指標を出していなかった → v11 へ

**次の仮説**
SDS のような単一指標でサイト間比較できる定量解析を作る。

---

### 分析 v9: 3軸統合(SWC × NDWI × VPD)

**目的**
「SWC 低 ∧ NDWI 高 ∧ VPD 高」の **同時成立日**を `deep_access` 群として抽出し、深根の直接証拠を探す。

**実施内容**
- NDWI を Sentinel-2 CSV から ±5日でマッチング
- SWC/VPD/NDWI の閾値を自動決定(KDE 谷検出 or unimodal 中央値)
- 6 群分類: deep_access / true_drought / atm_driven / compound_dry / wet_deep / wet_shallow
- GRACE 経度 0–360 変換

**結果**
- Oran deep_access: **0 日**
- Tarazona deep_access: **0 日**(NDWI 閾値=0.109 とした場合)
- Tarazona は v9 では deep_access=86日確認(NDWI 閾値が違うバージョン)

**解釈**
3軸同時成立が稀すぎて、`deep_access` 群がほぼ作れない。

**問題点**
- NDWI 絶対閾値(0.0 や 0.109)はサイト依存・季節依存
- 「SWC 低 ∧ NDWI 高」は自然界では同時に起きにくい(負の相関がある)

**次の仮説**
NDWI を anomaly(月別中央値からの偏差)で評価する v10 へ。

---

### 分析 v10: 交絡対策版

**目的**
v9 の根本問題を修正:
1. LE 絶対閾値 50 W/m²+ NDWI active-canopy filter
2. **Drought Response Ratio**: `LE_dry / LE_wet` の比をサイト間で比較
3. NDWI anomaly 化
4. ラグ相関を 30日移動平均からの偏差で計算
5. 月別 phenology-matched 比較

**結果**
Oran 雨養と Tarazona 灌漑で **絶対 LE 値が桁違い**(Oran ~3 W/m² vs Tarazona ~195 W/m²)を比較してしまう設計問題が露呈。

**解釈**
Cross-site の絶対比較は根本的に意味がない:
- 作物が違う(穀物 vs アーモンド)
- LAI が違う(0-4 vs 通年~4)
- 季節も違う(Oran 冬作 vs Tarazona 通年)
- "深根 vs 浅根" ではなく "種・気候の違い"を見ているだけ

**問題点**
DRR は形は良いが、cross-site でなく within-site で使うべき。

**次の仮説**
各サイト内で SDS を計算し、サイト固有の baseline で評価する。

---

### 分析 v11/v12: within-site SDS(SDS/VAS/DSO/CompoundDrop)

**目的**
v4 の 4 象限分類をそのまま使い、各サイト内で 4 つの干ばつ応答指標を bootstrap 5000回 + 95% CI で計算。

**指標定義**
```
SDS  = (LE_normal − LE_soil_dry) / LE_normal     (≈0 = 深根)
VAS  = (LE_atm_dry − LE_normal) / LE_normal      (高 = VPD駆動)
DSO  = (LE_compound − LE_soil_dry) / LE_soil_dry (高 = 土壌乾燥下VPD駆動)
CompoundDrop = (LE_normal − LE_compound) / LE_normal
```

**実施内容(v11)**
- parquet + json 入力
- 4 指標を LE_corr / EF_corr / Bowen_corr / ET の各変数で計算
- Mann-Whitney pairwise tests
- 8 枚 PNG 出力

**v11 で発覚したバグ**
- `DENOM_FLOOR = 5.0` を全変数一律に適用 → EF(0-1)/Bowen(~1)/ET(<5) で全て NaN
- Fix: 変数別 `DENOM_FLOORS = {"LE_corr": 5.0, "EF_corr": 0.05, "Bowen_corr": 0.1, "ET": 0.3}`
- Bowen は `|Bowen|<20` で極値クリップ追加

**v12 で改善**
- 8 枚の図を自己説明的に再設計
- 休眠期警告(LE_normal < 30 W/m² で注意)
- 自動判定ロジック(CI 重なり判定)

**結果(v12 Growing season, LE_corr)**
| 指標 | Oran | Tarazona |
|---|---|---|
| SDS | +0.312 [+0.21, +0.41] | +0.051 [-0.07, +0.15] |
| VAS | +2.498 | +1.167 |
| DSO | +1.581 | +1.015 |
| CompoundDrop | -0.776 | -0.912 |

→ 一見、深根仮説支持。CI が分離している。

**問題点**
- ユーザー指摘: 「Oran の生育期 11-6月、Tarazona 1-10月で全く違う」→ cross-site での主張は危うい
- さらに後で発覚: 季節間プーリングのアーティファクト(spring の高 SWC が "normal" 群を支配、summer の低 SWC が "soil_dry" を支配)
- v12 の +0.05 SDS は **アーティファクト**だった

**次の仮説**
衛星データで深根仮説を補強できないか。

---

### 分析 v13: 灌漑経過日数で深根 vs 灌漑を切り分け

**契機**
ユーザーから「Tarazona は灌漑あり」の情報。
→ 深根仮説と灌漑仮説が観測上区別不能。Irrig_mm column が CSV にある。

**目的**
「最後の灌漑からの経過日数(days_since_irrig)」でデータを階層化し、灌漑離脱後も LE 維持されるかを検定する。

**実施内容**
- Tarazona daily CSV から `Irrig_mm > 0.5` を灌漑イベント定義
- Days since last irrigation を計算(NaN for never-irrigated)
- Bucket: 0-1d, 2-3d, 4-7d, 8+d
- 各 bucket 内で SDS 計算

**結果**
| Bucket | LE SDS | CI | n_n, n_s |
|---|---|---|---|
| 0-1d | +0.207 | [-0.06, +0.42] | 46, 42 |
| 2-3d | +0.176 | [-0.07, +0.40] | 37, 34 |
| 4-7d | +0.310 | [-0.01, +0.44] | 43, 41 |
| **8+d** | **−0.151** | [-0.47, +0.13] | **249**, 69 |

灌漑統計:
- 138 イベント 2020-06 〜 2024-10
- 月別: 5月=18, 6月=17, 7月=47, 8月=35, 9月=14, 10月=7
- 11月–4月: 灌漑停止(完全に)

**解釈**
8+d バケットの SDS が負・サンプル数が異常に多い(249/69)。
傾き−0.094、相関−0.61。
→ 自動判定が「逆相関、不自然」と判定。

**問題点(Claude が気づく)**
**8+d バケットの 249 サンプルは大半が 1-4月**(灌漑停止期=春期)。
バケット番号 ≈ 季節 になっていた。経過日数の純粋効果を測れていない。

**次の仮説**
季節 × バケットの交絡を分離する v14 へ。

---

### 分析 v14: 季節分離 + 夏期内 bucket(現在の最終版)

**目的**
v13 の交絡を解決:
- Test A: 夏期(7-9月)限定で bucket 解析 → 夏期内で純粋な経過日数効果
- Test B: 季節別(spring 1-4月 / shoulder 5,6,10月 / summer 7-9月)で SDS → 灌漑停止期の純粋テスト

**実施内容**
- 季節定義: SEASONS = {"spring (1-4月)": [1,2,3,4], "shoulder (5,6,10月)": [5,6,10], "summer (7-9月)": [7,8,9]}
- 夏期内 bucket × 変数(LE_corr/EF_corr/ET)で SDS bootstrap
- 季節 × サイト で SDS bootstrap(Oran が Test B の control)
- 月 × bucket の heatmap で交絡可視化
- Verdict 自動判定: |SDS|≤0.15 OR CI が 0 を含む → 深根帯

**結果(v14)**

**Test A: 夏期 bucket(Tarazona LE_corr)**
| Bucket | SDS | CI | n_n, n_s |
|---|---|---|---|
| 0-1d | +0.051 | [-0.12, +0.43] | 28, 25 |
| 2-3d | +0.004 | [-0.29, +0.36] | 24, 21 |
| **4-7d** | **+0.384** | **[+0.11, +0.51]** | 27, 22 |
| **8+d** | **+0.283** | **[+0.14, +0.43]** | 21, 28 |

→ **明確な dose-response**: 灌漑直後 3 日は SDS≈0、4 日以降は SDS が立ち上がり CI が 0 を超える。

**Test B: 季節別 SDS (LE_corr)**
| 季節 | サイト | SDS | CI | n_n, n_s |
|---|---|---|---|---|
| spring | Oran | **+0.689** | [+0.62, +0.72] | 103, 86 |
| spring | Tarazona | +0.616 | [-0.98, +0.83] | 181, **15** (⚠) |
| shoulder | Tarazona | +0.429 | [+0.29, +0.55] | 94, 75 |
| summer | Tarazona | +0.272 | [+0.14, +0.35] | 100, 96 |

**解釈**
1. **Tarazona の夏期は灌漑依存パターン**: 0-3d で SDS≈0、4+d で SDS=0.28-0.38 → 灌漑効果の dose-response
2. **真の深根なら 8+d でも SDS≈0 のはず → 実際は +0.28 → 深根仮説棄却**
3. **Tarazona の shoulder/summer は SDS≈0.27-0.43 で表層感受性あり** → v12 の "+0.05" はアーティファクト
4. **Spring Tarazona は n_s=15 で判定不能**(冬の降雨で SWC>p40 が大半)

**重要な発見: 自動判定ロジックのバグ**
`in_band` 関数が「CI が 0 を含む → 深根」と判定するが、**n_s=15 で CI が広すぎる場合に偽陽性**。
Spring Tarazona の CI [-0.98, +0.83] は単に sample 不足の現れで、深根の証拠ではない。
→ v15 で fix 必要: CI 幅 > 0.5 で判定保留、n_s < 30 で判定保留など

**次の仮説**
- 真の深根仮説は棄却
- **灌漑による SWC-ET decoupling time-constant ~3-4日** が新しい主張
- 衛星 ET 推定モデル(SMAP-based)への応用が contribution として残る

---

## 4. コード変更履歴

### v4 (ユーザー側、既存)
- ファイル: ユーザー側 `Eddy Covariance Flux Analysis v4` script
- 変更箇所: 末尾に 2 行追加(Claude 提案)
  ```python
  import json
  daily.to_parquet(f"{OUT_DIR}daily_classified_v4.parquet")
  with open(f"{OUT_DIR}closure_slopes_v4.json", "w") as f:
      json.dump(closure_slopes, f)
  ```
- 理由: v11 以降の中間ファイル連携のため
- 現在の状態: ユーザーが反映済みの想定。Claude 側 repo にはなし。

### v9 (`analysis_A_v9.py`)
- 場所: `/home/user/bakanposs/analysis_A_v9.py`
- 内容: 3軸統合(SWC×NDWI×VPD)で deep_access 抽出を試行
- 状態: 完成、ただし deep_access=0 で目的未達

### v10 (`analysis_A_v10.py`)
- 場所: `/home/user/bakanposs/analysis_A_v10.py`
- 内容: DRR + NDWI anomaly + lag correlation
- 状態: 完成、ただし cross-site 比較の根本問題が残る

### v11 (`analysis_A_v11.py`)
- 場所: `/home/user/bakanposs/analysis_A_v11.py`
- 内容: 4-class SDS/VAS/DSO/CompoundDrop bootstrap
- 変更履歴:
  - 初版: DENOM_FLOOR=5.0 一律(バグ、EF/Bowen/ET 全 NaN)
  - 修正: `DENOM_FLOORS` を変数別に変更、`BOWEN_CLIP=20` 追加
- 状態: 完成

### v12 (`analysis_A_v12.py`)
- 場所: `/home/user/bakanposs/analysis_A_v12.py`
- 内容: 8 図を自己説明的に再設計、休眠期警告、自動判定
- 状態: 完成。**結果は今思えば季節プーリング artifact**

### v13 (`analysis_A_v13.py`)
- 場所: `/home/user/bakanposs/analysis_A_v13.py`
- 内容: 灌漑経過日数で階層化、Tarazona daily CSV `Irrig_mm` 利用
- 状態: 完成。8+d バケットの交絡が判明、v14 へ。

### v14 (`analysis_A_v14.py`)
- 場所: `/home/user/bakanposs/analysis_A_v14.py`
- 内容: 季節分離 + 夏期内 bucket、Oran spring を control
- 状態: 完成。**現時点の最終解析**
- 既知のバグ: `in_band` ロジックが「CI が 0 を含む」を無条件に深根と判定する偽陽性

### requirements.txt
- 場所: `/home/user/bakanposs/requirements.txt`
- 内容: v4-v14 全部の依存 (numpy, pandas, scipy, matplotlib, xarray, netCDF4, h5netcdf, dask, bottleneck, openpyxl, pyarrow)
- 状態: 完成・push 済

---

## 5. 試して失敗したこと(再試行禁止)

### F1: NDWI 絶対閾値での deep_access 抽出(v9)
- 試行: `NDWI > 0.0` または `> 0.109` を「植被が水を含む」基準にした
- 結果: deep_access 日が 0
- 理由: NDWI 絶対値はサイト・季節依存性が強すぎる。Oran は全体的に負、Tarazona は全体的に正で、固定閾値は意味を成さない
- 教訓: 必ず anomaly(月別中央値からの偏差)を使う

### F2: cross-site 絶対 LE 比較(v10)
- 試行: Oran 雨養日と Tarazona 灌漑日の絶対 LE を Mann-Whitney 検定 → p=1.5e-47
- 結果: 統計的に超有意だが、科学的に無意味
- 理由: Oran 雨養 cereal の dry SWC 日は収穫後/裸地、Tarazona は灌漑活動中で、種・LAI・季節すべてが異なる
- 教訓: cross-site の絶対比較はしない。各サイト内で baseline 比を取る

### F3: 季節間プーリング SDS(v12)
- 試行: Tarazona の生育期 1-10月を全部 1 群として SDS 計算 → SDS=+0.05(深根支持に見えた)
- 結果: v14 で季節分離すると spring/shoulder/summer 各 +0.27〜+0.62 で実は表層感受性あり
- 理由: spring が高 SWC × 低 LE で "normal" 群を支配、summer が低 SWC × 高 LE で "soil_dry" 群を支配 → "soil_dry" の方が見かけ上 LE が高くなる artifact
- 教訓: SWC 閾値や群定義は **季節別に**評価する

### F4: ラグ相関(v9)
- 試行: SWC と NDWI の lag correlation で深根を検出
- 結果: 季節周期に支配されてフラット、解釈不能
- 理由: SWC と NDWI 両方が年周期を持つので Pearson は周期成分で決まる
- 教訓: lag 解析は必ず detrend / deseasonalize してから

### F5: GRACE-FO を point sample で利用(v9)
- 試行: lat=38.82, lon=-1.86 で GRACE TWS を抽出
- 結果: 経度系の不一致で 0 ヶ月分(NaN)
- Fix: 経度 0-360 系に変換(-1.86 → 358.14)
- 教訓: GRACE は ~300km 解像度なのでサイト固有の解析には適さない(参考程度)

### F6: 自動判定の `in_band` 偽陽性(v12, v14 で再発)
- 試行: 「CI が 0 を含む → 深根帯」と単純判定
- 結果: n が小さく CI が広いケース(Tarazona spring n_s=15)で誤った "深根支持" 判定
- 改善案: CI 幅 > 0.5 で保留、n_s < 30 で保留
- 状態: v15 で fix 予定

### F7: DENOM_FLOOR 一律 5.0(v11 初版)
- 試行: 全変数で同じ最小分母閾値
- 結果: EF/Bowen/ET の指標が全 NaN
- Fix: 変数別 DENOM_FLOORS
- 教訓: スケールが違う物理量に同じ定数閾値を使わない

---

## 6. 現在の研究上の論点

### 論点 A: 深根仮説の最終判定
**現在のスタンス**: 棄却寄り
- 夏期 4+d bucket で SDS が有意に正 → 灌漑効果が消えると Tarazona も表層感受性を持つ
- 真の深根なら経過日数によらず SDS≈0 のはず
- ただし春期データ不足(n_s=15)が完全否定を弱める

### 論点 B: 春期の解釈
- Tarazona spring の SDS = +0.62 [-0.98, +0.83] は **判定不能**
- 解釈案 1: 冬の貯留水で SWC>p40 が大半 → 干ばつ条件成立せず → 深根仮説検定不能
- 解釈案 2: 春は花期/葉展開期で transpiration phenology が違う → SDS の意味自体が違う
- どちらにせよ春期では深根の決定的証拠にも反証にもならない

### 論点 C: Irrigation Decoupling Time-constant の意義
- 0-3d で SDS≈0、4+d で SDS が立ち上がる → time-constant ~3-4日
- これは **灌漑農地で衛星 SWC ET 推定モデルの系統誤差を生む現象**
- 論文化の新しい主軸候補

### 論点 D: Phenology vs Hydrology
- Oran と Tarazona の比較は phenology(作物の生育周期)も違う
- "Spring SDS" は両サイトとも生育期だが、Oran は cereal 中期で Tarazona は葉展開期
- 同じ "spring" の SDS=+0.69 vs +0.62 は両者似ているが、phenology 差が交絡しているか

---

## 7. 未解決課題(優先順)

### P1(最優先): 自動判定ロジックのバグ修正
- `in_band()` が CI 幅・サンプル数を考慮しない
- v14 の Tarazona spring 「深根支持」結論は誤り
- 修正案: `(lo <= 0 <= hi) AND (hi - lo < 0.5) AND n_s >= 30`

### P2: Recovery time analysis(灌漑応答時定数の精緻化)
- 現在の 0-1d / 2-3d / 4-7d / 8+d バケットは粗い
- 1日刻みで SDS をプロット → exponential fit → time-constant τ を抽出
- 論文の headline figure 候補

### P3: 春期データ不足の補正
- Tarazona spring で SWC<p40 が n=15 しかない
- Spring 内で別の SWC 閾値(例: p25)を使えば soil_dry サンプル増える可能性
- ただし定義変更による比較性損失とトレードオフ

### P4: Phenology 補正
- Oran vs Tarazona の比較は phenology 交絡あり
- LAI または NDVI で正規化した LE/EF 指標を併用検討
- データに NDVI あれば実装可能

### P5: 衛星補強(別データソース)
- MODIS LST × ERA5 Tair で ΔT(canopy stress)を独立検証
- ECOSTRESS 70m thermal ET で空間検証
- Sentinel-1 SAR で灌漑パルスを Irrig_mm と独立に検出
- どれも追加データ取得が必要

### P6: 多深度 SWC の必要性
- 5cm のみでは drip wetted bulb(20-50cm)を直接観測できない
- 既存タワーに 30cm/50cm SWC センサーあれば再解析価値大
- データ要確認

### P7: 論文骨子の準備
- 仮説変更後の新フレーミング: "Drip irrigation decouples SWC-ET on 3-4 day timescale"
- Introduction, Methods, Results, Discussion の骨組み
- v14 の図を中心に構成

---

## 8. 次セッション開始時にやること(順序)

1. **`SESSION_SUMMARY.md` を必ず最初に読む**(これ)
2. `/home/user/bakanposs/` の状態確認: `git log --oneline -20`, `ls`
3. 中間ファイル状態確認: `daily_classified_v4.parquet`, `closure_slopes_v4.json` が存在するか(v4 実行済みか)
4. v14 の output 確認: `output_analysis_A_v14/fig01-04.png` を視認
5. **ユーザーに優先課題を確認**:
   - P1(verdict bug fix)で v15 を書くか?
   - P2(recovery time analysis)を実装するか?
   - 論文骨子(Path 1)を書き始めるか?
   - 衛星補強(MODIS/ECOSTRESS, P5)を試すか?
6. **絶対に提案しないこと**:
   - cross-site の絶対 LE 比較(F2 で失敗済み)
   - NDWI 絶対閾値(F1 で失敗済み)
   - 季節プーリングなしの SDS(F3 で artifact 確認済み)

---

## 9. Claude への引き継ぎ指示

### 9.1 一行サマリー
"Mediterranean almond(Tarazona, 灌漑)と cereal(Oran, 雨養)の EC データで、当初「深根仮説」を検定したが、季節 + 灌漑経過日数の制御後、Tarazona は **灌漑による 3-4日の SWC-ET decoupling** を示すことが判明。深根仮説は棄却寄り、新しい論文方向性は灌漑による衛星 ET 推定モデルへの含意。"

### 9.2 重要な判断と制約
- **絶対にしてはいけない比較**: Oran と Tarazona の同日絶対 LE 比較(種・LAI・季節すべて違う)
- **必ず seasonal-stratified に**: 全季節プーリングは v12 で artifact を生んだ
- **n_s < 30 のときは結論保留**: spring Tarazona の偽陽性事例あり
- **SWC センサーは 5cm のみ**: drip wetted bulb 直接観測ではない、間接シグナルとして扱う

### 9.3 ユーザーの好み
- 図は **自己説明的**(タイトルに結論、軸ラベル明確、サンプル数表示)を強く好む
- 解析は **段階的に**(各バージョンが前の問題を1つ修正)
- バグ・問題点は**正直に指摘**することを歓迎(誤った "深根支持" の自動判定を指摘した経緯あり)
- コードは **コメントで "なぜそうしたか" を書く**ことを好む
- バージョン管理: `analysis_A_vN.py` と output ディレクトリ `output_analysis_A_vN/` で並走

### 9.4 当面の科学的方向性
- 論文化の方向: **"Drip irrigation decouples surface SWC from canopy ET on a 3-4 day timescale"**
- v14 の Test A 夏期 bucket dose-response が headline figure 候補
- 衛星 ET モデル(SMAP-based)への含意が contribution

### 9.5 注意事項
- v12 の "Tarazona SDS=+0.05 → 深根支持" は **artifact、引用しない**
- v13 の "SDS 逆相関、不自然" の自動判定は **間違い**、8+d が春期だった
- v14 の Tarazona spring "深根支持" 判定も **偽陽性**、n_s=15 が問題

---

## 10. 重要ファイル一覧

### コード
| パス | 内容 |
|---|---|
| `/home/user/bakanposs/analysis_A_v9.py` | NDWI 3軸統合(失敗例として保存) |
| `/home/user/bakanposs/analysis_A_v10.py` | DRR + cross-site(失敗例として保存) |
| `/home/user/bakanposs/analysis_A_v11.py` | within-site SDS(DENOM_FLOOR bug fix 済) |
| `/home/user/bakanposs/analysis_A_v12.py` | 自己説明的 figures(8 枚) |
| `/home/user/bakanposs/analysis_A_v13.py` | 灌漑 days_since_irrig 階層化 |
| `/home/user/bakanposs/analysis_A_v14.py` | **現最終版**: 季節分離 + 夏期 bucket |
| `/home/user/bakanposs/requirements.txt` | 全環境依存 |
| `/home/user/bakanposs/SESSION_SUMMARY.md` | **このドキュメント** |

### 中間ファイル(ユーザー側で生成、Claude 側 repo にはなし)
| パス | 内容 | 必要なバージョン |
|---|---|---|
| `daily_classified_v4.parquet` | v4 が出力する日次・分類済データ | v11-v14 が読込 |
| `closure_slopes_v4.json` | エネルギー閉合スロープ {Oran:0.7376, Tarazona:0.7098} | v11-v14 が読込 |

### 入力データ(ユーザー側 /home/shion-nagamine/, /mnt/hdd/)
| パス | 内容 |
|---|---|
| `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.xlsx` | Oran 半時間 EC(v4 用) |
| `/home/shion-nagamine/Dataset/Eddy data in Spain/EddyAlmond_Raw5years_withG.xlsx` | Tarazona 半時間 EC(v4 用) |
| `/home/shion-nagamine/Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv` | Tarazona 日次(`Irrig_mm`, `Rain_mm`, `IrrigRain_mm` 含む。v9/v10/v13/v14 用) |
| `/mnt/hdd/Dataset/Sentinel2_NDWI/Oran_NDWI_Export.csv` | Sentinel-2 NDWI(v9/v10 用) |
| `/mnt/hdd/Dataset/Sentinel2_NDWI/TzM_NDWI_Export.csv` | 同上 Tarazona |
| `/mnt/hdd/Dataset/ERA5_2m_Temperature/{year}.nc` (2018-2024) | ERA5 VPD 計算用 |
| `/mnt/hdd/Dataset/GRACE-FO_TWL/GRCTellus.JPL.200204_202602.GLO.RL06.3M.MSCNv04CRI.nc` | GRACE-FO TWS(参考) |

### 出力(ユーザー側)
| ディレクトリ | 内容 |
|---|---|
| `output_analysis_A_v9/` | v9 figures + CSV |
| `output_analysis_A_v10/` | v10 figures + CSV |
| `output_analysis_A_v11/` | v11 figures + CSV |
| `output_analysis_A_v12/` | v12 figures + CSV |
| `output_analysis_A_v13/` | v13 figures + CSV |
| `output_analysis_A_v14/` | **最終結果**: fig01-04 + CSV 2本 |

### 重要 git コミット
- `e9efc97` v9 added
- `2069be8` v10 + requirements
- `9ca8f73` v11 added
- `f547c3d` requirements consolidated
- `7d90947` v11 fix (DENOM_FLOORS)
- `0c12aea` v12 refactored
- `bcaf1a3` v13 irrigation analysis
- `1dc0ecb` v14 seasonal × bucket
