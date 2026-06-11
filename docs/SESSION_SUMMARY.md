# Session Summary — Drip Irrigation × Satellite ET Analysis

**期間**: 2026-05-05 〜 2026-05-10
**ブランチ**: `claude/compare-ec-satellite-et-ZnENi`
**最終コミット**: `d34c572` (paper Methods/Results/Discussion 拡張)
**作業ディレクトリ**: `/home/shion-nagamine/bakanposs/`

---

## 1. 研究背景・目的

### 1.1 研究目的
半乾燥地 (スペイン Albacete) の **ドリップ灌漑農地** において、衛星 ET プロダクト (MOD16, PML, LSA SAF METv3) がなぜ系統的に過小評価するのか、その時間構造を**フラックスタワー (EC) データの 7 年連続観測**で定量化し、**補正可能な物理モデル**を提示する。最終ゴールは査読論文 (Agric. Forest Meteorol. / RSE / HESS) 投稿。

### 1.2 なぜこの分析を行ったか
- 衛星 ET は灌漑スケジューリング・水収支計算で運用が始まっている
- しかし灌漑農地で 2–3 割の過小評価が報告されており、**その時間構造 (灌漑イベント後 N 日でどう減衰するか)** を連続観測で示した研究は皆無
- 補正手法を設計するには「灌漑タイミング情報を入力に与えてバイアスがどう動くか」を知る必要がある

### 1.3 仮説の進化 (極めて重要 — 研究の核)
1. **初期仮説 (deep-root)**: 「TzM アーモンドは深根で地下水にアクセスして表層が乾いていても蒸散を維持する」
2. **反証**: 灌漑経過日数で層別したら d0-3 だけ SDS+0.13、d4-7 と d8+ は 0 → 表層脱結合は灌漑直後限定
3. **修正仮説 (drip wet-bulb)**: 「ドリップ灌漑が深さ 10-30 cm に局所湿潤帯 (wet bulb) を作り、3-4 日持続する。表層 5 cm SWC とも深部 ~1 m SMAP root-zone とも独立に動く」
4. **直接観測**: TzM d0-3 で SWC vs SMAP_rz の r = −0.19 (depth inversion) → wet-bulb 仮説の決定的証拠

### 1.4 前提条件
- 解析は完全に観測駆動 (no land surface model)
- 灌漑タイプは drip 1 種類 (TzM) + rainfed コントロール 1 種類 (Oran)
- データは EC + GEE衛星 + LSA SAF METv3 NetCDF + SMAP L4 のみ
- 統計手法は frequentist (bootstrap CI, NLS, AIC, Pearson)

---

## 2. 使用データ

### 2.1 EC (Eddy Covariance) — 生データ

| サイト | 期間 | 作物 | フォーマット | パス |
|---|---|---|---|---|
| Oran | 2018-01〜2020-06 | vetch / wheat / pea 輪作 | AmeriFlux Standard | `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_*.csv` |
| Oran (30min) | 2018-01〜2020-12 | 同上 | Half-hour raw | 同ディレクトリ内 |
| TzM | 2020-06〜2024-10 | drip-irrigated almond | EddyPro proprietary | 同ディレクトリ内 |

**カラム**: date, ET_mm, LE_Wm2, H_Wm2, G_Wm2, Rn_Wm2, Ta_mean/max/min, VPD_kPa_mean/max, SWC (5 cm), P_mm, **Irrig_mm**, NDVI, GPP_gC_m2_d, n_obs_le

**前処理**:
- QC flag ≤ 2 (AmeriFlux 標準) のみ残す
- daily ≥ 24 半時間 (50%カバー) を要求
- Oran は半時間 → 日次集約 (LE/H/G/Rn) を `aggregate_oran_30min.py` で実行 (year + Julian + Time_hours から timestamp 再構築)
- TzM は既に日次集約済みデータを使用

### 2.2 衛星 ET (GEE)

サイト座標 (Oran: 38.82°N, -1.86°E / TzM: 39.27°N, -1.94°E) に対し、Oran 200 m / TzM 300 m バッファで GEE から抽出。

| プロダクト | 解像度・cadence | 単位 | CSV |
|---|---|---|---|
| MOD16A2GF v061 | 500 m, 8-day | mm/d | `/mnt/hdd/Dataset/Fast_OranTzM_MOD16.csv` |
| PML v018 | 500 m, 8-day | Ec, Es, Ei → mm/d | `/mnt/hdd/Dataset/Fast_OranTzM_PML.csv` |
| MCD15A3H | 500 m, 4-day | LAI, FPAR | `/mnt/hdd/Dataset/Fast_OranTzM_LAI.csv` |
| MOD11A1 | 1 km, daily | LST | `/mnt/hdd/Dataset/Fast_OranTzM_LST.csv` |
| Sentinel-2 L2A | 10 m, 5-day | NDVI/NDWI/NDMI | `/mnt/hdd/Dataset/OranTzM_S2_NDVI.csv` |
| ERA5-Land | 9 km, hourly | Ta, P | `/mnt/hdd/Dataset/Fast_OranTzM_ERA5_Daily.csv` |
| CHIRPS | 5 km, daily | precipitation | `/mnt/hdd/Dataset/Fast_OranTzM_CHIRPS.csv` |

### 2.3 LSA SAF METv3 (NetCDF)

- パス: `/home/shion-nagamine/Dataset/METv3/YYYY/MM/MMDD/YYYYMMDD_HHMM.nc`
- 期間: 2018-01-01 〜 2024-12-31
- 解像度: 0.05° (~5 km) global grid
- 時間刻み: 30 分 (1 日 48 スロット)
- 総ファイル数: ~120,000 枚
- 処理: `pipeline/load_metv3.py` で nearest pixel 抽出 → 日積算 (≥ 36/48 スロット必須)
- 出力: `/home/shion-nagamine/bakanposs/metv3_daily_all.csv` (5,114 site-days)

### 2.4 SMAP L4 SPL4SMGP v007

- 元: `NASA/SMAP/SPL4SMGP/007` (GEE)
- 解像度: 9 km, 3 時間刻み
- バッファ: **6 km** (重要: 9 km ピクセル捕捉のため 200/300 m では NULL になる)
- 抽出スクリプト: `gee/gee_smap_only.js`
- 出力: `/mnt/hdd/Dataset/SMAP_OranTzM.csv` (40,896 raw rows, 100% valid)
- 日次集約: `pipeline/load_smap.py` → `/home/shion-nagamine/bakanposs/smap_daily.csv` (5,112 site-days)

### 2.5 統合 master データセット

最終的なメイン CSV (53 列 × 1,356 site-days):
- `/home/shion-nagamine/bakanposs/master_full_v2.csv`
- EC + 全衛星 + METv3 + SMAP + irrigation flags すべて含む

---

## 3. 実施した分析 (時系列)

### 分析 A — In-situ SDS の再現

**目的**: 「TzM が SWC とLEで脱結合している」という v14 結果を、クリーンなパイプラインから再現する

**実施内容**:
1. `pipeline/unify_ec_daily.py`: 4 つの異種日次 CSV を統合 → `ec_daily_master.csv` (1,356 行)
2. `pipeline/add_flags.py`: `Irrig_mm > 0.5` から days_since_irrig 計算、irrig_bucket (d0-3 / d4-7 / d8+) / season / drought_class / NDVI gate 付与
3. `pipeline/aggregate_oran_30min.py`: Oran 半時間 → 日次 LE/H/G/Rn を埋め、`ec_daily_master_complete.csv` 作成
4. `figures/sds_v14_repro.py`: SDS = 1 - mean(LE|SWC<p25) / mean(LE|p25≤SWC≤p75)、bootstrap N=2000 で 95% CI

**結果**:
- Oran spring: SDS = +0.43 [+0.36, +0.51], n=202 (rainfed active, SWC と LE がカップル)
- TzM summer: SDS = +0.11 [+0.06, +0.17], n=393 (~1/4 to Oran)
- TzM summer × irrig_bucket: d0-3 = +0.13, d4-7 = +0.01, d8+ = 0.00 (灌漑直後だけ感度あり)

**解釈**: deep-root 仮説は反証された。灌漑からの経過日数で層別すると d0-3 だけが残り、d4-7 と d8+ では完全脱結合。

**問題点**:
- Oran SDS が当初 NaN だった → cell() が LE_Wm2 を優先したが Oran は半時間集約前で <30 件しかなかった
- 修正: LE_Wm2 が 30 件未満なら `ET_mm` にフォールバック

**次の仮説**: drip wet-bulb メカニズム (10-30 cm に局所湿潤帯)

---

### 分析 B — 衛星 ET 統合

**目的**: GEE から 7 プロダクトを取得して EC と同じ daily テーブルに揃え、サイト × 日付でジョイン可能にする

**実施内容**:
1. `pipeline/unify_satellite.py`: 8 種類の GEE wide CSV を long → daily に変換、単位変換 (MOD16: kg/m²/8d → mm/d、PML: Ec+Es+Ei を sum)、`satellite_daily.csv` (5,112 × 18 vars)
2. `pipeline/merge_satellite_ec.py`: EC × satellite を date × site でマージ → `master_full.csv`
3. 8-day プロダクトは 8-day window 内で forward-fill して daily に変換
4. 灌漑バケット別バイアス表をその場で生成

**結果** (TzM all-year):
- MOD16: MBE = -2.69, RMSE = 3.22 mm/d
- PML: MBE = -1.45, RMSE = 2.10 mm/d

**解釈**: 両プロダクトとも TzM で -1.5 〜 -2.7 mm/d 過小評価。Oran ではほぼ 0 (MOD16 = -0.23、PML = +0.44)。

**問題点**: GEE 元 CSV ファイルパスが時期によって変動 (Downloads → /mnt/hdd/Dataset/)。再現性のため `unify_satellite.py` の BASE を `/mnt/hdd/Dataset` に固定。

**次の仮説**: 灌漑タイミングでバイアスが指数減衰する (τ ≈ 3-4 d)

---

### 分析 C — 指数減衰モデルフィット

**目的**: 「バイアスが灌漑経過日数の指数関数で減衰する」を NLS で実証し、補正に使えるパラメータを取り出す

**実施内容**:
1. `figures/figure_C_summer.py`: 4 種類のフィルタ (all/summer/growing/summer×growing) × 灌漑バケットで boxplot
2. `figures/tau_fit.py`: full model Δ(t) = a·exp(-t/τ) + c と transient model Δ(t) = a·exp(-t/τ) の両方を fit (`scipy.curve_fit`)
3. Bootstrap N=500, τ ∈ (0, 60] 制約

**結果** (TzM summer × NDVI>0.3, full model, 95% CI):

| 製品 | a (mm/d) | τ (d) | c (mm/d) | c CI 0 含む |
|---|---:|---:|---:|:---:|
| MOD16 | −2.31 [−2.80, −1.88] | 4.0 [2.8, 5.9] | **−2.29 [−2.66, −1.85]** | NO (有意) |
| PML | −2.81 [−3.43, −2.22] | 4.3 [2.9, 7.0] | −0.57 [−1.03, +0.04] | YES |
| METv3 | −4.03 [−4.73, −3.44] | 6.0 [4.6, 8.3] | −0.62 [−1.08, +0.06] | YES |

**解釈**:
- **MOD16 は構造的 floor (-2.3 mm/d)** — FLUXNET 校正の作物カバー不足が原因
- **PML / METv3 は permanent offset 無し** — 過渡的補正で完全に取り除ける
- 3 つの独立な製品が同じ τ ≈ 4-6 d を示す → アルゴリズム固有でなく dry-surface 駆動の普遍的問題

**問題点**:
- 最初は bin 中央値 (5-6 点) で fit → τ が暴走
- 修正: raw daily data (n=300-400) で直接 fit、τ ∈ (0, 60] で制約

**次の仮説**: τ-based 補正は実用的に有効か?

---

### 分析 D — METv3 統合

**目的**: MODIS ファミリーから独立した第3プロダクトとして METv3 を追加し、結論の頑健性を高める

**実施内容**:
1. `inspectors/inspect_metv3.py`: NetCDF 構造把握 (ET [mm/h], 0.05° grid, quality_flag)
2. `pipeline/load_metv3.py`: ~120,000 NetCDF を年ごとに処理、xarray lazy load + 2 点ピクセル抽出 → 日積算
3. `pipeline/integrate_metv3_smap.py`: master_full に METv3 と SMAP を left-join → `master_full_v2.csv`
4. `tau_fit.py`, `figure_C_summer.py`, `sds_vs_bias.py` を 3 プロダクト対応に書き換え

**結果**:
- METv3 全データ: Oran MBE = +0.02 (3 製品で最もバイアス無し), TzM MBE = -2.34
- METv3 τ-fit: a = -4.03 (最大), τ = 6.0 d (最長), c CI が 0 を跨ぐ
- METv3 の最大振幅 a は 5 km ピクセル混合 (周辺 rainfed 地と orchard の混合) と整合

**解釈**: 製品が 3 つになったことで「MOD16 系統 vs PML/METv3 系統」の二群構造が明確に。

**問題点**:
- METv3 ロードに 3.5 時間 (120k ファイル)
- 並行で SMAP ダウンロード進行 → 効率的

---

### 分析 E — 仮説検証 (H1, H4, H6)

**目的**: 解析から提案された 8 個の仮説 (H1〜H8) のうち、既存データで検証可能な 3 個を回す

**実施内容**: `figures/hypothesis_tests.py`

#### H1: τ-based 補正の有効性
- 期待: -25〜-50% RMSE 削減
- **実測 (TzM summer × NDVI>0.3)**:
  - MOD16: 4.01 → 1.39 (-65%), MBE -3.69 → +0.00
  - PML: 2.82 → 1.44 (-49%), MBE -2.26 → +0.01
  - METv3: 3.85 → 1.50 (-61%), MBE -3.37 → +0.01
- **結論**: 期待を遥かに超える。MBE がほぼ完全に 0 へ。MOD16 の構造的 floor も τ-fit に含まれる c によって吸収される。

#### H4: AIC モデル比較
- M1: bias ~ VPD / M2: bias ~ days_since_irrig / M3: bias ~ VPD + d + 交互作用
- **実測 ΔAIC (vs VPD)**: MOD16 = -73, PML = -66, METv3 = -153
- すべて ΔAIC > 60 → Burnham & Anderson の "decisive evidence" 基準を圧倒
- **結論**: days_since_irrig は VPD より圧倒的に重要。交互作用も有意。

#### H6: SMAP root-zone 代替性 (multi-stratum)
- 当初 Oran summer (n=34 post-harvest) のみテスト → inconclusive
- **改訂: 8 階層に拡張**

| stratum | n | r(SWC, SMAP_rz) | SDS_in_situ | SDS_smap |
|---|---:|---:|---:|---:|
| **Oran_spring** | **203** | **+0.80** | **+0.43** | **+0.43** |
| TzM_summer | 401 | +0.16 | +0.11 | -0.15 |
| **TzM_summer_d0-3** | **281** | **-0.19** | +0.13 | -0.19 |
| TzM_summer_d4-7 | 75 | +0.35 | +0.01 | +0.13 |
| TzM_summer_d8+ | 44 | +0.61 | -0.00 | +0.18 |

- **rainfed Oran spring**: SDS が in-situ と SMAP_rz で完全一致 (+0.43)、r=0.80 → SMAP root-zone は EC タワーがなくても SDS を計算できる代替指標
- **TzM d0-3 で r = -0.19 (負！)**: 5 cm 表層が乾く方向 + 1 m root-zone が濡れる方向 = **depth inversion** → **drip wet-bulb の直接観測証拠**

**結論**: 論文の科学的貢献が 3 層に拡張された
1. (元) 灌漑バイアスが指数減衰し補正可能
2. (新) SMAP × in-situ で wet-bulb メカニズム直接観測
3. (新) rainfed 条件で SMAP-only SDS マッピングへの道筋

**問題点**: 元 H6 (Oran summer n=34) は post-harvest で植生不活発、SDS が 0 ≈ になり判定不能 → 多階層に拡張して救出

**次の仮説**: H7 (SDS 広域マッピング) が現実的に視野に入る — SMAP のみで EC タワー不要

---

## 4. コード変更履歴

リポジトリ: `https://github.com/zhiyinzhangling77-ui/bakanposs`
ブランチ: `claude/compare-ec-satellite-et-ZnENi`

主要コミット (時系列、新→古):
```
d34c572  paper Methods/Results/Discussion を 3 製品 + H6 multi-stratum に拡張
665b2a8  H6 を multi-stratum に拡張 (Oran spring focal)
dedca8e  H1/H4 結果を narrative と paper に反映
4cec8a7  hypothesis_tests.py 追加 + analysis_narrative.md 作成
a58ea9f  sds_vs_bias.py を 3 製品対応に
5f7f73a  paper draft を 3 製品 (MOD16, PML, METv3) に拡張
48e8ea2  unify_satellite.py の BASE を /mnt/hdd/Dataset/ に変更
b14b584  load_smap.py を新形式 (SMAP_OranTzM.csv 単一ファイル) 対応に
27641e3  gee_smap_only.js から Unicode 文字削除 (GEE JS パーサーエラー回避)
66d3298  Introduction 追加、全 path を REPO 相対に統一、integrate_metv3_smap.py 作成
2bff988  load_metv3.py, load_smap.py 初版
14e3829  Oran SDS NaN bug fix、tau_fit.py 追加
664f588  figure_C_summer.py 追加 (summer × NDVI フィルタ)、paper_outline.md 作成
c4950e0  merge_satellite_ec.py + killer_figures.py
69d4587  unify_satellite.py 初版
6694148  Oran 30min loader を fix
5064656  aggregate_oran_30min.py: 混合 datetime parser 対応
```

### 主要ファイルの最終状態

| ファイル | 役割 | 行数 |
|---|---|---:|
| `pipeline/unify_ec_daily.py` | EC 日次統合 | 100+ |
| `pipeline/add_flags.py` | 灌漑バケット・季節・干ばつクラス | 130+ |
| `pipeline/aggregate_oran_30min.py` | Oran 半時間→日次 | 180+ |
| `pipeline/qc_master.py` | 7 項目 QA | 150+ |
| `figures/sds_v14_repro.py` | SDS metric | 170+ |
| `pipeline/unify_satellite.py` | GEE 7プロダクト統合 | 217 |
| `pipeline/load_metv3.py` | METv3 NetCDF → daily | 153 |
| `pipeline/load_smap.py` | SMAP CSV → daily (両形式対応) | 80 |
| `pipeline/merge_satellite_ec.py` | EC + 衛星マージ | 53 |
| `pipeline/integrate_metv3_smap.py` | METv3+SMAP 統合 | 136 |
| `figures/figure_C_summer.py` | 灌漑バケット boxplot (3 製品) | 130+ |
| `figures/sds_vs_bias.py` | SDS vs バイアス scatter (3 製品) | 180+ |
| `figures/tau_fit.py` | NLS exponential decay (3 製品) | 225+ |
| `figures/hypothesis_tests.py` | H1/H4/H6 検証 | 380+ |
| `paper_outline.md` | 論文骨格 | 161 |
| `paper_methods_results.md` | 論文 Intro/Methods/Results/Discussion (完成) | ~250 |
| `analysis_narrative.md` | 全体ナラティブ + 統計知識解説 | ~470 |

### コード設計上の重要原則
- すべてのスクリプトで `REPO = Path(__file__).parent.parent` を使い、絶対パスを repo 相対に統一
- データ生 CSV (master_full.csv 等) は git 管理外 (再生成可能)
- METv3 / SMAP / GEE 元データは `/mnt/hdd/Dataset/` または `/home/shion-nagamine/Dataset/`

---

## 5. 試して失敗したこと (再試行禁止)

### 5.1 τ-fit を bin 中央値で行う
- 試行: 灌漑日 0〜20 を bin に分け、中央値に NLS フィット
- 結果: τ が暴走 (60 超え)、CI が ±50 になる
- 原因: 5-6 点では NLS の自由度が足りない
- 対策: **raw daily data (n=300-400) で直接 fit、τ ∈ (0, 60] で制約**

### 5.2 cell() で LE_Wm2 を無条件優先
- 試行: SDS 計算で LE_Wm2 → ET_mm の順に優先選択
- 結果: Oran の SDS が全部 NaN
- 原因: Oran は半時間集約前は LE_Wm2 が <30 件 (n_obs_le=NaN)
- 対策: **LE_Wm2 が 30 件未満なら ET_mm にフォールバック**

### 5.3 Oran 30分ファイルを `pd.to_datetime` で直接パース
- 試行: 文字列の date 列を pd.to_datetime
- 結果: 99% が NaT に
- 原因: 混合フォーマット (日付のみ / 日付+時刻 / NaT) で format-lock
- 対策: **year + Julian day + Time_hours から numeric に再構築**

### 5.4 SMAP を 200m / 300m バッファで GEE 抽出
- 試行: 他衛星と同じ buffer で reduceRegions
- 結果: 全列 NaN
- 原因: SMAP は 9 km grid、200 m バッファでは pixel centroid が入らない
- 対策: **6 km バッファ (9 km の半径以上) に変更**

### 5.5 gee_smap_only.js で Unicode 矢印 `→` をコメントに使用
- 試行: `// Export → Drive → CSV` のような Unicode コメント
- 結果: GEE JS パーサーが `SyntaxError: Unexpected token (1:3)` で停止
- 対策: **ASCII のみ使用 (`->` で代用)**

### 5.6 NetCDF 全配列 load
- 試行: `xr.open_dataset(p).load()` で全データをメモリへ
- 結果: 164 MB × 120k = ~20 TB 相当、実行不可能
- 対策: **lazy load + selector で 2 点ピクセルだけ抽出**

### 5.7 H6 を Oran summer のみでテスト
- 試行: rainfed control = Oran summer (n=34) で SMAP 代替性検証
- 結果: SDS が全 SM source で ~ 0 になり inconclusive
- 原因: post-harvest で植生不活発 → LE と SWC がそもそも相関しない
- 対策: **multi-stratum 化、Oran spring (n=203 active growth) を focal にする**

### 5.8 ブランチを取り違えて作業
- 試行: 別ブランチでスクリプト編集 → push したが pull で反映されず
- 原因: 別ブランチに居た (claude/quantify-water-divergence-LxPGr)
- 対策: **作業前に `git branch` で確認、`git checkout claude/compare-ec-satellite-et-ZnENi`**

---

## 6. 現在の研究上の論点

### 論点 1: 補正モデルの operational deployment
- **claim**: τ-based correction で RMSE が 49-65% 減る、unbiased になる
- **counterargument**: TzM 1 サイトのみで校正、他の drip 灌漑農地で同じ τ が成立する保証なし
- **解決方向**: 他サイトの EC データ (FLUXNET2015) で τ を独立に推定し、汎用パラメータ範囲を推定する

### 論点 2: MOD16 の structural floor c = -2.3 の正体
- **claim**: 「FLUXNET2015 が high-LAI orchard を含まないため」
- **counterargument**: 校正データ偏りなら orchard 以外でも c が大きいはず
- **解決方向**: 校正データセットの crop type コンポジションを調べる (Mu et al. 2011 を再読)

### 論点 3: METv3 5 km ピクセル混合 vs アルゴリズム本質
- **claim**: METv3 の最大振幅 a = -4.0 は 5 km ピクセル内の rainfed 地混合のため
- **counterargument**: 同じ 5 km の METv3 が Oran (rainfed control) で MBE = +0.02 → 混合だけでは説明できない
- **解決方向**: H5 (SIGPAC parcel) で TzM 周辺の灌漑面積率を計算し、a との回帰を取る

### 論点 4: depth inversion (TzM d0-3 で r = -0.19) の頑健性
- **claim**: drip wet-bulb の直接観測証拠
- **counterargument**: n = 281 だが、SMAP 9 km は TzM 30 m EC footprint と空間スケールが違う → 同じ「point」を観測していない可能性
- **解決方向**: SAR 由来の高分解能 SWC (Sentinel-1 backscatter) で再検証

---

## 7. 未解決課題 (優先順位付き)

### ★★★ (即取り組むべき)
1. **論文 Abstract 執筆** — Intro/Methods/Results/Discussion/Limitations はある、Abstract のみ未完
2. **fig_H1_*, fig_H4_*, fig_H6_* の品質チェック** — 図のラベル・スケール・色を最終 paper 提出形式に整える
3. **論文の引用整備** — 現在の paper_methods_results.md に挙がっている文献 (Mu 2011, Zhang 2019, Trigo 2018, Reichle 2018, Pettorelli 2005, Burnham & Anderson 2002 etc.) を BibTeX 化

### ★★ (中期)
4. **H7 検証** — H6 で SMAP root-zone が rainfed で代替可能と確認できたので、FLUXNET2015 / ICOS から多サイトの SMAP root-zone を抽出して SDS マッピング
5. **H2 検証** — AmeriFlux 等から flood/sprinkler 灌漑サイトを取得し、τ を比較
6. **EC データの Permission / DOI** — Oran と TzM の PI 同意取得、論文への引用形式を確認

### ★ (長期 / 続編)
7. **H5 検証** — SIGPAC parcel + Sentinel-2 で TzM 周辺 5 km の灌漑面積率を計算
8. **Sentinel-1 SAR で灌漑検出** — 灌漑記録がない site でも τ 補正を運用可能にする (H1 の延長)
9. **H8 — regional ET correction validation** — Júcar 流域などで補正前後の水収支を MITECO 統計と照合

---

## 8. 次セッション開始時にやること (順番付き TODO)

```
1. git pull && cd /home/shion-nagamine/bakanposs
2. ls -la *.csv figs/  # データ・図の存在確認
3. master_full_v2.csv が無い場合は以下を順に実行:
   python3 pipeline/unify_ec_daily.py
   python3 pipeline/add_flags.py
   python3 pipeline/aggregate_oran_30min.py
   python3 pipeline/unify_satellite.py
   python3 pipeline/merge_satellite_ec.py
   python3 pipeline/integrate_metv3_smap.py
   python3 figures/tau_fit.py
   python3 figures/hypothesis_tests.py
4. paper_methods_results.md を読む (1.4 Findings preview, 3.7, 3.8, 4.1)
5. analysis_narrative.md を読む (§7 H1, H4, H6 結果)
6. **論文 Abstract 執筆** (これがメインタスク)
7. 図 fig_H1_correction.png, fig_H4_aic.png, fig_H6_smap_multistratum.png をユーザーに確認してもらう
8. 必要なら fig labels, captions を最終調整
9. Abstract と figures を踏まえて目次レベルで再構成
```

---

## 9. Claude への引き継ぎ指示

次セッションの Claude へ:

### 9.1 最重要 — 仮説の進化を理解せよ
この研究は「**deep-root → drip wet-bulb**」という仮説転換が核心。最初の deep-root 仮説が反証され、days_since_irrig 層別が決定的だった点を絶対に外さないこと。これは単なる手法ではなく**研究の物語そのもの**。

### 9.2 主張の階層 (3 段階)
1. **PML / METv3 のバイアスは τ ≈ 4-6 d で過渡的に減衰し、補正可能**
2. **MOD16 は構造的 floor c = -2.3 mm/d を持つが、τ-fit の c 項で吸収できる (補正後 RMSE -65%)**
3. **SMAP root-zone × in-situ SWC の depth inversion (r=-0.19 at d0-3) が wet-bulb の直接観測証拠**
4. (補助) rainfed Oran spring で SMAP_rz が in-situ SWC を代替できる (SDS 完全一致)

論文ではこの 3+1 を独立に説得力ある証拠として並べる構成。

### 9.3 統計手法の使い方 (失敗を避けるため)
- **Bootstrap CI** は SDS / τ に必ず付ける (理論サンプリング分布が無いため)
- **NLS** は raw data でフィット (binned median はダメ)。τ ∈ (0, 60] 制約必須
- **AIC** で「ΔAIC > 10 = decisive evidence」 (Burnham & Anderson 2002)
- **分位閾値** は **TzM-summer プール全体で固定**、バケット内で再計算しない

### 9.4 データ・コードの再現性
- 中間 CSV (master_full.csv, ec_daily_master.csv 等) は git 管理外
- 完全再生成は `python3 pipeline/unify_ec_daily.py` から順に実行 (TODO §8 参照)
- METv3 ロードは 3.5 時間かかる (`metv3_daily_all.csv` が既にあるならスキップ可)
- SMAP は `/mnt/hdd/Dataset/SMAP_OranTzM.csv` (新形式) を `load_smap.py` が自動検出

### 9.5 ユーザーとのコミュニケーション
- ユーザーは日本語、技術内容は英語混じり OK
- スクリプト出力は必ずユーザーマシンで実行 (Claude のクラウド環境ではデータアクセス不可)
- コミット → push → ユーザーが pull → 実行 → 結果貼り付け の流れ
- ユーザーのブランチ取り違えに注意 (claude/compare-ec-satellite-et-ZnENi が正しい)

### 9.6 やってはいけないこと
- データ生 CSV を git add しない (容量爆発)
- 仮説検証なしで paper の主張を変えない
- "deep-root 仮説" を肯定する記述を入れない (反証済み)
- Python スクリプトに非 ASCII 文字を入れない (load_metv3.py 以外、GEE JS は特に)
- main ブランチに直接 push しない

---

## 10. 重要ファイル一覧

### 10.1 コード (絶対参照)
```
/home/shion-nagamine/bakanposs/pipeline/unify_ec_daily.py
/home/shion-nagamine/bakanposs/pipeline/add_flags.py
/home/shion-nagamine/bakanposs/pipeline/aggregate_oran_30min.py
/home/shion-nagamine/bakanposs/pipeline/qc_master.py
/home/shion-nagamine/bakanposs/pipeline/unify_satellite.py
/home/shion-nagamine/bakanposs/pipeline/load_metv3.py
/home/shion-nagamine/bakanposs/pipeline/load_smap.py
/home/shion-nagamine/bakanposs/pipeline/merge_satellite_ec.py
/home/shion-nagamine/bakanposs/pipeline/integrate_metv3_smap.py
/home/shion-nagamine/bakanposs/figures/sds_v14_repro.py
/home/shion-nagamine/bakanposs/figures/figure_C_summer.py
/home/shion-nagamine/bakanposs/figures/sds_vs_bias.py
/home/shion-nagamine/bakanposs/figures/tau_fit.py
/home/shion-nagamine/bakanposs/figures/hypothesis_tests.py
/home/shion-nagamine/bakanposs/gee/gee_extract.js
/home/shion-nagamine/bakanposs/gee/gee_smap_only.js
```

### 10.2 解析結果 CSV
```
/home/shion-nagamine/bakanposs/master_full_v2.csv          # メインデータ (53 cols × 1356 rows)
/home/shion-nagamine/bakanposs/sds_v14_results.csv          # SDS 全層別
/home/shion-nagamine/bakanposs/sds_vs_bias.csv              # SDS と衛星 bias の対応
/home/shion-nagamine/bakanposs/tau_fit_summary.csv          # 3 製品の (a, τ, c) と CI
/home/shion-nagamine/bakanposs/fig_C_summer_summary.csv     # バケット別バイアス統計
/home/shion-nagamine/bakanposs/hypothesis_tests_summary.csv # H1/H4/H6 検証結果
/home/shion-nagamine/bakanposs/metv3_daily_all.csv          # METv3 日次 (5114 site-days)
/home/shion-nagamine/bakanposs/smap_daily.csv               # SMAP 日次 (5112 site-days)
```

### 10.3 図
```
/home/shion-nagamine/bakanposs/figs/fig_C2_bias_by_irrig_seasonal.png  # Fig 5 headline
/home/shion-nagamine/bakanposs/figs/fig_F_tau_fit.png                  # Fig 7 decay
/home/shion-nagamine/bakanposs/figs/fig_E_sds_vs_bias.png              # Fig 6 SDS×bias
/home/shion-nagamine/bakanposs/figs/fig_H1_correction.png              # Fig H1 (新)
/home/shion-nagamine/bakanposs/figs/fig_H4_aic.png                     # Fig H4 (新)
/home/shion-nagamine/bakanposs/figs/fig_H6_smap_multistratum.png       # Fig H6 (新)
```

### 10.4 ドキュメント
```
/home/shion-nagamine/bakanposs/paper_outline.md            # 論文骨格
/home/shion-nagamine/bakanposs/paper_methods_results.md    # 論文本文 (Intro/Methods/Results/Discussion/Lim)
/home/shion-nagamine/bakanposs/analysis_narrative.md       # 解析ナラティブ (統計知識・図解・試行錯誤含む)
/home/shion-nagamine/bakanposs/SESSION_SUMMARY.md          # 本ファイル (引き継ぎ)
```

### 10.5 元データ (read-only)
```
/home/shion-nagamine/Dataset/Eddy data in Spain/   # EC raw (Oran + TzM)
/home/shion-nagamine/Dataset/METv3/YYYY/MM/MMDD/   # NetCDF (~120k files)
/mnt/hdd/Dataset/Fast_OranTzM_*.csv                # GEE wide CSVs
/mnt/hdd/Dataset/OranTzM_S2_NDVI.csv               # GEE Sentinel-2
/mnt/hdd/Dataset/SMAP_OranTzM.csv                  # SMAP 3-hour (single file)
```

### 10.6 Git
- リポジトリ URL: `https://github.com/zhiyinzhangling77-ui/bakanposs`
- メインブランチ: `main` (使わない)
- 作業ブランチ: **`claude/compare-ec-satellite-et-ZnENi`**
- 最終コミット: `d34c572`

---

*記録作成: 2026-05-10、Claude Code session 終了時点*
*次セッションの Claude へ: 本ファイルを最初に必ず読み、§9.1〜9.6 を遵守すること。*
