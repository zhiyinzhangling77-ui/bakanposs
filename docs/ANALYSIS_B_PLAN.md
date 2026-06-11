# Analysis B — Satellite ET Validation (PLAN)

> **Status**: ⬜ PLANNED(これから着手)
> **Role**: 解析A の結論を **独立データ(衛星 ET)** で検証
> **Why now**: 解析A は EC タワー単独。同じ τ ≈ 3-4d を独立データで再現できれば、最強の confirmation。再現できなければ、衛星 ET モデルの系統誤差を定量化(別の科学貢献)

---

## 1. 目的と研究全体での位置

### 1.1 解析B が答える問い
- **Q1**: 衛星 ET プロダクトは 解析A の τ ≈ 3-4 d を再現するか?
- **Q2**: 振幅 (LE_0 − LE_∞) の 4× scaling は衛星でも見えるか?
- **Q3**: 衛星 ET ←→ EC ET の bias は灌漑後何日で減衰するか?
- **Q4**: 周辺ピクセル(他作物・他管理)で τ はどう変わるか?

### 1.2 解析A への補強関係
```
解析A (EC タワー、点観測)
   ↓ τ ≈ 3-4d, 振幅 4× scaling
解析B (衛星 ET、空間観測)
   - τ 再現できる → 結論の独立確認 ★
   - τ 再現できない → 衛星モデル系統誤差の発見(別 contribution)
```

### 1.3 論文化への含意
- 解析A 主結果 + 解析B 衛星確認 → "**ground + satellite triangulation**" 強い主張
- もし衛星で bias が見えれば → 「**satellite ET retrieval in irrigated drylands**」の章を追加
- SMAP-based ET モデルへの直接的批評

---

## 2. 使用する衛星 ET プロダクト

| プロダクト | 解像度 | 期間 | 算法 | 推奨度 |
|---|---|---|---|---|
| **MOD16A2** | 500m / 8-day | 2000-現在 | PM-Mu モデル | ★★★ メイン |
| **ECOSTRESS** | 70m / 1-5 day | 2018-現在 | PT-JPL 熱赤外 | ★★★ サブメイン |
| **PML-V2** | 500m / 8-day | 2000-現在 | プロセス + ML | ★★ 補強 |
| **GLEAM** | 0.25° / daily | 1980-現在 | semi-empirical | ★ コース解像度参考 |
| **MODIS LST** (MOD11A1) | 1km / daily | 2000-現在 | 熱赤外 → ΔT 解析 | ★★ 補助 |

### データ取得方法
- **MOD16A2**: AppEEARS, GEE (Google Earth Engine)
- **ECOSTRESS**: AppEEARS, EarthData
- **PML-V2**: NESDC (中国データセンター) or GEE
- **GLEAM**: GLEAM.eu

→ **推奨**: GEE スクリプトで Tarazona/Oran 座標で時系列抽出 → CSV 出力

---

## 3. 設計図(コード設計)

### 3.1 想定ファイル構成
```
analysis_B_v1_satellite_ET.py
  → 衛星 ET 時系列抽出 + 解析A と同じ τ fit
analysis_B_v2_amplitude_check.py
  → 衛星振幅 vs EC 振幅 の照合
analysis_B_v3_bias_decay.py
  → (衛星 ET - EC ET) の時間構造
analysis_B_v4_spatial_extension.py
  → 周辺ピクセルでの τ パターン
```

### 3.2 必須機能(v1 で実装)

#### 入力
- 衛星 ET CSV(日次 or 8-day)
- 解析A の v4 parquet(EC 比較用)

#### 処理
1. 衛星 ET 時系列を読み込み
2. 解析A と同じ event 検出(雨/灌漑)
3. 解析A と同じ exp fit (`LE_inf + (LE_0 - LE_inf) × exp(-d/τ)`)
4. Bootstrap CI
5. EC ET と比較(scatter, time series)

#### 出力
- `v1_satellite_tau_comparison.csv`(衛星 τ vs EC τ)
- `fig01_satellite_recovery.png`(衛星 ET 回復曲線)
- `fig02_satellite_vs_ec.png`(scatter + 時系列)
- `fig03_tau_comparison.png`(衛星 τ vs EC τ pairwise)

### 3.3 重要な処理上の注意

#### 解像度差への対処
- MOD16 (500m) は **アーモンド orchard を1ピクセルでカバーしない**
- 周辺の混合(果樹園 + 裸地 + 他作物)が含まれる
- → ピクセル選定: タワー位置 + 周辺最大 3×3 grid を比較

#### 時間解像度差への対処
- MOD16 は 8-day composite → daily の解析A と直接比較不可
- 解決策:
  - 衛星 ET を daily に線形補間
  - OR 解析A も 8-day 集約版で比較
- ECOSTRESS は不定間隔(1-5 日) → match_asof で時間アライン

#### Bias 補正
- 衛星 ET と EC ET の絶対値はずれる(典型 10-30%)
- → 絶対値比較より **時間構造(τ)、相対変化(amplitude ratio)** を重視

---

## 4. 期待される結果のシナリオ

### Scenario A: 衛星 τ ≈ EC τ ≈ 3-4 d ★★★
**論文記述案**:
> "Independent assessment using MOD16A2 satellite ET yielded τ_MOD16 =
> X.X d (95% CI: ...), statistically indistinguishable from the EC-based
> estimate (Δτ < MDE). This cross-platform consistency provides strong
> support for the universal effective relaxation timescale in this
> Mediterranean system."

→ 解析A 結論が完全確証

### Scenario B: 衛星 τ ≠ EC τ ⚠
**論文記述案**:
> "MOD16A2 satellite ET yielded τ_MOD16 = X.X d, substantially longer/shorter
> than the EC-based estimate (3.36 d). This systematic discrepancy indicates
> that the MOD16 algorithm does not capture the irrigation pulse dynamics
> at the field scale. We attribute this to [LSM-based soil moisture
> parameterization / coarse resolution averaging / etc.]"

→ 衛星モデル系統誤差の発見(別の重要 finding)

### Scenario C: 振幅 scaling 確認 → SMAP モデルへの直接的含意
- 衛星 ET amplitude 比 ≈ 4× なら → "SMAP-driven ET retrievals must scale by ~4× in drip irrigated zones"

---

## 5. データ取得チェックリスト

### Phase 1: MOD16A2 抽出(必須、最初の解析)
- [ ] GEE script で Tarazona (lat 39.266, lon -1.9397) 8-day ET 時系列抽出
- [ ] GEE script で Oran (lat 38.82, lon -1.86) 8-day ET 時系列抽出
- [ ] 期間: 2018-2024(両サイトカバー)
- [ ] 単位確認:MOD16 は kg/m²/8day = mm/8day → 統一(mm/day or W/m²)

### Phase 2: ECOSTRESS 抽出(高解像度補強)
- [ ] EarthData アクセス
- [ ] 70m grid で Tarazona/Oran ピクセル抽出
- [ ] 時間アラインメント(不定間隔)

### Phase 3: 周辺ピクセル(空間拡張)
- [ ] Tarazona から ±5 km 範囲のピクセル
- [ ] 他作物(オリーブ、ぶどう等)タイプ判定
- [ ] LandCover map(Sentinel-2 derived)で作物区別

---

## 6. 解析B コード設計テンプレート(v1 で書くもの)

```python
"""
analysis_B_v1_satellite_ET.py

衛星 ET プロダクトで解析A 結論(τ ≈ 3-4d)を独立検証
"""

# === 1. 衛星 ET 読込 ===
def load_mod16(filepath, site):
    """MOD16A2 8-day CSV を読込み、daily 線形補間"""
    ...

# === 2. EC との時間アライン ===
def align_with_ec(satellite_df, ec_df):
    """衛星 ET と EC ET を共通 date grid で merge"""
    ...

# === 3. 解析A と同じ event detection + fit ===
from analysis_A_v27 import detect_events, fit_with_validation

def satellite_recovery_tau(sat_data, events, le_inf):
    """衛星 ET の event 後 τ を fit (解析A と同じロジック)"""
    ...

# === 4. EC vs Satellite 比較 ===
def compare_tau(tau_ec, tau_sat, ci_ec, ci_sat):
    """MDE-based 比較、CI 重なり check"""
    ...

# === MAIN ===
if __name__ == "__main__":
    # Tarazona
    sat_t = load_mod16("MOD16A2_Tarazona.csv", "Tarazona")
    # EC データは v27 と同じ parquet から
    ...
    # τ fit + comparison
    ...
    # 出力: 衛星 τ vs EC τ
```

### 想定する CLI:
```bash
python analysis_B_v1_satellite_ET.py \
    --mod16-tara MOD16A2_Tarazona.csv \
    --mod16-oran MOD16A2_Oran.csv \
    --ec-parquet daily_classified_v4.parquet \
    --out output_analysis_B_v1
```

---

## 7. データ準備手順(優先順)

### Step 1: GEE script 作成(まず必須)
```javascript
// Earth Engine script
var sites = [
  {name: 'Tarazona', lat: 39.266, lon: -1.9397},
  {name: 'Oran',     lat: 38.82,  lon: -1.86}
];
var collection = ee.ImageCollection("MODIS/061/MOD16A2GF")
    .filterDate('2018-01-01', '2025-01-01')
    .select('ET');
// 各 site で 8-day ET 抽出 → CSV export to Drive
```

### Step 2: ECOSTRESS 取得(可能なら)
- NASA EarthData login が必要
- AppEEARS で point extraction
- 70m 解像度なので Tarazona orchard をピクセルカバー可

### Step 3: PML-V2(中国データ、補助)
- GEE にも上がっている(`projects/google/PML_V2_V017`)
- 同じスクリプトで extract 可能

---

## 8. 統計検証戦略(解析A と同じロジックを再利用)

### 同じ方法論を採用
1. Event detection(Rain > 3mm / Irrig > 0.5mm)
2. Pooled fit with 2-param exp model
3. Bootstrap CI (Method A, 5000回)
4. BCa CI
5. Pairwise MDE comparison
6. Validation filter (R²≥0.3, τ in range, n≥3)

### 解析A と異なる注意点
- 衛星 ET は **絶対値が EC と異なる**(bias あり)
- → 絶対 τ 値の比較は妥当だが、絶対 LE_0 / LE_∞ の比較には注意
- → 振幅比(amp_sat / amp_ec)で management signal を再評価

---

## 9. 期待 timeline と難易度

| Step | 時間 | 難易度 |
|---|---|---|
| GEE script で MOD16 抽出 | 2-3h | 中(初心者なら半日) |
| analysis_B_v1 実装 | 2-3h | 易(解析A 流用) |
| Tarazona/Oran τ 再 fit | 1h | 易 |
| EC vs satellite 比較 | 2h | 中 |
| 結果可視化 | 1h | 易 |
| **合計 v1** | **約 1 日** | |
| ECOSTRESS 追加 | +1日 | 中 |
| 空間拡張(周辺ピクセル) | +1-2日 | 中-高 |

---

## 10. 解析B Success Criteria

### Minimum(これは必達)
- [ ] MOD16A2 で Tarazona/Oran の年単位 ET 時系列取得
- [ ] 解析A と同じロジックで τ_sat を計算
- [ ] τ_sat vs τ_EC 直接比較 + MDE check

### Stretch(あれば論文を強化)
- [ ] ECOSTRESS 70m での τ 再現
- [ ] 周辺ピクセル(他作物)で τ 多サンプル
- [ ] 衛星 ET bias の時間構造(灌漑後何日で平均化される?)

---

## 11. 次セッションで最初にやること

1. 本ファイル `ANALYSIS_B_PLAN.md` を読む
2. `ANALYSIS_A_FINAL.md` で解析A の最終 τ 値を確認
3. `RESEARCH_OVERVIEW.md` で全体図を確認
4. GEE script を作って Tarazona/Oran の MOD16A2 を CSV 化(または既に取得済みなら path 確認)
5. analysis_B_v1.py を書き始める(解析A v27 を base に)

質問あれば直接ユーザに確認:
- MOD16A2 はもう取得済みか? あればパスを教えてもらう
- ECOSTRESS の利用可能性
- GEE アカウントの有無
