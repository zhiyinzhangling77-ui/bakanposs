# `data_loaders.py` 採用ガイド — 解析 A/B 向け移行手順

このドキュメントは解析 A (`claude/quantify-water-divergence-LxPGr`) と
解析 B (`claude/compare-ec-satellite-et-ZnENi`) の作業者向け。
解析 C で発見・修正された **Oran ローダの致命的バグ** を共通モジュール
`data_loaders.py` で吸収しているので、各ブランチに取り込んで利用する。

---

## 0. なぜ移行が必要か（要約）

解析 A の `analysis_A_v9.load_oran_ec` (および v10–v14 で同じだった場合) は、
Oran 半時間値 CSV (52,606 行) のうち **真夜中行 1 つだけ (914/52,606 = 1.7%)
しか日次集計に使っていなかった**。原因:

- pandas の自動フォーマット推定が CSV 先頭行 `2018/01/01` (日付のみ) を
  `%Y/%m/%d` と確定し、残り 51,690 行 `2018/01/01 00:30:00` を全部 NaT 化
- `dropna(subset=['DateTime'])` で 914 行（各日真夜中の 1 行）が残る
- それを「日平均」と称して使う → Rn 中央値が **−63 W/m²** という物理的に
  不可能な値が出ていた

修正版 (`data_loaders.load_oran_ec_clean`):

| 指標 (Oran) | v9 (夜間値のみ) | C clean (正しい半時間→日平均) |
|---|---|---|
| 有効日数 | 914 | **922** |
| Rn 中央値 [W/m²] | −63.14 | **+92.82** |
| LE 中央値 [W/m²] | +2.01 | **+19.85** |
| ET 中央値 [mm/day] | 0.002 | **1.526** |
| EF 有効日数 | 28 | **899** |

→ **A/B の Oran 絶対値・サイト間比較は再評価が必要**

---

## 1. 移行手順（A ブランチで実行する例）

```bash
# 1. C ブランチから data_loaders.py を取り込む
git checkout claude/quantify-water-divergence-LxPGr
git checkout claude/ndvi-flux-analysis-iaS0b -- data_loaders.py
git add data_loaders.py
git commit -m "import data_loaders.py from analysis-C branch (TIMESTAMP/unit fixes)"
```

```bash
# 2. analysis_A_v9.py (or v14) のローダ呼び出しを差し替える
```

差し替え対象は通常以下 2 ヶ所:

```python
# Before:
from analysis_A_v9 import load_oran_ec, load_tarazona_ec
oran_ec = load_oran_ec(PATHS["oran_ec"])
tara_ec = load_tarazona_ec(PATHS["tara_ec"])

# After:
from data_loaders import load_oran_ec_clean, load_tarazona_ec_clean, normalize_swc
oran_ec = load_oran_ec_clean(PATHS["oran_ec"])
tara_ec = normalize_swc(load_tarazona_ec_clean(PATHS["tara_ec"]), "Tarazona")
```

注意点:

- `load_tarazona_ec_clean` は SWC を `m³/m³` で返す可能性がある
  (旧データソースの場合)。Oran と単位を揃えるため `normalize_swc` で
  `0–1` レンジを `%` に変換する。
- 新ローダは VPD を **kPa** に統一する (Oran は元 hPa から `÷10`)
- 新ローダは ET を `mm/day` に統一する (Tarazona は `ET_avg` ではなく `ET_sum`)

---

## 2. B ブランチも同様

```bash
git checkout claude/compare-ec-satellite-et-ZnENi
git checkout claude/ndvi-flux-analysis-iaS0b -- data_loaders.py
git add data_loaders.py
git commit -m "import data_loaders.py from analysis-C branch"
```

そして B 側のスクリプトでも `load_oran_ec_clean` / `load_tarazona_ec_clean`
を使うよう差し替え。

---

## 3. 移行後の影響

**期待される変化** (絶対値):

- Oran の LE, H, Rn, ET の絶対値が大きく変わる (上の表参照)
- Oran の EF が計算可能日数が 28 → 899 日に増える
- Oran VPD が kPa に統一されるので、サイト間 VPD 比較が初めて公平になる
- Tarazona ET が `mm/day` で正しく出る (旧 `ET_avg` は別単位だった)
- **半時間レコード使用率が 84% → 100%** (year/Julian/Time_hours フォールバックにより
  TIMESTAMP='nan' の 8517 行が復元される)

**期待される変化** (相関・サイト間比較):

- 解析 A の「Tarazona > Oran の蒸散」のサイト間比較は方向性は同じだが、
  Oran 側が修正で増える分、効果サイズは縮む可能性
- 解析 C の H1 結果 (灌漑日数で EF が急減衰) を踏まえると、
  解析 A の「深根アクセス仮説」は灌漑効果に置き換えるべき

---

## 4. テスト方法

`data_loaders.py` を取り込んだら、まず単体動作を確認:

```python
from data_loaders import load_oran_ec_clean, load_tarazona_ec_clean
oran = load_oran_ec_clean("/path/to/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv")
tara = load_tarazona_ec_clean("/path/to/Daily_Summary_Filtered_forPred_ActEne26.csv")

# 期待値
assert 900 <= len(oran) <= 950, f"Oran 日数異常: {len(oran)}"
assert 80 < oran['Rn'].median() < 110, f"Oran Rn med 異常: {oran['Rn'].median()}"
assert 700 <= len(tara) <= 750
assert tara['ET'].median() > 2, "Tarazona ET が小さすぎる (ET_avg を使ってる?)"
print("OK: data_loaders 健全")
```

---

## 5. data_loaders.py の API

```python
load_oran_ec_clean(filepath, verbose=True) -> pd.DataFrame
  columns: ['date', 'SWC', 'LE', 'H', 'G', 'Rn', 'VPD', 'ET', 'EF', 'site']
  units:   W/m² for fluxes, mm/day for ET, kPa for VPD, % for SWC

load_tarazona_ec_clean(filepath, verbose=True) -> pd.DataFrame
  columns: 上記に加えて Irrig_mm, Rain_mm, IrrigRain_mm, GPP_avg,
           NDVI_orig, NDVI_interp, Fv_interp (CSV にあれば)
  units:   同上

normalize_swc(df, site_name) -> pd.DataFrame
  SWC が 0-1 (m³/m³) なら % に変換

定数: EF_DENOM_MIN = 10.0, SENTINEL_THR = -9000.0, VPD_MAX_KPA = 10.0
```

---

## 6. 参照

- 解析 C レポート: `reports/analysis_C_report.md` §4 (試行錯誤のタイムライン)
- 解析 C 主スクリプト: `analysis_C_v1.py` (使用例)
- 解析 C ブランチ: `claude/ndvi-flux-analysis-iaS0b`
