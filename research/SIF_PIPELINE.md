# SIF 取得〜検証パイプライン（旗38）— 呼吸残差の記憶を独立シグナルで破る

**目的**：旗37 で「呼吸残差に約4日の記憶があり、気象では作れない／だが分割窓(4-7日)と交絡＝生物か
アルゴリズムか切れない」と分かった。**SIF は分割に依らない独立な光合成シグナル**なので、SIF が残差の
記憶を説明すれば＝基質供給が本物と証明でき、交絡を破れる（旗38 の検証エンジンで判定）。

## なぜ SIF（他の衛星でなく）か
旗37 で記憶の時間スケール＝**約4日**（e-folding 3〜7日）と特定。これは文献の**基質供給ラグ(1〜5日,
Stoy2007/Migliavacca2011)** と一致し、長い水分記憶(2〜10週, Cable2013)ではない。→ 狙う衛星は
**SMAP(週)でなく SIF(基質＝日〜数日スケール)**。時間スケールが product を決めた。

## プロダクト選定
| SIF | 解像度 | 期間 | 独立性 | 判断 |
|---|---|---|---|---|
| **CSIF**(Zhang) | 0.05°, **4-day** | 2000– | 中（MODIS 反射率で再構成＝完全独立でない）| **本命**：4日刻みが旗37の≈4日と一致し、フラックス長期記録を覆う |
| GOSIF(Li&Xiao) | 0.05°, 8-day | 2000– | 中（同上, 再構成）| 代替（8日はやや粗い）|
| TROPOMI SIF | ~7km, 日 | 2018– | 高（真の SIF）| 検証用：CSIF の結果を 2018+ で裏取り |

**方針**：主解析は CSIF(4-day)、頑健性確認に TROPOMI(2018+)。**再構成 SIF は「基質」か「緑度/フェノロジー」の
代理か曖昧**（LITERATURE_NOTES③の未解決点）なので、結論は控えめに＝「SIF が説明する」まで、断定は TROPOMI で。

## 手順

### 1. タワー座標を用意（前提）
```bash
python research/sif_coords.py --sites JP-Tak JP-Ta2 JP-Mse JP-BBY CN-HaM MN-Hst MN-Kbu MN-Nkh \
  JP-Fhk JP-Fjy RU-Ege TH-Mae --out site_coords.csv
```
BADM があれば lat/lon が入る。**空欄のサイトは AsiaFlux/FLUXNET のサイト情報から手入力**（座標は一次情報から、捏造不可）。

### 2. SIF をタワー画素で抽出

**推奨：GeoTIFF から抽出（GEE 認証不要・取得元非依存）** `sif_extract_geotiff.py`。
GOSIF/CSIF を公式配布から直接DL、または GEE でエクスポートした GeoTIFF、どれでも動く。
```bash
# 例: GOSIF 8-day（ファイル名 GOSIF_2018001.tif=年+通日, scale/nodata は要確認）
python research/sif_extract_geotiff.py --coords site_coords.csv --tifdir /path/to/GOSIF \
    --date-regex "(\d{4})(\d{3})" --date-fmt yyyyddd --scale 0.0001 --nodata 32767
# 例: CSIF 4-day（ファイル名に YYYYMMDD）
python research/sif_extract_geotiff.py --coords site_coords.csv --tifdir /path/to/CSIF \
    --date-regex "(\d{8})" --date-fmt yyyymmdd
```
出力：各サイト `<site>_sif.csv`（列 date, sif）。**--scale/--nodata は各プロダクト仕様を要確認**。

**代替：Google Earth Engine（アセットがあれば, ローカルで認証）**
GEE の Python API（`pip install earthengine-api`, `earthengine authenticate`）。CSIF/GOSIF が GEE アセットに
無い場合は (a) 各研究室の配布 GeoTIFF を直接ダウンロードしてタワー画素を抽出、(b) TROPOMI L2 SIF を GEE で。
下は **GOSIF(GEE コミュニティアセット例) or 任意の SIF ImageCollection** に対する雛形（アセット ID は各自確認）：

```python
# extract_sif_gee.py （ローカルで実行, 要 ee 認証）
import ee, pandas as pd
ee.Initialize()
coords = pd.read_csv("site_coords.csv").dropna(subset=["lat", "lon"])
SIF = ee.ImageCollection("<SIF_ASSET_ID>")          # ← CSIF/GOSIF/TROPOMI のアセットIDに置換
BAND = "<SIF_BAND>"                                  # ← 例 "SIF", "b1"
for _, r in coords.iterrows():
    pt = ee.Geometry.Point(float(r.lon), float(r.lat))
    fc = SIF.select(BAND).getRegion(pt, 5000).getInfo()   # 5km 近傍
    df = pd.DataFrame(fc[1:], columns=fc[0])
    df = df.rename(columns={BAND: "sif"})[["time", "sif"]].dropna()
    df["date"] = pd.to_datetime(df["time"], unit="ms")
    df[["date", "sif"]].to_csv(f"{r.site}_sif.csv", index=False)
    print(r.site, len(df))
```
出力：各サイト `<site>_sif.csv`（列 date, sif）。夏(7-8月)を含めば足りる。

### 2-TROPOMI. 真SIF（GOSIF が null だった後の、より公平なテスト）
旗38 で GOSIF(再構成8日)は記憶を説明せず、旗39 で記憶は「生物×分割窓の混在」と判明。**真のSIF(TROPOMI, 実測・
ほぼ日次)は GOSIF の2弱点(再構成・8日)を外す**＝「生物成分がどれだけ残るか」を測る。NetCDF 配布なので専用抽出器：
```bash
# 入手(ローカル): Caltech TROPOMI SIF  ftp://fluo.gps.caltech.edu/data/tropomi/  (2018-03〜2021-07, ほぼ日次)
# まず中身(変数名・座標名)を確認
python research/sif_extract_netcdf.py --ncdir /path/to/tropomi --list
# 抽出（--var は上の一覧から, 既定は名前に sif を含む変数を自動選択）
python research/sif_extract_netcdf.py --coords site_coords.csv --ncdir /path/to/tropomi --var <SIF変数>
```
**留保**：2018+ で夏の重なり~4年＝短いが、**ほぼ日次で4日記憶を分解できる**（GOSIF 8日はできなかった）。
画素~3.5×5.5km で footprint 不一致は残る。RTSIF 等の「再構成TROPOMI」は GOSIF 同様の再構成なので**避ける**（真SIFで）。

**重要（配布形態の実態）**：Caltech の無料公開は **ungridded L2（軌道ごとの sounding 点群, netCDF4）**。gridded は
リクエスト制、gridded TROPOSIF は 0.2°(≒22km)で GOSIF より粗い。**良 footprint の真SIF は ungridded L2 経路のみ**＝
格子でなく点群なので、タワー半径 N km 内の sounding を日次集約する専用抽出器 `sif_extract_ungridded.py` を使う：
```bash
# 入手(ローカル): Caltech ungridded TROPOMI SIF  ftp://fluo.gps.caltech.edu/data/tropomi/ (L2, 2018-03〜2021-07)
python research/sif_extract_ungridded.py --ncdir /path/to/L2 --list                 # 変数名を確認
python research/sif_extract_ungridded.py --coords site_coords.csv --ncdir /path/to/L2 \
    --var sif --lat-var lat --lon-var lon --time-var TIME --radius 20               # 半径20kmで日次集約
```
`sif_extract_netcdf.py`(格子用)は gridded 品にのみ有効。ungridded L2 は必ず `sif_extract_ungridded.py`。

### 3. 検証：SIF は残差の記憶を説明するか（旗38）
```bash
python research/sif_respiration_step38.py --site JP-Tak --sif JP-Tak_sif.csv --qc-max 1
```
- **★ SIF で残差 ACF が急落＋偏相関が高い** → 基質供給が記憶の正体＝分割窓アーティファクトでなく生物（交絡を破った）。
- ○ 一部説明 → 基質は寄与するが残りは別の未観測（深水分・微生物）。
- ・ 落ちない → SIF でも説明できず更に深い未観測（or 空間不一致/整列の問題）。

合成検証済み：SIF が隠れ基質を捉える→ACF 0.77→0.19・偏相関+0.94、無関係な対照→変化なし。

## 留保（正直に）
- **空間不一致**：SIF 画素(数km) vs タワー footprint(数百m)。不均一地形で残差（LITERATURE_NOTES③）。複数画素平均/
  ダウンスケールで緩和するが完全でない。
- **再構成 SIF の独立性**：CSIF/GOSIF は MODIS 反射率で再構成＝「基質」でなく「緑度/フェノロジー」を捉えている可能性。
  → 断定は TROPOMI(真の SIF, 2018+)で裏取り。
- **GER は分割派生量**：SIF が説明しても、SIF≈GPP なので「GPP と GER が実は同じ NEE 由来」の循環に注意
  （旗38 は SIF＝独立測定なので循環は回避されるが、SIF が GPP の代理である以上、光合成-呼吸の因果方向は別途要検討）。
- SIF は日中・晴天条件の制約、4/8日→日次補間の平滑。

## この検証が決める分岐
- SIF が記憶を説明 → **「呼吸の未観測駆動＝最近の光合成による基質供給」を独立シグナルで実証**（旗25→37→38 で環が閉じる）。
- SIF でも残らない → 基質でなく更に深い未観測（深土壌水分=SMAP L4、微生物動態）＝次の衛星/観測へ。
