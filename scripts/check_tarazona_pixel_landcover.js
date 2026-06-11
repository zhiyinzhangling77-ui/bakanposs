/**
 * scripts/check_tarazona_pixel_landcover.js
 * =============================================================
 *
 *  Tarazona EC tower 周辺の土地被覆チェック
 *  目的: H26 (SCATSAR-SWI 1 km) を fetch する価値があるかを判定する。
 *
 *  Tarazona orchard サイズ = 1 ha (0.01 km^2).
 *  各衛星 SM product の pixel に占める orchard 比:
 *    - H26 (1 km)          : 1 %     → pixel 平均は orchard 信号を見やすい
 *                                       「条件」を要する
 *    - ETv3 / MGPP (~5 km) : 0.04 %
 *    - ASCAT H141 (12.5 km): 0.006 %
 *
 *  → 1 km pixel に「灌漑活発な植生」がどれくらい入っているかで
 *    H26 fetch の費用対効果が決まる。
 *
 *  本スクリプトが返すもの:
 *    (1) 5 km tile の summer (Jun-Sep) 平均 NDVI ヒートマップ
 *    (2) 各バッファ (500 m / 1 km / 5 km / 12.5 km) ごとの:
 *        - 平均 summer NDVI
 *        - NDVI > 0.5 (active vegetation) のピクセル比率
 *        - ESA WorldCover 2021 のクラス頻度
 *    (3) NDVI コンポジットを GeoTIFF で Drive にエクスポート (任意)
 *
 *  判定の目安:
 *    - 1 km バッファで NDVI > 0.5 比率が 20% 以上    → H26 fetch 推奨
 *    - 1 km バッファで NDVI > 0.5 比率が 5–20%       → H26 fetch 価値あり
 *    - 1 km バッファで NDVI > 0.5 比率が 5% 未満     → H26 でも信号は
 *                                                        noise 以下、
 *                                                        candidate C 推奨
 *
 *  使い方:
 *    GEE Code Editor (https://code.earthengine.google.com/) に
 *    本スクリプト全文を貼り付けて Run.  Console に統計、Map に
 *    NDVI ヒートマップ + バッファ境界が出る。
 *
 *  Required permissions: 標準の GEE 公開データセットのみ使用。
 *                        個別認可は不要。
 */


// =============================================================
// 0. Site + buffers
// =============================================================

var TARAZONA = ee.Geometry.Point([-1.9397, 39.266]);
Map.centerObject(TARAZONA, 13);

// バッファ半径 = 各衛星 SM product の "半 pixel" に合わせる
// (= square pixel の中心からの距離.  buffer は circular だが
//  視覚的・統計的なスケール比較には十分)
//
// Tower footprint は典型的に 100-300 m (EC stability + wind 依存)
// → 200 m を採用.  以前 500 m にしていたが H26 と被って同じ統計が
//   出ていたバグを修正.
var buffers = {
  '200 m  (tower footprint)'   : TARAZONA.buffer(200),
  '1 km   (~H26 pixel)'        : TARAZONA.buffer(500),
  '5 km   (~ETv3/MGPP pixel)'  : TARAZONA.buffer(2500),
  '12.5 km (~ASCAT H141 pixel)': TARAZONA.buffer(6250)
};

// 表示用. 一番外側の 12.5 km バッファ全域で NDVI/LC を切り出す.
var roi_outer = buffers['12.5 km (~ASCAT H141 pixel)'];


// =============================================================
// 1. Sentinel-2 summer NDVI composite (2020-2024)
// =============================================================

// Tarazona EC tower の運用期間 (2020-2024) 全 5 年の Jun-Sep を平均.
// drip 灌漑の効果が NDVI に乗るとすれば summer.
var s2_summer = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi_outer)
  .filterDate('2020-06-01', '2024-09-30')
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 25))
  .map(function(img) {
    var ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
    return ndvi.copyProperties(img, ['system:time_start']);
  });

var ndvi_summer = s2_summer.mean().rename('NDVI_summer');

// "active summer vegetation" の binary mask (NDVI > 0.5)
//   - rainfed Mediterranean は summer に senescent (NDVI < 0.3)
//   - drip irrigated orchard / alfalfa などは summer も NDVI > 0.5
//   - つまり summer NDVI > 0.5 = "灌漑農地っぽい" の geographic proxy
var active_summer = ndvi_summer.gt(0.5).rename('active');


// =============================================================
// 2. ESA WorldCover 2021 (10 m)
// =============================================================

var worldcover = ee.ImageCollection('ESA/WorldCover/v200').first();

// WorldCover 主要クラス:
//   10 = Tree cover     (orchard 候補)
//   20 = Shrubland
//   30 = Grassland
//   40 = Cropland       (cereal / vineyard / 若い orchard)
//   50 = Built-up
//   60 = Bare / sparse vegetation
//   80 = Permanent water


// =============================================================
// 3. Map visualization
// =============================================================

Map.setOptions('SATELLITE');

// (a) Summer NDVI heatmap
Map.addLayer(
  ndvi_summer.clip(roi_outer),
  {min: 0.0, max: 0.8,
   palette: ['#f7f4ea', '#eddca8', '#cdd778', '#84b15a',
             '#3e8a3e', '#1a5e1a', '#073507']},
  '(a) Summer NDVI mean (Jun-Sep 2020-2024)'
);

// (b) Active summer vegetation mask
Map.addLayer(
  active_summer.selfMask().clip(roi_outer),
  {palette: ['#1D9E75']},
  '(b) Active in summer  (NDVI > 0.5)',
  false
);

// (c) ESA WorldCover
Map.addLayer(
  worldcover.clip(roi_outer),
  {},
  '(c) ESA WorldCover 2021',
  false
);

// (d) Buffer outlines
var outline = function(geom, color, name) {
  Map.addLayer(
    ee.Image().byte().paint({featureCollection:
      ee.FeatureCollection([ee.Feature(geom)]), color: 1, width: 3}),
    {palette: color},
    name
  );
};
outline(buffers['200 m  (tower footprint)'],    'red',    '200 m buffer');
outline(buffers['1 km   (~H26 pixel)'],         'orange', '1 km buffer (H26)');
outline(buffers['5 km   (~ETv3/MGPP pixel)'],   'blue',   '5 km buffer (ETv3)');
outline(buffers['12.5 km (~ASCAT H141 pixel)'], 'purple', '12.5 km buffer');

// (e) Tower point
Map.addLayer(
  ee.FeatureCollection([ee.Feature(TARAZONA)]),
  {color: 'black'},
  'Tarazona EC tower'
);


// =============================================================
// 4. Per-buffer statistics  → Console
// =============================================================

print('============================================================');
print('Tarazona pixel land-cover diagnostic');
print('  EC tower:  lat 39.266,  lon -1.9397');
print('  Orchard:   1 ha (0.01 km^2)');
print('============================================================');

Object.keys(buffers).forEach(function(label) {
  var roi = buffers[label];

  // mean summer NDVI
  var mean_ndvi = ndvi_summer.reduceRegion({
    reducer:    ee.Reducer.mean(),
    geometry:   roi,
    scale:      10,
    maxPixels:  1e9
  });

  // fraction NDVI > 0.5
  var frac_active = active_summer.reduceRegion({
    reducer:    ee.Reducer.mean(),
    geometry:   roi,
    scale:      10,
    maxPixels:  1e9
  });

  // WorldCover frequency
  var lc_histo = worldcover.reduceRegion({
    reducer:    ee.Reducer.frequencyHistogram(),
    geometry:   roi,
    scale:      10,
    maxPixels:  1e9
  });

  print('--- ' + label + ' ---');
  print('  mean summer NDVI    :', mean_ndvi);
  print('  fraction NDVI > 0.5 :', frac_active);
  print('  ESA WorldCover hist :', lc_histo);
});


// =============================================================
// 5. (Optional) Export NDVI composite for paper figure
// =============================================================
//
//  Drive にエクスポートしたい場合は下記コメントアウト解除:
//
// Export.image.toDrive({
//   image:        ndvi_summer.clip(buffers['5 km   (~ETv3/MGPP pixel)']),
//   description:  'Tarazona_summer_NDVI_5km_2020_2024',
//   folder:       'GEE_exports',
//   fileNamePrefix: 'Tarazona_summer_NDVI_5km',
//   scale:        10,
//   region:       buffers['5 km   (~ETv3/MGPP pixel)'],
//   maxPixels:    1e9
// });


// =============================================================
// 6. Decision rule  → 出力された数字を見てこれを当てる
// =============================================================
//
//  Console で 1 km buffer の "fraction NDVI > 0.5" を見る:
//
//   > 0.20  → 周辺に灌漑農地が集積.  H26 fetch を強く推奨.
//                pixel-mean SM 上昇が detectable.
//   0.05 – 0.20 → 一部に灌漑あり.  H26 fetch 価値あり (marginal).
//                  noise floor 近辺の解析になる.
//   < 0.05  → 周囲は rainfed / 乾燥地.  H26 でも検出できない.
//             → candidate C (衛星 SM を外して amp 比較 + tower SWC) を採用.
//
//  併せて WorldCover の class 40 (Cropland) vs 60 (Bare) の比率も
//  チェックすると、Sentinel-2 NDVI が薄雲で過小評価された場合の
//  二重チェックになる.
// =============================================================
