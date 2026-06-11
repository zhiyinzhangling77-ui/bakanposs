/**
 * scripts/check_both_sites_pixel_landcover.js
 * =============================================================
 *
 *  ★ Fig 1 SOURCE — Side-by-side pixel landscape diagnostic
 *
 *  ポスター/論文 Fig 1 のための GEE エクスポータ.
 *  両 EC tower (Oran rainfed cereal + Tarazona drip almond) の
 *  周辺 5 km の summer NDVI composite + ESA WorldCover を取得し、
 *  Python 側で 2-panel "pixel landscape comparison" 図に組む.
 *
 *  Narrative:
 *    - Oran     : uniform rainfed agricultural matrix
 *                 → satellite pixel content matches the site
 *    - Tarazona : isolated 1 ha orchard in a brown rainfed sea
 *                 → satellite pixel diluted, irrigation invisible
 *
 *  Period per site (matches EC observation period — site-specific):
 *    - Oran     : Jun–Sep  2018–2020  (3 yrs, cereal cycle 12 mo)
 *    - Tarazona : Jun–Sep  2020–2024  (5 yrs, almond active 9 mo)
 *
 *  Buffer scales (matching satellite SM/ET product pixel sizes):
 *    - 200 m   (tower footprint)
 *    - 1 km    (~H26 SCATSAR-SWI)
 *    - 5 km    (~Meteosat ETv3 / MGPP)
 *    - 12.5 km (~ASCAT H141/H142)
 *
 *  Outputs:
 *    1) Map layers per site (NDVI heatmap + WorldCover + buffers)
 *    2) Console: per-site × per-buffer NDVI stats + LC histogram
 *    3) Drive export: 5 km × 5 km NDVI GeoTIFF + WorldCover GeoTIFF
 *       per site, common color scale, ready for Python composition.
 *
 *  Decision rule (1 km buffer, fraction NDVI > 0.5):
 *    > 0.20  → pixel content uniform, satellite can resolve site
 *    0.05–0.20 → mixed, satellite signal diluted (marginal)
 *    < 0.05  → isolated patch, satellite blind (Tarazona: 0.045)
 *
 *  Confirmed earlier:
 *    Tarazona 1 km buffer: NDVI mean 0.31, fraction>0.5 = 4.5 %
 *    Tarazona pixel-mean ΔSM during irrigation ≈ 0.003 m³/m³
 *    → below noise floor → satellite SM cannot see drip events
 *
 *  Usage:
 *    GEE Code Editor → paste full script → Run.
 *    For Drive export: open "Tasks" tab → click Run on each export.
 *    Exported GeoTIFFs land in Drive folder "GEE_exports/".
 */


// =============================================================
// 0. Site definitions
// =============================================================

var SITES = {
  'Oran': {
    point      : ee.Geometry.Point([-1.86, 38.82]),
    date_start : '2018-06-01',
    date_end   : '2020-09-30',
    label      : 'Oran  (rainfed cereal, 2018–2020)',
    color      : 'orange'
  },
  'Tarazona': {
    point      : ee.Geometry.Point([-1.9397, 39.266]),
    date_start : '2020-06-01',
    date_end   : '2024-09-30',
    label      : 'Tarazona  (drip almond, 2020–2024)',
    color      : '#1D9E75'
  }
};

// buffer radius matches half-pixel of each product
var BUFFERS = {
  '200 m  (tower footprint)'   : 200,
  '1 km   (~H26 pixel)'        : 500,
  '5 km   (~ETv3/MGPP pixel)'  : 2500,
  '12.5 km (~ASCAT H141 pixel)': 6250
};

// Drive export tile = 5 km × 5 km centred on each tower
var EXPORT_BUFFER_M = 2500;


// =============================================================
// 1. Helpers
// =============================================================

function summer_ndvi(site_info) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(site_info.point.buffer(EXPORT_BUFFER_M + 500))
    .filterDate(site_info.date_start, site_info.date_end)
    .filter(ee.Filter.calendarRange(6, 9, 'month'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 25))
    .map(function(img) {
      return img.normalizedDifference(['B8', 'B4']).rename('NDVI')
                .copyProperties(img, ['system:time_start']);
    })
    .mean()
    .rename('NDVI_summer');
}


function summarize_site(name, info, ndvi, worldcover) {
  print('============================================================');
  print(info.label);
  print('  Sentinel-2 SR Harmonized, ' +
         info.date_start + '  →  ' + info.date_end);
  print('============================================================');

  var active = ndvi.gt(0.5).rename('active');
  Object.keys(BUFFERS).forEach(function(blabel) {
    var roi = info.point.buffer(BUFFERS[blabel]);
    var stats_ndvi = ndvi.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: roi,
      scale: 10, maxPixels: 1e9
    });
    var stats_active = active.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: roi,
      scale: 10, maxPixels: 1e9
    });
    var stats_lc = worldcover.reduceRegion({
      reducer: ee.Reducer.frequencyHistogram(), geometry: roi,
      scale: 10, maxPixels: 1e9
    });
    print('--- ' + blabel + ' ---');
    print('  mean summer NDVI    :', stats_ndvi);
    print('  fraction NDVI > 0.5 :', stats_active);
    print('  ESA WorldCover hist :', stats_lc);
  });
}


// =============================================================
// 2. Compute NDVI composites + WorldCover
// =============================================================

var WORLDCOVER = ee.ImageCollection('ESA/WorldCover/v200').first();
var NDVI = {};
Object.keys(SITES).forEach(function(name) {
  NDVI[name] = summer_ndvi(SITES[name]);
});


// =============================================================
// 3. Map visualisation (default centered on Tarazona — toggle in
//    Layer Manager to see Oran).
// =============================================================

Map.setOptions('SATELLITE');
Map.centerObject(SITES['Tarazona'].point, 13);

var NDVI_VIS = {
  min: 0.0, max: 0.8,
  palette: ['#f7f4ea', '#eddca8', '#cdd778', '#84b15a',
              '#3e8a3e', '#1a5e1a', '#073507']
};


function outline(geom, color, layer_name) {
  Map.addLayer(
    ee.Image().byte().paint({
      featureCollection: ee.FeatureCollection([ee.Feature(geom)]),
      color: 1, width: 3
    }),
    {palette: color},
    layer_name
  );
}


Object.keys(SITES).forEach(function(name) {
  var info = SITES[name];
  var roi_outer = info.point.buffer(EXPORT_BUFFER_M + 500);

  // (a) summer NDVI heatmap — Tarazona visible by default
  Map.addLayer(
    NDVI[name].clip(roi_outer), NDVI_VIS,
    '[' + name + ']  summer NDVI mean', name === 'Tarazona'
  );

  // (b) "active summer" mask (off by default)
  Map.addLayer(
    NDVI[name].gt(0.5).selfMask().clip(roi_outer),
    {palette: [SITES[name].color]},
    '[' + name + ']  active in summer (NDVI > 0.5)',
    false
  );

  // (c) ESA WorldCover (off by default)
  Map.addLayer(
    WORLDCOVER.clip(roi_outer), {},
    '[' + name + ']  ESA WorldCover 2021', false
  );

  // (d) buffer outlines
  Object.keys(BUFFERS).forEach(function(blabel) {
    outline(info.point.buffer(BUFFERS[blabel]),
              info.color,
              '[' + name + ']  buffer ' + blabel);
  });

  // (e) tower point
  Map.addLayer(
    ee.FeatureCollection([ee.Feature(info.point)]),
    {color: 'black'},
    '[' + name + ']  EC tower'
  );
});


// =============================================================
// 4. Console statistics (per site × per buffer)
// =============================================================

Object.keys(SITES).forEach(function(name) {
  summarize_site(name, SITES[name], NDVI[name], WORLDCOVER);
});


// =============================================================
// 5. Drive exports — for Python composition into Fig 1
// =============================================================
//
// Output filenames (all in Drive folder "GEE_exports/"):
//
//   pixel_landscape_Oran_summerNDVI_5km.tif
//   pixel_landscape_Tarazona_summerNDVI_5km.tif
//   pixel_landscape_Oran_worldcover_5km.tif
//   pixel_landscape_Tarazona_worldcover_5km.tif
//
// All centered on the EC tower, 5 km × 5 km extent, 10 m pixel.
// Same color scale (NDVI 0–0.8) → ready for direct Python composition.

Object.keys(SITES).forEach(function(name) {
  var info = SITES[name];
  var roi  = info.point.buffer(EXPORT_BUFFER_M);

  Export.image.toDrive({
    image:          NDVI[name].clip(roi),
    description:    'NDVI_5km_' + name,
    folder:         'GEE_exports',
    fileNamePrefix: 'pixel_landscape_' + name + '_summerNDVI_5km',
    scale:          10,
    region:         roi,
    maxPixels:      1e9
  });

  Export.image.toDrive({
    image:          WORLDCOVER.clip(roi),
    description:    'WorldCover_5km_' + name,
    folder:         'GEE_exports',
    fileNamePrefix: 'pixel_landscape_' + name + '_worldcover_5km',
    scale:          10,
    region:         roi,
    maxPixels:      1e9
  });
});


// =============================================================
// 6. Decision rule reminder
// =============================================================
//
//  Look at the 1 km buffer "fraction NDVI > 0.5" in console:
//
//    > 0.20         → pixel is uniform, satellite OK
//    0.05 – 0.20    → mixed, satellite signal diluted
//    < 0.05         → isolated patch, satellite blind
//
//  Expectation:
//    Oran     : ≫ 0.20  (rainfed cereal covers continuously)
//    Tarazona : ≈ 0.05  (already confirmed in earlier run; orchard
//                          is a 1 ha green dot, rest of pixel is
//                          rainfed dryland)
// =============================================================
