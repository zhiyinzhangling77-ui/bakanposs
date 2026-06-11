// GEE: SMAP L4 3-hourly soil moisture extraction at Oran and TzM.
// Uses 6 km buffer (> half of SMAP 9 km pixel) to ensure pixel capture.
// Export -> Google Drive -> SMAP_OranTzM.csv
//
// How to run:
//   1. Open code.earthengine.google.com
//   2. Paste this script, click Run
//   3. In Tasks tab, click Run next to "SMAP_OranTzM"
//   4. Download CSV from Google Drive

// --- Sites with large enough buffer for SMAP 9 km grid ---
var sites = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([-1.86, 38.82]).buffer(6000),
             {name: 'Oran'}),
  ee.Feature(ee.Geometry.Point([-1.9397, 39.2660]).buffer(6000),
             {name: 'TzM'})
]);

var START = '2018-01-01';
var END   = '2024-12-31';

// SMAP L4 Global 3-hourly 9 km
var smap = ee.ImageCollection('NASA/SMAP/SPL4SMGP/007')
  .filterDate(START, END)
  .select(['sm_surface', 'sm_rootzone']);

// Extract mean over 6 km buffer at each 3-hour timestep
var features = smap.map(function(img) {
  var t = img.date().format('YYYY-MM-dd HH:mm');
  var stats = img.reduceRegions({
    collection: sites,
    reducer: ee.Reducer.mean(),
    scale: 11000   // SMAP native ~9 km; use 11 km to match GEE tile
  });
  return stats.map(function(f) {
    return f.set({date_time: t});
  });
}).flatten();

print('sample', features.limit(5));

Export.table.toDrive({
  collection: features,
  description: 'SMAP_OranTzM',
  fileFormat: 'CSV',
  selectors: ['name', 'date_time', 'sm_surface', 'sm_rootzone']
});
