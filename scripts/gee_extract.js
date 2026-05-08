// GEE Code Editor script: extract daily/8-day satellite time series at 2 sites.
// Sites: Oran (winter cereal) and TzM (almond), Albacete/Murcia, Spain.
// Output: single CSV to Google Drive with columns
//   date, name, product, mean
// Pivot in pandas after download.

var sites = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([-1.86,    38.82]),  {name: 'Oran', buf: 200, crop: 'cereal'}),
  ee.Feature(ee.Geometry.Point([-1.9397,  39.2660]), {name: 'TzM',  buf: 300, crop: 'almond'})
]);

var regions = sites.map(function (f) {
  return f.buffer(ee.Number(f.get('buf')));
});

var START = '2018-01-01';
var END   = '2024-12-31';

var products = {
  MOD16: {
    ic: 'MODIS/061/MOD16A2GF',
    bands: ['ET', 'PET'],
    scale: 500
  },
  PML: {
    ic: 'CAS/IGSNRR/PML/V2_v018',
    bands: ['Ec', 'Es', 'Ei'],
    scale: 500
  },
  S2: {
    ic: 'COPERNICUS/S2_SR_HARMONIZED',
    bands: ['B3', 'B4', 'B8', 'B11'],
    scale: 10
  },
  LST: {
    ic: 'MODIS/061/MOD11A1',
    bands: ['LST_Day_1km', 'LST_Night_1km'],
    scale: 1000
  },
  LAI: {
    ic: 'MODIS/061/MCD15A3H',
    bands: ['Lai', 'Fpar'],
    scale: 500
  },
  SMAP: {
    ic: 'NASA/SMAP/SPL4SMGP/007',
    bands: ['sm_surface', 'sm_rootzone'],
    scale: 11000
  },
  ERA5: {
    ic: 'ECMWF/ERA5_LAND/HOURLY',
    bands: [
      'temperature_2m',
      'dewpoint_temperature_2m',
      'total_precipitation',
      'surface_net_solar_radiation'
    ],
    scale: 11132
  },
  CHIRPS: {
    ic: 'UCSB-CHG/CHIRPS/DAILY',
    bands: ['precipitation'],
    scale: 5566
  }
};

function extractProduct(p, key) {
  var ic = ee.ImageCollection(p.ic)
    .filterDate(START, END)
    .filterBounds(regions)
    .select(p.bands);

  return ic.map(function (img) {
    var stats = img.reduceRegions({
      collection: regions,
      reducer: ee.Reducer.mean(),
      scale: p.scale
    });
    return stats.map(function (f) {
      return f.set({
        date: img.date().format('YYYY-MM-dd'),
        product: key
      });
    });
  }).flatten();
}

var collections = Object.keys(products).map(function (k) {
  return extractProduct(products[k], k);
});

var all = ee.FeatureCollection(collections).flatten();

Export.table.toDrive({
  collection: all,
  description: 'OranTzM_satellite_timeseries',
  fileFormat: 'CSV',
  selectors: ['date', 'name', 'crop', 'product'].concat(
    // dump every band-mean column the reducer emits
    [].concat.apply([], Object.keys(products).map(function (k) {
      return products[k].bands.map(function (b) { return b + '_mean'; });
    }))
  )
});

// Sanity check in the Code Editor
Map.centerObject(sites, 9);
Map.addLayer(regions, {color: 'red'}, 'sites');
print('sites', sites);
