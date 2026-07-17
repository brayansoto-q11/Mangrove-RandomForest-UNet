// ============================================================
// Landsat Multi-Sensor Image Stack for U-Net Aquaculture Detection
// Tumbes Mangrove System, Northwestern Peru
//
// Purpose:
//   Generates annual multi-band image stacks used as input for
//   the U-Net deep learning model for aquaculture pond detection.
//   This script was used to produce stacks for both model training
//   (reference year) and annual inference (1993/1995–2024 time series).
//
// Supported sensors (auto-selected by year):
//   Landsat 5 TM   → 1993/1995–1999   [HARMONIZED to OLI, Roy et al. 2016]
//   Landsat 7 ETM+ → 2000–2012        [HARMONIZED to OLI, Roy et al. 2016]
//                    (SLC-off gap-filling via median composite)
//   Landsat 8 OLI  → 2013–2021        [native, no harmonization needed]
//   Landsat 9 OLI-2→ 2022–2024        [native, merged with L8 from 2022]
//
// Output: 23-band GeoTIFF stack (6 spectral + 11 static indices +
//         6 temporal metrics) at 30 m resolution, UTM Zone 17S
//
// Authors: [Brayan Soto Quispe, Fernando Alarcon Yllaconse,
//           Ulises Francisco Giraldo Malca, Pablito Marcelo López Serrano]

// GEE Project: ee-brayansotoquispe
// Last updated: 2024
//
// ============================================================


// ============================================================
// SECTION 0: MAIN PARAMETERS
//
// HOW TO USE:
//
// ── Option A: Single-year stack ─────────────────────────────
//    Set MESES_EXTRA = 0. Use when ≥ 5 cloud-free scenes exist.
//
//    Example:
//      var YEAR       = 2020;
//      var MESES_EXTRA = 0;
//      var MES_INI    = 1;
//
// ── Option B: Extended temporal window ──────────────────────
//    Set MESES_EXTRA > 0 to include scenes from months beyond
//    December of the target year. Recommended for years with
//    persistent cloud cover.
//    MESES_EXTRA = 12 extends the window by one full year.
//
//    Example — 18-month window starting January 1996:
//      var YEAR       = 1996;
//      var MESES_EXTRA = 6;
//      var MES_INI    = 1;
//
// Check col.size() in the GEE Console after running:
//   < 3 scenes  → increase MESES_EXTRA or raise CLOUD_THR
//   3–5 scenes  → extended window recommended
//   ≥ 5 scenes  → single-year stack is reliable
// ============================================================
var YEAR        = 1996;
var MESES_EXTRA = 12;   // 0 = target year only | 6, 12 = extended window
var MES_INI     = 1;    // start month of composite window (1 = January)
var CLOUD_THR   = 60;   // cloud cover threshold (%)


// ============================================================
// SECTION 1: STUDY AREA
//
// TO REPLICATE:
//   Upload asset_data/AREA_DEF.geojson to your GEE Assets and
//   replace the path below with your own GEE project path:
//   'projects/YOUR_GEE_PROJECT/assets/AREA_DEF'
// ============================================================
var studyArea = ee.FeatureCollection(
  'projects/YOUR_GEE_PROJECT/assets/AREA_DEF'
);

Map.centerObject(studyArea, 10);


// ============================================================
// SECTION 2: MULTI-SENSOR CONFIGURATION
//
// Each sensor entry defines:
//   collection — GEE image collection ID
//   bands      — native SR band names for that sensor
//   names      — standardized names used throughout the script
//   tcwCoefs   — Tasseled Cap Wetness coefficients (sensor-specific)
//
// TCW references for Surface Reflectance (SR):
//   L5/L7: Crist (1985), RSE 17:301-306 (Reflectance factor native)
//   L8/L9: Zhai et al. (2022), RSE 274:112992 (Derived for C2 L2 SR - 6-band model)
//
// harmonize — whether this sensor's SR reflectance must be converted
//             to OLI-equivalent reflectance before computing indices
//             (true for L5/L7, false for L8/L9 which are already OLI).
// ============================================================
var sensorConfig = {
  L5: {
    collection: 'LANDSAT/LT05/C02/T1_L2',
    bands:     ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
    names:     ['Blue',  'Green', 'Red',   'NIR',   'SWIR1', 'SWIR2'],
    tcwCoefs:  [0.0315,  0.2021,  0.3102,  0.1594,  -0.6806, -0.6109],
    harmonize: true
  },
  L7: {
    collection: 'LANDSAT/LE07/C02/T1_L2',
    bands:     ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
    names:     ['Blue',  'Green', 'Red',   'NIR',   'SWIR1', 'SWIR2'],
    tcwCoefs:  [0.0315,  0.2021,  0.3102,  0.1594,  -0.6806, -0.6109],
    harmonize: true
  },
  L8: {
    collection: 'LANDSAT/LC08/C02/T1_L2',
    bands:     ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
    names:     ['Blue',  'Green', 'Red',   'NIR',   'SWIR1', 'SWIR2'],
    tcwCoefs:  [0.0382,  0.2137,  0.3536,  0.2270,  -0.6108, -0.6351],
    harmonize: false
  },
  L9: {
    collection: 'LANDSAT/LC09/C02/T1_L2',
    bands:     ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
    names:     ['Blue',  'Green', 'Red',   'NIR',   'SWIR1', 'SWIR2'],
    tcwCoefs:  [0.0382,  0.2137,  0.3536,  0.2270,  -0.6108, -0.6351],
    harmonize: false
  }
};

// Automatic sensor selection by year
var sensor;
if      (YEAR >= 2022) { sensor = 'L9'; }
else if (YEAR >= 2013) { sensor = 'L8'; }
else if (YEAR >= 2000) { sensor = 'L7'; }
else                   { sensor = 'L5'; }

var cfg = sensorConfig[sensor];
print('Active sensor:', sensor);
print('Spectral harmonization to OLI (Roy et al. 2016) applied:', cfg.harmonize);


// ============================================================
// SECTION 2B: SPECTRAL HARMONIZATION — Roy et al. (2016)
//
// Reference:
//   Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L.,
//   Kumar, S.S., Egorov, A. (2016). Characterization of Landsat-7
//   to Landsat-8 reflective wavelength and normalized difference
//   vegetation index continuity. Remote Sensing of Environment,
//   185, 57-70.
//
// OLS (Ordinary Least Squares) transformation coefficients,
// Table 2 of Roy et al. (2016), form:
//
//     OLI_band = slope * ETM+_band + intercept
//
// These coefficients transform ETM+ (and, by extension, the
// spectrally near-identical TM/L5) surface reflectance into
// OLI-equivalent reflectance, so that L5/L7-derived reflectance
// and indices are radiometrically consistent with the native
// L8/L9 OLI/OLI-2 composites used from 2013 onward. This is the
// standard practice for building continuous, sensor-harmonized
// Landsat time series (e.g. Roy et al. 2016; USGS/CCDC workflows).
//
// Coefficients below (slope, intercept), in the order
// [Blue, Green, Red, NIR, SWIR1, SWIR2]:
// ============================================================
var roySlope     = [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071];
var royIntercept = [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172];
var harmonizeBandNames = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'];

// Applies OLI harmonization to an image already renamed to the
// standardized band names (Blue, Green, Red, NIR, SWIR1, SWIR2).
// Only used for L5/L7 (cfg.harmonize === true); L8/L9 pass through.
function harmonizeToOLI(img) {
  var harmonizedBands = harmonizeBandNames.map(function(name, i) {
    return img.select(name)
      .multiply(roySlope[i])
      .add(royIntercept[i])
      .rename(name);
  });
  return ee.Image(harmonizedBands.reduce(function(acc, b) {
    return ee.Image(acc).addBands(b);
  })).copyProperties(img, ['system:time_start', 'system:index']);
}

// ============================================================
// SECTION 3: TEMPORAL WINDOW
// ============================================================
var startDate = ee.Date.fromYMD(YEAR, MES_INI, 1);
var endDate   = startDate.advance(12 + MESES_EXTRA, 'month');

print('Composite window:', startDate.format('YYYY-MM-dd'),
      '→', endDate.format('YYYY-MM-dd'));

// ============================================================
// SECTION 4: PREPROCESSING — CLOUD MASKING, SCALING, HARMONIZATION
//
// QA_PIXEL bits used (identical across all C2 L2 sensors):
//   Bit 3 → cloud shadow
//   Bit 4 → snow/ice
//   Bit 5 → cloud
//
// For Landsat 7 SLC-off years (2003–2013), data gaps are
// masked natively in C2 and filled by the median composite
// across multiple acquisition dates.
//
// Scale factors (USGS Collection 2): × 0.0000275, offset −0.2
//
// Processing order per scene:
//   1. Cloud/shadow/snow mask (maskClouds)
//   2. Scale to reflectance + rename to standard bands (applyScale)
//   3. If sensor is L5/L7: harmonize to OLI-equivalent reflectance
//      (harmonizeToOLI), using Roy et al. (2016) coefficients
// ============================================================
function maskClouds(img) {
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)   // no cloud shadow
               .and(qa.bitwiseAnd(1 << 4).eq(0))  // no snow
               .and(qa.bitwiseAnd(1 << 5).eq(0));  // no cloud
  // L7 SLC-off: additional stripe mask via valid pixel check
  if (sensor === 'L7') {
    mask = mask.and(img.select('SR_B1').gt(0));
  }
  return img.updateMask(mask);
}

function applyScale(img) {
  var sr = img.select(cfg.bands)
    .rename(cfg.names)
    .multiply(0.0000275)
    .add(-0.2);
  sr = sr.updateMask(
    sr.reduce(ee.Reducer.min()).gt(-0.1)
      .and(sr.reduce(ee.Reducer.max()).lt(1.1))
  );
  return sr.copyProperties(img, ['system:time_start', 'system:index']);
}

function preprocessImg(img) {
  var scaled = applyScale(maskClouds(img));
  if (cfg.harmonize) {
    scaled = harmonizeToOLI(scaled);
  }
  return scaled;
}

// ============================================================
// SECTION 5: IMAGE COLLECTION
// For 2022+, Landsat 8 and 9 are merged into a single collection.
// Note: L8/L9 never pass through harmonizeToOLI (cfg.harmonize
// is false for both), since they are already native OLI/OLI-2.
// ============================================================
var filters = ee.Filter.and(
  ee.Filter.date(startDate, endDate),
  ee.Filter.bounds(studyArea),
  ee.Filter.lt('CLOUD_COVER', CLOUD_THR)
);

var col = ee.ImageCollection(cfg.collection)
  .filter(filters)
  .map(preprocessImg);

if (YEAR >= 2013) {
  var colL9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filter(filters)
    .map(preprocessImg);
  col = col.merge(colL9);
}

print('Available scenes — ' + YEAR + ' [' + sensor + '] '
      + '(' + MESES_EXTRA + ' extra months):', col.size());


// ============================================================
// SECTION 6: ANNUAL MEDIAN COMPOSITE
// ============================================================
var composite = col.median().clip(studyArea);

var Blue  = composite.select('Blue');
var Green = composite.select('Green');
var Red   = composite.select('Red');
var NIR   = composite.select('NIR');
var SWIR1 = composite.select('SWIR1');
var SWIR2 = composite.select('SWIR2');


// ============================================================
// SECTION 7: STATIC SPECTRAL INDICES (11 indices)
//
// Standardized band names (Blue, Green, Red, NIR, SWIR1, SWIR2)
// allow identical index formulas across all four sensors. Because
// L5/L7 reflectance has already been harmonized to OLI-equivalent
// values (Section 2B), these indices are computed on a radiometrically
// consistent basis across the full 1993/1995–2024 time series.
// ============================================================
var NDVI  = NIR.subtract(Red).divide(NIR.add(Red)).rename('NDVI');
var NDWI  = Green.subtract(NIR).divide(Green.add(NIR)).rename('NDWI');
var MNDWI = Green.subtract(SWIR1).divide(Green.add(SWIR1)).rename('MNDWI');

// AWEI — modified for sensors without coastal aerosol band
var AWEI = Green.multiply(4)
  .subtract(SWIR1.multiply(0.25).add(SWIR2.multiply(2.75)))
  .rename('AWEI');

var EVI = NIR.subtract(Red)
  .divide(NIR.add(Red.multiply(6)).subtract(Blue.multiply(7.5)).add(1))
  .multiply(2.5)
  .rename('EVI');

// Tasseled Cap Wetness — sensor-specific coefficients (Section 2)
var coefs = cfg.tcwCoefs;
var TCW = Blue.multiply(coefs[0])
  .add(Green.multiply(coefs[1]))
  .add(Red.multiply(coefs[2]))
  .add(NIR.multiply(coefs[3]))
  .add(SWIR1.multiply(coefs[4]))
  .add(SWIR2.multiply(coefs[5]))
  .rename('TCW');

var NDPI    = Green.subtract(SWIR1).divide(Green.add(SWIR1)).rename('NDPI');
var NDTI    = Red.subtract(Green).divide(Red.add(Green)).rename('NDTI');
var BSI     = SWIR1.add(Red).subtract(NIR.add(Blue))
                .divide(SWIR1.add(Red).add(NIR.add(Blue))).rename('BSI');
var SAVI    = NIR.subtract(Red).divide(NIR.add(Red).add(0.5)).multiply(1.5).rename('SAVI');
var RATIO_GN = Green.divide(NIR.add(1e-6)).rename('RATIO_GN');


// ============================================================
// SECTION 8: TEMPORAL METRICS (6 features)
//
// Computed from individual scenes before compositing. Since each
// scene has already passed through preprocessImg (mask → scale →
// harmonize), these metrics are also computed on OLI-harmonized
// reflectance for L5/L7 years.
//
// Capture intra-annual variability driven by tidal cycles,
// seasonal flooding, and ENSO precipitation anomalies.
//
//   MNDWI_amp / NDVI_amp  → tidal and seasonal amplitude
//   MNDWI_std / NDVI_std  → temporal variability
//   MNDWI_p50             → median water signal
//   AWEI_amp              → water extent variability
//   FREQ_AGUA             → fraction of time with surface water
// ============================================================
function addIndices(img) {
  var g  = img.select('Green');
  var r  = img.select('Red');
  var n  = img.select('NIR');
  var s1 = img.select('SWIR1');
  var s2 = img.select('SWIR2');
  return img.addBands([
    n.subtract(r).divide(n.add(r)).rename('NDVI'),
    g.subtract(s1).divide(g.add(s1)).rename('MNDWI'),
    g.multiply(4).subtract(s1.multiply(0.25).add(s2.multiply(2.75))).rename('AWEI')
  ]);
}

var colIdx = col.map(addIndices);

// MNDWI temporal metrics
var MNDWI_p10 = colIdx.select('MNDWI').reduce(ee.Reducer.percentile([10]));
var MNDWI_p50 = colIdx.select('MNDWI').reduce(ee.Reducer.percentile([50])).rename('MNDWI_p50');
var MNDWI_p90 = colIdx.select('MNDWI').reduce(ee.Reducer.percentile([90]));
var MNDWI_amp = MNDWI_p90.subtract(MNDWI_p10).rename('MNDWI_amp');
var MNDWI_std = colIdx.select('MNDWI').reduce(ee.Reducer.stdDev()).rename('MNDWI_std');

// NDVI temporal metrics
var NDVI_p10 = colIdx.select('NDVI').reduce(ee.Reducer.percentile([10]));
var NDVI_p90 = colIdx.select('NDVI').reduce(ee.Reducer.percentile([90]));
var NDVI_amp = NDVI_p90.subtract(NDVI_p10).rename('NDVI_amp');
var NDVI_std = colIdx.select('NDVI').reduce(ee.Reducer.stdDev()).rename('NDVI_std');

// AWEI temporal metrics
var AWEI_p10 = colIdx.select('AWEI').reduce(ee.Reducer.percentile([10]));
var AWEI_p90 = colIdx.select('AWEI').reduce(ee.Reducer.percentile([90]));
var AWEI_amp = AWEI_p90.subtract(AWEI_p10).rename('AWEI_amp');

// Annual water frequency: fraction of scenes where MNDWI > 0
var FREQ_AGUA = colIdx.select('MNDWI')
  .map(function(img) { return img.gt(0); })
  .mean()
  .rename('FREQ_AGUA');


// ============================================================
// SECTION 9: FINAL FEATURE STACK
//
// Composition:
//   6  spectral bands  (Blue, Green, Red, NIR, SWIR1, SWIR2)
//   11 static indices  (NDVI, NDWI, MNDWI, AWEI, EVI, TCW,
//                       NDPI, NDTI, BSI, SAVI, RATIO_GN)
//   6  temporal metrics(MNDWI_amp, MNDWI_std, MNDWI_p50,
//                       NDVI_amp, NDVI_std, AWEI_amp, FREQ_AGUA)
//
// This stack is the direct input to the U-Net model.
// Training stacks and inference stacks use identical band order.
// For L5/L7 years, all bands derive from OLI-harmonized reflectance
// (Roy et al. 2016), ensuring the model sees a radiometrically
// consistent input distribution across the whole time series.
// ============================================================
var unetStack = composite
  .select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
  .addBands(NDVI).addBands(NDWI).addBands(MNDWI)
  .addBands(AWEI).addBands(EVI).addBands(TCW).addBands(NDPI)
  .addBands(NDTI).addBands(BSI).addBands(SAVI).addBands(RATIO_GN)
  .addBands(MNDWI_amp).addBands(MNDWI_std).addBands(MNDWI_p50)
  .addBands(NDVI_amp).addBands(NDVI_std)
  .addBands(AWEI_amp).addBands(FREQ_AGUA);

print('Stack band names:', unetStack.bandNames());
print('Total bands (expected 23):', unetStack.bandNames().length());


// ============================================================
// SECTION 10: VISUALIZATION
// ============================================================
Map.addLayer(
  composite,
  {bands: ['Red', 'Green', 'Blue'], min: 0.0, max: 0.15, gamma: 1.4},
  sensor + ' ' + YEAR + ' — True color (window: +' + MESES_EXTRA + ' months)'
    + (cfg.harmonize ? ' [OLI-harmonized]' : '')
);

// ============================================================
// SECTION 11: EXPORT — U-NET INPUT STACK (GeoTIFF)
//
// Output naming convention:
//   stack_unet_tumbes_{YEAR}_{SENSOR}_ext{WINDOW}m[_OLIharm]
//   e.g. stack_unet_tumbes_1996_L5_ext12m_OLIharm
//
// Projection: UTM Zone 17S (EPSG:32717)
// Resolution: 30 m
// Folder:     Unet_Acuicultura_Datos (Google Drive)
// ============================================================
var exportName = 'stack_unet_tumbes_' + YEAR + '_' + sensor
                 + '_ext' + MESES_EXTRA + 'm'
                 + (cfg.harmonize ? '_OLIharm' : '');

Export.image.toDrive({
  image:       unetStack.toFloat(),
  description: exportName,
  folder:      'Unet_Acuicultura_Datos',
  region:      studyArea.geometry(),
  scale:       30,
  crs:         'EPSG:32717',
  fileFormat:  'GeoTIFF',
  maxPixels:   1e13
});
