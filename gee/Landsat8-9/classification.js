// ============================================================
// Supervised Land Cover Classification · Epoch 3 · 2013–2024
// Tumbes Mangrove System, Northwestern Peru
// Sensors: Landsat 8 OLI + Landsat 9 OLI-2 (USGS Collection 2 SR)
//
// Five land cover classes:
//   1 - Mangrove
//   2 - Water Bodies
//   3 - Bare Soil
//   4 - Mudflat
//   5 - Dry Forest
//
// Classifier: Random Forest (500 trees)
// Feature space: 22 spectral and temporal features
//
// Authors: [Brayan Soto Quispe, Fernando Alarcon Yllaconse,
//           Ulises Francisco Giraldo Malca]
// Institution: Universidad Peruana de Ciencias Aplicadas (UPC)
// GEE Project: ee-brayansotoquispe
// Last updated: 2024

// ============================================================


// ============================================================
// SECTION 0: MAIN PARAMETER
// Change only this value to reproduce classification for any year
// ============================================================
var YEAR = 2013;

// ============================================================
// SECTION 1: STUDY AREA BOUNDARIES
// ============================================================
// Import your own assets (from "asset_data") before running the script.

var areaManglar = ee.FeatureCollection('projects/your_username/assets/AREA_DEF');
var SNLMT = ee.FeatureCollection('projects/your_username/assets/SNLMT');
var Zona_Amortiguamiento = ee.FeatureCollection('projects/your_username/assets/Zona_Amortiguamiento');
Map.centerObject(studyArea, 10);


// ============================================================
// SECTION 2: SPECTRAL BANDS (Landsat 8/9 OLI)
// B2=Blue  B3=Green  B4=Red  B5=NIR  B6=SWIR1  B7=SWIR2
// ============================================================
var BANDS = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];

// ============================================================
// SECTION 3: PREPROCESSING — CLOUD MASKING AND SCALING
// Applies USGS Collection 2 scale factors (×0.0000275, offset −0.2)
// and masks cloud and cloud-shadow pixels using QA_PIXEL band
// ============================================================
function maskClouds(img) {
  var qa = img.select('QA_PIXEL');
  return img.updateMask(
    qa.bitwiseAnd(1 << 3).eq(0)
      .and(qa.bitwiseAnd(1 << 5).eq(0))
  );
}

function applyScaleFactor(img) {
  var sr = img.select(BANDS).multiply(0.0000275).add(-0.2);
  return sr
    .updateMask(
      sr.reduce(ee.Reducer.min()).gt(-0.1)
        .and(sr.reduce(ee.Reducer.max()).lt(1.1))
    )
    .copyProperties(img, ['system:time_start', 'system:index']);
}

function preprocessL89(img) {
  return applyScaleFactor(maskClouds(img));
}

// ============================================================
// SECTION 4: ANNUAL LANDSAT 8 + 9 IMAGE COLLECTION
// Both sensors are preprocessed identically and merged
// Cloud cover threshold: < 80%
// ============================================================
var startDate = ee.Date.fromYMD(YEAR, 1, 1);
var endDate   = ee.Date.fromYMD(YEAR, 12, 31);

var filters = ee.Filter.and(
  ee.Filter.date(startDate, endDate),
  ee.Filter.bounds(studyArea),
  ee.Filter.lt('CLOUD_COVER', 80)
);

var l8  = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filter(filters).map(preprocessL89);
var l9  = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filter(filters).map(preprocessL89);
var col = l8.merge(l9);

print('Available scenes for ' + YEAR + ':', col.size());


// ============================================================
// SECTION 5: ANNUAL MEDIAN COMPOSITE
// ============================================================
var composite = col.median().clip(studyArea);

// ============================================================
// SECTION 6: SPECTRAL INDICES (14 indices)
//
// Index selection rationale per target class:
// ─────────────────────────────────────────────────────────────
// MANGROVE:      NDVI, EVI (dense evergreen canopy),
//                CMRI (mangrove-specific; Jia et al. 2021),
//                NDMI / LSWI (high foliar water content)
//
// WATER BODIES:  AWEI_sh (turbid and shaded water; Feyisa et al. 2014),
//                MNDWI (water vs. soil/vegetation; Xu 2006),
//                NDWI (open water reference; McFeeters 1996)
//
// BARE SOIL:     BSI (exposed dry soil; Rikimaru et al. 2002),
//                RATIO_RS (RED/SWIR1 elevated in saline mudflats)
//
// MUDFLAT:       SAVI (sparse vegetation on wet soil; L=0.5; Huete 1988),
//                NDVI_p10 (captures unvegetated flood periods),
//                MNDWI_amp (tidal variability separates mudflat from dry soil)
//
// DRY FOREST:    EVI (low-LAI canopy more sensitive than NDVI),
//                SAVI (sparse vegetation on dry soil),
//                NDMI (low foliar moisture distinguishes from mangrove)
// ============================================================

// ── Vegetation indices ──────────────────────────────────────
var NDVI = composite.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

var EVI = composite.expression(
  '2.5 * ((NIR - RED) / (NIR + 6.0 * RED - 7.5 * BLUE + 1.0))',
  { NIR:  composite.select('SR_B5'),
    RED:  composite.select('SR_B4'),
    BLUE: composite.select('SR_B2') }
).rename('EVI');

var SAVI = composite.expression(
  '((NIR - RED) / (NIR + RED + L)) * (1 + L)',
  { NIR: composite.select('SR_B5'),
    RED: composite.select('SR_B4'),
    L: 0.5 }
).rename('SAVI');

var NDMI = composite.normalizedDifference(['SR_B5', 'SR_B6']).rename('NDMI');

var LSWI = composite.normalizedDifference(['SR_B5', 'SR_B7']).rename('LSWI');

// CMRI · Jia et al. 2021 — positive over mangrove, negative over mudflat
var CMRI = NDVI.subtract(
  composite.normalizedDifference(['SR_B3', 'SR_B5'])
).rename('CMRI');

// ── Water indices ───────────────────────────────────────────
var NDWI   = composite.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');
var MNDWI  = composite.normalizedDifference(['SR_B3', 'SR_B6']).rename('MNDWI');

// AWEI_sh = BLUE + 2.5×GREEN − 1.5×(NIR + SWIR1) − 0.25×SWIR2
var AWEI_sh = composite.expression(
  'BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2',
  { BLUE:  composite.select('SR_B2'),
    GREEN: composite.select('SR_B3'),
    NIR:   composite.select('SR_B5'),
    SWIR1: composite.select('SR_B6'),
    SWIR2: composite.select('SR_B7') }
).rename('AWEI_sh');

// ── Bare soil indices ───────────────────────────────────────
// BSI = ((SWIR1 + RED) − (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))
var BSI = composite.expression(
  '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))',
  { SWIR1: composite.select('SR_B6'),
    RED:   composite.select('SR_B4'),
    NIR:   composite.select('SR_B5'),
    BLUE:  composite.select('SR_B2') }
).rename('BSI');

var RATIO_RS = composite.select('SR_B4')
  .divide(composite.select('SR_B6').add(0.0001))
  .rename('RATIO_RS');


// ============================================================
// SECTION 7: INTRA-ANNUAL TEMPORAL METRICS (8 features)
//
//   Mangrove:      stable high NDVI year-round (evergreen)
//   Water Bodies:  stable high MNDWI year-round (permanent)
//   Mudflat:       high NDVI variability (seasonal inundation),
//                  high MNDWI amplitude (tidal cycle)
//   Dry Forest:    NDVI low in dry season, higher in wet season
//   Bare Soil:     NDVI ≈ 0 year-round, stable negative MNDWI
// ============================================================
var ndviCol = col.map(function(img) {
  return img.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')
    .set('system:time_start', img.get('system:time_start'));
});

var mndwiCol = col.map(function(img) {
  return img.normalizedDifference(['SR_B3', 'SR_B6']).rename('MNDWI')
    .set('system:time_start', img.get('system:time_start'));
});

var NDVI_p10  = ndviCol.reduce(ee.Reducer.percentile([10])).clip(studyArea).rename('NDVI_p10');
var NDVI_p90  = ndviCol.reduce(ee.Reducer.percentile([90])).clip(studyArea).rename('NDVI_p90');
var NDVI_std  = ndviCol.reduce(ee.Reducer.stdDev()).clip(studyArea).rename('NDVI_std');
var NDVI_amp  = NDVI_p90.subtract(NDVI_p10).rename('NDVI_amp');

var MNDWI_p10 = mndwiCol.reduce(ee.Reducer.percentile([10])).clip(studyArea).rename('MNDWI_p10');
var MNDWI_p90 = mndwiCol.reduce(ee.Reducer.percentile([90])).clip(studyArea).rename('MNDWI_p90');
var MNDWI_amp = MNDWI_p90.subtract(MNDWI_p10).rename('MNDWI_amp');
var MNDWI_std = mndwiCol.reduce(ee.Reducer.stdDev()).clip(studyArea).rename('MNDWI_std');

// Annual water frequency: permanent water ≈ 1.0 · mudflat ≈ 0.3–0.7 · bare soil ≈ 0.0
var FREQ_AGUA = mndwiCol
  .map(function(img) { return img.gt(0).rename('MNDWI'); })
  .mean().clip(studyArea).rename('FREQ_AGUA');

// ============================================================
// SECTION 8: FULL FEATURE STACK (22 features)
// 6 spectral bands + 8 static indices + 8 temporal metrics
// ============================================================
var features = composite.select(BANDS)
  .addBands(NDVI).addBands(EVI).addBands(SAVI)
  .addBands(NDMI).addBands(LSWI).addBands(CMRI)
  .addBands(NDWI).addBands(MNDWI).addBands(AWEI_sh)
  .addBands(BSI).addBands(RATIO_RS)
  .addBands(NDVI_p10).addBands(NDVI_p90)
  .addBands(NDVI_std).addBands(NDVI_amp)
  .addBands(MNDWI_p10).addBands(MNDWI_p90)
  .addBands(MNDWI_amp).addBands(MNDWI_std)
  .addBands(FREQ_AGUA);

var ALL_BANDS = features.bandNames();
print('Total features (expected: 22):', ALL_BANDS.size());
print('Feature stack bands:', ALL_BANDS);


// ============================================================
// SECTION 9: TRAINING POLYGONS — 5 CLASSES
//
// Class encoding (TC property):
//   TC = 1  Mangrove       dense evergreen canopy, interior pixels only
//   TC = 2  Water_Bodies   tidal channels, ponds, open water, aquaculture
//   TC = 3  Bare_Soil      tidal flats, dry sandy areas, unvegetated shores
//   TC = 4  Mudflat        flat expanses of mud or silt
//   TC = 5  Dry_Forest     inland dry-deciduous forest, no tidal influence
//
// Digitization guidelines:
//   · Polygons placed ≥ 1–2 pixels from class boundaries
//   · Minimum 50 polygons per class; target 80–120
//   · Reference image: false-color composite (SR_B5/B4/B3) at zoom 15–16
//
// Polygon counts:
//   Mangrove: 43 polys | Water_Bodies: 98 pts, 58 polys
//   Bare_Soil: 182 pts, 31 polys | Mudflat: 17 polys | Dry_Forest: 20 polys
// ============================================================
// ============================================================
//
// Two modes are available. Only ONE should be active at a time.
//
// MODE A (default): uses geometry imports drawn directly in the
//   GEE Code Editor. Works only with the authors' original script
//   where the imports Mangrove, Water_Bodies, Bare_soil, Mudflat,
//   and Dry_forest are defined as drawn geometries.
//
// MODE B (for replication): loads polygons from an uploaded GEE
//   Asset. Use this if you downloaded training_polygons_points_E3_L89.geojson
//   from the repository and uploaded it to your own GEE project.
//   To activate: delete the MODE A block and uncomment MODE B.
//
// Training polygons represent land cover conditions for one
// reference year per sensor epoch. For years with substantially
// different conditions, polygons should be reviewed and updated.
// ============================================================

// ── MODE A: AUTHORS' DIGITIZED IMPORTS (active by default) ──
// Delete this block and activate MODE B if replicating.
var mangrove  = Mangrove.map(function(f)     { return f.set('TC', 1); });
var water     = Water_Bodies.map(function(f) { return f.set('TC', 2); });
var bareSoil  = Bare_soil.map(function(f)    { return f.set('TC', 3); });
var mudflat   = Mudflat.map(function(f)      { return f.set('TC', 4); });
var dryForest = Dry_forest.map(function(f)  { return f.set('TC', 5); });

var trainingPolygons = mangrove
  .merge(water)
  .merge(bareSoil)
  .merge(mudflat)
  .merge(dryForest);

// ── MODE B: LOAD FROM ASSET (for replication) ───────────────
// To activate: delete MODE A above, then uncomment these lines.
// Replace YOUR_GEE_PROJECT with your GEE project ID.
// Change asset suffix: E1_L5 | E1_L7 | E3_L89
//
// var trainingAsset = ee.FeatureCollection(
//   'projects/YOUR_GEE_PROJECT/assets/training_polygons_points_E3_L89'
// );
// var mangrove  = trainingAsset.filter(ee.Filter.eq('TC', 1));
// var water     = trainingAsset.filter(ee.Filter.eq('TC', 2));
// var bareSoil  = trainingAsset.filter(ee.Filter.eq('TC', 3));
// var mudflat   = trainingAsset.filter(ee.Filter.eq('TC', 4));
// var dryForest = trainingAsset.filter(ee.Filter.eq('TC', 5));
// var trainingPolygons = trainingAsset;

// ============================================================
// SECTION 10: SPECTRAL SAMPLING
// ============================================================
var samples = features.sampleRegions({
  collection:  trainingPolygons,
  properties:  ['TC'],
  scale:        30,
  geometries:   true,
  tileScale:    4
}).filter(ee.Filter.notNull(ALL_BANDS.getInfo()));

print('Samples per class (TC 1–5):', samples.aggregate_histogram('TC'));
print('Total valid samples:', samples.size());


// ============================================================
// SECTION 11: POLYGON-BASED TRAIN / TEST SPLIT (70 / 30)
//
// All pixels from the same polygon are assigned exclusively to
// train OR test, preventing spatial autocorrelation leakage.
// Reference: Karasiak et al. (2022), Remote Sensing of Environment
// ============================================================
var polygonsWithRandom = trainingPolygons.randomColumn('random_poly', 789);

var samplesWithPolyRand = features.sampleRegions({
  collection:  polygonsWithRandom,
  properties:  ['TC', 'random_poly'],
  scale:        30,
  geometries:   true,
  tileScale:    4
}).filter(ee.Filter.notNull(ALL_BANDS.getInfo()));

var trainSet = samplesWithPolyRand.filter(ee.Filter.lt('random_poly',  0.7));
var testSet  = samplesWithPolyRand.filter(ee.Filter.gte('random_poly', 0.7));

print('▶ Training samples:', trainSet.size());
print('▶ Test samples:',     testSet.size());
print('▶ Class distribution — train:', trainSet.aggregate_histogram('TC'));
print('▶ Class distribution — test: ', testSet.aggregate_histogram('TC'));


// ============================================================
// SECTION 12: RANDOM FOREST CLASSIFIER
//
// Hyperparameter justification:
//   numberOfTrees = 500   stable variable importance estimates
//                         (Belgiu & Drăguţ 2016, ISPRS J. Photogramm.)
//   minLeafPopulation = 5 prevents overfitting on small classes
//                         (Gislason et al. 2006, Pattern Recognit. Lett.)
//   bagFraction = 0.7     standard for land cover mapping
//   variablesPerSplit = 5 ≈ sqrt(22 features); Breiman (2001)
//   seed = 456            fixed for full reproducibility
// ============================================================
var rfClassifier = ee.Classifier.smileRandomForest({
  numberOfTrees:     500,
  minLeafPopulation:   5,
  bagFraction:         0.7,
  variablesPerSplit:   5,
  seed:              456
}).train({
  features:        trainSet,
  classProperty:   'TC',
  inputProperties: ALL_BANDS
});

print('Variable importance:', ee.Dictionary(rfClassifier.explain().get('importance')));


// ============================================================
// SECTION 13: CLASSIFICATION
// ============================================================
var rawClassification = features.classify(rfClassifier);


// ============================================================
// SECTION 14: POST-CLASSIFICATION SPATIAL FILTERING
//
// Step 1: Modal filter (2 iterations) — removes salt-and-pepper noise
// Step 2: Morphological closing for Water Bodies (TC=2) —
//         regularizes channel and pond boundaries
// Step 3: Second modal pass — consolidates class edges
// Step 4: MMU filter — removes fragments < 1 pixel (< 900 m²)
// ============================================================
var smoothed1 = rawClassification.focal_mode({
  radius: 1, units: 'pixels', iterations: 2
});

var waterMask    = rawClassification.eq(2);
var waterClosed  = waterMask
  .focal_max({radius: 2, units: 'pixels', iterations: 1})
  .focal_min({radius: 2, units: 'pixels', iterations: 1});
var waterRefined = waterClosed
  .and(rawClassification.eq(2).focal_max({radius: 3, units: 'pixels', iterations: 1}));
var smoothed2 = smoothed1.where(waterRefined, 2);

var smoothed3 = smoothed2.focal_mode({
  radius: 1, units: 'pixels', iterations: 1
});

var finalClassification = smoothed3
  .where(smoothed3.connectedPixelCount(25, true).lt(1), 0)
  .selfMask()
  .rename('classification');


// ============================================================
// SECTION 15: VISUALIZATION
// ============================================================
var CLASS_PALETTE = [
  '#1A6B1A',  // TC=1 Mangrove      → dark forest green
  '#08519C',  // TC=2 Water Bodies  → deep blue
  '#D4A96A',  // TC=3 Bare Soil     → sandy ochre
  '#74C476',  // TC=4 Mudflat       → light yellow-green
  '#C9B96A'   // TC=5 Dry Forest    → golden brown
];
var visParams = {min: 1, max: 5, palette: CLASS_PALETTE};

Map.addLayer(composite,
  {bands: ['SR_B5', 'SR_B4', 'SR_B3'], min: 0, max: 0.3},
  'False color ' + YEAR);
Map.addLayer(composite,
  {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.3},
  'True color ' + YEAR, false);

Map.addLayer(MNDWI,     {min: -0.5, max: 0.5, palette: ['#8B4513','white','#0055E3']}, 'MNDWI', false);
Map.addLayer(AWEI_sh,   {min: -0.3, max: 0.3, palette: ['#8B4513','white','#0055E3']}, 'AWEI_sh (turbid water)', false);
Map.addLayer(BSI,       {min: -0.5, max: 0.5, palette: ['white','#D4A96A']},           'BSI (bare soil)', false);
Map.addLayer(CMRI,      {min: -0.3, max: 0.5, palette: ['white','#006400']},           'CMRI (mangrove-specific)', false);
Map.addLayer(NDVI_amp,  {min:    0, max: 0.5, palette: ['white','#FF6600']},           'NDVI amplitude (seasonality)', false);
Map.addLayer(MNDWI_amp, {min:    0, max: 0.4, palette: ['white','#0055E3']},           'MNDWI amplitude (tidal variability)', false);
Map.addLayer(FREQ_AGUA, {min:    0, max:   1, palette: ['white','#0055E3']},           'Annual water frequency', false);

Map.addLayer(rawClassification,   visParams, 'RF raw '         + YEAR, false);
Map.addLayer(smoothed1,           visParams, 'RF smoothed',            false);
Map.addLayer(finalClassification, visParams, 'RF final + MMU ' + YEAR);

Map.addLayer(SNLMT,      {color: 'yellow'}, 'SNLMT (national sanctuary)', false);
Map.addLayer(bufferZone, {color: 'orange'}, 'Buffer zone',                false);


// ============================================================
// SECTION 16: ACCURACY ASSESSMENT
//
// Metrics: OA (Overall Accuracy), Kappa, PA (Producer's Accuracy),
//          UA (User's Accuracy), F1 (harmonic mean of PA and UA)
// Reference: Congalton & Green (2019); Foody (2002, RSE)
//
// Note: errorMatrix rows/columns include a zero row for background;
// TC class indices (1–5) therefore map directly to list positions 1–5
// ============================================================
var testClassified = testSet.classify(rfClassifier);
var confMatrix     = testClassified.errorMatrix('TC', 'classification');

var OA    = ee.Number(confMatrix.accuracy());
var Kappa = ee.Number(confMatrix.kappa());
var PA_list = ee.List(confMatrix.producersAccuracy().toList()).flatten();
var UA_list = ee.List(confMatrix.consumersAccuracy().toList()).flatten();

print('════════════════════════════════════════════════════════');
print('ACCURACY ASSESSMENT · YEAR ' + YEAR + ' · 5 CLASSES');
print('Validation: polygon-based 70/30 split (no spatial leakage)');
print('════════════════════════════════════════════════════════');
print('Confusion matrix:', confMatrix);
print('Overall Accuracy (OA):', OA);
print('Kappa coefficient:', Kappa);
print('Producer\'s Accuracy (PA) by class:', PA_list);
print('User\'s Accuracy (UA) by class:', UA_list);


// ============================================================
// SECTION 17: CLASS LABEL DICTIONARY
// ============================================================
var CLASS_NAMES = ee.Dictionary({
  '1': 'Mangrove',
  '2': 'Water_Bodies',
  '3': 'Bare_Soil',
  '4': 'Mudflat',
  '5': 'Dry_Forest'
});


// ============================================================
// SECTION 18: EXPORT — ACCURACY TABLES (Google Drive)
// ============================================================
var TC_LIST = [1, 2, 3, 4, 5];

// File 1: Global accuracy metrics
Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, {
      Year:         YEAR,
      Epoch:        3,
      Sensor:       'L8+L9',
      N_Classes:    5,
      Split_Method: 'Polygon_70-30',
      N_Train:      trainSet.size(),
      N_Test:       testSet.size(),
      OA:           OA,
      Kappa:        Kappa
    })
  ]),
  description: 'Acc_Global_5C_E3_' + YEAR,
  folder:      'GEE_Manglares_Tumbes',
  fileFormat:  'CSV'
});

// File 2: Per-class PA, UA, F1
var perClassAccuracy = ee.FeatureCollection(
  ee.List(TC_LIST).map(function(tcNum) {
    var tc    = ee.Number(tcNum);
    var PA    = ee.Number(PA_list.get(tc));
    var UA    = ee.Number(UA_list.get(tc));
    var denom = PA.add(UA);
    var F1 = ee.Algorithms.If(
      denom.gt(0),
      ee.Number(2).multiply(PA).multiply(UA).divide(denom),
      ee.Number(0)
    );
    return ee.Feature(null, {
      Year:   YEAR,
      Epoch:  3,
      Sensor: 'L8+L9',
      TC:     tc,
      Class:  CLASS_NAMES.get(tc.format()),
      PA:     PA,
      UA:     UA,
      F1:     ee.Number(F1)
    });
  })
);

Export.table.toDrive({
  collection:  perClassAccuracy,
  description: 'Acc_PerClass_5C_E3_' + YEAR,
  folder:      'GEE_Manglares_Tumbes',
  fileFormat:  'CSV'
});

// File 3: Full confusion matrix
var matrixArray = ee.List(confMatrix.array().toList());

var confMatrixTable = ee.FeatureCollection(
  ee.List(TC_LIST).map(function(claseTC) {
    var tcNum  = ee.Number(claseTC);
    var rowRaw = ee.List(matrixArray.get(tcNum));
    return ee.Feature(null, {
      Year:           YEAR,
      Epoch:          3,
      Sensor:         'L8+L9',
      TC:             tcNum,
      Actual_Class:   CLASS_NAMES.get(tcNum.format()),
      Pred_Mangrove:  ee.Number(rowRaw.get(1)),
      Pred_Water:     ee.Number(rowRaw.get(2)),
      Pred_BareSoil:  ee.Number(rowRaw.get(3)),
      Pred_Mudflat:   ee.Number(rowRaw.get(4)),
      Pred_DryForest: ee.Number(rowRaw.get(5))
    });
  })
);

Export.table.toDrive({
  collection:  confMatrixTable,
  description: 'Acc_ConfMatrix_5C_E3_' + YEAR,
  folder:      'GEE_Manglares_Tumbes',
  fileFormat:  'CSV'
});


// ============================================================
// SECTION 19: AREA STATISTICS BY MANAGEMENT ZONE
// Output columns: Year · Zone · Class · TC · Area_ha
// ============================================================
var pixelAreaHa = ee.Image.pixelArea().divide(10000);

function areaByZone(image, zone, zoneName) {
  var stats = pixelAreaHa.addBands(image.rename('class'))
    .reduceRegion({
      reducer:   ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
      geometry:  zone.geometry(),
      scale:     30,
      maxPixels: 1e13,
      tileScale: 4
    });
  return ee.List(stats.get('groups')).map(function(item) {
    item = ee.Dictionary(item);
    var num = ee.Number(item.get('class')).toInt().format();
    return ee.Feature(null, {
      Year:    YEAR,
      Epoch:   3,
      Sensor:  'L8+L9',
      Zone:    zoneName,
      Class:   CLASS_NAMES.get(num),
      TC:      ee.Number(item.get('class')).toInt(),
      Area_ha: ee.Number(item.get('sum')).multiply(100).round().divide(100)
    });
  });
}

var areaStudyArea  = areaByZone(finalClassification, studyArea,  'Study_Area');
var areaSNLMT      = areaByZone(finalClassification, SNLMT,      'SNLMT');
var areaBufferZone = areaByZone(finalClassification, bufferZone, 'Buffer_Zone');

var areaTable = ee.FeatureCollection(areaStudyArea)
  .merge(ee.FeatureCollection(areaSNLMT))
  .merge(ee.FeatureCollection(areaBufferZone));

print('Area by class and management zone (ha):', areaTable);

Export.table.toDrive({
  collection:  areaTable,
  description: 'Areas_5C_ByZone_E3_' + YEAR,
  folder:      'GEE_Manglares_Tumbes',
  fileFormat:  'CSV'
});


// ============================================================
// SECTION 20: EXPORT — FINAL CLASSIFIED MAP (GeoTIFF)
// Projection: UTM Zone 17S (EPSG:32717) · Spatial resolution: 30 m
// ============================================================
Export.image.toDrive({
  image:          finalClassification.toByte(),
  description:    'Mangrove_5C_E3_L89' + YEAR,
  folder:         'GEE_Manglares_Tumbes',
  fileNamePrefix: 'mangrove_5c_' + YEAR,
  region:         studyArea.geometry(),
  scale:          30,
  crs:            'EPSG:32717',
  fileFormat:     'GeoTIFF',
  maxPixels:      1e13
});