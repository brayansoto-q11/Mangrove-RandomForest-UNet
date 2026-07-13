// ============================================================
// Supervised Land Cover Classification · Epoch 1 · 2000–2012
// Tumbes Mangrove System, Northwestern Peru
// Sensor: Landsat 7 ETM+ (USGS Collection 2 SR)
//
// Five land cover classes:
//   1 - Mangrove
//   2 - Water Bodies
//   3 - Bare Soil
//   4 - Mudflat (tidal mudflat / wet bare soil)
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
//
// ============================================================

// ============================================================
// LANDSAT 7 ETM+ — SCAN LINE CORRECTOR (SLC) FAILURE
// ============================================================
//
// On 31 May 2003, the Scan Line Corrector of Landsat 7 ETM+
// failed permanently. Post-failure images show parallel wedge-
// shaped data gaps covering approximately 22% of each scene.
// Gaps are widest at scene edges and nearly absent at nadir.
//
// Periods:
//   Epoch 1A (SLC-on):  1999–2002  → no gap problem
//   Epoch 1B (SLC-off): 2003–2013  → gap-filling required
//
// Gap-filling strategy implemented here:
//   Primary:  Annual median composite — valid pixels from
//             different acquisition dates fill gaps from other
//             dates. Effective when ≥ 5 cloud-free scenes are
//             available per year.
//   Fallback: 25th-percentile composite used to fill any
//             residual NoData pixels after the median composite.
//   Quality:  Pixels with fewer than 3 valid observations are
//             flagged and their temporal metrics are replaced
//             with neutral fill values to avoid biasing
//             the classifier.
//
// Reference: USGS (2021). Landsat 7 ETM+ SLC-off Products.
//            Scaramuzza et al. (2004). SLC-off Gap-Fill Algorithm.
// ============================================================


// ============================================================
// SECTION 0: MAIN PARAMETER
// Change only this value to reproduce classification for any year
// Valid range: 1999–2002 (SLC-on) · 2003–2013 (SLC-off)
// ============================================================
var YEAR = 2000;

// ============================================================
// SECTION 1: STUDY AREA BOUNDARIES
// ============================================================
var areaManglar = ee.FeatureCollection('projects/your_username/assets/AREA_DEF');
var SNLMT = ee.FeatureCollection('projects/your_username/assets/SNLMT');
var Zona_Amortiguamiento = ee.FeatureCollection('projects/your_username/assets/Zona_Amortiguamiento');

Map.centerObject(studyArea, 10);


// ============================================================
// SECTION 2: SPECTRAL BANDS — LANDSAT 7 ETM+
//
// Band equivalencies with Landsat 8/9 OLI:
//   L7 B1 (0.45–0.52 µm) → Blue   (≈ L8 B2)
//   L7 B2 (0.52–0.60 µm) → Green  (≈ L8 B3)
//   L7 B3 (0.63–0.69 µm) → Red    (≈ L8 B4)
//   L7 B4 (0.77–0.90 µm) → NIR    (≈ L8 B5)
//   L7 B5 (1.55–1.75 µm) → SWIR1  (≈ L8 B6)
//   L7 B7 (2.08–2.35 µm) → SWIR2  (≈ L8 B7)
//
// Bands are renamed to the L8/9 naming convention (SR_B2–SR_B7)
// after scaling to allow identical index formulas across epochs.
// ============================================================
var BANDS_L7    = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];
var BANDS_ALIAS = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];


// ============================================================
// SECTION 3: PREPROCESSING — CLOUD MASKING AND SCALING
// ============================================================

// 3A: Cloud and cloud-shadow masking using QA_PIXEL
//     Radiometric saturation flag also applied (QA_RADSAT)
function maskCloudsL7(img) {
  var qa     = img.select('QA_PIXEL');
  var cloud  = qa.bitwiseAnd(1 << 5).eq(0);
  var shadow = qa.bitwiseAnd(1 << 3).eq(0);
  var noSat  = img.select('QA_RADSAT').eq(0);
  return img.updateMask(cloud.and(shadow).and(noSat));
}

// 3B: Scale factor application and band renaming
//     USGS Collection 2 SR scale: × 0.0000275, offset −0.2
//     SLC-off gaps are already masked natively in C02 T1_L2;
//     the median composite fills them across acquisition dates.
function applyScaleL7(img) {
  var sr = img.select(BANDS_L7)
    .multiply(0.0000275)
    .add(-0.2)
    .rename(BANDS_ALIAS);
  return sr
    .updateMask(
      sr.reduce(ee.Reducer.min()).gt(-0.1)
        .and(sr.reduce(ee.Reducer.max()).lt(1.1))
    )
    .copyProperties(img, ['system:time_start', 'system:index']);
}

function preprocessL7(img) {
  return applyScaleL7(maskCloudsL7(img));
}

// ============================================================
// SECTION 4: ANNUAL LANDSAT 7 IMAGE COLLECTION
//
// The annual median composite is the primary gap-filling
// strategy for SLC-off years. Each pixel takes the median
// of all valid observations across the year. If a pixel
// falls in an SLC gap on date X but has valid data on
// dates X ± 16 days (orbital repeat cycle), the median
// recovers it automatically.
//
// Rule of thumb:
//   ≥ 5 valid scenes  → ~95% of SLC gaps filled
//   ≥ 8 valid scenes  → near-complete spatial coverage
// ============================================================
var startDate = ee.Date.fromYMD(YEAR, 1, 1);
var endDate   = ee.Date.fromYMD(YEAR, 12, 31);

var filters = ee.Filter.and(
  ee.Filter.date(startDate, endDate),
  ee.Filter.bounds(studyArea),
  ee.Filter.lt('CLOUD_COVER', 80)
);

var col = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filter(filters)
  .map(preprocessL7);

print('Available L7 scenes for ' + YEAR + ':', col.size());

// Observation density map — used to flag low-quality pixels
// Minimum recommended: ≥ 3 valid observations per pixel
var obsDensity = col.select('SR_B5').count().clip(studyArea);

Map.addLayer(obsDensity,
  {min: 0, max: 20, palette: ['red', 'yellow', 'green']},
  'Valid observations per pixel (SLC diagnostic)', false);

print('Minimum observations per pixel (target > 5):',
  obsDensity.reduceRegion({
    reducer:  ee.Reducer.min(),
    geometry: studyArea.geometry(),
    scale: 30, maxPixels: 1e9
  })
);

// Primary composite: annual median
var composite = col.median().clip(studyArea);

// Coverage check — fraction of pixels with valid data
print('Composite spatial coverage (1.0 = 100%, target > 0.95):',
  composite.select('SR_B4').mask().reduceRegion({
    reducer:  ee.Reducer.mean(),
    geometry: studyArea.geometry(),
    scale: 30, maxPixels: 1e9
  })
);

// Fallback composite: 25th percentile fills residual gaps
var compositeFallback = col.select(BANDS_ALIAS)
  .reduce(ee.Reducer.percentile([25]))
  .rename(BANDS_ALIAS);

var compositeFilled = composite.unmask(compositeFallback).clip(studyArea);

// ============================================================
// SECTION 5: SPECTRAL HARMONIZATION ETM+ → OLI
//
// Systematic inter-sensor differences between L7 ETM+ and
// L8/9 OLI are corrected using the OLS regression coefficients
// of Roy et al. (2016): OLI = a + b × ETM+ (surface reflectance)
//
// Coefficients (Roy et al. 2016, RSE 185:57-70, Table 2):
//   Blue  (B2): slope = 0.8474, intercept = 0.0003
//   Green (B3): slope = 0.8483, intercept = 0.0088
//   Red   (B4): slope = 0.9047, intercept = 0.0061
//   NIR   (B5): slope = 0.8462, intercept = 0.0412
//   SWIR1 (B6): slope = 0.8937, intercept = 0.0254
//   SWIR2 (B7): slope = 0.9071, intercept = 0.0172
//
// Harmonization is applied here to ensure spectral consistency
// when comparing Epoch 1 (L7) maps against Epoch 3 (L8/9) maps.
// For single-sensor, single-epoch classifications it is optional.
// ============================================================
var HARMONIZE_TO_OLI = true;

function harmonizeL7toL8(img) {
  img = ee.Image(img);
  return ee.Image.cat([
    img.select('SR_B2').multiply(0.8474).add(0.0003).rename('SR_B2'),  // Blue
    img.select('SR_B3').multiply(0.8483).add(0.0088).rename('SR_B3'),  // Green
    img.select('SR_B4').multiply(0.9047).add(0.0061).rename('SR_B4'),  // Red
    img.select('SR_B5').multiply(0.8462).add(0.0412).rename('SR_B5'),  // NIR
    img.select('SR_B6').multiply(0.8937).add(0.0254).rename('SR_B6'),  // SWIR1
    img.select('SR_B7').multiply(0.9071).add(0.0172).rename('SR_B7')   // SWIR2
  ]).copyProperties(img, img.propertyNames());
}

var compositeBase = ee.Image(
  HARMONIZE_TO_OLI
    ? harmonizeL7toL8(compositeFilled)
    : compositeFilled
);


// ============================================================
// SECTION 6: SPECTRAL INDICES (14 indices)
//
// Identical formulas to Epoch 3 (L8/9) — band renaming in
// Section 2 ensures full compatibility across epochs.
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
//                RATIO_RS (RED/SWIR1 elevated in saline tidal flats)
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
var NDVI = compositeBase.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

var EVI = compositeBase.expression(
  '2.5 * ((NIR - RED) / (NIR + 6.0 * RED - 7.5 * BLUE + 1.0))',
  { NIR:  compositeBase.select('SR_B5'),
    RED:  compositeBase.select('SR_B4'),
    BLUE: compositeBase.select('SR_B2') }
).rename('EVI');

var SAVI = compositeBase.expression(
  '((NIR - RED) / (NIR + RED + L)) * (1 + L)',
  { NIR: compositeBase.select('SR_B5'),
    RED: compositeBase.select('SR_B4'),
    L: 0.5 }
).rename('SAVI');

// NDMI · Gao 1996 — foliar moisture: high in mangrove, low in dry forest
var NDMI = compositeBase.normalizedDifference(['SR_B5', 'SR_B6']).rename('NDMI');

// LSWI · Xiao et al. 2004 — surface and leaf water content
var LSWI = compositeBase.normalizedDifference(['SR_B5', 'SR_B7']).rename('LSWI');

// CMRI · Jia et al. 2021 — positive over mangrove, negative over mudflat
var CMRI = NDVI.subtract(
  compositeBase.normalizedDifference(['SR_B3', 'SR_B5'])
).rename('CMRI');

// ── Water indices ───────────────────────────────────────────
var NDWI   = compositeBase.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');
var MNDWI  = compositeBase.normalizedDifference(['SR_B3', 'SR_B6']).rename('MNDWI');

// AWEI_sh = BLUE + 2.5×GREEN − 1.5×(NIR + SWIR1) − 0.25×SWIR2
var AWEI_sh = compositeBase.expression(
  'BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2',
  { BLUE:  compositeBase.select('SR_B2'),
    GREEN: compositeBase.select('SR_B3'),
    NIR:   compositeBase.select('SR_B5'),
    SWIR1: compositeBase.select('SR_B6'),
    SWIR2: compositeBase.select('SR_B7') }
).rename('AWEI_sh');

// ── Bare soil indices ───────────────────────────────────────
// BSI = ((SWIR1 + RED) − (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))
var BSI = compositeBase.expression(
  '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))',
  { SWIR1: compositeBase.select('SR_B6'),
    RED:   compositeBase.select('SR_B4'),
    NIR:   compositeBase.select('SR_B5'),
    BLUE:  compositeBase.select('SR_B2') }
).rename('BSI');

// RED/SWIR1 ratio — elevated in saline tidal flats
var RATIO_RS = compositeBase.select('SR_B4')
  .divide(compositeBase.select('SR_B6').add(0.0001))
  .rename('RATIO_RS');


// ============================================================
// SECTION 7: INTRA-ANNUAL TEMPORAL METRICS (8 features)
//
// Temporal metrics are computed from individual scene-level
// images before compositing, which increases robustness
// against SLC-off gaps relative to post-composite statistics.
//
// Interpretation per class:
//   Mangrove:      stable high NDVI year-round (evergreen)
//   Water Bodies:  stable high MNDWI year-round (permanent)
//   Mudflat:       high NDVI variability (seasonal inundation),
//                  high MNDWI amplitude (tidal flooding cycle)
//   Dry Forest:    NDVI low in dry season, higher in wet season
//   Bare Soil:     NDVI ≈ 0 year-round, stable negative MNDWI
//
// Pixels with < 3 valid observations produce unreliable
// temporal metrics and are filled with neutral constants
// (see Section 9) to prevent NoData propagation into sampling.
// ============================================================
var ndviCol = col.map(function(img) {
  return img.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')
    .set('system:time_start', img.get('system:time_start'));
});

var mndwiCol = col.map(function(img) {
  return img.normalizedDifference(['SR_B3', 'SR_B6']).rename('MNDWI')
    .set('system:time_start', img.get('system:time_start'));
});

// NDVI metrics
var NDVI_p10 = ndviCol.reduce(ee.Reducer.percentile([10])).clip(studyArea).rename('NDVI_p10');
var NDVI_p90 = ndviCol.reduce(ee.Reducer.percentile([90])).clip(studyArea).rename('NDVI_p90');
var NDVI_std = ndviCol.reduce(ee.Reducer.stdDev()).clip(studyArea).rename('NDVI_std');
var NDVI_amp = NDVI_p90.subtract(NDVI_p10).rename('NDVI_amp');

// MNDWI metrics
var MNDWI_p10 = mndwiCol.reduce(ee.Reducer.percentile([10])).clip(studyArea).rename('MNDWI_p10');
var MNDWI_p90 = mndwiCol.reduce(ee.Reducer.percentile([90])).clip(studyArea).rename('MNDWI_p90');
var MNDWI_amp = MNDWI_p90.subtract(MNDWI_p10).rename('MNDWI_amp');
var MNDWI_std = mndwiCol.reduce(ee.Reducer.stdDev()).clip(studyArea).rename('MNDWI_std');

// Annual water frequency: permanent water ≈ 1.0 · mudflat ≈ 0.3–0.7 · bare soil ≈ 0.0
var FREQ_AGUA = mndwiCol
  .map(function(img) { return img.gt(0).rename('MNDWI'); })
  .mean().clip(studyArea).rename('FREQ_AGUA');


// ============================================================
// SECTION 8: OBSERVATION QUALITY MASK
//
// Pixels with fewer than 3 valid observations are flagged.
// Their temporal metrics are replaced with neutral fill values
// that do not bias the Random Forest classifier.
//
// Fill values represent spectrally neutral conditions:
//   NDVI_p10  = 0.10  (sparse vegetation at seasonal minimum)
//   NDVI_p90  = 0.50  (moderate vegetation at seasonal maximum)
//   NDVI_std  = 0.05  (low temporal variability)
//   NDVI_amp  = 0.20  (moderate amplitude)
//   MNDWI_p10 = −0.30 (no surface water at seasonal minimum)
//   MNDWI_p90 = 0.10  (slight water signal at seasonal maximum)
//   MNDWI_amp = 0.20  (moderate tidal amplitude)
//   MNDWI_std = 0.08  (low variability)
//   FREQ_AGUA = 0.20  (water present ~20% of the year)
// ============================================================
var sufficientObs = obsDensity.gte(3);

Map.addLayer(sufficientObs,
  {min: 0, max: 1, palette: ['red', 'green']},
  'Pixels with ≥ 3 valid observations', false);

var NDVI_p10_f  = NDVI_p10.updateMask(sufficientObs).unmask(0.10);
var NDVI_p90_f  = NDVI_p90.updateMask(sufficientObs).unmask(0.50);
var NDVI_std_f  = NDVI_std.updateMask(sufficientObs).unmask(0.05);
var NDVI_amp_f  = NDVI_amp.updateMask(sufficientObs).unmask(0.20);
var MNDWI_p10_f = MNDWI_p10.updateMask(sufficientObs).unmask(-0.30);
var MNDWI_p90_f = MNDWI_p90.updateMask(sufficientObs).unmask(0.10);
var MNDWI_amp_f = MNDWI_amp.updateMask(sufficientObs).unmask(0.20);
var MNDWI_std_f = MNDWI_std.updateMask(sufficientObs).unmask(0.08);
var FREQ_AGUA_f = FREQ_AGUA.updateMask(sufficientObs).unmask(0.20);


// ============================================================
// SECTION 9: FULL FEATURE STACK (22 features)
// 6 spectral bands + 8 static indices + 8 temporal metrics
// ============================================================
var features = compositeBase.select(BANDS_ALIAS)
  .addBands(NDVI).addBands(EVI).addBands(SAVI)
  .addBands(NDMI).addBands(LSWI).addBands(CMRI)
  .addBands(NDWI).addBands(MNDWI).addBands(AWEI_sh)
  .addBands(BSI).addBands(RATIO_RS)
  .addBands(NDVI_p10_f).addBands(NDVI_p90_f)
  .addBands(NDVI_std_f).addBands(NDVI_amp_f)
  .addBands(MNDWI_p10_f).addBands(MNDWI_p90_f)
  .addBands(MNDWI_amp_f).addBands(MNDWI_std_f)
  .addBands(FREQ_AGUA_f);

// Explicit JS array avoids server-side resolution issues
// when passing band names to Random Forest inputProperties
var ALL_BANDS = [
  'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',
  'NDVI', 'EVI', 'SAVI', 'NDMI', 'LSWI', 'CMRI',
  'NDWI', 'MNDWI', 'AWEI_sh', 'BSI', 'RATIO_RS',
  'NDVI_p10', 'NDVI_p90', 'NDVI_std', 'NDVI_amp',
  'MNDWI_p10', 'MNDWI_p90', 'MNDWI_amp', 'MNDWI_std',
  'FREQ_AGUA'
];

print('Total features (expected: 22):', ALL_BANDS.length);
print('Feature stack bands:', ALL_BANDS);

// ============================================================
// SECTION 10: TRAINING POLYGONS — 5 CLASSES
//
// Class encoding (TC property):
//   TC = 1  Mangrove       dense evergreen canopy, interior pixels only
//   TC = 2  Water_Bodies   tidal channels, ponds, open water
//   TC = 3  Bare_Soil      dry tidal flats, unvegetated sandy shores
//   TC = 4  Mudflat        tidal mudflat / wet bare soil; seasonally
//                          flooded, sparse or absent vegetation;
//                          high MNDWI amplitude due to tidal cycle
//   TC = 5  Dry_Forest     inland dry-deciduous forest, no tidal influence
//
// Digitization guidelines:
//   · Avoid areas where observation density (Section 4) < 3
//   · Place polygons ≥ 1–2 pixels from class boundaries
//   · Reference image: false-color composite (SR_B5/B4/B3) at zoom 15–16
//   · Mudflat: low-lying intertidal zones, dark wet surface tone,
//              often adjacent to mangrove fringe
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
//   Asset. Use this if you downloaded training_polygons_points_E2_L7.geojson
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
// Change asset suffix: E1_L5 | E2_L7 | E3_L89
//
// var trainingAsset = ee.FeatureCollection(
//   'projects/YOUR_GEE_PROJECT/assets/training_polygons_points_E2_L7'
// );
// var mangrove  = trainingAsset.filter(ee.Filter.eq('TC', 1));
// var water     = trainingAsset.filter(ee.Filter.eq('TC', 2));
// var bareSoil  = trainingAsset.filter(ee.Filter.eq('TC', 3));
// var mudflat   = trainingAsset.filter(ee.Filter.eq('TC', 4));
// var dryForest = trainingAsset.filter(ee.Filter.eq('TC', 5));
// var trainingPolygons = trainingAsset;

// ============================================================
// SECTION 11: POLYGON-BASED TRAIN / TEST SPLIT (70 / 30)
//
// All pixels from the same polygon are assigned exclusively to
// train OR test, preventing spatial autocorrelation leakage.
// Reference: Karasiak et al. (2022), Remote Sensing of Environment
// ============================================================
var polygonsWithRandom = trainingPolygons.randomColumn('random_poly', 42);

var samplesWithPolyRand = features.sampleRegions({
  collection:  polygonsWithRandom,
  properties:  ['TC', 'random_poly'],
  scale:        30,
  geometries:   false,
  tileScale:    4
});

print('Samples per class (TC 1–5):', samplesWithPolyRand.aggregate_histogram('TC'));
print('Total valid samples:', samplesWithPolyRand.size());

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
//   seed = 789            fixed for full reproducibility
// ============================================================
var rfClassifier = ee.Classifier.smileRandomForest({
  numberOfTrees:     500,
  minLeafPopulation:   5,
  bagFraction:         0.7,
  variablesPerSplit:   5,
  seed:              789
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

var waterClosed = rawClassification.eq(2)
  .focal_max({radius: 2, units: 'pixels', iterations: 1})
  .focal_min({radius: 2, units: 'pixels', iterations: 1});
var waterRefined = waterClosed.and(
  rawClassification.eq(2).focal_max({radius: 3, units: 'pixels', iterations: 1})
);
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

Map.addLayer(compositeBase,
  {bands: ['SR_B5', 'SR_B4', 'SR_B3'], min: 0, max: 0.3},
  'False color L7 ' + YEAR);
Map.addLayer(compositeBase,
  {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.3},
  'True color L7 ' + YEAR, false);

Map.addLayer(obsDensity,
  {min: 0, max: 20, palette: ['red', 'yellow', 'green']},
  'Valid observations per pixel (SLC diagnostic)', false);

Map.addLayer(MNDWI,    {min: -0.5, max: 0.5, palette: ['#8B4513','white','#0055E3']}, 'MNDWI', false);
Map.addLayer(AWEI_sh,  {min: -0.3, max: 0.3, palette: ['#8B4513','white','#0055E3']}, 'AWEI_sh (turbid water)', false);
Map.addLayer(BSI,      {min: -0.5, max: 0.5, palette: ['white','#D4A96A']},           'BSI (bare soil)', false);
Map.addLayer(CMRI,     {min: -0.3, max: 0.5, palette: ['white','#006400']},           'CMRI (mangrove-specific)', false);
Map.addLayer(NDVI_amp, {min:    0, max: 0.5, palette: ['white','#FF6600']},           'NDVI amplitude (seasonality)', false);

Map.addLayer(rawClassification,   visParams, 'RF raw L7 '      + YEAR, false);
Map.addLayer(finalClassification, visParams, 'RF final L7 '    + YEAR);

Map.addLayer(SNLMT,      {color: 'yellow'}, 'SNLMT (national sanctuary)', false);
Map.addLayer(bufferZone, {color: 'orange'}, 'Buffer zone',                false);


// ============================================================
// SECTION 16: ACCURACY ASSESSMENT
//
// Metrics: OA (Overall Accuracy), Kappa, PA (Producer's Accuracy),
//          UA (User's Accuracy), F1 (harmonic mean of PA and UA)
// Reference: Congalton & Green (2019); Foody (2002, RSE)
//
// Note: errorMatrix includes a zero row/column for background;
// TC indices (1–5) map directly to list positions 1–5
// ============================================================
var testClassified = testSet.classify(rfClassifier);
var confMatrix     = testClassified.errorMatrix('TC', 'classification');

var OA      = ee.Number(confMatrix.accuracy());
var Kappa   = ee.Number(confMatrix.kappa());
var PA_list = ee.List(confMatrix.producersAccuracy().toList()).flatten();
var UA_list = ee.List(confMatrix.consumersAccuracy().toList()).flatten();

print('════════════════════════════════════════════════════════');
print('ACCURACY ASSESSMENT · YEAR ' + YEAR + ' · L7 ETM+ · 5 CLASSES');
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
var TC_LIST     = [1, 2, 3, 4, 5];
var matrixArray = ee.List(confMatrix.array().toList());

// File 1: Global accuracy metrics
Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, {
      Year:         YEAR,
      Epoch:        1,
      Sensor:       'L7_ETM+',
      N_Classes:    5,
      Split_Method: 'Polygon_70-30',
      N_Train:      trainSet.size(),
      N_Test:       testSet.size(),
      OA:           OA,
      Kappa:        Kappa
    })
  ]),
  description: 'Acc_Global_5C_E1_L7_' + YEAR,
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
      Epoch:  1,
      Sensor: 'L7_ETM+',
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
  description: 'Acc_PerClass_5C_E1_L7_' + YEAR,
  folder:      'GEE_Manglares_Tumbes',
  fileFormat:  'CSV'
});

// File 3: Full confusion matrix
var confMatrixTable = ee.FeatureCollection(
  ee.List(TC_LIST).map(function(claseTC) {
    var tcNum  = ee.Number(claseTC);
    var rowRaw = ee.List(matrixArray.get(tcNum));
    return ee.Feature(null, {
      Year:           YEAR,
      Epoch:          1,
      Sensor:         'L7_ETM+',
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
  description: 'Acc_ConfMatrix_5C_E1_L7_' + YEAR,
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
      Epoch:   1,
      Sensor:  'L7_ETM+',
      Zone:    zoneName,
      Class:   CLASS_NAMES.get(num),
      TC:      ee.Number(item.get('class')).toInt(),
      Area_ha: ee.Number(item.get('sum')).multiply(100).round().divide(100)
    });
  });
}

var areaTable = ee.FeatureCollection(areaByZone(finalClassification, studyArea,  'Study_Area'))
  .merge(ee.FeatureCollection(areaByZone(finalClassification, SNLMT,      'SNLMT')))
  .merge(ee.FeatureCollection(areaByZone(finalClassification, bufferZone, 'Buffer_Zone')));

print('Area by class and management zone (ha):', areaTable);

Export.table.toDrive({
  collection:  areaTable,
  description: 'Areas_5C_ByZone_E1_L7_' + YEAR,
  folder:      'GEE_Manglares_Tumbes',
  fileFormat:  'CSV'
});


// ============================================================
// SECTION 20: EXPORT — FINAL CLASSIFIED MAP (GeoTIFF)
// Projection: UTM Zone 17S (EPSG:32717) · Spatial resolution: 30 m
// ============================================================
Export.image.toDrive({
  image:          finalClassification.toByte(),
  description:    'Mangrove_5C_E2_L7_' + YEAR,
  folder:         'GEE_Manglares_Tumbes',
  fileNamePrefix: 'mangrove_L7_5c_' + YEAR,
  region:         studyArea.geometry(),
  scale:          30,
  crs:            'EPSG:32717',
  fileFormat:     'GeoTIFF',
  maxPixels:      1e13
});
