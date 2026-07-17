// ============================================================
// Supervised Land Cover Classification · Epoch 1 · 1993/1995–1999
// Tumbes Mangrove System, Northwestern Peru
// Sensor: Landsat 5 TM (USGS Collection 2 SR)
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
//           Ulises Francisco Giraldo Malca, Pablito Marcelo López Serrano]
// GEE Project: ee-brayansotoquispe
// Last updated: 2024
//
// ============================================================


// ============================================================
// LANDSAT 5 TM — KEY DIFFERENCES FROM LANDSAT 8/9 OLI
// ============================================================
//
// Band equivalencies (TM → OLI naming convention):
//   TM B1 (0.45–0.52 µm) → Blue   (≈ OLI B2)
//   TM B2 (0.52–0.60 µm) → Green  (≈ OLI B3)
//   TM B3 (0.63–0.69 µm) → Red    (≈ OLI B4)
//   TM B4 (0.76–0.90 µm) → NIR    (≈ OLI B5)
//   TM B5 (1.55–1.75 µm) → SWIR1  (≈ OLI B6)
//   TM B7 (2.08–2.35 µm) → SWIR2  (≈ OLI B7)
//
// Important notes:
//   · TM has no coastal aerosol band (OLI B1 equivalent)
//   · TM B6 is the thermal band — excluded from spectral stack
//   · AWEI_nsh cannot be computed (requires coastal band);
//     AWEI_sh is used instead (requires only GREEN, NIR, SWIR1, SWIR2)
//   · Collection 2 Level-2 scale factors are identical to L8/9
//   · Spatial resolution: 30 m (same as L8/9)
//   · Last fully reliable year: 2011 (gyroscope failure in 2012)
//   · Revisit cycle: 16 days (~22 potential scenes per year)
// ============================================================

// ============================================================
// SECTION 0: MAIN PARAMETER
//
// HOW TO USE:
//
// ── Option A: Single-year composite ─────────────────────────
//    Set all three variables to the same year.
//    Use this when the target year has enough cloud-free scenes
//    (typically ≥ 5 scenes after cloud filtering).
//
//    Example for 2006:
//      var YEAR       = 1996;
//      var YEAR_START = 1996;
//      var YEAR_END   = 1996;
//
// ── Option B: Multi-year composite ──────────────────────────
//    Set YEAR_START and YEAR_END to bracket the window.
//    YEAR is the label used in export filenames and does NOT
//    need to equal YEAR_START — set it to the ecological year
//    of interest (e.g. the central year of the window).
//    Use this when the target year has few usable scenes due to
//    persistent cloud cover. A ±1 year window (3 years total) is recommended
//    as a first attempt; expand to ±2 years only if needed.
//
//    Example — 3-year window centered on 1993:
//      var YEAR       = 1995;   // label for exports
//      var YEAR_START = 1993;   // composite starts here
//      var YEAR_END   = 1995;   // composite ends here
//
//    Example — 2-year window (e.g. before/after gap):
//      var YEAR       = 1997;
//      var YEAR_START = 1996;
//      var YEAR_END   = 1997;
//
// RULE OF THUMB:
//   Check col.size() in the GEE Console after running.
//   < 3 scenes  → widen the window or raise CLOUD_COVER to 80%
//   3–5 scenes  → acceptable; multi-year composite recommended
//   ≥ 5 scenes  → single-year composite is reliable
// ============================================================
var YEAR       = 1995;   // Target / label year (used in export filenames)
var YEAR_START = 1993;   // Start of composite window  ← set = YEAR for single-year
var YEAR_END   = 1995;   // End of composite window    ← set = YEAR for single-year
// ============================================================
// SECTION 1: STUDY AREA BOUNDARIES
// ============================================================
var areaManglar = ee.FeatureCollection('projects/your_username/assets/AREA_DEF');
var SNLMT = ee.FeatureCollection('projects/your_username/assets/SNLMT');
var Zona_Amortiguamiento = ee.FeatureCollection('projects/your_username/assets/Zona_Amortiguamiento');

Map.centerObject(studyArea, 10);

// ============================================================
// SECTION 2: SPECTRAL BANDS — LANDSAT 5 TM
// Native TM band names used for collection filtering;
// indices computed using original TM band designations
// ============================================================
var BANDS_TM = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];

// ============================================================
// SECTION 3: PREPROCESSING — CLOUD MASKING AND SCALING
// Collection 2 Level-2 scale factors identical to L8/9:
//   scale × 0.0000275, offset −0.2
// QA_PIXEL bit flags identical to L8/9 in Collection 2
// ============================================================
function maskCloudsL5(img) {
  var qa = img.select('QA_PIXEL');
  return img.updateMask(
    qa.bitwiseAnd(1 << 3).eq(0)          // no cloud shadow
      .and(qa.bitwiseAnd(1 << 5).eq(0))  // no cloud
  );
}

function applyScaleL5(img) {
  var sr = img.select(BANDS_TM).multiply(0.0000275).add(-0.2);
  return sr
    .updateMask(
      sr.reduce(ee.Reducer.min()).gt(-0.1)
        .and(sr.reduce(ee.Reducer.max()).lt(1.1))
    )
    .copyProperties(img, ['system:time_start', 'system:index']);
}

function preprocessL5(img) {
  return applyScaleL5(maskCloudsL5(img));
}

// ============================================================
// SECTION 4: ANNUAL LANDSAT 5 IMAGE COLLECTION
//
// Tumbes has a rainy season (Jan–Apr, high cloud cover) and a
// dry season (May–Dec, better image quality). Temporal metrics
// in Section 7 capture this intra-annual seasonality.
// If fewer than 3 scenes are available, consider raising the
// cloud cover threshold to 80% or using an adjacent year.
// ============================================================

var startDate = ee.Date.fromYMD(YEAR_START, 1, 1);
var endDate   = ee.Date.fromYMD(YEAR_END,   12, 31);

var filters = ee.Filter.and(
  ee.Filter.date(startDate, endDate),
  ee.Filter.bounds(studyArea),
  ee.Filter.lt('CLOUD_COVER', 80)
);

var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filter(filters)
  .map(preprocessL5);

print('Scenes used (' + YEAR_START + '–' + YEAR_END + '):', col.size());
print('⚠ If count < 5, widen the year window or raise CLOUD_COVER threshold');

// ============================================================
// SECTION 5: ANNUAL MEDIAN COMPOSITE
// ============================================================
var composite = col.median().clip(studyArea);

// ============================================================
// SECTION 6: SPECTRAL HARMONIZATION TM → OLI
//
// Systematic inter-sensor differences between L5 TM and
// L8/9 OLI are corrected using the OLS regression coefficients
// of Roy et al. (2016): OLI = a + b × ETM+ (surface reflectance)
//
// The same coefficients apply to both TM and ETM+, as both
// sensors share nearly identical spectral response functions
// relative to OLI.
//
// Coefficients (Roy et al. 2016, RSE 185:57-70, Table 2):
//   Blue  (B1→B2): slope = 0.8474, intercept = 0.0003
//   Green (B2→B3): slope = 0.8483, intercept = 0.0088
//   Red   (B3→B4): slope = 0.9047, intercept = 0.0061
//   NIR   (B4→B5): slope = 0.8462, intercept = 0.0412
//   SWIR1 (B5→B6): slope = 0.8937, intercept = 0.0254
//   SWIR2 (B7→B7): slope = 0.9071, intercept = 0.0172
//
// Applied here to ensure spectral consistency when comparing
// Epoch 1 (L5 TM) maps against Epoch 3 (L8/9 OLI) maps.
// For single-sensor classifications it is optional.
//
// Reference: Roy et al. (2016). Remote Sensing of Environment,
//            185, 57–70.
// ============================================================
var HARMONIZE_TO_OLI = true;

function harmonizeTMtoOLI(img) {
  img = ee.Image(img);
  return ee.Image.cat([
    img.select('SR_B1').multiply(0.8474).add(0.0003).rename('SR_B1'),  // Blue
    img.select('SR_B2').multiply(0.8483).add(0.0088).rename('SR_B2'),  // Green
    img.select('SR_B3').multiply(0.9047).add(0.0061).rename('SR_B3'),  // Red
    img.select('SR_B4').multiply(0.8462).add(0.0412).rename('SR_B4'),  // NIR
    img.select('SR_B5').multiply(0.8937).add(0.0254).rename('SR_B5'),  // SWIR1
    img.select('SR_B7').multiply(0.9071).add(0.0172).rename('SR_B7')   // SWIR2
  ]).copyProperties(img, img.propertyNames());
}

// Composite
var compositeBase = ee.Image(
  HARMONIZE_TO_OLI
    ? harmonizeTMtoOLI(composite)
    : composite
);

// ============================================================
// SECTION 7: SPECTRAL INDICES (14 indices)
//
// All indices use native TM band designations:
//   Blue  = SR_B1 | Green = SR_B2 | Red   = SR_B3
//   NIR   = SR_B4 | SWIR1 = SR_B5 | SWIR2 = SR_B7
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
// NDVI · Rouse et al. 1974
var NDVI = compositeBase.normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI');
//                                              NIR       Red (TM)

// EVI · Huete et al. 2002
var EVI = compositeBase.expression(
  '2.5 * ((NIR - RED) / (NIR + 6.0 * RED - 7.5 * BLUE + 1.0))',
  { NIR:  compositeBase.select('SR_B4'),  // TM NIR
    RED:  compositeBase.select('SR_B3'),  // TM Red
    BLUE: compositeBase.select('SR_B1')   // TM Blue
  }
).rename('EVI');

// SAVI · Huete 1988 — soil-adjusted index (L=0.5) for sparse vegetation
var SAVI = compositeBase.expression(
  '((NIR - RED) / (NIR + RED + L)) * (1 + L)',
  { NIR: compositeBase.select('SR_B4'),
    RED: compositeBase.select('SR_B3'),
    L: 0.5 }
).rename('SAVI');

// NDMI · Gao 1996 — foliar moisture: high in mangrove, low in dry forest
var NDMI = compositeBase.normalizedDifference(['SR_B4', 'SR_B5']).rename('NDMI');
//                                              NIR       SWIR1 (TM)

// LSWI · Xiao et al. 2004 — surface and leaf water content
var LSWI = compositeBase.normalizedDifference(['SR_B4', 'SR_B7']).rename('LSWI');
//                                              NIR       SWIR2 (TM)

// CMRI · Jia et al. 2021 — mangrove-specific index
// CMRI = NDVI − NDWI; positive over mangrove, negative over mudflat
var NDWI_tmp = compositeBase.normalizedDifference(['SR_B2', 'SR_B4']);
var CMRI = NDVI.subtract(NDWI_tmp).rename('CMRI');

// ── Water indices ───────────────────────────────────────────
// NDWI · McFeeters 1996
var NDWI = compositeBase.normalizedDifference(['SR_B2', 'SR_B4']).rename('NDWI');
//                                              Green     NIR (TM)

// MNDWI · Xu 2006 — suppresses vegetation and soil better than NDWI
var MNDWI = compositeBase.normalizedDifference(['SR_B2', 'SR_B5']).rename('MNDWI');
//                                               Green     SWIR1 (TM)

// AWEI_sh · Feyisa et al. 2014 — optimized for turbid and shaded water
// Does not require the coastal aerosol band → valid for TM
// AWEI_sh = BLUE + 2.5×GREEN − 1.5×(NIR + SWIR1) − 0.25×SWIR2
var AWEI_sh = compositeBase.expression(
  'BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2',
  { BLUE:  compositeBase.select('SR_B1'),  // TM B1
    GREEN: compositeBase.select('SR_B2'),  // TM B2
    NIR:   compositeBase.select('SR_B4'),  // TM B4
    SWIR1: compositeBase.select('SR_B5'),  // TM B5
    SWIR2: compositeBase.select('SR_B7')   // TM B7
  }
).rename('AWEI_sh');

// ── Bare soil indices ───────────────────────────────────────
// BSI · Rikimaru et al. 2002
// BSI = ((SWIR1 + RED) − (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))
var BSI = compositeBase.expression(
  '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))',
  { SWIR1: compositeBase.select('SR_B5'),
    RED:   compositeBase.select('SR_B3'),
    NIR:   compositeBase.select('SR_B4'),
    BLUE:  compositeBase.select('SR_B1')
  }
).rename('BSI');

// RED/SWIR1 ratio — elevated in saline tidal flats and exposed mudflats
var RATIO_RS = compositeBase.select('SR_B3')
  .divide(compositeBase.select('SR_B5').add(0.0001))
  .rename('RATIO_RS');
//               Red (TM)              SWIR1 (TM)


// ============================================================
// SECTION 8: INTRA-ANNUAL TEMPORAL METRICS (8 features)
//
// Computed from individual scenes before compositing.
// With a 16-day revisit cycle, ~22 potential scenes exist per
// year; typically 5–15 usable scenes after cloud filtering.
//
// Interpretation per class:
//   Mangrove:      stable high NDVI year-round (evergreen)
//   Water Bodies:  stable high MNDWI year-round (permanent)
//   Mudflat:       high NDVI variability (seasonal inundation),
//                  high MNDWI amplitude (tidal flooding cycle)
//   Dry Forest:    NDVI low in dry season, higher in wet season
//   Bare Soil:     NDVI ≈ 0 year-round, stable negative MNDWI
// ============================================================
var ndviCol = col.map(function(img) {
  return img.normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI')
    .set('system:time_start', img.get('system:time_start'));
});

var mndwiCol = col.map(function(img) {
  return img.normalizedDifference(['SR_B2', 'SR_B5']).rename('MNDWI')
    .set('system:time_start', img.get('system:time_start'));
});

// NDVI metrics
var NDVI_p10 = ndviCol.reduce(ee.Reducer.percentile([10])).clip(studyArea).rename('NDVI_p10');
var NDVI_p90 = ndviCol.reduce(ee.Reducer.percentile([90])).clip(studyArea).rename('NDVI_p90');
var NDVI_std = ndviCol.reduce(ee.Reducer.stdDev()).clip(studyArea).rename('NDVI_std');
// High amplitude → seasonally deciduous dry forest or tidally flooded mudflat
// Low amplitude  → evergreen mangrove or permanently unvegetated bare soil
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
// SECTION 9: FULL FEATURE STACK (22 features)
// 6 TM spectral bands + 8 static indices + 8 temporal metrics
//
// Band count matches Epoch 3 (L8/9) stack exactly, ensuring
// a methodologically consistent feature space across all epochs
// ============================================================
var features = compositeBase.select(BANDS_TM)
  .addBands(NDVI).addBands(EVI).addBands(SAVI)
  .addBands(NDMI).addBands(LSWI).addBands(CMRI)
  .addBands(NDWI).addBands(MNDWI).addBands(AWEI_sh)
  .addBands(BSI).addBands(RATIO_RS)
  .addBands(NDVI_p10).addBands(NDVI_p90)
  .addBands(NDVI_std).addBands(NDVI_amp)
  .addBands(MNDWI_p10).addBands(MNDWI_p90)
  .addBands(MNDWI_amp).addBands(MNDWI_std)
  .addBands(FREQ_AGUA);

// Explicit JS array avoids server-side resolution issues
// when passing band names to Random Forest inputProperties
var ALL_BANDS = [
  'SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7',
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
// Digitization guidelines for TM imagery:
//   · Use false-color composite (SR_B4/B3/B2 = NIR/Red/Green) at zoom 15–16
//   · TM has slightly higher radiometric noise than OLI —
//     use more conservative (interior) polygon placement
//   · For historical years (pre-2000), use the L5 false-color
//     image itself as the primary reference; supplement with
//     available INRENA / SERNANP / ANA land cover maps
//   · Mudflat: dark wet tidal surface, often adjacent to
//     mangrove fringe; differentiate from dry bare soil by
//     tone and proximity to tidal channels
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
//   Asset. Use this if you downloaded training_polygons_points_E1_L5.geojson
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
//   'projects/YOUR_GEE_PROJECT/assets/training_polygons_points_E1_L5'
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
}).filter(ee.Filter.notNull(ALL_BANDS));

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
// Hyperparameters are identical across all epochs (L5, L7, L8/9)
// to ensure methodological consistency and cross-epoch comparability.
//
//   numberOfTrees = 500   stable variable importance estimates
//                         (Belgiu & Drăguţ 2016, ISPRS J. Photogramm.)
//   minLeafPopulation = 5 prevents overfitting on small classes
//                         (Gislason et al. 2006, Pattern Recognit. Lett.)
//   bagFraction = 0.7     standard for land cover mapping
//   variablesPerSplit = 5 ≈ sqrt(22 features); Breiman (2001)
//   seed = 123            fixed for full reproducibility (Epoch 1 L5)
// ============================================================
var rfClassifier = ee.Classifier.smileRandomForest({
  numberOfTrees:     500,
  minLeafPopulation:   5,
  bagFraction:         0.7,
  variablesPerSplit:   5,
  seed:              123
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

// False color TM: B4/B3/B2 = NIR/Red/Green (equivalent to B5/B4/B3 in L8/9)
Map.addLayer(compositeBase,
  {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.3},
  'False color TM ' + YEAR);
Map.addLayer(compositeBase,
  {bands: ['SR_B3', 'SR_B2', 'SR_B1'], min: 0, max: 0.3},
  'True color TM ' + YEAR, false);

Map.addLayer(MNDWI,     {min: -0.5, max: 0.5, palette: ['#8B4513','white','#0055E3']}, 'MNDWI (TM)', false);
Map.addLayer(AWEI_sh,   {min: -0.3, max: 0.3, palette: ['#8B4513','white','#0055E3']}, 'AWEI_sh turbid water (TM)', false);
Map.addLayer(BSI,       {min: -0.5, max: 0.5, palette: ['white','#D4A96A']},           'BSI bare soil (TM)', false);
Map.addLayer(CMRI,      {min: -0.3, max: 0.5, palette: ['white','#006400']},           'CMRI mangrove-specific (TM)', false);
Map.addLayer(NDVI_amp,  {min:    0, max: 0.5, palette: ['white','#FF6600']},           'NDVI amplitude seasonality (TM)', false);
Map.addLayer(MNDWI_amp, {min:    0, max: 0.4, palette: ['white','#0055E3']},           'MNDWI amplitude tidal variability (TM)', false);
Map.addLayer(FREQ_AGUA, {min:    0, max:   1, palette: ['white','#0055E3']},           'Annual water frequency (TM)', false);

Map.addLayer(rawClassification,   visParams, 'RF raw TM '      + YEAR, false);
Map.addLayer(smoothed1,           visParams, 'RF smoothed (TM)',        false);
Map.addLayer(finalClassification, visParams, 'RF final TM '    + YEAR);

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
print('ACCURACY ASSESSMENT · YEAR ' + YEAR + ' · L5 TM · 5 CLASSES');
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
      Sensor:       'L5_TM',
      N_Classes:    5,
      Split_Method: 'Polygon_70-30',
      N_Train:      trainSet.size(),
      N_Test:       testSet.size(),
      OA:           OA,
      Kappa:        Kappa
    })
  ]),
  description: 'Acc_Global_5C_E1_L5_' + YEAR,
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
      Sensor: 'L5_TM',
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
  description: 'Acc_PerClass_5C_E1_L5_' + YEAR,
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
      Sensor:         'L5_TM',
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
  description: 'Acc_ConfMatrix_5C_E1_L5_' + YEAR,
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
      Sensor:  'L5_TM',
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
  description: 'Areas_5C_ByZone_E1_L5_' + YEAR,
  folder:      'GEE_Manglares_Tumbes',
  fileFormat:  'CSV'
});


// ============================================================
// SECTION 20: EXPORT — FINAL CLASSIFIED MAP (GeoTIFF)
// Projection: UTM Zone 17S (EPSG:32717) · Spatial resolution: 30 m
// ============================================================
Export.image.toDrive({
  image:          finalClassification.toByte(),
  description:    'Mangrove_5C_E1_L5_' + YEAR,
  folder:         'GEE_Manglares_Tumbes',
  fileNamePrefix: 'mangrove_L5_5c_' + YEAR,
  region:         studyArea.geometry(),
  scale:          30,
  crs:            'EPSG:32717',
  fileFormat:     'GeoTIFF',
  maxPixels:      1e13
});
