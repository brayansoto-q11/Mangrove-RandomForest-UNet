# Google Earth Engine (GEE) Scripts

This folder contains the Google Earth Engine (GEE) scripts used to classify mangrove forests using the Random Forest algorithm and Landsat imagery.

The workflow was developed to generate Landsat-based land cover classifications for different temporal periods using Landsat 5 TM, Landsat 7 ETM+, and Landsat 8/9 OLI imagery.

---

## Folder structure

```
gee/
│
├── asset_data/
│   ├── AREA_DEF.geojson
│   ├── SNLMT.geojson
│   └── Zona_Amortiguamiento.geojson
│
├── Landsat5/
│   ├── RF_Landsat5.js
│   └── training_polygons_points.geojson
│
├── Landsat7/
│   ├── RF_Landsat7.js
│   └── training_polygons_points.geojson
│
└── Landsat8-9/
    ├── RF_Landsat8_9.js
    └── training_polygons_points.geojson
```

---

# Description

Each Landsat folder contains:

- **classification.js**: Google Earth Engine scripts implementing the Random Forest classification workflow.
- **training_polygons_points.geojson**: Training polygons and points manually digitized for the corresponding Landsat temporal period.

The **asset_data** folder contains the geographic boundaries required by all scripts:

- `AREA_DEF.geojson` — complete study area.
- `SNLMT.geojson` — Santuario Nacional Los Manglares de Tumbes (protected core area).
- `Zona_Amortiguamiento.geojson` — buffer zone surrounding the protected area.

---

# Training polygons — important note

The training polygons and points included in each Landsat folder were manually digitized in the Google Earth Engine Code Editor using false-color median composites as visual reference.

These samples represent land cover conditions for a reference period and are provided to facilitate reproducibility of the classification workflow.

They should not be considered a universal training dataset.

For years with substantial land cover changes, unusual climatic conditions, or major disturbances, the training samples should be reviewed and updated before classification.

Particularly important periods include strong ENSO events:

- 1997–1998
- 2016–2017

Users can:

1. Add new polygons directly in the GEE Code Editor.
2. Modify the GeoJSON files in GIS software (e.g., QGIS).
3. Re-upload the updated samples as GEE Assets.

---

# Training samples by Landsat period

| File | Reference period | Sensor | Recommended use |
|------|-----------------|--------|----------------|
| Landsat5/training_polygons_points_E1_L5.geojson | Epoch 1 | Landsat 5 TM | 1993/1995–1999 |
| Landsat7/training_polygons_points_E2_L7.geojson | Epoch 2 | Landsat 7 ETM+ | 2000–2012 |
| Landsat8-9/training_polygons_points_E2_L89.geojson | Epoch 3 | Landsat 8/9 OLI | 2013–2024 |

---

# Temporal coverage

The classification periods were defined according to Landsat data availability:

| Sensor | Period | Processing strategy |
|--------|--------|--------------------|
| Landsat 5 TM | 1993–1999 | Annual or multi-year median composites depending on image availability. Multi-year composites may be required for years with limited cloud-free observations. |
| Landsat 7 ETM+ | 2000–2012 | Annual median composites. Landsat 7 SLC-off effects were considered after 2003. |
| Landsat 8/9 OLI | 2013–2024 | Annual median composites combining Landsat 8 and Landsat 9 imagery. |

---

# Required Assets

Before running any script, upload the following files as Google Earth Engine Assets:

1. `AREA_DEF.geojson`
2. `SNLMT.geojson`
3. `Zona_Amortiguamiento.geojson`
4. `training_polygons_points.geojson` from the corresponding Landsat folder.

After uploading the files, update the asset paths at the beginning of each script:

```javascript
var studyArea = ee.FeatureCollection(
  'projects/YOUR_GEE_PROJECT/assets/AREA_DEF'
);
```

Replace:

```
YOUR_GEE_PROJECT
```

with your own Google Earth Engine project ID.

---

# Workflow

## Option A — Use provided training polygons (recommended for reproduction)

1. Upload all files from **asset_data** as GEE Assets.
2. Upload the corresponding `training_polygons_points.geojson`.
3. In **SECTION 10** of the script:
   - Uncomment the asset-loading block.
   - Comment out the original digitized imports.
4. Update all Asset paths.
5. Define the processing period in **SECTION 0**.
6. Run the script.
7. Submit the export tasks.

---

## Option B — Digitize new training polygons

Users can create their own training samples:

1. Open the script in the GEE Code Editor.
2. Draw polygons or points using the Geometry Imports panel.
3. Name each class exactly as required:

```
Mangrove
Water_Bodies
Bare_soil
Mudflat
Dry_forest
```

4. Keep the default import configuration.
5. Run the classification workflow.

---

# Single-year vs multi-year composites

The scripts support both single-year and multi-year median composites (ONLY FOR LANDSAT 5).

These options are controlled in **SECTION 0**:

```javascript
// Single-year composite
var YEAR       = 1996;
var YEAR_START = 1996;
var YEAR_END   = 1996;


// Multi-year composite
var YEAR       = 1995;
var YEAR_START = 1993;
var YEAR_END   = 1995;
```

Multi-year composites are recommended when the number of available cloud-free observations is insufficient.

After execution, check:

```javascript
col.size()
```

in the GEE Console:

- **< 3 scenes** → expand the temporal window.
- **3–5 scenes** → multi-year composite recommended.
- **≥ 5 scenes** → single-year composite generally acceptable.

---

# Software

- Google Earth Engine JavaScript API.
- Random Forest classifier:
  - `ee.Classifier.smileRandomForest`
- Landsat Collection 2 Level-2 Surface Reflectance products (USGS).
