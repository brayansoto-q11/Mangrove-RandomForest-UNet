# U-Net Deep Learning Workflow

This folder contains the complete workflow used to train and apply a U-Net semantic segmentation model for aquaculture pond mapping in the Tumbes mangrove ecosystem (Peru).

---

# Folder structure

```
unet/
│
├── README.md
├── UNet_Aquaculture_Tumbes.ipynb
├── UNet_Inference_TimeSeries_Tumbes.ipynb
│
├── training_data/
│   ├── AREA_DEF.geojson
│   ├── download_training_image.js
│   ├── pozas_acuicultura_1995.geojson
│   ├── pozas_acuicultura_2003.geojson
│   └── pozas_acuicultura_2024.geojson
│
└── model_weight/
    ├── mejor_modelo_unet.pth
    └── metadata_unet.pkl
```

---

# Workflow

The U-Net workflow consists of two independent stages:

1. Model training
2. Time-series inference

---

# 1. Model training

The notebook

```
UNet_Aquaculture_Tumbes.ipynb
```

is used to train the U-Net model.

## Required training data

The **training_data** folder contains all datasets required for model training.

### Study area

```
AREA_DEF.geojson
```

Defines the study area used to generate the training imagery.

### Training polygons

Three manually digitized polygon datasets are provided:

- `pozas_acuicultura_1995.geojson`
- `pozas_acuicultura_2003.geojson`
- `pozas_acuicultura_2024.geojson`

These polygons represent aquaculture ponds for three reference years and are used to generate the training masks.

### Training imagery

Training images are **not included** in this repository.

Instead, they can be reproduced using:

```
download_training_image.js
```

This Google Earth Engine script generates the spectral-index images used as input for the U-Net model.

---

# 2. Time-series inference

After training, the notebook

```
UNet_Inference_TimeSeries_Tumbes.ipynb
```

is used to apply the trained model to Landsat imagery from multiple years (applied for all the imagen per year download from download_training_image.js).

The notebook performs inference for the complete temporal series.

---

# Required model files

The **model_weight** folder contains the trained model and its metadata:

- `mejor_modelo_unet.pth`
- `metadata_unet.pkl`

Inside the inference notebook there are sections where the paths to these files must be specified.

Update these paths according to the location of your local repository before running the notebook.

---

# Running the workflow

## Step 1

Generate the training images using

```
download_training_image.js
```

in Google Earth Engine.

---

## Step 2

Run

```
UNet_Aquaculture_Tumbes_RSASE.ipynb
```

to train the model.

The notebook will generate:

- trained model weights
- training history
- performance metrics

---

## Step 3

Run

```
UNet_Inference_TimeSeries_Tumbes.ipynb
```

to predict aquaculture ponds for the complete Landsat time series.

Before execution, update the file paths to:

- `mejor_modelo_unet.pth`
- `metadata_unet.pkl`

---

# Important note

The provided model was trained specifically for the mangrove ecosystem of Tumbes, Peru.

Although the workflow can be adapted to other study areas, the trained weights should **not** be assumed to generalize to different geographic regions without additional training and validation.

Users applying this workflow outside the Tumbes region are encouraged to create new training samples and retrain or fine-tune the model.

The released model represents the first version developed for this study and should be considered a baseline implementation. Future improvements in training data, model architecture, or transfer learning strategies may further enhance its performance.

---

# Software

- Python 3
- PyTorch
- Rasterio
- NumPy
- OpenCV
- GDAL
- Scikit-learn
- Matplotlib