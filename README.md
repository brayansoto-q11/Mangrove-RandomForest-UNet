# Mangrove Mapping using Random Forest and U-Net

This repository contains the complete workflow developed for mangrove ecosystem mapping in the Tumbes region (Peru) using Landsat imagery, Google Earth Engine (GEE), Random Forest classification, and U-Net deep learning.

---

# Repository structure

```
Mangrove-RandomForest-UNet/
│
├── README.md
├── LICENSE
├── CITATION.cff
│
├── gee/
│   ├── README.md
│   ├── asset_data/
│   ├── Landsat5/
│   ├── Landsat7/
│   └── Landsat8-9/
│
├── unet/
    ├── README.md
    ├── UNet_training.ipynb
    ├── requirements.txt
    ├── training_data/
    ├── model_weights/
    └── temporal_models/
```

---

# Project overview

The workflow consists of two complementary stages.

## Stage 1 — Random Forest (Google Earth Engine)

The first stage uses Google Earth Engine to generate annual Random Forest land-cover classifications from Landsat imagery.

The workflow includes:

- Image preprocessing
- Cloud masking
- Spectral index calculation
- Random Forest training
- Classification
- Accuracy assessment

Detailed instructions are available in:

```
gee/README.md
```

---

## Stage 2 — U-Net Deep Learning

The second stage uses the Random Forest classifications as reference labels to train a U-Net semantic segmentation model.

The U-Net workflow includes:

- Model weight
- Trainign data
- Model evaluation
- Export of trained weights
- Temporal prediction using the trained model

Detailed instructions are available in:

```
unet/README.md
```

---

# Study area

Mangrove ecosystem of Tumbes, Peru.

---

# Software

- Google Earth Engine
- Python
- TensorFlow / Keras
- NumPy
- Rasterio
- GDAL
- Scikit-learn

---

# Citation

If you use this repository in your research, please cite both the accompanying publication and this repository.

See:

```
CITATION.cff
```

for citation information.

---

# License

This project is distributed under the MIT License.

See:

```
LICENSE
```
