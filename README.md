
# Project Overview

This repository contains various folders and scripts related to the Master thesis "Northern Boreal Alluvial Meadows (Svamängar) spatio-temporal changes detection using satellite data, AI and Google Earth Engine". Below is a brief description of each folder:

## Folders Description

**lulc_rf_models**  
scripts for Random Forest models for land use and land cover (LULC) classification.
   - `rfmodel_approach_1a.ipynb` & `rfmodel_approach_1b.ipynb` approach 1 in the thesis report. 
   - `rfmodel_approach_2a.ipynb`, `rfmodel_approach_1c.ipynb` & `rfmodel_approach_1d.ipynb` refers to approach 2 in the thesis report. 
   - `rfmodel_approach_3.ipynb` refers to same approach as approach 2d, but with the use of ESA world cover as LULC product for the reference data instead of Google Dynamiv. 
   - `labeling.ipynb`, refers to a notebook used to label some data manually
   - `lulc_ts_gw.ipynb` refers to a notebook to generate time series directly from the lulc maps product.
   - `csv_charging.ipynb` is for the display of the data
   - `create_composite_l8.py` & `simple_random_forest.py` are not used atm but will be useful for deployment of the method. 
   - `spatial_changes_detection_ml.ipynb` is a notebook to detect and display spatial changes btw two images for lulc method.

   As some models were built to see the impact of the training data (temporal agreggation and reference data) on the model performance, all models were built simply, using GEE RF model & platform to handle huge amount of multitemporal & multispectral satellite data. Some preliminary results on hyperparameters optimization can be found in rf_results.xlsx. 

**data_fusion_landsat_modis/**  
   - landsat_modis folder for data fusion to generate new images with same resolution as Landsat using Modis images for prediction. 
   - linear_transformation_l7_l8_s2.ipynb for radiometric comparison btw landsat and sentinel 2 sensors.

**global-canopy-height-model-main**  
   models from https://github.com/langnico/global-canopy-height-model. Nico Lang DL model to predict vegetation height with S2 images, and trained with LiDAR data. The adaptation finally didn't work. 

**wetlands_detection/**
scripts for **wetlands detection** using Sentinel-2 optical data:
- `water_time_series_l8.ipynb`: Analyzes water time series for Landsat 8.
- `water_time_series_s2.ipynb`: Analyzes water time series for Sentinel-2.
- `water_time_series_month_s2.ipynb`: Monthly water time series analysis for Sentinel-2.
- `wetlands_mask_th.py`: Generates wetland masks using thresholding techniques.
- `wetlands_uns_clustering.py`: Performs unsupervised clustering for wetlands identification.

**wetlands_detection_radar/**
scripts for **wetlands detection** using Sentinel-1 radar data:
- `radar_preprocessing.py`: Prepares Sentinel-1 radar data.
- `radar_change_detection.ipynb`: Detects changes using radar-based analysis.
- `radar_wetlands_mask.py`: Creates wetland masks using radar data.
- `radar_unsupervised_clustering.py`: Performs unsupervised clustering on radar data for wetland detection.

**spectral_indexing/**
scripts for spectral index analysis, focusing on time series analysis and correlation:
- `ACP_indexing.ipynb`: Principal Component Analysis (PCA) for spectral indices.
- `dem_ndvi_density_correlation.ipynb`: Analyzes the correlation between NDVI and elevation (DEM).
- `landsat_indices_time_series.ipynb`: Processes spectral indices for Landsat data, and generate some time series of it, as well as time series from climate data.
- `s2_indices_time_series.ipynb`: Processes spectral indices for Sentinel-2 data, and generate some time series of it,.
- `s2_indices_ts_correlation.ipynb`: Correlates spectral indices over Sentinel-2 time series to estimate the best vegetation spectral index to use.
- `changes_btw_2_dates.py`: Process the spatial changes detection between two temporal spectral index bands images.
- `spatial_changes_detection.ipynb`: Use the file `changes_btw_2_dates.py` to display spatial changes. 

**preprocessing**  
scripts for part 3 (Preprocessing) of the thesis report.
   - `s2_preprocessing.py` is used for Sentinel 2 data preprocessing 
   - `s1_preprocessing.py` is used for Sentinel 1 radar data preprocessing 
   - `l_preprocessing.py` is used for Landsat data preprocessing 
   - `modis_preprocessing.py` is used for MODIS data preprocessing 
   - `data_avalaible.ipynb` generates the frequency of observations available for 

**mad_transfo**  
script for Multivariate Alteration Detection (MAD) transformation, used for change detection.

**ndvi_unsup_d_clustering**  
scripts for part 6 (NDVI density clustering) of the thesis report.
   - `kmeans_clustering.py` compute the unsupervised clustering on a Sentinel 2 image and define the land cover type of each cluster, using different spectral indices median value for each cluster. 
   - `ndvi_density_clustering.ipynb` uses the `kmeans_clustering.py` to display and create time series of land cover in wetlands. 
   - `spatial_changes_detection_unsupervised_clustering.ipynb` is a notebook to detect and display spatial changes btw two images using `kmeans_clustering.py` for land cover identification. 
---

## How to Use

### Installation

### Prerequisites
The scripts require the following Python packages:
- **Core libraries**: `earth-engine-api`, `pandas`, `numpy`, `matplotlib`, `plotly`.
- **Optional libraries**: Check the `import` cells in each notebook for specific dependencies.
 
Those scripts can be utilized as a solid fundation for future work on LULC classification in Bredforsen region, wetlands monitoring, and future Svämängar precise detection. The use of those scripts & notebooks can be made generally downloading the preprocessing folder, as well as the wetlands detection folder and the notebook or folder for the method chosen. Some notebooks doesn't require any of them. 
If needed, some additional comments can be added in the different scripts regarding the needs of datascientists of Vattenfall. If you wish so, please contact me at morgane.magnier@vattenfall.com, or later, morgane.magnier@grenoble-inp.pro 