# PRISM-MS: Mass-Guided Single-Cell MALDI Imaging of Low-Mass Metabolites

DOI: [10.1002/advs.202410506](https://doi.org/10.1002/advs.202410506)

This repository hosts the PRISM-MS application and related resources for the study described in our publication:
**Mass‐Guided Single‐Cell MALDI Imaging of Low‐Mass Metabolites Reveals Cellular Activation Markers.**

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Usage Instructions](#usage-instructions)
- [Supplementary Data](#supplementary-data)

## Overview
PRISM-MS is a Shiny application designed for MALDI imaging analysis workflows. The app simplifies the selection of regions of interest in MALDI datasets and enables guided DeepScan measurements. PRISM-MS is fully self-contained, requiring no additional installation of R or other software.

## Installation
Clone or download this repository.
   ```bash
   git clone https://github.com/CeMOS-Mannheim/PRISM-MS.git
   ```
   
Quick Start (without cloning the repository) :
1. Download `PRISM-MS.rar` [Direct Download newest Release](https://github.com/CeMOS-Mannheim/PRISM-MS/releases/tag/v1.0.0)

   
2.  Unzip the downloaded `PRISM-MS.rar` to a local directory, such as your Desktop.

## Usage Instructions

### Step 1: Run the Application
1. Navigate to the unzipped PRISM-MS folder.
2. Double-click on the `run.bat` file.
   - This will open a command prompt and launch the PRISM-MS Shiny application.
   - No additional installations (e.g., R) are required.

### Step 2: Load Your Data
1. Load an `.imzML` file and the corresponding `.idb` file from your measurements.
   - Sample data is available in the folder: `sample data\231027_GUVs_A`.
2. Load the corresponding `.mis` file.
   - A sample file is also included in the same directory.

### Step 3: Configure Analysis Settings
1. Select a normalization method (usually "none" is sufficient).
2. Set your mass range and tolerance values.
   - Preset values work well for the sample dataset.

### Step 4: Process Data
1. If all data is loaded correctly, a `Load` button will appear. Click it to proceed.
2. Specify your target m/z value (e.g., `786.6` or `678.5` for the GUV dataset).

### Step 5: Adjust Visualization
1. Manually adjust threshold levels to refine the selected mask.
   - Use the Otsu thresholding method until the desired mask is achieved.
2. Enter the acquisition path and adjust spatial resolution and spot dilation as needed.
   - Adjustments will affect the resulting `.mis` file.

### Step 6: Export Data
1. The thresholded image (bottom right of the app) will be converted into a new `.mis` file upon pressing the `Write MIS File` button.
2. Save the `.mis` file to a designated folder for further measurements.

### Step 7: Measurement
1. Use the generated `.mis` file to guide your DeepScan measurement.
2. Example output: `231027_GUVs_A_5_1.mis` (from the sample dataset) demonstrates a DeepScan with 5-micron spatial resolution and a spot dilation factor of 1.

## Supplementary Data
Supplementary datasets used in the manuscript are included in this repository. Notably:
- **Supplementary Dataset S1**: Available as a zip file from [figshare](https://figshare.com/articles/dataset/Supplementary_Dataset_S1_zip/27951516?file=50940864).
- **Measurement Data**: Available as .imzML and .ibd files https://metaspace2020.eu/project/PRISM-MS 
