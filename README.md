# ITDM - Modal Identification in the Time Domain with MATLAB

This repository contains a MATLAB implementation of system identification techniques using time-domain methods, particularly focusing on Ibrahim Time Domain (ITD) methods. The code allows for the processing of experimental data to identify dynamic characteristics, such as natural frequencies, damping ratios, and mode shapes of a system.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Data Preparation](#data-preparation)
  - [Running the Analysis](#running-the-analysis)
- [Functions Overview](#functions-overview)
- [Visualization](#visualization)
- [Examples](#examples)
- [License](#license)

## Introduction

This MATLAB project provides tools to analyze time-domain data for system identification purposes. The main focus is on applying the Ibrahim Time Domain (ITD) method, which is widely used for modal analysis in structural dynamics. The repository includes functions for computing auto and cross-correlations, plotting power spectral densities (PSD), and performing stabilization analysis to identify stable poles.

## Features

- **Auto and Cross Correlation:** Calculate the auto-correlation and cross-correlation matrices of the input data.
- **Power Spectral Density (PSD):** Compute and visualize the PSD of the input signal.
- **ITD Method:** Apply the Ibrahim Time Domain method to extract modal parameters.
- **Stability Analysis:** Perform stabilization analysis to identify and filter stable poles.
- **Mode Shape Visualization:** Visualize the identified mode shapes.
- **Graphical User Interface (GUI):** Select relevant parameters for analysis through an intuitive GUI.

## Installation

1. Clone the repository to your local machine:

   ```bash
   git clone https://github.com/shahin1009/ITDM.git
   cd ITDMwithmatlab
   ```

2. Open MATLAB and navigate to the cloned directory.

3. Make sure all dependencies (e.g., Signal Processing Toolbox) are installed in MATLAB.

## Usage

### Data Preparation

Before running the analysis, ensure your data is in the correct format:

- The data should be stored in `.mat` files within the same directory.
- Each `.mat` file should contain the variables `y` (data) and `fsamp` (sampling frequency).

### Running the Analysis

To run the main analysis, execute the following commands in MATLAB:

```matlab
main_ITDM.m
```

This script will perform the following:

1. Load the data from `.mat` files.
2. Calculate auto and cross-correlation matrices.
3. Plot the PSD of the signal.
4. Apply the ITDM method to identify modal parameters.
5. Visualize the stability of the identified poles and mode shapes.

## Functions Overview

- **`Autocorr(Data)`**: Computes the auto-correlation matrix of the input data.
- **`AutocorrVisual(G, lag, limit_x, name)`**: Visualizes the auto-correlation matrix.
- **`PowerSpecPlot(Data, fsamp, splitsec, overlap, name)`**: Plots the power spectral density of the input signal.
- **`ITDM2svd(G, fsamp, order, p, name, Gyy, freq_of_frf)`**: Applies the ITD method incorporating Singular Value Decomposition (SVD) to reduce the system order.
- **`ModeVisual(params_matrix1, 'data1')`**: Visualizes the identified mode shapes.

## Visualization

The repository includes scripts for visualizing various aspects of the analysis:

- **Auto-Correlation**: Plots auto-correlation of the input data.
- **Power Spectral Density (PSD)**: Displays the frequency content of the data.
- **Stability Diagrams**: Shows the stable poles across different model orders.
- **Mode Shapes**: Visualizes the mode shapes of the system.

