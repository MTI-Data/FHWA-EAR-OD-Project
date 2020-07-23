# FHWA-EAR-OD-Project
This repository contains source codes developed by the Maryland Transportation Institute for a Federal Highway Administration Exploratory Advanced Research Program Project titled "Data Analytics and Modeling Methods for Tracking and Predicting Origin-Destination Travel Trends based on Mobile Device Data".


This repository contains the source codes used to produce passenger origin-destination tables from anonymized location-based services data and truck origin-destination tables from GPS data. The repository only includes the algorithms that are developed by the University of Maryland. Other algorithms used in the project were developed by the sub-contractors and are proprietary to the sub-contractors,

The Codes directory structure is as follows:

* Person OD: codes that were developed to produce OD tables for person trips.
  - Data Cleaning: data cleaning R script reads daily raw data and perform the cleaning procedures to raw observations to ensure the quality of data.
  - Mode Imputation: a set of R and Python scripts to impute travel mode for trips identified from mobile data location data.
    1. Attribute Extraction: R script to produce attributes for mode imputation.
    2. Distance_Calculation: Python script to calculate distance attributes for mode imputation.
    3. Prediction: Python script to predict the mode based on a trained model.
  - Sociodemographic Imputation: R markdown codes for assigning an age, gender, and income category to each device ID, based on imputed home Census Block Group.
  - Trip Identification: R script to identify which device points form a trip together.

* Truck OD: codes that were developed to produce OD tables for truck trips.
  - Trip_Chaining: Python script to chain truck trips based on truck parking lots and gas stations.

