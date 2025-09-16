# Sample Code for Machine Learning

This directory contains code for evaluating the prediction accuracy of a machine learning model on **leaked** (non-novel) vs. **unleaked** (novel) test data in the context of materials property prediction.  

## Overview
Two types of test splits are compared:

- **Leaked test data**: test structures whose graph IDs also appear in the training set (structurally non-novel).  
- **Unleaked test data**: test structures with graph IDs not present in the training set (structurally novel).  

This comparison illustrates how data leakage can distort performance evaluation in materials informatics.

## Dataset
**Materials Project (2019.04.01 release)**  
Download: [Graphs of Materials Project (20190401)](https://figshare.com/articles/dataset/Graphs_of_Materials_Project_20190401/8097992)

- Place the JSON file (`mp.2019.04.01.json`) in the working directory.  
- The script will automatically generate and cache a processed pickle file (`mp.2019.04.01.pickle`) containing featurized data and graph IDs.
