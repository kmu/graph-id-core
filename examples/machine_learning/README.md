# Sample code for Machine Learning

This directory contains the code used to compare the prediction accuracy of a machine learning model on **leaked** (non-novel) and **unleaked** (novel) test data in the context of materials property prediction.  

## Overview
Two types of test splits are evaluated:

- **Leaked test data**: test structures whose graph IDs also appear in the training set (i.e., structurally non-novel).
- **Unleaked test data**: test structures with graph IDs not seen in the training set (structurally novel).

This comparison highlights the impact of data leakage on performance evaluation in materials informatics.


## Dataset
Materials Project (2019.04.01 release)  
  Download: [Graphs of Materials Project (20190401)](https://figshare.com/articles/dataset/Graphs_of_Materials_Project_20190401/8097992)

- The JSON file (`mp.2019.04.01.json`) should be placed in the working directory.  
  The script will automatically generate and cache a processed pickle file (`mp.2019.04.01.pickle`) with featurized data and graph IDs.


