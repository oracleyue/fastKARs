# k-ARs: a fast time-series clustering method for large-scale problems

This project implements the mixture AR method and the k-ARs method to
perform time-series clustering. They both utilize the EM algorithm. The
k-ARs method is proposed in [Yue & Solo, ICASSP, 2019] to improve
computation speed considerably for large-scale problems.

## Functions

The following functions in MATLAB implement the methods proposed in
[ICASSP, 2020]:

- `simARs.m`: simulates mixture AR models to generate benchmark data;

- `kARs.m`: applies k-ARs method to perform time-series clustering;

- `mixARs.m`: applies mixARs method to perform time-series clustering;

- `labelmatch.m`: matches the labels of ground truth clusters and the
  labels of clustering results.

## Demo Scripts

A demo is provided in the script `main.m`.

The following scripts are used to generate benchmark results:

- `main_speed.m`: benchmark computational time costs over data with
  different number of signals per clusters;
  
- `main_accuray.m`: benchmark clustering accuracy over data with different
  number of clusters.
