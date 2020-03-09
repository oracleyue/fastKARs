# k-ARs: a fast method for large-scale time-series clustering

This project implements the mixture AR method and the k-ARs method to
perform time-series clustering. They both utilize the EM algorithm. The
k-ARs method is proposed in [Yue & Solo, 2020] to improve computation speed
considerably for large-scale problems.

*Zuogong Yue and Victor Solo (2020). Large-Scale Time Series Clustering
with k-ARs. In 2020 IEEE International Conference on Acoustics, Speech and
Signal Processing (ICASSP). IEEE, pp. (to be released online)*.

## Functions

The following functions in MATLAB implement the methods proposed in
[ICASSP, 2020]:

- `simARs.m`: simulates mixture AR models to generate benchmark data;

- `kARs.m`: applies k-ARs method to perform time-series clustering;

- `mixARs.m`: applies mixARs method to perform time-series clustering;

- `minModelGaps.m` and `minNumMismatch.m`: matches the labels of ground
  truth clusters and the labels of clustering results.

## Demo Scripts

A demo is provided in the script `main.m`.

The following scripts are used to generate benchmark results:

- `main_speed.m`: benchmark computational time costs over data with
  different number of signals per clusters;
  
- `main_precision.m`: benchmark clustering accuracy over data with different
  number of clusters.
