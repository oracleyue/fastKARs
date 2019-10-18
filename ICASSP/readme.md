# Reproducing Numerical Examples in the ICASSP Paper

## Dependencies

These benchmarks depend on the following function in MATLAB: 

- `simARs.m`: simulates mixture AR models to generate benchmark data;

- `kARs.m`: applies k-ARs method to perform time-series clustering;

- `mixARs.m`: applies mixARs method to perform time-series clustering;

- `minModelGaps.m` and `minNumMismatch.m`: matches the labels of ground
  truth clusters and the labels of clustering results.

## Benchmark results

One needs to run the following scripts in the project root to generate
benchmark results:

- `main_speed.m`: benchmark computational time costs over data with
  different number of signals per clusters;
  
- `main_precision.m`: benchmark clustering accuracy over data with different
  number of clusters.
  
## Figures

The following scripts to generate figures in the ICASSP paper:

- `plot_speed.m`: script to generate the figure on time cost benchmark;

- `plot_precision.m`: script the generate the figure on clustering accuracy.