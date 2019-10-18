% This scripts perform clustering for multiple datasets and save the
% clustering results.
% Last modified: 18 Oct 2019


clear all; close all
diary main_precision.log
rng(4)

% model specification
p = 5;               % AR(p)
T = 250;             % time length
nGroupList = 2:2:50; % list of numbers of clusters
nLim = 40;           % number of signals per clusters
Var.init = .5;       % variance for initial states
Var.noise = 1;       % variance for noise
tol = 1e-8;          % tolerance for EM

% declare spaces
nRuns = length(nGroupList);
nSamples = nGroupList * nLim;
% ground truth
gtLabels = cell(nRuns, 1);
gtModels = cell(nRuns, 1);
% results: column 1: kARs; column 2: mixARs
estModels = cell(nRuns, 2);
estLabels = cell(nRuns, 2);

% benchmarking
for n = 1:nRuns
    nGroup = nGroupList(n);

    % data generation
    [data, gtLabel, gtModel] = simARs(p, T, nGroup, nLim, Var);
    fprintf('... %d-th dataset generated\n', n)

    % time-series clustering
    karTimer = tic;
    [estLabelsKARs, estModelsKARs, estLkdKARs] = kARs(data, nGroup, p, tol);
    fprintf('... %d-th dataset: k-ARs done\n', n)
    toc(karTimer)

    marTimer = tic;
    [estLabelsMARs, estModelsMARs, estLkdMARs] = mixARs(data, nGroup, p, tol);
    fprintf('... %d-th dataset: mixARs done\n', n)
    toc(marTimer)

    % save results
    gtModels{n} = gtModel;
    gtLabels{n} = gtLabel;
    estModels{n,1} = estModelsKARs; estModels{n,2} = estModelsMARs;
    estLabels{n,1} = estLabelsKARs; estLabels{n,2} = estLabelsMARs;
end

% save results
garbageVars = {'data','gtLabel','gtModel', ...
               'estLabelsKARs', 'estModelsKARs', 'estLkdKARs', ...
               'estLabelsMARs', 'estModelsMARs', 'estLkdMARs'};
clear(garbageVars{:});
save('benchmark_precision.mat');
diary off; system(['mv *.log ./ICASSP/']);