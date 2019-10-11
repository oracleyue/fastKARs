clear all; close all
rng(4)

% model specification
p = 5;               % AR(p)
T = 250;             % time length
nGroupList = 2:2:50; % list of numbers of clusters
nLim = 40;           % number of signals per group
Var.init = .5;       % variance for initial states
Var.noise = 1;       % variance for noise
tol = 1e-8;          % tolerance for EM

% declare spaces
nRuns = length(nGroupList);
numTSeries = zeros(nRuns, 2);
% column 1: kARs; column 2: mixARs
numWrong = zeros(nRuns, 2);
numMissed = zeros(nRuns, 2);

% benchmarking
for n = 1:nRuns
    nGroup = nGroupList(n);

    % data generation
    [data, gtLabels, gtModels] = simARs(p, T, nGroup, nLim, Var);

    % time-series clustering
    karTimer = tic;
    [estLabelsKARs, estModelsKARs, estLkdKARs] = kARs(data, nGroup, p, tol, 'fast');
    fprintf('... %d-th dataset: k-ARs done', n)
    toc(karTimer)

    marTimer = tic;
    [estLabelsMARs, estModelsMARs, estLkdMARs] = mixARs(data, nGroup, p, tol);
    fprintf('... %d-th dataset: mixARs done', n)
    toc(marTimer)

    % performance
    [estLabels, labelMissed, idxPerm] = labelmatch(gtModels, estModelsKARs, estLabelsKARs);
    idxLabelWrong = find(estLabelsKARs-gtLabels);
    numWrong(n,1) = length(idxLabelWrong);
    numMissed(n,1) = length(labelMissed);
    numTSeries(n) = nLim * nGroup;
end

% save results
garbageVars = {'data','gtLabels','gtModels', ...
               'estLabelsKARs', 'estModelsKARs', 'estLkdKARs', ...
               'estLabelsMARs', 'estModelsMARs', 'estLkdMARs', ...
               'labelMissed', 'idxPerm', 'idxLabelWrong'};
save('benchmark_accuracy.mat');