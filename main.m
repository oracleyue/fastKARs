clear all; close all
rng(4)

% model specification
p = 5;            % AR(p)
T = 250;          % time length
K = 6;            % number of clusters
nLim = 40;        % [min,max] or n number of signals per group
N = K*nLim;       % total number of time series
Var.init = .5;    % Var for initial states
Var.noise = 1;    % Var for noise
tol = 1e-8;       % tolerance

% data generation
[data, gtLabels, gtModels] = simARs(p, T, K, nLim, Var);

% time-series clustering
tic
% [estLabels, estModels, estLkd] = mixARs(data, K, p, tol);
[estLabels, estModels, estLkd] = kARs(data, K, p, tol);
toc

% result summary
[updateLabels, idxmap, labelMissed] = minModelGaps(gtModels, estModels, estLabels);
% [updateLabels, idxmap, labelMissed] = minNumMismatch(gtLabels, estLabels);
idxLabelWrong = find(updateLabels-gtLabels);
if isempty(idxLabelWrong)
    disp('All have been successfully clustered.')
else
    fprintf(['Summary of Clustering Results\n' ...
             '-----------------------------\n']);
    fprintf('> %d (%d) time series misclustered.\n', ...
            length(idxLabelWrong), N);
    fprintf('> %d groups missed:\n  ', length(labelMissed));
    fprintf('%d ', labelMissed);
    fprintf('\n\n');
end