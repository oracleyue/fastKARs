clear all; close all
rng(4)

% model specification
p = 5;            % AR(p)
T = 250;          % time length
nGroup = 160;     % number of clusters
nLim = 20;       % [min,max] or n number of signals per group
Var.init = .5;    % Var for initial states
Var.noise = 1;    % Var for noise
tol = 1e-8;       % tolerance

% data generation
[data, gtLabels, gtModels] = simARs(p, T, nGroup, nLim, Var);

% time-series clustering
tic
[estLabels, estModels, estLkd] = mixARs(data, nGroup, p, tol);
% [estLabels, estModels, estLkd] = kARs(data, nGroup, p, tol, 'fast');
toc

% result summary
[estLabels, labelMissed, idxPerm] = labelmatch(gtModels, estModels, estLabels);
idxLabelWrong = find(estLabels-gtLabels);
if isempty(idxLabelWrong)
    disp('All have been successfully clustered.')
else
    fprintf(['Summary of Clustering Results\n' ...
             '-----------------------------\n']);
    fprintf('\n> %d time series misclustered:\n  ', length(idxLabelWrong));
    % fprintf('%d ', idxLabelWrong); fprintf('\n');

    fprintf('\n> %d groups missed:\n  ', length(labelMissed));
    fprintf('%d ', labelMissed);
    fprintf('\n');

%     disp('Labels: ground truth vs. clustered')
%     [gtLabels, estLabels]'
end