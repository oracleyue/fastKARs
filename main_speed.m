clear all; close all
rng(4)

% model specification
p = 5;                % AR(p)
T = 250;              % time length
nGroup = 160;         % number of clusters
nLimList = 10:10:200; % list of numbers of signals per group
Var.init = .5;        % variance for initial states
Var.noise = 1;        % variance for noise
tol = 1e-8;           % tolerance for EM
method = 'kARs';
% method = 'mixARs';

switch method
  case 'kARs'
    diary main_speed-kARs.log
  case 'mixARs'
    diary main_speed-mixARs.log
end

% declare spaces
nRuns = length(nLimList);
switch method
  case 'kARs'
    eTimeKARs = zeros(nRuns,1);
    labelMissedKARs = cell(nRuns,1);
    wrongIndexKARs = cell(nRuns,1);
  case 'mixARs'
    eTimeMARs = zeros(nRuns,1);
    labelMissedMARs = cell(nRuns,1);
    wrongIndexMARs = cell(nRuns,1);
end

% benchmarking
for n = 1:nRuns
    nLim = nLimList(n);

    % data generation
    [data, gtLabels, gtModels] = simARs(p, T, nGroup, nLim, Var);

    % time-series clustering
    switch method
      case 'kARs'
        karTimer = tic;
        [estLabels, estModels, estLkd] = kARs(data, nGroup, p, tol);
        eTimeKARs(n) = toc(karTimer);
        fprintf('... %d-th dataset (#signals per group = %d): k-ARs done\n', ...
                n, nLim)
        fprintf('    elapsed time: %f sec.\n', eTimeKARs(n));

      case 'mixARs'
        marTimer = tic;
        [estLabels, estModels, estLkd] = mixARs(data, nGroup, p, tol);
        eTimeMARs(n) = toc(marTimer);
        fprintf('... %d-th dataset (#signals per group = %d): mixARs done\n', ...
                n, nLim)
        fprintf('    elapsed time: %f sec.\n', eTimeMARs(n));
    end

    % performance
    [estLabels, labelMissed, idxPerm] = labelmatch(gtModels, estModels, estLabels);
    idxLabelWrong = find(estLabels-gtLabels);
    switch method
      case 'kARs'
        labelMissedKARs{n} = labelMissed;
        wrongIndexKARs{n} = idxLabelWrong;
      case 'mixARs'
        labelMissedMARs{n} = labelMissed;
        wrongIndexMARs{n} = idxLabelWrong;
    end
end

% save results
garbageVars = {'data','gtLabels','gtModels', ...
               'estLabels', 'estModels', 'estLkd', ...
               'labelMissed', 'idxPerm', 'idxLabelWrong'};
switch method
  case 'kARs'
    clear(garbageVars{:})
    save('benchmark_speed_kARs.mat');
  case 'mixARs'
    clear(garbageVars{:})
    save('benchmark_speed_mixARs2.mat');
end
diary off; system(['mv *.log ./ICASSP/']);