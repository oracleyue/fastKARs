function [labels, models, ldf] = kARs(tsData, nGroup, p, tol, type)
% KARS applies the k-ARs approach to clustering time-series data, which
% is an extended version of k-means for time series.
 %
% INPUT:
%   tsData   :   (N x T) matrix
%                N - number of samples, each of which is a univariate
%                    time series
%                T - length of time series
%   nGroup     :   positive integer
%                number of groups/clusters
%   p        :   positive integer; the order of AR(p)
%   tol      :   postive float; tolerance, close to 0
%   type     :   char: 'basic', 'fast' (default)
%                'basic' - basic algorithm
%                'fast'  - fast implementation
%
% OUTPUT:
%   labels   :   (N x 1) vector of integers in {1, ..., nGroup}
%                the i-th value tells the group label of the i-th sample
%   models   :   struct
%                models.Pi: (1 x nGroup) matrix of doubles in (0,1)
%                models.A: (p x nGroup) matrix for AR(p) models
%                models.noiseVar: (1 x nGroup) vector for sigma^2 of nGroup models
%   ldf      :   (N x nGroup) matrix of likelihoods being each groups
%                the row indexes refer to samples
%                the column indexes refer to N groups
%
% Examples:
%   labels = kARs(tsData, nGroup)
%   labels = kARs(tsData, nGroup, 'fast')

% Copyright (c) 2019, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 09 Oct 2019


% flags
debugFlag = 0;
warnFlag = 0;

% Parsing arguments
[N, Tstar] = size(tsData);
T = Tstar - p;
assert(nGroup < N, ...
       'Argument "nGroup" should be set smaller than #samples.')
if nargin == 5
    if ~any(strcmpi(type, {'basic', 'fast'}))
        error('Argument "type" is NOT valid.')
    end
else
    type = 'fast';
end

% Declare storage variables
XCell = cell(1, N);        % Xi
yMat  = zeros(T, N);       % yi, for yi = Xi*as for AR(p)
switch type
  case 'basic'
    astarMat = zeros(p, N);    % astar for N signals
  case 'fast'
    yqMat = zeros(p, N);
    RCell = cell(1,N);
end
errMat = zeros(T, N);       % pred error
aMat   = zeros(p, nGroup);  % a for each group/model
sigVec = zeros(1, nGroup);  % sigma^2 for each group/model
piVec  = zeros(1, nGroup);  % pi for each group/model

% Intermediate variables
rMat = ones(N, nGroup)*.5;
eNormMat = zeros(N, nGroup);

% Precompution
idxInit = sort(randi(N, [1, nGroup]));
for i = 1:N
    [X, y] = arToepl(tsData(i,:), p);
    XCell{i} = X;
    switch type
      case 'basic'
        yMat(:,i) = y;
        astarMat(:,i) = X\y;

      case 'fast'
        % use QR decomposition
        [Q, R] = qr(X, 0);  % Q: Txp, R: pxp
        yq = Q'*y;
        yqMat(:,i) = yq;
        RCell{i} = R;
        % initialization
        [~, idxes] = find(idxInit == i);
        if ~isempty(idxes)
            idxStat = nGroup-length(idxInit) + 1;
            idxEnd = idxStat + length(idxes) - 1;
            aMat(:,idxStat:idxEnd) = bsxfun(@times, R\yq, ones(1,length(idxes)));
            idxInit(idxes) = [];
        end
    end
end

% Initialization
if strcmp(type, 'basic')
    aMat = astarMat(:, idxInit);
end
sigVec = ones(1, nGroup);
piVec = ones(1, nGroup)*(1/nGroup);

% debug: initialize convergence plot
if debugFlag
    iter = -100:1:0;
    aList = zeros(size(iter));
    aList(end) = norm(aMat, 'fro');
    sigList = zeros(size(iter));
    sigList(end) = norm(sigVec, 2);
    piList = zeros(size(iter));
    piList(end) = norm(piVec, 2);

    fig_hl = figure(1);
    set(fig_hl, 'units', 'inches', ...
                'position', [5.2500 7.9861 14.1944 4.3750]);
    subplot(1,3,1)
    plH_a = plot(iter, aList, '*-');
    xlabel('Iteration');
    ylabel('Norm of A of groups');
    subplot(1,3,2)
    plH_sig = plot(iter, sigList, '*-');
    xlabel('Iteration');
    ylabel('Norm of sigma^2 of groups');
    subplot(1,3,3)
    plH_pi = plot(iter, piList, '*-');
    xlabel('Iteration');
    ylabel('Norm of Pi of groups');
end

% Expectation Maximization
while 1
    % parameter backup
    aMatPrev = aMat;
    sigVecPrev = sigVec;
    piVecPrev = piVec;

    % E-step
    for i = 1:N
        X_i = XCell{i};
        y_i = yMat(:,i);
        if strcmp(type, 'fast')
            yq_i = yqMat(:,i);
            R_i = RCell{i};
            estarNorm = norm(y_i,2)^2 - norm(yq_i,2)^2;
        end

        for k = 1:nGroup
            a_k = aMat(:,k);
            switch type
              case 'basic'
                e_ik = y_i - X_i * a_k;
                eNormMat(i,k) = norm(e_ik,2)^2;
              case 'fast'
                eNormMat(i,k) = estarNorm + norm(yq_i - R_i*a_k, 2)^2;
            end
        end
        [~, idx] = min(eNormMat(i,:));
        rMat(i,:) = double((1:nGroup) == idx);
    end

    % M-step
    piVec = mean(rMat);
    for k = 1:nGroup
        XX = zeros(p, p);
        Xa = zeros(p, 1);
        idxList = find(rMat(:,k));
        if isempty(idxList)
            msg = sprintf(['No time series contributes to the %d-th ' ...
                           'model.\n Model update has been skipped.'],...
                          k);
            if warnFlag
                warning(msg);
            end
        else
            for i = idxList'
                XX = XX + XCell{i}' * XCell{i};
                switch type
                  case 'basic'
                    Xa = Xa + XCell{i}' * XCell{i} * astarMat(:,i);
                  case 'fast'
                    Xa = Xa + RCell{i}' * yqMat(:,i);
                end
            end
            aMat(:,k) = XX\Xa;
        end
        sigVec(k) = sum(eNormMat(idxList, k))/T/length(idxList);
    end

    % debug: plot convergence curves
    if debugFlag
        iter = iter + 1;
        aList = circshift(aList, -1);
        aList(end) = norm(aMat, 'fro');
        sigList = circshift(sigList, -1);
        sigList(end) = norm(sigVec, 2);
        piList = circshift(piList, -1);
        piList(end) = norm(piVec, 2);
        set(plH_a,   'Xdata', iter, 'Ydata', aList);
        set(plH_sig, 'Xdata', iter, 'Ydata', sigList);
        set(plH_pi,  'Xdata', iter, 'Ydata', piList);
        pause(.1);
    end

    % stopping criteria
    incPi = norm(piVec-piVecPrev, 2);
    incA = norm(aMat-aMatPrev, 'fro');
    incSigma = norm(sigVec-sigVecPrev, 2);
    if max([incPi, incA, incSigma]) < tol
        labels = labelGroups(rMat);
        ldf = rMat;
        models.Pi = piVec;
        models.A = aMat;
        models.noiseVar = sigVec;

        break
    end
end

end % END of kARs

% ================================================================
% Local Functions
% ================================================================

function [X, y] = arToepl(ts, p)
% Build the matrix X for AR(p) estimation min(|y-Xa|).
% Notes: "ts" is a (1xT) row vector.

    T = length(ts);
    rawX = toeplitz(ts);
    X = rawX(p:T-1, 1:p);
    y = ts(p+1:end)';

end % END of arToepl

function labels = labelGroups(rMat)
% Determine group labels from variable rMat.

    N = size(rMat, 1);
    labels = zeros(N, 1);
    for n = 1:N
        [~, labels(n)] = max(rMat(n,:));
    end

end % END of labelGroups
