function [labels, models, ldf] = mixARs(tsData, nGroup, p, tol)
% MIXARS applies the mixAR-EM approach to clustering time-series data.
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
%   labels = mixARs(tsData, nGroup)
%   labels = mixARs(tsData, nGroup, 'fast')

% Copyright (c) 2019, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 03 Oct 2019


% flags
debugFlag = 0;

% Parsing arguments
[N, Tstar] = size(tsData);
T = Tstar - p;
assert(nGroup < N, ...
       'Argument "nGroup" should be set smaller than #samples.')

% Declare storage variables
XCell = cell(1, N);        % Xi
yMat  = zeros(T, N);       % yi, for yi = Xi*as for AR(p)
astarMat = zeros(p, N);    % astar for N signals
errMat = zeros(T, N);      % pred error
aMat  = zeros(p, nGroup);  % a for each group/model
sigVec = zeros(1, nGroup); % sigma^2 for each group/model
piVec = zeros(1, nGroup);  % pi for each group/model

% Intermediate variables
rMat = ones(N, nGroup)*.5;
phiMat = zeros(N, nGroup);
psiMat = zeros(N, nGroup);
eNormMat = zeros(N, nGroup);

% Precompution
for i = 1:N
    [X, y] = arToepl(tsData(i,:), p);
    astar = X\y;
    XCell{i} = X;
    yMat(:,i) = y;
    astarMat(:,i) = astar;
end
[iX, jX] = size(X);

% Initialization
aMat = astarMat(:, randi(N, [1, nGroup]));
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

        gamma = min(sigVec);
        for k = 1:nGroup
            a_k = aMat(:,k);
            e_ik = y_i - X_i * a_k;
            eNormMat(i,k) = norm(e_ik,2)^2;

            % phiMat(i,k) = evalProb(e_ik, sigVec(k), T);
            % scaling formula to avoid underflow
            psiMat(i,k) = evalPsi(e_ik, sigVec(k), T, gamma);
        end
        % rMat(i,:) = (piVec .* phiMat(i,:)) / sum(piVec .* phiMat(i,:));
        % scaling formula to avoid underflow
        rMat(i,:) = weightedAverage(piVec, psiMat(i,:));
    end

    % M-step
    piVec = mean(rMat);
    for k = 1:nGroup
        XX = zeros(jX, jX);
        Xa = zeros(jX, 1);
        re = 0;
        for i = 1:N
            XX = XX + rMat(i,k) * XCell{i}' * XCell{i};
            Xa = Xa + rMat(i,k) * XCell{i}' * XCell{i} * astarMat(:,i);
            re = re + rMat(i,k) * eNormMat(i,k);
        end
        aMat(:,k) = XX\Xa;

        sigVec(k) = re/T/sum(rMat(:,k));
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

function pval = evalProb(e, sigma, T)
% Evaluate the Gaussian-like function for phi_ik.

    pval = exp(-norm(e, 2)^2 /2/sigma) / sigma^(T/2);

end % END of evalProb

function labels = labelGroups(rMat)
% Determine group labels from variable rMat.

    N = size(rMat, 1);
    labels = zeros(N, 1);
    for n = 1:N
        [~, labels(n)] = max(rMat(n,:));
    end

end % END of labelGroups

function psi = evalPsi(e, sigma, T, gamma)
% Calculate the psi for avoid underflow of computing probability.

    psi = norm(e,2)^2 /2/sigma + T/2*log(sigma/gamma);

end % END of evalPsi

function riVec = weightedAverage(piVec, psiVec)
% Calculate the ldf of being k-th group (k = 1,..., nGroup) by
% the nGroup number of weighted average.

    minPsi = min(psiVec);
    expPsi = exp(-(psiVec - minPsi));
    riVec = (piVec .* expPsi) / sum(piVec .* expPsi);

end % END weightedAverage