function [data, labels, models] = simARs(p, T, nGroup, nLim, var)
% SIMARS simulates a prefixed number of AR(p) models for time series
% of a prefixed length.
%
% INPUT:
%   p        :   integer; AR(p)
%   T   	 :   integer; length of scalar time seires
%   nGroup   :   integer; number of clusters, i.e. AR models
%   nLim     :   vector "[nMin, nMax]" OR an integer "n"
%                [nMin, nMax] - the lower and upper bound of #data for
%                each group
%                n - the number of time series per group
%   var      :   struct; variances
%                var.noise - vector/scalar of noise variances for "nGroup"
%                models
%                var.init - scalar/vector of variances of initial states
%
% OUTPUT:
%   data     :   (NxT) matrix; the N number of time series of length T
%   labels   :   (Nx1) vector of group labels {1, ..., nGroup}
%   models   :   (nGroup x 1) cell of struct
%                models{m}.A : AR parameters of model m
%                models{m}.initVar: initial state variance of model m
%                models{m}.noiseVar: noise variance of model m
%
% Examples:
%   var.noise = 1; var.init = .5;
%   data = simARs(4, 100, 3, [10, 20], var);

% Copyright (c) 2019, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 07 Oct 2019


% debug flags
debugPlot = 0;

% Parsing arguments
if length(nLim) == 1
    nSignal = ones(nGroup, 1) * nLim;
elseif length(nLim) == 2
    nSignal = randi(nLim, nGroup, 1);
end
N = sum(nSignal);
assert(T > p, "The argument T must be larger than p for AR(p).");
initVar = var.init;
noiseVar = var.noise;
if isscalar(initVar)
    initVar = ones(1,nGroup) * initVar;
end
if isscalar(noiseVar)
    noiseVar = ones(1,nGroup) * noiseVar;
end

% declare space
data = zeros(N, T);
labels = zeros(N, 1);
models.Pi = nSignal/sum(nSignal);
models.A = zeros(p, nGroup);
models.noiseVar = noiseVar;
models.initVar = initVar;

n = 0;  % indexing samples
for m = 1:nGroup  % indexing models/groups

    A = stabAR(p)';  % coefficents: xt = A(1)*x(t-1) + ... + A(p)*x(t-p)
    c = 0;           % fixed offset, if needed
    models.A(:,m) = A;

    for k = 1:nSignal(m)
        n = n + 1;
        x = zeros(1, T);

        % initial states
        % NOTE: add vairance of initial states
        x(1:p) = randn(1,p) * sqrt(initVar(m));

        % simulation
        for t = p+1:T
            et = randn(1) * sqrt(noiseVar(m));
            x(t) = x(t-p:t-1) * flipud(A) + c + et;
        end

        % recording
        data(n, :) = x;
        labels(n) = m;
    end
end

% permute row orders
rowIdx = randperm(N);
data = data(rowIdx, :);
labels = labels(rowIdx);

% sample plots of time seris
if debugPlot
    figure;
    for k = 1:10
        n = randi(size(data,1));
        plot(data(n,:), '-o');
        hold on
    end
    hold off
    xlabel('time')
    str = sprintf('response of AR(%d) with var(e) = %f', p, noiseVar(1));
    ylabel(str);
end

end % END of simARs

% ================================================================
% Local Functions
% ================================================================

function coeff = stabAR(n)
% STABAR create an AR(n) model, i.e.
%    x(t) = a1*x(t-1) + a2*x(t-2) + ... + an*x(t-n) + e(t),
% where Phi(z) = 1 - a1*z - ... - an*z^n has all poles with norm larger
% than 1. The output is
%    coeff = [a1, a2, ..., an].

    % prob of an integrator is 0.10 for the first and 0.01 for all others
    nint = (rand(1,1)<0.10)+sum(rand(n-1,1)<0.01);

    % prob of repeated roots is 0.05
    nrepeated = floor(sum(rand(n-nint,1)<0.05)/2);

    % prob of complex roots is 0.5
    ncomplex = floor(sum(rand(n-nint-2*nrepeated,1)<0.5)/2);
    nreal = n-nint-2*nrepeated-2*ncomplex;

    % determine random poles
    rep = 2*rand(nrepeated,1)-1;
    mag = rand(ncomplex,1);
    ang = pi*rand(ncomplex,1);
    jay = sqrt(-1);
    complex = mag.*exp(jay*ang);
    re = real(complex); im = imag(complex);
    shrink = .9;
    poles = [re+jay*im; re-jay*im;...
             ones(nint,1)*shrink; ...
             rep;rep;2*rand(nreal,1)*shrink-1];

    % coeff polynomial
    pvec = poly(poles);
    coeff = -pvec(2:end);

end % END of stabAR