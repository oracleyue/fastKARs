function [labels, idxmap, labelMissed] = minNumMismatch(gtLabels, estLabels)
% MINNUMMISMATCH matches the cluster labels of the clustering results to the
% ground-truth cluster labels, by minimizing the number of misclustered
% time series in each cluster.
%
% INPUT:
%   gtLabels   :   (Nx1) vector of {1,...,K}
%                  ground-truth cluster labels of N samples
%   estLabels  :   (Nx1) vector of {1,...,K}
%                  labels of clustering results of N samples
% where
%   K - number of clusters
%   N - number of time series (samples).
%
% OUTPUT:
%   labels      :   (Nx1) vector of {1,...,K}
%                   labelling the clustering results using the ground-truth
%                   cluster names.
%   idxmap      :   (Kx1) vector; a permutation of {1,...,K}
%                   e.g., idxmap(i) = j implies that label "i" in the
%                   clustering results matches label "j" in the ground truth.
%   labelMissed :   vector of {1,...,K}
%                   clusters not appearing in results after matching.

% Copyright (c) 2019, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
% Last update on 18 Oct 2019


% parse arguments
N = length(gtLabels);
K = max(gtLabels);

% declare variables
labels = zeros(N, 1);
idxmap = zeros(K, 1);

% label matching
for k = 1:K
    kIndexes = find(gtLabels == k);  % indexes of time series
    kLabels = estLabels(kIndexes);   % clustering labels for k-th gt cluster

    % count frequencies
    uniqLabels = sort(unique(kLabels));
    labelFreq = histcounts(uniqLabels, [uniqLabels; uniqLabels(end)*1.1]);

    % choose the label that maximizes the correct clustering
    [freq, idx] = max(labelFreq);
    klabel = uniqLabels(idx);
    idxmap(k) = klabel;
end

% correct labels for clustering results
for n = 1:N
    labels(n) = idxmap(estLabels(n));
end

% check if theres is any cluster label missed
idxuniq = unique(idxmap);
if length(idxuniq) < length(idxmap)
    warning('There exist several clusters missed during matching.')
end
labelMissed = [];
for k = 1:K
    if isempty(find(idxuniq == k))
        labelMissed = [labelMissed k];
    end
end