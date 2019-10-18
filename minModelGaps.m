function [labels, idxmap, labelMissed] = minModelGaps(gtModels, estModels, estLabels)
% MINMODELGAPS matches the cluster labels of the clustering results to the
% ground-truth cluster labels, by minimizing the gaps between estimated
% models and the ground-truth.
%
% INPUT:
%   gtModels    :   struct; ground truth models
%   estModels   :   struct; estimated models
%
% OUTPUT:
%   labels      :   a vector of {1,...,nGroup}, which has been matched
%                   to share the same group indexing as the ground truth
%   idxmap      :   a permutation of vector 1:nGroup
%                   e.g., idxmap(i) = j implies the label "i" in estLabels
%                   (estModels, correspondingly) refers to the label "j" in
%                   the ground truth.
%   labelMissed :   a vector of integers: missed labels detected during
%                   label matching based model differences
%
% Notes:
%   It is possible that "idxmap" has duplicate integers, since we might
%   miss some models due to local optimality. In a result, we force to
%   have "nGroup" models, however, some of them might refers to the same
%   one in the ground truth. That is, we missed certain groups/models,
%   meanwhile, we mistakenly treate the model with minor difference as
%   ones for different groups. It becomes critical when the number of
%   groups increases.

% Copyright (c) 2019, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
% Last update on 18 Oct 2019


nGroup = size(gtModels.A, 2);
modelGaps = - ones(nGroup, nGroup);
for m = 1:nGroup
    target = estModels.A(:, m);
    for k = 1:nGroup
        modelGaps(k, m) = norm(gtModels.A(:,k) - target, 2);
    end
end
[~, idxmap] = min(modelGaps, [], 1);

% check if idxmap is a permutation map
idxuniq = unique(idxmap);
if length(idxuniq) < length(idxmap)
    warning(['There exist several clusters missed, due to some ground-truth ' ...
             'model are too close, or estimation errors are too large']);
end

% re-labeling according to the matching map
labels = estLabels;
for k = 1:nGroup
    idx = find(estLabels == k);
    labels(idx) = idxmap(k)/10;
end
labels = labels*10;

% missed labels
labelMissed = [];
for k = 1:nGroup
    if isempty(find(idxuniq == k))
        labelMissed = [labelMissed k];
    end
end