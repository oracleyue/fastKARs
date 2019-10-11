function [labels, labelMissed, idxPerm] = labelmatch(gtModels, estModels, estLabels)
% LABELMATCH computes the permutation map that matches the group labels
% of the estimated models to the ground truth labels.
%
% INPUT:
%   gtModels    :   struct; ground truth models
%   estModels   :   struct; estimated models
%
% OUTPUT:
%   labels      :   a vector of {1,...,nGroup}, which has been matched
%                   to share the same group indexing as the ground truth
%   labelMissed :   a vector of integers: missed labels detected during
%                   label matching based model differences
%   idxPerm     :   a permutation of vector 1:nGroup
%                   e.g., idxPerm(i) = j implies the label "i" in estLabels
%                   (estModels, correspondingly) refers to the label "j" in
%                   the ground truth.
%
% Examples:
%   [estLabels, idxPerm] = matching(gtModels, estModels, estLabels)
%   % Now the integers in estLabels and gtLabels refer to the same groups.
%
% Notes:
%   It is possible that "idxPerm" has duplicate integers, since we might
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
%
% Last update on 09 Oct 2019


nGroup = size(gtModels.A, 2);
modelGaps = - ones(nGroup, nGroup);
for m = 1:nGroup
    target = estModels.A(:, m);
    for k = 1:nGroup
        modelGaps(k, m) = norm(gtModels.A(:,k) - target, 2);
    end
end
[~, idxPerm] = min(modelGaps, [], 1);

% check if idxPerm is a permutation map
idxUniq = unique(idxPerm);
if length(idxUniq) < length(idxPerm)
    warning(['The matching failed due to certain ground truth models are too' ...
             'close, or estimation errors are too large.']);
end

% re-labeling according to the matching map
labels = estLabels;
for k = 1:nGroup
    idx = find(estLabels == k);
    labels(idx) = idxPerm(k)/10;
end
labels = labels*10;

% missed labels
labelMissed = [];
for k = 1:nGroup
    if isempty(find(idxUniq == k))
        labelMissed = [labelMissed k];
    end
end