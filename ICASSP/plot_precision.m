% This script loads data from `main_precision.m` and performs the
% cluster matching to calculate the clustering precision.
% Last modified on 18 Oct 2019

clear all; close all
addpath('../');

% load data
load benchmark_precision.mat

% choose matching methods
type = 'minModelGaps';
% type = 'minNumMismatch';

% declare variables (row 1: k-ARs; row 2: MixARs)
numWrong = zeros(2, nRuns);
numMissed = zeros(2, nRuns);

% match cluster labels
for n = 1:nRuns
    % unwrap data
    K = nGroupList(n);
    N = nSamples(n);
    gtModel = gtModels{n};
    gtLabel = gtLabels{n};
    estModelKARs = estModels{n,1};
    estLabelKARs = estLabels{n,1};
    estModelMARs = estModels{n,2};
    estLabelMARs = estLabels{n,2};

    % matching
    switch type
      case 'minModelGaps'
        % k-ARs
        [updateLabelKARs, idxmapKARs, labelMissedKARs] = minModelGaps(gtModel, estModelKARs, estLabelKARs);
        % mixARs
        [updateLabelMARs, idxmapMARs, labelMissedMARs] = minModelGaps(gtModel, estModelMARs, estLabelMARs);

      case 'minNumMismatch'
        % ...
    end

    % collect statistics
    numWrong(1,n)  = length(find(updateLabelKARs - gtLabel));
    numMissed(1,n) = length(labelMissedKARs);
    numWrong(2,n)  = length(find(updateLabelMARs - gtLabel));
    numMissed(2,n) = length(labelMissedMARs);
end

%% plot and export
fig_hl = figure(2);
pctWrong  = numWrong(1,:) ./ nSamples * 100;
plot(nSamples, 100-pctWrong, '-o');
hold on
pctWrong  = numWrong(2,:) ./ nSamples * 100;
plot(nSamples, 100-pctWrong, '-*');
xlabel('total number of samples $N$', 'interpreter', 'latex')
ylabel('precision (\%)', 'interpreter', 'latex')
legend('k-ARs', 'MixARs', 'interpreter', 'latex')
grid on

filename = ['benchmark-precision', '.pdf'];
pos = [6.8472 6.1111 5.6528 2];
set(fig_hl,'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, filename, '-dpdf', '-r0')


% % plot of clusters missed
% pctMissed = numMissed(1,:) ./ nGroupList * 100;
% plot(nGroupList, pctMissed, '-o');
% hold on
% pctMissed = numMissed(2,:) ./ nGroupList * 100;
% plot(nGroupList, pctMissed, '-o');
% xlabel('number of clusters K')
% ylabel('clusters missed (%)')
% legend('k-ARs', 'MixARs')
