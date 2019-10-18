% This script loads speed benchmark data and generates the figure used
% in (Yue & Solo, ICASSP 2020).


% data for kARs (nLimList = 10:10:200)
load benchmark_speed_kARs.mat
nLimList_KARs = nLimList;
% eTimeKARs

% data for mixARs (nLimList = 10:10:100)
load benchmark_speed_mixARs.mat  % nLimList = 10:10:50
nLimList_MARs = nLimList;
% Note: if not running separately, skip the following merging operations.
load benchmark_speed_mixARs2.mat % nLimList2 = 60:10:100
nLimList_MARs = [nLimList nLimList2];
eTimeMARs = [eTimeMARs; eTimeMARs2];

% calculate #samples (or #time-series) in dataset
nSampleKARs = nLimList_KARs * nGroup;
nSampleMARs = nLimList_MARs * nGroup;

%% plot and export
fig_hl = figure(1);
% cmap = colormap('lines');
semilogy(nSampleKARs, eTimeKARs, '-o', 'linewidth', .5);
hold on
semilogy(nSampleMARs, eTimeMARs, '-*', 'linewidth', .5);
line_1h  = 3600*1 * ones(length(nLimList_KARs), 1);
line_10h = 3600*10 * ones(length(nLimList_KARs), 1);
semilogy(nSampleKARs, line_10h, 'r-.', 'linewidth', 1);
text(1.8e4, 1e5, '($\geq$ 10 hours)', 'interpreter', 'latex');
semilogy(nSampleKARs, line_1h, 'k-.', 'linewidth', 1);
text(1.8e4, 1e4, '($\geq$ 1 hours)', 'interpreter', 'latex');
hold off
xlim([1.6e3 3.2e4]);
ylim([1 3e5]);
yticks(logspace(0,5,6));
yticklabels({'10^0', '10^1', '10^2', '10^3', '10^4', '10^5'});
xlabel('total number of time series $N$', 'interpreter','latex');
ylabel('CPU time (s)', 'interpreter','latex');
leg_hl = legend('k-ARs', 'mixARs', ...
                'location', 'southeast', 'interpreter','latex');
grid on

% export in PDF
filename = ['benchmark-speed', '.pdf'];
pos = [5.3646 5.2812 4.5833 1.6771];
set(fig_hl,'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, filename, '-dpdf', '-r0');