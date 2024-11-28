%% Script to generate various scatter summary plots
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.stim_mode, summary.short_stim, ...
    summary.long_stim, summary.response_p, summary.depth,...
    summary.fr_r, summary.fr_z, summary.delta_byEye, ...
    summary.gamma_byEye, summary.delta_IRASA, ...
    summary.gamma_IRASA] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end

fbands = {[2 5], [6 10], [12 28], [30 55]};
c_list = {'red', 'blue','magenta', 'cyan'};

% Load the list of final opto cells
load('data\FinalOptoCells.mat');
dStr_mask = (contains(summary.labels, dStr_opto) &  summary.depth < 3.5);
vStr_mask = (contains(summary.labels, vStr_opto) &  summary.depth >= 3.5);

% Another exclusion criterion is if the difference between the 2 trial counts for a given phase binning is greater than 5% of the trials
% max_dif = 0.05; % Subject to change
% dif_mask = summary.trialbin_dif <= max_dif;

% Significance mask
z_thresh = 2;
sig_mask = (summary.fr_z > z_thresh);
%% Figure7C: Create CSV for UPSET Plot
delta_bool = [0;1;0;0;0;1;1;1;0;0;0;1;1;1;0;1];
theta_bool = [0;0;1;0;0;1;0;0;1;1;0;1;1;0;1;1];
beta_bool  = [0;0;0;1;0;0;1;0;1;0;1;1;0;1;1;1];
gamma_bool = [0;0;0;0;1;0;0;1;0;1;1;0;1;1;1;1];
delta_bool = delta_bool == 1;
theta_bool = theta_bool == 1;
beta_bool = beta_bool == 1;
gamma_bool = gamma_bool == 1;

keep = sig_mask;
[dStr_sig, vStr_sig,all_sig] = deal(zeros(size(delta_bool)));

dStr_sig(1) = sum(~keep(dStr_mask,1) & ~keep(dStr_mask,2) & ~keep(dStr_mask,3) & ~keep(dStr_mask,4));
dStr_sig(2) = sum(keep(dStr_mask,1) & ~keep(dStr_mask,2) & ~keep(dStr_mask,3) & ~keep(dStr_mask,4));
dStr_sig(3) = sum(~keep(dStr_mask,1) & keep(dStr_mask,2) & ~keep(dStr_mask,3) & ~keep(dStr_mask,4));
dStr_sig(4) = sum(~keep(dStr_mask,1) & ~keep(dStr_mask,2) & keep(dStr_mask,3) & ~keep(dStr_mask,4));
dStr_sig(5) = sum(~keep(dStr_mask,1) & ~keep(dStr_mask,2) & ~keep(dStr_mask,3) & keep(dStr_mask,4));
dStr_sig(6) = sum(keep(dStr_mask,1) & keep(dStr_mask,2) & ~keep(dStr_mask,3) & ~keep(dStr_mask,4));
dStr_sig(7) = sum(keep(dStr_mask,1) & ~keep(dStr_mask,2) & keep(dStr_mask,3) & ~keep(dStr_mask,4));
dStr_sig(8) = sum(keep(dStr_mask,1) & ~keep(dStr_mask,2) & ~keep(dStr_mask,3) & keep(dStr_mask,4));
dStr_sig(9) = sum(~keep(dStr_mask,1) & keep(dStr_mask,2) & keep(dStr_mask,3) & ~keep(dStr_mask,4));
dStr_sig(10) = sum(~keep(dStr_mask,1) & keep(dStr_mask,2) & ~keep(dStr_mask,3) & keep(dStr_mask,4));
dStr_sig(11) = sum(~keep(dStr_mask,1) & ~keep(dStr_mask,2) & keep(dStr_mask,3) & keep(dStr_mask,4));
dStr_sig(12) = sum(keep(dStr_mask,1) & keep(dStr_mask,2) & keep(dStr_mask,3) & ~keep(dStr_mask,4));
dStr_sig(13) = sum(keep(dStr_mask,1) & keep(dStr_mask,2) & ~keep(dStr_mask,3) & keep(dStr_mask,4));
dStr_sig(14) = sum(keep(dStr_mask,1) & ~keep(dStr_mask,2) & keep(dStr_mask,3) & keep(dStr_mask,4));
dStr_sig(15) = sum(~keep(dStr_mask,1) & keep(dStr_mask,2) & keep(dStr_mask,3) & keep(dStr_mask,4));
dStr_sig(16) = sum(keep(dStr_mask,1) & keep(dStr_mask,2) & keep(dStr_mask,3) & keep(dStr_mask,4));

vStr_sig(1) = sum(~keep(vStr_mask,1) & ~keep(vStr_mask,2) & ~keep(vStr_mask,3) & ~keep(vStr_mask,4));
vStr_sig(2) = sum(keep(vStr_mask,1) & ~keep(vStr_mask,2) & ~keep(vStr_mask,3) & ~keep(vStr_mask,4));
vStr_sig(3) = sum(~keep(vStr_mask,1) & keep(vStr_mask,2) & ~keep(vStr_mask,3) & ~keep(vStr_mask,4));
vStr_sig(4) = sum(~keep(vStr_mask,1) & ~keep(vStr_mask,2) & keep(vStr_mask,3) & ~keep(vStr_mask,4));
vStr_sig(5) = sum(~keep(vStr_mask,1) & ~keep(vStr_mask,2) & ~keep(vStr_mask,3) & keep(vStr_mask,4));
vStr_sig(6) = sum(keep(vStr_mask,1) & keep(vStr_mask,2) & ~keep(vStr_mask,3) & ~keep(vStr_mask,4));
vStr_sig(7) = sum(keep(vStr_mask,1) & ~keep(vStr_mask,2) & keep(vStr_mask,3) & ~keep(vStr_mask,4));
vStr_sig(8) = sum(keep(vStr_mask,1) & ~keep(vStr_mask,2) & ~keep(vStr_mask,3) & keep(vStr_mask,4));
vStr_sig(9) = sum(~keep(vStr_mask,1) & keep(vStr_mask,2) & keep(vStr_mask,3) & ~keep(vStr_mask,4));
vStr_sig(10) = sum(~keep(vStr_mask,1) & keep(vStr_mask,2) & ~keep(vStr_mask,3) & keep(vStr_mask,4));
vStr_sig(11) = sum(~keep(vStr_mask,1) & ~keep(vStr_mask,2) & keep(vStr_mask,3) & keep(vStr_mask,4));
vStr_sig(12) = sum(keep(vStr_mask,1) & keep(vStr_mask,2) & keep(vStr_mask,3) & ~keep(vStr_mask,4));
vStr_sig(13) = sum(keep(vStr_mask,1) & keep(vStr_mask,2) & ~keep(vStr_mask,3) & keep(vStr_mask,4));
vStr_sig(14) = sum(keep(vStr_mask,1) & ~keep(vStr_mask,2) & keep(vStr_mask,3) & keep(vStr_mask,4));
vStr_sig(15) = sum(~keep(vStr_mask,1) & keep(vStr_mask,2) & keep(vStr_mask,3) & keep(vStr_mask,4));
vStr_sig(16) = sum(keep(vStr_mask,1) & keep(vStr_mask,2) & keep(vStr_mask,3) & keep(vStr_mask,4));

all_sig(1) = sum(~keep(dStr_mask | vStr_mask,1) & ~keep(dStr_mask | vStr_mask,2) & ~keep(dStr_mask | vStr_mask,3) & ~keep(dStr_mask | vStr_mask,4));
all_sig(2) = sum(keep(dStr_mask | vStr_mask,1) & ~keep(dStr_mask | vStr_mask,2) & ~keep(dStr_mask | vStr_mask,3) & ~keep(dStr_mask | vStr_mask,4));
all_sig(3) = sum(~keep(dStr_mask | vStr_mask,1) & keep(dStr_mask | vStr_mask,2) & ~keep(dStr_mask | vStr_mask,3) & ~keep(dStr_mask | vStr_mask,4));
all_sig(4) = sum(~keep(dStr_mask | vStr_mask,1) & ~keep(dStr_mask | vStr_mask,2) & keep(dStr_mask | vStr_mask,3) & ~keep(dStr_mask | vStr_mask,4));
all_sig(5) = sum(~keep(dStr_mask | vStr_mask,1) & ~keep(dStr_mask | vStr_mask,2) & ~keep(dStr_mask | vStr_mask,3) & keep(dStr_mask | vStr_mask,4));
all_sig(6) = sum(keep(dStr_mask | vStr_mask,1) & keep(dStr_mask | vStr_mask,2) & ~keep(dStr_mask | vStr_mask,3) & ~keep(dStr_mask | vStr_mask,4));
all_sig(7) = sum(keep(dStr_mask | vStr_mask,1) & ~keep(dStr_mask | vStr_mask,2) & keep(dStr_mask | vStr_mask,3) & ~keep(dStr_mask | vStr_mask,4));
all_sig(8) = sum(keep(dStr_mask | vStr_mask,1) & ~keep(dStr_mask | vStr_mask,2) & ~keep(dStr_mask | vStr_mask,3) & keep(dStr_mask | vStr_mask,4));
all_sig(9) = sum(~keep(dStr_mask | vStr_mask,1) & keep(dStr_mask | vStr_mask,2) & keep(dStr_mask | vStr_mask,3) & ~keep(dStr_mask | vStr_mask,4));
all_sig(10) = sum(~keep(dStr_mask | vStr_mask,1) & keep(dStr_mask | vStr_mask,2) & ~keep(dStr_mask | vStr_mask,3) & keep(dStr_mask | vStr_mask,4));
all_sig(11) = sum(~keep(dStr_mask | vStr_mask,1) & ~keep(dStr_mask | vStr_mask,2) & keep(dStr_mask | vStr_mask,3) & keep(dStr_mask | vStr_mask,4));
all_sig(12) = sum(keep(dStr_mask | vStr_mask,1) & keep(dStr_mask | vStr_mask,2) & keep(dStr_mask | vStr_mask,3) & ~keep(dStr_mask | vStr_mask,4));
all_sig(13) = sum(keep(dStr_mask | vStr_mask,1) & keep(dStr_mask | vStr_mask,2) & ~keep(dStr_mask | vStr_mask,3) & keep(dStr_mask | vStr_mask,4));
all_sig(14) = sum(keep(dStr_mask | vStr_mask,1) & ~keep(dStr_mask | vStr_mask,2) & keep(dStr_mask | vStr_mask,3) & keep(dStr_mask | vStr_mask,4));
all_sig(15) = sum(~keep(dStr_mask | vStr_mask,1) & keep(dStr_mask | vStr_mask,2) & keep(dStr_mask | vStr_mask,3) & keep(dStr_mask | vStr_mask,4));
all_sig(16) = sum(keep(dStr_mask | vStr_mask,1) & keep(dStr_mask | vStr_mask,2) & keep(dStr_mask | vStr_mask,3) & keep(dStr_mask | vStr_mask,4));

dStr_table = table(delta_bool, theta_bool, beta_bool, gamma_bool, dStr_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz', '12-28 Hz', '30-55 Hz', 'Count'});
vStr_table = table(delta_bool, theta_bool, beta_bool, gamma_bool, vStr_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz', '12-28 Hz', '30-55 Hz', 'Count'});
all_table = table(delta_bool, theta_bool, beta_bool, gamma_bool, all_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz', '12-28 Hz', '30-55 Hz', 'Count'});
writetable(dStr_table, 'C:\Users\mvdmlab\Desktop\dStr_sig.csv');
writetable(vStr_table, 'C:\Users\mvdmlab\Desktop\vStr_sig.csv');
writetable(all_table, 'C:\Users\mvdmlab\Desktop\all_sig.csv');


%% Figure7A and 7B: Plot depth of modulation and Z-Score side-by side
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ax = subplot(2,4,iF);
    hold on
    scatter(summary.depth(dStr_mask), summary.fr_r(dStr_mask,iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
    scatter(summary.depth(dStr_mask & sig_mask(:,iF)), summary.fr_r(dStr_mask & sig_mask(:, iF),iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 50);
    scatter(summary.depth(vStr_mask), summary.fr_r(vStr_mask,iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    scatter(summary.depth(vStr_mask & sig_mask(:,iF)), summary.fr_r(vStr_mask & sig_mask(:,iF),iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'Marker', 'd', 'SizeData', 50);
    [r,p] = corr(summary.fr_r(:,iF), summary.depth);
    legend({sprintf('R^{2} = %.2f, p = %.2f', r,p)})
    ylim([0 0.6])
    xlim([2 5])
    xticks([2 3.5 5])
    yticks([0 0.3 0.6])
    xlabel('Recording depth (mm)')
    ylabel('Modulation strength')
    axis square
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';

    ax = subplot(2,4,iF+4);
    hold on
    scatter(summary.depth(dStr_mask), summary.fr_z(dStr_mask,iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
    scatter(summary.depth(vStr_mask), summary.fr_z(vStr_mask,iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    yline(2, '--black')
    [r,p] = corr(summary.fr_z(:,iF), summary.depth);
    legend({sprintf('R^{2} = %.2f, p = %.2f', r,p), ''})
    ylim([-2 10])
    xlim([2 5])
    xticks([2 3.5 5])
    yticks([-2 2 6 10])
    xlabel('Recording depth (mm)')
    ylabel('z-score')
    axis square;
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
end
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts


%% Diagnostic figure

keep = dStr_mask | vStr_mask;

keep2 = keep & (summary.delta_byEye==1); % High delta sessions checked by eye
keep3 = keep & (summary.delta_byEye==0); % High delta sessions checked by eye
subplot(2,4,1)
title('2 - 5 Hz (by eye)')
hold on
scatter(zeros(1, sum(keep2)), summary.fr_r(keep2,1), 'red', 'filled', 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_r(keep3,1), 'red', 'SizeData', 200);
scatter(0, mean(summary.fr_r(keep2,1)), 'black', 'filled', 'SizeData', 100);
scatter(1, mean(summary.fr_r(keep3,1)), 'black', 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({sprintf('Good delta sessions, n = %d', sum(keep2)), ...
    sprintf('Other sessions, n = %d', sum(keep3))});
ylim([0 0.8])
ylabel('Raw Mod strength')

subplot(2,4,5)
hold on
scatter(zeros(1, sum(keep2)), summary.fr_z(keep2,1), 'red', 'filled', 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_z(keep3,1), 'red', 'SizeData', 200);
scatter(0, mean(summary.fr_z(keep2,1)), 'black', 'filled', 'SizeData', 100);
scatter(1, mean(summary.fr_z(keep3,1)), 'black', 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({sprintf('Good delta sessions, n = %d', sum(keep2)), ...
    sprintf('Other sessions, n = %d', sum(keep3))});
ylim([-3 8])
ylabel('z-scored Mod strength')

keep2 = keep & (summary.gamma_byEye==1); % High gamma sessions checked by eye
keep3 = keep & (summary.gamma_byEye==0); % High gamma sessions checked by eye
subplot(2,4,2)
title('30 - 55 Hz (by eye)')
hold on
scatter(zeros(1, sum(keep2)), summary.fr_r(keep2,4), 'cyan', 'filled', 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_r(keep3,4), 'cyan', 'SizeData', 200);
scatter(0, mean(summary.fr_r(keep2,4)), 'black', 'filled', 'SizeData', 100);
scatter(1, mean(summary.fr_r(keep3,4)), 'black', 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({sprintf('Good gamma sessions, n = %d', sum(keep2)), ...
    sprintf('Other sessions, n = %d', sum(keep3))});
ylim([0 0.8])
ylabel('Raw Mod strength')

subplot(2,4,6)
hold on
scatter(zeros(1, sum(keep2)), summary.fr_z(keep2,4), 'cyan', 'filled', 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_z(keep3,4), 'cyan', 'SizeData', 200);
scatter(0, mean(summary.fr_z(keep2,4)), 'black', 'filled', 'SizeData', 100);
scatter(1, mean(summary.fr_z(keep3,4)), 'black', 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({sprintf('Good gamma sessions, n = %d', sum(keep2)), ...
    sprintf('Other sessions, n = %d', sum(keep3))});
ylim([-3 8])
ylabel('z-scored Mod strength')

keep = dStr_mask | vStr_mask;

keep2 = keep & (summary.delta_IRASA==1); % High delta sessions as per PSD
keep3 = keep & (summary.delta_IRASA==0); % High delta sessions as per PSD
subplot(2,4,3)
title('2 - 5 Hz (by psd)')
hold on
scatter(zeros(1, sum(keep2)), summary.fr_r(keep2,1), 'red', 'filled', 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_r(keep3,1), 'red', 'SizeData', 200);
scatter(0, mean(summary.fr_r(keep2,1)), 'black', 'filled', 'SizeData', 100);
scatter(1, mean(summary.fr_r(keep3,1)), 'black', 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({sprintf('High delta sessions, n = %d', sum(keep2)), ...
    sprintf('Low delta sessions, n = %d', sum(keep3))});
ylim([0 0.8])
ylabel('Raw Mod strength')

subplot(2,4,7)
hold on
scatter(zeros(1, sum(keep2)), summary.fr_z(keep2,1), 'red', 'filled', 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_z(keep3,1), 'red', 'SizeData', 200);
scatter(0, mean(summary.fr_z(keep2,1)), 'black', 'filled', 'SizeData', 100);
scatter(1, mean(summary.fr_z(keep3,1)), 'black', 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({sprintf('High delta sessions, n = %d', sum(keep2)), ...
    sprintf('Low delta sessions, n = %d', sum(keep3))});
ylim([-3 8])
ylabel('z-scored Mod strength')

keep2 = keep & (summary.gamma_IRASA==1); % High gamma sessions as per PSD
keep3 = keep & (summary.gamma_IRASA==0); % High gamma sessions as per PSD
subplot(2,4,4)
title('30 - 55 Hz (by psd)')
hold on
scatter(zeros(1, sum(keep2)), summary.fr_r(keep2,4), 'cyan', 'filled', 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_r(keep3,4), 'cyan', 'SizeData', 200);
scatter(0, mean(summary.fr_r(keep2,4)), 'black', 'filled', 'SizeData', 100);
scatter(1, mean(summary.fr_r(keep3,4)), 'black', 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({sprintf('High gamma sessions, n = %d', sum(keep2)), ...
    sprintf('Low gamma sessions, n = %d', sum(keep3))});
ylim([0 0.8])
ylabel('Raw Mod strength')

subplot(2,4,8)
hold on
scatter(zeros(1, sum(keep2)), summary.fr_z(keep2,4), 'cyan', 'filled', 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_z(keep3,4), 'cyan', 'SizeData', 200);
scatter(0, mean(summary.fr_z(keep2,4)), 'black', 'filled', 'SizeData', 100);
scatter(1, mean(summary.fr_z(keep3,4)), 'black', 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({sprintf('High gamma sessions, n = %d', sum(keep2)), ...
    sprintf('Low Gamma sessions, n = %d', sum(keep3))});
ylim([-3 8])
ylabel('z-scored Mod strength')
%% Make Rebuttal figure (Only show delta by eye)
fig = figure('WindowState', 'maximized');
keep = dStr_mask | vStr_mask;

keep2 = keep & (summary.delta_byEye==1); % High delta sessions checked by eye
keep3 = keep & (summary.delta_byEye==0); % Low delta sessions checked by eye
ax = subplot(2,2,1);
title('2 - 5 Hz', 'FontSize', 35)
hold on
scatter(zeros(1, sum(keep2)), summary.fr_r(keep2,1), 'MarkerFaceColor', 'red', ...
    'MarkerEdgeColor', 'red', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_r(keep3,1), 'MarkerFaceColor', 'red', ...
    'MarkerEdgeColor', 'red', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
scatter(0, mean(summary.fr_r(keep2,1)), 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 100);
scatter(1, mean(summary.fr_r(keep3,1)), 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({'High delta-power sessions', 'Low delta-power sessions'});
% xticklabels({sprintf('High delta-power sessions, n = %d', sum(keep2)), ...
%     sprintf('Low delta-power sessions, n = %d', sum(keep3))});
ylim([0 0.6])
ylabel('Raw Mod strength')
box off;
ax.TickDir = 'out';

ax = subplot(2,2,3);
hold on
scatter(zeros(1, sum(keep2)), summary.fr_z(keep2,1), 'MarkerFaceColor', 'red', ...
    'MarkerEdgeColor', 'red', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
scatter(ones(1, sum(keep3)), summary.fr_z(keep3,1),  'MarkerFaceColor', 'red', ...
    'MarkerEdgeColor', 'red', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
scatter(0, mean(summary.fr_z(keep2,1)), 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 100);
scatter(1, mean(summary.fr_z(keep3,1)), 'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 100);
xlim([-0.25 1.25])
xticks([0 1])
xticklabels({'High delta-power sessions', 'Low delta-power sessions'});
% xticklabels({sprintf('High delta-power sessions, n = %d', sum(keep2)), ...
%     sprintf('Low delta-power sessions, n = %d', sum(keep3))});
ylim([-3 6])
ylabel('z-scored Mod strength')
box off;
ax.TickDir = 'out';
fontname(fig, 'Helvetica')

fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts
exportgraphics(fig,'data\rebuttal.eps','BackgroundColor','none','ContentType','vector')
%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    fbands = {[2 5], [6 10], [12 28], [30 55]};
    delta_eye = load('E:\goodDeltaByEye.mat');
    gamma_eye = load('E:\goodGammaByEye.mat');
    delta_psd = load('E:\goodDeltaIRASA.mat');
    gamma_psd = load('E:\goodDeltaIRASA.mat');
    temp = split(pwd, '\');
    temp = string(temp{end});

    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');

        % Load the phase_response
        load(strcat(fn_prefix, '_phase_response_5_bins.mat')); % Change this to what is decided to be the best binning option
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];
        s_out.response_p = [s_out.response_p; out.overall_response];
        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.stim_mode = [s_out.stim_mode; string(ExpKeys.stim_mode)];
        s_out.fr_z = [s_out.fr_z; out.fr.zscore];
        s_out.fr_r = [s_out.fr_r; out.fr.ratio];

        s_out.delta_byEye = [s_out.delta_byEye; any(strcmp(temp, delta_eye.temp))]; % Good manually picked delta session
        s_out.gamma_byEye = [s_out.gamma_byEye; any(strcmp(temp, gamma_eye.temp))]; % Good manually picked gamma session
        s_out.delta_IRASA = [s_out.delta_IRASA; any(strcmp(temp, delta_psd.temp))]; 
        s_out.gamma_IRASA = [s_out.gamma_IRASA; any(strcmp(temp, gamma_psd.temp))]; 

    end
end

