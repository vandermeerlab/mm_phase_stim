%% Script to compare the effect of amplitude based thresholding of trials
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.stim_mode, summary.short_stim, ...
    summary.long_stim, summary.depth,...
    summary.fr_r100, summary.fr_z100, ...
    summary.fr_r75, summary.fr_z75, ...
    summary.fr_r50, summary.fr_z50] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end

fbands = {[2 5], [6 10], [12 28] [30 55]};
c_list = {'red', 'blue', 'magenta', 'cyan'};

% Load the list of final opto cells
load('data\FinalOptoCells.mat');
dStr_mask = (contains(summary.labels, dStr_opto) &  summary.depth < 3.5);
vStr_mask = (contains(summary.labels, vStr_opto) &  summary.depth >= 3.5);

%%
% Significance mask
z_thresh = 2;
sig_mask100 = (summary.fr_z100 > z_thresh);
sig_mask75 = (summary.fr_z75 > z_thresh);
sig_mask50 = (summary.fr_z50 > z_thresh);

%%
fig = figure('WindowState', 'maximized');
for iF =  1:length(fbands)
    xcoords = [50, 75, 100];
    ycoords1 = [summary.fr_r50(:,iF), summary.fr_r75(:,iF), summary.fr_r100(:,iF)];
    ycoords2 = [summary.fr_z50(:,iF), summary.fr_z75(:,iF), summary.fr_z100(:,iF)];
    ycoords3 = [sum(sig_mask50(:,iF)), sum(sig_mask75(:,iF)), sum(sig_mask100(:,iF))];
    ax = subplot(3,4,iF);
    hold on
    for iC = 1:size(summary.fr_r50,1)
        plot(xcoords, ycoords1(iC,:), 'Color', c_list{iF});
        scatter(xcoords, ycoords1(iC,:), 'MarkerFaceColor', c_list{iF}, ...
            'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
    end
    plot(xcoords, mean(ycoords1), 'Color', 'black', 'LineWidth', 4);
    xticks(xcoords);
    yticks([0 0.5 1]);
    xlim([45,105]);
    ylim([-0.25 1]);
    xlabel('Percentage of trials', 'FontSize', 14);
    ylabel('FR mod. raw', 'FontSize', 14);
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 16);
    
    ax = subplot(3,4,iF+4);
    hold on
    for iC = 1:size(summary.fr_r50,1)
        plot(xcoords, ycoords2(iC,:), 'Color', c_list{iF});
        scatter(xcoords, ycoords2(iC,:), 'MarkerFaceColor', c_list{iF}, ...
            'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
    end
    plot(xcoords, nanmean(ycoords2), 'Color', 'black', 'LineWidth', 4);
    xticks(xcoords);
    yticks([-2,2,6,10]);
    xlim([45,105]);
    ylim([-3,10]);
    xlabel('Percentage of trials', 'FontSize', 14);
    ylabel('FR mod. (z-scored)', 'FontSize', 14);
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 16);
    dummy =1;

    ax = subplot(3,4,iF+8);
    hold on
    plot(xcoords, ycoords3, 'Color', c_list{iF});
    scatter(xcoords, ycoords3, 'MarkerFaceColor', c_list{iF}, ...
            'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
    xticks(xcoords);
    xlim([45,105]);
    xlabel('Percentage of trials', 'FontSize', 14);
    ylabel('Significant cell count', 'FontSize', 14);
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 18);
end


%% Figure: Create CSV for UPSET Plot
delta_bool = [0;1;0;0;1;1;0;1];
theta_bool = [0;0;1;0;1;0;1;1];
gamma_bool = [0;0;0;1;0;1;1;1];
delta_bool = delta_bool == 1;
theta_bool = theta_bool == 1;
gamma_bool = gamma_bool == 1;

keep = sig_mask;
[dStr_sig, vStr_sig,all_sig] = deal(zeros(size(delta_bool)));


dStr_sig(1) = sum(~keep(dStr_mask,1) & ~keep(dStr_mask,2) & ~keep(dStr_mask,3));
dStr_sig(2) = sum(keep(dStr_mask,1) & ~keep(dStr_mask,2) & ~keep(dStr_mask,3));
dStr_sig(3) = sum(~keep(dStr_mask,1) & keep(dStr_mask,2) & ~keep(dStr_mask,3));
dStr_sig(4) = sum(~keep(dStr_mask,1) & ~keep(dStr_mask,2) & keep(dStr_mask,3));
dStr_sig(5) = sum(keep(dStr_mask,1) & keep(dStr_mask,2) & ~keep(dStr_mask,3));
dStr_sig(6) = sum(keep(dStr_mask,1) & ~keep(dStr_mask,2) & keep(dStr_mask,3));
dStr_sig(7) = sum(~keep(dStr_mask,1) & keep(dStr_mask,2) & keep(dStr_mask,3));
dStr_sig(8) = sum(keep(dStr_mask,1) & keep(dStr_mask,2) & keep(dStr_mask,3));

vStr_sig(1) = sum(~keep(vStr_mask,1) & ~keep(vStr_mask,2) & ~keep(vStr_mask,3));
vStr_sig(2) = sum(keep(vStr_mask,1) & ~keep(vStr_mask,2) & ~keep(vStr_mask,3));
vStr_sig(3) = sum(~keep(vStr_mask,1) & keep(vStr_mask,2) & ~keep(vStr_mask,3));
vStr_sig(4) = sum(~keep(vStr_mask,1) & ~keep(vStr_mask,2) & keep(vStr_mask,3));
vStr_sig(5) = sum(keep(vStr_mask,1) & keep(vStr_mask,2) & ~keep(vStr_mask,3));
vStr_sig(6) = sum(keep(vStr_mask,1) & ~keep(vStr_mask,2) & keep(vStr_mask,3));
vStr_sig(7) = sum(~keep(vStr_mask,1) & keep(vStr_mask,2) & keep(vStr_mask,3));
vStr_sig(8) = sum(keep(vStr_mask,1) & keep(vStr_mask,2) & keep(vStr_mask,3));

all_sig(1) = sum(~keep((dStr_mask | vStr_mask),1) & ~keep((dStr_mask | vStr_mask),2) & ~keep((dStr_mask | vStr_mask),3));
all_sig(2) = sum(keep((dStr_mask | vStr_mask),1) & ~keep((dStr_mask | vStr_mask),2) & ~keep((dStr_mask | vStr_mask),3));
all_sig(3) = sum(~keep((dStr_mask | vStr_mask),1) & keep((dStr_mask | vStr_mask),2) & ~keep((dStr_mask | vStr_mask),3));
all_sig(4) = sum(~keep((dStr_mask | vStr_mask),1) & ~keep((dStr_mask | vStr_mask),2) & keep((dStr_mask | vStr_mask),3));
all_sig(5) = sum(keep((dStr_mask | vStr_mask),1) & keep((dStr_mask | vStr_mask),2) & ~keep((dStr_mask | vStr_mask),3));
all_sig(6) = sum(keep((dStr_mask | vStr_mask),1) & ~keep((dStr_mask | vStr_mask),2) & keep((dStr_mask | vStr_mask),3));
all_sig(7) = sum(~keep((dStr_mask | vStr_mask),1) & keep((dStr_mask | vStr_mask),2) & keep((dStr_mask | vStr_mask),3));
all_sig(8) = sum(keep((dStr_mask | vStr_mask),1) & keep((dStr_mask | vStr_mask),2) & keep((dStr_mask | vStr_mask),3));

dStr_table = table(delta_bool, theta_bool, gamma_bool, dStr_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz','30-55 Hz', 'Count'});
vStr_table = table(delta_bool, theta_bool, gamma_bool, vStr_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz','30-55 Hz', 'Count'});
all_table = table(delta_bool, theta_bool, gamma_bool, all_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz','30-55 Hz', 'Count'});
writetable(dStr_table, 'C:\Users\mvdmlab\Desktop\dStr_sig.csv');
writetable(vStr_table, 'C:\Users\mvdmlab\Desktop\vStr_sig.csv');
writetable(all_table, 'C:\Users\mvdmlab\Desktop\all_sig.csv');


%% Figure3: Plot depth of modulation and Z-Score side-by side
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ax = subplot(2,3,iF);
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

    ax = subplot(2,3,iF+3);
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

%%
% Get the correlation values between these
[r1,p1] = corr(summary.fr_r(:,1), summary.depth);
[r2,p2] = corr(summary.fr_r(:,2), summary.depth);
[r3,p3] = corr(summary.fr_r(:,3), summary.depth);

[r4,p4] = corr(summary.fr_z(:,1), summary.depth);
[r5,p5] = corr(summary.fr_z(:,2), summary.depth);
[r6,p6] = corr(summary.fr_z(:,3), summary.depth);

fprintf("%d - %d Hz, modCorr : %.2f, p-value: %.3f,\t\t z_corr : %.2f, p-value: %.3f" + ...
    "\n", fbands{1}(1), fbands{1}(2), r1, p1, r4, p4);
fprintf("%d - %d Hz, modCorr : %.2f, p-value: %.3f,\t\t z_corr : %.2f, p-value: %.3f" + ...
    "\n", fbands{2}(1), fbands{2}(2), r2, p2, r5, p5);
fprintf("%d - %d Hz, modCorr : %.2f, p-value: %.3f,\t\t z_corr : %.2f, p-value: %.3f" + ...
    "\n", fbands{3}(1), fbands{3}(2), r3, p3, r6, p6);

%% Get the population stats
delta_mod = summary.fr_r(sig_mask(:,1));
delta_z = summary.fr_z(sig_mask(:,1));
theta_mod = summary.fr_r(sig_mask(:,2));
theta_z = summary.fr_z(sig_mask(:,2));
gamma_mod = summary.fr_r(sig_mask(:,3));
gamma_z = summary.fr_z(sig_mask(:,3));

%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');

        % Load the phase_response
        top100 = load(strcat(fn_prefix, '_top_100pct_trials_phase_response_5_bins.mat'));
        top75 = load(strcat(fn_prefix, '_top_75pct_trials_phase_response_5_bins.mat'));
        top50 = load(strcat(fn_prefix, '_top_50pct_trials_phase_response_5_bins.mat'));

        s_out.depth = [s_out.depth; ExpKeys.probeDepth];
        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.stim_mode = [s_out.stim_mode; string(ExpKeys.stim_mode)];
        s_out.fr_z100 = [s_out.fr_z100; top100.out.fr.zscore];
        s_out.fr_r100 = [s_out.fr_r100; top100.out.fr.ratio];
        s_out.fr_z75 = [s_out.fr_z75; top75.out.fr.zscore];
        s_out.fr_r75 = [s_out.fr_r75; top75.out.fr.ratio];
        s_out.fr_z50 = [s_out.fr_z50; top50.out.fr.zscore];
        s_out.fr_r50 = [s_out.fr_r50; top50.out.fr.ratio];

    end
end

