%% Script to generate various scatter summary plots
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.bfr] = deal([]);
[summary.labels, summary.stim_mode, summary.short_stim, ...
    summary.long_stim, summary.response_p, summary.depth, ...
    summary.fr_r, summary.fr_z, summary.c_fr_r, summary.c_fr_z, summary.c_fr_z_old, ...
    summary.phaselock_plv, summary.phaselock_mean_phase, ...
    summary.phaselock_pct, summary.phaselock_z, summary.phaselock_max_shufplv, ...
    summary.phaselock_circ_pct, summary.phaselock_circ_z, ...
    summary.phaselock_max_circ_shufplv, summary.excitable_phase, ...
    summary.ntrials, summary.corrected_ex_phase] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end

fbands = {[2 5], [6 10], [30 55]};
c_list = {'red', 'blue', 'green'};

% Load the list of final opto cells and keep the results from only those
load('E:\Dropbox (Dartmouth College)\AnalysisResults\phase_stim_results\FinalOptoCells.mat');
keep = contains(summary.labels, dStr_opto) | contains(summary.labels, vStr_opto);
fn = fieldnames(summary);
for i = 1:numel(fn)
    temp = summary.(fn{i});
    summary.(fn{i}) = temp(keep,:);
end
clear fn temp

dStr_mask = (contains(summary.labels, dStr_opto) &  summary.depth < 3.5);
vStr_mask = (contains(summary.labels, vStr_opto) &  summary.depth >= 3.5);

% Significance mask
z_thresh = 2;
sig_mask = (summary.fr_z > z_thresh);
%% Figure6 (Optional/Extended): Scatter plot of depth vs PLV
pl_thresh = 0.99;
pl_mask = summary.phaselock_pct >= pl_thresh;

fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ax = subplot(1,3,iF);
    hold on
    scatter(summary.depth(dStr_mask), summary.phaselock_plv(dStr_mask,iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
    scatter(summary.depth(dStr_mask & pl_mask(:,iF)), summary.phaselock_plv(dStr_mask & pl_mask(:, iF),iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 50);
    scatter(summary.depth(vStr_mask), summary.phaselock_plv(vStr_mask,iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    scatter(summary.depth(vStr_mask & pl_mask(:,iF)), summary.phaselock_plv(vStr_mask & pl_mask(:,iF),iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'Marker', 'd', 'SizeData', 50);
    ylim([0 0.8])
    xlim([2 5])
    xlabel('Recording depth (mm)')
    ylabel('Phase Locking Value')
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
end
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts

% Get the correlation values between these
[r1,p1] = corr(summary.depth(~isnan(summary.phaselock_plv(:,1)),1), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,1)),1));
[r2,p2] = corr(summary.depth(~isnan(summary.phaselock_plv(:,2))), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,2)),2));
[r3,p3] = corr(summary.depth(~isnan(summary.phaselock_plv(:,3))), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,3)),3));

fprintf("%d - %d Hz, Corr: %.2f, p-value: %.3f\n", fbands{1}(1), fbands{1}(2), r1, p1);
fprintf("%d - %d Hz, Corr: %.2f, p-value: %.3f\n", fbands{2}(1), fbands{2}(2), r2, p2);
fprintf("%d - %d Hz, Corr: %.2f, p-value: %.3f\n", fbands{3}(1), fbands{3}(2), r3, p3);

%% Figure6A: Scatter plot of neurons with significant PLV and mod_strength as well as Z_score (version 1)
pl_thresh = 0.99;
pl_mask = summary.phaselock_pct >= pl_thresh;
z_thresh = 2;
sig_mask = (summary.fr_z > z_thresh);
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ax = subplot(2,3,iF);
    hold on
    keep = dStr_mask & pl_mask(:,iF);
    keep_sig = dStr_mask & pl_mask(:,iF) & sig_mask(:,iF);
    scatter(summary.fr_r(keep,iF), summary.phaselock_plv(keep), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'SizeData', 200);
    scatter(summary.fr_r(keep_sig,iF), summary.phaselock_plv(keep_sig), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 50);
    keep = vStr_mask & pl_mask(:,iF);
    keep_sig = vStr_mask & pl_mask(:,iF) & sig_mask(:,iF);
    scatter(summary.fr_r(keep,iF), summary.phaselock_plv(keep), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'Marker', 'd', 'SizeData', 200);
    scatter(summary.fr_r(keep_sig,iF), summary.phaselock_plv(keep_sig), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'Marker', 'd', 'SizeData', 50);
    xlim([0 0.6])
    ylim([0 0.8])
    xticks([0 0.6])
    yticks([0 0.8])
    ylabel('Phase Locking Value')
    xlabel('Modulation strength')
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';

    ax = subplot(2,3,iF+3);
    hold on
    keep = dStr_mask & pl_mask(:,iF);
    keep_sig = dStr_mask & pl_mask(:,iF) & sig_mask(:,iF);
    scatter(summary.fr_z(keep,iF), summary.phaselock_plv(keep), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'SizeData', 200);
%     scatter(summary.phaselock_plv(keep_sig), summary.fr_z(keep_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
%         'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 50);
    keep = vStr_mask & pl_mask(:,iF);
    keep_sig = vStr_mask & pl_mask(:,iF) & sig_mask(:,iF);
    scatter(summary.fr_z(keep,iF), summary.phaselock_plv(keep), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'Marker', 'd', 'SizeData', 200);
%     scatter(summary.phaselock_plv(keep_sig), summary.fr_z(keep_sig,iF) , 'MarkerFaceColor', c_list{iF}, ...
%         'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'Marker', 'd', 'SizeData', 50);
    xlim([-2.5 8])
    ylim([0 0.8])
    xticks([-2 2 8])
    yticks([0 0.8])
    xline(2, '--black')
    ylabel('Phase Locking Value')
    xlabel('Z-scored modulation strength')
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
end

fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts

% Get the correlation values between these
[r1,p1] = corr(summary.fr_r(~isnan(summary.phaselock_plv(:,1)),1), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,1)),1));
[r2,p2] = corr(summary.fr_r(~isnan(summary.phaselock_plv(:,2))), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,2)),2));
[r3,p3] = corr(summary.fr_r(~isnan(summary.phaselock_plv(:,3))), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,3)),3));
fprintf("%d - %d Hz, Corr b/w mod-strength and PLV: %.2f, p-value: %.3f\n", fbands{1}(1), fbands{1}(2), r1, p1);
fprintf("%d - %d Hz, Corr b/w mod-strength and PLV: %.2f, p-value: %.3f\n", fbands{2}(1), fbands{2}(2), r2, p2);
fprintf("%d - %d Hz, Corr b/w mod-strength and PLV: %.2f, p-value: %.3f\n", fbands{3}(1), fbands{3}(2), r3, p3);

[r1,p1] = corr(summary.fr_z(~isnan(summary.phaselock_plv(:,1)),1), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,1)),1));
[r2,p2] = corr(summary.fr_z(~isnan(summary.phaselock_plv(:,2))), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,2)),2));
[r3,p3] = corr(summary.fr_z(~isnan(summary.phaselock_plv(:,3))), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,3)),3));
fprintf("%d - %d Hz, Corr b/w z-score and PLV: %.2f, p-value: %.3f\n", fbands{1}(1), fbands{1}(2), r1, p1);
fprintf("%d - %d Hz, Corr b/w z-score and PLV: %.2f, p-value: %.3f\n", fbands{2}(1), fbands{2}(2), r2, p2);
fprintf("%d - %d Hz, Corr b/w z-score and PLV: %.2f, p-value: %.3f\n", fbands{3}(1), fbands{3}(2), r3, p3);

%% Figure6A: Scatter plot of neurons with significant PLV and mod_strength (version 2)
pct_thresh = 0.99;
pl_mask = summary.phaselock_pct >= pct_thresh;
z_thresh = 2;
sig_mask = (summary.fr_z > z_thresh);
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ax = subplot(1,3,iF);
    hold on
    keep = dStr_mask & pl_mask(:,iF);
    keep_sig = dStr_mask & pl_mask(:,iF) & sig_mask(:,iF);
    scatter(summary.fr_r(keep,iF), summary.phaselock_plv(keep), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'SizeData', 200);
    scatter(summary.fr_r(keep_sig,iF), summary.phaselock_plv(keep_sig), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0, 'SizeData', 50);
    keep = vStr_mask & pl_mask(:,iF);
    keep_sig = vStr_mask & pl_mask(:,iF) & sig_mask(:,iF);
    scatter(summary.fr_r(keep,iF), summary.phaselock_plv(keep), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'Marker', 'd', 'SizeData', 200);
    scatter(summary.fr_r(keep_sig,iF), summary.phaselock_plv(keep_sig), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0, 'Marker', 'd', 'SizeData', 50);
    xlim([0 0.6])
    ylim([0 0.8])
    xticks([0 0.6])
    yticks([0 0.8])
    ylabel('Phase Locking Value')
    xlabel('Modulation strength')
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
end

fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts

% Get the correlation values between these
[r1,p1] = corr(summary.fr_r(~isnan(summary.phaselock_plv(:,1)),1), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,1)),1));
[r2,p2] = corr(summary.fr_r(~isnan(summary.phaselock_plv(:,2))), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,2)),2));
[r3,p3] = corr(summary.fr_r(~isnan(summary.phaselock_plv(:,3))), ...
    summary.phaselock_plv(~isnan(summary.phaselock_plv(:,3)),3));
fprintf("%d - %d Hz, Corr b/w mod-strength and PLV: %.2f, p-value: %.3f\n", fbands{1}(1), fbands{1}(2), r1, p1);
fprintf("%d - %d Hz, Corr b/w mod-strength and PLV: %.2f, p-value: %.3f\n", fbands{2}(1), fbands{2}(2), r2, p2);
fprintf("%d - %d Hz, Corr b/w mod-strength and PLV: %.2f, p-value: %.3f\n", fbands{3}(1), fbands{3}(2), r3, p3);

%% Figure6A: Scatter plot of neurons with significant PLV and mod_strength (version 3)
pct_thresh = 0.99;
pl_mask = summary.phaselock_pct >= pct_thresh;
z_thresh = 2;
sig_mask = (summary.fr_z > z_thresh);
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ax = subplot(1,3,iF);
    hold on
    keep = pl_mask(:,iF) &  sig_mask(:,iF);
    scatter(summary.fr_r(keep,iF), summary.phaselock_plv(keep), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'SizeData', 800);
    plot([-1 1], [-1 1], '--black')
    [r1,p1] = corr(summary.fr_r(keep,iF), summary.phaselock_plv(keep,iF));
    legend({sprintf('R^{2} = %.2f, p = %.3f', r1,p1), ''}, 'Location', 'northwest', 'FontSize', 25)
    xlim([0 0.8])
    ylim([0 0.8])
    xticks([0 0.8])
    yticks([0 0.8])
    ylabel('Phase Locking Value')
    xlabel('Modulation strength')
    axis square;
%     title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
    ax.XAxis.FontSize = 40;
    ax.YAxis.FontSize = 40;
end

fontname(fig, 'Helvetica')
% fontsize(fig, 30, 'points')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts


%% Input to the Venn diagrams to show phase dependent excitability and phase locking
pct_thresh = 0.99;
z_thresh = 2;
sig_mask = summary.fr_z >= z_thresh;
pl_mask = summary.phaselock_pct >= pct_thresh;
clc
[sum(dStr_mask), sum(dStr_mask & sig_mask(:,1)), sum(dStr_mask & pl_mask(:,1)), ...
    sum(dStr_mask & sig_mask(:,1) & pl_mask(:,1))]
[sum(dStr_mask), sum(dStr_mask & sig_mask(:,2)), sum(dStr_mask & pl_mask(:,2)), ...
    sum(dStr_mask & sig_mask(:,2) & pl_mask(:,2))]
[sum(dStr_mask), sum(dStr_mask & sig_mask(:,3)), sum(dStr_mask & pl_mask(:,3)), ...
    sum(dStr_mask & sig_mask(:,3) & pl_mask(:,3))]
[sum(vStr_mask), sum(vStr_mask & sig_mask(:,1)), sum(vStr_mask & pl_mask(:,1)), ...
    sum(vStr_mask & sig_mask(:,1) & pl_mask(:,1))]
[sum(vStr_mask), sum(vStr_mask & sig_mask(:,2)), sum(vStr_mask & pl_mask(:,2)), ...
    sum(vStr_mask & sig_mask(:,2) & pl_mask(:,2))]
[sum(vStr_mask), sum(vStr_mask & sig_mask(:,3)), sum(vStr_mask & pl_mask(:,3)), ...
    sum(vStr_mask & sig_mask(:,3) & pl_mask(:,3))]

%% Figure 6C: Polar Plots (version1: All, but weighted sum and unweighted sum)
phase_bins = [-pi:2*pi/5:pi];
bin_centers = 0.5*(phase_bins(1:5)+phase_bins(2:6));
bin_marks = unique(sort(mod(rad2deg(phase_bins +2 * pi),360)));

z_thresh = 2;
pct_thresh = 0.99;
sig_mask = summary.fr_z >= z_thresh;
pl_mask = summary.phaselock_pct >= pct_thresh;
clear i

fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    % Look at the circular description of most excitable phases
    keep = find(sig_mask(:,iF));
    this_phase = summary.excitable_phase(keep);
    this_theta = bin_centers(this_phase);
    this_rho = summary.fr_r(keep,iF);
    ax = subplot(2,3,iF);
    for iBin = 1:5
        temp = find(this_phase==iBin);
        for iC = 1:length(temp)  
            circ_jit = power(-1, iC)*iC*pi/18; % for better visualization        
            polarplot([this_theta(temp(iC))+circ_jit, this_theta(temp(iC))+circ_jit], ...
                [0 this_rho(temp(iC))], 'red', 'LineWidth', 1);
            hold on;
        end
    end
    
    % Unweighted sum won't have a length
    mean_angle = circmean(this_theta);
    polarplot([mean_angle, mean_angle], [0, 5], 'Color', 'cyan', 'LineStyle', '--', 'LineWidth', 1)    
    % weighted sum will have an angle and a length
    this_vec = [];
    for iC = 1:length(this_rho)
        this_vec(iC) = this_rho(iC) * (cos(this_theta(iC)) + i*sin(this_theta(iC)));
    end
    sum_vec = sum(this_vec);
    w_rho = abs(sum_vec)/length(this_vec);
    w_theta = angle(sum_vec);
    polarplot([w_theta, w_theta], [0, 5], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1)
    polarplot([w_theta, w_theta], [0, w_rho], 'Color', 'black', 'LineWidth', 3)

    rlim([0, 0.6])
    rticks([0 0.3 0.6])
    rticklabels([])
    thetaticks(bin_marks);
    thetaticklabels({})
%     title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 25)

    % Look at the circular distribution of phase-locked neurons
    keep = find(pl_mask(:,iF));
    this_phase = summary.phaselock_mean_phase(keep,iF);
    this_rho = summary.phaselock_plv(keep,iF);
    ax = subplot(2,3,iF+3);
    for iC = 1:length(this_phase)
        polarplot([this_phase(iC) this_phase(iC)], [0 this_rho(iC)], 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        hold on;
    end
     
    % Unweighted sum won't have a length
    mean_angle = circmean(this_phase);
    polarplot([mean_angle, mean_angle], [0, 5], 'Color', 'cyan', 'LineStyle', '--', 'LineWidth', 1)    
    % weighted sum will have an angle and a length
    this_vec = [];
    for iC = 1:length(this_rho)
        this_vec(iC) = this_rho(iC) * (cos(this_phase(iC)) + i*sin(this_phase(iC)));
    end
    sum_vec = sum(this_vec);
    w_rho = abs(sum_vec)/length(this_vec);
    w_theta = angle(sum_vec);
    polarplot([w_theta, w_theta], [0, 5], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1)
    polarplot([w_theta, w_theta], [0, w_rho], 'Color', 'black', 'LineWidth', 3)
    rlim([0, 0.8])
    rticks([0 0.4 0.8])
    rticklabels([])
    thetaticks(bin_marks);
    thetaticklabels({})
%     title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 25)
end
fontname(fig, 'Helvetica')
fontsize(fig, 45, 'points')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts


%% Figure 6C: Polar Plots (version2: All as well as BOTH, but weighted sum and unweighted sum)
phase_bins = [-pi:2*pi/5:pi];
bin_centers = 0.5*(phase_bins(1:5)+phase_bins(2:6));
bin_marks = unique(sort(mod(rad2deg(phase_bins +2 * pi),360)));

z_thresh = 2;
pct_thresh = 0.99;
sig_mask = summary.fr_z >= z_thresh;
pl_mask = summary.phaselock_pct >= pct_thresh;
clear i

fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    % Look at the circular description of most excitable phases
    keep = find(sig_mask(:,iF));
    keep2 = find(sig_mask(:,iF) & pl_mask(:,iF));
    keep3 = [];
    for iC = 1:length(keep)
        keep3(iC) = sum(keep(iC) == keep2);
    end
    this_phase = summary.excitable_phase(keep);
    this_theta = bin_centers(this_phase);
    this_rho = summary.fr_r(keep,iF);
    ax = subplot(2,3,iF);
    for iBin = 1:5
        temp = find(this_phase==iBin);
        temp2 = find(this_phase(keep3==1));
        temp3 = [];
        for iC = 1:length(temp)
            temp3(iC) = sum(temp(iC) == temp2);
        end
        for iC = 1:length(temp)  
            circ_jit = power(-1, iC)*iC*pi/18; % for better visualization        
            if temp3(iC) == 1
                polarplot([this_theta(temp(iC))+circ_jit, this_theta(temp(iC))+circ_jit], ...
                [0 this_rho(temp(iC))], 'green', 'LineWidth', 1);
                hold on;
            else
                polarplot([this_theta(temp(iC))+circ_jit, this_theta(temp(iC))+circ_jit], ...
                [0 this_rho(temp(iC))], 'red', 'LineWidth', 1);
                hold on;
            end
        end
    end
    
    % Unweighted sum won't have a length
    mean_angle = circmean(this_theta);
    polarplot([mean_angle, mean_angle], [0, 5], 'Color', 'cyan', 'LineStyle', '--', 'LineWidth', 1)
    % For both
    mean_angle = circmean(this_theta(keep3 == 1));
    polarplot([mean_angle, mean_angle], [0, 5], 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 1)
    % weighted sum will have an angle and a length
    this_vec = [];
    for iC = 1:length(this_rho)
        this_vec(iC) = this_rho(iC) * (cos(this_theta(iC)) + i*sin(this_theta(iC)));
    end
    sum_vec = sum(this_vec);
    w_rho = abs(sum_vec)/length(this_vec);
    w_theta = angle(sum_vec);
    polarplot([w_theta, w_theta], [0, 5], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1)
    polarplot([w_theta, w_theta], [0, w_rho], 'Color', 'black', 'LineWidth', 2)
    % For both
    sum_vec = sum(this_vec(keep3 == 1));
    w_rho = abs(sum_vec)/length(this_vec(keep3 ==1));
    w_theta = angle(sum_vec);
    polarplot([w_theta, w_theta], [0, 5], 'Color', 'magenta', 'LineStyle', '--', 'LineWidth', 1)
    polarplot([w_theta, w_theta], [0, w_rho], 'Color', 'magenta', 'LineWidth', 2)
    rlim([0, 0.8])
    rticks([0 0.4 0.8])
    rticklabels([])
    thetaticks(bin_marks);
    thetaticklabels({})
%     title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 25)
       

    % Look at the circular distribution of phase-locked neurons
    keep = find(pl_mask(:,iF));
    this_phase = summary.phaselock_mean_phase(keep,iF);
    this_rho = summary.phaselock_plv(keep,iF);
    ax = subplot(2,3,iF+3);
    for iC = 1:length(this_phase)
        polarplot([this_phase(iC) this_phase(iC)], [0 this_rho(iC)], 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        hold on;
    end
     
    % Unweighted sum won't have a length
    mean_angle = circmean(this_phase);
    polarplot([mean_angle, mean_angle], [0, 5], 'Color', 'cyan', 'LineStyle', '--', 'LineWidth', 1)    
    % weighted sum will have an angle and a length
    this_vec = [];
    for iC = 1:length(this_rho)
        this_vec(iC) = this_rho(iC) * (cos(this_phase(iC)) + i*sin(this_phase(iC)));
    end
    sum_vec = sum(this_vec);
    w_rho = abs(sum_vec)/length(this_vec);
    w_theta = angle(sum_vec);
    polarplot([w_theta, w_theta], [0, 5], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1)
    polarplot([w_theta, w_theta], [0, w_rho], 'Color', 'black', 'LineWidth', 3)
    rticklabels([])
    thetaticks(bin_marks);
    thetaticklabels({})
%     title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 25)

% Look at the circular distribution of phase-locked neurons (BOTH)
    keep = find(pl_mask(:,iF) & sig_mask(:,iF));
    this_phase = summary.phaselock_mean_phase(keep,iF);
    this_rho = summary.phaselock_plv(keep,iF);
    for iC = 1:length(this_phase)
        polarplot([this_phase(iC) this_phase(iC)], [0 this_rho(iC)], 'Color', 'green', 'LineWidth', 1);
        hold on;
    end
     
    % Unweighted sum won't have a length
    mean_angle = circmean(this_phase);
    polarplot([mean_angle, mean_angle], [0, 5], 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 1)    
    % weighted sum will have an angle and a length
    this_vec = [];
    for iC = 1:length(this_rho)
        this_vec(iC) = this_rho(iC) * (cos(this_phase(iC)) + i*sin(this_phase(iC)));
    end
    sum_vec = sum(this_vec);
    w_rho = abs(sum_vec)/length(this_vec);
    w_theta = angle(sum_vec);
    polarplot([w_theta, w_theta], [0, 5], 'Color', 'magenta', 'LineStyle', '--', 'LineWidth', 1)
    polarplot([w_theta, w_theta], [0, w_rho], 'Color', 'magenta', 'LineWidth', 3)
    rlim([0, 0.8])
    rticks([0 0.4 0.8])
    rticklabels([])
    thetaticks(bin_marks);
    thetaticklabels({})
%     title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 25)
end
fontname(fig, 'Helvetica')
fontsize(fig, 45, 'points')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts


%% Figure 6D Scatter between mean_phase and most_excitable phase
[mean_angles, max_ex_angles, foi, locs] = deal([]);
z_thresh = 2;
pct_thresh = 0.99;
sig_mask = summary.fr_z >= z_thresh;
pl_mask = summary.phaselock_pct >= pct_thresh;
phase_bins = [-pi:2*pi/5:pi];
bin_centers = 0.5*(phase_bins(1:5)+phase_bins(2:6));
for iF = 1:3
    keep = sig_mask(:,iF) & pl_mask(:,iF);
    mean_angles = [mean_angles, summary.phaselock_mean_phase(keep,iF)'];
    max_ex_angles = [max_ex_angles, bin_centers(summary.excitable_phase(keep,iF)')];
    foi = [foi, repmat(iF,1, sum(keep))];
    locs = [locs, dStr_mask(keep)'];
end

fig = figure('WindowState', 'maximized');
for iF = 1:3
    idx = find(foi == iF);
    scatter(max_ex_angles(idx), mean_angles(idx), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'SizeData', 200);
    hold on

end
xlim([-3.25 3.25])
ylim([-3.25 3.25])
xticks([-pi,0,pi])
yticks([-pi,0,pi])
xticklabels({'-{\pi}','0','{\pi}'})
yticklabels({'-{\pi}','0','{\pi}'})
ylabel('Mean phase')
xlabel('Most excitable phase')
ax = gca;
ax.TickDir = 'out';
ax.TickLength(1) = 0.03;
ax.Box = 'off';
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts
% [r, p] = circ_corrcc(mean_angles, max_ex_angles);

%% Figure 6D Scatter between mean_phase and most_excitable phase (Version 2 : different bands in different columnss)
z_thresh = 2;
pct_thresh = 0.99;
sig_mask = summary.fr_z >= z_thresh;
pl_mask = summary.phaselock_pct >= pct_thresh;
phase_bins = -pi:2*pi/5:pi;
bin_centers = 0.5*(phase_bins(1:5)+phase_bins(2:6));
fig = figure('WindowState', 'maximized');
for iF = 1:3
    ax = subplot(1,3,iF);
    hold on;
    keep = sig_mask(:,iF) & pl_mask(:,iF);
    mean_angles = summary.phaselock_mean_phase(keep,iF)';
    max_ex_angles = bin_centers(summary.excitable_phase(keep,iF)');
    scatter(max_ex_angles, mean_angles, 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'SizeData', 800);
    plot([-5 5], [-5 5], '--black')
    xlim([-pi pi])
    ylim([-pi pi])
    xticks([-pi,0,pi])
    yticks([-pi,0,pi])
    xticklabels({'-{\pi}','0','{\pi}'})
    yticklabels({'-{\pi}','0','{\pi}'})
    ylabel('Mean phase')
    xlabel('Most excitable phase')
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
    axis square;
    ax.XAxis.FontSize = 45;
    ax.YAxis.FontSize = 45;
end
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts

fig = figure('WindowState', 'maximized');
num_shufs = 1000;
% Plot shuffle histograms for 
for iF = 1:3
    ax = subplot(1,3,iF);
    hold on;
    keep = sig_mask(:,iF) & pl_mask(:,iF);
    mean_angles = summary.phaselock_mean_phase(keep,iF)';
    max_ex_angles = bin_centers(summary.excitable_phase(keep,iF)');
    % absolute
    norm_diff = mod(mean_angles - max_ex_angles, 2*pi);
    min_diff = min(2*pi-norm_diff, norm_diff);
    mean_diff = circ_mean(min_diff');
    sdiff = zeros(1,num_shufs);
    for iShuf = 1:num_shufs
        sidx = randperm(length(max_ex_angles));
        temp_diff = mod(mean_angles - max_ex_angles(sidx), 2*pi);
        temp_diff = min(2*pi-temp_diff, temp_diff);
        sdiff(iShuf) = circ_mean(temp_diff');
    end
    histogram(sdiff, 0:pi/10:pi, 'FaceColor', c_list{iF});
    xline(mean_diff, '--black', 'LineWidth', 3);
%     [this_sum, this_bin] = histcounts(sdiff, 0:pi/10:pi, 'Normalization','cdf');
%     this_bin = 0.5*(this_bin(1:end-1) + this_bin(2:end));
%     this_cdf = [this_bin', this_sum'];
%     [h, p] = kstest(mean_diff, 'CDF', this_cdf);
%     legend({'',sprintf('%.2f',p)});
    legend({'',sprintf('%.2f', 1 - sum(mean_diff<sdiff)/num_shufs)});
    xlim([0 pi])
%     ylim([0 pi])
    xticks([0,pi])
%     yticks([-pi,0,pi])
    xticklabels({'0','{\pi}'})
%     yticklabels({'-{\pi}','0','{\pi}'})
    ylabel('Count')
    xlabel('Mean phase difference')
    t = ax.YTick;
    yticks([t(1),t(end)]);
    ylim([t(1),t(end)]);
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
    ax.XAxis.FontSize = 45;
    ax.YAxis.FontSize = 45;
end
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts
% [r, p] = circ_corrcc(mean_angles, max_ex_angles);
%% Figure 6 Extended 1
%% Figure 6D Scatter between mean_phase and most_excitable phase (Version 2 : different bands in different columnss)
z_thresh = 2;
pct_thresh = 0.99;
sig_mask = summary.fr_z >= z_thresh;
pl_mask = summary.phaselock_pct >= pct_thresh;
phase_bins = -pi:2*pi/5:pi;
bin_centers = 0.5*(phase_bins(1:5)+phase_bins(2:6));
fig = figure('WindowState', 'maximized');
for iF = 1:3
    ax = subplot(1,3,iF);
    hold on;
    keep = sig_mask(:,iF) & pl_mask(:,iF);
    mean_angles = summary.phaselock_mean_phase(keep,iF)';
    max_ex_angles = bin_centers(summary.excitable_phase(keep,iF)');
    scatter(max_ex_angles, mean_angles, 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0, 'SizeData', 800);
    plot([-5 5], [-5 5], '--black')
    xlim([-pi pi])
    ylim([-pi pi])
    xticks([-pi,0,pi])
    yticks([-pi,0,pi])
    xticklabels({'-{\pi}','0','{\pi}'})
    yticklabels({'-{\pi}','0','{\pi}'})
    ylabel('Mean phase')
    xlabel('Most excitable phase')
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
    axis square;
    ax.XAxis.FontSize = 45;
    ax.YAxis.FontSize = 45;
end
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts

%% Statistical test for above
fig = figure('WindowState', 'maximized');
num_shufs = 1000;
% Plot shuffle histograms for 
for iF = 1:3
    ax = subplot(1,3,iF);
    hold on;
    keep = sig_mask(:,iF) & pl_mask(:,iF);
    mean_angles = summary.phaselock_mean_phase(keep,iF)';
    max_ex_angles = bin_centers(summary.excitable_phase(keep,iF)');
    % absolute
    norm_diff = mod(mean_angles - max_ex_angles, 2*pi);
    min_diff = min(2*pi-norm_diff, norm_diff);
    mean_diff = circ_mean(min_diff');
    sdiff = zeros(1,num_shufs);
    for iShuf = 1:num_shufs
        sidx = randperm(length(max_ex_angles));
        temp_diff = mod(mean_angles - max_ex_angles(sidx), 2*pi);
        temp_diff = min(2*pi-temp_diff, temp_diff);
        sdiff(iShuf) = circ_mean(temp_diff');
    end
    histogram(sdiff, 0:pi/10:pi, 'FaceColor', c_list{iF});
    xline(mean_diff, '--black', 'LineWidth', 3);
%     [this_sum, this_bin] = histcounts(sdiff, 0:pi/10:pi, 'Normalization','cdf');
%     this_bin = 0.5*(this_bin(1:end-1) + this_bin(2:end));
%     this_cdf = [this_bin', this_sum'];
%     [h, p] = kstest(mean_diff, 'CDF', this_cdf);
%     legend({'',sprintf('%.2f',p)});
    legend({'',sprintf('%.2f', 1 - sum(mean_diff<sdiff)/num_shufs)});
    xlim([0 pi])
%     ylim([0 pi])
    xticks([0,pi])
%     yticks([-pi,0,pi])
    xticklabels({'0','{\pi}'})
%     yticklabels({'-{\pi}','0','{\pi}'})
    ylabel('Count')
    xlabel('Mean phase difference')
    t = ax.YTick;
    yticks([t(1),t(end)]);
    ylim([t(1),t(end)]);
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
    ax.XAxis.FontSize = 45;
    ax.YAxis.FontSize = 45;
end
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts
%% Figure 6D alternate version (Plot the histogram of circular distance)
fig = figure('WindowState', 'maximized');
phase_difs = [];
z_thresh = 2;
pct_thresh = 0.99;
sig_mask = summary.fr_z >= z_thresh;
pl_mask = summary.phaselock_pct >= pct_thresh;
phase_bins = [-pi:2*pi/5:pi];
bin_centers = 0.5*(phase_bins(1:5)+phase_bins(2:6));
for iF = 1:3
    keep = sig_mask(:,iF) & pl_mask(:,iF);
    mean_angles = summary.phaselock_mean_phase(keep,iF)';
    max_ex_angles = bin_centers(summary.excitable_phase(keep,iF)');
    phase_difs = [phase_difs, abs(circ_dist(mean_angles, max_ex_angles))];
end
histogram(phase_difs, 0:pi/20:pi)
hold on 
xticks([0, pi/3, 2*pi/3, pi])
yticks([0 6])
ylim([0 6])
xticklabels({'0','{\pi}/3', '2{\pi}/3', '{\pi}'})
xline(mean(phase_difs), 'LineWidth', 2, 'LineStyle', '--');
ylabel('Count')
xlabel('Absolute circular distance between mean phase and most-excitable phase')
ax = gca;
ax.TickDir = 'out';
ax.TickLength(1) = 0.03;
ax.Box = 'off';
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts
%% Figure 6D alternate version 2 (Plot the histogram of circular distance)
fig = figure('WindowState', 'maximized');
z_thresh = 2;
pct_thresh = 0.99;
sig_mask = summary.fr_z >= z_thresh;
pl_mask = summary.phaselock_pct >= pct_thresh;
phase_bins = [-pi:2*pi/5:pi];
bin_centers = 0.5*(phase_bins(1:5)+phase_bins(2:6));
for iF = 1:3
    ax = subplot(1,3,iF);
    keep = sig_mask(:,iF) & pl_mask(:,iF);
    mean_angles = summary.phaselock_mean_phase(keep,iF)';
    max_ex_angles = bin_centers(summary.excitable_phase(keep,iF)');
    phase_difs = abs(circ_dist(mean_angles, max_ex_angles));
    histogram(phase_difs, 0:pi/20:pi)
    hold on 
    xticks([0, pi/3, 2*pi/3, pi])
    yticks([0 3])
    ylim([0 3])
    xticklabels({'0','{\pi}/3', '2{\pi}/3', '{\pi}'})
    xline(mean(phase_difs), 'LineWidth', 2, 'LineStyle', '--');
    ylabel('Count')
    xlabel('Absolute circular distance between mean phase and most-excitable phase')
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
end
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts


%% Figure 6-Extended: For neurons that are both phase locked and significantly phase modulated, how do they show up in the non_stim test
fig = figure('WindowState', 'maximized');
for iF = 1:3
%     subplot(2,3,iF)
%     scatter(summary.fr_r(:,iF), summary.c_fr_r(:,iF), 'SizeData', 100, ...
%         'MarkerFaceColor', c_list{iF}, 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0)
%     xlabel('Original Mod. Strength')
%     ylabel('Corrected Mod. Strength')
%     xlim([0 0.8])
%     ylim([0 0.8])
%     title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))

    ax = subplot(1,3,iF);
    scatter(summary.fr_z(:,iF), summary.c_fr_z(:,iF), 'SizeData', 200, ...
        'MarkerFaceColor', c_list{iF}, 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0)
    xlabel('Original Z score')
    ylabel('Corrected Z score')
    hold on
    plot([2,2], [-3,2], '--black')
    plot([2,12], [2,2], '--black')
    plot([-3,2], [2,2], '--black')
    plot([2,2], [2,12], '--black')
    xlim([-2.5 10])
    ylim([-2.5 10])
    xticks([-2 2 10])
    yticks([-2 2 10])
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
end

fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts
%% Version 1
z_thresh= 2;
pct_thresh = 0.99;
sig_mask = summary.fr_z > z_thresh;
c_sig = summary.c_fr_z > z_thresh;
pl_mask = summary.phaselock_pct > pct_thresh;
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ax = subplot(1,3,iF);
    hold on
    plot([1 2], [summary.fr_r(dStr_mask,iF), summary.c_fr_r(dStr_mask,iF)], ...
    c_list{iF}, 'LineStyle', '--')
    plot([1 2], [summary.fr_r(vStr_mask,iF), summary.c_fr_r(vStr_mask,iF)], c_list{iF})
    ylim([0 0.6])
    xlim([0.9 2.1])
    xticks([1 2])
    xticklabels({'Original', 'Corrected'})
    ylabel('Modulation Strength')
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)));
    yticks([0 0.3 0.6])
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
end
%% Version 2

z_thresh= 2;
pct_thresh = 0.99;
sig_mask = summary.fr_z > z_thresh;
c_sig = summary.c_fr_z > z_thresh;
pl_mask = summary.phaselock_pct > pct_thresh;

fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ax = subplot(1,3,iF);
    hold onv 
    plot([1 2], [summary.fr_r(dStr_mask,iF), summary.c_fr_r(dStr_mask,iF)], ...
    c_list{iF}, 'LineStyle', '--')
    plot([1 2], [summary.fr_r(vStr_mask,iF), summary.c_fr_r(vStr_mask,iF)], c_list{iF})
    ylim([0 0.6])
    xlim([0.9 2.1])
    xticks([1 2])
    xticklabels({'Original', 'Corrected'})
    ylabel('Modulation Strength')
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)));
    yticks([0 0.3 0.6])
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
end
%% Diagnostic plot: Look at proportions of neurons
z_thresh = 2;
pct_thresh = 0.99;
sig_ns = summary.ns_fr_z >= z_thresh;
pl_mask = summary.phaselock_pct >= pct_thresh;

fprintf('Dorsal Striatum\n');
for iF = 1:length(fbands)
    keep = find(sig_mask(:,iF) & dStr_mask);
    fprintf('%d - %d Hz\n', fbands{iF}(1),  fbands{iF}(2));
    fprintf('Label\t\t\t\t\tExPhase\t\tNS_Sig\t\tNS_ExPhase\t\tPhaseLocked\n');
    for iC = 1:length(keep)
        fprintf('%s\t %d\t\t\t %d\t\t\t %d\t\t\t\t %d\n', summary.labels(keep(iC)), ...
           summary.excitable_phase(keep(iC),iF), sig_ns(keep(iC),iF), ...
           summary.ns_excitable_phase(keep(iC),iF),pl_mask(keep(iC),iF));
    end
end

fprintf('Ventral Striatum\n')
for iF = 1:length(fbands)
    keep = find(sig_mask(:,iF) & vStr_mask); 
    fprintf('%d - %d Hz\n', fbands{iF}(1),  fbands{iF}(2));
    fprintf('Label\t\t\t\t\tExPhase\t\tNS_Sig\t\tNS_ExPhase\t\tPhaseLocked\n');
    for iC = 1:length(keep)
        fprintf('%s\t %d\t\t\t %d\t\t\t %d\t\t\t\t %d\n', summary.labels(keep(iC)), ...
           summary.excitable_phase(keep(iC),iF), sig_ns(keep(iC),iF), ...
           summary.ns_excitable_phase(keep(iC),iF),pl_mask(keep(iC),iF));
    end
end

%% Diagnostic Plot: Look at proportions of significant phase-locking depending on what thresholding/shuffle method was used 
pl_thresh = 0.99;
z_thresh = 3;
pl_mask = summary.phaselock_pct >= pl_thresh;
pl_zmask = summary.phaselock_z >= z_thresh;
pl_circ_mask = summary.phaselock_circ_pct >= pl_thresh;
pl_circ_zmask = summary.phaselock_circ_z >= z_thresh;
total = ~isnan(summary.phaselock_plv);

fprintf("Fbands:\t\t\t\t\t\t\t\t2 - 5 Hz\t6 - 10 Hz\t 30 - 55Hz\n");

fprintf("dStr Sig Phase Locked:\t\t\t\t %d / %d\t %d / %d\t %d / %d\n", ...
    sum(pl_mask(dStr_mask,1)), sum(total(dStr_mask,1)), sum(pl_mask(dStr_mask,2)), sum(total(dStr_mask,2)), ...
        sum(pl_mask(dStr_mask,3)), sum(total(dStr_mask,3)));
fprintf("vStr Sig Phase Locked:\t\t\t\t %d / %d\t %d / %d\t %d / %d\n", ...
    sum(pl_mask(vStr_mask,1)), sum(total(vStr_mask,1)), sum(pl_mask(vStr_mask,2)), sum(total(vStr_mask,2)), ...
        sum(pl_mask(vStr_mask,3)), sum(total(vStr_mask,3)));

fprintf("dStr Sig Phase Locked(z):\t\t\t %d / %d\t %d / %d\t %d / %d\n", ...
    sum(pl_zmask(dStr_mask,1)), sum(total(dStr_mask,1)), sum(pl_zmask(dStr_mask,2)), sum(total(dStr_mask,2)), ...
        sum(pl_zmask(dStr_mask,3)), sum(total(dStr_mask,3)));
fprintf("vStr Sig Phase Locked(z):\t\t\t %d / %d\t %d / %d\t %d / %d\n", ...
    sum(pl_zmask(vStr_mask,1)), sum(total(vStr_mask,1)), sum(pl_zmask(vStr_mask,2)), sum(total(vStr_mask,2)), ...
        sum(pl_zmask(vStr_mask,3)), sum(total(vStr_mask,3)));

fprintf("dStr Sig Phase Locked(circ):\t\t %d / %d\t\t %d / %d\t\t %d / %d\n", ...
    sum(pl_circ_mask(dStr_mask,1)), sum(total(dStr_mask,1)), sum(pl_circ_mask(dStr_mask,2)), sum(total(dStr_mask,2)), ...
        sum(pl_circ_mask(dStr_mask,3)), sum(total(dStr_mask,3)));
fprintf("vStr Sig Phase Locked(circ):\t\t %d / %d\t\t %d / %d\t\t %d / %d\n", ...
    sum(pl_circ_mask(vStr_mask,1)), sum(total(vStr_mask,1)), sum(pl_circ_mask(vStr_mask,2)), sum(total(vStr_mask,2)), ...
        sum(pl_circ_mask(vStr_mask,3)), sum(total(vStr_mask,3)));

fprintf("dStr Sig Phase Locked(circz):\t\t %d / %d\t\t %d / %d\t\t %d / %d\n", ...
    sum(pl_circ_zmask(dStr_mask,1)), sum(total(dStr_mask,1)), sum(pl_circ_zmask(dStr_mask,2)), sum(total(dStr_mask,2)), ...
        sum(pl_circ_zmask(dStr_mask,3)), sum(total(dStr_mask,3)));
fprintf("vStr Sig Phase Locked(circz):\t\t %d / %d\t\t %d / %d\t\t %d / %d\n", ...
    sum(pl_circ_zmask(vStr_mask,1)), sum(total(vStr_mask,1)), sum(pl_circ_zmask(vStr_mask,2)), sum(total(vStr_mask,2)), ...
        sum(pl_circ_zmask(vStr_mask,3)), sum(total(vStr_mask,3)));
%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
%     fbands = {[2 5], [6 10], [12 30], [30 55]};
    fbands = {[2 5], [6 10], [30 55]};
    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');

        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];
        s_out.stim_mode = [s_out.stim_mode; string(ExpKeys.stim_mode)];
        s_out.short_stim = [s_out.short_stim; ExpKeys.short_stim_pulse_width];
        s_out.long_stim = [s_out.long_stim; ExpKeys.long_stim_pulse_width];

        % Load the stim_responses
        load('stim_phases.mat');
        goodTrials = ExpKeys.goodTrials(iC,:);
        s_out.ntrials = [s_out.ntrials; goodTrials(2) + 1 - goodTrials(1)];  

        % Load the stim-phase responses
        load(strcat(fn_prefix, '_phase_response_5_bins.mat'));
        s_out.response_p = [s_out.response_p; out.overall_response];
        s_out.bfr = [s_out.bfr; {out.bfr}];
        s_out.fr_z = [s_out.fr_z; out.fr.zscore];
        s_out.fr_r = [s_out.fr_r; out.fr.ratio];

        % Load the corrected-phase responses
        load(strcat(fn_prefix, '_corrected_phase_response_5_bins.mat'));       
        s_out.c_fr_z_old = [s_out.c_fr_z_old; corrected_fr_z_oldshufs'];
        s_out.c_fr_z = [s_out.c_fr_z; corrected_fr_z'];
        s_out.c_fr_r = [s_out.c_fr_r; corrected_fr_r'];
        [~, corrected_ex_phase] = max(corrected_fr_bin, [] ,2);
        s_out.corrected_ex_phase = [s_out.corrected_ex_phase; corrected_ex_phase'];

        % Load the phase locking stuff
        load(strcat(fn_prefix, '_spike_phaselock_causal_plv.mat'));
        fn_prefix = strrep(fn_prefix, '_', '-');
        load(strcat(fn_prefix, '_shuf_spec_circ_plv.mat')); % First circularly shifted and then subsampled
%         load(strcat(fn_prefix, '_shuf_spec_circ2_plv.mat')); % Uniform
        load(strcat(fn_prefix, '_shuf_spec_plv.mat'));  % Uniformly distributed fake spikes

        % Get rid of all the 3rd band stuff IF there are 4 bands
        if size(causal_phase, 1) == 4 causal_phase(3,:) = []; end
        if size(shuf_circ_plv, 2) == 4 shuf_circ_plv(:,3) = []; end
        if size(shuf_plv, 2) == 4 shuf_plv(:,3) = []; end     
%         if length(all_spk_phase) == 4 all_spk_phase(3) = []; end
%         if length(all_subsampled_mean_phase) == 4 all_subsampled_mean_phase(3) = []; end
%         if length(all_unsampled_mean_phase) == 4 all_unsampled_mean_phase(3) = []; end
%         if length(all_subsampled_plv) == 4 all_subsampled_plv(3) = []; end
%         if length(all_unsampled_plv) == 4 all_unsampled_plv(3) = []; end
%         if length(trial_spk_phase) == 4 trial_spk_phase(3) = []; end
%         if length(trial_subsampled_mean_phase) == 4 trial_subsampled_mean_phase(3) = []; end
%         if length(trial_unsampled_mean_phase) == 4 trial_unsampled_mean_phase(3) = []; end
%         if length(trial_subsampled_plv) == 4 trial_subsampled_plv(3) = []; end
%         if length(trial_unsampled_plv) == 4 trial_unsampled_plv(3) = []; end
%         if length(pre_spk_phase) == 4 pre_spk_phase(3) = []; end
%         if length(pre_subsampled_mean_phase) == 4 pre_subsampled_mean_phase(3) = []; end
%         if length(pre_unsampled_mean_phase) == 4 pre_unsampled_mean_phase(3) = []; end
%         if length(pre_subsampled_plv) == 4 pre_subsampled_plv(3) = []; end
%         if length(pre_unsampled_plv) == 4 pre_unsampled_plv(3) = []; end
%         if length(post_spk_phase) == 4 post_spk_phase(3) = []; end
%         if length(post_subsampled_mean_phase) == 4 post_subsampled_mean_phase(3) = []; end
%         if length(post_unsampled_mean_phase) == 4 post_unsampled_mean_phase(3) = []; end
%         if length(post_subsampled_plv) == 4 post_subsampled_plv(3) = []; end
%         if length(post_unsampled_plv) == 4 post_unsampled_plv(3) = []; end

        [this_pct, this_circ_pct, this_circ2_pct, this_ex_phase, ...
            this_z, this_circ_z, this_circ2_z, this_ns_ex_phase,] = deal(zeros(1,length(fbands)));
        for iF = 1:length(fbands)
            if isnan(trial_subsampled_plv(iF))
                this_pct(iF) = nan;
                this_z(iF) = nan;
                this_circ_pct(iF) = nan;
                this_circ_z(iF) = nan;
%                 fprintf("There are only %d spikes in %s!\n", ...
%                     size(trial_spk_phase,2), fn_prefix);
%                 this_circ2_pct(iF) = nan;
%                 this_circ2_z(iF) = nan;
            else
                this_pct(iF) = sum(trial_subsampled_plv(iF) > shuf_plv(:,iF))/length(shuf_plv);
                this_z(iF) = (trial_subsampled_plv(iF) - mean(shuf_plv(:,iF)))/std(shuf_plv(:,iF));
                this_circ_pct(iF) = sum(trial_subsampled_plv(iF) > shuf_circ_plv(:,iF))/length(shuf_circ_plv);
                this_circ_z(iF) = (trial_subsampled_plv(iF) - mean(shuf_circ_plv(:,iF)))/std(shuf_circ_plv(:,iF));
%                 this_circ2_pct(iF) = sum(trial_subsampled_plv(iF) > shuf_circ2_plv(:,iF))/length(shuf_circ2_plv);
%                 this_circ2_z(iF) = (trial_subsampled_plv(iF) - mean(shuf_circ2_plv(:,iF)))/std(shuf_circ2_plv(:,iF));
            end
            % The maximally excitable phase
            [~, midx] =  max(out.fr.bin(iF,:));%max(out.fr.bin{iF});
            this_ex_phase(iF) = midx;
        end
      
        if isnan(trial_subsampled_plv(iF)) %this_pct is all nans in this case
            s_out.phaselock_z = [s_out.phaselock_z; this_pct];
            s_out.phaselock_pct = [s_out.phaselock_pct; this_pct];
            s_out.phaselock_max_shufplv = [s_out.phaselock_max_shufplv; this_pct];
            
            s_out.phaselock_circ_z = [s_out.phaselock_circ_z; this_pct];
            s_out.phaselock_circ_pct = [s_out.phaselock_circ_pct; this_pct];
            s_out.phaselock_max_circ_shufplv = [s_out.phaselock_max_circ_shufplv; this_pct];

%             s_out.phaselock_cir2c_z = [s_out.phaselock_circ2_z; this_pct];
%             s_out.phaselock_circ2_pct = [s_out.phaselock_circ2_pct; this_pct];
%             s_out.phaselock_max_circ2_shufplv = [s_out.phaselock_max_circ2_shufplv; this_pct];
            
            s_out.phaselock_plv = [s_out.phaselock_plv; this_pct];
            s_out.phaselock_mean_phase = [s_out.phaselock_mean_phase; this_pct];
        else
            s_out.phaselock_z = [s_out.phaselock_z; this_z];
            s_out.phaselock_pct = [s_out.phaselock_pct; this_pct];
            s_out.phaselock_max_shufplv = [s_out.phaselock_max_shufplv; max(shuf_plv, [], 1)];
            
            s_out.phaselock_circ_z = [s_out.phaselock_circ_z; this_circ_z];
            s_out.phaselock_circ_pct = [s_out.phaselock_circ_pct; this_circ_pct];
            s_out.phaselock_max_circ_shufplv = [s_out.phaselock_max_circ_shufplv; max(shuf_circ_plv, [], 1)];

%             s_out.phaselock_circ2_z = [s_out.phaselock_circ2_z; this_circ2_z];
%             s_out.phaselock_circ2_pct = [s_out.phaselock_circ2_pct; this_circ2_pct];
%             s_out.phaselock_max_circ2_shufplv = [s_out.phaselock_max_circ2_shufplv; max(shuf_circ2_plv, [], 1)];
           
            s_out.phaselock_plv = [s_out.phaselock_plv; trial_subsampled_plv];
            s_out.phaselock_mean_phase = [s_out.phaselock_mean_phase; trial_subsampled_mean_phase];
        end
        s_out.excitable_phase = [s_out.excitable_phase; this_ex_phase];
    end
end

