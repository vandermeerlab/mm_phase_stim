%% Script to generate various scatter summary plots
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.bfr] = deal({});
[summary.labels, summary.response_p, summary.depth,...
    summary.fr_r, summary.fr_z,  summary.phaselock_plv, ...
    summary.phaselock_mean_phase, summary.phaselock_pct, ...
    summary.phaselock_max_shufplv, summary.excitable_phase, ...
    summary.ntrials, summary.trialbin_dif, summary.nbins] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end
% fbands = {[2 5], [6 10], [12 30], [30 55]};
fbands = {[2 5], [6 10], [30 55]};
% c_list = {'red', 'blue','magenta', 'green'};
c_list = {'red', 'blue', 'green'};
% Load the list of final opto cells
load('E:\Dropbox (Dartmouth College)\AnalysisResults\phase_stim_results\FinalOptoCells.mat');
dStr_mask = (contains(summary.labels, dStr_opto) &  summary.depth < 3.5);
vStr_mask = (contains(summary.labels, vStr_opto) &  summary.depth >= 3.5);

% Another exclusion criterion is if the difference between the 2 trial counts for a given phase binning is greater than 5% of the trials
max_dif = 0.05; % Subject to change
dif_mask = summary.trialbin_dif <= max_dif;

% Significance mask
z_thresh = 2;
sig_mask = (summary.fr_z > z_thresh);
%%  Create CSV for UPSET Plot
delta_bool = [0;1;0;0;1;1;0;1];
theta_bool = [0;0;1;0;1;0;1;1];
gamma_bool = [0;0;0;1;0;1;1;1];
delta_bool = delta_bool == 1;
theta_bool = theta_bool == 1;
gamma_bool = gamma_bool == 1;

keep = sig_mask & dif_mask;
[dStr_sig, vStr_sig] = deal(zeros(size(delta_bool)));


dStr_sig(1) = sum(~keep(dStr_mask,1) & ~keep(dStr_mask,2) & ~keep(dStr_mask,4));
dStr_sig(2) = sum(keep(dStr_mask,1) & ~keep(dStr_mask,2) & ~keep(dStr_mask,4));
dStr_sig(3) = sum(~keep(dStr_mask,1) & keep(dStr_mask,2) & ~keep(dStr_mask,4));
dStr_sig(4) = sum(~keep(dStr_mask,1) & ~keep(dStr_mask,2) & keep(dStr_mask,4));
dStr_sig(5) = sum(keep(dStr_mask,1) & keep(dStr_mask,2) & ~keep(dStr_mask,4));
dStr_sig(6) = sum(keep(dStr_mask,1) & ~keep(dStr_mask,2) & keep(dStr_mask,4));
dStr_sig(7) = sum(~keep(dStr_mask,1) & keep(dStr_mask,2) & keep(dStr_mask,4));
dStr_sig(8) = sum(keep(dStr_mask,1) & keep(dStr_mask,2) & keep(dStr_mask,4));

vStr_sig(1) = sum(~keep(vStr_mask,1) & ~keep(vStr_mask,2) & ~keep(vStr_mask,4));
vStr_sig(2) = sum(keep(vStr_mask,1) & ~keep(vStr_mask,2) & ~keep(vStr_mask,4));
vStr_sig(3) = sum(~keep(vStr_mask,1) & keep(vStr_mask,2) & ~keep(vStr_mask,4));
vStr_sig(4) = sum(~keep(vStr_mask,1) & ~keep(vStr_mask,2) & keep(vStr_mask,4));
vStr_sig(5) = sum(keep(vStr_mask,1) & keep(vStr_mask,2) & ~keep(vStr_mask,4));
vStr_sig(6) = sum(keep(vStr_mask,1) & ~keep(vStr_mask,2) & keep(vStr_mask,4));
vStr_sig(7) = sum(~keep(vStr_mask,1) & keep(vStr_mask,2) & keep(vStr_mask,4));
vStr_sig(8) = sum(keep(vStr_mask,1) & keep(vStr_mask,2) & keep(vStr_mask,4));

dStr_table = table(delta_bool, theta_bool, gamma_bool, dStr_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz','30-55 Hz', 'Count'});
vStr_table = table(delta_bool, theta_bool, gamma_bool, vStr_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz','30-55 Hz', 'Count'});
writetable(dStr_table, 'C:\Users\mvdmlab\Desktop\dStr_sig.csv');
writetable(vStr_table, 'C:\Users\mvdmlab\Desktop\vStr_sig.csv');

%% Sort by recording depth and fr_modulation_depth
[~, depth_sorted] = sort(summary.depth);

%% Plot depth of modulation vs Z-score
fig = figure('WindowState', 'maximized');
ax1 = subplot(2,1,1);
hold on;
ax2 = subplot(2,1,2);
hold on;
for iF = 1:length(fbands)
    scatter(ax1, summary.fr_r(dStr_mask & dif_mask(:,iF),iF), summary.fr_z(dStr_mask & dif_mask(:,iF),iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    scatter(ax2, summary.fr_r(vStr_mask & dif_mask(:,iF),iF), summary.fr_z(vStr_mask & dif_mask(:,iF),iF) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
%     legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
end
yline(ax1, 2, '--black');
yline(ax2, 2, '--black');
legend(ax1, {sprintf('%d / %d', sum(dStr_mask & dif_mask(:,1) & (summary.fr_z(:,1)> 2)), sum(dStr_mask & dif_mask(:,1))), ...
    sprintf('%d / %d', sum(dStr_mask & dif_mask(:,2) & (summary.fr_z(:,2)> 2)), sum(dStr_mask & dif_mask(:,2))), ...
    sprintf('%d / %d', sum(dStr_mask & dif_mask(:,3) & (summary.fr_z(:,3)> 2)), sum(dStr_mask & dif_mask(:,3)))}, 'Location', 'southeast');
legend(ax2, {sprintf('%d / %d', sum(vStr_mask & dif_mask(:,1) & (summary.fr_z(:,1)> 2)), sum(vStr_mask & dif_mask(:,1))), ...
    sprintf('%d / %d', sum(vStr_mask & dif_mask(:,2) & (summary.fr_z(:,2)> 2)), sum(vStr_mask & dif_mask(:,2))), ...
    sprintf('%d / %d', sum(vStr_mask & dif_mask(:,3) & (summary.fr_z(:,3)> 2)), sum(vStr_mask & dif_mask(:,3)))}, 'Location', 'southeast');
ax1.Title.String = 'dStr';
ax1.XLim = [0 0.5];
ax1.YLim = [-3 12];
ax1.YLabel.String = 'Z-score';
ax1.XLabel.String = 'Depth of Modulation';
ax1.TickDir = 'out';
ax1.YTick = [-2:2:12];
ax2.Title.String = 'vStr';
ax2.XLim = [0 0.5];
ax2.YLim = [-3 12];
ax2.YLabel.String = 'Z-score';
ax2.XLabel.String = 'Depth of Modulation';
ax2.TickDir = 'out';
ax2.YTick = [-2:2:12];

%% Look at proportions of phase-locking and phase-modulation
pl_thresh = 0.95;
pl_mask = summary.phaselock_pct >= pl_thresh;
all_clean_mask = (pl_mask & sig_mask & dif_mask);


%% Scatter plot of depth of modulation vs Zscore

% dStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
   
    scatter(summary.fr_r(dStr_mask,iF), repmat(iF, sum(dStr_mask),1) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on;
    scatter(summary.fr_r(dStr_mask & this_ex_sig,iF), repmat(iF, sum(dStr_mask & this_ex_sig),1) ,  'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 50);
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
end

yticklabels({'Delta', 'Theta', 'Beta', 'Gamma'}) 
xlabel('Depth of Modulation', 'FontSize', 12)
ylim([0 5])
xlim([0 0.3])
ax = gca;
ax.TickDir = 'out';
sgtitle('Dorsal Striatum')
%
% vStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
   
    scatter(summary.fr_r(vStr_mask,iF), repmat(iF, sum(vStr_mask),1) , 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on;
    scatter(summary.fr_r(dStr_mask & this_ex_sig,iF), repmat(iF, sum(dStr_mask & this_ex_sig),1) ,  'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 50);
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
end

yticklabels({'Delta', 'Theta', 'Beta', 'Gamma'}) 
xlabel('Depth of Modulation', 'FontSize', 12)
ylim([0 5])
xlim([0 0.3])
ax = gca;
ax.TickDir = 'out';
sgtitle('Dorsal Striatum')
%% Scatter plot of baseline firing rate vs depth of modulation
% dStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(cellfun(@mean, summary.bfr(dStr_mask)), summary.fr_r(dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on; 
    scatter(cellfun(@mean, summary.bfr(dStr_mask  & this_ex_sig)), summary.fr_r(dStr_mask & this_ex_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 50);
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    xlabel('Baseline Firing Rate', 'FontSize', 12) 
    ylabel('Depth of Modulation', 'FontSize', 12)
    xlim([-5 35])
    ylim([-0.1 0.5])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
    sgtitle('Dorsal Striatum')
end

% vStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(cellfun(@mean, summary.bfr(vStr_mask)), summary.fr_r(vStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on; 
    scatter(cellfun(@mean, summary.bfr(vStr_mask  & this_ex_sig)), summary.fr_r(vStr_mask & this_ex_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 50);
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    xlabel('Baseline Firing Rate', 'FontSize', 12) 
    ylabel('Depth of Modulation', 'FontSize', 12)
    xlim([-5 35])
    ylim([-0.1 0.5])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
    sgtitle('Ventral Striatum')
end

%% Scatter plot of excitable phase vs intrinsic phase
%Plot for Dorsal Striatum
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_pl_sig = summary.phaselock_sig(:,iF) == 1;
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.phaselock_phase(dStr_mask,iF), summary.excitable_phase(dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on; 
    scatter(summary.phaselock_phase(dStr_mask & this_pl_sig,iF), summary.excitable_phase(dStr_mask & this_pl_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 25);
    scatter(summary.phaselock_phase(dStr_mask & this_ex_sig,iF), summary.excitable_phase(dStr_mask & this_ex_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0.5, 'Marker', '+','SizeData', 200);
    plot([-5 5],[-5 5], 'Color', c_list{iF})
    legend({'all ', 'Sig. Phase Locked', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    xlabel('Intrinsic Phase', 'FontSize', 12) 
    ylabel('Most Excitable Phase', 'FontSize', 12)
    xlim([-3.2 3.2])
    ylim([-3.2 3.2])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Dorsal Striatum')

%% Plot for Ventral Striatum
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_pl_sig = summary.phaselock_sig(:,iF) == 1;
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.phaselock_phase(vStr_mask,iF), summary.excitable_phase(vStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
    hold  on; 
    scatter(summary.phaselock_phase(vStr_mask & this_pl_sig,iF), summary.excitable_phase(vStr_mask & this_pl_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SizeData', 25);
    scatter(summary.phaselock_phase(vStr_mask & this_ex_sig,iF), summary.excitable_phase(vStr_mask & this_ex_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0.5, 'Marker', '+','SizeData', 200);
    plot([-5 5],[-5 5], 'Color', c_list{iF})
    legend({'all ', 'Sig. Phase Locked', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    xlabel('Intrinsic Phase', 'FontSize', 12) 
    ylabel('Most Excitable Phase', 'FontSize', 12)
    xlim([-3.2 3.2])
    ylim([-3.2 3.2])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Ventral Striatum')
%% Scatter of most excitable phases  vs depth of modulationfor significant results
% dStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.excitable_phase(dStr_mask,iF), summary.fr_r(dStr_mask,iF),'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on; 
    scatter(summary.excitable_phase(dStr_mask & this_ex_sig,iF), summary.fr_r(dStr_mask & this_ex_sig,iF),  'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 25);
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    ylabel('Depth of modulation', 'FontSize', 12) 
    xlabel('Most Excitable Phase', 'FontSize', 12)
    ylim([-0.1 0.5])
    xlim([-3.2 3.2])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Dorsal Striatum')

%% vStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.excitable_phase(vStr_mask,iF), summary.fr_r(vStr_mask,iF),'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on; 
    scatter(summary.excitable_phase(vStr_mask & this_ex_sig,iF), summary.fr_r(vStr_mask & this_ex_sig,iF),  'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 25);
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    ylabel('Depth of modulation', 'FontSize', 12) 
    xlabel('Most Excitable Phase', 'FontSize', 12)
    ylim([-0.1 0.5])
    xlim([-3.2 3.2])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Ventral Striatum')

%% Look at why high depth of modulation is not necessarily significant

% dStr
fig = figure('WindowState', 'maximized');
for iF = 1%1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2; 
    ax = subplot(2,2,iF);
    sel1 = find(~this_ex_sig & dStr_mask);
    sel2 = find(this_ex_sig & dStr_mask);
    for iS = 1:length(sel1)
%         text(summary.excitable_phase(sel1(iS),iF), summary.fr_r(sel1(iS),iF), ...
%             string(summary.response_p(sel1(iS))));
        text(3*summary.excitable_phase(sel1(iS),iF), 5*summary.fr_r(sel1(iS),iF), ...
            summary.labels(sel1(iS)), 'Interpreter', 'none');
    end
    for iS = 1:length(sel2)
%         text(summary.excitable_phase(sel2(iS),iF), summary.fr_r(sel2(iS),iF), ...
%             string(summary.response_p(sel1(iS))), 'FontWeight', 'bold');
        text(3*summary.excitable_phase(sel2(iS),iF), 5*summary.fr_r(sel2(iS),iF), ...
            summary.labels(sel2(iS)), 'FontWeight', 'bold', 'Interpreter', 'none');
    end
    ylabel('Depth of modulation', 'FontSize', 12) 
    xlabel('Most Excitable Phase', 'FontSize', 12)
    ylim([-0.1 2.5])
    xlim([-9.6 9.6])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Dorsal Striatum')

%% vStr
fig = figure('WindowState', 'maximized');
for iF = 4%1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
    ax = subplot(2,2,iF);
    sel1 = find(~this_ex_sig & vStr_mask);
    sel2 = find(this_ex_sig & vStr_mask);
    for iS = 1:length(sel1)
%         text(summary.excitable_phase(sel1(iS),iF), summary.fr_r(sel1(iS),iF), ...
%             string(summary.response_p(sel1(iS))));
        text(3*summary.excitable_phase(sel1(iS),iF), 5*summary.fr_r(sel1(iS),iF), ...
            summary.labels(sel1(iS)), 'Interpreter', 'none');
    end
    for iS = 1:length(sel2)
%         text(summary.excitable_phase(sel2(iS),iF), summary.fr_r(sel2(iS),iF), ...
%             string(summary.response_p(sel1(iS))), 'FontWeight', 'bold');
        text(3*summary.excitable_phase(sel2(iS),iF), 5*summary.fr_r(sel2(iS),iF), ...
            summary.labels(sel2(iS)), 'FontWeight', 'bold', 'Interpreter', 'none');
    end
    ylabel('Depth of modulation', 'FontSize', 12) 
    xlabel('Most Excitable Phase', 'FontSize', 12)
    ylim([-0.1 2.5])
    xlim([-9.6 9.6])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Ventral Striatum')


%% Plot for Ventral Striatum
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_pl_sig = summary.phaselock_sig(:,iF) == 1;
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.phaselock_phase(vStr_mask,iF), summary.excitable_phase(vStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
    hold  on; 
    scatter(summary.phaselock_phase(vStr_mask & this_pl_sig,iF), summary.excitable_phase(vStr_mask & this_pl_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SizeData', 25);
    scatter(summary.phaselock_phase(vStr_mask & this_ex_sig,iF), summary.excitable_phase(vStr_mask & this_ex_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0.5, 'Marker', '+','SizeData', 200);
    plot([-5 5],[-5 5], 'Color', c_list{iF})
    legend({'all ', 'Sig. Phase Locked', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    xlabel('Intrinsic Phase', 'FontSize', 12) 
    ylabel('Most Excitable Phase', 'FontSize', 12)
    xlim([-3.2 3.2])
    ylim([-3.2 3.2])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Ventral Striatum')

%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    p_thresh = 0.99;
    fbands = {[2 5], [6 10], [12 30], [30 55]};
    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        % Load the phase_response
        load(strcat(fn_prefix, '_phase_response_5_bins.mat')); % Change this to what is decided to be the best binning option
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];
        s_out.response_p = [s_out.response_p; out.overall_response];
        s_out.bfr{length(s_out.bfr)+1} = out.bfr;
        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.fr_z = [s_out.fr_z; out.fr.zscore];
        s_out.fr_r = [s_out.fr_r; out.fr.ratio];

        % Load the stim_responses
        load('stim_phases.mat');
        goodTrials = ExpKeys.goodTrials(iC,:);

        s_out.ntrials = [s_out.ntrials; goodTrials(2) + 1 - goodTrials(1)];
        
        % Load phase_lock and shuf_spec
        fn_prefix = strrep(fn_prefix, '_', '-');
        load(strcat(fn_prefix, '_spike_phaselock_plv.mat'));
        load(strcat(fn_prefix, '_shuf_spec_plv.mat'));  
        [this_pct, this_ex_phase, this_trialbin_dif] = deal(zeros(1,length(fbands)));
        for iF = 1:length(fbands)
            if isempty(trial_subsampled_plv)
                this_pct(iF) = nan;
            else
                this_pct(iF) = sum(trial_subsampled_plv(iF) > shuf_plv(:,iF))/length(shuf_plv);
            end
            % Change accordingly
            this_nbins = 5; %length(out.fr.bin{iF}); %5;
            phase_bins = -pi:2*pi/this_nbins:pi;
            this_phase = causal_phase(iF,goodTrials(1):goodTrials(2));
            [this_count, ~, ~] = histcounts(this_phase, phase_bins);
            this_trialbin_dif(iF) = (max(this_count) - min(this_count))/(goodTrials(2) + 1 - goodTrials(1));
            % The maximally excitable phase
            [~, midx] =  max(out.fr.bin(iF,:));%max(out.fr.bin{iF});
            this_ex_phase(iF) = mean(phase_bins(midx:midx+1));
        end
        if isempty(trial_subsampled_plv)
            s_out.phaselock_pct = [s_out.phaselock_pct; this_pct];
            s_out.phaselock_plv = [s_out.phaselock_plv; this_pct];
            s_out.phaselock_max_shufplv = [s_out.phaselock_max_shufplv; this_pct];
            s_out.phaselock_mean_phase = [s_out.phaselock_mean_phase; this_pct];
        else
            s_out.phaselock_pct = [s_out.phaselock_pct; this_pct];
            s_out.phaselock_plv = [s_out.phaselock_plv; trial_subsampled_plv'];
            s_out.phaselock_max_shufplv = [s_out.phaselock_max_shufplv; max(shuf_plv, [], 1)];
            s_out.phaselock_mean_phase = [s_out.phaselock_mean_phase; trial_subsampled_mean_phase'];
        end
        s_out.excitable_phase = [s_out.excitable_phase; this_ex_phase];
        s_out.nbins = [s_out.nbins; this_nbins];
        s_out.trialbin_dif = [s_out.trialbin_dif; this_trialbin_dif];
    end
end
