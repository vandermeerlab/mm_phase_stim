%% Script to generate various scatter summary plots
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.bfr, summary.bfr10, summary.bfr250] = deal({});
[summary.labels, summary.stim_mode, summary.short_stim, ...
    summary.long_stim, summary.response_p, summary.depth,...
    summary.fr_r, summary.fr_z,  summary.phaselock_plv, ...
    summary.fr_r10, summary.fr_z10, summary.fr_r250, summary.fr_z250, ...
    summary.phaselock_mean_phase, summary.phaselock_pct, summary.phaselock_z, ...
    summary.phaselock_max_shufplv, summary.phaselock_circ_pct, summary.phaselock_circ_z, ...
    summary.phaselock_max_circ_shufplv, summary.excitable_phase, summary.ntrials, ...
    summary.excitable_phase10, summary.excitable_phase250, summary.trialbin_dif, summary.nbins] = deal([]);
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

% This step is important because originally there were 4 frequency bands, but now we
% are skipping the 3rd
fn = fieldnames(summary);
for i = 1:numel(fn)
   if size(summary.(fn{i}),2) == 4
        temp = summary.(fn{i});
        temp(:,3) = temp(:,4);
        summary.(fn{i}) = temp(:,1:3);
   end
end



% Load the list of final opto cells
load('E:\Dropbox (Dartmouth College)\AnalysisResults\phase_stim_results\FinalOptoCells.mat');
dStr_mask = (contains(summary.labels, dStr_opto) &  summary.depth < 3.5);
vStr_mask = (contains(summary.labels, vStr_opto) &  summary.depth >= 3.5);

% Another exclusion criterion is if the difference between the 2 trial counts for a given phase binning is greater than 5% of the trials
% max_dif = 0.05; % Subject to change
% dif_mask = summary.trialbin_dif <= max_dif;

% Significance mask
z_thresh = 2;
sig_mask = (summary.fr_z > z_thresh);
%% Figure3: Create CSV for UPSET Plot
delta_bool = [0;1;0;0;1;1;0;1];
theta_bool = [0;0;1;0;1;0;1;1];
gamma_bool = [0;0;0;1;0;1;1;1];
delta_bool = delta_bool == 1;
theta_bool = theta_bool == 1;
gamma_bool = gamma_bool == 1;

keep = sig_mask;
[dStr_sig, vStr_sig] = deal(zeros(size(delta_bool)));


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

dStr_table = table(delta_bool, theta_bool, gamma_bool, dStr_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz','30-55 Hz', 'Count'});
vStr_table = table(delta_bool, theta_bool, gamma_bool, vStr_sig, 'VariableNames', {'2-5 Hz', '6-10 Hz','30-55 Hz', 'Count'});
writetable(dStr_table, 'C:\Users\mvdmlab\Desktop\dStr_sig.csv');
writetable(vStr_table, 'C:\Users\mvdmlab\Desktop\vStr_sig.csv');



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
    ylim([-0.05 0.5])
    xlim([2 5])
    xlabel('Recording depth (mm)')
    ylabel('Modulation strength')
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
    ylim([-3 12])
    xlim([2 5])
    xlabel('Recording depth (mm)')
    ylabel('z-score')
    yline(2, '--black')
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
end
fontname(fig, 'Helvetica')
fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts


%% Figure 4: Input to the Venn diagrams to show phase dependent excitability and phase locking
pl_thresh = 0.95;
z_thresh = 2;
pl_mask = summary.phaselock_pct >= pl_thresh;
pl_zmask = summary.phaselock_z >= z_thresh;
pl_circ_mask = summary.phaselock_circ_pct >= pl_thresh;
pl_circ_zmask = summary.phaselock_circ_z >= z_thresh;

[sum(dStr_mask), sum(dStr_mask & sig_mask(:,1)), sum(dStr_mask & pl_zmask(:,1))],...
[sum(dStr_mask & sig_mask(:,1)), sum(dStr_mask & pl_zmask(:,1)), sum(pl_zmask(dStr_mask,1) & sig_mask(dStr_mask,1)), ...
sum(dStr_mask & sig_mask(:,1) & pl_zmask(:,1))]

[sum(dStr_mask), sum(dStr_mask & sig_mask(:,2)), sum(dStr_mask & pl_zmask(:,2))],...
[sum(dStr_mask & sig_mask(:,2)), sum(dStr_mask & pl_zmask(:,2)), sum(pl_zmask(dStr_mask,2) & sig_mask(dStr_mask,2)), ...
sum(dStr_mask & sig_mask(:,2) & pl_zmask(:,2))]

[sum(dStr_mask), sum(dStr_mask & sig_mask(:,3)), sum(dStr_mask & pl_zmask(:,3))],...
[sum(dStr_mask & sig_mask(:,3)), sum(dStr_mask & pl_zmask(:,3)), sum(pl_zmask(dStr_mask,3) & sig_mask(dStr_mask,3)), ...
sum(dStr_mask & sig_mask(:,3) & pl_zmask(:,3))]

[sum(vStr_mask), sum(vStr_mask & sig_mask(:,1)), sum(vStr_mask & pl_zmask(:,1))],...
[sum(vStr_mask & sig_mask(:,1)), sum(vStr_mask & pl_zmask(:,1)), sum(pl_zmask(vStr_mask,1) & sig_mask(vStr_mask,1)), ...
sum(vStr_mask & sig_mask(:,1) & pl_zmask(:,1))]

[sum(vStr_mask), sum(vStr_mask & sig_mask(:,2)), sum(vStr_mask & pl_zmask(:,2))],...
[sum(vStr_mask & sig_mask(:,2)), sum(vStr_mask & pl_zmask(:,2)), sum(pl_zmask(vStr_mask,2) & sig_mask(vStr_mask,2)), ...
sum(vStr_mask & sig_mask(:,2) & pl_zmask(:,2))]

[sum(vStr_mask), sum(vStr_mask & sig_mask(:,3)), sum(vStr_mask & pl_zmask(:,3))],...
[sum(vStr_mask & sig_mask(:,3)), sum(vStr_mask & pl_zmask(:,3)), sum(pl_zmask(vStr_mask,3) & sig_mask(vStr_mask,3)), ...
sum(vStr_mask & sig_mask(:,3) & pl_zmask(:,3))]

%% Figure 4: For neurons that are both phase locked and significantly phase modulated, how do they show up in the control_stim terst
z_thresh = 2;
sig_10 = summary.fr_z10 >= z_thresh;
sig_250 = summary.fr_z250 >= z_thresh;
pl_mask = summary.phaselock_z >= z_thresh;

fprintf('Dorsal Striatum\n');
for iF = 1:length(fbands)
    keep = find(sig_mask(:,iF) & pl_mask(:,iF) & dStr_mask); 
    fprintf('%d - %d Hz\n', fbands{iF}(1),  fbands{iF}(2));
    fprintf('Label\t\t\t\t\tExPhase\t\tSig10\t\tExPhase10\t\tSig250\t\tExPhase250\n');
    for iC = 1:length(keep)
        fprintf('%s\t %d\t\t\t %d\t\t\t %d\t\t\t\t %d\t\t\t\t %d\n', summary.labels(keep(iC)), ...
           summary.excitable_phase(keep(iC),iF), sig_10(keep(iC),iF), summary.excitable_phase10(keep(iC),iF), ...
           sig_250(keep(iC),iF), summary.excitable_phase250(keep(iC),iF));
    end
end

fprintf('Ventral Striatum\n')
for iF = 1:length(fbands)
    keep = find(sig_mask(:,iF) & pl_mask(:,iF) & vStr_mask); 
    fprintf('%d - %d Hz\n', fbands{iF}(1),  fbands{iF}(2));
    fprintf('Label\t\t\t\t\tExPhase\t\tSig10\t\tExPhase10\t\tSig250\t\tExPhase250\n');
    for iC = 1:length(keep)
        fprintf('%s\t %d\t\t\t %d\t\t\t %d\t\t\t\t %d\t\t\t\t %d\n', summary.labels(keep(iC)), ...
           summary.excitable_phase(keep(iC),iF), sig_10(keep(iC),iF), summary.excitable_phase10(keep(iC),iF), ...
           sig_250(keep(iC),iF), summary.excitable_phase250(keep(iC),iF));
    end
end



%% Diagnostic Plot: Look at proportions of significant phase-locking depending on what thresholding/shuffle method was used 
pl_thresh = 0.95;
z_thresh = 2;
pl_mask = summary.phaselock_pct >= pl_thresh;
pl_zmask = summary.phaselock_z >= z_thresh;
pl_circ_mask = summary.phaselock_circ_pct >= pl_thresh;
pl_circ_zmask = summary.phaselock_circ_z >= z_thresh;

fprintf("Fbands:\t\t\t\t\t\t\t\t2 - 5 Hz\t6 - 10 Hz\t 30 - 55Hz\n");
fprintf("dStr Sig Phase Locked:\t\t\t\t %d / %d\t %d / %d\t %d / %d\n", ...
    sum(pl_mask(dStr_mask,1)), sum(dStr_mask), sum(pl_mask(dStr_mask,2)), sum(dStr_mask), ...
        sum(pl_mask(dStr_mask,3)), sum(dStr_mask));
fprintf("vStr Sig Phase Locked:\t\t\t\t %d / %d\t %d / %d\t %d / %d\n", ...
    sum(pl_mask(vStr_mask,1)), sum(vStr_mask), sum(pl_mask(vStr_mask,2)), sum(vStr_mask), ...
        sum(pl_mask(vStr_mask,3)), sum(vStr_mask));

fprintf("dStr Sig Phase Locked(z):\t\t\t %d / %d\t %d / %d\t %d / %d\n", ...
    sum(pl_zmask(dStr_mask,1)), sum(dStr_mask), sum(pl_zmask(dStr_mask,2)), sum(dStr_mask), ...
        sum(pl_zmask(dStr_mask,3)), sum(dStr_mask));
fprintf("vStr Sig Phase Locked(z):\t\t\t %d / %d\t %d / %d\t %d / %d\n", ...
    sum(pl_zmask(vStr_mask,1)), sum(vStr_mask), sum(pl_zmask(vStr_mask,2)), sum(vStr_mask), ...
        sum(pl_zmask(vStr_mask,3)), sum(vStr_mask));

fprintf("dStr Sig Phase Locked(circ):\t\t %d / %d\t\t %d / %d\t\t %d / %d\n", ...
    sum(pl_circ_mask(dStr_mask,1)), sum(dStr_mask), sum(pl_circ_mask(dStr_mask,2)), sum(dStr_mask), ...
        sum(pl_circ_mask(dStr_mask,3)), sum(dStr_mask));
fprintf("vStr Sig Phase Locked(circ):\t\t %d / %d\t\t %d / %d\t\t %d / %d\n", ...
    sum(pl_circ_mask(vStr_mask,1)), sum(vStr_mask), sum(pl_circ_mask(vStr_mask,2)), sum(vStr_mask), ...
        sum(pl_circ_mask(vStr_mask,3)), sum(vStr_mask));

fprintf("dStr Sig Phase Locked(circz):\t\t %d / %d\t\t %d / %d\t\t %d / %d\n", ...
    sum(pl_circ_zmask(dStr_mask,1)), sum(dStr_mask), sum(pl_circ_zmask(dStr_mask,2)), sum(dStr_mask), ...
        sum(pl_circ_zmask(dStr_mask,3)), sum(dStr_mask));
fprintf("vStr Sig Phase Locked(circz):\t\t %d / %d\t\t %d / %d\t\t %d / %d\n", ...
    sum(pl_circ_zmask(vStr_mask,1)), sum(vStr_mask), sum(pl_circ_zmask(vStr_mask,2)), sum(vStr_mask), ...
        sum(pl_circ_zmask(vStr_mask,3)), sum(vStr_mask));

%% Diagnostic Plot: Z-score and Mod-depth plot for neurons that are both
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ax = subplot(2,3,iF);
    hold on
    scatter(summary.depth(dStr_mask & sig_mask(:,iF)),  summary.fr_r(dStr_mask & sig_mask(:,iF),iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
    scatter(summary.depth(dStr_mask & sig_mask(:,iF) & pl_mask(:,iF)), summary.fr_r(dStr_mask & sig_mask(:, iF) & pl_mask(:,iF),iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 50);
    scatter(summary.depth(vStr_mask & sig_mask(:,iF)), summary.fr_r(vStr_mask & sig_mask(:,iF),iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    scatter(summary.depth(vStr_mask & sig_mask(:,iF) & pl_mask(:,iF)), summary.fr_r(vStr_mask & sig_mask(:,iF) & pl_mask(:,iF),iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'Marker', 'd', 'SizeData', 50);
    ylim([-0.05 0.5])
    xlim([2 5])
    xlabel('Recording depth (mm)')
    ylabel('Modulation strength')
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';

    ax = subplot(2,3,iF+3);
    hold on
     scatter(summary.depth(dStr_mask & sig_mask(:,iF)),  summary.fr_z(dStr_mask & sig_mask(:,iF),iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
    scatter(summary.depth(dStr_mask & sig_mask(:,iF) & pl_mask(:,iF)), summary.fr_z(dStr_mask & sig_mask(:, iF) & pl_mask(:,iF),iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'SizeData', 50);
    scatter(summary.depth(vStr_mask & sig_mask(:,iF)), summary.fr_z(vStr_mask & sig_mask(:,iF),iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    scatter(summary.depth(vStr_mask & sig_mask(:,iF) & pl_mask(:,iF)), summary.fr_z(vStr_mask & sig_mask(:,iF) & pl_mask(:,iF),iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'Marker', 'd', 'SizeData', 50);
    ylim([-3 12])
    xlim([2 5])
    xlabel('Recording depth (mm)')
    ylabel('z-score')
    ax.TickDir = 'out';
    ax.TickLength(1) = 0.03;
    ax.Box = 'off';
end

%% Diagnostic plot: Scatter plot of baseline firing rate vs depth of modulation
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

%% Diagnostic Plot: Scatter of most excitable phases vs depth of modulationfor significant results
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

% vStr
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

%% Diagnostic Plot: Look at why high depth of modulation is not necessarily significant

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

% vStr
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


%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    fbands = {[2 5], [6 10], [12 30], [30 55]};
    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');

        % Load the phase_response
        load(strcat(fn_prefix, '_phase_response_5_bins.mat')); % Change this to what is decided to be the best binning option
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];
        s_out.response_p = [s_out.response_p; out.overall_response];
        s_out.bfr{length(s_out.bfr)+1} = out.bfr;
        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.stim_mode = [s_out.stim_mode; string(ExpKeys.stim_mode)];
        s_out.fr_z = [s_out.fr_z; out.fr.zscore];
        s_out.fr_r = [s_out.fr_r; out.fr.ratio];
        s_out.short_stim = [s_out.short_stim; ExpKeys.short_stim_pulse_width];
        s_out.long_stim = [s_out.long_stim; ExpKeys.long_stim_pulse_width];
        

        % Load the control phase_response
        load(strcat(fn_prefix, '_control_phase_response_5_bins.mat'));
        s_out.bfr10{length(s_out.bfr10)+1} = out10.bfr;
        s_out.fr_z10 = [s_out.fr_z10; out10.fr.zscore];
        s_out.fr_r10 = [s_out.fr_r10; out10.fr.ratio];
        s_out.bfr250{length(s_out.bfr250)+1} = out250.bfr;
        s_out.fr_z250 = [s_out.fr_z250; out250.fr.zscore];
        s_out.fr_r250 = [s_out.fr_r250; out250.fr.ratio];
        
        % Load the stim_responses
        load('stim_phases.mat');
        goodTrials = ExpKeys.goodTrials(iC,:);

        s_out.ntrials = [s_out.ntrials; goodTrials(2) + 1 - goodTrials(1)];
        
        % Load phase_lock and shuf_spec
        fn_prefix = strrep(fn_prefix, '_', '-');
        load(strcat(fn_prefix, '_spike_phaselock_plv.mat'));
        load(strcat(fn_prefix, '_shuf_spec_circ_plv.mat')); % First circularly shifted and then subsampled
%         load(strcat(fn_prefix, '_shuf_spec_circ2_plv.mat')); % Uniform
        load(strcat(fn_prefix, '_shuf_spec_plv.mat'));  % Uniformly distributed fake spikes
        [this_pct, this_circ_pct, this_circ2_pct, this_ex_phase, this_trialbin_dif, ...
            this_z, this_circ_z, this_circ2_z, this_ex_phase10, this_ex_phase250] = deal(zeros(1,length(fbands)));
        for iF = 1:length(fbands)
            if isempty(trial_subsampled_plv)
                this_pct(iF) = nan;
                this_z(iF) = nan;
                this_circ_pct(iF) = nan;
                this_circ_z(iF) = nan;
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
            % Change accordingly
            this_nbins = 5; %length(out.fr.bin{iF}); %5;
            phase_bins = -pi:2*pi/this_nbins:pi;
            this_phase = causal_phase(iF,goodTrials(1):goodTrials(2));
            [this_count, ~, ~] = histcounts(this_phase, phase_bins);
            this_trialbin_dif(iF) = (max(this_count) - min(this_count))/(goodTrials(2) + 1 - goodTrials(1));
            % The maximally excitable phase
            [~, midx] =  max(out.fr.bin(iF,:));%max(out.fr.bin{iF});
            this_ex_phase(iF) = midx;

            % The maximally excitable phase in the control conditions
            [~, midx] =  max(out10.fr.bin(iF,:));%max(out.fr.bin{iF});
            this_ex_phase10(iF) = midx;
            [~, midx] =  max(out250.fr.bin(iF,:));%max(out.fr.bin{iF});
            this_ex_phase250(iF) = midx;
            
        end
        if isempty(trial_subsampled_plv) %this_pct is all nans in this case
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
           
            s_out.phaselock_plv = [s_out.phaselock_plv; trial_subsampled_plv'];
            s_out.phaselock_mean_phase = [s_out.phaselock_mean_phase; trial_subsampled_mean_phase'];
        end
        s_out.excitable_phase = [s_out.excitable_phase; this_ex_phase];
        s_out.excitable_phase10 = [s_out.excitable_phase10; this_ex_phase10];
        s_out.excitable_phase250 = [s_out.excitable_phase250; this_ex_phase250];

        s_out.nbins = [s_out.nbins; this_nbins];
        s_out.trialbin_dif = [s_out.trialbin_dif; this_trialbin_dif];
    end
end

