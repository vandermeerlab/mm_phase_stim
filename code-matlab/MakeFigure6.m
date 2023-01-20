%% Script to generate various scatter summary plots
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.stim_mode, summary.short_stim, ...
    summary.long_stim, summary.response_p, summary.depth,...
    summary.fr_r, summary.fr_z] = deal([]);
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
%     fbands = {[2 5], [6 10], [12 30], [30 55]};
    fbands = {[2 5], [6 10], [30 55]};
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

    end
end

