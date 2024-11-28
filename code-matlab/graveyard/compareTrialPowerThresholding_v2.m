%% Script to compare the effect of amplitude based thresholding of trials
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.stim_mode, summary.depth,...
    summary.fr_r, summary.fr_z, ...
    summary.top_fr_r, summary.top_fr_z, ...
    summary.bottom_fr_r, summary.bottom_fr_z, ...
    summary.ratio_pct, summary.zscore_pct] = deal([]);
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
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    ycoords1 = summary.top_fr_r(vStr_mask | dStr_mask,iF);
    ycoords2 = summary.bottom_fr_r(vStr_mask | dStr_mask,iF);
    ax = subplot(2,4,iF);
    hold on
%     bar(1, mean(ycoords1), 0.5, 'FaceColor', c_list{iF}, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    scatter(ones(size(ycoords1)), ycoords1, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0);
%     bar(2, mean(ycoords2), 0.5,  'FaceColor', c_list{iF}, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    scatter(2*ones(size(ycoords1)), ycoords2, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0);
    for iC = 1:length(ycoords1)
        plot([1,2],[ycoords1(iC),ycoords2(iC)], c_list{iF},  'LineWidth', 0.2)
    end
    plot([1,2], [mean(ycoords1), mean(ycoords2)], 'black', 'LineWidth', 4)
    ylim([0 1])
    xticks([1 2])
    xlim([0.6 2.4])
    xticklabels({'Top-50', 'Bottom-50'})
    ylabel('FR mod. raw', 'FontSize', 14);
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 16);
    ax.XAxis.FontSize = 14;
    fprintf('Mean raw FR in %d- %d Hz for top-50 = %.3f, for bottom-50 = %.3f\n', ...
        fbands{iF}(1), fbands{iF}(2), mean(ycoords1), mean(ycoords2));
    [h,p] = ttest2(ycoords1,ycoords2);
    fprintf('p-val = %.3f\n', p);
    % T-test result
    
    ycoords1 = summary.top_fr_z(vStr_mask | dStr_mask,iF);
    ycoords2 = summary.bottom_fr_z(vStr_mask | dStr_mask,iF);
    ax = subplot(2,4,iF+4);
    hold on
%     bar(1, mean(ycoords1), 0.5, 'FaceColor', c_list{iF}, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    scatter(ones(size(ycoords1)), ycoords1, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0);
%     bar(2, mean(ycoords2), 0.5,  'FaceColor', c_list{iF}, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    scatter(2*ones(size(ycoords1)), ycoords2, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0);
    for iC = 1:length(ycoords1)
        plot([1,2],[ycoords1(iC),ycoords2(iC)], c_list{iF},  'LineWidth', 0.2)
    end
    plot([1,2], [mean(ycoords1), mean(ycoords2)], 'black', 'LineWidth', 4)
    ylim([-3 5])
    xticks([1 2])
    xlim([0.6 2.4])
    xticklabels({'Top-50', 'Bottom-50'})
    ylabel('FR mod. z-scored', 'FontSize', 14);
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 16);
    ax.XAxis.FontSize = 14;
    fprintf('Mean z-scored FR in %d- %d Hz for top-50 = %.3f, for bottom-50 = %.3f\n', ...
        fbands{iF}(1), fbands{iF}(2), mean(ycoords1), mean(ycoords2));
    [h,p] = ttest2(ycoords1,ycoords2);
    fprintf('p-val = %.3f\n', p);

%     ycoords3 = 
%     ax = subplot(3,4,iF+8);
%     hold on
%     plot(xcoords, ycoords3, 'Color', c_list{iF});
%     scatter(xcoords, ycoords3, 'MarkerFaceColor', c_list{iF}, ...
%             'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
%     xticks(xcoords);
%     xlim([45,105]);
%     xlabel('Percentage of trials', 'FontSize', 14);
%     ylabel('Significant cell count', 'FontSize', 14);
%     title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)), 'FontSize', 18);
end
%%
% For neurons that were significant, was the difference more than expected
% by chance

% Significance mask
z_thresh = 2;
sig_mask = (summary.fr_z > z_thresh);

sum(sig_mask & (dStr_mask | vStr_mask)  & summary.zscore_pct < 0.05)


%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');

        % Load the thresholding based results
        load(strcat(fn_prefix, '_split_by_power_phase_response_5_bins.mat'));
        out1 = load(strcat(fn_prefix, '_phase_response_5_bins.mat'));
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];
        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.stim_mode = [s_out.stim_mode; string(ExpKeys.stim_mode)];

        s_out.fr_z = [s_out.fr_z; out1.out.fr.zscore];
        s_out.fr_r = [s_out.fr_r; out1.out.fr.ratio];
        s_out.top_fr_z = [s_out.top_fr_z; out_top.fr.zscore];
        s_out.top_fr_r = [s_out.top_fr_r; out_top.fr.ratio];
        s_out.bottom_fr_z = [s_out.bottom_fr_z; out_bottom.fr.zscore];
        s_out.bottom_fr_r = [s_out.bottom_fr_r; out_bottom.fr.ratio];

        shuf_ratio_diff = out_top_control.fr.ratio - out_bottom_control.fr.ratio;
        shuf_z_diff = out_top_control.fr.zscore - out_bottom_control.fr.zscore;
        ratio_diff = out_top.fr.ratio - out_bottom.fr.ratio;
        z_diff = out_top.fr.zscore - out_bottom.fr.zscore;

        % Probability that difference was greater than expected by chance
        s_out.ratio_pct = [s_out.ratio_pct; sum(ratio_diff < shuf_ratio_diff)/1000];
        s_out.zscore_pct = [s_out.zscore_pct; sum(z_diff < shuf_z_diff)/1000];

    end
end

