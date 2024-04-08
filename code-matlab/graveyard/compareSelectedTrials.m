%% Script to generate various scatter summary plots
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.depth, summary.fr_r, summary.fr_z, ...
    summary.fr_r65, summary.fr_z65, summary.fr_r60, ...
    summary.fr_z60, summary.num_trials, summary.num_trials65, ...
    summary.num_trials60, summary.ratio_chance65, summary.ratio_chance60, ...
    summary.zscore_chance65, summary.zscore_chance60] = deal([]);

for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end

fbands = {[2 5], [6 10], [12,28], [30 55]};
c_list = {'red', 'blue','magenta','cyan'};

% Load the list of final opto cells and keep the results from only those
load('data\FinalOptoCells.mat');
keep = contains(summary.labels, dStr_opto) | contains(summary.labels, vStr_opto);
fn = fieldnames(summary);
for i = 1:numel(fn)
    temp = summary.(fn{i});
    summary.(fn{i}) = temp(keep,:);
end
clear fn temp

dStr_mask = (contains(summary.labels, dStr_opto) &  summary.depth < 3.5);
vStr_mask = (contains(summary.labels, vStr_opto) &  summary.depth >= 3.5);

%% Plot stuff for 60 pct comparison
fig = figure('WindowState', 'maximized');
min_trials = 600; 
for iF = 1:length(fbands)
    keep = (vStr_mask | dStr_mask) & (summary.num_trials60(:,iF) >= min_trials);

    % Plot raw score
    ax = subplot(3,4,iF);
    hold on
    y_og = summary.fr_r(keep,iF);
    y_new = summary.fr_r60(keep,iF);
    scatter(ones(size(y_og)), y_og, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0);
    scatter(2*ones(size(y_new)), y_new, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0);
    for iC = 1:length(y_og)
        plot([1,2],[y_og(iC),y_new(iC)], c_list{iF},  'LineWidth', 0.2)
    end
    plot([1,2], [mean(y_og), mean(y_new)], 'black', 'LineWidth', 2)
%     ylim([0 0.6])
    xticks([1 2])
    xlim([0.8 2.2])
    title(sprintf('%d - %d Hz (n = %d)', fbands{iF}(1), fbands{iF}(2), length(y_og)), 'FontSize', 18);
    ylabel('FR mod. raw', 'FontSize', 14);
    xticklabels({'Original', 'Top60'});
    ax.XAxis.FontSize = 14;
    fprintf('Mean raw FR in %d- %d Hz for all = %.3f, for top-60 = %.3f\n', ...
        fbands{iF}(1), fbands{iF}(2), mean(y_og), mean(y_new));
    [h,p] = kstest2(y_og,y_new,'Tail','smaller');
    fprintf('p-val = %.3f\n', p);

    % Plot Z-scores
    ax = subplot(3,4,iF+4);
    hold on
    y_og = summary.fr_z(keep,iF);
    y_new = summary.fr_z60(keep,iF);
    scatter(ones(size(y_og)), y_og, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0);
    scatter(2*ones(size(y_new)), y_new, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0);
    for iC = 1:length(y_og)
        plot([1,2],[y_og(iC),y_new(iC)], c_list{iF},  'LineWidth', 0.2)
    end
    plot([1,2], [mean(y_og), mean(y_new)], 'black', 'LineWidth', 2)
%     ylim([0 0.6])
    xticks([1 2])
    xlim([0.8 2.2])
    ylabel('FR mod. z-scored', 'FontSize', 14);
    xticklabels({'Original', 'Top60'});
    ax.XAxis.FontSize = 14;
    fprintf('Mean z-scored FR in %d- %d Hz for all = %.3f, for top-60 = %.3f\n', ...
        fbands{iF}(1), fbands{iF}(2), mean(y_og), mean(y_new));
    [h,p] = kstest2(y_og,y_new,'Tail','smaller');
    fprintf('p-val = %.3f\n', p)

    % Plot proportion of neurons significant in each case
    ax = subplot(3,4,iF+8);
    hold on
    this_proportion = [sum(summary.fr_z((vStr_mask|dStr_mask),iF)>=2)/sum(vStr_mask|dStr_mask), ...
        sum(summary.fr_z60(keep,iF)>=2)/sum(keep)];
    hb = bar([1 2],100*this_proportion, 'FaceColor', c_list{iF}, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    text(hb.XData(1) + hb.XOffset, hb.YData(1), sprintf('%d/%d', ...
        sum(summary.fr_z((vStr_mask|dStr_mask),iF)>=2),sum(vStr_mask|dStr_mask)), ...
        'VerticalAlignment','bottom', 'HorizontalAlignment', 'center', 'FontSize', 16);
    text(hb.XData(2) + hb.XOffset, hb.YData(2), sprintf('%d/%d', ...
        sum(summary.fr_z60(keep,iF)>=2),sum(keep)), ...
        'VerticalAlignment','bottom', 'HorizontalAlignment', 'center', 'FontSize', 16);
    xticks([1 2])
    xticklabels({'Original', 'Top60'});
    ylabel('% significant', 'FontSize', 14);
    ax.XAxis.FontSize = 14;
end

sgtitle('Comparing trials within top 60 percentile')

%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    evs = LoadEvents([]);
    cfg = [];
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.min_cluster_quality = 3;
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
    cfg.fc = ExpKeys.goodCell;

    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');

        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];

        % Load the og results and selected trial results
        load(strcat(fn_prefix, '_phase_response_5_bins.mat'));
        load(strcat(fn_prefix, '_selected_trials_phase_response_5_bins.mat'));

        s_out.fr_z = [s_out.fr_z; out.fr.zscore];
        s_out.fr_r = [s_out.fr_r; out.fr.ratio];

        s_out.fr_z60 = [s_out.fr_z60; out_top.fr.zscore(1,:)];
        s_out.fr_r60 = [s_out.fr_r60; out_top.fr.ratio(1,:)];

        s_out.fr_z65 = [s_out.fr_z65; out_top.fr.zscore(2,:)];
        s_out.fr_r65 = [s_out.fr_r65; out_top.fr.ratio(2,:)];

        s_out.num_trials = [s_out.num_trials; diff(ExpKeys.goodTrials(iC,:))+1];
        s_out.num_trials60 = [s_out.num_trials60; out_top.num_trials(1,:)];
        s_out.num_trials65 = [s_out.num_trials65; out_top.num_trials(2,:)];
        
        % Calculating if the difference between the all trials and selected
        % trials is more than expected by chance
        s_out.ratio_chance60 = [s_out.ratio_chance60; ...
            sum(out.fr.ratio >= squeeze(out_top.control.fr.ratio(1,:,:)))];
        s_out.ratio_chance65 = [s_out.ratio_chance65; ...
            sum(out.fr.ratio >= squeeze(out_top.control.fr.ratio(2,:,:)))];
        s_out.zscore_chance60 = [s_out.zscore_chance60; ...
            sum(out.fr.zscore >= squeeze(out_top.control.fr.zscore(1,:,:)))];
        s_out.zscore_chance65 = [s_out.zscore_chance65; ...
            sum(out.fr.zscore >= squeeze(out_top.control.fr.zscore(2,:,:)))];

    end
end

