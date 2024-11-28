%% Script to generate various scatter summary plots
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.depth, summary.stim_50_pct, summary.stim_55_pct, ...
    summary.stim_60_pct, summary.stim_65_pct, summary.stim_70_pct, ...
    summary.stim_75_pct, summary.stim_90_pct] = deal([]);

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

%%

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

        % Load the good stim results
        load('good_stim.mat')
        this_stim = cellfun(@(x) length(x), good_stim);
        s_out.stim_50_pct = [s_out.stim_50_pct; this_stim(:,1)'];
        s_out.stim_55_pct = [s_out.stim_55_pct; this_stim(:,2)'];
        s_out.stim_60_pct = [s_out.stim_60_pct; this_stim(:,3)'];
        s_out.stim_65_pct = [s_out.stim_65_pct; this_stim(:,4)'];
        s_out.stim_70_pct = [s_out.stim_70_pct; this_stim(:,5)'];
        s_out.stim_75_pct = [s_out.stim_75_pct; this_stim(:,6)'];
        s_out.stim_90_pct = [s_out.stim_90_pct; this_stim(:,7)'];
    end
end

