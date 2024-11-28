%% Script to generate various scatter summary plots
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.stim_mode, summary.short_sw, ...
    summary.long_sw, summary.depth, ...
    summary.pre_stim, summary.trial_stim, summary.post_stim, ...
    summary.long_stim, summary.ntrials] = deal([]);
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

% Load the list of final non-opto cells and keep the results from only those
load('data\FinalNonOptoCells.mat');
keep = contains(summary.labels, dStr_others) | contains(summary.labels, vStr_others);
fn = fieldnames(summary);
for i = 1:numel(fn)
    temp = summary.(fn{i});
    summary.(fn{i}) = temp(keep,:);
end
clear fn temp

dStr_mask = (contains(summary.labels, dStr_others) &  summary.depth < 3.5);
vStr_mask = (contains(summary.labels, vStr_others) &  summary.depth >= 3.5);

%% Figure without separating into dStr and vStr
keep = dStr_mask | vStr_mask;

fig = figure('WindowState', 'maximized');
ax = subplot(1,4,1);
imagesc([-20:1:20], 1:sum(keep), summary.pre_stim(keep,:));
yticks([]);
xticks([-20:10:20]);
ylabel('Non Opto cells')
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Pre-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(1,4,2);
imagesc([-20:1:20], 1:sum(keep), summary.trial_stim(keep,:));
yticks([])
xticks([-20:10:20]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Trial-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(1,4,3);
imagesc([-20:1:20], 1:sum(keep), summary.post_stim(keep,:));
yticks([]);
xticks([-20:10:20]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Post-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(1,4,4);
imagesc([-200:10:200], 1:sum(keep), summary.long_stim(keep,:));
yticks([]);
xticks([-200:100:200]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Long-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

%%
fig = figure('WindowState', 'maximized');

% dStr stuff
ax = subplot(2,4,1);
imagesc([-20:1:20], 1:sum(dStr_mask), summary.pre_stim(dStr_mask,:));
yticks([]);
xticks([-20:10:20]);
ylabel('dStr')
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Pre-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(2,4,2);
imagesc([-20:1:20], 1:sum(dStr_mask), summary.trial_stim(dStr_mask,:));
yticks([])
xticks([-20:10:20]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Trial-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(2,4,3);
imagesc([-20:1:20], 1:sum(dStr_mask), summary.post_stim(dStr_mask,:));
yticks([]);
xticks([-20:10:20]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Post-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(2,4,4);
imagesc([-200:10:200], 1:sum(dStr_mask), summary.long_stim(dStr_mask,:));
yticks([]);
xticks([-200:100:200]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Long-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

% vStr stuff
ax = subplot(2,4,5);
imagesc([-20:1:20], 1:sum(vStr_mask), summary.pre_stim(vStr_mask,:));
yticks([]);
xticks([-20:10:20]);
ylabel('vStr')
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Pre-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(2,4,6);
imagesc([-20:1:20], 1:sum(vStr_mask), summary.trial_stim(vStr_mask,:));
yticks([]);
xticks([-20:10:20]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Trial-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(2,4,7);
imagesc([-20:1:20], 1:sum(vStr_mask), summary.post_stim(vStr_mask,:));
yticks([]);
xticks([-20:10:20]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Post-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(2,4,8);
imagesc([-200:10:200], 1:sum(vStr_mask), summary.long_stim(vStr_mask,:));
yticks([]);
xticks([-200:100:200]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Long-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

fontname(fig, 'Helvetica');
fig.Renderer = 'painters';
%% Diagnosis figure

ax = subplot(1,2,1);
imagesc([-20:1:20], 1:sum(dStr_mask), summary.trial_stim(dStr_mask,:));
% yticks([]);
yticks(1:sum(dStr_mask))
yticklabels(summary.labels(dStr_mask));
xticks([-20:10:20]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Trial-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

ax = subplot(1,2,2);
imagesc([-20:1:20], 1:sum(vStr_mask), summary.trial_stim(vStr_mask,:));
% yticks([]);
yticks(1:sum(vStr_mask))
yticklabels(summary.labels(vStr_mask));
xticks([-20:10:20]);
xlabel('Time (ms)')
xline(0, '--red', 'LineWidth', 2)
title('Trial-stim')
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.03 0.02];

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
%     cfg.fc = ExpKeys.goodCell;
   
    S = LoadSpikes(cfg);
    % Reject goodCell and continue if other cells in the session
    S = SelectTS([],S,~ismember(S.label, ExpKeys.goodCell));
    if length(S.label) ~= 0
        % Set Variables
        max_delay = 0.01; % sec (window for the first response since stimulus)
        max_long_delay = 0.4; % sec
        short_win = 0.02; % sec
        long_win = 0.2; % sec
        short_dt = 0.001; % sec
        long_dt = 0.01; % sec
        fbands = {[2 5], [6 10], [30 55]};
    
        if contains(ExpKeys.light_source, 'LASER')
            start_delay = 0.0011;
            stop_delay = 0.0012; %0.0022 is the spike width in the case of 1 msec laser pulse
        else
            start_delay = 0;
            stop_delay = 0;
        end
    
        pre_stim_on = evs.t{strcmp(evs.label, ExpKeys.pre_trial_stim_on)} + start_delay;
        if ~isempty(pre_stim_on) && ~isempty(ExpKeys.pre_stim_times)
            pre_stim_on = pre_stim_on(pre_stim_on >= ExpKeys.pre_stim_times(1) & ...
                                      pre_stim_on <= ExpKeys.pre_stim_times(2));
        else
            pre_stim_on = [];
        end
    
        stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
        if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
            stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                              stim_on <= ExpKeys.stim_times(2));
        else
            stim_on = [];
        end
    
        post_stim_on = evs.t{strcmp(evs.label, ExpKeys.post_trial_stim_on)} + start_delay;
        if ~isempty(post_stim_on) && ~isempty(ExpKeys.post_stim_times)
            post_stim_on = post_stim_on(post_stim_on >= ExpKeys.post_stim_times(1) & ...
                                        post_stim_on <= ExpKeys.post_stim_times(2));
        else
            post_stim_on = [];
        end
    
        if sum(strcmp(evs.label, ExpKeys.long_stim_on)) ~= 0
            long_stim_on = evs.t{strcmp(evs.label, ExpKeys.long_stim_on)} + start_delay;
        else
            long_stim_on = [];
        end
        if ~isempty(long_stim_on) && ~isempty(ExpKeys.long_stim_times)
            long_stim_on = long_stim_on(long_stim_on >= ExpKeys.long_stim_times(1) & ...
                    long_stim_on <= ExpKeys.long_stim_times(2));
            % Take only the first 25 long stim if special stim present
            if strcmp(ExpKeys.hasSpecialStim, 'Yes')
                long_stim_on = long_stim_on(1:25); 
            end
        end
    
        for iC = 1:length(S.t)
            fn_prefix = extractBefore(S.label{iC}, '.t');
            this_cell = SelectTS([], S, iC);
            goodTrials = [min(min(ExpKeys.goodTrials)), max(max(ExpKeys.goodTrials))]; % Take all the trials for this
            
            % Bin the pre-stim stuff
            this_on_events = pre_stim_on;
            edges = -short_win:short_dt:short_win;
            bin_spks = zeros(length(this_on_events), length(edges)-1);
            for iT = 1:length(this_on_events)
                temp = restrict(this_cell, iv(this_on_events(iT) - short_win, ...
                    this_on_events(iT) + short_win));
                temp = temp.t{1} - this_on_events(iT); % center it around 0
                [count, ~, ~] = histcounts(temp, edges);
                bin_spks(iT,:) = count;
            end
            pre_bin_spks = mean(bin_spks);
            pre_bin_spks = (pre_bin_spks - min(pre_bin_spks))/...
                (max(pre_bin_spks) - min(pre_bin_spks));
            
            % Bin the trial stim stuff
            this_on_events = stim_on(goodTrials(1):goodTrials(2));
            edges = -short_win:short_dt:short_win;
            bin_spks = zeros(length(this_on_events), length(edges)-1);
            for iT = 1:length(this_on_events)
                temp = restrict(this_cell, iv(this_on_events(iT) - short_win, ...
                    this_on_events(iT) + short_win));
                temp = temp.t{1} - this_on_events(iT); % center it around 0
                [count, ~, ~] = histcounts(temp, edges);
                bin_spks(iT,:) = count;
            end
            trial_bin_spks = mean(bin_spks);
            trial_bin_spks = (trial_bin_spks - min(trial_bin_spks))/...
                (max(trial_bin_spks) - min(trial_bin_spks));
    
            % Bin the post-stim stuff
            if ~isempty(post_stim_on)
                this_on_events = post_stim_on;
                edges = -short_win:short_dt:short_win;
                bin_spks = zeros(length(this_on_events), length(edges)-1);
                for iT = 1:length(this_on_events)
                    temp = restrict(this_cell, iv(this_on_events(iT) - short_win, ...
                        this_on_events(iT) + short_win));
                    temp = temp.t{1} - this_on_events(iT); % center it around 0
                    [count, ~, ~] = histcounts(temp, edges);
                    bin_spks(iT,:) = count;
                end
                post_bin_spks = mean(bin_spks);
                post_bin_spks = (post_bin_spks - min(post_bin_spks))/...
                    (max(post_bin_spks) - min(post_bin_spks));
            else
                post_bin_spks = zeros(size(trial_bin_spks));
            end
    
            % Bin the long-stim stuff
            if ~isempty(long_stim_on)
                this_on_events = long_stim_on;
                edges = -long_win:long_dt:long_win;
                bin_spks = zeros(length(this_on_events), length(edges)-1);
                for iT = 1:length(this_on_events)
                    temp = restrict(this_cell, iv(this_on_events(iT) - long_win, ...
                        this_on_events(iT) + long_win));
                    temp = temp.t{1} - this_on_events(iT); % center it around 0
                    [count, ~, ~] = histcounts(temp, edges);
                    bin_spks(iT,:) = count;
                end
                long_bin_spks = mean(bin_spks);
                long_bin_spks = (long_bin_spks - min(long_bin_spks))/...
                    (max(long_bin_spks) - min(long_bin_spks));
            else
                long_bin_spks = zeros(size(trial_bin_spks));
            end
            
            s_out.pre_stim = [s_out.pre_stim; pre_bin_spks];
            s_out.trial_stim = [s_out.trial_stim; trial_bin_spks];
            s_out.post_stim = [s_out.post_stim; post_bin_spks];
            s_out.long_stim = [s_out.long_stim; long_bin_spks];
    
            s_out.labels = [s_out.labels; string(fn_prefix)];
            s_out.depth = [s_out.depth; ExpKeys.probeDepth];
            s_out.stim_mode = [s_out.stim_mode; string(ExpKeys.stim_mode)];
            s_out.short_sw = [s_out.short_sw; ExpKeys.short_stim_pulse_width];
            s_out.long_sw = [s_out.long_sw; ExpKeys.long_stim_pulse_width];
    
            % Load the stim_responses
            load('stim_phases.mat');
            s_out.ntrials = [s_out.ntrials; goodTrials(2) + 1 - goodTrials(1)];  
        end
    end
end

