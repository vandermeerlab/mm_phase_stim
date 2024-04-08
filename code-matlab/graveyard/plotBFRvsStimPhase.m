%% Script to generate the relationships between spiking and LFP Phase
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.depth] = deal([]);
[summary.bfr, summary.phase] = deal({});
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end
fbands = {[2 5], [6 10], [12 30], [30 55]};
c_list = {'red', 'blue','magenta', 'green'};
dStr_mask = (summary.depth < 3.5)';


%% Plot histogram of corrleation between baseline firing rate and stim_phase
fig = figure('WindowState', 'maximized');

sel = find(dStr_mask);
dStr_corr = zeros(length(sel), 4);
for i = 1:length(sel)
    for j = 1:4
        this_corr = corrcoef(summary.bfr{sel(i)}', summary.phase{sel(i)}(j,:));
        dStr_corr(i,j) = this_corr(1,2);
    end
end
for i = 1:4
    subplot(2,4,i)
    histogram(dStr_corr(:,i),[-0.1:0.01:0.1], 'FaceColor', c_list{i});
    xlabel('Correlation b/w stim-phase and pre-stim firing-rate', 'FontSize', 12);
    title(sprintf('dStr: %d Hz - %d Hz', fbands{i}(1), fbands{i}(2)))
end

sel = find(~dStr_mask);
vStr_corr = zeros(length(sel), 4);
for i = 1:length(sel)
    for j = 1:4
        this_corr = corrcoef(summary.bfr{sel(i)}', summary.phase{sel(i)}(j,:));
        dStr_corr(i,j) = this_corr(1,2);
    end
end
for i = 1:4
    subplot(2,4,4+i)
    histogram(vStr_corr(:,i),[-0.1:0.01:0.1], 'FaceColor', c_list{i});
    xlabel('Correlation b/w stim-phase and pre-stim firing-rate', 'FontSize', 12);
    title(sprintf('vStr: %d Hz - %d Hz', fbands{i}(1), fbands{i}(2)))
end


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
        % Load stim phases
        load('stim_phases.mat');
        % Load the stim responses
        load(strcat(fn_prefix, '_stim_response.mat'));
        
        sham_dfr = od.sham_stim.fr' - od.sham_stim.bfr';
        opto_dfr = od.trial_stim.fr(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2)) - ...
            od.trial_stim.bfr(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2));
        [h,~,~] = kstest2(sham_dfr, opto_dfr, 'Alpha',0.01);
        
        % Proceed only if kstest2 passes 
        if h == 1
            s_out.depth = [s_out.depth; ExpKeys.probeDepth];
            s_out.bfr{size(s_out.bfr,2)+1} = ...
                od.trial_stim.bfr(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2));
            s_out.labels    = [s_out.labels; string(fn_prefix)];
            s_out.phase{size(s_out.phase,2)+1} = ...
                causal_phase(:,ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2));

        end
    end
end
