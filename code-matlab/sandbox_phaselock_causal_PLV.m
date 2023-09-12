%% Script to generate spike-LFP phase locking
rng(491994); % Setting up seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess)); 
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff
    end
end

function doStuff
    % Declare parameters and variables
    fbands = {[2 5], [6 10], [30 55]};
    wsz = 0.5; % seconds, window to use for causal phase estimation
    min_spikes = 200;
    num_subamples = 1000;
 
    LoadExpKeys;
    evs = LoadEvents([]);
    if isempty(ExpKeys.goodCell)
        return
    end
    
    cfg = []; cfg.fc = ExpKeys.goodLFP;
    if contains(cfg.fc, '-')
        temp = split(cfg.fc,'-');
        cfg.fc = {cat(2,temp{1},'.ncs')};
        csc = LoadCSC(cfg);
        cfg_temp.fc = {cat(2,temp{2},'.ncs')};
        ref = LoadCSC(cfg_temp);
        csc.data = csc.data - ref.data;
        clear temp ref;
    else
        csc = LoadCSC(cfg);
    end
    % restrict csc to stim_times
    csc = restrict(csc, iv(ExpKeys.stim_times));
    Fs = 1/median(diff(csc.tvec));

    cfg = []; cfg.fc = ExpKeys.goodCell;
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
    S = LoadSpikes(cfg);
    
    % Remove spikes during short stim times
    if contains(ExpKeys.light_source, 'LASER')
        start_delay = 0.0011;
        stop_delay = 0.0012; %0.0022 is the spike width in the case of 1 msec laser pulse
    else
        start_delay = 0;
        stop_delay = 0;
    end
    
    stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
    if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
        stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                          stim_on <= ExpKeys.stim_times(2));
    else
        stim_on = [];
    end

    all_stim = [stim_on-0.01, stim_on+0.25];
    clean_iv = MergeIV([], iv(all_stim));
    clean_iv = InvertIV(clean_iv, ExpKeys.stim_times(1), ExpKeys.stim_times(2));
    S = restrict(S, clean_iv);
    
    for iC = 1:length(S.t)
        this_cell = SelectTS([], S, iC);
        all_ends = this_cell.t{1};
        % Get rid of spikes that are too close to the border
        all_ends = all_ends(all_ends - csc.tvec(1) >= wsz);
        [trial_unsampled_plv, trial_subsampled_plv, trial_spk_phase, ...
            trial_unsampled_mean_phase, trial_subsampled_mean_phase] = deal([]); 
        if length(all_ends) >= min_spikes
            [trial_unsampled_mean_phase, trial_unsampled_plv, ...
                trial_subsampled_mean_phase, trial_subsampled_plv] = deal(nan(size(fbands)));
            trial_spk_phase = nan(length(fbands), length(all_ends));
            all_starts = nearest_idx3(all_ends - wsz, csc.tvec);
            all_ends = nearest_idx3(all_ends, csc.tvec);
            for iF = 1:length(fbands)
                for iS = 1:length(all_ends)
                    this_echt = echt(csc.data(all_starts(iS):all_ends(iS)), ...
                    fbands{iF}(1), fbands{iF}(2), Fs);
                    this_phase = angle(this_echt);
                    trial_spk_phase(iF, iS) = this_phase(end);
                end
                trial_unsampled_plv(iF) = resultantlength(trial_spk_phase(iF,:));
                trial_unsampled_mean_phase(iF) = circ_mean(trial_spk_phase(iF,:), [], 2);
                [temp_mean_phase, temp_plv] = deal(nan(1,num_subamples));
                for iSub = 1:num_subamples
                    this_idx = randperm(length(all_ends));
                    temp_mean_phase(iSub) = circ_mean(trial_spk_phase(iF,this_idx(1:200)), [], 2);
                    temp_plv(iSub) = resultantlength(trial_spk_phase(iF,this_idx(1:200)));   
                end
                trial_subsampled_mean_phase(iF) = circ_mean(temp_mean_phase, [], 2); 
                trial_subsampled_plv(iF) = mean(temp_plv);
                clear temp_plv temp_mean_phase
            end
        end
        % Save variables
        fn_prefix = extractBefore(this_cell.label{1}, '.t');
        fn_prefix = strrep(fn_prefix, '_', '-');
        save(strcat(fn_prefix, '_spike_phaselock_causal_plv'), ...
            'trial_unsampled_plv', 'trial_unsampled_mean_phase', ...
            'trial_subsampled_plv', 'trial_subsampled_mean_phase', ...
            'trial_spk_phase');
    end
end

%% Helper functions 
% Borrowed from FieldTrip
function [resLen] = resultantlength(angles)

n = sum(~isnan(angles));
resLen = abs(nansum(angles))./n;

end
