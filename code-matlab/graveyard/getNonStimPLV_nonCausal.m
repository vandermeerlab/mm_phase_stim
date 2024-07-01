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
    fbands = {[2 5], [6 10], [12 30], [30 55]};
    stim_win = 0.25; % seconds around stim to ignore spikes
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

    % Downsample CSC if Sampling Frequency >2.6 kHz for faster computation
    Fs = 1/median(diff(csc.tvec));
    if  Fs > 3000
        csc.data = decimate(csc.data,12);
        csc.tvec = csc.tvec(1:12:end);
        csc.cfg.hdr{1}.SamplingFrequency = csc.cfg.hdr{1}.SamplingFrequency/12;
    end

    % Filter the CSCs in each of the bands, and obtain the hilbert
    % transform phases in each band
    filt_phase = cell(length(fbands),1);
    for iB = 1:length(fbands)
        cfg_filt.type = 'fdesign'; 
        cfg_filt.f  = fbands{iB};
        filt_lfp = FilterLFP(cfg_filt, csc);
        filt_phase{iB} = hilbert(filt_lfp.data);
        % Normalizing for operations in complex plane
        filt_phase{iB} = filt_phase{iB}./abs(filt_phase{iB});
    end
    

    % Find gaps > 2 msec
    gaps = find(diff(csc.tvec) > 0.002);

    overall_max = csc.cfg.hdr{1}.InputRange * 1e-6;
    overall_min = -overall_max;
    all_saturated = (csc.data(1,:) == overall_max | csc.data(1,:) == overall_min);
    all_saturated = csc.tvec(all_saturated);

    % clean_spikes
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
    
    post_stim_on = evs.t{strcmp(evs.label, ExpKeys.post_trial_stim_on)};
    if ~isempty(post_stim_on) && ~isempty(ExpKeys.post_stim_times)
        post_stim_on = post_stim_on + start_delay;
        post_stim_on = post_stim_on(post_stim_on >= ExpKeys.post_stim_times(1) & ...
                                    post_stim_on <= ExpKeys.post_stim_times(2));
    else
        post_stim_on = [];
    end
    
    if ~isempty(ExpKeys.long_stim_times)
        long_stim_on = evs.t{strcmp(evs.label, ExpKeys.long_stim_on)};
        if ~isempty(long_stim_on)
            long_stim_on = long_stim_on + start_delay;
            long_stim_on = long_stim_on(long_stim_on >= ExpKeys.long_stim_times(1) & ...
                                        long_stim_on <= ExpKeys.long_stim_times(2));
        end
    else
        long_stim_on = [];
    end
    
    all_stim = [pre_stim_on-stim_win, pre_stim_on+stim_win; stim_on-stim_win,stim_on+stim_win; ...
        post_stim_on-stim_win,post_stim_on+stim_win; long_stim_on-stim_win,long_stim_on+stim_win];
    
    clean_iv = MergeIV([], iv(all_stim));
    clean_iv = InvertIV(clean_iv, ExpKeys.recording_times(1), ExpKeys.recording_times(2));
    S = restrict(S, clean_iv);
    
    for iC = 1:length(S.t)
        this_cell = SelectTS([], S, iC);

        % Set all_start and all_stop to the largest uninterrupted, unsaturated epoch
        this_start = ExpKeys.recording_times(1);
        this_stop = ExpKeys.recording_times(2);
        this_gaps = gaps(csc.tvec(gaps) <= this_stop & csc.tvec(gaps) >= this_start);
        this_saturated = all_saturated(all_saturated <= this_stop & all_saturated >= this_start);
        this_bounds = [this_start this_stop];
        for iG = 1:length(this_gaps)
            this_bounds = [this_bounds csc.tvec(this_gaps(iG))];
            if csc.tvec(this_gaps(iG)+1) <= this_stop
                this_bounds = [this_bounds csc.tvec(this_gaps(iG)+1)];
            end
        end
        this_bounds = [this_bounds this_saturated];
        this_bounds = sort(this_bounds);
        [max_len, midx] = max(diff(this_bounds));
        fprintf('Taking the longest unsaturated segment of %.2f sec out of total %.2f sec\n', ...
                        max_len, this_stop-this_start);
        this_start = this_bounds(midx);
        this_stop = this_bounds(midx+1);
       
        this_S = restrict(this_cell, iv([this_start, this_stop]));
        spk_idx = nearest_idx3(this_S.t{1}, csc.tvec);
        all_spk_count = length(spk_idx);
        if all_spk_count > min_spikes
            all_spk_phase = cellfun(@(x) x(spk_idx), filt_phase, 'UniformOutput', false);
            all_unsampled_mean_phase = cellfun(@(x) circ_mean(angle(x), [], 2), all_spk_phase);
            all_unsampled_plv = cellfun(@(x) resultantlength(x), all_spk_phase);
            % Subsample PLV to de-bias
            [temp_mean_phase, temp_plv] = deal(nan(length(fbands),num_subamples));
            for iSub = 1:num_subamples
                this_idx = randperm(all_spk_count);
                temp_mean_phase(:,iSub) = cellfun(@(x) circ_mean(angle(x(this_idx(1:min_spikes))), [], 2), all_spk_phase);
                temp_plv(:,iSub) = cellfun(@(x) resultantlength(x(this_idx(1:min_spikes))), all_spk_phase);    
            end
            all_subsampled_mean_phase = circ_mean(temp_mean_phase, [], 2); 
            all_subsampled_plv = mean(temp_plv, 2);
            clear temp_mean_phase temp_plv

        else % Save all results as null
            all_spk_phase = [];
            all_unsampled_mean_phase = [];
            all_unsampled_plv = [];
            all_subsampled_plv = [];
            all_subsampled_mean_phase = [];
        end
       

        % Set pre_start and pre_stop to the largest uninterrupted, unsaturated epoch
        this_start = ExpKeys.recording_times(1);
        this_stop = ExpKeys.pre_stim_times(2);
        this_gaps = gaps(csc.tvec(gaps) <= this_stop & csc.tvec(gaps) >= this_start);
        this_saturated = all_saturated(all_saturated <= this_stop & all_saturated >= this_start);
        this_bounds = [this_start this_stop];
        for iG = 1:length(this_gaps)
            this_bounds = [this_bounds csc.tvec(this_gaps(iG))];
            if csc.tvec(this_gaps(iG)+1) <= this_stop
                this_bounds = [this_bounds csc.tvec(this_gaps(iG)+1)];
            end
        end
        this_bounds = [this_bounds this_saturated];
        this_bounds = sort(this_bounds);
        [max_len, midx] = max(diff(this_bounds));
        fprintf('Taking the longest unsaturated segment of %.2f sec out of total %.2f sec\n', ...
                        max_len, this_stop-this_start);
        this_start = this_bounds(midx);
        this_stop = this_bounds(midx+1);
       
        this_S = restrict(this_cell, iv([this_start, this_stop]));
        spk_idx = nearest_idx3(this_S.t{1}, csc.tvec);
        pre_spk_count = length(spk_idx);
        if pre_spk_count > min_spikes
            pre_spk_phase = cellfun(@(x) x(spk_idx), filt_phase, 'UniformOutput', false);
            pre_unsampled_mean_phase = cellfun(@(x) circ_mean(angle(x), [], 2), pre_spk_phase);
            pre_unsampled_plv = cellfun(@(x) resultantlength(x), pre_spk_phase);
            % Subsample PLV to de-bias
            [temp_mean_phase, temp_plv] = deal(nan(length(fbands),num_subamples));
            for iSub = 1:num_subamples
                this_idx = randperm(pre_spk_count);
                temp_mean_phase(:,iSub) = cellfun(@(x) circ_mean(angle(x(this_idx(1:min_spikes))), [], 2), pre_spk_phase);
                temp_plv(:,iSub) = cellfun(@(x) resultantlength(x(this_idx(1:min_spikes))), pre_spk_phase);    
            end
            pre_subsampled_mean_phase = circ_mean(temp_mean_phase, [], 2); 
            pre_subsampled_plv = mean(temp_plv, 2);
            clear temp_mean_phase temp_plv
        else % Save subsampled results as null
            pre_spk_phase = [];
            pre_unsampled_plv = [];
            pre_unsampled_mean_phase = [];
            pre_subsampled_plv = [];
            pre_subsampled_mean_phase = [];
        end

        % Set trial_start and trial_stop to the largest uninterrupted, unsaturated epoch
        this_start = ExpKeys.stim_times(1);
        this_stop = ExpKeys.stim_times(2);
        this_gaps = gaps(csc.tvec(gaps) <= this_stop & csc.tvec(gaps) >= this_start);
        this_saturated = all_saturated(all_saturated <= this_stop & all_saturated >= this_start);
        this_bounds = [this_start this_stop];
        for iG = 1:length(this_gaps)
            this_bounds = [this_bounds csc.tvec(this_gaps(iG))];
            if csc.tvec(this_gaps(iG)+1) <= this_stop
                this_bounds = [this_bounds csc.tvec(this_gaps(iG)+1)];
            end
        end
        this_bounds = [this_bounds this_saturated];
        this_bounds = sort(this_bounds);
        [max_len, midx] = max(diff(this_bounds));
        fprintf('Taking the longest unsaturated segment of %.2f sec out of total %.2f sec\n', ...
                        max_len, this_stop-this_start);
        this_start = this_bounds(midx);
        this_stop = this_bounds(midx+1);
       
        this_S = restrict(this_cell, iv([this_start, this_stop]));
        spk_idx = nearest_idx3(this_S.t{1}, csc.tvec);
        trial_spk_count = length(spk_idx);
        if trial_spk_count > min_spikes
            trial_spk_phase = cellfun(@(x) x(spk_idx), filt_phase, 'UniformOutput', false);
            trial_unsampled_mean_phase = cellfun(@(x) circ_mean(angle(x), [], 2), trial_spk_phase);
            trial_unsampled_plv = cellfun(@(x) resultantlength(x), trial_spk_phase);
            % Subsample PLV to de-bias
            [temp_mean_phase, temp_plv] = deal(nan(length(fbands),num_subamples));
            for iSub = 1:num_subamples
                this_idx = randperm(trial_spk_count);
                temp_mean_phase(:,iSub) = cellfun(@(x) circ_mean(angle(x(this_idx(1:min_spikes))), [], 2), trial_spk_phase);
                temp_plv(:,iSub) = cellfun(@(x) resultantlength(x(this_idx(1:min_spikes))), trial_spk_phase);    
            end
            trial_subsampled_mean_phase = circ_mean(temp_mean_phase, [], 2); 
            trial_subsampled_plv = mean(temp_plv, 2);
            clear temp_mean_phase temp_plv
        else % Save subsampled results as null
            trial_spk_phase = [];
            trial_unsampled_plv = [];
            trial_unsampled_mean_phase = [];
            trial_subsampled_plv = [];
            trial_subsampled_mean_phase = [];
        end
        
        if ~isempty(ExpKeys.post_baseline_times)
            % Set post_start and post_stop to the largest uninterrupted, unsaturated epoch
            this_start = ExpKeys.post_baseline_times(1);
            this_stop = ExpKeys.recording_times(2);
            this_gaps = gaps(csc.tvec(gaps) <= this_stop & csc.tvec(gaps) >= this_start);
            this_saturated = all_saturated(all_saturated <= this_stop & all_saturated >= this_start);
            this_bounds = [this_start this_stop];
            for iG = 1:length(this_gaps)
                this_bounds = [this_bounds csc.tvec(this_gaps(iG))];
                if csc.tvec(this_gaps(iG)+1) <= this_stop
                    this_bounds = [this_bounds csc.tvec(this_gaps(iG)+1)];
                end
            end
            this_bounds = [this_bounds this_saturated];
            this_bounds = sort(this_bounds);
            [max_len, midx] = max(diff(this_bounds));
            fprintf('Taking the longest unsaturated segment of %.2f sec out of total %.2f sec\n', ...
                            max_len, this_stop-this_start);
            this_start = this_bounds(midx);
            this_stop = this_bounds(midx+1);
            this_S = restrict(this_cell, iv([this_start, this_stop]));
            spk_idx = nearest_idx3(this_S.t{1}, csc.tvec);
            post_spk_count = length(spk_idx);
            if post_spk_count > min_spikes
                post_spk_phase = cellfun(@(x) x(spk_idx), filt_phase, 'UniformOutput', false);
                post_unsampled_mean_phase = cellfun(@(x) circ_mean(angle(x), [], 2), post_spk_phase);
                post_unsampled_plv = cellfun(@(x) resultantlength(x), post_spk_phase); 
                % Subsample PLV to de-bias
                [temp_mean_phase, temp_plv] = deal(nan(length(fbands),num_subamples));
                for iSub = 1:num_subamples
                    this_idx = randperm(post_spk_count);
                    temp_mean_phase(:,iSub) = cellfun(@(x) circ_mean(angle(x(this_idx(1:min_spikes))), [], 2), post_spk_phase);
                    temp_plv(:,iSub) = cellfun(@(x) resultantlength(x(this_idx(1:min_spikes))), post_spk_phase);    
                end
                post_subsampled_mean_phase = circ_mean(temp_mean_phase, [], 2); 
                post_subsampled_plv = mean(temp_plv, 2);
                clear temp_mean_phase temp_plv
            else % Save subsampled results as null
                post_spk_phase = [];
                post_unsampled_plv = [];
                post_unsampled_mean_phase = [];
                post_subsampled_plv = [];
                post_subsampled_mean_phase = [];
            end
        else
                post_spk_count = 0;
                post_spk_phase = [];
                post_unsampled_plv = [];
                post_unsampled_mean_phase = [];
                post_subsampled_plv = [];
                post_subsampled_mean_phase = [];
        end
        % Save variables
        fn_prefix = extractBefore(this_cell.label{1}, '.t');
        fn_prefix = strrep(fn_prefix, '_', '-');
        save(strcat(fn_prefix, '_spike_phaselock_plv'), ...
            'all_spk_count', 'all_unsampled_plv', 'all_unsampled_mean_phase', 'all_subsampled_plv', 'all_subsampled_mean_phase', ...
            'pre_spk_count', 'pre_unsampled_plv', 'pre_unsampled_mean_phase', 'pre_subsampled_plv', 'pre_subsampled_mean_phase', ...
            'trial_spk_count', 'trial_unsampled_plv', 'trial_unsampled_mean_phase', 'trial_subsampled_plv', 'trial_subsampled_mean_phase', ...
            'post_spk_count', 'post_unsampled_plv', 'post_unsampled_mean_phase', 'post_subsampled_plv', 'post_subsampled_mean_phase', ...
            'all_spk_phase', 'pre_spk_phase', 'trial_spk_phase', 'post_spk_phase');
    end
end

%% Helper functions 
% Borrowed from FieldTrip
function [resLen] = resultantlength(angles)

n = sum(~isnan(angles));
resLen = abs(nansum(angles))./n;

end
