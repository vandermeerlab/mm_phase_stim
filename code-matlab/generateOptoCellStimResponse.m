%% Assumes that spike-sorting has been done

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
%%
function doStuff
    LoadExpKeys;
    evs = LoadEvents([]);
    cfg_spk = [];
    cfg_spk.fc = ExpKeys.goodCell;
    if isempty(cfg_spk.fc)
        return
    end
    % cfg_spk.min_cluster_quality = 3;
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg_spk.getRatings = 1;
        cfg_spk.uint = '64';
    end
    S = LoadSpikes(cfg_spk);
    
    %% Set variables
    max_delay = 0.01; % sec (window for the first response since stimulus)
    max_long_delay = 0.25; % sec
    
    %% Remove spikes during short stim times
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

    all_short_stim = [pre_stim_on; stim_on; post_stim_on];

    to_be_removed = [all_short_stim,  all_short_stim + ExpKeys.short_stim_pulse_width + stop_delay];

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
        to_be_removed = [to_be_removed; long_stim_on, long_stim_on + ExpKeys.long_stim_pulse_width + stop_delay];
    end

    clean_iv = InvertIV(iv(to_be_removed), ExpKeys.recording_times(1), ...
        ExpKeys.recording_times(2));
    restricted_S = restrict(S, clean_iv);
    
    
    %% Set variables and parameters
    % snippet for autocorrelation
    for iC = 1:length(restricted_S.label)
        fn_prefix = extractBefore(restricted_S.label{iC}, '.t');
        fn_prefix = strrep(fn_prefix, '_', '-');    
        this_cell = SelectTS([], S, iC);
        restricted_cell = SelectTS([], restricted_S, iC);
        od = [];
        od.pre_stim = [];
        od.trial_stim = [];
        od.post_stim = [];
        od.long_stim = [];
        
        % Pre-stim response
        if ~isempty(ExpKeys.pre_stim_times)
            this_on_events = pre_stim_on;
            latency = nan(size(this_on_events));
            latency_wo_stim = nan(size(this_on_events));
            fr = zeros(size(this_on_events));
            fr_wo_stim = zeros(size(this_on_events));
            bfr = zeros(size(this_on_events));
            for iStim = 1:length(latency_wo_stim)
                st = restrict(this_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
                st_wo_stim = restrict(restricted_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
                baseline = restrict(this_cell, iv(this_on_events(iStim)-max_delay, this_on_events(iStim)));
                if ~isempty(st.t{1})
                    latency(iStim) = st.t{1}(1) - this_on_events(iStim);
                    fr(iStim) = length(st.t{1})/max_delay;
                    bfr(iStim) = length(baseline.t{1})/max_delay;
                end
                if ~isempty(st_wo_stim.t{1})
                    latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                    assert(latency_wo_stim(iStim) >= ExpKeys.short_stim_pulse_width, ...
                        sprintf('Latency for pre-trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
                    fr_wo_stim(iStim) = length(st_wo_stim.t{1})/max_delay;
                end
            end
            od.pre_stim.latency = latency;
            od.pre_stim.latency_wo_stim = latency_wo_stim;
            od.pre_stim.fr = fr;
            od.pre_stim.fr_wo_stim = fr_wo_stim;
            od.pre_stim.bfr = bfr;
        end
       
        % Trial-stim response
        if ~isempty(ExpKeys.stim_times)
            this_on_events = stim_on;
            latency = nan(size(this_on_events));
            latency_wo_stim = nan(size(this_on_events));
            fr = zeros(size(this_on_events));
            fr_wo_stim = zeros(size(this_on_events));
            bfr = zeros(size(this_on_events));
            for iStim = 1:length(latency_wo_stim)
                st = restrict(this_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
                st_wo_stim = restrict(restricted_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
                baseline = restrict(this_cell, iv(this_on_events(iStim)-max_delay, this_on_events(iStim)));
                if ~isempty(st.t{1})
                    latency(iStim) = st.t{1}(1) - this_on_events(iStim);
                    fr(iStim) = length(st.t{1})/max_delay;
                    bfr(iStim) = length(baseline.t{1})/max_delay;
                end
                if ~isempty(st_wo_stim.t{1})
                    latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                    assert(latency_wo_stim(iStim) >= ExpKeys.short_stim_pulse_width, ...
                        sprintf('Latency for Trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
                    fr_wo_stim(iStim) = length(st_wo_stim.t{1})/max_delay;
                end
            end
            od.trial_stim.latency = latency;
            od.trial_stim.latency_wo_stim = latency_wo_stim;
            od.trial_stim.fr = fr;
            od.trial_stim.fr_wo_stim = fr_wo_stim;
            od.trial_stim.bfr = bfr;
        end

        % Post-stim response
        if ~isempty(ExpKeys.post_stim_times)
            this_on_events = post_stim_on;
            latency = nan(size(this_on_events));
            latency_wo_stim = nan(size(this_on_events));
            fr = zeros(size(this_on_events));
            fr_wo_stim = zeros(size(this_on_events));
            bfr = zeros(size(this_on_events));
            for iStim = 1:length(latency_wo_stim)
                st = restrict(this_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
                st_wo_stim = restrict(restricted_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
                baseline = restrict(this_cell, iv(this_on_events(iStim)-max_delay, this_on_events(iStim)));
                if ~isempty(st.t{1})
                    latency(iStim) = st.t{1}(1) - this_on_events(iStim);
                    fr(iStim) = length(st.t{1})/max_delay;
                    bfr(iStim) = length(baseline.t{1})/max_delay;
                end
                if ~isempty(st_wo_stim.t{1})
                    latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                    assert(latency_wo_stim(iStim) >= ExpKeys.short_stim_pulse_width, ...
                        sprintf('Latency for post-trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
                    fr_wo_stim(iStim) = length(st_wo_stim.t{1})/max_delay;
                end
            end
            od.post_stim.latency = latency;
            od.post_stim.latency_wo_stim = latency_wo_stim;
            od.post_stim.fr = fr;
            od.post_stim.fr_wo_stim = fr_wo_stim;
            od.post_stim.bfr = bfr;
        end

        % Long-stim response
        if ~isempty(ExpKeys.long_stim_times)
            this_on_events = long_stim_on;
            latency = nan(size(this_on_events));
            latency_wo_stim = nan(size(this_on_events));
            fr = zeros(size(this_on_events));
            fr_wo_stim = zeros(size(this_on_events));
            bfr = zeros(size(this_on_events));
            for iStim = 1:length(latency_wo_stim)
                st = restrict(this_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_long_delay));
                st_wo_stim = restrict(restricted_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_long_delay));
                baseline = restrict(this_cell, iv(this_on_events(iStim)-max_delay, this_on_events(iStim)));
                if ~isempty(st.t{1})
                    latency(iStim) = st.t{1}(1) - this_on_events(iStim);
                    fr(iStim) = length(st.t{1})/max_long_delay;
                    bfr(iStim) = length(baseline.t{1})/max_delay;
                end
                if ~isempty(st_wo_stim.t{1})
                    latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                    assert(latency_wo_stim(iStim) >= ExpKeys.long_stim_pulse_width, ...
                        sprintf('Latency for post-trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
                    fr_wo_stim(iStim) = length(st_wo_stim.t{1})/max_long_delay;
                end
            end
            od.long_stim.latency = latency;
            od.long_stim.latency_wo_stim = latency_wo_stim;
            od.long_stim.fr = fr;
            od.long_stim.fr_wo_stim = fr_wo_stim;
            od.long_stim.bfr = bfr;
        end       
      
        % Save variables
        fn_prefix = extractBefore(restricted_S.label{iC}, '.t');
        save(strcat(fn_prefix, '_stim_response'), 'od');
        close;
    end
end