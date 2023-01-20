%% Assumes that spike-sorting has been done

top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1 :length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        rng(491994) % Random seed set for each session
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff
    end

end
%%
function doStuff
    LoadExpKeys;
    evs = LoadEvents([]);
    cfg = [];
    if isempty(ExpKeys.goodCell)
        return
    end
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.min_cluster_quality = 3;
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
%     cfg.fc = ExpKeys.goodCell;
    S = LoadSpikes(cfg);
    
    %% Set variables
    max_delay = 0.01; % sec (window for the first response since stimulus)
    max_long_delay = 0.4; % sec
    num_sham = 1000000; % Number of sham stim used
    num_nonstim = 10000; % Number of spoints to be selected for non_stim
    
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
        this_cell = SelectTS([], S, iC);
        restricted_cell = SelectTS([], restricted_S, iC);
        od = [];
        od.pre_stim = [];
        od.trial_stim = [];
        od.post_stim = [];
        od.long_stim = [];
        od.sham_stim = [];
%         od.trial_nonstim = [];
%         od.matched_nonstim = [];
        
        % Pre-stim response
        if ~isempty(ExpKeys.pre_stim_times)
            this_on_events = pre_stim_on;
            latency = nan(size(this_on_events));
            latency_wo_stim = nan(size(this_on_events));
            fr = zeros(size(this_on_events));
            fr_wo_stim = zeros(size(this_on_events));
            bfr = zeros(size(this_on_events));
            for iStim = 1:length(latency_wo_stim)
                st.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)) & ...
                    (this_cell.t{1} < this_on_events(iStim)+max_delay));
                st_wo_stim.t{1} = restricted_cell.t{1}((restricted_cell.t{1} >= this_on_events(iStim)) & ...
                    (restricted_cell.t{1} < this_on_events(iStim)+max_delay));
                baseline.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)-max_delay) & ...
                    (this_cell.t{1} < this_on_events(iStim)));

                if ~isempty(st.t{1})
                    latency(iStim) = st.t{1}(1) - this_on_events(iStim);
                end
                fr(iStim) = length(st.t{1})/max_delay;
                bfr(iStim) = length(baseline.t{1})/max_delay;
                if ~isempty(st_wo_stim.t{1})
                    latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                    assert(latency_wo_stim(iStim) >= ExpKeys.short_stim_pulse_width, ...
                        sprintf('Latency for pre-trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
                end
                fr_wo_stim(iStim) = length(st_wo_stim.t{1})/max_delay;
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
                st.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)) & ...
                    (this_cell.t{1} < this_on_events(iStim)+max_delay));
                st_wo_stim.t{1} = restricted_cell.t{1}((restricted_cell.t{1} >= this_on_events(iStim)) & ...
                    (restricted_cell.t{1} < this_on_events(iStim)+max_delay));
                baseline.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)-max_delay) & ...
                    (this_cell.t{1} < this_on_events(iStim)));
                if ~isempty(st.t{1})
                    latency(iStim) = st.t{1}(1) - this_on_events(iStim);
                end
                fr(iStim) = length(st.t{1})/max_delay;
                bfr(iStim) = length(baseline.t{1})/max_delay;
                if ~isempty(st_wo_stim.t{1})
                    latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                    assert(latency_wo_stim(iStim) >= ExpKeys.short_stim_pulse_width, ...
                        sprintf('Latency for Trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
                end
                fr_wo_stim(iStim) = length(st_wo_stim.t{1})/max_delay;
            end
            od.trial_stim.latency = latency;
            od.trial_stim.latency_wo_stim = latency_wo_stim;
            od.trial_stim.fr = fr;
            od.trial_stim.fr_wo_stim = fr_wo_stim;
            od.trial_stim.bfr = bfr;
            % Saving mean firing rate to generate poisson spikes later
            temp = restrict(this_cell, iv(ExpKeys.stim_times(1), ExpKeys.stim_times(2)));
            od.trial_stim.mfr = length(temp.t{1})/diff(ExpKeys.stim_times);
        end


%         % Non-stim evoked response (Choose 1 point each from each ISI to
%         % match the characteristics of the stim evoked response
%         if ~isempty(ExpKeys.stim_times)
%             dt = 0.001; % in seconds
%             this_on_events = zeros(size(stim_on));
%             for iStim = 1:length(stim_on)-1
%                 tvec = [stim_on(iStim)+0.25:dt:stim_on(iStim+1)-0.01];
%                 this_on_events(iStim) = tvec(randi(length(tvec)));           
%             end
%             tvec = stim_on(end)+0.25:dt:stim_on(end)+1;
%             this_on_events(end) = tvec(randi(length(tvec)));
%             clear tvec
%             latency = nan(size(this_on_events));
%             fr = zeros(size(this_on_events));
%             bfr = zeros(size(this_on_events));
%             for iStim = 1:length(this_on_events)
%                 st.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)) & ...
%                     (this_cell.t{1} < this_on_events(iStim)+max_delay));
%                 baseline.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)-max_delay) & ...
%                     (this_cell.t{1} < this_on_events(iStim)));
%                 if ~isempty(st.t{1})
%                     latency(iStim) = st.t{1}(1) - this_on_events(iStim);
%                 end
%                 fr(iStim) = length(st.t{1})/max_delay;
%                 bfr(iStim) = length(baseline.t{1})/max_delay;
%             end
%             od.trial_nonstim.latency = latency;
%             od.trial_nonstim.fr = fr;
%             od.trial_nonstim.bfr = bfr;
%             od.trial_nonstim.on_events = this_on_events;
%         end

        % Post-stim response
        if ~isempty(ExpKeys.post_stim_times)
            this_on_events = post_stim_on;
            latency = nan(size(this_on_events));
            latency_wo_stim = nan(size(this_on_events));
            fr = zeros(size(this_on_events));
            fr_wo_stim = zeros(size(this_on_events));
            bfr = zeros(size(this_on_events));
            for iStim = 1:length(latency_wo_stim)
                st.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)) & ...
                    (this_cell.t{1} < this_on_events(iStim)+max_delay));
                st_wo_stim.t{1} = restricted_cell.t{1}((restricted_cell.t{1} >= this_on_events(iStim)) & ...
                    (restricted_cell.t{1} < this_on_events(iStim)+max_delay));
                baseline.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)-max_delay) & ...
                    (this_cell.t{1} < this_on_events(iStim)));
                if ~isempty(st.t{1})
                    latency(iStim) = st.t{1}(1) - this_on_events(iStim);
                end
                fr(iStim) = length(st.t{1})/max_delay;
                bfr(iStim) = length(baseline.t{1})/max_delay;
                if ~isempty(st_wo_stim.t{1})
                    latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                    assert(latency_wo_stim(iStim) >= ExpKeys.short_stim_pulse_width, ...
                        sprintf('Latency for post-trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
                end
                fr_wo_stim(iStim) = length(st_wo_stim.t{1})/max_delay;
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
                st.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)) & ...
                    (this_cell.t{1} < this_on_events(iStim)+max_long_delay));
                st_wo_stim.t{1} = restricted_cell.t{1}((restricted_cell.t{1} >= this_on_events(iStim)) & ...
                    (restricted_cell.t{1} < this_on_events(iStim)+max_long_delay));
                baseline.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)-max_delay) & ...
                    (this_cell.t{1} < this_on_events(iStim)));
                if ~isempty(st.t{1})
                    latency(iStim) = st.t{1}(1) - this_on_events(iStim);
                end
                fr(iStim) = length(st.t{1})/max_long_delay;
                bfr(iStim) = length(baseline.t{1})/max_delay;
                if ~isempty(st_wo_stim.t{1})
                    latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                    assert(latency_wo_stim(iStim) >= ExpKeys.long_stim_pulse_width, ...
                        sprintf('Latency for post-trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
                end
                fr_wo_stim(iStim) = length(st_wo_stim.t{1})/max_long_delay;
            end
            od.long_stim.latency = latency;
            od.long_stim.latency_wo_stim = latency_wo_stim;
            od.long_stim.fr = fr;
            od.long_stim.fr_wo_stim = fr_wo_stim;
            od.long_stim.bfr = bfr;
        end

        % Sham-stim from real data as well as poisson train, response to confirm opto cell
        if ~isempty(ExpKeys.stim_times)
            this_time = ExpKeys.stim_times(1):1/32000:ExpKeys.stim_times(2);
            % Generate num_random sham_stim fo
            this_on_events = sort(randsample(this_time, num_sham));
            fr = zeros(size(this_on_events));
            bfr = zeros(size(this_on_events));
%             p_fr = zeros(size(this_on_events));
%             p_bfr = zeros(size(this_on_events));

            
%             % Generate a poisson spike train in this period using the mfr
%             % generate dfr based on poisson spike train with pre_stim_baseline mfr
%             dt = 1/32000;
%             pspike = dt*od.trial_stim.mfr; % based on actual mean_firing_rate
%             spk_poiss = rand(size(this_time)); % random numbers between 0 and 1
%             spk_poiss_t = this_time(spk_poiss < pspike)'; % generate spike train
%             p_cell = this_cell;
%             p_cell.t{1} = spk_poiss_t;

            for iStim = 1:length(this_on_events)
                st.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)) & ...
                    (this_cell.t{1} < this_on_events(iStim)+max_delay));
                baseline.t{1} = this_cell.t{1}((this_cell.t{1} >= this_on_events(iStim)-max_delay) & ...
                    (this_cell.t{1} < this_on_events(iStim)));
%                 p_st.t{1} = p_cell.t{1}((p_cell.t{1} >= this_on_events(iStim)) & ...
%                     (p_cell.t{1} < this_on_events(iStim)+max_delay));
%                 p_baseline.t{1} = p_cell.t{1}((p_cell.t{1} >= this_on_events(iStim)-max_delay) & ...
%                     (p_cell.t{1} < this_on_events(iStim)));

                fr(iStim) = length(st.t{1})/max_delay;
                bfr(iStim) = length(baseline.t{1})/max_delay;
%                 p_fr(iStim) = length(p_st.t{1})/max_delay;
%                 p_bfr(iStim) = length(p_baseline.t{1})/max_delay;
            end
            od.sham_stim.fr = fr;
            od.sham_stim.bfr = bfr;
%             od.sham_stim.p_fr = p_fr;
%             od.sham_stim.p_bfr = p_bfr;
        end
  
        % Save variables
        fn_prefix = extractBefore(restricted_S.label{iC}, '.t');
        save(strcat(fn_prefix, '_stim_response'), 'od');
        close;
    end
end
