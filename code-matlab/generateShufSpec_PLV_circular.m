%% Script to generate circularly shifted shuffles for signifcance testing of spike-phase locking (PLV version)
% Assumes that spike_phaselock_plv.mat exists in each session's folder
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff;
    end
end
%%
function doStuff
    rng(491994); % Setting the seed for reproducibility (needs to be set for each session at the start)
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    evs  = LoadEvents([]);
    nshufs = 1000;
    num_subsamples = 1000;
    min_spikes = 200;
    stim_win = 0.25; % seconds
    bin_width = 0.01; % seconds
    
    cfg = [];
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.min_cluster_quality = 3;
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
    cfg.fc = ExpKeys.goodCell;
    S = LoadSpikes(cfg);

    cfg = [];
    cfg.fc = ExpKeys.goodLFP;
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
    Fs = 1/median(diff(csc.tvec));

    for iC = 1:length(ExpKeys.goodCell)
        % Load the spikes
        this_cell = SelectTS([], S, iC);
        this_cell = restrict(this_cell, iv(ExpKeys.stim_times));

        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
%         fn_prefix = strrep(fn_prefix, '_', '-');
        % Load the spike_phaselock data
        load(strcat(fn_prefix, '_nonstim_spk_phases'),'phase_out');
        spk_train = false(size(phase_out.all_phase{1}));
        
        % Get rid of all spikes that happen within 0.25 around stim
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
        
        last_idx = 0;
        for iStim = 2:length(stim_on)
            end_t = stim_on(iStim) - 0.01;
            start_t = stim_on(iStim-1) + stim_win;
            if end_t - start_t > 0.5
                this_spks = this_cell.t{1}((this_cell.t{1} >= start_t) & ...
                    (this_cell.t{1} < end_t));
                end_idx = nearest_idx3(end_t, csc.tvec);
                start_idx = nearest_idx3(start_t, csc.tvec);
                actual_t = csc.tvec(start_idx:end_idx);
                real_idx = nearest_idx3(this_spks, actual_t);
                spk_train(last_idx+real_idx) = true;
                last_idx = last_idx + length(actual_t);
            end
        end
        assert(last_idx == length(spk_train), "Some logical error");
        assert(sum(spk_train) == length(phase_out.all_spk_phase{1}), "Some other logical error");
        shuf_circ_plv = nan(length(phase_out.all_spk_phase),nshufs);
        trial_spk_count = phase_out.spk_count;
        if trial_spk_count >  min_spikes
            to_shift = randperm(last_idx);
            to_shift = to_shift(1:nshufs);
            for iShuf = 1:nshufs
                % circularly shift spikes
                this_shuf_train = circshift(spk_train, to_shift(iShuf));
                this_shuf_phase = cellfun(@(x) x(this_shuf_train), phase_out.all_phase, 'UniformOutput', false);
                % Subsample PLV to de-bias
                temp_plv = nan(length(this_shuf_phase),num_subsamples);
                for iSub = 1:num_subsamples
                    this_idx = randperm(trial_spk_count);
                    temp_plv(:,iSub) = cellfun(@(x) resultantlength(x(this_idx(1:min_spikes))), this_shuf_phase);
                end
                shuf_circ_plv(:,iShuf) = (mean(temp_plv,2))';
            end
        end
        save(strcat(fn_prefix, '_shuf_spec_circ_plv_reworked'), 'shuf_circ_plv');
    end
end

%% Helper functions 
% Borrowed from FieldTrip
function [resLen] = resultantlength(angles)
    n = sum(~isnan(angles));
    resLen = abs(nansum(angles))./n;
end
