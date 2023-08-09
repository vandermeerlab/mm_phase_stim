%% Script to generate fake spike-phase locking data for significance testing (PLV version)

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
    spk_dt = 0.0125; % interspike interval for fake spikes to be generated for this session
 
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
    
    % Find gaps > 2 msec
    gaps = find(diff(csc.tvec) > 0.002);

    overall_max = csc.cfg.hdr{1}.InputRange * 1e-6;
    overall_min = -overall_max;
    all_saturated = (csc.data(1,:) == overall_max | csc.data(1,:) == overall_min);
    all_saturated = csc.tvec(all_saturated);

    % code to separate out trial-stim duration
    % Set this_start and this_stop to the largest uninterrupted, unsaturated epoch
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
    csc = restrict(csc, iv([this_start, this_stop]));

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
    
    %  Generate fake spikes
    cfg = []; cfg.fc = ExpKeys.goodCell(1);
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
    S = LoadSpikes(cfg);
    % Make fake spikes
    S.t{1} = csc.tvec(1):spk_dt:csc.tvec(end);
    
    % Get rid of spikes around the stim_times
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
 
    all_stim = [stim_on-stim_win,stim_on+stim_win];
    clean_iv = MergeIV([], iv(all_stim));
    clean_iv = InvertIV(clean_iv, ExpKeys.recording_times(1), ExpKeys.recording_times(2));
    S = restrict(S, clean_iv);
    spk_idx = nearest_idx3(S.t{1}, csc.tvec);
    all_spk_count = length(spk_idx);
    fprintf("Total number of fake spikes is %d\n", all_spk_count);
    fake_spk_phase = cellfun(@(x) x(spk_idx), filt_phase, 'UniformOutput', false);

    save('surrogate_plv','fake_spk_phase');
end
