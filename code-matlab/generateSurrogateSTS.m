%% Script to generate fake spike-phase locking data for significance testing

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
    this_csc = restrict(csc, iv([this_start, this_stop]));
    
    % Convert this_csc into fieldtrip format
    ft_csc = convert_tsd_to_ft(this_csc);
    ft_csc.time{1} = ft_csc.time{1}*1e6;
    ft_csc.hdr.FirstTimeStamp = ft_csc.hdr.FirstTimeStamp*1e6;
    ft_csc.hdr.LastTimeStamp = ft_csc.hdr.LastTimeStamp*1e6;
    ft_csc.hdr.TimeStampPerSample = ft_csc.hdr.TimeStampPerSample*1e6;

    %  Generate fake spikes
    cfg = []; cfg.fc = ExpKeys.goodCell(1);
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
    S = LoadSpikes(cfg);
    % Make fake spikes
    S.t{1} = this_csc.tvec(1):spk_dt:this_csc.tvec(end);
    
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
 
    all_stim = [stim_on-0.5,stim_on+0.5];
    clean_iv = MergeIV([], iv(all_stim));
    clean_iv = InvertIV(clean_iv, ExpKeys.recording_times(1), ExpKeys.recording_times(2));
    S = restrict(S, clean_iv);
    if ~strcmp(ExpKeys.experimenter, 'EC')
        S.ft_spikes =  ft_read_spike(S.label{1}, 'encoding', 64);
    else
        S.ft_spikes =  ft_read_spike(S.label{1}, 'encoding', 32);
    end
        
    S.ft_spikes.timestamp{1} = S.t{1}*1e6;
    fake_data = ft_appendspike([],ft_csc, S.ft_spikes);
    cfg=[];
    cfg.begsample = ft_csc.sampleinfo(1);
    cfg.endsample = ft_csc.sampleinfo(2);
    fake_data = ft_redefinetrial(cfg, fake_data);
    fprintf("Number of fake spikes is %d\n", sum(fake_data.trial{1}(2,:)));

    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.foi = 1:1:100;
    cfg.t_ftimwin = 5./cfg.foi;
    cfg.taper = 'hanning';
    cfg.spikechannel =  S.ft_spikes.label{1};
    cfg.channel = fake_data.label{1};
    pool_sts = ft_spiketriggeredspectrum(cfg, fake_data);  
    save('surrogate_sts','pool_sts');
end
