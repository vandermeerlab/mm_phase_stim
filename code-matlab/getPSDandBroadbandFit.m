%% Assumes that good LFPs have been picked out

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

function doStuff
    % Setting up parameters
    LoadExpKeys;
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

%     % Downsample CSC if Sampling Frequency >2.6 kHz for faster computation
%         Fs = 1/median(diff(csc.tvec));
%     if  Fs > 3000
%         csc.data = decimate(csc.data,12);
%         csc.tvec = csc.tvec(1:12:end);
%         csc.cfg.hdr{1}.SamplingFrequency = csc.cfg.hdr{1}.SamplingFrequency/12;
%     end

    % Find gaps > 2 msec
    gaps = find(diff(csc.tvec) > 0.002);

    overall_max = csc.cfg.hdr{1}.InputRange * 1e-6;
    overall_min = -overall_max;
    all_saturated = (csc.data(1,:) == overall_max | csc.data(1,:) == overall_min);
    all_saturated = csc.tvec(all_saturated);

    % Set this_start and this_stop to the largest uninterrupted, unsaturated epoch
    this_start = ExpKeys.stim_times(1);
    this_stop = ExpKeys.stim_times(2);
    this_gaps = gaps(csc.tvec(gaps) <= this_stop & csc.tvec(gaps) >= this_start);
    this_saturated = all_saturated(all_saturated <= this_stop & all_saturated >= this_start)';
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
    
    % Convert this_csc into fieldtrip format for trying out IRASA stuff
    ft_csc = convert_tsd_to_ft(csc);

    cfg = [];
    cfg.length = 2;
    cfg.overlap = 0.5;
    data = ft_redefinetrial(cfg, ft_csc);
    cfg               = [];
    cfg.foi        = 1:1:120;
    cfg.pad           = 'nextpow2';
    cfg.method        = 'irasa';
    cfg.output        = 'original';
    original = ft_freqanalysis(cfg, data);
    cfg.output        = 'fractal';
    fractal = ft_freqanalysis(cfg, data);

    cfg = [];
    cfg.freq_range = [original.freq(1) original.freq(end)];
    cfg.power_line = '60';
    cfg.peak_width_limits = [0.5,12];
    cfg.max_peaks = 4;
    cfg.min_peak_height = 0.3;
    cfg.aperiodic_mode = 'knee'; %Check with 'fixed' first
    cfg.peak_threshold = 2;
    cfg.return_spectrum = 1;
    cfg.border_threshold = 1;
    cfg.peak_type = 'best'; %There is an error in documenation where it says 'both'
    cfg.proximity_threshold = 2;
    cfg.guess_weight = 'none';
    cfg.thresh_after = 1;
    cfg.sort_type = 'param';
    cfg.sort_param = 'frequency';
    cfg.sort_bands = {{'delta'}, {'2', '5'}; {'theta'}, {'6', '10'}; {'beta'},{'12', '28'}
                      {'gamma1'}, {'30',' 55'}};
    [fs, fg] = process_fooof('FOOOF_matlab', reshape(original.powspctrm, 1, 1, length(original.powspctrm)), ...
        original.freq, cfg, 1);
    powspctrm_f = cat(1, fg.ap_fit);
    aperiodic_P = [];
    for k = 1:size(powspctrm_f,1)
        aperiodic_P(k,:) = interp1(fs, powspctrm_f(k,:), original.freq, 'linear', nan);
    end

    % Save variables for easy retrieval later
    psd = [];
    psd.freq = original.freq;
    psd.original = original.powspctrm;
    psd.irasa  = fractal.powspctrm;
    psd.fooof = aperiodic_P;

    save('trial_psd', 'psd');
  
end