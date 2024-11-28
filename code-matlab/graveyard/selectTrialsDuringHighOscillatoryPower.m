%% Script to select trials that have high oscillatory power in each of the frequency bands

top_dir = 'data\';
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

%
function doStuff
    % Declaring variables
    % Setting up parameters
    fbands = {[2 5], [6 10], [12 28], [30 55]};

    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    evs = LoadEvents([]);

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
    
    % Sometimes csc.tvec can be have weird elements because of gaps in recording
    to_remove = find(diff(csc.tvec)<=0);
    while (~isempty(to_remove))
        csc.tvec(to_remove+1) = [];
        csc.data(to_remove+1) = [];
        to_remove = find(diff(csc.tvec)<=0);
    end
   
    Fs = 1/median(diff(csc.tvec));
    if contains(ExpKeys.light_source, 'LASER')
        start_delay = 0.0011;
    else
        start_delay = 0;
    end

    stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
    if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
        stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                          stim_on <= ExpKeys.stim_times(2));
    else
        stim_on = [];
    end

    pct_thresh = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.9]; % Percentile based thresholding of oscillatory power
    ncycles = 3; 
    
    good_stim = cell(length(fbands), length(pct_thresh));

    for iB = 1:length(fbands)
        cfg_filter = [];
        cfg_filter.f = fbands{iB};
        cfg_filter.type = 'fdesign';
        cfg_filter.order = 4;

        for iT = 1:length (pct_thresh)
            % Event Detection
            cfg_evt = [];
            cfg_evt.filter_cfg = cfg_filter;
            cfg_evt.minlen = ncycles./mean(fbands{iB});
            cfg_evt.smooth = 0.05; % convolve with Gaussian of this SD
            cfg_evt.threshold = pct_thresh(iT);
            cfg_evt.method = 'percentile';
            [evt,evt_thr] = detectOscillatoryEvents(cfg_evt,csc,ExpKeys);

            % Track which stim happened within these oscillatory events
            keep = [];
            for iS = 1:length(stim_on)
                if any((stim_on(iS) >= evt.tstart) & (stim_on(iS) <= evt.tend))
                    keep = [keep iS];
                end
            end
            good_stim{iB,iT} = keep;
        end
    end

    % Assume you are in the correct folder
    save('good_stim','good_stim'); % should add option to save in specified output dir
end