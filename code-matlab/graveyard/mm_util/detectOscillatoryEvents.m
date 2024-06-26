function [evt_iv,evt_thr] = detectOscillatoryEvents(cfg_in,csc,ExpKeys)
% function [evt_iv,evt_thr] = detectOscillatortEvents(cfg_in,csc,ExpKeys) %
% Modified from Julien_DetectEvents
%
% detects LFP events in specified frequency band exceeding specified
% threshold
%
% identifies threshold based on z-scoring within a specific epoch (returned
% as evt_thr)
%
% 
% manishm 2024-06-13

% input processing
cfg_def = [];
cfg_def.verbose = 1;
cfg_def.cscLabel = []; % need to specify this if csc has more than one signal
cfg_def.filter_cfg = []; % input config for FilterLFP()
cfg_def.smooth = []; % if specified, convolve with gaussian of this SD (in SECONDS)
cfg_def.signalType = 'amplitude'; % {'amplitude','power'}
cfg_def.threshold = 3; % threshold for detecting events
cfg_def.method = 'zscore'; % what threshold means
cfg_def.minlen = 0.05; % minimum event length
cfg_def.verbose = 0; 
cfg = ProcessConfig(cfg_def,cfg_in);

if isempty(cfg.filter_cfg)
    error('Must specify cfg_in.filter_cfg!');
end

% select specified channel if multidimensional data
nInputSignals = size(csc.data,1);
if nInputSignals == 1
    this_csc = csc;
else
    if isempty(cfg.cscLabel)
       error('Multiple signals present but cfg.cscLabel is not defined.'); 
    end

    this_csc = TSD_SelectChannel(csc,cfg.cscLabel);
end

% filter
cscF = FilterLFP(cfg.filter_cfg,this_csc); % note FilterLFP handles NaNs in the data, but puts them back!

% replace NaNs temporarily
nan_idx = find(isnan(cscF.data));
cscF.data(nan_idx) = 0;

% extract envelope (or power)
cscF.data = abs(hilbert(cscF.data));

switch cfg.signalType
    case 'amplitude'
        % do nothing
    case 'power'
        cscF.data = cscF.data.^2;
end

% optionally, smooth
if ~isempty(cfg.smooth) % hmm, this should be farmed out to a nice standardized function
    conv_window_len = 1; % seconds
    dt = median(diff(cscF.tvec));
    fprintf('Julien_DetectEvents.m: smoothing with %.2f sec window, dt %.5f...\n',conv_window_len,dt);
    
    gauss_window = conv_window_len./dt; % 1 second window
    gauss_SD = cfg.smooth./dt;
    
    gk = gausskernel(gauss_window,gauss_SD); gk = gk./max(gk);
    cscF.data = conv2(cscF.data,gk,'same');
end

% Use Pre-stim and post-stim baseline to determine z-score
if isempty(ExpKeys.post_baseline_times)
    cscF_r = restrict(csc, iv(ExpKeys.pre_baseline_times));
else
    cscF_r = restrict(csc, MergeIV([], iv([ExpKeys.pre_baseline_times; ExpKeys.post_baseline_times])));
end

% z-score and threshold
cfg_det = [];
cfg_det.operation = '>';
cfg_det.method = cfg.method;
cfg_det.threshold = cfg.threshold;
cfg_det.minlen = cfg.minlen;

[~,evt_thr] = TSDtoIV(cfg_det,cscF_r);

% now apply to all data
cfg_sel = [];
cfg_sel.operation = '>';
cfg_sel.method = 'raw';
cfg_sel.threshold = evt_thr;
cfg_sel.minlen = cfg.minlen;

evt_iv = TSDtoIV(cfg_sel,cscF);