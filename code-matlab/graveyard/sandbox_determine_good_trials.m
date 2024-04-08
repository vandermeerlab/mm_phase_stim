%% Script to play around with power_thresholds to determine what are good trials
cd('data\M019\M019-2019-04-14_vStr_4p2_light_cells_TT5 _min')
LoadExpKeys;
if isempty(ExpKeys.goodCell)
        return
end
fbands = {[2 5], [6 10], [12,28], [30 55]};
c_list = {'red', 'blue','magenta','cyan'};

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
% csc = restrict(csc, iv(ExpKeys.stim_times));

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

%%

% set up filter
cfg_filter = [];
cfg_filter.f = fbands{1};
cfg_filter.type = 'fdesign';%'cheby1';
cfg_filter.order = 4; %5;

ncycles = 3;
% Event Detection
cfg_evt = [];
% cfg_evt.epoch = 'all';
cfg_evt.epochLength = 5*60; % was 5 * 60
cfg_evt.filter_cfg = cfg_filter;
cfg_evt.minlen = ncycles./mean(fbands{1}); % or, 0.05
%cfg_evt.minlen = 0.05;
cfg_evt.smooth = 0.05; % convolve with Gaussian of this SD
cfg_evt.threshold = 0.75; %cfg.detect_thr(iFreq);
cfg_evt.method = 'percentile';
[evt,evt_thr] = detectOscillatoryEvents(cfg_evt,csc,ExpKeys);

% Check if any stim is within this limit
goodStim = [];
for iS = 1:length(stim_on)
    if any((stim_on(iS) >= evt.tstart) & (stim_on(iS) <= evt.tend))
        goodStim = [goodStim iS];
    end
end
disp(length(goodStim))

%%

% Plot data and see
figure;
plot(csc.tvec, csc.data, 'black');
hold on;
evt_start = nearest_idx3(evt.tstart, csc.tvec);
evt_end = nearest_idx3(evt.tend, csc.tvec);
for iEvt = 1:length(evt.tstart)
    plot(csc.tvec(evt_start(iEvt):evt_end(iEvt)), ...
        csc.data(evt_start(iEvt):evt_end(iEvt)), 'red');
end
 xline(stim_on, 'blue')



%%
causal_phase = nan(length(fbands), length(stim_on));
causal_power = cell(length(fbands), length(stim_on));
nEnds = nearest_idx3(stim_on, csc.tvec);
for iB = 1:length(fbands)
    win_length  = 1; % sec
    nStarts = nearest_idx3(stim_on - win_length, csc.tvec);
    for iS = 1:length(stim_on)
        this_echt = echt(csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
        this_phase = angle(this_echt);
        causal_phase(iB,iS) = this_phase(end); % The last sample's phase
        causal_power{iB,iS} = abs(this_echt); 
% %             Diagnostic to check how are angles assined (Uncomment to run)
%             diag_fig = figure;
%             ax1 = subplot(2,1,1);
%             plot(abs(this_echt))
%             ax2 = subplot(2,1,2);
%             plot(angle(this_echt))
%             hold on
%             yline(0, 'red')
%             yline(pi/2, 'green')
%             yline(pi, 'black')
%             linkaxes([ax1,ax2],'x')
%             close(diag_fig)
    end
end
% Inverting phases (because LFP was recorded upside down)
causal_phase = -1*causal_phase;

% 