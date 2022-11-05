% STATE_Evt_Spec



%% load data
cd('/Volumes/Fenrir/State_dep/M13-2018-12-05_dStr_2p2_light_min')

% get csc
cfg_csc = [];
cfg_csc.fc = {'CSC22.ncs'};
cfg_csc.decimateByFactor = 15;

csc = LoadCSC(cfg_csc);

% get spikes
cfg_s = [];
cfg_s.getTTnumbers = 0;
S = LoadSpikes(cfg_s);

%% again with FT

addpath(PARAMS.ft_dir)
ft_defaults

fc = {'CSC22.ncs'};
data = ft_read_neuralynx_interp(fc);


%% get and append spikes
t_id = FindFile_str(cd, '.t');
% hack around above finding '.txt files
keep_idx = [];
for ii =1:length(t_id)
    if strcmp(t_id{ii}(end-1:end), '.t')
        keep_idx(ii) = 1;
    else
        keep_idx(ii) = 0;
    end
end
t_id(keep_idx==0) = [];
t_id

S_list = {};
for iS = 1:length(t_id)
    spike = ft_read_spike(t_id{iS}); % needs fixed read_mclust_t.m
    spike.label{1} = t_id{iS}(1:end-2);
    disp([spike.label{1} ' Contained: ' num2str(length(spike.timestamp{1})) ' spikes'])

    S_list{end+1} = spike.label{1};
    data = ft_appendspike([],data, spike);
end



%% get events for the main recording only

evt = LoadEvents([]);
for ii = 1:length(evt.label)
    if strcmp(evt.label{ii}(end-4:end), '00A).'); % this is the NLX event for the laser on TS
        laser_on_idx = ii;
    end
end

% find the longest recording
for ii = 1:length(evt.t{1})
    rec_times(ii) = evt.t{2}(ii)-evt.t{1}(ii);
end
[duration, main_rec_idx] = max(rec_times);
disp(['Longest Recording interval is ' num2str(duration/60) 'minutes in recording slot number ' num2str(main_rec_idx)])

laser_on_keep = evt.t{6}(evt.t{1}(main_rec_idx)<evt.t{6} & evt.t{6}<evt.t{2}(main_rec_idx));

% make an events trial set 
cfg = [];
cfg.t = laser_on_keep(1:100);
cfg.mode = 'nlx';
cfg.hdr = data.hdr;
cfg.twin = [-0.5 1];
trl = ft_maketrl(cfg);

% redefine the data to laser pulses.  
cfg = [];
cfg.trl = trl;
data_trl = ft_redefinetrial(cfg, data)

%% Make some Event locked spectrograms
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = 'CSC22';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies of interest
cfg.t_ftimwin    = ones(size(cfg.foi)).*0.5;  % window size: fixed at 0.5s
cfg.toi          = -0.5:0.05:1; % times of interest
 
TFR = ft_freqanalysis(cfg, data_trl);
 
figure
cfg = []; cfg.channel = 'CSC22';
cfg.trials =[1:10];
ft_singleplotTFR(cfg, TFR);
caxis([0 30])
