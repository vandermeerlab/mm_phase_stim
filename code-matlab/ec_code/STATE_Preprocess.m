function [STATE_data, fnames] = STATE_Preprocess(cfg_in)
%% STATE_Preprocess: this function will load the 'main CSC', all spikes, and any additional experimental data (such as wheel turns)
%
%
%
% Inputs:
%
%
% Outputs: STATE_data: [struct] contains fields: 
%                                                   S - Spike
%                                                   CSC - lfp data
%                                                   Evt - events
%                                                   Dist - distance date from the wheel 
%                                                   Hdr -  header 
%
%
%
% EC 2018-12-17

%% set up defults
 global PARAMS
 
 cfg_def = [];
 cfg_def.csc_chan = {'CSC22.ncs'}; % default csc channel which corresponds to Site # 1 on the NN Buz32
 cfg_def.S_chan = []; % keep empty if you want all the .t files in the current folder.  
 cfg_def.font_size = 18;
 
 cfg = ProcessConfig(cfg_def, cfg_in);

%% get the information from the current folder

if isunix
    fname = strsplit(cd, '/');
else
    fname = strsplit(cd, '/');
end
fname = fname{end}; 
fname = fname(1:strfind(fname,'p')+1);
fname = strrep(fname, '-', '_');

%% make a header 'hdr' with general information taken from the fname.  

hdr = [];
hdr.subject = fname(1:3);
hdr.date = fname(5:14);
hdr.depth = str2double([fname(strfind(fname, 'p')-1) '.' fname(strfind(fname, 'p')-+1)]); % convert the depth from the title p = . in file names used here

if hdr.depth < 2; % if less than 2mm define as cortex
    hdr.target = 'crtx';
elseif hdr.depth > 2 && hdr.depth < 4 % define dorsal striatum
    hdr.target = 'dStr';
elseif hdr.depth > 4 % define as vStr
    hdr.depth = 'vStr';
end
    

% check for proper names
if ~strcmp(hdr.subject(1), 'M') || ~strcmp(hdr.date(1:3), '201') || ~isnumeric(hdr.depth)
    error('This does not use the standard naming convention')
end

%% Load a CSC

cfg_lfp = [];
cfg_lfp.fc = cfg.csc_chan;

CSC = LoadCSC(cfg_lfp);


hdr.csc_Fs = CSC.cfg.hdr{1}.SamplingFrequency; % get the sampling freq. 

if hdr.csc_Fs < 32000
    warning(['CSC Fs is ' num2str(hdr.csc_Fs) '... Expecting 32kHz'])
end


%% Load the Spikes

cfg_S = [];
cfg_S.getTTnumbers = 0;  % workaround for TT number identifier.  
if ~isempty(cfg.S_chan);
cfg_S.fc = cfg.S_chan;
end

% SHOULD ADD SOMETHING FOR ONLY LOADING CERTAIN TYPES: 'OK', 'Good', 'Art'

S = LoadSpikes(cfg_S);


%% Load the Events

cfg_Evt = [];
cfg_Evt.eventList = {'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000).'};
evt_off = LoadEvents(cfg_Evt); % these are the 'off' events corresponding to all the other actual event starts

cfg_Evt = [];
cfg_Evt.eventList = {'TTL Input on AcqSystem1_0 board 0 port 2 value (0x000A).'};
cfg_Evt.eventLabel = {'laser on'};
laser_on = LoadEvents(cfg_Evt);

% get the laser off times
laser_off = [];
laser_off.type = 'ts'; 
laser_off.cfg = laser_on.cfg;
laser_off.label = {'laser off'};
keep_idx = nearest_idx3(laser_on.t{1}, evt_off.t{1},1);
laser_off.t = {evt_off.t{1}(keep_idx)};

cfg_Evt = [];
cfg_Evt.eventList = {'Starting Recording', 'Stopping Recording'};
cfg_Evt.eventLabel = {'start', 'stop'};
Evt.start_stop = LoadEvents(cfg_Evt);

% find the longest recording
for ii = 1:length(Evt.start_stop.t{1})
    rec_times(ii) = Evt.start_stop.t{2}(ii)-Evt.start_stop.t{1}(ii);
end
[duration, main_rec_idx] = max(rec_times);
disp(['Longest Recording interval is ' num2str(duration/60) ' minutes in recording slot number ' num2str(main_rec_idx)])


Evt.laser_on = restrict(laser_on, Evt.start_stop.t{1}(main_rec_idx), Evt.start_stop.t{2}(main_rec_idx)); 
Evt.laser_off = restrict(laser_off, Evt.start_stop.t{1}(main_rec_idx), Evt.start_stop.t{2}(main_rec_idx)); 

% Evt.laser_on = restrict(laser_on, start_stop.t{1}(main_rec_idx+2), start_stop.t{2}(main_rec_idx+2)); 


% check number of pulses. 
if length(Evt.laser_on.t{1}) ~= 1000 && length(Evt.laser_on.t{1}) ~= 600
   warning('Wrong number of laser pulses. %0.2f when it should be 1000 or in early sessions 600',length(laser_on.t{1}))  
end

% get the wide bursts from the post section of the session.  This should be
% 25 x 100ms pulses with a 100ms isi.  
cfg_Evt = [];
cfg_Evt.eventList = {'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0012).'};
cfg_Evt.eventLabel = {'Wide laser on'};
Evt.Wide_laser_on = LoadEvents(cfg_Evt);


Evt.Wide_laser_off = [];
Evt.Wide_laser_off.type = 'ts'; 
Evt.Wide_laser_off.cfg = laser_off.cfg;
Evt.Wide_laser_off.label = {'laser off'};
keep_idx = nearest_idx3(Evt.Wide_laser_on.t{1}, evt_off.t{1},1);
Evt.Wide_laser_off.t = {evt_off.t{1}(keep_idx)};

% get the wide laser off


%% Load the wheel data

csv_id = FindFile('*.csv');

if ~isempty(csv_id)
    csv_data = readtable(csv_id);
    if ~strcmp(csv_id(end-8:end-5), 'all')
        warning('Wheel data might not contain the entire session')
    end
    
    Dist.all.tvec = table2array(csv_data(:,1));
    Dist.all.data = table2array(csv_data(:,3));
    
    % restrict to the different epochs.   DOESN'T work since the time is
    % not consistent between NLX and Wheel. 
%     Wheel_Fs = 1 ./ median(diff(Dist.all.tvec)) % This doesn't really work.
%     Dist.main = restrict(Dist.all
    
    
    
    
else
    warning('no wheel data found')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Get some basic information for a summary plot %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Loop through cells %%%%%

    for iS = 1:length(S.t)
        % restric to specific cell
    S_r = [];
    S_r.type = 'ts';
    S_r.cfg = S.cfg;
    S_r.t = S.t(iS);
    S_r.label = S.label(iS);

%% fname and location
figure(iS)
p1 = subplot(4,4,1);
% h(1) = plot(1:10,nan(1,10), 'color', [1 1 1]);
text(0,0.8, [hdr.subject ' ' strrep(hdr.date, '_', '-') ' '], 'FontSize', cfg.font_size);
text(0, 0.4, [hdr.target ' depth: ' num2str(hdr.depth) 'mm'], 'FontSize', cfg.font_size);
text(0, 0.0, ['Cell: ' S_r.label{1}], 'FontSize', cfg.font_size)

axis off

%% PSD 
% generate a PSD
wSize = 4098;
[Pxx,F] = pwelch(CSC.data, rectwin(wSize), wSize/2, [], hdr.csc_Fs);

% plot it
subplot(4,3,2)
plot(F, 10*log10(Pxx), 'k', 'LineWidth', 2);
    set(gca, 'XLim', [0 120], 'FontSize', cfg.font_size); grid on;
    xlabel('Frequency (Hz)');


    
%% Spike Peth for Wide band
% NEED TO FIX POSSIBLE TIMING ISSUE

cfg_peth = [];
cfg_peth.window = [-.25 2];
cfg_peth.plot = 'off';
cfg_plot.evt_color = [4,172,218]./255;
[outputS, outputT, outputGau] =SpikePETH(cfg_peth, S_r, [Evt.Wide_laser_on.t{1}; Evt.Wide_laser_off.t{1}]);


subplot(4,3,3);
plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
xlabel('peri-event (sec)');
ylabel('Event #');
ylim([1 length(Evt.Wide_laser_on.t{1})])
xlim(cfg_peth.window);
hold on
if size(Evt.Wide_laser_on.t{1},2) > 1
    rectangle('position', [0 1 abs(mode(Evt.Wide_laser_on.t{1}-Evt.Wide_laser_off.t{1})) length(Evt.Wide_laser_on.t{1})], 'facecolor', [cfg_plot.evt_color 0.5], 'edgecolor', [cfg_plot.evt_color 0.5])
else
    rectangle('position', [0 1 0.001  length(Evt.Wide_laser_on.t{1})], 'facecolor', [cfg_plot.evt_color 0.5], 'edgecolor', [cfg_plot.evt_color 0.5])
end


%% inspect stim artifact
this_stim_binned = zeros(size(CSC.tvec));
idx = nearest_idx3(Evt.laser_on.t{1}, CSC.tvec);
this_spk_binned(idx) = 1;
[xc, tvec] = xcorr(this_csc.data, this_spk_binned, 500);
tvec = tvec ./ Fs;

figure;
plot(tvec, xc, 'k', 'LineWidth', 2);

%% Get the Phase-spike relationships

cfg_phase_peth = [];
cfg_phase_peth.mode = 'phase'; % can be 'phase' or 'power'


STATE_state_peth(cfg_phase_peth, S_r, CSC, Evt.laser_on);







%% get the Power-spike relationships




    end % Cell loop 'iS'