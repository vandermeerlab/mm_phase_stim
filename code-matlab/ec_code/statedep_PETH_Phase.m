function statedep_latency_phase_sandbox
% restoredefaultpath;

if isunix
    addpath(genpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'));
    addpath('/Users/jericcarmichael/Documents/GitHub/EC_state/Basic_functions');
    
    all_fig_dir = '/Volumes/Fenrir/State_dep/all_checks/';
    %     all_lat_dir = '/Volumes/Fenrir/State_dep/all_lat/';
    
else
    %     addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
    %addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
    %     addpath('C:\Users\mvdm\Documents\GitHub\EC_state\Basic_functions');
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\EC_State'))
    %addpath('D:\My_Documents\GitHub\EC_state\Basic_functions');
    all_fig_dir = 'G:\State_data\all_checks\';
    %     all_lat_dir = 'G:\State_data\all_lat\';
    
    %     %cd('D:\data\EC_state\M14_2018-12-01_vStr_light_min');
    %     cd('C:\data\state-dep\M14_2018-12-01_vStr_light_min');
    %     %cd('D:\data\EC_state\M13-2018-12-05_dStr_2p2_light_min');
    %     cd('C:\data\state-dep\M13-2018-12-05_dStr_2p2_light_min');
    
end

global PARAMS


mkdir(all_fig_dir);
%% defaults
font_size = 18;
LoadExpKeys
%% load CSC
cfg = [];
%cfg.decimateByFactor = 30;
cfg.fc = {ExpKeys.goodCSC};
this_csc = LoadCSC(cfg);
Fs = 1 ./ median(diff(this_csc.tvec));

%% load events
cfg = [];
cfg.eventList = {ExpKeys.laser_on};
cfg.eventLabel = {'laser on'};
laser_on = LoadEvents(cfg);

cfg = [];
cfg.eventList = {'Starting Recording', 'Stopping Recording'};
cfg.eventLabel = {'start', 'stop'};
start_stop = LoadEvents(cfg);

% find the longest recording
for ii = 1:length(start_stop.t{1})
    rec_times(ii) = start_stop.t{2}(ii)-start_stop.t{1}(ii);
end
[duration, main_rec_idx] = max(rec_times);
disp(['Longest Recording interval is ' num2str(duration/60) ' minutes in recording slot number ' num2str(main_rec_idx)])


laser_on = restrict(laser_on, start_stop.t{1}(main_rec_idx), start_stop.t{2}(main_rec_idx));

% check number of pulses.
if length(laser_on.t{1}) ~= 1000 && length(laser_on.t{1}) ~= 600
    warning('Wrong number of laser pulses. %0.2f when it should be 1000 or in early sessions 600',length(laser_on.t{1}))
    
end

%% load spikes
cfg = []; cfg.getTTnumbers = 0;
if isfield(ExpKeys, 'goodCell')
    cfg.fc = ExpKeys.goodCell;
end
S = LoadSpikes(cfg);

for iC = 1:length(S.label)
    
    
    %pick this cell
    this_S = SelectTS([], S, iC);
    cell_id = this_S.label{1}(1:end-2);
    cell_id = strrep(cell_id, '_SS', '');
    
    if strcmpi(this_S.label{1}(end-4:end-2), 'Art')  %|| ~ismember([ExpKeys.subject '_' ExpKeys.date '_' cell_id], good_sess_list) % don't bother processing artifacts
        continue
    end
    %% get some LFP phases (filtfilt)
    f_list = {[3 5], [7 10],[15 25], [30 40],[40 60], [60 80]};
    f_list_label = {'3 - 5', '7 - 10', '15 - 25', '30 - 40', '40 - 60', '60 - 80'};

    
    for iF = 1:length(f_list) % loop across freqs
        %%  filter & hilbertize
        cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
        csc_f = FilterLFP(cfg_filt, this_csc);
        
        % get the phase
        csc_f.data = angle(hilbert(csc_f.data));
        
        % get phase for each laser stim (actual)
        stim_phase_idx = nearest_idx3(laser_on.t{1}, csc_f.tvec);
        stim_phase = csc_f.data(stim_phase_idx);
        
        
        % convert to histogram of spikes relative to laser onset (based on ccf function by MvdM
        cfg_ccf =[];
        cfg_ccf.smooth = 0; % get raw values
        cfg_ccf.binsize = 0.0001;
        cfg_ccf.max_t = 0.015;
        %         [ccf_raw, tvec] = ccf(cfg_ccf, laser_on.t{1}, this_S.t{1});
        
        xbin_centers = -cfg_ccf.max_t-cfg_ccf.binsize:cfg_ccf.binsize:cfg_ccf.max_t+cfg_ccf.binsize; % first and last bins are to be deleted later
        
        
        % verison with phase bin loops
        n_phases = 5; 
        c_ord = linspecer(n_phases); % get colors
        fprintf('\nPhase split into %1d bins\n', n_phases)
        [~, edges, ~] = histcounts(-pi:pi, n_phases, 'BinLimits', [-pi, pi]);
        
        %% make a PETH with color coding
        subplot(ceil(length(f_list)/2),2, iF)        % Get the phase labels
                    evt_count = 0;

        for iPhase = 1:n_phases
            [~, this_phase_idx] = find(stim_phase > edges(iPhase) & stim_phase < edges(iPhase+1));

            hold on
            % plot(linspace(0,0.001),-1*ones(1,length(linspace(0,0.001))))
            for iEvt = this_phase_idx
                evt_count = evt_count+1;
                this_event_spikes = restrict(this_S, laser_on.t{1}(iEvt), laser_on.t{1}(iEvt)+0.025);
                
                if ~isempty(this_event_spikes.t{1}) && length(this_event_spikes.t{1}) >2
                    plot([(this_event_spikes.t{1}-laser_on.t{1}(iEvt))*1000 (this_event_spikes.t{1}-laser_on.t{1}(iEvt))*1000],[iEvt-0.5 iEvt+0.5],'Color',c_ord(iPhase, :), 'linewidth', 3);
                elseif ~isempty(this_event_spikes.t{1}) && length(this_event_spikes.t{1}) ==2 % for odd behaviour that plots two values on a diagonal
                    for ii = 1:2
                        plot([(this_event_spikes.t{1}(ii)-laser_on.t{1}(iEvt))*1000 (this_event_spikes.t{1}(ii)-laser_on.t{1}(iEvt))*1000],[iEvt-0.5 iEvt+0.5],'Color',c_ord(iPhase, :), 'linewidth', 3);
                    end
                elseif ~isempty(this_event_spikes.t{1}) && length(this_event_spikes.t{1}) ==1
                    plot([(this_event_spikes.t{1}(1)-laser_on.t{1}(iEvt))*1000 (this_event_spikes.t{1}(1)-laser_on.t{1}(iEvt))*1000],[iEvt-0.5 iEvt+0.5],'Color',c_ord(iPhase, :), 'linewidth', 3);
                else
                    disp([num2str(iEvt) 'n spikes ' num2str(length(this_event_spikes.t{1}))])
                end


            end
            xlim([-1 10])
            ylim([1 length(laser_on.t{1})])
            set(gca, 'xtick',  -1:1:10)
            xlabel('time (ms)')
            ylabel('stim #')
        end
        rectangle('position', [0, 0, 1, length(laser_on.t{1})], 'facecolor',[([4,172,218]./255) 0.5], 'edgecolor',[([4,172,218]./255) 0.5] )
    end
    SetFigure([], gcf)
    text(0.35, 0.98,[ExpKeys.subject '_' ExpKeys.date '  Cell ' cell_id], 'fontsize', font_size)
    saveas(gcf, [all_fig_dir ExpKeys.subject '_' ExpKeys.date '_' cell_id(1:end-3) '_peth.png']);
    saveas_eps([ExpKeys.subject '_' ExpKeys.date '_' cell_id(1:end-3) '_peth'], all_fig_dir)
close all
end
