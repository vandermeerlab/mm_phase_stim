function statedep_latency_phase_manish
% restoredefaultpath;

% if isunix
%     addpath(genpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'));
%     addpath('/Users/jericcarmichael/Documents/GitHub/EC_state/Basic_functions');
%     
%     all_fig_dir = '/Volumes/Fenrir/State_dep/all_checks/';
%     all_lat_dir = '/Volumes/Fenrir/State_dep/all_lat/';
%     
% else
%     %     addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
%     addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
%     %addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
%     %     addpath('C:\Users\mvdm\Documents\GitHub\EC_state\Basic_functions');
%     addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\EC_State'))
%     %addpath('D:\My_Documents\GitHub\EC_state\Basic_functions');
%     all_fig_dir = 'G:\State_data\all_checks\';
%     all_lat_dir = 'G:\State_data\all_lat\';good
%     
%     %     %cd('D:\data\EC_state\M14_2018-12-01_vStr_light_min');
%     %     cd('C:\data\state-dep\M14_2018-12-01_vStr_light_min');
%     %     %cd('D:\data\EC_state\M13-2018-12-05_dStr_2p2_light_min');
%     %     cd('C:\data\state-dep\M13-2018-12-05_dStr_2p2_light_min');
%     
% end


% 
% global PARAMS
% 
% all_fig_dir = [PARAMS.inter_dir  'all_checks' filesep];
% all_lat_dir = [PARAMS.inter_dir  'all_lat' filesep];
% 
all_lat_dir = 'latencies/';
all_fig_dir = 'figures/';

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
%%
%expkeys duplicate check
these_files = dir;
for iF = 1:length(these_files)
    if length(these_files(iF).name) <=2
        continue
    else
        if strcmp(these_files(iF).name(1:2), '._')
            delete(these_files(iF).name)
        end
    end
end
%%

mkdir(all_lat_dir); mkdir(all_fig_dir);
%% defaults
font_size = 18;
LoadExpKeys
%% set up some baselines
cfg_def = [];
cfg_def.baseline_cor = 'on';


%% load CSC
cfg = [];
cfg.fc = ExpKeys.goodLFP;
if contains(cfg.fc, '-')
    temp = split(cfg.fc,'-');
    cfg.fc = {cat(2,temp{1},'.ncs')};
    this_csc = LoadCSC(cfg);
    cfg_temp.fc = {cat(2,temp{2},'.ncs')};
    ref = LoadCSC(cfg_temp);
    this_csc.data = this_csc.data - ref.data;
    clear temp ref;
else
    this_csc = LoadCSC(cfg);
end
Fs = 1 ./ median(diff(this_csc.tvec));
% Downsample the LFP
if this_csc.cfg.hdr{1}.SamplingFrequency >2000
    d_fac = this_csc.cfg.hdr{1}.SamplingFrequency/2000;
    this_csc.data = decimate(this_csc.data,d_fac);
    this_csc.tvec = this_csc.tvec(1:d_fac:end);
    this_csc.cfg.hdr{1}.SamplingFrequency = this_csc.cfg.hdr{1}.SamplingFrequency./d_fac;
end

%% load events
evs = LoadEvents([]);
% cfg = [];
% cfg.eventList = {'Starting Recording', 'Stopping Recording'};
% cfg.eventLabel = {'start', 'stop'};
% start_stop = LoadEvents(cfg);

% % find the longest recording
% for ii = 1:length(start_stop.t{1})
%     rec_times(ii) = start_stop.t{2}(ii)-start_stop.t{1}(ii);
% end
% [duration, main_rec_idx] = max(rec_times);
% disp(['Longest Recording interval is ' num2str(duration/60) ' minutes in recording slot number ' num2str(main_rec_idx)])
laser_on = evs.t{strcmp(evs.label,ExpKeys.trial_stim_on)};



%% load spikes
cfg = []; cfg.getTTnumbers = 0;
cfg.uint = '64';
S = LoadSpikes(cfg);

%% select a cell
for iC = 1:length(S.label)
    
    %pick this cell
    cfg_select.verbose = 0; 
    this_S = SelectTS(cfg_select, S, iC);
    cell_id = this_S.label{1}(1:end-2);
    
    % check if this is a 'good cell' and not an artifact or non-light mod
    % cell
    if ~ismember(this_S.label, ExpKeys.goodCell)
    continue
    end
    
%     if strcmpi(this_S.label{1}(end-4:end-2), 'Art')  || ~ismember([ExpKeys.subject_id '_' ExpKeys.date '_' cell_id], PARAMS.Good_cells) % don't bother processing artifacts
%         continue
%     end
    
    
    %% get some LFP phases (filtfilt)
    f_list = {[2 5], [6 10],[25 55], [65 100]};
    f_list_label = {'2 - 5', '6 - 10', '25 - 55', '65 - 100'};
    fstop_list = {[1.5 5.5], [5.5 10.5],[24.5 55.5], [64.5 100.5]};

    nShuf = 100;
    
    for iF = 1:length(f_list) % loop across freqs
        %  filter & hilbertize
                cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
                csc_f = FilterLFP(cfg_filt, this_csc);
        
                csc_f.data_amp = abs(hilbert(csc_f.data));
                
                % get the phase
                csc_f.data = angle(hilbert(csc_f.data));
                
                % get the amplitude
                
        
                % get phase for each laser stim (actual)
                stim_phase_idx = nearest_idx3(laser_on, csc_f.tvec);
                stim_phase = csc_f.data(stim_phase_idx);
                % get the amp for each laser stim
                stim_amp = csc_f.data_amp(stim_phase_idx);
                
                
                % get the stim phase histogram (maybe make this a pie?
                % shuffle?
                figure(1); subplot(5, 3, 2);
                [n_val, bin_x] = hist(stim_phase, 5);
%                 title(sprintf('stim phase histo (%.1f-%.1f Hz)', f_list{iF}(1), f_list{iF}(2)));
%                 text(bin_x(1), median(n_val), num2str(circ_rtest(stim_phase),3))
%                 hist(stim_phase, 5);
%                 subplot(4, 3, iF+3);
%                 hist(stim_phase, 36);
                pie(n_val)
%                 saveas(gcf, [all_fig_dir ExpKeys.subject_id '_' ExpKeys.date '_' cell_id(1:end-3) '_f' num2str(floor(f_list{iF}(1))) '_' num2str(floor(f_list{iF}(2))) '_phase_dist.fig']);
   

        % convert to histogram of spikes relative to laser onset (based on ccf function by MvdM
        cfg_ccf =[];
        cfg_ccf.smooth = 0; % get raw values
        cfg_ccf.binsize = 0.0001;
        cfg_ccf.max_t = 0.015;
        
        xbin_centers = -cfg_ccf.max_t-cfg_ccf.binsize:cfg_ccf.binsize:cfg_ccf.max_t+cfg_ccf.binsize; % first and last bins are to be deleted later
        
        
        % verison with phase bin loops
        n_phases = 5;
        fprintf('\nPhase split into %1d bins\n', n_phases)
        [~, edges, ~] = histcounts(-pi:pi, n_phases, 'BinLimits', [-pi, pi]);
        
        % allocate some space
        all_lat =  NaN(2, length(laser_on));
        all_lat_shuf =  NaN(2, length(laser_on),nShuf);
        all_count = NaN(2,length(laser_on));
        all_count_shuf = NaN(2,length(laser_on), nShuf);
        all_phase_vals = NaN(1,length(laser_on));
        all_resp = NaN(2, length(laser_on));
        all_resp_shuf = NaN(2,length(laser_on), nShuf);
        all_amp =  NaN(2, length(laser_on)); % get the phase at each stim
        
        % track the amplitude at each stim
        all_amp(2, :) = stim_amp; % get the amplitude at each stim for the given frequency

        % get the raw phase value for each stim
        all_phase_vals = stim_phase;
        
        % Get the phase labels
        for iPhase = 1:n_phases
            [~, this_phase_idx] = find(stim_phase > edges(iPhase) & stim_phase < edges(iPhase+1));
            
            all_lat(1,this_phase_idx) = iPhase;
            all_lat_shuf(1,this_phase_idx) = iPhase;
            all_count(1,this_phase_idx) = iPhase;
            all_count_shuf(1,this_phase_idx) = iPhase;
            all_amp(1,this_phase_idx) = iPhase;

        end
        

        
        % get the latencies
        [~, cell_idx] = ismember(this_S.label, ExpKeys.goodCell);
        for iEvt = ExpKeys.goodTrials(cell_idx,1):1:ExpKeys.goodTrials(cell_idx,2)
            
            if this_S.t{1}(nearest_idx3(laser_on(iEvt), this_S.t{1}, 1)) > laser_on(iEvt)
                all_lat(2, iEvt) = this_S.t{1}(nearest_idx3(laser_on(iEvt), this_S.t{1}, 1)) - laser_on(iEvt);
            else
                all_lat(2, iEvt) = NaN;
            end
            if all_lat(2, iEvt) >cfg_ccf.max_t
                all_lat(2, iEvt) = NaN;
            end
        end
        
        
        bin_counts = histc(all_lat(2,:),xbin_centers); %get a hist of the
        [~,peak_idx] = max(bin_counts);
        hist_lim = [(xbin_centers(peak_idx)-0.001) (xbin_centers(peak_idx)+0.001)];
        
        n_drops = zeros(n_phases, length(all_lat));
        
        for iEvt = ExpKeys.goodTrials(cell_idx,1):1 :ExpKeys.goodTrials(cell_idx,2)
            if isnan(all_lat(2,iEvt))
                continue
            else % if there is a spike within the 'pre' window remove that event.
                if all_lat(2,iEvt) > (xbin_centers(peak_idx)+0.001) || all_lat(2,iEvt) < (xbin_centers(peak_idx)-0.001)
                    n_drops(iPhase, iEvt) = 1;
                    all_lat(2,iEvt) = NaN;
                end
                % get the number of spikes that happen within the
                % 'main' spike window as determined by the max in
                % the histogram of spike latencies
                if all_lat(2,iEvt) < (xbin_centers(peak_idx)+0.001) && all_lat(2,iEvt) > (xbin_centers(peak_idx)-0.001)
                    % get the spike count
                    this_event_spikes = restrict(this_S, (xbin_centers(peak_idx)-0.001)+laser_on(iEvt), (xbin_centers(peak_idx)+0.001)+laser_on(iEvt));
                    all_count(2,iEvt) = length(this_event_spikes.t{1});
                    % baseline correct (take the same window before the
                    % event and subtract that from the post-event
                    if strcmp(cfg_def.baseline_cor, 'on') || cfg_def.baseline_cor ==1
                        
                        this_event_spikes_pre = restrict(this_S, laser_on(iEvt)-0.002,laser_on(iEvt));
                        all_base_cout(2,iEvt) = length(this_event_spikes.t{1}) - length(this_event_spikes_pre.t{1});
                    end
                end
            end
        end
        fprintf('\n %.0f# events dropped\n', sum(n_drops(iPhase,:)))
        
        % convert lags from S to ms
        
        all_lat(2,:)  = all_lat(2,:)*1000;
        xbin_centers_ms = xbin_centers*1000;
        
        % get the response 1 or nan
        for iEvt = ExpKeys.goodTrials(cell_idx,1): ExpKeys.goodTrials(cell_idx,2)
            if ~isnan(all_count(2, iEvt)) && all_count(2,iEvt) >=1
                all_resp(2,iEvt) = 1;
            else
                all_resp(2,iEvt) = 0;
            end
        end
        
        
        % get the shuffle
        for iS = nShuf:-1:1
            mix = randperm(length(all_lat(2,:)));
            all_lat_shuf(1,:,iS) = all_lat(1,mix);
            all_count_shuf(1,:,iS) = all_lat(1,mix);
            all_resp_shuf(1,:,iS) = all_lat(1,mix);
            
            all_lat_shuf(2,:,iS) = all_lat(2,mix);
            all_count_shuf(2,:,iS)  = all_count(2,mix);
            all_resp_shuf(2,:,iS) = all_resp(2,mix);
            
        end
        % average over shuffles to get a temporary shuffle array
        for iEvt = 1:length(laser_on)
            temp_lat_shuf(iEvt) = nanmean(all_lat_shuf(2,iEvt,:));
            temp_count_shuf(iEvt) = nanmean(all_count_shuf(2,iEvt,:));
            
            all_resp_shuf_mean(iEvt) = nanmean(all_resp_shuf(2,iEvt,:));
        end
        
        % get the means for each phase bin
        for iPhase = 1:n_phases
            [~, this_phase_idx] = find(all_lat(1,:)==iPhase); % get the index of values for this phase
            
            mean_lat(iPhase) = nanmean(all_lat(2,this_phase_idx));
            std_lat(iPhase) = nanstd(all_lat(2,this_phase_idx));
            mean_count(iPhase) = nanmean(all_count(2,this_phase_idx));
            std_count(iPhase) = nanstd(all_count(2,this_phase_idx));
            mean_response(iPhase) = nanmean(all_resp(2,this_phase_idx));
            
            % get the means for shuffles
            mean_shuf(iPhase) = nanmean(temp_lat_shuf(this_phase_idx));
            std_shuf(iPhase) = nanstd(temp_lat_shuf(this_phase_idx));
        end
        
        % make phase labels
        for ii = unique(all_lat(1,:))
            phase_labels{ii} = ['p' num2str(ii)];
        end
        
        %% make a figure
        c_ord = linspecer(n_phases);
        
        % try to make a colored sine wave to match the phases.
        subplot(5, 3, 3)
        hold on
        x = -pi:pi/(50*n_phases):pi;
        wave_phase = sin(x);
        for iPhase = 1:n_phases
            
            plot(-99+(100*iPhase):(100*iPhase), wave_phase(-99+(100*iPhase):(100*iPhase)), 'color', c_ord(iPhase,:), 'linewidth', 5)
            
        end
        axis off
        title(sprintf('%.1f to %.1f Hz', f_list{iF}(1), f_list{iF}(2)), 'fontsize', font_size)
        
        
        for iPhase = 1:n_phases
            [~, this_phase_idx] = find(all_lat(1,:)==iPhase); % get the index of values for this phase
            
            subplot(5, 3, 3*iPhase - 2)
            %             h = histogram(all_lat(iPhase,:), xbin_centers);
            h = histogram(all_lat(2,this_phase_idx), xbin_centers_ms);
            h.FaceColor = c_ord(iPhase,:);
            h.EdgeColor = c_ord(iPhase,:);
            %             ylim([0 20])
            xlim(hist_lim*1000)
            xlabel('latency (ms)')
            set(gca, 'fontsize', font_size)

            
            subplot(5, 3, [5 6 8 9 11 12 14 15])
            hold on
            he = errorbar(iPhase, mean_lat(iPhase), 1.96*std_lat(iPhase), 'o', 'color', c_ord(iPhase,:),'MarkerFaceColor', c_ord(iPhase,:), 'markersize', 10);
            %shuf for this phase bin
            he = errorbar(iPhase, mean_shuf(iPhase), 1.96*std_shuf(iPhase), 'o', 'color', 'k','MarkerFaceColor', 'k', 'markersize', 10);
            
            set(gca, 'fontsize', font_size);
            %             xlabel('Phase')
            %             xlim([-0.5 n_phases+0.5])
            ylabel('Mean latency (ms)')
            set(gca, 'xtick', 1:n_phases, 'xticklabel', phase_labels)
            %             set(gca, 'xticklabel', '')
        end
        
                    % put in the phase of the LFP at stim as a pie chart
%                     subplot(4, ceil(n_phases/2), [7 10])
        
        %         SetFigure([], gcf)
        %         set(gcf, 'position', [435, 50, 949, 697]);
        %
        ha = axes('Position',[0 0 1 1],'Xlim',[-0.1 1.1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.35, 0.98,cell_id, 'fontsize', font_size)
        saveas(gcf, [all_fig_dir cell_id '_f' num2str(floor(f_list{iF}(1))) '_' num2str(floor(f_list{iF}(2))) '_latency.fig']);
%         saveas_eps([cell_id '_f' num2str(floor(f_list{iF}(1))) '_' num2str(floor(f_list{iF}(2))) '_latency'], all_fig_dir)
        close
        %%
        %         % make a count figure
        %         subplot(6, ceil(n_phases/2), [n_phases+8  6*ceil(n_phases/2)])
        %         hold on
        %         he = errorbar(0, all_count_shuf_mean, 2*all_count_shuf_std, 'o','color', 'k','MarkerFaceColor', 'k', 'markersize', 10);
        %
        %         for iPhase = 1:n_phases
        %             [~, this_phase_idx] = find(all_count(1,:)==iPhase); % get the index of values for this phase
        %
        %             hold on
        %             he = errorbar(iPhase, nanmean(all_count(2,this_phase_idx), 2), 2*nanstd(all_count(2,this_phase_idx), 0,2), 'o', 'color', c_ord(iPhase,:),'MarkerFaceColor', c_ord(iPhase,:), 'markersize', 10);
        %             set(gca, 'fontsize', font_size);
        %             xlabel('Phase')
        %             xlim([-0.5 n_phases+0.5])
        %             ylabel('Spike count')
        %             set(gca, 'xticklabel', ['shuf' phase_labels])
        %         end
        %
        %         SetFigure([], gcf)
        %         set(gcf, 'position', [440    37   930   761]);
        %
        %         ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        %         text(0.35, 0.98,[ExpKeys.subject_id '_' ExpKeys.date '  Cell ' cell_id], 'fontsize', font_size)
        %         saveas(gcf, [all_fig_dir ExpKeys.subject_id '_' ExpKeys.date '_' cell_id(1:end-3) '_f' num2str(floor(f_list{iF}(1))) '_' num2str(floor(f_list{iF}(2))) '_latency.png']);
        %         saveas_eps([ExpKeys.subject_id '_' ExpKeys.date '_' cell_id(1:end-3) '_f' num2str(floor(f_list{iF}(1))) '_' num2str(floor(f_list{iF}(2))) '_latency'], all_fig_dir)
        %         %
        %         close all
        %         %                 if strfind(f_list_label{iF}), '.')
        freq_label = ['f_' strrep(strrep(f_list_label{iF}, ' ', ''), '-', '_')];
        %         %                 else
        %
        %                 end
        out.(strrep(cell_id,'-','_')).(freq_label).latency = all_lat;
        out.(strrep(cell_id,'-','_')).(freq_label).latency_shuf = temp_lat_shuf;
        out.(strrep(cell_id,'-','_')).(freq_label).count = all_count;
        out.(strrep(cell_id,'-','_')).(freq_label).resp = all_resp;
        out.(strrep(cell_id,'-','_')).(freq_label).resp_shuf = all_resp_shuf_mean;
        out.(strrep(cell_id,'-','_')).(freq_label).amp = all_amp;
        out.(strrep(cell_id,'-','_')).(freq_label).phase_vals = all_phase_vals;


        close all
    end % end freq loop
    % put the cell info in here for later
    out.(strrep(cell_id,'-','_')).ExpKeys = ExpKeys;
end % end cell loop
if exist('out', 'var')
    save(strcat(all_lat_dir, ExpKeys.subject_id, '-', ExpKeys.date, '_lat.mat'), 'out', '-v7.3')
end

end % end of function
