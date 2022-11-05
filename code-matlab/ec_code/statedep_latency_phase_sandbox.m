function statedep_latency_phase_sandbox
% restoredefaultpath;

if isunix
    addpath(genpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'));
    addpath('/Users/jericcarmichael/Documents/GitHub/EC_state/Basic_functions');
    
    all_fig_dir = '/Volumes/Fenrir/State_dep/all_checks/';
    all_lat_dir = '/Volumes/Fenrir/State_dep/all_lat/';
    
else
    %     addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
    %addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
    %     addpath('C:\Users\mvdm\Documents\GitHub\EC_state\Basic_functions');
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\EC_State'))
    %addpath('D:\My_Documents\GitHub\EC_state\Basic_functions');
    all_fig_dir = 'G:\State_data\all_checks\';
    all_lat_dir = 'G:\State_data\all_lat\';
    
    %     %cd('D:\data\EC_state\M14_2018-12-01_vStr_light_min');
    %     cd('C:\data\state-dep\M14_2018-12-01_vStr_light_min');
    %     %cd('D:\data\EC_state\M13-2018-12-05_dStr_2p2_light_min');
    %     cd('C:\data\state-dep\M13-2018-12-05_dStr_2p2_light_min');
    
end

mkdir(all_lat_dir); mkdir(all_fig_dir);
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
S = LoadSpikes(cfg);

%% select a cell
for iC = 1:length(S.label)
    
    %pick this cell
    this_S = SelectTS([], S, iC);
    cell_id = this_S.label{1}(1:end-2);
    cell_id = strrep(cell_id, '_SS', '');
    
    if strcmpi(this_S.label{1}(end-4:end-2), 'Art') % don't bother processing artifacts
        continue
    end
    %% get some LFP phases (filtfilt)
    f_list = {[3 5], [6.5 9.5],[15 25], [30 40],[40 60], [60 80]};
    f_list_label = {'3 - 5', '7 - 10', '15 - 25', '30 - 40', '40 - 60', '60 - 80'};
    nShuf = 100;
    
    for iF = 1:length(f_list) % loop across freqs
        %%  filter & hilbertize
        cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
        csc_f = FilterLFP(cfg_filt, this_csc);
        
        % get the phase
        csc_f.data = angle(hilbert(csc_f.data));
        
        % get phase for each laser stim (actual)
        stim_phase_idx = nearest_idx3(laser_on.t{1}, csc_f.tvec);
        stim_phase = csc_f.data(stim_phase_idx);
        
        %same but for shuffle
        for iS = nShuf:-1:1
            
            csc_shuf = csc_f;
            if strcmp(version('-release'), '2014b') % version differences in how circshift is handled.
                csc_shuf.data = circshift(csc_shuf.data, round(rand(1) .* 0.5*length(csc_shuf.data)), 2);
            else
                csc_shuf.data = circshift(csc_shuf.data, round(rand(1) .* 0.5*length(csc_shuf.data)));
            end
            stim_phase_idx = nearest_idx3(laser_on.t{1}, csc_shuf.tvec);
            stim_phase_shuf(iS, :) = csc_shuf.data(stim_phase_idx);
        end % of shuffles
        
        
        %% convert to histogram of spikes relative to laser onset (based on ccf function by MvdM
        cfg_ccf =[];
        cfg_ccf.smooth = 0; % get raw values
        cfg_ccf.binsize = 0.0001;
        cfg_ccf.max_t = 0.015;
        %         [ccf_raw, tvec] = ccf(cfg_ccf, laser_on.t{1}, this_S.t{1});
        
        xbin_centers = -cfg_ccf.max_t-cfg_ccf.binsize:cfg_ccf.binsize:cfg_ccf.max_t+cfg_ccf.binsize; % first and last bins are to be deleted later
        this_ccf = zeros(size(xbin_centers));
        
        for iEvt = 1:length(laser_on.t{1})
            
            relative_spk_t = this_S.t{1} - laser_on.t{1}(iEvt);
            
            this_ccf = this_ccf + hist(relative_spk_t,xbin_centers); % note that histc puts all spikes outside the bin centers in the first and last bins! delete later.
            
            % get the latancy
            all_lat(iEvt) = this_S.t{1}(nearest_idx3(laser_on.t{1}(iEvt), this_S.t{1}, 1)) - laser_on.t{1}(iEvt);
            
        end
        
        this_ccf = this_ccf(2:end-1); % correct for everything being lumped into the first and last bins
        tvec = xbin_centers(2:end-1); % remove unwanted bins
        
        %         zero_idx = find(tvec == 0);
        %         % align to zero
        %         this_ccf = this_ccf(zero_idx:end);
        %         tvec = tvec(zero_idx:end);
        
        
        
        %% verison with phase bin loops
        n_phases = 5;
        fprintf('\nPhase split into %1d bins\n', n_phases)
        [~, edges, ~] = histcounts(-pi:pi, n_phases, 'BinLimits', [-pi, pi]);
        %         all_lat = NaN(n_phases,length(laser_on.t{1}));
        all_lat_2 =  NaN(2, length(laser_on.t{1}));
        all_lat_shuf2 =  NaN(2, length(laser_on.t{1}),nShuf);
        all_count = NaN(2,length(laser_on.t{1}));
        all_count_shuf = NaN(2,length(laser_on.t{1}), nShuf);
        %         shuf_lat  = NaN(nShuf, length(laser_on.t{1}), n_phases);
        all_phase_labels = NaN(1,length(laser_on.t{1}));
        
        for iPhase = 1:n_phases
            [~, this_phase_idx] = find(stim_phase > edges(iPhase) & stim_phase < edges(iPhase+1));
            
            all_phase_labels(this_phase_idx) = iPhase;
                        
            phase_labels{iPhase} = ['p' num2str(iPhase)]; % used for the xlabels
            
            all_lat_2(1,this_phase_idx) = iPhase;
            
            for iEvt = this_phase_idx
                %                 all_lat(iPhase, iEvt) = this_S.t{1}(nearest_idx3(laser_on.t{1}(iEvt), this_S.t{1}, 1)) - laser_on.t{1}(iEvt);
                % correct for small chance of the last spike occuring
                % before the last laser stim.
                if this_S.t{1}(nearest_idx3(laser_on.t{1}(iEvt), this_S.t{1}, 1)) > laser_on.t{1}(iEvt)
                    all_lat_2(2, iEvt) = this_S.t{1}(nearest_idx3(laser_on.t{1}(iEvt), this_S.t{1}, 1)) - laser_on.t{1}(iEvt);
                else
                    all_lat_2(2, iEvt) = NaN;
                end
                %                 if all_lat(iPhase, iEvt) >cfg_ccf.max_t
                if all_lat_2(2, iEvt) >cfg_ccf.max_t
                    %                     all_lat(iPhase, iEvt) = NaN;
                    all_lat_2(2, iEvt) = NaN;
                end
            end
            
            
            
            % get shuffle for each phase
            for iS = nShuf:-1:1
                %                 all_lat_shuf2(1,this_phase_idx,iS) = iPhase;
                all_lat_shuf2(1,:,iS) = all_lat_2(1,randperm(length(all_lat_2(1,:))));
                %                 [vals, this_phase_idx] = find(stim_phase_shuf(iS,:) > edges(iPhase) & stim_phase_shuf(iS,:) < edges(iPhase+1));
                [~, this_phase_idx] = find(all_lat_shuf2(1,:,iS)==iPhase);
                
                for iEvt = this_phase_idx
                    %                     shuf_lat(iS, iEvt, iPhase) = this_S.t{1}(nearest_idx3(laser_on.t{1}(iEvt), this_S.t{1}, 1)) - laser_on.t{1}(iEvt);
                    if this_S.t{1}(nearest_idx3(laser_on.t{1}(iEvt), this_S.t{1}, 1)) > laser_on.t{1}(iEvt)
                        all_lat_shuf2(2, iEvt, iS) = this_S.t{1}(nearest_idx3(laser_on.t{1}(iEvt), this_S.t{1}, 1)) - laser_on.t{1}(iEvt);
                    else%                     if shuf_lat(iS, iEvt, iPhase) >cfg_ccf.max_t
                        all_lat_shuf2(2, iEvt, iS) = NaN;
                    end
                    if all_lat_shuf2(2, iEvt, iS) >cfg_ccf.max_t
                        
                        %                         shuf_lat(iS, iEvt, iPhase) = NaN;
                        all_lat_shuf2(2, iEvt, iS) = NaN;
                    end
                end
            end
        end
        
        %% center on the 'main peak'
        n_drops = zeros(n_phases, length(all_lat_2));
        
        for iPhase = 1:n_phases
            [~, this_phase_idx] = find(all_lat_2(1,:)==iPhase); % get the index of values for this phase
            bin_counts = histc(all_lat_2(2,this_phase_idx),xbin_centers); %get a hist of the
            [~,peak_idx] = max(bin_counts);
            this_peak_time{iPhase}  =[(xbin_centers(peak_idx)-0.001), (xbin_centers(peak_idx)+0.05)];
            hist_lim{iPhase} = [(xbin_centers(peak_idx)-0.001) (xbin_centers(peak_idx)+0.001)];
            all_count(1,this_phase_idx) = iPhase; % keep track of phase assignments for counts of cells per stim
            all_count_shuf(1,this_phase_idx) = iPhase; % keep track of phase assignments for counts of cells per stim for the shuffle
            
            % cycle all of this phase segment and remove values outside of
            % the main peak
            for iEvt = this_phase_idx
                if isnan(all_lat_2(2,iEvt))
                    continue
                else
                    if all_lat_2(2,iEvt) > (xbin_centers(peak_idx)+0.001) || all_lat_2(2,iEvt) < (xbin_centers(peak_idx)-0.001)
                        n_drops(iPhase, iEvt) = 1;
                        all_lat_2(2,iEvt) = NaN;
                    end
                    % get the number of spikes that happen within the
                    % 'main' spike window as determined by the max in
                    % the histogram of spike latencies
                    if all_lat_2(2,iEvt) < (xbin_centers(peak_idx)+0.05) && all_lat_2(2,iEvt) > (xbin_centers(peak_idx)-0.001)
                        % get the spike count
                        this_event_spikes = restrict(this_S, this_peak_time{iPhase}(1)+laser_on.t{1}(iEvt), this_peak_time{iPhase}(2)+laser_on.t{1}(iEvt));
                        all_count(2,iEvt) = length(this_event_spikes.t{1});
                    end
                end
            end
            fprintf('\n %.1f events dropped in phase #', sum(n_drops(iPhase,:)))
            
            % same for shuffle
            for iS = nShuf:-1:1
                [~, this_phase_idx] = find(all_lat_shuf2(1,:,iS)==iPhase);
                
                for iEvt = this_phase_idx
                    if isnan(all_lat_shuf2(2,iEvt,iS))
                        continue
                    else
                        if all_lat_shuf2(2,iEvt,iS) > (xbin_centers(peak_idx)+0.001) || all_lat_shuf2(2,iEvt,iS) < (xbin_centers(peak_idx)-0.001)
                            all_lat_shuf2(2,iEvt,iS) = NaN;
                        end
                        
                        %get the number of spikes
                        if all_lat_shuf2(2,iEvt,iS) < (xbin_centers(peak_idx)+0.05) && all_lat_shuf2(2,iEvt,iS) > (xbin_centers(peak_idx)-0.001)
                            % get the spike count
                            this_event_spikes = restrict(this_S, this_peak_time{iPhase}(1)+laser_on.t{1}(iEvt), this_peak_time{iPhase}(2)+laser_on.t{1}(iEvt));
                            all_count_shuf(2,iEvt,iS) = length(this_event_spikes.t{1});
                        end
                    end
                end
                
            end
            
        end
        
        %         %% center on the 'main peak'
        %         n_drops = zeros(size(all_lat));
        %
        %         for iPhase = 1:n_phases
        %             bin_counts = histc(all_lat(iPhase,:),xbin_centers);
        %             [~,peak_idx] = max(bin_counts);
        %             this_peak_time{iPhase}  =[(xbin_centers(peak_idx)-0.001), (xbin_centers(peak_idx)+0.001)];
        %
        %             for iEvt = 1:length(all_lat(iPhase,:))
        %                 if isnan(all_lat(iPhase, iEvt))
        %                     continue
        %                 else
        %                     if all_lat(iPhase, iEvt) > (xbin_centers(peak_idx)+0.001) || all_lat(iPhase, iEvt) < (xbin_centers(peak_idx)-0.001)
        %                         n_drops(iPhase, iEvt) = 1;
        %                         all_lat(iPhase, iEvt) = NaN;
        %                     else
        %                         %                         % get the spike count
        %                         %                         this_event_spikes = restrict(this_S, this_peak_time{iPhase}(1)+laser_on.t{1}(iEvt), this_peak_time{iPhase}(2)+laser_on.t{1}(iEvt));
        %                         %                         all_count(iPhase, iEvt) = length(this_event_spikes.t{1});
        %                     end
        %                 end
        %             end
        %             fprintf('\n %.1f events dropped in phase #', sum(n_drops(iPhase,:)))
        %
        %             % same for shuffle
        %             for iS = nShuf:-1:1
        %                 for iEvt = 1:length(shuf_lat(iS,:,iPhase))
        %                     if isnan(shuf_lat(iS,iEvt,iPhase))
        %                         continue
        %                     else
        %                         if shuf_lat(iS,iEvt,iPhase) > (xbin_centers(peak_idx)+0.001) || shuf_lat(iS,iEvt,iPhase) < (xbin_centers(peak_idx)-0.001)
        %                             shuf_lat(iS,iEvt,iPhase) = NaN;
        %                         end
        %                     end
        %                 end
        %                 %                 % get the spike count
        %                 %                 this_event_spikes = restrict(this_S, this_peak_time{iPhase}(1)+laser_on.t{1}(iEvt), this_peak_time{iPhase}(2)+laser_on.t{1}(iEvt));
        %                 %                 all_count_shuf(iPhase, iEvt) = length(this_event_spikes.t{1});
        %             end
        %
        %         end
        
        
        %% get the stats for the shuffles
        
        all_shuf_mean = nanmean(reshape(all_lat_shuf2(2,:,:), 1, numel(all_lat_shuf2(2,:,:))));
        all_shuf_std = nanstd(reshape(all_lat_shuf2(2,:,:), 1, numel(all_lat_shuf2(2,:,:))));
        
        
        all_count_shuf_mean = nanmean(reshape(all_count_shuf(2,:,:), 1, numel(all_count_shuf(2,:,:))));
        all_count_shuf_std = nanstd(reshape(all_count_shuf(2,:,:), 1, numel(all_count_shuf(2,:,:))));
        
        %         for iPhase = 1:n_phases
        %             [~, this_phase_idx] = find(all_lat_2(1,:,:)==iPhase); % get the index of values for this phase
        %
        %             this_shuf = shuf_lat(:,:,iPhase);
        %             all_shuf_phase_mean(iPhase) = nanmean(reshape(this_shuf, 1, numel(this_shuf)));
        %             all_shuf_phase_std(iPhase) = nanstd(reshape(this_shuf, 1, numel(this_shuf)));
        %
        %         end
        
        
        %                 Old version
        %         all_shuf_mean = nanmean(reshape(shuf_lat, 1, numel(shuf_lat)));
        %         all_shuf_std = nanstd(reshape(shuf_lat, 1, numel(shuf_lat)));
        %
        %         for iPhase = 1:n_phases
        %             this_shuf = shuf_lat(:,:,iPhase);
        %             all_shuf_phase_mean(iPhase) = nanmean(reshape(this_shuf, 1, numel(this_shuf)));
        %             all_shuf_phase_std(iPhase) = nanstd(reshape(this_shuf, 1, numel(this_shuf)));
        %
        %         end
        
        %% convert daata to ms
        all_lat_2(2,:) = all_lat_2(2,:)*1000;
        all_lat_shuf2(2,:,:) = all_lat_shuf2(2,:,:)*1000;
        xbin_centers = xbin_centers*1000;
        all_shuf_mean = all_shuf_mean*1000;
        all_shuf_std = all_shuf_std*1000;
        %         all_lat = all_lat*1000;
        %         all_shuf_mean = all_shuf_mean*1000;
        %         all_shuf_std = all_shuf_std*1000;
        %         xbin_centers = xbin_centers*1000;
        %% make a figure
        figure(1);
        c_ord = linspecer(n_phases);
        
        % try to make a colored sine wave to match the phases.
        subplot(6, ceil(n_phases/2), 1)
        hold on
        x = -pi:pi/(50*n_phases):pi;
        wave_phase = sin(x);
        for iPhase = 1:n_phases
            
            plot(-99+(100*iPhase):(100*iPhase), wave_phase(-99+(100*iPhase):(100*iPhase)), 'color', c_ord(iPhase,:), 'linewidth', 5)
            
        end
        axis off
        title(sprintf('%.1f to %.1f Hz', f_list{iF}(1), f_list{iF}(2)), 'fontsize', font_size)
        
        
        subplot(6, ceil(n_phases/2), [n_phases+2  4*ceil(n_phases/2)])
        hold on
        he = errorbar(0, all_shuf_mean, 2*all_shuf_std, 'o','color', 'k','MarkerFaceColor', 'k', 'markersize', 10);
        
        for iPhase = 1:n_phases
            [~, this_phase_idx] = find(all_lat_2(1,:)==iPhase); % get the index of values for this phase
            
            subplot(6, ceil(n_phases/2), iPhase+1)
            %             h = histogram(all_lat(iPhase,:), xbin_centers);
            h = histogram(all_lat_2(2,this_phase_idx), xbin_centers);
            h.FaceColor = c_ord(iPhase,:);
            h.EdgeColor = c_ord(iPhase,:);
            ylim([0 20])
            xlim(hist_lim{iPhase}*1000)
            xlabel('latency (ms)')
            %             title(sprintf('N stims: %.1f', N_stims{iPhase}), 'fontsize', font_size)
            %             legend boxoff
            set(gca, 'fontsize', font_size)
            
            subplot(6, ceil(n_phases/2), [n_phases+2  4*ceil(n_phases/2)])
            hold on
            %             he = errorbar(iPhase, nanmean(all_lat(iPhase,:), 2), 2*nanstd(all_lat(iPhase,:), 0,2), 'o', 'color', c_ord(iPhase,:),'MarkerFaceColor', c_ord(iPhase,:), 'markersize', 10);
            he = errorbar(iPhase, nanmean(all_lat_2(2,this_phase_idx), 2), 2*nanstd(all_lat_2(2,this_phase_idx), 0,2), 'o', 'color', c_ord(iPhase,:),'MarkerFaceColor', c_ord(iPhase,:), 'markersize', 10);
            
            set(gca, 'fontsize', font_size);
            %             xlabel('Phase')
            xlim([-0.5 n_phases+0.5])
            ylabel('Mean latency (ms)')
            %             set(gca, 'xticklabel', ['shuf' phase_labels])
            set(gca, 'xticklabel', '')
        end
        
        %         SetFigure([], gcf)
        %         set(gcf, 'position', [435, 50, 949, 697]);
        %
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.35, 0.98,[ExpKeys.subject '_' ExpKeys.date '  Cell ' cell_id], 'fontsize', font_size)
        %         saveas(gcf, [all_fig_dir ExpKeys.subject '_' ExpKeys.date '_' cell_id(1:end-3) '_f' num2str(floor(f_list{iF}(1))) '_' num2str(floor(f_list{iF}(2))) '_latency.png']);
        %         saveas_eps([ExpKeys.subject '_' ExpKeys.date '_' cell_id(1:end-3) '_f' num2str(floor(f_list{iF}(1))) '_' num2str(floor(f_list{iF}(2))) '_latency'], all_fig_dir)
        
        % make a count figure
        subplot(6, ceil(n_phases/2), [n_phases+8  6*ceil(n_phases/2)])
        hold on
        he = errorbar(0, all_count_shuf_mean, 2*all_count_shuf_std, 'o','color', 'k','MarkerFaceColor', 'k', 'markersize', 10);
        
        for iPhase = 1:n_phases
            [~, this_phase_idx] = find(all_count(1,:)==iPhase); % get the index of values for this phase
            
            hold on
            he = errorbar(iPhase, nanmean(all_count(2,this_phase_idx), 2), 2*nanstd(all_count(2,this_phase_idx), 0,2), 'o', 'color', c_ord(iPhase,:),'MarkerFaceColor', c_ord(iPhase,:), 'markersize', 10);
            set(gca, 'fontsize', font_size);
            xlabel('Phase')
            xlim([-0.5 n_phases+0.5])
            ylabel('Spike count')
            set(gca, 'xticklabel', ['shuf' phase_labels])
        end
        
        SetFigure([], gcf)
        set(gcf, 'position', [440    37   930   761]);

        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.35, 0.98,[ExpKeys.subject '_' ExpKeys.date '  Cell ' cell_id], 'fontsize', font_size)
        saveas(gcf, [all_fig_dir ExpKeys.subject '_' ExpKeys.date '_' cell_id(1:end-3) '_f' num2str(floor(f_list{iF}(1))) '_' num2str(floor(f_list{iF}(2))) '_latency.png']);
        saveas_eps([ExpKeys.subject '_' ExpKeys.date '_' cell_id(1:end-3) '_f' num2str(floor(f_list{iF}(1))) '_' num2str(floor(f_list{iF}(2))) '_latency'], all_fig_dir)
        %
                close all
%                 if strfind(f_list_label{iF}), '.')
                    freq_label = ['f_' strrep(strrep(f_list_label{iF}, ' ', ''), '-', '_')];
%                 else
                    
%                 end
    out.(cell_id).(freq_label).latency = all_lat_2;
    out.(cell_id).(freq_label).latency_shuf = all_lat_shuf2;
    out.(cell_id).(freq_label).count = all_count;
    out.(cell_id).(freq_label).count_shuf = all_count_shuf;
    end % end freq loop

end % end cell loop

    save([ all_lat_dir  ExpKeys.subject '_' ExpKeys.date '_lat.mat'], 'out', '-v7.3')


end % end of function
