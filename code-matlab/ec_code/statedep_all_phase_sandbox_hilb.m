function statedep_all_phase_sandbox_hilb
% restoredefaultpath;

if isunix
    addpath(genpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'));
    addpath('/Users/jericcarmichael/Documents/GitHub/EC_state/Basic_functions');
    
    all_fig_dir = '/Volumes/Fenrir/State_dep/all_checks/';
    all_ccf_dir = '/Volumes/Fenrir/State_dep/all_ccf/';
    
    
else
    %     addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
    %addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
    %     addpath('C:\Users\mvdm\Documents\GitHub\EC_state\Basic_functions');
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\EC_State'))
    %addpath('D:\My_Documents\GitHub\EC_state\Basic_functions');
    all_fig_dir = 'G:\State_data\all_checks\';
    all_ccf_dir = 'G:\State_data\all_ccf\';
    
    %     %cd('D:\data\EC_state\M14_2018-12-01_vStr_light_min');
    %     cd('C:\data\state-dep\M14_2018-12-01_vStr_light_min');
    %     %cd('D:\data\EC_state\M13-2018-12-05_dStr_2p2_light_min');
    %     cd('C:\data\state-dep\M13-2018-12-05_dStr_2p2_light_min');
    
end

mkdir(all_ccf_dir); mkdir(all_fig_dir); mkdir([all_fig_dir '/CCFs'])
%% defaults
font_size = 18;
% for naming figures;
if isunix
    fname = strsplit(cd, '/');
else
    fname = strsplit(cd, '\');
end
fname = fname{end};
fname = fname(1:strfind(fname,'p')+1);
fname = strrep(fname, '-', '_');
%% load CSC
cfg = [];
%cfg.decimateByFactor = 30;
cfg.fc = {'CSC22.ncs'};
this_csc = LoadCSC(cfg);
Fs = 1 ./ median(diff(this_csc.tvec));

%% make psd
if Fs >=30000
    wSize =  32768;%16384;
else
    wSize = 1024;
end

[Pxx,F] = pwelch(diff(this_csc.data), rectwin(wSize), wSize/8, [], Fs);
    %% look at the psd
    figure(3)
    subplot(3,2,1);
    plot(F, 10*log10(Pxx), 'k', 'LineWidth', 2);
    set(gca, 'XLim', [0 100], 'FontSize', font_size); grid on;
    xlabel('Frequency (Hz)');
    title(this_csc.label{1},'fontsize', font_size);
    SetFigure([], gcf)
    %% load expkeys
    LoadExpKeys
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

%% restric the CSC to the run portion
this_csc = restrict(this_csc, start_stop.t{1}(main_rec_idx), start_stop.t{2}(main_rec_idx));


%% inspect stim artifact
this_stim_binned = zeros(size(this_csc.tvec));
idx = nearest_idx3(laser_on.t{1}, this_csc.tvec);
this_spk_binned(idx) = 1;
[xc, tvec] = xcorr(this_csc.data, this_spk_binned, 500);
tvec = tvec ./ Fs;

figure(101)
plot(tvec, xc, 'k', 'LineWidth', 2);

%% load spikes
cfg = []; cfg.getTTnumbers = 0;
S = LoadSpikes(cfg);

%% select a cell
for iC = 1:length(S.label)
    
    this_S = SelectTS([], S, iC);
    
    %% ID plot
    figure(3)
    subplot(3,3,1);
    % h(1) = plot(1:10,nan(1,10), 'color', [1 1 1]);
    
    cell_id = this_S.label{1}(1:end-2);
    cell_id = strrep(cell_id, '_SS', '');
    
    hdr = [];
    hdr.subject = fname(1:3);
    hdr.date = fname(5:14);
    hdr.depth = str2double([fname(strfind(fname, 'p')-1) '.' fname(strfind(fname, 'p')-+1)]); % convert the depth from the title p = . in file names used here
    if hdr.depth < 2; % if less than 2mm define as cortex
        hdr.target = 'crtx';
    elseif hdr.depth > 2 && hdr.depth < 4 % define dorsal striatum
        hdr.target = 'dStr';
    elseif hdr.depth > 4 % define as vStr
        hdr.target = 'vStr';
    end
    
    text(0,0.8, [hdr.subject ' ' strrep(hdr.date, '_', '-') ' '], 'FontSize', font_size);
    text(0, 0.4, [hdr.target ' depth: ' num2str(hdr.depth) 'mm'], 'FontSize', font_size);
    text(0, 0.0, ['Cell: ' strrep(this_S.label{1}, '_', '-')], 'FontSize', font_size)
    
    
    axis off
    %% overall PETH
    cfg = [];
    cfg.binsize = 0.0005;
    cfg.max_t = 0.02;
    [this_ccf, tvec] = ccf(cfg, laser_on.t{1}, this_S.t{1});
    
    figure(3);
    subplot(332);
    plot(tvec, this_ccf, 'k', 'LineWidth', 2);
    h = title('all stim PETH', 'fontsize', font_size);
    set(h, 'Interpreter', 'none');
    set(gca, 'FontSize', font_size); xlabel('time (s)'); ylabel('spike count');
    

    %% get some LFP phases (filtfilt)
    f_list = {[3 5], [6.5 9.5],[15 25], [30 40],[40 60], [60 80]};
    f_list_label = {'3 - 5', '7 - 10', '15 - 25', '30 - 40', '40 - 60', '60 - 80'};
    nShuf = 100;
    sub4_id = 1; % allows for colum looping in Fig 4
    
    for iF = 1:length(f_list) % loop across freqs
        %%  filter & hilbertize
        cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
        csc_f = FilterLFP(cfg_filt, this_csc);
        
        
        %         %% overall PETH
        %         cfg = [];
        %         cfg.binsize = 0.0005;
        %         cfg.max_t = 0.02;
        %         [this_ccf, tvec] = ccf(cfg, laser_on.t{1}, this_S.t{1});
        %
        %         figure(3);
        %         subplot(332);
        %         plot(tvec, this_ccf, 'k', 'LineWidth', 2); h = title(sprintf('all stim PETH %s', this_S.label{1})); set(h, 'Interpreter', 'none');
        %         set(gca, 'FontSize', font_size); xlabel('time (s)'); ylabel('spike count');
        %         title(strcat('Cell id:    ', this_S.label));
        %
        
        
        %% SHUFFLE WOULD DO SOME KIND OF RANDOM CIRCSHIFT HERE %%%
        
        csc_f.data = angle(hilbert(csc_f.data));
        hold on
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
        
        % get phase for each laser stim (actual)
        stim_phase_idx = nearest_idx3(laser_on.t{1}, csc_f.tvec);
        stim_phase = csc_f.data(stim_phase_idx);
        
        % STIM PHASE HISTO THIS IS IMPORTANT
        figure(2); subplot(2, 3, iF);
        hist(stim_phase, 36); title(sprintf('stim phase histo (%.1f-%.1f Hz)', f_list{iF}(1), f_list{iF}(2)), 'fontsize', font_size);
        
        
        %%  phase split
        % instead use Histc to split the phases by some number of bins.
        % Then take all the events that occur with thise bin edges and loop
        % over those bins.  In the loop I will run the CCF on the events
        % within that phase bin and get a ccf which will make an events x
        % phase bin array.  This will then become the imagesc.
        
        % first use a simple version for checks
        phase_low_idx = find(stim_phase < 0);
        phase_high_idx = find(stim_phase >= 0);
        
        [this_ccf_low, ~] = ccf(cfg, laser_on.t{1}(phase_low_idx), this_S.t{1});
        [this_ccf_high, ~] = ccf(cfg, laser_on.t{1}(phase_high_idx), this_S.t{1});
        
        for iS = nShuf:-1:1
            
            phase_low_idx = find(stim_phase_shuf(iS,:) < 0);
            phase_high_idx = find(stim_phase_shuf(iS,:) >= 0);
            
            [ccf_low_shuf(iS,:), ~] = ccf(cfg, laser_on.t{1}(phase_low_idx), this_S.t{1});
            [ccf_high_shuf(iS,:), ~] = ccf(cfg, laser_on.t{1}(phase_high_idx), this_S.t{1});
            
        end % of shuffles
        
        % complex version with phase splits
        
        %get some division of phases
        n_phases = 5;
        c_ord = linspecer(n_phases);
        fprintf('\nPhase split into %1d bins\n', n_phases)
        [~, edges, ~] = histcounts(-pi:pi, n_phases, 'BinLimits', [-pi, pi]);
        
        for iPhase = 1:length(edges)-1
            
            [vals, this_phase_idx] = find(stim_phase > edges(iPhase) & stim_phase < edges(iPhase+1));
            N_stims{iF} = sum(vals);
            
            [all_ccf(iPhase,:), tvec] = ccf(cfg, laser_on.t{1}(this_phase_idx), this_S.t{1});
            
            % same for shuffles
            for iS = nShuf:-1:1
                
                this_phase_idx = find(stim_phase_shuf(iS,:) > edges(iPhase) & stim_phase_shuf(iS,:) < edges(iPhase+1));
                %                 phase_high_idx = find(stim_phase_shuf(iS,:) >= 0);
                
                [all_ccf_shuf(iPhase,:,iS), ~] = ccf(cfg, laser_on.t{1}(this_phase_idx), this_S.t{1});
                %                 [ccf_high_shuf(iS,:), ~] = ccf(cfg, laser_on.t{1}(phase_high_idx), this_S.t{1});
                
            end % of shuffles
            temp_all_ccf(iPhase,:) = all_ccf(iPhase,:);
            temp_all_ccf(iPhase,temp_all_ccf(iPhase,:) <= squeeze(nanmean(all_ccf_shuf(iPhase,:,:),3) + 1.96*nanstd(all_ccf_shuf(iPhase,:,:), [], 3))) = NaN;
            all_ccf_vShuff(iPhase,:) = temp_all_ccf(iPhase, :);
        end
        %% multi-phase breakdown
        figure(111)
        subplot(ceil(length(f_list)/2),2, iF);
        ms_tvec = tvec*1000; %convert to ms
        for iPhase = 1:length(edges)-1
%             h(3)= plot(tvec, nanmean(all_ccf_shuf(iPhase,:,:),3), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
%             h(4)= plot(tvec, nanmean(all_ccf_shuf(iPhase,:,:),3) + 1.96*nanstd(all_ccf_shuf(iPhase,:,:),[],3), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
%             h(5) = plot(tvec, nanmean(all_ccf_shuf(iPhase,:,:),3) - 1.96*nanstd(all_ccf_shuf(iPhase,:,:),[],3), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
                         h(1) = plot(ms_tvec, all_ccf(iPhase,:), '-', 'color', c_ord(iPhase,:), 'LineWidth', 1); hold on;

            max_val(iPhase) = max(all_ccf(iPhase,:));
            max_shuf_val(iPhase) = max(nanmean(all_ccf_shuf(iPhase,:,:),3));
            
        end
         set(gca, 'FontSize', font_size); xlabel('time (ms)'); ylabel('response %');
        title(sprintf('%.1f - %.1f Hz', f_list{iF}(1), f_list{iF}(2)), 'fontsize', font_size);
        % get a zoom in
        [~, max_idx] = max(diff(all_ccf(iPhase,:)));
%         xlim([tvec(max_idx)-0.0005 tvec(max_idx)+0.002]);
        val_max = max(max_val);
        val_min = min(max_val);
%         ylim([val_max-0.1 val_max+0.02])
SetFigure([], gcf)
    set(gcf, 'position', [600 50 640*2 420*2]);
        xlim([-5 15])
        ymax = ylim;
                rectangle('position', [0, 0, 1, ymax(2)], 'facecolor',[([4,172,218]./255) 0.5], 'edgecolor',[([4,172,218]./255) 0.5] )
        set(gca, 'ylim', ymax);
        MagInset(gcf, gca, [ms_tvec(max_idx)-0.5,ms_tvec(max_idx)+2,val_min-0.05, val_max+0.005],[ms_tvec(max_idx)+5,ms_tvec(max_idx)+10,val_max-0.2, val_max+0.005], {'NE','NW'; 'SE','SW'}) 

% Square_subplots
        %%
        
        
        % hold on to values across frequencies
        ccf_out.(cell_id).all_f_ccf{iF} = all_ccf;
        %         ccf_out.all_f_ccf_Z{iF} = zscore(all_ccf, 1);
        ccf_out.(cell_id).all_f_ccf_vShuff{iF} = all_ccf_vShuff;
        ccf_out.(cell_id).all_f_ccf_shuf{iF} = all_ccf_shuf;
        %         all_f_ccf_VShuff{iF} =
        
        figure(4)
        subplot( length(f_list)/2,4, sub4_id); sub4_id = sub4_id+1;
        
        plot(tvec, this_ccf_low - this_ccf_high, 'g', 'LineWidth', 2); hold on;
        
        plot(tvec, nanmean(ccf_low_shuf - ccf_high_shuf), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        plot(tvec, nanmean(ccf_low_shuf - ccf_high_shuf) + 1.96*nanstd(ccf_low_shuf - ccf_high_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        plot(tvec, nanmean(ccf_low_shuf - ccf_high_shuf) - 1.96*nanstd(ccf_low_shuf - ccf_high_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        
        set(gca, 'FontSize', font_size); xlabel('time (s)'); ylabel('spike count');
        title(sprintf('phase split %.1f-%.1f Hz', f_list{iF}(1), f_list{iF}(2)), 'fontsize', font_size);
        
        subplot(length(f_list)/2,4, sub4_id); sub4_id = sub4_id+1;
        imagesc( tvec, 1:length(edges), all_ccf)
        axis xy
        caxis
        set(gca, 'ytick', [1 length(edges)/2 length(edges)], 'yticklabel', {'-\pi' '0' '\pi'}, 'fontname', 'helvetica',  'FontSize', font_size);
        xlabel('time (s)'); ylabel('phase');
        
        %% same thing but relative to the shuffle
        %         [all_ccf_zscore, mu, sigma] = zscore(all_ccf, 0, 2);
        
        figure(6)
        subplot(ceil(length(f_list)/2),2, iF);
        imagesc( tvec, 1:length(edges), all_ccf_vShuff)
        axis xy
        %         caxis([nanmean(ccf__shuf) + 1.96*nanstd(ccf_low_shuf)])
        set(gca, 'ytick', [1 length(edges)/2 length(edges)], 'yticklabel', {'-\pi' '0' '\pi'}, 'fontname', 'helvetica',  'FontSize', font_size);
        xlabel('time (s)'); ylabel('phase');
        title(sprintf('ccf exceeding shuffle (%.1f-%.1f Hz) Stims = %.1f', f_list{iF}(1), f_list{iF}(2), N_stims{iF}), 'fontsize', font_size);
        
        
        
        %%
        figure(3);
        subplot(3, 3, 3 + iF);
        h(1) = plot(tvec, this_ccf_low, 'b', 'LineWidth', 2); hold on;
        h(2) = plot(tvec, this_ccf_high, 'r', 'LineWidth', 2);
        
        h(3)= plot(tvec, nanmean(ccf_low_shuf), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h(4)= plot(tvec, nanmean(ccf_low_shuf) + 1.96*nanstd(ccf_low_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h(5) = plot(tvec, nanmean(ccf_low_shuf) - 1.96*nanstd(ccf_low_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        
        legend(h, {'phase < 0', 'phase >= 0'}, 'Location', 'Northwest'); legend boxoff;
        set(gca, 'FontSize', font_size); xlabel('time (s)'); ylabel('spike count');
        title(sprintf('%.1f-%.1f Hz nStims: %.1f', f_list{iF}(1), f_list{iF}(2),  N_stims{iF}), 'fontsize', font_size);
        
        %         figure(3);
        %         subplot(3, 3, 3 + iF);
        %         plot(tvec, this_ccf_low - this_ccf_high, 'g', 'LineWidth', 2); hold on;
        %
        %         plot(tvec, nanmean(ccf_low_shuf - ccf_high_shuf), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        %         plot(tvec, nanmean(ccf_low_shuf - ccf_high_shuf) + 1.96*nanstd(ccf_low_shuf - ccf_high_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        %         plot(tvec, nanmean(ccf_low_shuf - ccf_high_shuf) - 1.96*nanstd(ccf_low_shuf - ccf_high_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        %
        %         set(gca, 'FontSize', font_size); xlabel('time (s)'); ylabel('spike count');
        %         title(sprintf('phase split diff %.1f-%.1f Hz', f_list{iF}(1), f_list{iF}(2)));
        
    end % of freq loop
    
    %% try a stacked imagesc of the phase binned ccf
    for iPhase =1:length(ccf_out.(cell_id).all_f_ccf)
        M(:,:,iPhase) = ccf_out.(cell_id).all_f_ccf{iPhase};
    end
    
    figure(5)
    subplot(1,2,1)
    hs = slice(M,[],[],1:length(ccf_out.(cell_id).all_f_ccf)) ;
    colormap('Parula')
    %     shading interp
    set(hs,'FaceAlpha',1)
    x_val = get(gca, 'xtick');
    title('Raw', 'fontsize', font_size)
    set(gca,  'xtick', [x_val(1) x_val(ceil(length(x_val)/2)) x_val(end)],'xticklabel', [tvec(1) tvec(ceil(length(tvec)/2)) tvec(end)],...
        'ytick', [1 length(edges)/2 length(edges)], 'Yticklabel', {'-\pi' '0' '\pi'}, ...
        'ztick', [1:length(ccf_out.(cell_id).all_f_ccf)], 'zticklabel', f_list_label, 'fontsize', font_size, 'fontname', 'helvetica');
    xlabel('time (s)');
    ylabel('phase');
    zlabel('frequency band (Hz)');
    view(-4.5, 14)
    
    % same thing for exceeding shuffle
    for iPhase =1:length(ccf_out.(cell_id).all_f_ccf_vShuff)
        this_phase_ccf_vShuf = ccf_out.(cell_id).all_f_ccf_vShuff{iPhase};
        this_phase_ccf_vShuf(isnan(this_phase_ccf_vShuf)) = 0;
        M(:,:,iPhase) = this_phase_ccf_vShuf;
    end
    
    figure(5)
    subplot(1,2,2)
    hs = slice(M,[],[],1:length(ccf_out.(cell_id).all_f_ccf)) ;
    colormap('Parula')
    %     shading interp
    set(hs,'FaceAlpha',1)
    title('Exceeding 2SD of shuffle', 'fontsize', font_size)
    x_val = get(gca, 'xtick');
    set(gca,  'xtick', [x_val(1) x_val(ceil(length(x_val)/2)) x_val(end)],'xticklabel', [tvec(1) tvec(ceil(length(tvec)/2)) tvec(end)],...
        'ytick', [1 length(edges)/2 length(edges)], 'Yticklabel', {'-\pi' '0' '\pi'}, ...
        'ztick', [1:length(ccf_out.(cell_id).all_f_ccf_vShuff)], 'zticklabel', f_list_label, 'fontsize', font_size, 'fontname', 'helvetica');
    xlabel('time (s)');
    ylabel('phase');
    zlabel('frequency band (Hz)');
    view(-4.5, 14)
    
    
    %% save the split diff plot as a summary
    figure(2)
    cfg_fig.ft_size = font_size;
    SetFigure(cfg_fig, gcf)
    tightfig
    set(gcf, 'position', [235 158  1121 528]);
    %     saveas(gcf, [fname '_' cell_id '_hist.png']);
    saveas(gcf, [all_fig_dir fname '_' cell_id '_hist.png']);
    %     saveas_eps([fname '_' cell_id '_hist'], cd)
    saveas_eps([fname '_' cell_id '_hist'], all_fig_dir)
    
    figure(3)
    cfg_fig.ft_size = font_size;
    SetFigure(cfg_fig, gcf)
    set(gcf, 'position', [600 50 640*2 420*2]);
    %     saveas(gcf, [fname '_' cell_id '_ccf.png']);
    saveas(gcf, [all_fig_dir fname '_' cell_id '_ccf.png']);
    saveas(gcf, [all_fig_dir '/CCFs/' fname '_' cell_id '_ccf.png']);
    %     saveas_eps([fname '_' cell_id '_ccf'], cd)
    saveas_eps([fname '_' cell_id '_ccf'], all_fig_dir)
    
    figure(4)
    tightfig
    cfg_fig.ft_size = font_size;
    SetFigure(cfg_fig, gcf)
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.25, 0.98,[fname '  Cell ' cell_id], 'fontsize', font_size)
    set(gcf, 'position', [600 50 640*2 420*2]);
    %     saveas(gcf, [fname '_' cell_id '_both.png']);
    saveas(gcf, [all_fig_dir fname '_' cell_id '_both.png']);

    %     saveas_eps([fname '_' cell_id '_both'], cd)
    saveas_eps([fname '_' cell_id '_both'], all_fig_dir)
    
    
    % stack of raw imagescs
    figure(5)
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.35, 0.98,[fname '  Cell ' cell_id], 'fontsize', font_size)
    SetFigure([], gcf)
    set(gcf, 'position', [100    38   1340   751]);
    %     saveas(gcf, [fname '_' cell_id '_stack.png']);
    saveas(gcf, [all_fig_dir fname '_' cell_id '_stack.png']);
    %     saveas_eps([fname '_' cell_id '_stack'], cd)
    saveas_eps([fname '_' cell_id '_stack'], all_fig_dir)
    
    %only imagesc for ccf exceeding shuffle 1.96
    figure(6)
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.35, 0.98,[fname '  Cell ' cell_id], 'fontsize', font_size)
    cfg_fig.ft_size = font_size;
    SetFigure(cfg_fig, gcf)
    set(gcf, 'position', [600 50 640*2 420*2]);
    %     saveas(gcf, [fname '_' cell_id '_all_phase.png']);
    saveas(gcf, [all_fig_dir fname '_' cell_id '_all_phase.png']);
    %     saveas_eps([fname '_' cell_id '_all_phase'], cd)
    saveas_eps([fname '_' cell_id '_all_phase'], all_fig_dir)
    
    
        figure(111)
    cfg_fig.ft_size = font_size;
%     SetFigure(cfg_fig, gcf)
    set(gcf, 'position', [600 50 640*2 420*2]);
    %     saveas(gcf, [fname '_' cell_id '_ccf.png']);
    saveas(gcf, [all_fig_dir fname '_' cell_id '_ccf_5.png']);
    saveas(gcf, [all_fig_dir '/CCFs/' fname '_' cell_id '_ccf_5.png']);
    %     saveas_eps([fname '_' cell_id '_ccf'], cd)
    saveas_eps([fname '_' cell_id '_ccf_5'], all_fig_dir)
    
        close all

    
    %% save the ccf info
%     save([ all_ccf_dir  hdr.subject '_' hdr.date '.mat'], 'ccf_out', '-v7.3')
    
    
    
    fprintf(['Complete Cell #' cell_id '...\n'])
    
end % end of Spike loop

    %% save the ccf info
        save([ all_ccf_dir  hdr.subject '_' hdr.date '.mat'], 'ccf_out', '-v7.3')

end

