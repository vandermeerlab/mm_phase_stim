function STATE_state_peth(cfg_in, S, CSC, evt_in)
%% STATE_state_peth: takes in a spike 'ts' and a csc 'tsd' creates PETHS
% The CSC input will be band-pass filtered between pairs specified in
% cfg.freq (default cfg.freq = [1, 4; 7, 11; 20, 30; 45, 65; 70, 90]).

% Two cfg.modes for either 'power' or 'phase' mode.
%
% In 'power' mode (cfg.mode = 'power') only requires frequencies and either
% cfg.power_mode = 'zscore' or cfg.power_mode = 'raw'
%
% In 'Phase' mode (cfg.mode = 'phase') you can specify how many
% devisions of pi you would like (default cfg.phase_div = 2.
%
%
% Inputs:
%   cfg_in: [struct] = contains config specifications.  Will override
%   cfg_def below.
%   S: ['ts' struct] = contains spike(s)
%   CSC: ['tsd' struct] = contains LFP
%
% EC 2018-12-31

%% set up configs
%
cfg_def = [];
cfg_def.mode = 'phase';
cfg_def.phase_div = 2; % splits phases into >pi or <=pi
cfg_def.freq = {[1, 4], [7, 11], [20, 30], [45, 65], [70, 90]}; % band pass frequency pairs
cfg_def.power_mode = 'zscore';
cfg_def.ft_size = 18;
cfg_def.nShuf = 100;
cfg_def.c_ord = linspecer(4);

cfg = ProcessConfig2(cfg_def, cfg_in);

%%

for iF = 1:length(cfg.freq) % loop across freqs
    
    % filter & hilbertize
    cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = cfg.freq{iF};
    csc_f = FilterLFP(cfg_filt, CSC);
    
    
    fprintf(['Processing Cell #' S.label{1} '...\n'])
    %% overall PETH
    cfg_peth = [];
    cfg_peth.binsize = 0.0005;
    cfg_peth.max_t = 0.02;
    [this_ccf, tvec] = ccf(cfg_peth, evt_in.t{1}, S.t{1});
    
    figure(3);
    ax_1 = subplot(321);
    plot(tvec, this_ccf, 'k', 'LineWidth', 2); h = title(sprintf('all stim PETH %s', S.label{1})); set(h, 'Interpreter', 'none');
    set(gca, 'FontSize', cfg_peth.ft_size); xlabel('time (s)'); ylabel('spike count');
    title(strcat('Cell id:    ', S.label));
    
    
    %% SHUFFLE WOULD DO SOME KIND OF RANDOM CIRCSHIFT HERE %%%
    switch cfg.mode
        
        
        case strcmp(cfg.mode, 'phase')
            disp(['Entering ' cfg.mode ' mode!'])
            csc_f.data = angle(hilbert(csc_f.data));
            
        case strcmp(cfg.mode, 'power')
            disp(['Entering ' cfg.mode ' mode!'])
            csc_f.data = abs(hilbert(csc_f.data));
    end
    hold on
    for iS = cfg.nShuf:-1:1
        
        csc_shuf = csc_f;
        if strcmp(version('-release'), '2014b') % version differences in how circshift is handled.
            csc_shuf.data = circshift(csc_shuf.data, round(rand(1) .* 0.5*length(csc_shuf.data)), 2);
        else
            csc_shuf.data = circshift(csc_shuf.data, round(rand(1) .* 0.5*length(csc_shuf.data)));
        end
        stim_phase_idx = nearest_idx3(evt_in.t{1}, csc_shuf.tvec);
        stim_phase_shuf(iS, :) = csc_shuf.data(stim_phase_idx);
    end % of shuffles
    
    % get phase for each laser stim (actual)
    stim_phase_idx = nearest_idx3(evt_in.t{1}, csc_f.tvec);
    stim_phase = csc_f.data(stim_phase_idx);
    
    % STIM PHASE HISTO THIS IS IMPORTANT
    figure(2); subplot(2, 2, iF);
    hist(stim_phase, 36); title(sprintf('stim phase histo (%.1f-%.1f Hz)', cfg.freq{iF}(1), cfg.freq{iF}(2)));
    
    
    % example phase split
      phase_low_idx = find(-pi <= stim_phase & stim_phase <= -pi/2);
      phase_low_mid_idx = find(-pi/2 < stim_phase & stim_phase <= 0);
      phase_mid_high_idx = find(0 < stim_phase & stim_phase <= pi/2);
      phase_high_idx = find(pi/2 < stim_phase & stim_phase <= pi);

 
      
      [this_ccf_low, tvec] = ccf(cfg, evt_in.t{1}(phase_low_idx), S.t{1});
      [this_ccf_low_mid, ~] = ccf(cfg, evt_in.t{1}(phase_low_mid_idx), S.t{1});
      [this_ccf_mid_high, ~] = ccf(cfg, evt_in.t{1}(phase_mid_high_idx), S.t{1});
      [this_ccf_high, ~] = ccf(cfg, evt_in.t{1}(phase_high_idx), S.t{1});
    
    % same for shuffles
    for iS = cfg.nShuf:-1:1
        
      phase_low_idx = find(-pi <= stim_phase_shuf(iS,:) & stim_phase_shuf(iS,:) <= -pi/2);
      phase_low_mid_idx = find(-pi/2 < stim_phase_shuf(iS,:) & stim_phase_shuf(iS,:) <= 0);
      phase_mid_high_idx = find(0 < stim_phase_shuf(iS,:) & stim_phase_shuf(iS,:) <= pi/2);
      phase_high_idx = find(pi/2 < stim_phase_shuf(iS,:) & stim_phase_shuf(iS,:) <= pi);
        
        [ccf_low_shuf(iS,:), ~] = ccf(cfg, evt_in.t{1}(phase_low_idx), S.t{1});
        [ccf_low_mid_shuf(iS,:), ~] = ccf(cfg, evt_in.t{1}(phase_low_mid_idx), S.t{1});
        [ccf_mid_high_shuf(iS,:), ~] = ccf(cfg, evt_in.t{1}(phase_mid_high_idx), S.t{1});
        [ccf_high_shuf(iS,:), ~] = ccf(cfg, evt_in.t{1}(phase_high_idx), S.t{1});
        
    end % of shuffles
    
    figure(1);
%     subplot(3, 2, 2 + iF);
    h(1) = plot(tvec, this_ccf_low, 'color', cfg.c_ord(1,:) , 'LineWidth', 2); hold on;
    h(2) = plot(tvec, this_ccf_low_mid, 'color', cfg.c_ord(2,:) , 'LineWidth', 2);
    h(3) = plot(tvec, this_ccf_mid_high, 'color', cfg.c_ord(3,:) , 'LineWidth', 2); 
    h(4) = plot(tvec, this_ccf_high, 'color', cfg.c_ord(4,:) , 'LineWidth', 2);
    
    h(5)= plot(tvec, nanmean(ccf_low_shuf), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    h(6)= plot(tvec, nanmean(ccf_low_shuf) + 1.96*nanstd(ccf_low_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    h(7) = plot(tvec, nanmean(ccf_low_shuf) - 1.96*nanstd(ccf_low_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
    legend(h, {'-pi to -pi/2', '-pi/2 to 0', '0 to pi/2','pi/2 to pi'}, 'Location', 'Northwest'); legend boxoff;
    set(gca, 'FontSize', cfg.ft_size); xlabel('time (s)'); ylabel('spike count');
    title(sprintf('phase split %.1f-%.1f Hz', cfg.freq{iF}(1), cfg.freq{iF}(2)));
    
    
    %% got to here on JAN 3rd
    figure(3);
    subplot(3, 2, 2 + iF);
    plot(tvec, this_ccf_low - this_ccf_high, 'g', 'LineWidth', 2); hold on;
    
    plot(tvec, nanmean(ccf_low_shuf - ccf_high_shuf), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    plot(tvec, nanmean(ccf_low_shuf - ccf_high_shuf) + 1.96*nanstd(ccf_low_shuf - ccf_high_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    plot(tvec, nanmean(ccf_low_shuf - ccf_high_shuf) - 1.96*nanstd(ccf_low_shuf - ccf_high_shuf), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
    set(gca, 'FontSize', cfg.ft_size); xlabel('time (s)'); ylabel('spike count');
    title(sprintf('phase split diff %.1f-%.1f Hz', cfg.freq{iF}(1), cfg.freq{iF}(2)));
    
end % of freq loop

%save the split diff plot as a summary
figure(3)
cfg_fig.ft_size = cfg.ft_size;
SetFigure(cfg_fig, gcf)
saveas_eps([fname '_' S.label{iC}(1:end-2)], cd)
close all
fprintf(['Complete Cell #' S.label{iC} '...\n'])
