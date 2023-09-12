%% Assumes that good LFPs have been picked out

top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
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

%%
function doStuff
    % Declaring variables
    % Setting up parameters
    fbands = {[2 5], [6 10],[30 55]};
    c_list = {'cyan', 'red','magenta', 'green'};

    LoadExpKeys;
    evs = LoadEvents([]);
    cfg_spk = [];
    cfg_spk.fc = ExpKeys.goodCell;
    if isempty(cfg_spk.fc)
        return
    end
%     % For debugging
%     if ~ExpKeys.has60Hz
%         return
%     end
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
    % Keep the pre-stim baseline csc for obtaining optimal parameters for phase estimation
    eval_csc = restrict(csc, iv(ExpKeys.pre_baseline_times));
    csc = restrict(csc, iv(ExpKeys.stim_times));
    
    % Uncomment below if things are too slow;
%     if csc.cfg.hdr{1}.SamplingFrequency > 30000
%         cfg = []; cfg.decimateFactor = 12;
%         csc = decimate_tsd(cfg, csc);
%     end
    
    
    fig = figure('WindowState', 'maximized');
%     ax = subplot(4, 12, [1 2 3 4 13 14 15 16 25 26 27 28]);
%     % Plot PSD of TrialCSC
%     Fs = 1/median(diff(csc.tvec));
%     wsize = pow2(floor(log2(4*Fs)));  % arbitrary decision on Window Size
%     [Pxx, F] = pwelch(csc.data, hanning(wsize), wsize/2, [], Fs);
%     %Restricting Frequency to >0 and <200 Hz for FOOOF to work
%     f_idx = find(F>0 & F<120);
%     Pxx = Pxx(f_idx); F = F(f_idx);
%     P = 10*log10(Pxx); % Converting into decibels
%     plot(F, P, 'black');
%     hold on;
%     xlabel('Frequency (Hz)');
%     
%     % Testing FOOOF
%     F_fooof = F';
%     P_fooof = (Pxx(1:length(F)))';
%     reshaped_P = reshape(P_fooof,1 ,1, length(P_fooof));
%     % All these parameters are borrowed from "process_fooof.m"
%     opt.freq_range = [F(1) F(end)];
%     opt.power_line = '60';
%     opt.peak_width_limits = [0.5,12];
%     opt.max_peaks = 5;
%     opt.min_peak_height = 0.3;
%     opt.aperiodic_mode = 'knee'; %Check with 'fixed' first
%     opt.peak_threshold = 2;
%     opt.return_spectrum = 1;
%     opt.border_threshold = 1;
%     opt.peak_type = 'best'; %There is an error in documenation where it says 'both'
%     opt.proximity_threshold = 2;
%     opt.guess_weight = 'none';
%     opt.thresh_after = 1;
%     opt.sort_type = 'param';
%     opt.sort_param = 'frequency';
%     opt.sort_bands = {{'delta'}, {'2', '5'}; {'theta'}, {'6', '10'}; {'beta'},{'12', '30'}
%                       {'gamma1'}, {'30',' 55'};{'gamma2'}, {'65','90'}};
%     [fs, fg] = process_fooof('FOOOF_matlab', reshaped_P, F_fooof, opt, 1);
%     powspctrm_f = cat(1, fg.ap_fit);
%     for k = 1:size(powspctrm_f,1)
%         aperiodic_P(k,:) = interp1(fs, powspctrm_f(k,:), F, 'linear', nan);
%     end
%     aperiodic_P = 10*log10(aperiodic_P); % Converting to decibels
%     plot(F, aperiodic_P, '--black')
%     for iF = 1:length(fbands)
%         % Fill is tricky, change this later if need be
%         f_idx = find(round(F) >= fbands{iF}(1) & round(F) <= fbands{iF}(2));
% %         curve2 = P(f_idx)';
% %         curve1 = aperiodic_P(f_idx);
% %         x2 = [F(f_idx)', fliplr(F(f_idx)')];
% %         inBetween = [curve1, curve2];
% %         fill(x2, inBetween,c_list{iF})
%         area(F(f_idx), P(f_idx), 'FaceColor', c_list{iF}, 'FaceAlpha', 0.5, 'BaseValue', min(P))                    
%     end
%     xlim([0 120])

%     if ExpKeys.has60Hz
%         d = designfilt('bandstopiir','FilterOrder',2, ...
%             'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5, ...
%             'DesignMethod','butter','SampleRate',Fs);
%         filt_csc = filtfilt(d, csc.data);
%         [Pyy, ~] = pwelch(filt_csc, hanning(wsize), wsize/2, [], Fs);
%         plot(F, 10*log10(Pyy), 'blue');
%         legend({'Raw CSC','Notched CSC'})
%     end
    title('PSD');
    ax.FontSize = 12;
    %%
    % The filtering is a time taking step and is best done outside the optimization step
    winRange = [0.5 1.5]; % in seconds
    winInc = 0.1; % in seconds
    nSamples = 1000;
    od = [];
    od.win = cell(1,length(fbands));
    od.fval = cell(1,length(fbands));
    od.time_taken = zeros(1,length(fbands));
    od.best_win = zeros(1,length(fbands));
    ax0 = subplot(4, 12, [37 38 39]);
    hold on
    title('Optimization Results')
    ylabel('Ang. var. diff')
    xlabel('Window Size (sec)')
    xlim([winRange(1)-0.1, winRange(2)+0.1])
    for iB = 1:length(fbands)
        cfg_filt.type = 'fdesign'; 
        cfg_filt.f  = fbands{iB};
        filt_lfp = FilterLFP(cfg_filt, eval_csc);
        filt_phase = angle(hilbert(filt_lfp.data));
        opt_run = get_optimal_echt(eval_csc, filt_phase, fbands{iB}, Fs, nSamples, winRange, winInc);
        od.win{iB} = opt_run.wsz;
        od.fval{iB} = opt_run.fval;
        od.time_taken(iB) = opt_run.time_taken;
        [~, bidx] = min(od.fval{iB});
        od.best_win(iB) = od.win{iB}(bidx);
        plot(ax0, od.win{iB}, od.fval{iB}, 'Color', c_list{iB})
        xline(ax0, od.best_win(iB), 'Color', c_list{iB});
        ax0.FontSize = 12;
        win_length = od.best_win(iB);
        % Choose nEnds in a way such that the smallest sample is win_length long
        min_start = ceil(win_length*Fs);
        nEnds = randi(length(eval_csc.data) - min_start, nSamples, 1) + min_start;
        nStarts = nearest_idx3(eval_csc.tvec(nEnds) - win_length, eval_csc.tvec);
        estimated_phase = zeros(1,nSamples);
        true_phase = filt_phase(nEnds);
        for iS = 1:nSamples
            this_echt = echt(eval_csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
            this_phase = angle(this_echt);
            estimated_phase(iS) = this_phase(end); % The last sample's phase
        end
        ax1 = subplot(4,12, [(12*(iB-1))+5:(12*(iB-1))+8]);
        scatter(estimated_phase, true_phase, c_list{iB});
        hold on;
        plot([-pi pi], [-pi pi], 'k'); 
        
        axis tight; grid on; set(gca, 'XTick', -pi:pi/2:pi, 'YTick', -pi:pi/2:pi, 'XTickLabel', l, 'YTickLabel', l);
        title(sprintf("%d Hz - %d Hz", fbands{iB}(1), fbands{iB}(2))); xlabel('estimated phase'); ylabel('true phase');
        ax1.FontSize = 12;
        ax2 = subplot(4,12, (12*(iB-1))+9);
        histogram(true_phase, -pi:2*pi/5:pi, 'FaceColor', 'blue');
        title('True phases')
        ax2.FontSize = 12;
        ax3 = subplot(4,12, (12*(iB-1))+10);
        histogram(estimated_phase, -pi:2*pi/5:pi, 'FaceColor', [0.8500 0.3250 0.0980]);
        title('Causal Phases')
        ax3.FontSize = 12;
    end  
    % Plot phase distribution of stim before and after rejecting plots
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
    ISI = [100 diff(stim_on)'];
    causal_phase = zeros(length(fbands), length(stim_on));
    nEnds = nearest_idx3(stim_on, csc.tvec);

    for iB = 1:length(fbands)
        win_length  = od.best_win(iB);
        nStarts = nearest_idx3(stim_on - win_length, csc.tvec);
        for iS = 1:length(stim_on)
            this_echt = echt(csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
            this_phase = angle(this_echt);
            causal_phase(iB,iS) = this_phase(end); % The last sample's phase
        end
        ax4 = subplot(4,12, (12*(iB-1))+11);
        histogram(causal_phase(iB,:), -pi:2*pi/5:pi, 'FaceColor', 'blue');
        title('All Trials');
        ax4.FontSize = 12;
        
        keep = ISI >= od.best_win(iB);
        ax5 = subplot(4,12, (12*(iB-1))+12);
        histogram(causal_phase(iB,keep), -pi:2*pi/5:pi, 'FaceColor',  [0.8500 0.3250 0.0980]);
        title(sprintf('Trials > ITI: %.1f sec',od.best_win(iB)));
        ax5.FontSize = 12;
    end
    od.causal_phase = causal_phase;
    % Assume you are in the correct folder
    save('good_lfp_phases','od'); % should add option to save in specified output dir
    fn_prefix = split(pwd, '\');
    fn_prefix = fn_prefix{end};
    print(fig, '-dpng',  strcat('E:\Dropbox (Dartmouth College)\AnalysisResults\phase_stim_results\', fn_prefix, '-LFP-Report'));
    close;
end

%% Function to return optimal parameters for echt
function [opt_run] = get_optimal_echt(csc, filt_phase, fband, Fs, nSamples, winRange, winInc)
    all_win = winRange(1):winInc:winRange(2);
    if all_win(end) ~= winRange(2)
        all_win = [all_win, winRange(2)];
    end
    opt_run = [];
    opt_run.wsz = all_win;
    opt_run.fval = zeros(size(all_win));
    tic;
    for iW = 1:length(all_win)
        wsz = all_win(iW);
        min_start = ceil(wsz*Fs);
        nEnds = randi(length(csc.data) - min_start, nSamples, 1) + min_start;
        nStarts = nearest_idx3(csc.tvec(nEnds) - wsz, csc.tvec);  
        true_phase = filt_phase(nEnds);
        output_phase = zeros(1, nSamples);
        for iS = 1:nSamples
            this_echt = echt(csc.data(nStarts(iS):nEnds(iS)), fband(1), fband(2), Fs);
            this_phase = angle(this_echt);
            output_phase(iS) = this_phase(end); % The last sample's phase
        end
        % Using angular variance of difference as the cost function
        opt_run.fval(iW) = 1 - abs(mean(exp(1i*true_phase')./exp(1i*output_phase')));
    end
    opt_run.time_taken = toc;
end