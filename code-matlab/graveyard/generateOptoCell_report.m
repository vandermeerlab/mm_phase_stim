%% Assumes that spike-sorting has been done

top_dir = 'data\';
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
    LoadExpKeys;
    evs = LoadEvents([]);
    cfg_spk = [];
    cfg_spk.fc = ExpKeys.goodCell;
    if isempty(cfg_spk.fc)
        return
    end
    % cfg_spk.min_cluster_quality = 3;
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg_spk.getRatings = 1;
        cfg_spk.uint = '64';
    end
    S = LoadSpikes(cfg_spk);
    
    %% Set variables
    text_fs = 11; % Text font size
    small_label_fs = 12; 
    big_label_fs = 16; 
    max_delay = 0.01; % sec (window for the first response since stimulus)
    
    %% Remove spikes during short stim times
    if contains(ExpKeys.light_source, 'LASER')
        start_delay = 0.0011;
        stop_delay = 0.0012; %0.0022 is the spike width in the case of 1 msec laser pulse
    else
        start_delay = 0;
        stop_delay = 0;
    end

    pre_stim_on = evs.t{strcmp(evs.label, ExpKeys.pre_trial_stim_on)} + start_delay;
    if ~isempty(pre_stim_on) && ~isempty(ExpKeys.pre_stim_times)
        pre_stim_on = pre_stim_on(pre_stim_on >= ExpKeys.pre_stim_times(1) & ...
                                  pre_stim_on <= ExpKeys.pre_stim_times(2));
    else
        pre_stim_on = [];
    end

    stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
    if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
        stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                          stim_on <= ExpKeys.stim_times(2));
    else
        stim_on = [];
    end

    post_stim_on = evs.t{strcmp(evs.label, ExpKeys.post_trial_stim_on)} + start_delay;
    if ~isempty(post_stim_on) && ~isempty(ExpKeys.post_stim_times)
        post_stim_on = post_stim_on(post_stim_on >= ExpKeys.post_stim_times(1) & ...
                                    post_stim_on <= ExpKeys.post_stim_times(2));
    else
        post_stim_on = [];
    end

    all_short_stim = [pre_stim_on; stim_on; post_stim_on];
    
    all_long_stim = [evs.t{strcmp(evs.label, ExpKeys.long_stim_on)}];
    
    to_be_removed = [all_short_stim,  all_short_stim + ExpKeys.short_stim_pulse_width + stop_delay];
    
    
    clean_iv = InvertIV(iv(to_be_removed), ExpKeys.recording_times(1), ...
        ExpKeys.recording_times(2));
    restricted_S = restrict(S, clean_iv);
    
    
    %% Set variables and parameters
    % snippet for autocorrelation
    for iC = 1:length(restricted_S.label)
        fn_prefix = extractBefore(restricted_S.label{iC}, '.t');
        fn_prefix = strrep(fn_prefix, '_', '-');
        
        %Extract AD2BitVoltConversioFactor from corresponding .ntt file
        temp_idx = find(fn_prefix == '-');
        ntt_fn = strcat(extractBefore(fn_prefix, temp_idx(end)), '.ntt');
        if isfile(ntt_fn)
            [~, ~, ~, ~, ~, hdr] = Nlx2MatSpike(ntt_fn, [1 1 1 1 1], 1, 1, []);
            temp_idx = find(contains(hdr, 'ADBitVolts'));
            if ~ isempty(temp_idx)
                hdr = hdr{temp_idx};
                hdr = split(hdr, ' ');
                conv_factors = cellfun(@str2double, hdr(2:end));
                %Convert to microvolts
                conv_factors = conv_factors*1e6;
            else
               fprintf("Something wrong with the header of %s, Skipping\n", ntt_fn);
               continue;
            end
        else
            fprintf("%s not found in the directory, Skipping\n", ntt_fn);
            continue;
        end
    
        this_fig = figure('WindowState', 'maximized');
        
        ax = subplot(5,4,1);
        % Load ClusterQual
        temp_cq = load(strcat(fn_prefix, '-ClusterQual.mat'));
        text(0, 1, fn_prefix, 'FontSize', text_fs, 'Interpreter', 'none');
        text(0, 0.8, sprintf("LRatio: %.2f ", temp_cq.CluSep.Lratio), ...
            'FontSize', text_fs);
        text(0, 0.6, sprintf("Isolation Distance: %.2f ", ...
            temp_cq.CluSep.IsolationDist), 'FontSize', text_fs);
        text(0, 0.2, sprintf("Mean Firing Rate: %.2f Hz", ...
            length(restricted_S.t{iC})/sum(clean_iv.tend - clean_iv.tstart)), ...
            'FontSize', text_fs); 
        if ~strcmp(ExpKeys.experimenter, 'EC')
            text(0, 0.4, sprintf("Cluster Rating: %d", ...
                restricted_S.usr.rating(iC)), 'FontSize', text_fs);
%             text(0, 0, sprintf("Stimulus amplitude: %.2f mW", ...
%                 interp1(ExpKeys.dial_values_calib, ExpKeys.output_power_calib, ...
%                 ExpKeys.dial_value, 'spline')), 'FontSize', text_fs);
            text(0, 0, sprintf("Stimulus amplitude: %.2f mW", ...
                ExpKeys.light_power), 'FontSize', text_fs);
        end
        axis off;
        
        % Plot autocorrelation
        ax = subplot(5,4,[5,6]);
        cfg_xc.smooth = 1;
        cfg_xc.max_t = 1;
        [test_xc, test_t] = ccf(cfg_xc,restricted_S.t{iC},restricted_S.t{iC});
        plot(test_t, test_xc);
        xlim([test_t(find(test_t > 0, 1)), test_t(end)])
        xlabel('Time (sec)')
        title('Autocorrelation')
        ax.XAxis.FontSize = small_label_fs;
        ax.YAxis.FontSize = small_label_fs;
        ax.XAxis.Label.FontSize = small_label_fs;
        ax.YAxis.Label.FontSize = small_label_fs;
        ax.Title.FontSize = big_label_fs;

        
        % Plot logISI
        ax = subplot(5,4,[7,8]);
        test_isi = diff(restricted_S.t{iC});
        histogram(log10(test_isi), -3.75:0.5:2.75);
        xline(log10(0.002), 'r');
        this_str = sprintf(', %d/%d ISIs less than 2 msec', ...
            sum(test_isi<0.002), length(test_isi));
        title(strcat('LogISI', this_str))
        ax.XAxis.FontSize = small_label_fs;
        ax.YAxis.FontSize = small_label_fs;
        ax.XAxis.Label.FontSize = small_label_fs;
        ax.YAxis.Label.FontSize = small_label_fs;
        ax.Title.FontSize = big_label_fs;
        
        % Load average waveformfile and plot it
        ax = subplot(5,4,2);
        temp_wv = load(strcat(fn_prefix,'-wv.mat'));
        %Convert_Values to microvolts
        mWV = temp_wv.mWV .* conv_factors';
        plot(mWV);
        ylabel('micro-Volts');
        title('Average SpikeWaveform')
        ax.XAxis.FontSize = small_label_fs;
        ax.YAxis.FontSize = small_label_fs;
        ax.XAxis.Label.FontSize = small_label_fs;
        ax.YAxis.Label.FontSize = small_label_fs;
        ax.Title.FontSize = big_label_fs;
        
        % Extract Spike amplitudes from coresponding .awv file
        awv_fn = strcat(fn_prefix, '-awv.mat');
        if isfile(awv_fn)
            ax = subplot(5,4,[3,4]);
            temp_awv = load(awv_fn);
            temp_awv = temp_awv.WV;
            max_awv = max(temp_awv,[],3);
            max_awv = max_awv .* conv_factors';
            plot(S.t{iC}, max_awv, '.', 'MarkerSize', 10)
            % Uncomment the following section to figure out overlay trial_stim
            hold on;
            trial_stim = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)};
            % Overlay horizontal lines for every 50 stim
            topY = max(max(max(max_awv))) - 50;
            for iStim = 1:50:length(trial_stim)
                xline(trial_stim(iStim), '--black')
                text(trial_stim(iStim), topY, num2str(iStim), 'Rotation', 90, 'FontSize', text_fs)
            end
            ylabel('micro-Volts');
            ax.XAxis.FontSize = small_label_fs;
            ax.YAxis.FontSize = small_label_fs;
            ax.XAxis.Label.FontSize = small_label_fs;
            ax.YAxis.Label.FontSize = small_label_fs;
            ax.Title.FontSize = big_label_fs;
        else
           fprintf("All Waveforms file missing for %s, Skipping\n", fn_prefix);
           continue
        end
        
        this_cell = SelectTS([], S, iC);
        restricted_cell = SelectTS([], restricted_S, iC);
        
        % Pre-stim raster
        if ~isempty(ExpKeys.pre_stim_times)
            ax = subplot(5,4,[9,13]);
            this_on_events = pre_stim_on;
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm([], ...
                this_cell, this_on_events, '', 0.5);
            hold on
            plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
            plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            plot([ExpKeys.short_stim_pulse_width+stop_delay ExpKeys.short_stim_pulse_width+stop_delay], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            ylabel('Pre-Stim #')
            ylim([1 length(this_on_events)])
            xlim([-0.01 0.01]);
            xlabel("Time (sec)")
            ax.XAxis.FontSize = small_label_fs;
            ax.YAxis.FontSize = small_label_fs;
            ax.XAxis.Label.FontSize = small_label_fs;
            ax.YAxis.Label.FontSize = big_label_fs;
            ax.YAxis.Label.FontWeight = 'bold';
        end

        % Pre-stim response hist
        latency = nan(size(this_on_events));
        latency_wo_stim = nan(size(this_on_events));
        for iStim = 1:length(latency_wo_stim)
            st = restrict(this_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
            st_wo_stim = restrict(restricted_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
            if ~isempty(st.t{1})
                latency(iStim) = st.t{1}(1) - this_on_events(iStim);
            end
            if ~isempty(st_wo_stim.t{1})
                latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                assert(latency_wo_stim(iStim) >= ExpKeys.short_stim_pulse_width, ...
                    sprintf('Latency for pre-trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
            end
        end
        ax = subplot(5,4,17);
        hold on;
        histogram(latency*1000, [0:0.5:10], 'FaceAlpha', 0.5, 'FaceColor', 'Green')
        histogram(latency_wo_stim*1000, [0:0.5:10], 'FaceAlpha', 0.5, 'FaceColor', 'Red')
        leg = legend({sprintf('\\color{green}Prob: %.2f', sum(~isnan(latency))/length(latency)), ...
            sprintf('\\color{red}Prob w/o stim: %.2f', sum(~isnan(latency_wo_stim))/length(latency_wo_stim))});
        leg.FontSize = text_fs;
        leg.Location = 'best';
        xlabel('Delay since stim (msec)')
        ylabel('Count')
        leg.Box = 'off';
        ax.XAxis.FontSize = small_label_fs;
        ax.YAxis.FontSize = small_label_fs;
        ax.XAxis.Label.FontSize = small_label_fs;
        ax.YAxis.Label.FontSize = big_label_fs;
        ax.YAxis.Label.FontWeight = 'bold';
        
%         ax.Title.Interpreter = 'tex';
%         ax.Title.String = sprintf("{\\color{green}prob: %.2f}, {\\color{red}prob w/o stim: %.2f}" , ...
%             sum(~isnan(latency))/length(latency), sum(~isnan(latency_wo_stim))/length(latency_wo_stim));
% 
%         ax.Title.FontSize = small_label_fs;
       
        % Trial-stim raster
        if ~isempty(ExpKeys.stim_times)
            ax = subplot(5,4,[10,14]);
            this_on_events = stim_on;
            if ~isempty(ExpKeys.goodTrials) % Only keep the good trials
                this_on_events = this_on_events(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2));
            end
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm([], ...
                this_cell, this_on_events, '', 0.5);
            hold on
            plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
            plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            plot([ExpKeys.short_stim_pulse_width+stop_delay ExpKeys.short_stim_pulse_width+stop_delay], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            ylabel(' Trial Stim #');
            ylim([1 length(this_on_events)])
            xlim([-0.01 0.01]);
            xlabel("Time (sec)", 'FontSize', small_label_fs)
            ax.XAxis.FontSize = small_label_fs;
            ax.YAxis.FontSize = small_label_fs;
            ax.XAxis.Label.FontSize = small_label_fs;
            ax.YAxis.Label.FontSize = big_label_fs;
            ax.YAxis.Label.FontWeight = 'bold';
        end

        % Trial-stim response hist
        latency = nan(size(this_on_events));
        latency_wo_stim = nan(size(this_on_events));
        for iStim = 1:length(latency_wo_stim)
            st = restrict(this_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
            st_wo_stim = restrict(restricted_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
            if ~isempty(st.t{1})
                latency(iStim) = st.t{1}(1) - this_on_events(iStim);
            end
            if ~isempty(st_wo_stim.t{1})
                latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                assert(latency_wo_stim(iStim) >= ExpKeys.short_stim_pulse_width, ...
                    sprintf('Latency for Trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
            end
        end
        ax = subplot(5,4,18);
        hold on;
        histogram(latency*1000, [0:0.5:10], 'FaceAlpha', 0.5, 'FaceColor', 'Green')
        histogram(latency_wo_stim*1000, [0:0.5:10], 'FaceAlpha', 0.5, 'FaceColor', 'Red')
        leg = legend({sprintf('\\color{green}Prob: %.2f', sum(~isnan(latency))/length(latency)), ...
            sprintf('\\color{red}Prob w/o stim: %.2f', sum(~isnan(latency_wo_stim))/length(latency_wo_stim))});
        leg.FontSize = text_fs;
        leg.Location = 'best';
        xlabel('Delay since stim (msec)')
        ylabel('Count')
        leg.Box = 'off';
        ax.XAxis.FontSize = small_label_fs;
        ax.YAxis.FontSize = small_label_fs;
        ax.XAxis.Label.FontSize = small_label_fs;
        ax.YAxis.Label.FontSize = big_label_fs;
        ax.YAxis.Label.FontWeight = 'bold';

        % Post-stim raster
        if ~isempty(ExpKeys.post_stim_times)
            ax = subplot(5,4,[11,15]);
            this_on_events = post_stim_on;
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm([], ...
                this_cell, this_on_events, '', 0.5);
            hold on
            plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
            plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            plot([ExpKeys.short_stim_pulse_width+stop_delay ExpKeys.short_stim_pulse_width+stop_delay], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            ylabel('Post Stim #');
            ylim([1 length(this_on_events)])
            xlim([-0.01 0.01]);
            xlabel("Time (sec)", 'FontSize', small_label_fs)
            ax.XAxis.FontSize = small_label_fs;
            ax.YAxis.FontSize = small_label_fs;
            ax.XAxis.Label.FontSize = small_label_fs;
            ax.YAxis.Label.FontSize = big_label_fs;
            ax.YAxis.Label.FontWeight = 'bold';
        end

        % Post-stim response hist
        latency = nan(size(this_on_events));
        latency_wo_stim = nan(size(this_on_events));
        for iStim = 1:length(latency_wo_stim)
            st = restrict(this_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
            st_wo_stim = restrict(restricted_cell, iv(this_on_events(iStim), this_on_events(iStim)+max_delay));
            if ~isempty(st.t{1})
                latency(iStim) = st.t{1}(1) - this_on_events(iStim);
            end
            if ~isempty(st_wo_stim.t{1})
                latency_wo_stim(iStim) = st_wo_stim.t{1}(1) - this_on_events(iStim);
                assert(latency_wo_stim(iStim) >= ExpKeys.short_stim_pulse_width, ...
                    sprintf('Latency for post-trial stim #%d is %.2f msec', iStim, latency_wo_stim(iStim)*1000));
            end
        end
        ax = subplot(5,4,19);
        hold on;
        histogram(latency*1000, [0:0.5:10], 'FaceAlpha', 0.5, 'FaceColor', 'Green')
        histogram(latency_wo_stim*1000, [0:0.5:10], 'FaceAlpha', 0.5, 'FaceColor', 'Red')
        leg = legend({sprintf('\\color{green}Prob: %.2f', sum(~isnan(latency))/length(latency)), ...
            sprintf('\\color{red}Prob w/o stim: %.2f', sum(~isnan(latency_wo_stim))/length(latency_wo_stim))});
        leg.FontSize = text_fs;
        leg.Location = 'best';
        xlabel('Delay since stim (msec)')
        ylabel('Count')
        leg.Box = 'off';
        ax.XAxis.FontSize = small_label_fs;
        ax.YAxis.FontSize = small_label_fs;
        ax.XAxis.Label.FontSize = small_label_fs;
        ax.YAxis.Label.FontSize = big_label_fs;
        ax.YAxis.Label.FontWeight = 'bold';
        
        % Long-stim raster
        if ~isempty(ExpKeys.long_stim_times)
            ax = subplot(5,4,[12,16]);
            this_on_events = evs.t{strcmp(evs.label, ExpKeys.long_stim_on)} + start_delay;
            this_on_events = this_on_events(this_on_events >= ExpKeys.long_stim_times(1) & ...
                this_on_events <= ExpKeys.long_stim_times(2));
            %Take only the first 25 long stim if special stim present
            if strcmp(ExpKeys.hasSpecialStim, 'Yes')
                this_on_events = this_on_events(1:25); 
            end
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm([], ...
                this_cell, this_on_events, '', 0.5);
            hold on
            plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
            plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            plot([ExpKeys.long_stim_pulse_width ExpKeys.long_stim_pulse_width], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            ylabel('Long Stim #');
            ylim([1 length(this_on_events)])
            xlim([-0.1 0.2]);
            xlabel("Time (sec)", 'FontSize', small_label_fs)
            ax.XAxis.FontSize = small_label_fs;
            ax.YAxis.FontSize = small_label_fs;
            ax.XAxis.Label.FontSize = small_label_fs;
            ax.YAxis.Label.FontSize = big_label_fs;
            ax.YAxis.Label.FontWeight = 'bold';
        end
        savefig(this_fig, strcat(fn_prefix,'-CellReport'));
%         print(this_fig, '-dpdf', '-fillpage', strcat(fn_prefix,'-CellReport'));
%         print(this_fig, '-dpng',  strcat('data\', fn_prefix, '-CellReport'));
        close;
    end
end
