% Assumes that spike-sorting has been done
cd('E:\Dropbox (Dartmouth College)\manish_data\M320\M320-2022-05-26');
LoadExpKeys;
evs = LoadEvents([]);
cfg_spk.min_cluster_quality = 3;
cfg_spk.uint = '64';
S = LoadSpikes(cfg_spk);

%% Set variables
fsz1 = 7; % Text font size
fsz2 = 8; % 

%% Remove spikes during stim times
all_short_stim = [evs.t{strcmp(evs.label, ExpKeys.pre_trial_stim_on)}; ...
    evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)}; ...
    evs.t{strcmp(evs.label, ExpKeys.post_trial_stim_on)}];

all_long_stim = [evs.t{strcmp(evs.label, ExpKeys.long_stim_on)}];
to_be_removed = [all_short_stim,  all_short_stim + ExpKeys.short_stim_pulse_width; ...
    all_long_stim, all_long_stim + ExpKeys.long_stim_pulse_width];


clean_iv = InvertIV(iv(to_be_removed), ExpKeys.recording_times(1), ...
    ExpKeys.recording_times(2));
restricted_S = restrict(S, clean_iv);


%% Set variables and parameters
% snippet for autocorrelation
for iC = 1:length(restricted_S.label)
    fn_prefix = extractBefore(restricted_S.label{iC}, '.t');
    temp_idx = find(fn_prefix == '_');
    fn_prefix(temp_idx) = '-';
    
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
    subplot(4,4,1)
    % Load ClusterQual
    temp_cq = load(strcat(fn_prefix, '-ClusterQual.mat'));
    text(0, 1, fn_prefix, 'FontSize', fsz1, 'Interpreter', 'none');
    text(0, 0.8, sprintf("LRatio: %.2f ", temp_cq.CluSep.Lratio), ...
        'FontSize', fsz1);
    text(0, 0.6, sprintf("Isolation Distance: %.2f ", ...
        temp_cq.CluSep.IsolationDist), 'FontSize', fsz1);
    text(0, 0.4, sprintf("Cluster Rating: %d", ...
        restricted_S.usr.rating(iC)), 'FontSize', fsz1);
    text(0, 0.2, sprintf("Mean Firing Rate: %.2f Hz", ...
        length(restricted_S.t{iC})/sum(clean_iv.tend - clean_iv.tstart)), ...
        'FontSize', fsz1); 
    text(0, 0, sprintf("Stimulus amplitude: %.2f mW", ...
        interp1(ExpKeys.dial_values_calib, ExpKeys.output_power_calib, ...
        ExpKeys.dial_value, 'spline')), 'FontSize', fsz1);
    axis off;
    
    % Plot autocorrelation
    subplot(4,4,2)
    cfg_xc.smooth = 1;
    [test_xc, test_t] = ccf(cfg_xc,restricted_S.t{iC},restricted_S.t{iC});
    plot(test_t, test_xc);
    xlim([test_t(find(test_t > 0, 1)), test_t(end)])
    xlabel('Time in Seconds', 'FontSize', fsz2)
    title('Autocorrelation', 'FontSize', fsz2)
    
    % Plot logISI
    subplot(4,4,3)
    test_isi = diff(restricted_S.t{iC});
    histogram(log10(test_isi), -3.75:0.5:2.75);
    xline(log10(0.002), 'r');
    this_str = sprintf(', %d/%d ISIs less than 2 msec', ...
        sum(test_isi<0.002), length(test_isi));
    title(strcat('LogISI', this_str), 'FontSize', fsz2)
    
    % Load average waveformfile and plot it
    subplot(3,4,4)
    temp_wv = load(strcat(fn_prefix,'-wv.mat'));
    %Convert_Values to microvolts
    mWV = temp_wv.mWV .* conv_factors';
    plot(mWV);
    ylabel('micro-Volts');
    title('Average SpikeWaveform', 'FontSize', fsz2)
    
    %Extract Spike amplitudes from coresponding .awv file
    awv_fn = strcat(fn_prefix, '-awv.mat');
    if isfile(awv_fn)
        subplot(4,4,[5,6,7,8])
        temp_awv = load(awv_fn);
        temp_awv = temp_awv.WV;
        max_awv = max(temp_awv,[],3);
        max_awv = max_awv .* conv_factors';
        plot(S.t{iC}, max_awv, '.', 'MarkerSize', 10)
        ylabel('micro-Volts');
        dummy2 = 0;
    else
       fprintf("All Waveforms file missing for %s, Skipping\n", fn_prefix);
       continue
    end
    
    this_cell = SelectTS([], S, iC);
    
    %Pre-stim raster
    subplot(4,4,[9,13])
    this_on_events = evs.t{strcmp(evs.label, ExpKeys.pre_trial_stim_on)};
    [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm([], ...
        this_cell, this_on_events, '', 0.5);
    hold on
    plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
    plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
    plot([ExpKeys.short_stim_pulse_width ExpKeys.short_stim_pulse_width], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
    ylabel('Stim #');
    ylim([1 length(this_on_events)])
    xlim([-0.01 0.01]);
    xlabel("Time in seconds", 'FontSize', fsz2)
    title('Pre-stim')
   
    %Trial-stim raster
    subplot(4,4,[10,14])
    this_on_events = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)};
    [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm([], ...
        this_cell, this_on_events, '', 0.5);
    hold on
    plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
    plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
    plot([ExpKeys.short_stim_pulse_width ExpKeys.short_stim_pulse_width], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
    ylabel('Stim #');
    ylim([1 length(this_on_events)])
    xlim([-0.01 0.01]);
    xlabel("Time in seconds", 'FontSize', fsz2)
    title('Trial-stim')
    
    %Post-stim raster
    subplot(4,4,[11,15])
    this_on_events = evs.t{strcmp(evs.label, ExpKeys.post_trial_stim_on)};
    [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm([], ...
        this_cell, this_on_events, '', 0.5);
    hold on
    plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
    plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
    plot([ExpKeys.short_stim_pulse_width ExpKeys.short_stim_pulse_width], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
    ylabel('Stim #');
    ylim([1 length(this_on_events)])
    xlim([-0.01 0.01]);
    xlabel("Time in seconds", 'FontSize', fsz2)
    title('Post-stim')
    
    %Long-stim raster
    subplot(4,4,[12,16])
    this_on_events = evs.t{strcmp(evs.label, ExpKeys.long_stim_on)};
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
    ylabel('Stim #');
    ylim([1 length(this_on_events)])
    xlim([-0.1 0.1]);
    xlabel("Time in seconds", 'FontSize', fsz2)
    title('Long-stim')
    
    print(this_fig, '-dpdf', '-fillpage', strcat(fn_prefix,'-CellReport'));
    close;
end



