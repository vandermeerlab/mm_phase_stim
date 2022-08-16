% Assumes that spike-sorting has been done
cd('E:\Dropbox (Dartmouth College)\manish_data\M321\M321-2022-07-09');
LoadExpKeys;
evs = LoadEvents([]);
cfg_spk.min_cluster_quality = 3;
cfg_spk.uint = '64';
S = LoadSpikes(cfg_spk);

%% Set variables
fsz1 = 8; % Text font size
fsz2 = 10; % 

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
    

    this_fig = figure('WindowState', 'maximized');
    subplot(3,4,1)
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
    subplot(3,4,2)
    cfg_xc.smooth = 1;
    [test_xc, test_t] = ccf(cfg_xc,restricted_S.t{iC},restricted_S.t{iC});
    plot(test_t, test_xc);
    xlabel('Time in Seconds', 'FontSize', fsz2)
    title('Autocorrelation', 'FontSize', fsz2)
    
    % Plot logISI
    subplot(3,4,3)
    test_isi = diff(restricted_S.t{iC});
    histogram(log10(test_isi), -3.75:0.5:2.75);
    xline(-3, 'r');
    title('LogISI', 'FontSize', fsz2)
    
    % Load average waveformfile and plot it
    subplot(3,4,4)
    temp_wv = load(strcat(fn_prefix,'-wv.mat'));
    plot(temp_wv.mWV);
    title('Average SpikeWaveform', 'FontSize', fsz2)
    
    this_cell = SelectTS([], S, iC);
    
    %Pre-stim raster
    subplot(3,4,[5,9])
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
    subplot(3,4,[6,10])
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
    subplot(3,4,[7,11])
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
    subplot(3,4,[8,12])
    this_on_events = evs.t{strcmp(evs.label, ExpKeys.long_stim_on)};
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



