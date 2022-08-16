% Assumes that references and all LFP channels have the same sampling
% frequency
cd('E:\Dropbox (Dartmouth College)\manish_data\M325\M325-2022-08-03');
LoadExpKeys;
evs = LoadEvents([]);


%% Set variables and parameters
title_font_size = 18;
subplot_font_size = 14;
text_font_size = 10;
sta_calculation_wsize = [-0.5 0.5]; % in seconds
sta_plot_wsize = [-0.1 0.1]; % in seconds
lfp_plot_wsize = [-0.5 0.5]; % in seconds
num_snips = 10; %How many snippets of LFP to create;


%%
cfg_lfp.fc = ExpKeys.LFP_channels;
all_lfp = LoadCSC(cfg_lfp);
% Modify the all_lfp to now have referenced LFPs
% TODO: Modify headers accordingly
if strcmp(ExpKeys.isReferenceRecordedSeparately,'Yes')
    cfg_refs.fc = ExpKeys.ref_channels;
    all_refs = LoadCSC(cfg_refs);
    num_lfp = length(cfg_lfp.fc);
    num_ref = length(cfg_refs.fc);
    % Create reference subtracted versions 
    referenced_labels = cell(1,num_lfp*num_ref);
    stacked_ref_data = zeros(length(referenced_labels), ...
        length(all_refs.data));
    for iR = 1:num_ref
       stacked_ref_data(((iR-1)*num_lfp)+1:iR*num_lfp,:) = ...
           repmat(all_refs.data(iR,:),num_lfp,1);
    end
    referenced_data = repmat(all_lfp.data,num_ref, 1) - stacked_ref_data;
    for iR = 1:num_ref
        for iL = 1:num_lfp
           referenced_labels{(iR-1)*num_lfp + iL} = ...
               cat(2, extractBefore(all_lfp.label{iL},'.ncs'), '-', ...
               extractBefore(all_refs.label{iR}, '.ncs'));
        end
    end
    og_lfp = all_lfp; %Make a copy for comparison later
    all_lfp.data = referenced_data;
    all_lfp.label = referenced_labels; 
else
    all_lfp.label = cellfun(@(x) extractBefore(x,'.ncs'), ...
        all_lfp.label,'UniformOutput',false);
end


%% Generate overall PSDs
Fs = all_lfp.cfg.hdr{1}.SamplingFrequency;
wsize = round(4*Fs);  % arbitrary decision on Window Size
fig1 = figure;
for i = 1:length(all_lfp.label)
    [Pxx, F] = pwelch(all_lfp.data(i,:), rectwin(wsize), wsize/2, [], Fs);
    plot(F, 10*log10(Pxx));
    hold on
end
legend(all_lfp.label)
xlabel('Frequency in Hertz');
xlim([0,120]);
title('All LFP PSD', 'FontSize', title_font_size);
%Save as PDF and .fig
print(fig1, '-dpdf', '-fillpage', strcat(ExpKeys.subject_id,'-', ...
    ExpKeys.date,'-','All_LFP_PSD'))
% saveas(fig1, strcat(ExpKeys.subject_id,'-', ExpKeys.date,'-','All_LFP_PSD'));
if strcmp(ExpKeys.isReferenceRecordedSeparately,'Yes')
    Fs = og_lfp.cfg.hdr{1}.SamplingFrequency;
    wsize = round(4*Fs);  % arbitrary decision on Window Size
    fig2 = figure;
    for i = 1:length(og_lfp.label)
        [Pxx, F] = pwelch(og_lfp.data(i,:), rectwin(wsize), wsize/2, [], Fs);
        plot(F, 10*log10(Pxx));
        hold on
    end
    legend(og_lfp.label)
    xlabel('Frequency in Hertz');
    xlim([0,120]);
    title('OG LFP PSD', 'FontSize', title_font_size);
    %Save as PDF and .fig
    print(fig2, '-dpdf', '-fillpage', strcat(ExpKeys.subject_id,'-', ...
        ExpKeys.date,'-','OG_LFP_PSD'))
%     saveas(fig2, strcat(ExpKeys.subject_id,'-', ExpKeys.date,'-','OG_LFP_PSD'));
    
    Fs = all_refs.cfg.hdr{1}.SamplingFrequency;
    wsize = round(4*Fs);  % arbitrary decision on Window Size
    fig3 = figure;
    for i = 1:length(all_refs.label)
        [Pxx, F] = pwelch(all_refs.data(i,:), rectwin(wsize), wsize/2, [], Fs);
        plot(F, 10*log10(Pxx));
        hold on
    end
    legend(all_refs.label)
    xlabel('Frequency in Hertz');
    xlim([0,120]);
    title('Reference Channel PSD', 'FontSize', title_font_size);
    print(fig3, '-dpdf', '-fillpage', strcat(ExpKeys.subject_id,'-', ...
        ExpKeys.date,'-','REF_PSD'))
%     saveas(fig3, strcat(ExpKeys.subject_id,'-', ExpKeys.date,'-','REF_PSD'));
end

%% Generate Epoch-wise reports for each candidate LFP
close all;
cfg_filter.f = [1 500];
cfg_filter.order = 3;
cfg_filter.type = 'butter';

% Pre-stim
cfg_pre.lfp = restrict(all_lfp, iv(ExpKeys.pre_stim_times));
%Filter the lfp for plotting traces
cfg_pre.filtered_lfp = FilterLFP(cfg_filter, cfg_pre.lfp);
cfg_pre.subplot_font_size = subplot_font_size;
cfg_pre.text_font_size = text_font_size;
cfg_pre.sta_calculation_wsize = sta_calculation_wsize;
cfg_pre.sta_plot_wsize = sta_plot_wsize;
cfg_pre.lfp_plot_wsize = lfp_plot_wsize;
cfg_pre.stim_on = evs.t{strcmp(evs.label, ExpKeys.pre_trial_stim_on)};
cfg_pre.snip_idx = sort(randsample(1:length(cfg_pre.stim_on),num_snips));
cfg_pre.plot_prefix = "Pre-stim";
cfg_pre.probe_depth = ExpKeys.probeDepth;
cfg_pre.hemi = ExpKeys.hemisphere;
cfg_pre.stim_width = ExpKeys.short_stim_pulse_width;
for iL = 1:length(all_lfp.label)
    cfg_pre.idx = iL;
    this_fig = createFigure(cfg_pre);
    this_fn = strcat(ExpKeys.subject_id,'-', ExpKeys.date, ...
        '-Pre-Stim-',all_lfp.label{iL});
    print(this_fig, '-dpdf', '-fillpage', this_fn);
%     saveas(this_fig, this_fn);
    close;
end

% Trial stim
cfg_trial.lfp = restrict(all_lfp, iv(ExpKeys.stim_times));
%Filter the lfp for plotting traces
cfg_trial.filtered_lfp = FilterLFP(cfg_filter, cfg_trial.lfp);
cfg_trial.subplot_font_size = subplot_font_size;
cfg_trial.text_font_size = text_font_size;
cfg_trial.sta_calculation_wsize = sta_calculation_wsize;
cfg_trial.sta_plot_wsize = sta_plot_wsize;
cfg_trial.lfp_plot_wsize = lfp_plot_wsize;
cfg_trial.stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)};
cfg_trial.snip_idx = sort(randsample(1:length(cfg_trial.stim_on),num_snips));
cfg_trial.plot_prefix = "Trial-stim";
cfg_trial.probe_depth = ExpKeys.probeDepth;
cfg_trial.hemi = ExpKeys.hemisphere;
cfg_trial.stim_width = ExpKeys.short_stim_pulse_width;
for iL = 1:length(all_lfp.label)
    cfg_trial.idx = iL;
    this_fig = createFigure(cfg_trial);
    this_fn = strcat(ExpKeys.subject_id,'-', ExpKeys.date, ...
        '-Trial-Stim-',all_lfp.label{iL});
    print(this_fig, '-dpdf', '-fillpage', this_fn);
%     saveas(this_fig, this_fn);
    close;
end

% Post-stim
cfg_post.lfp = restrict(all_lfp, iv(ExpKeys.post_stim_times));
%Filter the lfp for plotting traces
cfg_post.filtered_lfp = FilterLFP(cfg_filter, cfg_post.lfp);
cfg_post.text_font_size = text_font_size;
cfg_post.subplot_font_size = subplot_font_size;
cfg_post.sta_calculation_wsize = sta_calculation_wsize;
cfg_post.sta_plot_wsize = sta_plot_wsize;
cfg_post.lfp_plot_wsize = lfp_plot_wsize;
cfg_post.stim_on = evs.t{strcmp(evs.label, ExpKeys.post_trial_stim_on)};
cfg_post.snip_idx = sort(randsample(1:length(cfg_post.stim_on),num_snips));
cfg_post.plot_prefix = "Post-stim";
cfg_post.probe_depth = ExpKeys.probeDepth;
cfg_post.hemi = ExpKeys.hemisphere;
cfg_post.stim_width = ExpKeys.short_stim_pulse_width;
for iL = 1:length(all_lfp.label)
    cfg_post.idx = iL;
    this_fig = createFigure(cfg_post);
    this_fn = strcat(ExpKeys.subject_id,'-', ExpKeys.date, ...
        '-Post-Stim-',all_lfp.label{iL});
    print(this_fig, '-dpdf', '-fillpage', this_fn);
%     saveas(this_fig, this_fn);
    close;
end

% Long-stim
cfg_long.lfp = restrict(all_lfp, iv(ExpKeys.long_stim_times));
%Filter the lfp for plotting traces
cfg_long.filtered_lfp = FilterLFP(cfg_filter, cfg_long.lfp);
cfg_long.text_font_size = text_font_size;
cfg_long.subplot_font_size = subplot_font_size;
cfg_long.sta_calculation_wsize = sta_calculation_wsize;
cfg_long.sta_plot_wsize = sta_plot_wsize;
cfg_long.lfp_plot_wsize = lfp_plot_wsize;
cfg_long.stim_on = evs.t{strcmp(evs.label, ExpKeys.long_stim_on)};
cfg_long.snip_idx = sort(randsample(1:length(cfg_long.stim_on),num_snips));
cfg_long.plot_prefix = "Long-stim";
cfg_long.probe_depth = ExpKeys.probeDepth;
cfg_long.hemi = ExpKeys.hemisphere;
cfg_long.stim_width = ExpKeys.long_stim_pulse_width;
for iL = 1:length(all_lfp.label)
    cfg_long.idx = iL;
    this_fig = createFigure(cfg_long);
    this_fn = strcat(ExpKeys.subject_id,'-', ExpKeys.date, ...
        '-Long-Stim-',all_lfp.label{iL});
    print(this_fig, '-dpdf', '-fillpage', this_fn);
%     saveas(this_fig, this_fn);
    close;
end

%% Helper functions
function fh = createFigure(cfg_in)
    fh = figure;
    this_lfp = cfg_in.lfp;
    filtered_lfp = cfg_in.filtered_lfp;
    idx = cfg_in.idx;
    Fs = this_lfp.cfg.hdr{1}.SamplingFrequency;
    wsize = 4*Fs;
    fsz1 = cfg_in.subplot_font_size;
    fsz2 = cfg_in.text_font_size;
    wsz1 = cfg_in.sta_calculation_wsize;
    wsz2 = cfg_in.sta_plot_wsize;
    wsz3 = cfg_in.lfp_plot_wsize;
    stim_on = cfg_in.stim_on;
    snip_idx = cfg_in.snip_idx;
    stim_width = cfg_in.stim_width;
    depth = cfg_in.probe_depth;
    hemi = cfg_in.hemi;
    nRows = ceil(((length(snip_idx))+4)/2);
    nCols = 2;
    
    % Info about session
    subplot(nRows,nCols,1)
    hold on;
    text(0, 1, sprintf("Recording Depth: %.2f mm", depth), 'FontSize', fsz2);
    text(0, 0.5, sprintf("Hemisphere: %s", hemi), 'FontSize', fsz2);
    text(0, 0, sprintf("Stimulus width: %.2f msec", stim_width*1000), 'FontSize', fsz2);
    axis off;
    
    %Plot raw LFP trace to see the absurd jumps
    subplot(nRows, nCols, 2)
    hold on
    plot(this_lfp.tvec,this_lfp.data(idx,:), 'green');
    for iS = 1:length(snip_idx)
        xline(stim_on(snip_idx(iS)), 'red');
    end 
    xticks([]);
    title(sprintf('%s LFP trace for %s',cfg_in.plot_prefix, this_lfp.label{idx}), 'FontSize', fsz1);
    
    
    % Plot PSD
    subplot(nRows,nCols,3)
    hold on
    [Pxx, F] = pwelch(this_lfp.data(idx,:), rectwin(wsize), wsize/2, [], Fs);
    plot(F, 10*log10(Pxx), 'black');
    xlabel('Frequency in Hertz')
    xlim([0 120])
    title(sprintf('%s PSD for %s',cfg_in.plot_prefix, this_lfp.label{idx}), 'FontSize', fsz1);

    %Plot STA
    subplot(nRows,nCols,4)
    hold on
    num_stim = length(stim_on);
    this_tvec = wsz1(1):1/Fs:wsz1(2); % time axis for STA
    this_on_sta = zeros(num_stim, length(this_tvec));
    filtered_snips = zeros(size(this_on_sta));
    for iEvt = 1:num_stim % for each stim ...
        on_sta_t = stim_on(iEvt)+wsz1(1);
        this_on_sta_idx = (nearest_idx3(on_sta_t,this_lfp.tvec));
        % grab LFP snippet for this window
        this_on_toAdd = this_lfp.data(idx,this_on_sta_idx:this_on_sta_idx+ ...
            length(this_tvec)-1);
        this_on_sta(iEvt,:) = this_on_toAdd';
        %also grab Filtered LFP Snippets
        filtered_toAdd = filtered_lfp.data(idx,this_on_sta_idx:this_on_sta_idx+ ...
            length(this_tvec)-1);
        filtered_snips(iEvt,:) = filtered_toAdd';
    end
    on_sta = mean(this_on_sta,1);
    plot(this_tvec, on_sta);
    xline(0, 'red');
    if stim_width > 0.001
        xline(stim_width, 'red');
    end
    xlim(wsz2);
    xlabel('Time in seconds');
    title('STA aligned to STIM ON', 'FontSize', fsz1);

    %Plot num_snip filtered LFP snippets
    for iN = 1:length(snip_idx)
        subplot(nRows,nCols,(5+iN-1))
        hold on;
%         plot(this_tvec, this_on_sta(snip_idx(iN),:));
        plot(this_tvec, filtered_snips(snip_idx(iN),:), 'green');
        xline(0, 'red');
        if stim_width > 0.001
            xline(stim_width, 'red');
        end
    end

end

%%

