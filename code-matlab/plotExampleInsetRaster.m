%% Script to generate part of Figure 4
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        this_label = 'M019-2019-04-14-TT05_1.t';
        doStuff(this_label)
    end

end
%%
function doStuff(label)
    LoadExpKeys;
    evs = LoadEvents([]);
    cfg_spk = [];
    cfg_spk.fc = ExpKeys.goodCell;
    if isempty(cfg_spk.fc)
        return
    end
    % if not the required cell, skip
    if ~strcmp(cfg_spk.fc, label)
        return
    end
    % cfg_spk.min_cluster_quality = 3;
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg_spk.getRatings = 1;
        cfg_spk.uint = '64';
    end
    S = LoadSpikes(cfg_spk);
    
    % Set variables
    axisLabelFS = 30;
    tickLabelFS = 20;
    nbins = 5;
%     fbands = {[2 5], [6 10], [12 30], [30 55]};
    fbands = {[2 5], [6 10], [30 55]};
    c_list = {'red', 'cyan','magenta', 'green', 'blue'}; % One color for each phase bin
    temp = linspecer(nbins);
    c_rgb = {};
    for iB = 1:nbins
        c_rgb{iB} = temp(iB,:);
    end
%     c_rgb = {[1 0 0], [0 1 1], [1 0 1], [0 1 0], [0 0 1]};
    clear temp
    
    % Remove spikes during short stim times
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

    this_fig = figure('WindowState','maximized');
    % Load phase at stim
    load('stim_phases.mat');
    
    for iC = 1:length(restricted_S.label)
        fn_prefix = extractBefore(restricted_S.label{iC}, '.t');
        % Load 
        load(strcat(fn_prefix, '_phase_response_', string(nbins), '_bins.mat'));
        delta_fr = out.fr.bin;
          
        this_cell = SelectTS([], S, iC);
        restricted_cell = SelectTS([], restricted_S, iC);
        goodTrials = ExpKeys.goodTrials(iC,:);

        % Trial-stim raster
        if ~isempty(ExpKeys.stim_times)
            ax = subplot(1,1,1);
            this_on_events = stim_on;
            if ~isempty(ExpKeys.goodTrials) % Only keep the good trials
                this_on_events = this_on_events(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2));
            end
            hold on
            start_idx = 810;
            max_trials = 25;
            % Use delta phase for illustration   
            this_phase = causal_phase(1,goodTrials(1):goodTrials(2));
            phase_bins = -pi:2*pi/nbins:pi;
            [~, ~, this_bin] = histcounts(this_phase, phase_bins);
            for iT = 1:max_trials
                   fill([-10,-10,10,10],[iT-0.5,iT+0.5,iT+0.5,iT-0.5], ...
                       c_rgb{this_bin(start_idx+iT)}, 'FaceAlpha', 0.25, 'LineStyle', 'none')
            end
            % Use MultiRaster to do this, because why not
            fake_S = this_cell;
            for iT = 1:max_trials
                temp_S = restrict(S, iv(this_on_events(start_idx+iT) - 0.01, this_on_events(start_idx+iT) + 0.01));
                fake_S.t{iT} = 1000*(temp_S.t{1} - this_on_events(start_idx+iT));
                fake_S.label{iT} = temp_S.label{1};
            end
            cfg = [];
            cfg.SpikeHeight = 0.48;
            cfg.openNewFig = 0;
            MultiRaster(cfg,fake_S);
          
            xline(0, '--cyan', 'LineWidth', 1);
            xline(ExpKeys.short_stim_pulse_width+stop_delay*1000, '--cyan', 'LineWidth', 1)
            ylabel('Trials colored by LFP Phase at stim');
            xticks([-10 0 10])
            ylim([0.5 max_trials+0.5])
            xlim([-10 10]);
            xlabel("Time (ms)")
            ax.Box = 'off';
            ax.TickDir = 'out';
            ax.XAxis.FontSize = tickLabelFS;
            ax.YAxis.FontSize = tickLabelFS;
            ax.XLabel.FontSize = axisLabelFS;
            ax.YLabel.FontSize = axisLabelFS;
        end
        fontname(this_fig, 'Helvetica')
        this_fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts
        exportgraphics(this_fig, strcat('C:\Users\mvdmlab\Desktop\', fn_prefix,sprintf('-Inset-%d-%d-.eps',start_idx, max_trials)))
        close;
    end
end
