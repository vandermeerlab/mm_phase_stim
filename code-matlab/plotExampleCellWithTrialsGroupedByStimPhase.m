%% Assumes that spike-sorting has been done

top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
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
    text_fs = 11; % Text font size
    small_label_fs = 12; 
    big_label_fs = 16; 
    max_delay = 0.01; % sec (window for the first response since stimulus)
    nbins = 5;
    fbands = {[2 5], [6 10], [12 30], [30 55]};
    c_list = {'red', 'cyan','magenta', 'green', 'blue'}; % One color for each phase bin
    c_rgb = {[1 0 0], [0 1 1], [1 0 1], [0 1 0], [0 0 1]};
    
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
            ax = subplot(4,5,[1,6,11]);
            this_on_events = stim_on;
            if ~isempty(ExpKeys.goodTrials) % Only keep the good trials
                this_on_events = this_on_events(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2));
            end
            [outputS, outputT, ~, ~, ~] = SpikePETHvdm([], ...
                this_cell, this_on_events, '', 0.5);
            hold on
            plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
            plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            plot([ExpKeys.short_stim_pulse_width+stop_delay ExpKeys.short_stim_pulse_width+stop_delay], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            ylabel(' Trial Stim #');
            ylim([1 length(this_on_events)])
            xlim([-0.01 0.01]);
            xlabel("Time (sec)")
        end
%%
        % outputT has the trial numbers, and outputS has the spiketiming within
        % that given trial. If a trial has no spikes, it doesn't exist in outputT
        % For each frequency band, divide the trials into phase bins and group them
        % by that
        phase_bins = -pi:2*pi/nbins:pi;
        x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));
        for iF = 1:length(fbands)
            ax = subplot(4,5,[iF+1,iF+6,iF+11]);
            if iF == 1
                ax.YLabel.String = 'Trials grouped by LFP phase-bin';
            end
            hold on
            this_phase = causal_phase(iF,goodTrials(1):goodTrials(2));
            [this_count, ~, this_bin] = histcounts(this_phase, phase_bins);
            outputB = this_bin(outputT); % This works because outputT has trial numbers as the members
            group_end = zeros(1,nbins);
            for iBin = 1:nbins
                sel = (outputB == iBin);
                bin_trials = outputT(sel);
                groupS = outputS(sel);
                groupT = zeros(size(groupS));
                groupT(1) = 1;
                for iT = 2:length(bin_trials)
                    if bin_trials(iT) == bin_trials(iT-1)
                        groupT(iT) = groupT(iT-1);
                    else
                       groupT(iT) = groupT(iT-1)+1;
                    end
                end
                if iBin == 1
                   group_end(iBin) =  groupT(end);
                   fill([-0.1,-0.1,0.1,0.1],[0,group_end(iBin),group_end(iBin),0], ...
                       c_list{iBin}, 'FaceAlpha', 0.25, 'LineStyle', 'none')
                else
                    group_end(iBin) = group_end(iBin-1) + groupT(end);
                    groupT = groupT + group_end(iBin-1); 
                    fill([-0.1,-0.1,0.1,0.1],[group_end(iBin-1),group_end(iBin),group_end(iBin),group_end(iBin-1)], ...
                       c_list{iBin}, 'FaceAlpha', 0.25, 'LineStyle', 'none')
                end
                plot(groupS, groupT+0.5, 'k.', 'MarkerSize', 5);
            end
            plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            plot([ExpKeys.short_stim_pulse_width+stop_delay ExpKeys.short_stim_pulse_width+stop_delay], ...
                [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
            ax.XLim = [-0.01 0.01];
            ax.YLim = goodTrials;
            ax.YTick = 0.5*([0,group_end(1:end-1)]+group_end(1:end));
%             ax.YTickLabel = {'1', '2', '3', '4', '5'};
            ax.YTickLabel = {};
            ax.TickDir = 'out';
            ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));

            ax = subplot(4,5,iF+16);
            b = bar(delta_fr(iF,:),1);
            b.FaceColor = 'flat';
            b.FaceAlpha = 0.25;
            for iBin = 1:nbins
                b.CData(iBin,:) = c_rgb{iBin};
            end
            ax.TickDir = 'out';
            ax.YLabel.String = '{\Delta} FR';
            ax.XLabel.String = 'Phase Bins';
            ax.YLim = [0 60];
            dummy = 1;
        end



        savefig(this_fig, strcat('C:\Users\mvdmlab\Desktop\', fn_prefix,'-TrialsGroupedByStimPhase'));
%         print(this_fig, '-dpdf', '-fillpage', strcat(fn_prefix,'-CellReport'));
%         print(this_fig, '-dpng',  strcat('E:\Dropbox (Dartmouth College)\EC_State_inProcess\', fn_prefix, '-CellReport'));
        close;
    end
end
