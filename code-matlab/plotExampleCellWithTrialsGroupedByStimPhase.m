%% Script to generate Figures 4 (partially), 5 and 8 (partially)
top_dir = 'data\'
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
%         this_label = 'M019-2019-04-14-TT05_1.t'; % Fig 4A vStr Didactic example; choose trials 810-835 for inset
%         this_label = 'M016-2019-02-18-TT04_1.t'; % Fig 5A vStr
%         this_label = 'M018-2019-04-14-TT03_1.t'; % Fig 5B vStr
%         this_label = 'M017-2019-02-19-TT05_1.t'; % Fig 5C dStr
%          this_label = 'M295-2022-01-06-TT06_4.t'; % dStr clearly inhibited non-opto cell
        this_label = 'M016-2019-02-15-TT06_2.t'; % 2nd vStr clearly inhibited non-opto cell
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
%     % if not the required cell, skip
%     if ~strcmp(cfg_spk.fc, label)
%         return 
%     end
    % cfg_spk.min_cluster_quality = 3;
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg_spk.getRatings = 1;
        cfg_spk.uint = '64';
    end
    cfg_spk.fc = {label}; % Hacky insert
    S = LoadSpikes(cfg_spk);
    
    % Set variables
    axisLabelFS = 30;
    tickLabelFS = 30;
    nbins = 5;
    fbands = {[2 5], [6 10], [12 28], [30 55]};
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
    axs = [];
    % Load phase at stim
    load('stim_phases.mat');
    
    for iC = 1:length(restricted_S.label)
        fn_prefix = extractBefore(restricted_S.label{iC}, '.t');
        % Load 
        load(strcat(fn_prefix, '_phase_response20ms_', string(nbins), '_bins.mat'));
        delta_fr = out.fr.bin;
          
        this_cell = SelectTS([], S, iC);
        restricted_cell = SelectTS([], restricted_S, iC);
        goodTrials = ExpKeys.goodTrials(iC,:);

%         % Trial-stim raster
%         if ~isempty(ExpKeys.stim_times)
%             ax = subplot(4,4,[1,5,9]);
%             ax = subplot(3,4,[1,5]);
%             this_on_events = stim_on;
%             if ~isempty(ExpKeys.goodTrials) % Only keep the good trials
%                 this_on_events = this_on_events(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2));
%             end
%             % Use MultiRaster to do this, because why not
%             fake_S = this_cell;
%             for iT = 1:length(this_on_events)
%                 temp_S = restrict(S, iv(this_on_events(iT) - 0.01, this_on_events(iT) + 0.01));
%                 fake_S.t{iT} = 1000*(temp_S.t{1} - this_on_events(iT));
%                 fake_S.label{iT} = temp_S.label{1};
%             end
%             cfg = [];
%             cfg.SpikeHeight = 0.48;
%             cfg.openNewFig = 0;
%             MultiRaster(cfg,fake_S);
%             hold on
%             xline(0, '--cyan', 'LineWidth', 1);
%             xline(ExpKeys.short_stim_pulse_width+stop_delay*1000, '--cyan', 'LineWidth', 1)
%             yticks([0 length(this_on_events)])
%             xticks([-10 0 10])
%             ylim([0 length(this_on_events)])
%             xlim([-10 10]);
%             % Plot the box around the inset trials (only for the didactic example)
%             plot([-10 10], [810 810], 'Color', 'black', 'LineWidth', 2);
%             plot([-10 10], [835 835], 'Color', 'black', 'LineWidth', 2);
%             plot([-10 -10], [810 835], 'Color', 'black', 'LineWidth', 2);
%             plot([10 10], [810 835], 'Color', 'black', 'LineWidth', 2);
%             xlabel("Time (ms)")
%             ylabel('Trial #');
%             ax.Box = 'off';
%             ax.TickDir = 'out';
%             ax.TickLength = [0.03 0.02];
%             ax.XAxis.FontSize = tickLabelFS;
%             ax.YAxis.FontSize = tickLabelFS;
%             ax.XLabel.FontSize = axisLabelFS;
%             ax.YLabel.FontSize = axisLabelFS;
%             axs = [axs, ax];
%         end
% 
%         % Plot PSD on the bottom left
%         load('trial_psd.mat');
% %         ax = subplot(4,4,13);
%         ax = subplot(3,4,9);
% %         ax.PositionConstraint = "innerposition";
%         hold on;
%         plot(psd.freq, 10*log10(psd.original), 'black', 'LineWidth', 1.5);
%         for iF = 1:length(fbands)
%             f_idx = find(round(psd.freq) >= fbands{iF}(1) & round(psd.freq) <= fbands{iF}(2));
% %             area(psd.freq(f_idx), 10*log10(psd.original(f_idx)), 'FaceColor', c_list{iF}, ...
% %                 'FaceAlpha', 0.5, 'BaseValue', min(10*log10(psd.original)))
%             area(psd.freq(f_idx), 10*log10(psd.original(f_idx)), 'FaceColor', [0.7 0.7 0.7], ...
%                 'FaceAlpha', 0.5, 'BaseValue', min(10*log10(psd.original)))
%         end
%         plot(psd.freq, 10*log10(psd.irasa), '--black', 'LineWidth', 1.5);
%         xlim([0 100])
% %         yticklabels([])
% %         yticks([0 length(this_on_events)])
%         ylabel('PSD')
%         xlabel('Frequency (Hz)')
%         ax.YLim(1) =  min(10*log10(psd.original));
%         ax.Box = 'off';
%         ax.TickDir = 'out';
%         ax.TickLength = [0.06 0.02];
%         ax.XAxis.FontSize = tickLabelFS;
%         ax.YAxis.FontSize = tickLabelFS;
%         ax.XLabel.FontSize = axisLabelFS;
%         ax.YLabel.FontSize = axisLabelFS;
%         axs = [axs, ax];

        this_on_events = stim_on;
        if ~isempty(ExpKeys.goodTrials) % Only keep the good trials
            this_on_events = this_on_events(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2));
        end
        fake_S = this_cell;
        for iT = 1:length(this_on_events)
            temp_S = restrict(S, iv(this_on_events(iT) - 0.01, this_on_events(iT) + 0.01));
            fake_S.t{iT} = 1000*(temp_S.t{1} - this_on_events(iT));
            fake_S.label{iT} = temp_S.label{1};
        end

        % outputT has the trial numbers, and outputS has the spiketiming within
        % that given trial. If a trial has no spikes, it doesn't exist in outputT
        % For each frequency band, divide the trials into phase bins and group them
        % by that
        phase_bins = -pi:2*pi/nbins:pi;
        for iF = 1:length(fbands)
%             ax = subplot(4,4,[iF+1,iF+5,iF+9]);
            ax = subplot(3,4,[iF,iF+4]);
            hold on
            this_phase = causal_phase(iF,goodTrials(1):goodTrials(2));
            [this_count, ~, this_bin] = histcounts(this_phase, phase_bins);
            % Group the trials according to the bin assigned in this_bin
            fake_S2 = fake_S;
            tc = 0;
            for iBin = 1:nbins
                bin_trials = find(this_bin == iBin);
                for iT = 1:this_count(iBin)
                    tc = tc+1;
                    fake_S2.t{tc} = fake_S.t{bin_trials(iT)};
                end
                if iBin == 1
                   fill([-10,-10,10,10],[0,sum(this_count(1:iBin)),sum(this_count(1:iBin)),0], ...
                       c_rgb{iBin}, 'FaceAlpha', 0.5, 'LineStyle', 'none')
                else
                    fill([-10,-10,10,10],[sum(this_count(1:iBin-1)),sum(this_count(1:iBin)), ...
                        sum(this_count(1:iBin)),sum(this_count(1:iBin-1))], ...
                       c_rgb{iBin}, 'FaceAlpha', 0.5, 'LineStyle', 'none')
                end
            end
            assert(tc == goodTrials(2) - goodTrials(1) + 1, "All trials are not regrouped");
            cfg = [];
            cfg.SpikeHeight = 0.48;
            cfg.openNewFig = 0;
            MultiRaster(cfg,fake_S2);
            xline(0, '--cyan', 'LineWidth', 1);
            xline(ExpKeys.short_stim_pulse_width+stop_delay*1000, '--cyan', 'LineWidth',1)
            hold on
            xlabel("Time (ms)")
            yticks([0 length(this_on_events)])
            xticks([-10 0 10])
            ax.XLim = [-10 10];
            ax.YLim = goodTrials;
            if iF == 1
                ax.YLabel.String = 'Trials grouped by LFP phase-bin';
            else
                ax.YLabel.String = {};
            end
            ax.TickDir = 'out';
            ax.Title.String = sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2));
            ax.Title.FontSize = 45;
            ax.Box = 'off';
            ax.TickLength = [0.03 0.02];
            ax.XAxis.FontSize = tickLabelFS;
            ax.YAxis.FontSize = tickLabelFS;
            ax.XLabel.FontSize = axisLabelFS;
            ax.YLabel.FontSize = axisLabelFS;
            axs = [axs, ax];

%             ax = subplot(4,4,iF+13);
            ax = subplot(3,4,iF+8);
%             ax.PositionConstraint = "outerposition";
            hold on
            b = bar(delta_fr(iF,:),1, 'EdgeColor', 'none');
            b.FaceColor = 'flat';
            b.FaceAlpha = 0.5;
            for iBin = 1:nbins
                b.CData(iBin,:) = c_rgb{iBin};
            end
%             errorbar(1:5, delta_fr(iF,:), out.fr.sd(iF,:), 'black', 'LineStyle', 'none');
%             xticklabels([]);
            q0 = squeeze(out.fr.shufs(iF,:,:));
            q1 = max(q0,[],2);
            q2 =  min(q0,[],2);
            q3 = (q1 - q2)./(q2 + q1);
            qz2 = 2*std(q3)+ mean(q3);
            [qmin, min_idx] = min(out.fr.bin(iF,:));
            [qmax, max_idx] = max(out.fr.bin(iF,:));
            qoff = qmin*(qz2 + 1)/(1 - qz2);
            plot([min_idx, 5.65], [qmin, qmin], '--black')   
            plot([max_idx, 5.65], [qmax, qmax], '--black')
            errorbar(5.75, 0.5*(qmax+qmin), 0.5*(qmax-qmin), 'black', 'LineStyle', 'none');
            errorbar(6, 0.5*(qoff+qmin), 0.5*(qoff-qmin), 'red', 'LineStyle', 'none');
%             plot([6, 6], [qmin, qmax], '--black', 'LineWidth', 1)
%             plot([6.25,6.25], [qmin, qmin+qoff], '--red', 'LineWidth', 1)
            clear q0 q1 q2 q3 qz2 qoff qmin qmax min_idx max_idx
            xticks(1:5)
            ax.XLim = [0.5 6.5];
            ax.TickDir = 'out';
            if iF == 1
                ax.YLabel.String = '{\Delta} FR (Hz)';
            else
                ax.YLabel.String = {};
            end
            ylim([-10 0]); % Need to change this in a case by case basis
            yticks([-10 0]);
            ax.XLabel.String = 'Phase bin';
            ax.Box = 'off';
            ax.TickLength = [0.06 0.02];
            ax.XAxis.FontSize = tickLabelFS;
            ax.YAxis.FontSize = tickLabelFS;
            ax.XLabel.FontSize = axisLabelFS;
            ax.YLabel.FontSize = axisLabelFS;
            axs = [axs, ax];
        end

        fontname(this_fig, 'Helvetica');
        this_fig.Renderer = 'painters'; % makes sure tht the figure is exported with customizable parts
%   
%         % Uncomment this only for debugging
%         fig2 = figure('WindowState', 'maximized');
%         clist2 = {'red', 'blue', 'green'};
%         for iF2 = 1:3
%             subplot(1,3,iF2)
%             hold on
%             q0 = squeeze(out.fr.shufs(iF2,:,:));
%             q1 = max(q0,[],2);
%             q2 =  min(q0,[],2);
%             q3 = (q1 - q2)./(q2 + q1);
%             qz2 = 2*std(q3)+ mean(q3);
%             histogram(q3, 'FaceColor',clist2{iF2});
%             xline(out.fr.ratio((iF2)), '--black', 'LineWidth', 2);
%             xline(qz2, '--red', 'LineWidth', 2);
%             xlabel('Modulation Stength')
%             title(sprintf('%d - %d Hz', fbands{iF2}(1), fbands{iF2}(2)));
%         end
%         sgtitle('Shuffle Histograms')
%         clear q0 q1 q2 q3 clist2 iF2
%         close(fig2)

        % Need to pause here for this to work
%         pause(2)
%         % Manual plot manipulation to make the final figure looks pretty
%         % Make a copy of original positions
%         daxs_in = [];
%         daxs_out = [];
%         for i = 1:length(axs)
%             daxs_in = [daxs_in; axs(i).Position];
%             daxs_out = [daxs_out; axs(i).OuterPosition];
%         end
% 
%         axs(2).OuterPosition(1) = axs(1).OuterPosition(1);
%         axs(1).OuterPosition(2) = axs(1).OuterPosition(2)+0.05;
%         axs(2).Position(1) = axs(1).Position(1);
%         axs(2).Position(3) = axs(1).Position(3);
% 
%         axs(3).OuterPosition(1) = axs(1).OuterPosition(1) + axs(1).OuterPosition(3);
%         axs(3).OuterPosition(2) = axs(3).OuterPosition(2)+0.05;
%         axs(4).Position(1) = axs(3).Position(1);
%         axs(4).Position(3) = axs(3).Position(3)+0.03;
% 
%         axs(5).OuterPosition(1) = axs(3).OuterPosition(1) + axs(3).OuterPosition(3);
%         axs(5).OuterPosition(2) = axs(5).OuterPosition(2)+0.05;
%         axs(6).Position(1) = axs(5).Position(1);
%         axs(6).Position(3) = axs(5).Position(3)+0.03;
% 
%         axs(7).OuterPosition(1) = axs(5).OuterPosition(1) + axs(5).OuterPosition(3);
%         axs(7).OuterPosition(2) = axs(7).OuterPosition(2)+0.05;
%         axs(8).Position(1) = axs(7).Position(1);
%         axs(8).Position(3) = axs(7).Position(3)+0.03;

        exportgraphics(this_fig, strcat('output\', fn_prefix,'-TrialsGroupedByStimPhase.eps'))
%         savefig(this_fig, strcat('data\', fn_prefix,'-TrialsGroupedByStimPhase'));
%         print(this_fig, '-dpdf', '-fillpage', strcat(fn_prefix,'-CellReport'));
%         print(this_fig, '-dpng',  strcat('data\', fn_prefix, '-CellReport'));
        close;
    end
end
