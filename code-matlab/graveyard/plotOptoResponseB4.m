cd('data\M019\M019-2019-04-14_vStr_4p2_light_cells_TT5 _min\');
%     cd('data\M018\M018-2019-04-12_dStr_3p8_light_cells_TT7_min\'); % Select Cell
%     cd('data\M319\M319-2022-06-28\'); % Select Cell
% cd('data\M019\M019-2019-04-12_vStr_4p2_light_cells_TT7_min');
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

%% Remove spikes during short stim times
if contains(ExpKeys.light_source, 'LASER')
    start_delay = 0.0011;
    stop_delay = 0.0012; %0.0022 is the spike width in the case of 1 msec laser pulse
else
    start_delay = 0;
    stop_delay = 0;
end

stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
    stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                      stim_on <= ExpKeys.stim_times(2));
else
    stim_on = [];
end

for iC = 1:length(S.label)
    % Set variables and parameters
    fbands = {[2 5], [6 10], [30 55]};
%     c_list = {'red', 'cyan','magenta', 'green', 'blue'};
    c_list = {'red', 'blue','magenta', 'green'};
    nbins = 5;
    phase_bins = -pi:2*pi/nbins:pi;
    x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));

    fn_prefix = extractBefore(S.label{iC}, '.t');
    fn_prefix(1:end-2) = strrep(fn_prefix(1:end-2), '_', '-');

    goodTrials = ExpKeys.goodTrials(iC,:);

    % Load the phases at stim_on in various frequency bands
    load('stim_phases.mat');
    causal_phase = causal_phase(:,goodTrials(1):goodTrials(2));
    stim_response = load(strcat(fn_prefix, '_stim_response.mat'));
    load(strcat(fn_prefix, '_phase_response.mat'));
    fr = stim_response.od.trial_stim.fr(goodTrials(1):goodTrials(2));
    this_cell = SelectTS([], S, iC);
    % Trial-stim raster
        this_on_events = stim_on;
        if ~isempty(ExpKeys.goodTrials) % Only keep the good trials
            this_on_events = this_on_events(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2));
        end
        for iStim = length(this_on_events)
%         for iStim = 1:10
            fig = figure('WindowState','maximized');
            ax = subplot(4,2,[1,3,5,7]);
            hold off;
            temp_on = this_on_events(1:iStim); % testing
            temp_fr = fr(1:iStim);
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm([], ...
                this_cell, temp_on, '', 0.5);
            hold on
            plot(outputS, outputT, 'k.', 'MarkerSize', 20);
            plot([0 0], [0 1600], 'color', 'red', 'linewidth', 1);
            plot([ExpKeys.short_stim_pulse_width+stop_delay ExpKeys.short_stim_pulse_width+stop_delay], ...
                [0 1600], 'color', 'red', 'linewidth', 1);
            ylabel('Stim #');
%             ylim([0 100])
                ylim([0 1500])
            xlim([-0.01 0.01]);
%             ax.YTick = [1 10 50 100];
            ax.YTick = [1, 100, 500, iStim];
            ax.TickDir = 'out';
            area([0, ExpKeys.short_stim_pulse_width+stop_delay], [1600 1600], ...
                'FaceColor', 'cyan', 'FaceAlpha', 0.1);
            xlabel('Time (sec)');
            ax.YDir = 'reverse';
            ax.YAxis.Label.FontWeight = 'bold';
            
            % Plot Frequency Hists
            for iF = 1:4
                ax = subplot(4,2,2*iF);
                hold off
                this_phase = causal_phase(iF,1:iStim);                
                [this_count, ~, this_bin] = histcounts(this_phase, phase_bins);
                this_fr = zeros(size(this_count));
                for iB = 1:nbins
                    if this_count(iB) ~= 0
                        this_fr(iB) = sum(temp_fr(this_bin==iB))/this_count(iB);
                    end
                
                end
                bar(ax,x_ticks,this_fr,1,c_list{iF});
                ax.Title.String = sprintf('%d Hz - %d Hz, mod-depth: %.2f, zscore: %.2f ', ...
                    fbands{iF}(1), fbands{iF}(2), out.fr.ratio(iF), out.fr.zscore(iF));
%                 ax.Title.FontSize = 14;
                ax.XLim = [-3.25 3.25];
                xlabel('Phase bins', 'FontSize', 15);
                ylabel('Firing Rate', 'FontSize',15);
            end
            sgtitle(fn_prefix);
            % Save this in some folder
            print(fig, '-dpng', '-r300', strcat('data\', ...
                fn_prefix, '_stim', num2str(iStim)));
            close;
        end

        dummy = 1;

end
