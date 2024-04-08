%% Script to generate the relationships between spiking and LFP Phase separately for trials grouped into two halves based on the oscillatory power
% Assumes that stim_phases.mat and *stim_response.mat already exist in each folder
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        rng(491994); % Setting the seed for reproducibility
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff;
    end
end
%%
function doStuff
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    fbands = {[2 5], [6 10], [12 28], [30 55]};
    c_list = {'red', 'blue','magenta', 'cyan'};
    nshufs = 1000;
    nsplits = 1000;
    nbins = 5;
    all_pct = [0.5,0.55,0.6,0.65,0.7,0.75,0.9];
    sel_pct = [0.6, 0.65]; % Selected based on a good number of cells having trials
    min_trials = 600; % Skip calculating stuff if the number of 'good stim' is less than this threshold

    % Load the phases at stim_on in various frequency bands
    load('stim_phases.mat');
%     mean_power = cellfun(@(x) mean(x), causal_power);

    % Load list of trials thresholded by oscillatory power
    load('good_stim.mat');

    cfg = [];
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.min_cluster_quality = 3;
        cfg.getRatings = 1;
        cfg.uint = '64';
    end   
    S = LoadSpikes(cfg);

    for iC = 1:length(S.label)
        fn_prefix = extractBefore(S.label{iC}, '.t');
        % Load the stim_responses
        load(strcat(fn_prefix, '_stim_response.mat'));
        iGC = find(strcmp(S.label{iC}, ExpKeys.goodCell));
        if isempty(iGC)
            continue; % Don't calculate stuff for cells that are not 'good'
        else
            goodTrials = ExpKeys.goodTrials(iGC,:);
        end

        % Load the actual results
        load(strcat(fn_prefix, '_phase_response_5_bins.mat'));
%  
        phase_bins = -pi:2*pi/nbins:pi;
        x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));

        % Results field to be saved at the end
        [out_top, out_top.fr, out_top.control] = deal([]);
        out_top.nshufs = nshufs;
        out_top.nsplits = nsplits;
        out_top.bin_count = nbins;
        out_top.pct = sel_pct;
        out_top.num_trials  = zeros(length(sel_pct), length(fbands));
        out_top.min_trials = min_trials;
        [out_top.fr.bin, out_top.fr.sd] = deal(nan(length(sel_pct),length(fbands),nbins));
        [out_top.fr.ratio, out_top.fr.zscore] = deal(nan(length(sel_pct),length(fbands)));
        [out_top.fr.shufs] = nan(length(sel_pct),length(fbands),nshufs,nbins);
        
        [out_top.control.fr.ratio, out_top.control.fr.zscore] = ...
            deal(nan(length(sel_pct),nsplits,length(fbands))); 

        % Since all calcluations need to happen differently for each
        % frequency band and each threshold
        for iF = 1:length(fbands)
            for iP = 1:length(sel_pct)
                pidx = find(sel_pct(iP) == all_pct);
                this_stim = good_stim{iF, pidx};
                keep = this_stim((this_stim >= goodTrials(1)) & ...
                    (this_stim <= goodTrials(2)));
                out_top.num_trials(iP,iF) = length(keep);
                if length(keep) < min_trials 
                    % Do nothing, let everything stay nan
                    continue;
                elseif length(keep) == diff(goodTrials)+1 % All trials occur during good oscillations, so just copy over the pre-calculated results
                    out_top.fr.bin(iP,iF,:) = out.fr.bin(iF,:);
                    out_top.fr.sd(iP,iF,:) = out.fr.sd(iF,:);
                    out_top.fr.ratio(iP,iF) = out.fr.ratio(iF);
                    out_top.fr.zscore(iP,iF) = out.fr.zscore(iF);
                else % Calculate stuff
                    this_phase = causal_phase(iF,keep);
                    [this_count, ~, this_bin] = histcounts(this_phase, phase_bins);
                    all_bfr = od.trial_stim.bfr(keep);
                    all_fr = od.trial_stim.fr(keep) - all_bfr;
                    [this_fr, this_fr_sd] = deal(zeros(size(this_count)));
                    for iBin = 1:nbins
                        this_fr(iBin) = mean(all_fr(this_bin==iBin));
                        this_fr_sd(iBin) = std(all_fr(this_bin==iBin));
                    end
                    % Using abs(max-min)/(abs(max)+abs(min)) as the ratio
                    out_top.fr.bin(iP,iF,:) = this_fr;
                    out_top.fr.sd(iP,iF,:) = this_fr_sd;
                    out_top.fr.ratio(iP,iF) = abs(max(this_fr) - min(this_fr))/(abs(max(this_fr)) + abs(min(this_fr)));
                    
                    % Generate indices for shuffling
                    shuf_idx = zeros(nshufs, length(keep));
                    for iShuf = 1:nshufs
                        shuf_idx(iShuf,:) = randperm(length(keep));
                    end
                    shuf_fr_ratio = zeros(nshufs,1);
                    for iShuf = 1:nshufs
                        shuf_bin = this_bin(shuf_idx(iShuf,:));
                        shuf_fr = zeros(1,nbins);
                        for iBin = 1:nbins
                            shuf_fr(iBin) = sum(all_fr(shuf_bin==iBin))/this_count(iBin);
                        end
                        shuf_fr_ratio(iShuf) = abs(max(shuf_fr) - min(shuf_fr))/...
                            (abs(max(shuf_fr)) + abs(min(shuf_fr)));
                        out_top.fr.shufs(iP,iF,iShuf,:) = shuf_fr';
                    end
                    out_top.fr.zscore(iP,iF) = (out_top.fr.ratio(iP,iF) - ...
                        mean(shuf_fr_ratio))/std(shuf_fr_ratio); 
                
                    % clear some variables from workspace just to be safe
                    clear iShuf this_fr this_fr_sd this_count this_bin this_phase all_bfr all_fr shuf_idx

                    % generate random splits (instead of choosing the events, choose randomly)
                    for iSplit = 1:nsplits
                        this_keep = sort(randsample(goodTrials(1):goodTrials(2),length(keep)));
                        this_phase = causal_phase(iF,this_keep);
                        [this_count, ~, this_bin] = histcounts(this_phase, phase_bins);
                        all_bfr = od.trial_stim.bfr(this_keep);
                        all_fr = od.trial_stim.fr(this_keep) - all_bfr;
                        [this_fr, this_fr_sd] = deal(zeros(size(this_count)));
                        for iBin = 1:nbins
                            this_fr(iBin) = mean(all_fr(this_bin==iBin));
                            this_fr_sd(iBin) = std(all_fr(this_bin==iBin));
                        end
                        % Using abs(max-min)/(abs(max)+abs(min)) as the ratio
                        out_top.control.fr.ratio(iP,iSplit,iF) = ...
                            abs(max(this_fr) - min(this_fr))/(abs(max(this_fr)) + abs(min(this_fr)));
                    
                        % Generate indices for shuffling
                        shuf_idx = zeros(nshufs, length(this_keep));
                        for iShuf = 1:nshufs
                            shuf_idx(iShuf,:) = randperm(length(this_keep));
                        end
                        shuf_fr_ratio = zeros(nshufs,1);
                        for iShuf = 1:nshufs
                            shuf_bin = this_bin(shuf_idx(iShuf,:));
                            shuf_fr = zeros(1,nbins);
                            for iBin = 1:nbins
                                shuf_fr(iBin) = sum(all_fr(shuf_bin==iBin))/this_count(iBin);
                            end
                            shuf_fr_ratio(iShuf) = abs(max(shuf_fr) - min(shuf_fr))/...
                                (abs(max(shuf_fr)) + abs(min(shuf_fr)));
                        end
                        out_top.control.fr.zscore(iP,iSplit,iF) = (out_top.control.fr.ratio(iP,iSplit,iF) - ...
                            mean(shuf_fr_ratio))/std(shuf_fr_ratio);
                    end
                end
            end
        end
        
        % Plot stuff
        this_fig = figure('WindowState', 'maximized');
        for iF = 1:length(fbands)
            % The first 2 columns always exist 
            this_phase = causal_phase(iF,goodTrials(1):goodTrials(2));
            [this_count, ~, ~] = histcounts(this_phase, phase_bins);
            
            ax = subplot(5,6,(iF-1)*6+1);
            bar(ax,x_ticks,this_count/sum(this_count),1,c_list{iF});          
            ax.XLim = [-3.25 3.25];
            if iF == 1
                ax.Title.String = sprintf('Stim Phase - all trials');
%                 ax.Title.FontSize = 12;;
            end
            ax.YLabel.String = sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2));
            ax.XAxis.TickDirection = 'out';

            ax = subplot(5,6,(iF-1)*6+2);
            bar(ax,x_ticks,out.fr.bin(iF,:),1,c_list{iF})          
            ax.XLim = [-3.25 3.25];
            if iF == 1
                ax.Title.String = sprintf('FR mod - all trials');
%                 ax.Title.FontSize = 12;
            end
            ax.XAxis.TickDirection = 'out';

            % 3rd and 4th column exist only if the trials in this set exceed min_trials
            if out_top.num_trials(1,iF) >= out_top.min_trials
                pidx = find(sel_pct(1) == all_pct);
                this_stim = good_stim{iF, pidx};
                keep = this_stim((this_stim >= goodTrials(1)) & ...
                    (this_stim <= goodTrials(2)));
                this_phase = causal_phase(iF,keep);
                [this_count, ~, ~] = histcounts(this_phase, phase_bins);

                ax = subplot(5,6,(iF-1)*6+3);
                bar(ax,x_ticks,this_count/sum(this_count),1,c_list{iF});          
                ax.XLim = [-3.25 3.25];
                if iF == 1
                    ax.Title.String = sprintf('Stim Phase - top %d %% trials', sel_pct(1)*100);
    %                 ax.Title.FontSize = 12;;
                end
                ax.XAxis.TickDirection = 'out';
    
                ax = subplot(5,6,(iF-1)*6+4);
                bar(ax,x_ticks,squeeze(out_top.fr.bin(1,iF,:)),1,c_list{iF})          
                ax.XLim = [-3.25 3.25];
                if iF == 1
                    ax.Title.String = sprintf('FR mod - top %d %% trials', sel_pct(1)*100);
    %                 ax.Title.FontSize = 12;
                end
                ax.XAxis.TickDirection = 'out';
            end

            % 5th and 6th column exist only if the trials in this set exceed min_trials
            if out_top.num_trials(2,iF) >= out_top.min_trials
                pidx = find(sel_pct(2) == all_pct);
                this_stim = good_stim{iF, pidx};
                keep = this_stim((this_stim >= goodTrials(1)) & ...
                    (this_stim <= goodTrials(2)));
                this_phase = causal_phase(iF,keep);
                [this_count, ~, ~] = histcounts(this_phase, phase_bins);

                ax = subplot(5,6,(iF-1)*6+5);
                bar(ax,x_ticks,this_count/sum(this_count),1,c_list{iF});          
                ax.XLim = [-3.25 3.25];
                if iF == 1
                    ax.Title.String = sprintf('Stim Phase - top %d %% trials', sel_pct(2)*100);
    %                 ax.Title.FontSize = 12;;
                end
                ax.XAxis.TickDirection = 'out';
    
                ax = subplot(5,6,(iF-1)*6+6);
                bar(ax,x_ticks,squeeze(out_top.fr.bin(2,iF,:)),1,c_list{iF})          
                ax.XLim = [-3.25 3.25];
                if iF == 1
                    ax.Title.String = sprintf('FR mod - top %d %% trials', sel_pct(2)*100);
    %                 ax.Title.FontSize = 12;
                end
                ax.XAxis.TickDirection = 'out';
            end
        end
       
        % Last row is text and scatter plots
        ax = subplot(5,6,25);
        hold off
        text(0, 0.75, sprintf('Total trials: %d', diff(goodTrials)+1), 'Color', 'black', 'FontSize', 15);
        ax.Box = 'off';
        ax.Visible = 'off';

        ax = subplot(5,6,26);
        ax.FontSize = 12;
        hold on
        for iF = 1:length(fbands)
            scatter(mean(fbands{iF}), out.fr.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                'MarkerEdgeColor', c_list{iF})
        end
        ax.XLim = [0 60];
        ax.YLim = [min([out.fr.zscore,out_top.fr.zscore(1,:),out_top.fr.zscore(2,:)])-0.25, ...
            max([out.fr.zscore,out_top.fr.zscore(1,:),out_top.fr.zscore(2,:)])+0.25];
        ax.XAxis.Label.String = 'Freq (Hz)';
        ax.YAxis.Label.String = 'Z-score';

        ax = subplot(5,6,27);
        hold off
        for iF=1:length(fbands)
            text(0, 1-(iF*0.25), sprintf('Trials in %d - %d Hz: %d', ...
                fbands{iF}(1), fbands{iF}(2), out_top.num_trials(1,iF)), 'Color', c_list{iF}, 'FontSize', 15);
        end
        ax.Box = 'off';
        ax.Visible = 'off';

        ax = subplot(5,6,28);
        ax.FontSize = 12;
        hold on
        for iF = 1:length(fbands)
            scatter(mean(fbands{iF}), out_top.fr.zscore(1,iF), 'MarkerFaceColor', c_list{iF}, ...
                'MarkerEdgeColor', c_list{iF})
        end
        ax.XLim = [0 60];
        ax.YLim = [min([out.fr.zscore,out_top.fr.zscore(1,:),out_top.fr.zscore(2,:)])-0.25, ...
            max([out.fr.zscore,out_top.fr.zscore(1,:),out_top.fr.zscore(2,:)])+0.25];
        ax.XAxis.Label.String = 'Freq (Hz)';
        ax.YAxis.Label.String = 'Z-score';

        ax = subplot(5,6,29);
        hold off
        for iF=1:length(fbands)
            text(0, 1-(iF*0.25), sprintf('Trials in %d - %d Hz: %d', ...
                fbands{iF}(1), fbands{iF}(2), out_top.num_trials(2,iF)), 'Color', c_list{iF}, 'FontSize', 15);
        end
        ax.Box = 'off';
        ax.Visible = 'off';

        ax = subplot(5,6,30);
        ax.FontSize = 12;
        hold on
        for iF = 1:length(fbands)
            scatter(mean(fbands{iF}), out_top.fr.zscore(2,iF), 'MarkerFaceColor', c_list{iF}, ...
                'MarkerEdgeColor', c_list{iF})
        end
        ax.XLim = [0 60];
        ax.YLim = [min([out.fr.zscore,out_top.fr.zscore(1,:),out_top.fr.zscore(2,:)])-0.25, ...
            max([out.fr.zscore,out_top.fr.zscore(1,:),out_top.fr.zscore(2,:)])+0.25];
        ax.XAxis.Label.String = 'Freq (Hz)';
        ax.YAxis.Label.String = 'Z-score';
 
        fn_prefix = strcat(fn_prefix, '_selected_trials');
        sgtitle(fn_prefix, 'Interpreter', 'none');
        % Save variables and figure
        save(strcat(fn_prefix, '_phase_response_5_bins'), 'out_top');
        print(this_fig, '-dpng', strcat(fn_prefix, '_phase_response_hist_5_bins'));
        close all;
    end
end
