%% Script to generate the relationships between spiking and LFP Phase
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
    bin_counts = 5; % [5, 7, 10, 15, 20];
    % This parameter (pow_thresh) will decide how many trials to keep
    % based on the power in the frequency band in the last 0.25 seconds
    % prior to the stim. 0 is all trials, 50 is trials in the top 50 
    % percentile power, 75 is trials in the top 75 percentile power and
    % so on
    trial_thresh = 0;

    % Load the phases at stim_on in various frequency bands
    load('stim_phases.mat');
    mean_power = cellfun(@(x) mean(x), causal_power);
    power_thresh = prctile(mean_power,trial_thresh,2);

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
            goodTrials = [min(min(ExpKeys.goodTrials)), max(max(ExpKeys.goodTrials))];
        else
            goodTrials = ExpKeys.goodTrials(iGC,:);
        end

        % Generate indices for shuffling
            
        % Treating trial_thresh = 0 as a separate case to keep the results
        % from before the same as before
        [freqwise_trials, freqwise_shuf_idx] = deal(cell(size(fbands)));
        if trial_thresh == 0 % include all trials
            shuf_idx = zeros(nshufs, diff(goodTrials)+1);
            for iShuf = 1:nshufs
                shuf_idx(iShuf,:) = randperm(diff(goodTrials)+1);
            end
            for iF = 1:length(fbands)
                freqwise_trials{iF} = true(1,diff(goodTrials)+1);
                freqwise_shuf_idx{iF} = shuf_idx;
            end
            clear shuf_idx
        else
            for iF = 1:length(fbands)
                keep_trials = mean_power(iF,:) >= power_thresh(iF);
                freqwise_trials{iF} = keep_trials(goodTrials(1):goodTrials(2));
                shuf_idx = zeros(nshufs, sum(freqwise_trials{iF}));
                for iShuf = 1:nshufs
                    shuf_idx(iShuf,:) = randperm(sum(freqwise_trials{iF}));
                end
                freqwise_shuf_idx{iF} = shuf_idx;
                clear keep_trials shuf_idx
            end
        end
        
        % Repeat for each bin count
        for iN = 1:length(bin_counts)
            this_fig = figure('WindowState', 'maximized');
            nbins = bin_counts(iN);
            phase_bins = -pi:2*pi/nbins:pi;
            x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));
            % Results field to be saved at the end
            [out, out.lat, out.lat_ws, out.fr, out.fr_ws] = deal([]);
            out.nshufs = nshufs;
            out.bin_count = nbins;
            [out.lat.bin, out.lat_ws.bin, out.fr.bin, out.fr_ws.bin, out.fr.sd, out.fr_ws.sd] = deal(nan(length(fbands),nbins));
            [out.lat.ratio, out.lat.zscore, out.lat_ws.ratio, out.lat_ws.zscore, ...
                out.fr.ratio, out.fr.zscore, out.fr_ws.ratio, out.fr_ws.zscore] = deal(nan(1,length(fbands)));

            % Since all calcluations need to happen differently for each
            % frequency band separately, start a loop here
            for iF = 1:length(fbands)
                this_trial_idx = (goodTrials(1):goodTrials(2)) & freqwise_trials{iF};
                all_lat = od.trial_stim.latency(this_trial_idx);
                all_lat_ws = od.trial_stim.latency_wo_stim(this_trial_idx);
                % Now subtracting baseline firing rate
                all_bfr = od.trial_stim.bfr(this_trial_idx);
                all_fr = od.trial_stim.fr(this_trial_idx) - all_bfr;
                all_fr_ws = od.trial_stim.fr_wo_stim(this_trial_idx) - all_bfr;

                this_phase = causal_phase(iF,this_trial_idx);
                [this_count, ~, this_bin] = histcounts(this_phase, phase_bins);

                [this_lat, this_lat_ws, this_fr, this_fr_ws, this_fr_sd, this_fr_ws_sd] = deal(zeros(size(this_count)));
                for iB = 1:nbins
                   this_lat(iB) = sum(~isnan(all_lat(this_bin==iB)))/this_count(iB);                   
                   this_lat_ws(iB) = sum(~isnan(all_lat_ws(this_bin==iB)))/this_count(iB);
                   this_fr(iB) = mean(all_fr(this_bin==iB));
                   this_fr_sd(iB) = std(all_fr(this_bin==iB));
                   this_fr_ws(iB) = mean(all_fr_ws(this_bin==iB));
                   this_fr_ws_sd(iB) = std(all_fr_ws(this_bin==iB));
                end

                % Using (max-min)/(max+min) as the ratio
                out.lat.bin(iF,:)= this_lat;
                out.lat.ratio(iF) = (max(this_lat) - min(this_lat))/(max(this_lat) + min(this_lat));
                out.lat_ws.bin(iF,:) = this_lat_ws;
                out.lat_ws.ratio(iF) = (max(this_lat_ws) - min(this_lat_ws))/(max(this_lat_ws) + min(this_lat_ws));
                out.fr.bin(iF,:) = this_fr;
                out.fr.ratio(iF) = abs(max(this_fr) - min(this_fr))/(abs(max(this_fr)) + abs(min(this_fr)));
                out.fr_ws.bin(iF,:) = this_fr_ws;
                out.fr_ws.ratio(iF) = abs(max(this_fr_ws) - min(this_fr_ws))/(abs(max(this_fr_ws)) + abs(min(this_fr_ws)));
                out.fr.sd(iF,:) = this_fr_sd;
                out.fr_ws.sd(iF,:) = this_fr_ws_sd;

                % Generate shuffles
                [shuf_lat_ratio, shuf_lat_ws_ratio, shuf_fr_ratio, shuf_fr_ws_ratio] = deal(zeros(nshufs, 1));
                for iShuf = 1:nshufs
                    % shuffle the bin indices
                    shuf_bin = this_bin( freqwise_shuf_idx{iF}(iShuf,:));
                    [shuf_lat, shuf_lat_ws, shuf_fr, shuf_fr_ws] = deal(zeros(1,nbins));
                    for iB = 1:nbins
                       shuf_lat(iB) = sum(~isnan(all_lat(shuf_bin==iB)))/this_count(iB);
                       shuf_lat_ws(iB) = sum(~isnan(all_lat_ws(shuf_bin==iB)))/this_count(iB);
                       shuf_fr(iB) = sum(all_fr(shuf_bin==iB))/this_count(iB);
                       shuf_fr_ws(iB) = sum(all_fr_ws(shuf_bin==iB))/this_count(iB);
                    end
                    shuf_lat_ratio(iShuf) = (max(shuf_lat) - min(shuf_lat))/(max(shuf_lat) + min(shuf_lat));
                    shuf_lat_ws_ratio(iShuf) = (max(shuf_lat_ws) - min(shuf_lat_ws))/(max(shuf_lat_ws) + min(shuf_lat_ws));
                    shuf_fr_ratio(iShuf) = abs(max(shuf_fr) - min(shuf_fr))/(abs(max(shuf_fr)) + abs(min(shuf_fr)));
                    shuf_fr_ws_ratio(iShuf) = abs(max(shuf_fr_ws) - min(shuf_fr_ws))/(abs(max(shuf_fr_ws)) + abs(min(shuf_fr_ws)));
                    % Also save shufs for significance calculation and description in figures later
                    out.lat.shufs(iF,iShuf,:) = shuf_lat'; 
                    out.lat_ws.shufs(iF,iShuf,:) = shuf_lat_ws';
                    out.fr.shufs(iF,iShuf,:) = shuf_fr';
                    out.fr_ws.shufs(iF,iShuf,:) = shuf_fr_ws';
                end
                % Saving Z-scored values
                out.lat.zscore(iF) = (out.lat.ratio(iF) - mean(shuf_lat_ratio))/std(shuf_lat_ratio);
                out.lat_ws.zscore(iF) = (out.lat_ws.ratio(iF) - mean(shuf_lat_ws_ratio))/std(shuf_lat_ws_ratio);
                out.fr.zscore(iF) = (out.fr.ratio(iF) - mean(shuf_fr_ratio))/std(shuf_fr_ratio);
                out.fr_ws.zscore(iF) = (out.fr_ws.ratio(iF) - mean(shuf_fr_ws_ratio))/std(shuf_fr_ws_ratio);

                % Plot
                ax = subplot(5,5,(iF-1)*5+1);
                bar(ax,x_ticks,this_count/sum(this_count),1,c_list{iF});
                ax.Title.String = sprintf('Stim Phase distribution');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
                ax.YLabel.String = sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2));
    
                ax = subplot(5,5,(iF-1)*5+2);
                bar(ax,x_ticks,this_lat,1,c_list{iF});
                ax.Title.String = sprintf('Response prob.');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
    
                ax = subplot(5,5,(iF-1)*5+3);
                bar(ax,x_ticks,this_lat_ws,1,c_list{iF});
                ax.Title.String = sprintf('Response prob. w/o stim');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
    
                ax = subplot(5,5,(iF-1)*5+4);
                bar(ax,x_ticks,this_fr,1,c_list{iF});
                ax.Title.String = sprintf('Firing rate');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
    
                ax = subplot(5,5,(iF-1)*5+5);
                bar(ax,x_ticks,this_fr_ws,1,c_list{iF});
                ax.Title.String = sprintf('Firing rate w/o stim');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];       
            end
                
%             ax = subplot(5,5,21);
%             ax.FontSize = 25;
%             hold off
%             text(0, 0.75, sprintf('Overall response ratio: % .2f', out.overall_response))
%             text(0, 0.25, sprintf('Overall response w/o stim ratio: % .2f', out.overall_response_ws));
%             ax.Box = 'off';
%             ax.Visible = 'off';
    
            ax = subplot(5,5,22);
            ax.FontSize = 12;
            hold on
            for iF = 1:length(fbands)
                scatter(mean(fbands{iF}), out.lat.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                    'MarkerEdgeColor', c_list{iF})
            end
            ax.XLim = [0 60];
    %         ax.YLim = [-3 10];
            ax.XAxis.Label.String = 'Freq (Hz)';
            ax.YAxis.Label.String = 'Z-score';
    
            ax = subplot(5,5,23);
            ax.FontSize = 12;
            hold on
            for iF = 1:length(fbands)
                scatter(mean(fbands{iF}), out.lat_ws.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                    'MarkerEdgeColor', c_list{iF})
            end
            ax.XLim = [0 60];
    %         ax.YLim = [-3 10];
            ax.XAxis.Label.String = 'Freq (Hz)';
            ax.YAxis.Label.String = 'Z-score';
                
            ax = subplot(5,5,24);
            ax.FontSize = 12;
            hold on
            for iF = 1:length(fbands)
                scatter(mean(fbands{iF}), out.fr.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                    'MarkerEdgeColor', c_list{iF})
            end
            ax.XLim = [0 60];
    %         ax.YLim = [-3 10];
            ax.XAxis.Label.String = 'Freq (Hz)';
                ax.YAxis.Label.String = 'Z-score';
        
            ax = subplot(5,5,25);
            ax.FontSize = 12;
            hold on
            for iF = 1:length(fbands)
                scatter(mean(fbands{iF}), out.fr_ws.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                    'MarkerEdgeColor', c_list{iF})
            end
            ax.XLim = [0 60];
    %         ax.YLim = [-3 10];
            ax.XAxis.Label.String = 'Freq (Hz)';
            ax.YAxis.Label.String = 'Z-score';
            
            fn_prefix = strcat(fn_prefix, '_top_', num2str(100-trial_thresh), 'pct_trials');
            sgtitle(fn_prefix, 'Interpreter', 'none');
            % Save variables and figure
            save(strcat(fn_prefix, '_phase_response_',string(bin_counts(iN)),'_bins'), 'out');
            print(this_fig, '-dpng', strcat(fn_prefix, '_phase_response_hist_', string(bin_counts(iN)),'_bins'));
    %             savefig(this_fig, strcat(fn_prefix, '_phase_response_hist_', string(nbins),'_bins'));
            close;
        end
    end
end
