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
    bin_counts = 5; % [5, 7, 10, 15, 20];
    % This parameter (pow_thresh) will decide how many trials to keep
    % based on the power in the frequency band in the last 0.25 seconds
    % prior to the stim. 0 is all trials, 50 is trials in the top 50 
    % percentile power, 75 is trials in the top 75 percentile power and
    % so on
    trial_thresh = 50;

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
 
        % Repeat for each bin count
        for iN = 1:length(bin_counts)
            this_fig = figure('WindowState', 'maximized');
            nbins = bin_counts(iN);
            phase_bins = -pi:2*pi/nbins:pi;
            x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));
            
            % Results field to be saved at the end
            [out_top, out_top.lat, out_top.lat_ws, out_top.fr, out_top.fr_ws] = deal([]);
            out_top.nshufs = nshufs;
            out_top.bin_count = nbins;
            [out_top.lat.bin, out_top.lat_ws.bin, out_top.fr.bin, ...
                out_top.fr_ws.bin, out_top.fr.sd, out_top.fr_ws.sd] = deal(nan(length(fbands),nbins));
            [out_top.lat.ratio, out_top.lat.zscore, out_top.lat_ws.ratio, out_top.lat_ws.zscore, ...
                out_top.fr.ratio, out_top.fr.zscore, out_top.fr_ws.ratio, out_top.fr_ws.zscore] = deal(nan(1,length(fbands)));
            
            [out_bottom, out_bottom.lat, out_bottom.lat_ws, out_bottom.fr, out_bottom.fr_ws] = deal([]);
            out_bottom.nshufs = nshufs;
            out_bottom.bin_count = nbins;
            [out_bottom.lat.bin, out_bottom.lat_ws.bin, out_bottom.fr.bin, ...
                out_bottom.fr_ws.bin, out_bottom.fr.sd, out_bottom.fr_ws.sd] = deal(nan(length(fbands),nbins));
            [out_bottom.lat.ratio, out_bottom.lat.zscore, out_bottom.lat_ws.ratio, out_bottom.lat_ws.zscore, ...
                out_bottom.fr.ratio, out_bottom.fr.zscore, out_bottom.fr_ws.ratio, out_bottom.fr_ws.zscore] = deal(nan(1,length(fbands)));

            [out_binwise_top, out_binwise_top.lat, out_binwise_top.lat_ws, out_binwise_top.fr, out_binwise_top.fr_ws] = deal([]);
            out_binwise_top.nshufs = nshufs;
            out_binwise_top.bin_count = nbins;
            [out_binwise_top.lat.bin, out_binwise_top.lat_ws.bin, ...
                out_binwise_top.fr.bin, out_binwise_top.fr_ws.bin, ...
                out_binwise_top.fr.sd, out_binwise_top.fr_ws.sd] = deal(nan(length(fbands),nbins));
            [out_binwise_top.lat.ratio, out_binwise_top.lat.zscore, ...
                out_binwise_top.lat_ws.ratio, out_binwise_top.lat_ws.zscore, ...
                out_binwise_top.fr.ratio, out_binwise_top.fr.zscore, out_binwise_top.fr_ws.ratio, out_binwise_top.fr_ws.zscore] = deal(nan(1,length(fbands)));
            
            [out_binwise_bottom, out_binwise_bottom.lat, out_binwise_bottom.lat_ws, ...
                out_binwise_bottom.fr, out_binwise_bottom.fr_ws] = deal([]);
            out_binwise_bottom.nshufs = nshufs;
            out_binwise_bottom.bin_count = nbins;
            [out_binwise_bottom.lat.bin, out_binwise_bottom.lat_ws.bin, ...
                out_binwise_bottom.fr.bin, out_binwise_bottom.fr_ws.bin, ...
                out_binwise_bottom.fr.sd, out_binwise_bottom.fr_ws.sd] = deal(nan(length(fbands),nbins));
            [out_binwise_bottom.lat.ratio, out_binwise_bottom.lat.zscore, ...
                out_binwise_bottom.lat_ws.ratio, out_binwise_bottom.lat_ws.zscore, ...
                out_binwise_bottom.fr.ratio, out_binwise_bottom.fr.zscore, ...
                out_binwise_bottom.fr_ws.ratio, out_binwise_bottom.fr_ws.zscore] = deal(nan(1,length(fbands)));

            [out_top_control, out_top_control.lat, out_top_control.lat_ws, ...
                out_top_control.fr, out_top_control.fr_ws] = deal([]);
            out_top_control.nsplits = nsplits;
            out_top_control.nshufs = nshufs;
            out_top_control.bin_count = nbins;
            [out_top_control.lat.bin, out_top_control.lat_ws.bin, ...
                out_top_control.fr.bin, out_top_control.fr_ws.bin, ...
                out_top_control.fr.sd, out_top_control.fr_ws.sd] = deal(nan(nsplits,length(fbands),nbins));
            [out_top_control.lat.ratio, out_top_control.lat.zscore, ...
                out_top_control.lat_ws.ratio, out_top_control.lat_ws.zscore, ...
                out_top_control.fr.ratio, out_top_control.fr.zscore, ...
                out_top_control.fr_ws.ratio, out_top_control.fr_ws.zscore] = deal(nan(nsplits,length(fbands)));

            [out_bottom_control, out_bottom_control.lat, out_bottom_control.lat_ws, ...
                out_bottom_control.fr, out_bottom_control.fr_ws] = deal([]);
            out_bottom_control.nsplits = nsplits;
            out_bottom_control.nshufs = nshufs;
            out_bottom_control.bin_count = nbins;
            [out_bottom_control.lat.bin, out_bottom_control.lat_ws.bin, ...
                out_bottom_control.fr.bin, out_bottom_control.fr_ws.bin, ...
                out_bottom_control.fr.sd, out_bottom_control.fr_ws.sd] = deal(nan(nsplits,length(fbands),nbins));
            [out_bottom_control.lat.ratio, out_bottom_control.lat.zscore, ...
                out_bottom_control.lat_ws.ratio, out_bottom_control.lat_ws.zscore, ...
                out_bottom_control.fr.ratio, out_bottom_control.fr.zscore, ...
                out_bottom_control.fr_ws.ratio, out_bottom_control.fr_ws.zscore] = deal(nan(nsplits,length(fbands)));

            % Since all calcluations need to happen differently for each
            % frequency band separately, start a loop here
            for iF = 1:length(fbands)
                % Divide trials into two groups based on power (top half vs
                % bottom half)
                keep_trials = mean_power(iF,:) >= power_thresh(iF);
                top_trials = keep_trials(goodTrials(1):goodTrials(2));
                bottom_trials = ~keep_trials(goodTrials(1):goodTrials(2));

                % Divide trials into two groups based on power in each bin separately
                % (top half vs bottom half)
                this_trial_idx = (goodTrials(1):goodTrials(2));
                this_phase = causal_phase(iF,this_trial_idx);
                [~, ~, this_bin] = histcounts(this_phase, phase_bins);
                clear this_trial_idx

                [top_binwise_trials, bottom_binwise_trials] = deal(false(size(top_trials)));
                [top_control_trials, bottom_control_trials] = deal(false(nsplits,length(top_trials)));
                for iBin = 1:nbins
                    this_pow_thresh = prctile(mean_power(iF,this_bin==iBin),trial_thresh);
                    top_binwise_trials(this_bin==iBin) = mean_power(iF,this_bin==iBin) >= this_pow_thresh;
                    bottom_binwise_trials(this_bin==iBin) = mean_power(iF,this_bin==iBin) < this_pow_thresh;
                    % Split each bin into two halves randomly instead on
                    % the basis of power
                    this_bin_trial = find(this_bin==iBin);
                    mid = ceil(length(this_bin_trial)/2);
                    for iSplit = 1:nsplits
                        this_perm = randperm(length(this_bin_trial));
                        % choose half of these
                        top_control_trials(iSplit,this_bin_trial(this_perm(1:mid))) = true;
                        bottom_control_trials(iSplit,this_bin_trial(this_perm(mid+1:end))) = true;
                    end
                end

                top_phase = this_phase(top_trials);
                bottom_phase = this_phase(bottom_trials);
                top_binwise_phase = this_phase(top_binwise_trials);
                bottom_binwise_phase = this_phase(bottom_binwise_trials);    

                % generate shuffle indices
                top_shuf_idx = zeros(nshufs, sum(top_trials));
                bottom_shuf_idx = zeros(nshufs, sum(bottom_trials));
                top_binwise_shuf_idx = zeros(nshufs, sum(top_binwise_trials));
                bottom_binwise_shuf_idx = zeros(nshufs, sum(bottom_binwise_trials));
                top_control_shuf_idx = zeros(nshufs, sum(top_control_trials(1,:)));
                bottom_control_shuf_idx = zeros(nshufs, sum(bottom_control_trials(1,:)));

                for iShuf = 1:nshufs
                    top_shuf_idx(iShuf,:) = randperm(sum(top_trials));
                    bottom_shuf_idx(iShuf,:) = randperm(sum(bottom_trials));
                    top_binwise_shuf_idx(iShuf,:) = randperm(sum(top_binwise_trials));
                    bottom_binwise_shuf_idx(iShuf,:) = randperm(sum(bottom_binwise_trials));
                    top_control_shuf_idx(iShuf,:) = randperm(sum(top_control_trials(1,:)));
                    bottom_control_shuf_idx(iShuf,:) = randperm(sum(bottom_control_trials(1,:)));
                end

                for iSplit = 1:nsplits
                    % Calculate stuff for top control
                    this_trial_idx = (goodTrials(1):goodTrials(2)) & top_control_trials(iSplit,:);
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
                    out_top_control.lat.bin(iSplit,iF,:)= this_lat;
                    out_top_control.lat.ratio(iSplit,iF) = (max(this_lat) - min(this_lat))/(max(this_lat) + min(this_lat));
                    out_top_control.lat_ws.bin(iSplit,iF,:) = this_lat_ws;
                    out_top_control.lat_ws.ratio(iSplit,iF) = (max(this_lat_ws) - min(this_lat_ws))/(max(this_lat_ws) + min(this_lat_ws));
                    out_top_control.fr.bin(iSplit,iF,:) = this_fr;
                    out_top_control.fr.ratio(iSplit,iF) = abs(max(this_fr) - min(this_fr))/(abs(max(this_fr)) + abs(min(this_fr)));
                    out_top_control.fr_ws.bin(iSplit,iF,:) = this_fr_ws;
                    out_top_control.fr_ws.ratio(iSplit,iF) = abs(max(this_fr_ws) - min(this_fr_ws))/(abs(max(this_fr_ws)) + abs(min(this_fr_ws)));
                    out_top_control.fr.sd(iSplit,iF,:) = this_fr_sd;
                    out_top_control.fr_ws.sd(iSplit,iF,:) = this_fr_ws_sd;
    
                    % Generate shuffles
                    [shuf_lat_ratio, shuf_lat_ws_ratio, shuf_fr_ratio, shuf_fr_ws_ratio] = deal(zeros(nshufs, 1));
                    for iShuf = 1:nshufs
                        % shuffle the bin indices
                        shuf_bin = this_bin(top_control_shuf_idx(iShuf,:));
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
                        out_top_control.lat.shufs(iSplit,iF,iShuf,:) = shuf_lat'; 
                        out_top_control.lat_ws.shufs(iSplit,iF,iShuf,:) = shuf_lat_ws';
                        out_top_control.fr.shufs(iSplit,iF,iShuf,:) = shuf_fr';
                        out_top_control.fr_ws.shufs(iSplit,iF,iShuf,:) = shuf_fr_ws';
                    end
                    % Saving Z-scored values
                    out_top_control.lat.zscore(iSplit,iF) = (out_top_control.lat.ratio(iSplit,iF) - mean(shuf_lat_ratio))/std(shuf_lat_ratio);
                    out_top_control.lat_ws.zscore(iSplit,iF) = (out_top_control.lat_ws.ratio(iSplit,iF) - mean(shuf_lat_ws_ratio))/std(shuf_lat_ws_ratio);
                    out_top_control.fr.zscore(iSplit,iF) = (out_top_control.fr.ratio(iSplit,iF) - mean(shuf_fr_ratio))/std(shuf_fr_ratio);
                    out_top_control.fr_ws.zscore(iSplit,iF) = (out_top_control.fr_ws.ratio(iSplit,iF) - mean(shuf_fr_ws_ratio))/std(shuf_fr_ws_ratio);

                    % Calculate stuff for bottom control
                    this_trial_idx = (goodTrials(1):goodTrials(2)) & bottom_control_trials(iSplit,:);
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
                    out_bottom_control.lat.bin(iSplit,iF,:)= this_lat;
                    out_bottom_control.lat.ratio(iSplit,iF) = (max(this_lat) - min(this_lat))/(max(this_lat) + min(this_lat));
                    out_bottom_control.lat_ws.bin(iSplit,iF,:) = this_lat_ws;
                    out_bottom_control.lat_ws.ratio(iSplit,iF) = (max(this_lat_ws) - min(this_lat_ws))/(max(this_lat_ws) + min(this_lat_ws));
                    out_bottom_control.fr.bin(iSplit,iF,:) = this_fr;
                    out_bottom_control.fr.ratio(iSplit,iF) = abs(max(this_fr) - min(this_fr))/(abs(max(this_fr)) + abs(min(this_fr)));
                    out_bottom_control.fr_ws.bin(iSplit,iF,:) = this_fr_ws;
                    out_bottom_control.fr_ws.ratio(iSplit,iF) = abs(max(this_fr_ws) - min(this_fr_ws))/(abs(max(this_fr_ws)) + abs(min(this_fr_ws)));
                    out_bottom_control.fr.sd(iSplit,iF,:) = this_fr_sd;
                    out_bottom_control.fr_ws.sd(iSplit,iF,:) = this_fr_ws_sd;
    
                    % Generate shuffles
                    [shuf_lat_ratio, shuf_lat_ws_ratio, shuf_fr_ratio, shuf_fr_ws_ratio] = deal(zeros(nshufs, 1));
                    for iShuf = 1:nshufs
                        % shuffle the bin indices
                        shuf_bin = this_bin(bottom_control_shuf_idx(iShuf,:));
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
                        out_bottom_control.lat.shufs(iSplit,iF,iShuf,:) = shuf_lat'; 
                        out_bottom_control.lat_ws.shufs(iSplit,iF,iShuf,:) = shuf_lat_ws';
                        out_bottom_control.fr.shufs(iSplit,iF,iShuf,:) = shuf_fr';
                        out_bottom_control.fr_ws.shufs(iSplit,iF,iShuf,:) = shuf_fr_ws';
                    end
                    % Saving Z-scored values
                    out_bottom_control.lat.zscore(iSplit,iF) = (out_bottom_control.lat.ratio(iSplit,iF) - mean(shuf_lat_ratio))/std(shuf_lat_ratio);
                    out_bottom_control.lat_ws.zscore(iSplit,iF) = (out_bottom_control.lat_ws.ratio(iSplit,iF) - mean(shuf_lat_ws_ratio))/std(shuf_lat_ws_ratio);
                    out_bottom_control.fr.zscore(iSplit,iF) = (out_bottom_control.fr.ratio(iSplit,iF) - mean(shuf_fr_ratio))/std(shuf_fr_ratio);
                    out_bottom_control.fr_ws.zscore(iSplit,iF) = (out_bottom_control.fr_ws.ratio(iSplit,iF) - mean(shuf_fr_ws_ratio))/std(shuf_fr_ws_ratio);
                end

                clear this_count

                % Calculate stuff for top half
                this_trial_idx = (goodTrials(1):goodTrials(2)) & top_trials;
                this_phase = causal_phase(iF,this_trial_idx);
                [top_count, out_top] = doStuff2(out_top, od, iF, this_trial_idx, ...
                    top_shuf_idx, this_phase, phase_bins, nshufs, nbins);
 
                % Calculate stuff for bottom half
                this_trial_idx = (goodTrials(1):goodTrials(2)) & bottom_trials;
                this_phase = causal_phase(iF,this_trial_idx);
                [bottom_count, out_bottom] = doStuff2(out_bottom, od, iF, this_trial_idx, ...
                    bottom_shuf_idx, this_phase, phase_bins, nshufs, nbins);

                % Calculate stuff for binwise top half
                this_trial_idx = (goodTrials(1):goodTrials(2)) & top_binwise_trials;
                this_phase = causal_phase(iF,this_trial_idx);
                [top_binwise_count, out_binwise_top] = doStuff2(out_binwise_top, od, iF, this_trial_idx, ...
                    top_binwise_shuf_idx, this_phase, phase_bins, nshufs, nbins);

                % Calculate stuff for binwise bottom half
                this_trial_idx = (goodTrials(1):goodTrials(2)) & bottom_binwise_trials;
                this_phase = causal_phase(iF,this_trial_idx);
                [bottom_binwise_count, out_binwise_bottom] = doStuff2(out_binwise_bottom, od, iF, this_trial_idx, ...
                    bottom_binwise_shuf_idx, this_phase, phase_bins, nshufs, nbins);

                clear this_trial_idx this_perm this_bin_trial
 
                ax = subplot(4,4,(iF-1)*4+1);
                bar(ax,x_ticks,[top_count/sum(top_count);bottom_count/sum(bottom_count)],2)
                ax.Title.String = sprintf('Stim Phase distribution for overall split');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
                ax.YLabel.String = sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2));
                ax.XAxis.TickDirection = 'out';
                legend({'Top', 'Bottom'},'FontSize', 12, 'Location', 'bestoutside')
                
                ax = subplot(4,4,(iF-1)*4+2);
                bar(ax,x_ticks,[out_top.fr.bin(iF,:);out_bottom.fr.bin(iF,:)],2);
                ax.Title.String = sprintf('Firing Rate for overall split');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
    
                ax = subplot(4,4,(iF-1)*4+3);
                bar(ax,x_ticks,[top_binwise_count/sum(top_binwise_count);bottom_binwise_count/sum(bottom_binwise_count)],2)
                ax.Title.String = sprintf('Stim Phase distribution for binwise split');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
                ax.XAxis.TickDirection = 'out';
                legend({'Top', 'Bottom'},'FontSize', 12, 'Location', 'bestoutside')
                
                ax = subplot(4,4,(iF-1)*4+4);
                bar(ax,x_ticks,[out_binwise_top.fr.bin(iF,:);out_binwise_bottom.fr.bin(iF,:)],2);
                ax.Title.String = sprintf('Firing Rate for binwise split');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];   
            end

            fn_prefix = strcat(fn_prefix, '_split_by_power');
            sgtitle(fn_prefix, 'Interpreter', 'none');
            % Save variables and figure
            save(strcat(fn_prefix, '_phase_response_',string(bin_counts(iN)),'_bins'), 'out_top', 'out_bottom', ...
                'out_binwise_top', 'out_binwise_bottom', 'out_top_control', 'out_bottom_control');
            print(this_fig, '-dpng', strcat(fn_prefix, '_phase_response_hist_', string(bin_counts(iN)),'_bins'));
            close all;
        end
    end
end

function [out_count, out_struct] = doStuff2(in_struct1, in_struct2, iF, ...
    this_trial_idx, this_shuf_idx, this_phase, phase_bins, nshufs, nbins)
    out_struct = in_struct1;

    all_lat = in_struct2.trial_stim.latency(this_trial_idx);
    all_lat_ws = in_struct2.trial_stim.latency_wo_stim(this_trial_idx);
    % Now subtracting baseline firing rate
    all_bfr = in_struct2.trial_stim.bfr(this_trial_idx);
    all_fr = in_struct2.trial_stim.fr(this_trial_idx) - all_bfr;
    all_fr_ws = in_struct2.trial_stim.fr_wo_stim(this_trial_idx) - all_bfr;

%     this_phase = causal_phase(iF,this_trial_idx);
    [out_count, ~, this_bin] = histcounts(this_phase, phase_bins);
    [this_lat, this_lat_ws, this_fr, this_fr_ws, this_fr_sd, this_fr_ws_sd] = deal(zeros(size(out_count)));
    for iB = 1:nbins
       this_lat(iB) = sum(~isnan(all_lat(this_bin==iB)))/out_count(iB);                   
       this_lat_ws(iB) = sum(~isnan(all_lat_ws(this_bin==iB)))/out_count(iB);
       this_fr(iB) = mean(all_fr(this_bin==iB));
       this_fr_sd(iB) = std(all_fr(this_bin==iB));
       this_fr_ws(iB) = mean(all_fr_ws(this_bin==iB));
       this_fr_ws_sd(iB) = std(all_fr_ws(this_bin==iB));
    end

    % Using (max-min)/(max+min) as the ratio
    out_struct.lat.bin(iF,:)= this_lat;
    out_struct.lat.ratio(iF) = (max(this_lat) - min(this_lat))/(max(this_lat) + min(this_lat));
    out_struct.lat_ws.bin(iF,:) = this_lat_ws;
    out_struct.lat_ws.ratio(iF) = (max(this_lat_ws) - min(this_lat_ws))/(max(this_lat_ws) + min(this_lat_ws));
    out_struct.fr.bin(iF,:) = this_fr;
    out_struct.fr.ratio(iF) = abs(max(this_fr) - min(this_fr))/(abs(max(this_fr)) + abs(min(this_fr)));
    out_struct.fr_ws.bin(iF,:) = this_fr_ws;
    out_struct.fr_ws.ratio(iF) = abs(max(this_fr_ws) - min(this_fr_ws))/(abs(max(this_fr_ws)) + abs(min(this_fr_ws)));
    out_struct.fr.sd(iF,:) = this_fr_sd;
    out_struct.fr_ws.sd(iF,:) = this_fr_ws_sd;

    % Generate shuffles
    [shuf_lat_ratio, shuf_lat_ws_ratio, shuf_fr_ratio, shuf_fr_ws_ratio] = deal(zeros(nshufs, 1));
    for iShuf = 1:nshufs
        % shuffle the bin indices
        shuf_bin = this_bin(this_shuf_idx(iShuf,:));
        [shuf_lat, shuf_lat_ws, shuf_fr, shuf_fr_ws] = deal(zeros(1,nbins));
        for iB = 1:nbins
           shuf_lat(iB) = sum(~isnan(all_lat(shuf_bin==iB)))/out_count(iB);
           shuf_lat_ws(iB) = sum(~isnan(all_lat_ws(shuf_bin==iB)))/out_count(iB);
           shuf_fr(iB) = sum(all_fr(shuf_bin==iB))/out_count(iB);
           shuf_fr_ws(iB) = sum(all_fr_ws(shuf_bin==iB))/out_count(iB);
        end
        shuf_lat_ratio(iShuf) = (max(shuf_lat) - min(shuf_lat))/(max(shuf_lat) + min(shuf_lat));
        shuf_lat_ws_ratio(iShuf) = (max(shuf_lat_ws) - min(shuf_lat_ws))/(max(shuf_lat_ws) + min(shuf_lat_ws));
        shuf_fr_ratio(iShuf) = abs(max(shuf_fr) - min(shuf_fr))/(abs(max(shuf_fr)) + abs(min(shuf_fr)));
        shuf_fr_ws_ratio(iShuf) = abs(max(shuf_fr_ws) - min(shuf_fr_ws))/(abs(max(shuf_fr_ws)) + abs(min(shuf_fr_ws)));
        % Also save shufs for significance calculation and description in figures later
        out_struct.lat.shufs(iF,iShuf,:) = shuf_lat'; 
        out_struct.lat_ws.shufs(iF,iShuf,:) = shuf_lat_ws';
        out_struct.fr.shufs(iF,iShuf,:) = shuf_fr';
        out_struct.fr_ws.shufs(iF,iShuf,:) = shuf_fr_ws';
    end
    % Saving Z-scored values
    out_struct.lat.zscore(iF) = (out_struct.lat.ratio(iF) - mean(shuf_lat_ratio))/std(shuf_lat_ratio);
    out_struct.lat_ws.zscore(iF) = (out_struct.lat_ws.ratio(iF) - mean(shuf_lat_ws_ratio))/std(shuf_lat_ws_ratio);
    out_struct.fr.zscore(iF) = (out_struct.fr.ratio(iF) - mean(shuf_fr_ratio))/std(shuf_fr_ratio);
    out_struct.fr_ws.zscore(iF) = (out_struct.fr_ws.ratio(iF) - mean(shuf_fr_ws_ratio))/std(shuf_fr_ws_ratio);
end