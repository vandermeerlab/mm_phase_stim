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

        % Generate indices for shuffling
            
        % Treating trial_thresh = 0 as a separate case to keep the results
        % from before the same as before
        [freqwise_trials, freqwise_bottom] = deal(cell(size(fbands)));
        if trial_thresh == 0 % include all trials
            for iF = 1:length(fbands)
                freqwise_trials{iF} = true(1,diff(goodTrials)+1);
            end
        else
            for iF = 1:length(fbands)
                keep_trials = mean_power(iF,:) >= power_thresh(iF);
                freqwise_trials{iF} = keep_trials(goodTrials(1):goodTrials(2));
                freqwise_bottom{iF} = ~keep_trials(goodTrials(1):goodTrials(2));
                clear keep_trials
            end
        end
        
        % Repeat for each bin count
        for iN = 1:length(bin_counts)
            this_fig = figure('WindowState', 'maximized');
            nbins = bin_counts(iN);
            phase_bins = -pi:2*pi/nbins:pi;
            x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));

            % Since all calcluations need to happen differently for each
            % frequency band separately, start a loop here
            for iF = 1:length(fbands)
                this_trial_idx = (goodTrials(1):goodTrials(2)) & freqwise_trials{iF};
                this_bottom_idx = (goodTrials(1):goodTrials(2)) & freqwise_bottom{iF};

                this_phase = causal_phase(iF,this_trial_idx);
                [this_count, ~, this_bin] = histcounts(this_phase, phase_bins);

                this_bphase = causal_phase(iF,this_bottom_idx);
                [this_bcount, ~, this_b_bin] = histcounts(this_bphase, phase_bins);

     
                % Plot
                ax = subplot(4,2,(iF-1)*2+1);
                bar(ax,x_ticks,this_count/sum(this_count),1,c_list{iF});
                ax.Title.String = sprintf('Stim Phase distribution for top 50% trials');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
                ax.YLabel.String = sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2));

                ax = subplot(4,2,(iF-1)*2+2);
                bar(ax,x_ticks,this_bcount/sum(this_bcount),1,c_list{iF});
                ax.Title.String = sprintf('Stim Phase distribution for bottom 50% trials');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
                ax.YLabel.String = sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2));     
            end
                 
            fn_prefix = strcat(fn_prefix, '_sanityCheck');
            sgtitle(fn_prefix, 'Interpreter', 'none');
            % Save variables and figure
            print(this_fig, '-dpng', strcat(fn_prefix, '_', string(bin_counts(iN)),'_bins'));
            close;
        end
    end
end
