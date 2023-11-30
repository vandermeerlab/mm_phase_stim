%% Script to generate the relationships between spiking and LFP Phase
% Assumes that stim_phases.mat and *stim_response.mat already exist in each folder
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
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
    
    fbands = {[2 5], [6 10], [30 55]};
    c_list = {'red', 'blue','green'};
    nshufs = 1000;
    stim_binc = 5;
    ns_binc = 25;
    phase_bins = -pi:2*pi/stim_binc:pi;
    ns_bins =  -pi:2*pi/ns_binc:pi;
    x_ticks = 0.5*(phase_bins(1:end-1) + phase_bins(2:end));
    ns_ticks = 0.5*(ns_bins(1:end-1) + ns_bins(2:end));

    % Load the phases at stim_on in various frequency bands
    load('stim_phases.mat');
    % Get rid of all the 3rd band stuff IF there are 4 bands
    if size(causal_phase, 1) == 4 causal_phase(3,:) = []; end

    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        % Load the stim_responses
        load(strcat(fn_prefix, '_stim_response.mat'));
        goodTrials = ExpKeys.goodTrials(iC,:);
        
        % Now subtracting baseline firing rate
        all_bfr = od.trial_stim.bfr(goodTrials(1):goodTrials(2));
        all_fr = od.trial_stim.fr(goodTrials(1):goodTrials(2)) - all_bfr;
        
        % Load the uncorrected stuff
        load(strcat(fn_prefix, '_phase_response_5_bins.mat'));

        % Load non_stim lookup at various frequency bands
        load(strcat(fn_prefix, '_nonstim_spk_phases.mat'));

        ns_lookup = phase_out.ns_lookup;
        ns_fr = phase_out.binned_fr;
        ns_spk_count = phase_out.binned_spk_count;
        corrected_fr_bin = zeros(length(fbands), stim_binc);
        corrected_fr_shufs = zeros(length(fbands), nshufs, stim_binc);
        [corrected_fr_r, corrected_fr_z, corrected_fr_z_oldshufs] =  ...
            deal(zeros(length(fbands),1));

        this_fig = figure('WindowState', 'maximized');
        for iF = 1:length(fbands)
            
            this_phase = causal_phase(iF,goodTrials(1):goodTrials(2));
            [~,~,this_lookup_bin] = histcounts(this_phase, ns_bins);
            
            % Apply corrections
            for iB = 1:ns_binc
                idx = find(this_lookup_bin == iB);
                all_fr(idx) = all_fr(idx) - ns_lookup(iF,iB);
            end
            
            [this_count, ~, this_bin] = histcounts(this_phase, phase_bins);
            for iB = 1:stim_binc
               corrected_fr_bin(iF,iB) = mean(all_fr(this_bin==iB));
            end
            corrected_fr_r(iF) = (max(corrected_fr_bin(iF,:)) - min(corrected_fr_bin(iF,:)))/...
                (max(corrected_fr_bin(iF,:)) + min(corrected_fr_bin(iF,:)));

            % Generate indices for shuffling
            shuf_idx = zeros(nshufs, length(all_fr));
            for iShuf = 1:nshufs
                shuf_idx(iShuf,:) = randperm(length(all_fr));
            end

            % Generate new shuffles
            shuf_fr_ratio = zeros(nshufs, 1);
            for iShuf = 1:nshufs
                % shuffle the bin indices
                shuf_bin = this_bin(shuf_idx(iShuf,:));
                shuf_fr = deal(zeros(1,stim_binc));
                for iB = 1:stim_binc
                   shuf_fr(iB) = sum(all_fr(shuf_bin==iB))/this_count(iB);
                end
                shuf_fr_ratio(iShuf) = (max(shuf_fr) - min(shuf_fr))/(max(shuf_fr) + min(shuf_fr));
                corrected_fr_shufs(iF,iShuf,:) = shuf_fr;
            end
            corrected_fr_z(iF) = (corrected_fr_r(iF) - mean(shuf_fr_ratio))/std(shuf_fr_ratio);

            % Calculate shuf ratio based on older shuffles
            q0 = squeeze(out.fr.shufs(iF,:,:));
            q1 = max(q0, [] , 2);
            q2 = min(q0, [], 2);
            q3 = (q1 - q2)./(q1 + q2); 
            corrected_fr_z_oldshufs(iF) = (corrected_fr_r(iF) - mean(q3))/std(q3);
            clear q0 q1 q2 q3 this_bin this_count

            % Plot the Spike-phase distribution
            ax = subplot(4,5,(iF-1)*5+1);
            bar(ax,ns_ticks,ns_spk_count(iF,:),1,c_list{iF});
            xlabel('Phase Bin')
            ylabel('Spike Count')
            ax.Title.String = sprintf('%d - %d Hz', ...
                fbands{iF}(1), fbands{iF}(2));
            ax.XLim = [-3.25 3.25];

            % Plot the FR_phase distribution
            ax = subplot(4,5,(iF-1)*5+2);
            bar(ax,ns_ticks,ns_fr(iF,:),1,c_list{iF});
            xlabel('Phase Bin')
            ylabel('FR')
            ax.XLim = [-3.25 3.25];

            % Plot the lookup
            ax = subplot(4,5,(iF-1)*5+3);
            bar(ax,ns_ticks,ns_lookup(iF,:),1,c_list{iF});
            ax.Title.String = 'NonStim Lookup';
            xlabel('Phase Bin')
            ylabel('{\Delta} FR')
            ax.XLim = [-3.25 3.25];
            
            ax = subplot(4,5,(iF-1)*5+4);
            xlabel('Phase Bin')
            ylabel('{\Delta} FR')
            bar(ax,x_ticks,out.fr.bin(iF,:),1,c_list{iF});
            ax.Title.String = 'Original {\Delta FR}';
            ax.XLim = [-3.25 3.25];

            ax = subplot(4,5,(iF-1)*5+5);
            xlabel('Phase Bin')
            ylabel('{\Delta} FR')
            bar(ax,x_ticks,corrected_fr_bin(iF,:),1,c_list{iF});
            ax.Title.String = 'Corrected {\Delta FR}';
            ax.XLim = [-3.25 3.25];  
        end

        ax = subplot(4,5,19);
        hold on
        for iF = 1:length(fbands)
            scatter(mean(fbands{iF}), out.fr.zscore(iF), 'MarkerFaceColor', ...
                c_list{iF}, 'MarkerEdgeColor', c_list{iF});
        end
        ax.XLim = [0 60];
        ax.YLim = [-3 12];
        ax.XAxis.Label.String = 'Freq (Hz)';
        ax.YAxis.Label.String = 'Z-score';

        ax = subplot(4,5,20);
        hold on
        for iF = 1:length(fbands)
            scatter(mean(fbands{iF}), corrected_fr_z(iF), 'MarkerFaceColor', ...
                c_list{iF}, 'MarkerEdgeColor', c_list{iF});
        end
        ax.XLim = [0 60];
        ax.YLim = [-3 12];
        ax.XAxis.Label.String = 'Freq (Hz)';
        ax.YAxis.Label.String = 'Z-score';

        sgtitle(fn_prefix, 'Interpreter', 'none');
        % Save variables and figure
        save(strcat(fn_prefix, '_corrected_phase_response_5_bins'), ...
            'corrected_fr_bin', 'corrected_fr_r', 'corrected_fr_z', 'corrected_fr_shufs', 'corrected_fr_z_oldshufs');
        print(this_fig, '-dpng', strcat(fn_prefix, '_corrected_phase_response_hist_5_bins'));
        close;
    end
end