%% Script to generate the relationships between spiking and LFP Phase
% Assumes that stim_phases.mat and *stim_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
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
    fbands = {[2 5], [6 10], [12 30], [30 55]};
    c_list = {'cyan', 'red','magenta', 'green'};
    nshufs = 100;
    bin_counts = 5; %[5, 7, 10, 15, 20]; % in msec
    
    % Load the phases at stim_on in various frequency bands
    load('control_stim_phases.mat');

    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        % Load the stim_responses
        load(strcat(fn_prefix, '_stim_response.mat'));
        goodTrials = ExpKeys.goodTrials(iC,:);

        all_lat1 = od.stim_pre10.latency(goodTrials(1):goodTrials(2));
        all_lat2 = od.stim_pre250.latency(goodTrials(1):goodTrials(2));
        % Now subtracting baseline firing rate
        all_bfr1 = od.stim_pre10.bfr(goodTrials(1):goodTrials(2));
        all_bfr2 = od.stim_pre250.bfr(goodTrials(1):goodTrials(2));
        all_fr1 = od.stim_pre10.fr(goodTrials(1):goodTrials(2)) - all_bfr1;
        all_fr2 = od.stim_pre250.fr(goodTrials(1):goodTrials(2)) - all_bfr2;
        
        % Repeat for each bin count
        for iN = 1:length(bin_counts)   
            nbins = bin_counts(iN);
            phase_bins = -pi:2*pi/nbins:pi;
            x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));
            % Results field to be saved at the end
            [out10, out10.lat, out10.fr, out250, out250.lat, out250.fr] = deal([]);
            out10.nshufs = nshufs;
            out10.bin_count = nbins;
            out10.overall_response = sum(~isnan(all_lat1))/length(all_lat1);
            out10.bfr = all_bfr1;
            out250.nshufs = nshufs;
            out250.bin_count = nbins;
            out250.overall_response = sum(~isnan(all_lat2))/length(all_lat2);
            out250.bfr = all_bfr2;

            [out10.lat.bin, out10.fr.bin, out250.lat.bin, out250.fr.bin,] = deal(nan(4,nbins));
            [out10.lat.ratio, out10.lat.zscore, out10.fr.ratio, out10.fr.zscore, ...
                out250.lat.ratio, out250.lat.zscore, out250.fr.ratio, out250.fr.zscore,] = deal(nan(1,4));
    
            % Generate indices for shuffling
            shuf_idx  = zeros(nshufs, length(all_lat1));
            for iShuf = 1:nshufs
                shuf_idx(iShuf,:) = randperm(length(all_lat1));
            end
    
            this_fig = figure('WindowState', 'maximized'); 
            for iF = 1:length(fbands)
                this_phase1 = causal_phase_pre10(iF,goodTrials(1):goodTrials(2));
                this_phase2 = causal_phase_pre250(iF,goodTrials(1):goodTrials(2));
                [this_count1, ~, this_bin1] = histcounts(this_phase1, phase_bins);
                [this_count2, ~, this_bin2] = histcounts(this_phase2, phase_bins);
    
                [this_lat1, this_fr1] = deal(zeros(size(this_count1)));
                [this_lat2, this_fr2] = deal(zeros(size(this_count2)));
                for iB = 1:nbins
                   this_lat1(iB) = sum(~isnan(all_lat1(this_bin1==iB)))/this_count1(iB);
                   this_fr1(iB) = sum(all_fr1(this_bin1==iB))/this_count1(iB);
                   this_lat2(iB) = sum(~isnan(all_lat2(this_bin2==iB)))/this_count2(iB);
                   this_fr2(iB) = sum(all_fr2(this_bin2==iB))/this_count2(iB);
                end
                
                % Using (max-min)/(max+min) as the ratio
                out10.lat.bin(iF,:)= this_lat1;
                out10.lat.ratio(iF) = (max(this_lat1) - min(this_lat1))/(max(this_lat1) + min(this_lat1));
                out10.fr.bin(iF,:) = this_fr1;
                out10.fr.ratio(iF) = (max(this_fr1) - min(this_fr1))/(max(this_fr1) + min(this_fr1));
                out250.lat.bin(iF,:)= this_lat2;
                out250.lat.ratio(iF) = (max(this_lat2) - min(this_lat2))/(max(this_lat2) + min(this_lat2));
                out250.fr.bin(iF,:) = this_fr2;
                out250.fr.ratio(iF) = (max(this_fr2) - min(this_fr2))/(max(this_fr2) + min(this_fr2));
    
                % Generate shuffles
                [shuf_lat_ratio1, shuf_fr_ratio1, shuf_lat_ratio2, shuf_fr_ratio2] = deal(zeros(nshufs, 1));
                for iShuf = 1:nshufs
                    % shuffle the bin indices
                    shuf_bin1 = this_bin1(shuf_idx(iShuf,:));
                    shuf_bin2 = this_bin2(shuf_idx(iShuf,:));
                    [shuf_lat1, shuf_fr1, shuf_lat2, shuf_fr2] = deal(zeros(1,nbins));
                    for iB = 1:nbins
                       shuf_lat1(iB) = sum(~isnan(all_lat1(shuf_bin1==iB)))/this_count1(iB);
                       shuf_fr1(iB) = sum(all_fr1(shuf_bin1==iB))/this_count1(iB);
                       shuf_lat2(iB) = sum(~isnan(all_lat2(shuf_bin2==iB)))/this_count2(iB);
                       shuf_fr2(iB) = sum(all_fr2(shuf_bin2==iB))/this_count2(iB);
                    end
                    shuf_lat_ratio1(iShuf) = (max(shuf_lat1) - min(shuf_lat1))/(max(shuf_lat1) + min(shuf_lat1));
                    shuf_fr_ratio1(iShuf) = (max(shuf_fr1) - min(shuf_fr1))/(max(shuf_fr1) + min(shuf_fr1));
                    shuf_lat_ratio2(iShuf) = (max(shuf_lat2) - min(shuf_lat2))/(max(shuf_lat2) + min(shuf_lat2));
                    shuf_fr_ratio2(iShuf) = (max(shuf_fr2) - min(shuf_fr2))/(max(shuf_fr2) + min(shuf_fr2));
                end
                 
                % Saving Z-scored values
                out10.lat.zscore(iF) = (out10.lat.ratio(iF) - mean(shuf_lat_ratio1))/std(shuf_lat_ratio1);
                out10.fr.zscore(iF) = (out10.fr.ratio(iF) - mean(shuf_fr_ratio1))/std(shuf_fr_ratio1);
                out250.lat.zscore(iF) = (out250.lat.ratio(iF) - mean(shuf_lat_ratio2))/std(shuf_lat_ratio2);
                out250.fr.zscore(iF) = (out250.fr.ratio(iF) - mean(shuf_fr_ratio2))/std(shuf_fr_ratio2);
    
                ax = subplot(5,4,(iF-1)*4+1);
                bar(ax,x_ticks,this_count1/sum(this_count1),1,c_list{iF});
                ax.Title.String = sprintf('Stim-10 Phase distribution for %d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
    
                ax = subplot(5,4,(iF-1)*4+2);
                bar(ax,x_ticks,this_fr1,1,c_list{iF});
                ax.Title.String = sprintf('Firing rate (Stim - 10)');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];

                ax = subplot(5,4,(iF-1)*4+3);
                bar(ax,x_ticks,this_count2/sum(this_count2),1,c_list{iF});
                ax.Title.String = sprintf('Stim-250 Phase distribution for %d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
    
                ax = subplot(5,4,(iF-1)*4+4);
                bar(ax,x_ticks,this_fr2,1,c_list{iF});
                ax.Title.String = sprintf('Firing rate (Stim - 250)');
                ax.Title.FontSize = 12;
                ax.XLim = [-3.25 3.25];
            end
                     
            ax = subplot(5,4,18);
            ax.FontSize = 12;
            hold on
            for iF = 1:length(fbands)
                scatter(mean(fbands{iF}), out10.fr.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                    'MarkerEdgeColor', c_list{iF})
            end
            ax.XLim = [0 60];
    %         ax.YLim = [-3 10];
            ax.XAxis.Label.String = 'Freq (Hz)';
            ax.YAxis.Label.String = 'Z-score';
    
            ax = subplot(5,4,20);
            ax.FontSize = 12;
            hold on
            for iF = 1:length(fbands)
                scatter(mean(fbands{iF}), out250.fr.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                    'MarkerEdgeColor', c_list{iF})
            end
            ax.XLim = [0 60];
    %         ax.YLim = [-3 10];
            ax.XAxis.Label.String = 'Freq (Hz)';
            ax.YAxis.Label.String = 'Z-score';
            
            sgtitle(strcat(fn_prefix, ' Control'), 'Interpreter', 'none');
            % Save variables and figure
            save(strcat(fn_prefix, '_control_phase_response_',string(bin_counts(iN)),'_bins'), 'out10', 'out250');
            print(this_fig, '-dpng', strcat(fn_prefix, '_control_phase_response_hist_', string(bin_counts(iN)),'_bins'));
%             savefig(this_fig, strcat(fn_prefix, '_phase_response_hist_', string(nbins),'_bins'));
            close;
        end
    end
end