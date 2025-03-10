%% Script to generate the relationships between spiking and LFP Phase
% Assumes that stim_phases.mat and *stim_response.mat already exist in each folder
top_dir = 'data\'
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 12%:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1%:length(sid)
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
    fbands = {[2 5], [6 10], [12,28], [30 55]};
    c_list = {'red', 'blue','magenta','cyan'};
    nshufs = 1000;
    bin_counts = 5; % [5, 7, 10, 15, 20];
    
    % Load the phases at stim_on in various frequency bands
    load('stim_phases.mat');

    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        fn_prefix = 'M295-2022-01-06-TT06_4'
        % Load the stim_responses
        load(strcat(fn_prefix, '_stim_response.mat'));
        goodTrials = ExpKeys.goodTrials(iC,:);

        all_lat = od.trial_stim.latency(goodTrials(1):goodTrials(2));
        all_lat_ws = od.trial_stim.latency_wo_stim(goodTrials(1):goodTrials(2));
        % Now subtracting baseline firing rate
        all_bfr = od.trial_stim.bfr(goodTrials(1):goodTrials(2));
        all_fr = od.trial_stim.fr(goodTrials(1):goodTrials(2)) - all_bfr;
        all_fr_ws = od.trial_stim.fr_wo_stim(goodTrials(1):goodTrials(2)) - all_bfr;
        
        % Repeat for each bin count
        for iN = 1:length(bin_counts)   
            nbins = bin_counts(iN);
            phase_bins = -pi:2*pi/nbins:pi;
            x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));
            % Results field to be saved at the end
            [out, out.lat, out.lat_ws, out.fr, out.fr_ws] = deal([]);
            out.nshufs = nshufs;
            out.bin_count = nbins;
            out.overall_response = sum(~isnan(all_lat))/length(all_lat);
            out.overall_response_ws = sum(~isnan(all_lat_ws))/length(all_lat_ws);
            out.bfr = all_bfr;

            [out.lat.bin, out.lat_ws.bin, out.fr.bin, out.fr_ws.bin, out.fr.sd, out.fr_ws.sd] = deal(nan(length(fbands),nbins));
            [out.lat.ratio, out.lat.zscore, out.lat_ws.ratio, out.lat_ws.zscore, ...
                out.fr.ratio, out.fr.zscore, out.fr_ws.ratio, out.fr_ws.zscore] = deal(nan(1,length(fbands)));
            % Saving all the shufs for later
            [out.lat.shufs, out.lat_ws.shufs, out.fr.shufs, out.fr_ws.shufs] = deal(nan(length(fbands),nshufs,nbins)); 
    
            % Generate indices for shuffling
            shuf_idx = zeros(nshufs, length(all_lat));
            for iShuf = 1:nshufs
                shuf_idx(iShuf,:) = randperm(length(all_lat));
            end
    
            this_fig = figure('WindowState', 'maximized'); 
            for iF = 1:length(fbands)
                this_phase = causal_phase(iF,goodTrials(1):goodTrials(2));
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

                out.lat.bin(iF,:)= this_lat;
                out.lat.ratio(iF) = (max(this_lat) - min(this_lat))/(max(this_lat) + min(this_lat));
                out.lat_ws.bin(iF,:) = this_lat_ws;
                out.lat_ws.ratio(iF) = (max(this_lat_ws) - min(this_lat_ws))/(max(this_lat_ws) + min(this_lat_ws));
                out.fr.bin(iF,:) = this_fr;
                out.fr.ratio(iF) = abs((max(this_fr) - min(this_fr))/(max(this_fr) + min(this_fr)));
                out.fr_ws.bin(iF,:) = this_fr_ws;
                out.fr_ws.ratio(iF) = abs((max(this_fr_ws) - min(this_fr_ws))/(max(this_fr_ws) + min(this_fr_ws)));
                out.fr.sd(iF,:) = this_fr_sd;
                out.fr_ws.sd(iF,:) = this_fr_ws_sd;
                % Generate shuffles
                [shuf_lat_ratio, shuf_lat_ws_ratio, shuf_fr_ratio, shuf_fr_ws_ratio] = deal(zeros(nshufs, 1));
                for iShuf = 1:nshufs
                    % shuffle the bin indices
                    shuf_bin = this_bin(shuf_idx(iShuf,:));
                    [shuf_lat, shuf_lat_ws, shuf_fr, shuf_fr_ws] = deal(zeros(1,nbins));
                    for iB = 1:nbins
                       shuf_lat(iB) = sum(~isnan(all_lat(shuf_bin==iB)))/this_count(iB);
                       shuf_lat_ws(iB) = sum(~isnan(all_lat_ws(shuf_bin==iB)))/this_count(iB);
                       shuf_fr(iB) = sum(all_fr(shuf_bin==iB))/this_count(iB);
                       shuf_fr_ws(iB) = sum(all_fr_ws(shuf_bin==iB))/this_count(iB);
                    end
                    shuf_lat_ratio(iShuf) = (max(shuf_lat) - min(shuf_lat))/(max(shuf_lat) + min(shuf_lat));
                    shuf_lat_ws_ratio(iShuf) = (max(shuf_lat_ws) - min(shuf_lat_ws))/(max(shuf_lat_ws) + min(shuf_lat_ws));
                    shuf_fr_ratio(iShuf) = abs((max(shuf_fr) - min(shuf_fr))/(max(shuf_fr) + min(shuf_fr)));
                    shuf_fr_ws_ratio(iShuf) = abs((max(shuf_fr_ws) - min(shuf_fr_ws))/(max(shuf_fr_ws) + min(shuf_fr_ws)));
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
            
            ax = subplot(5,5,21);
            ax.FontSize = 25;
            hold off
            text(0, 0.75, sprintf('Overall response ratio: % .2f', out.overall_response))
            text(0, 0.25, sprintf('Overall response w/o stim ratio: % .2f', out.overall_response_ws));
            ax.Box = 'off';
            ax.Visible = 'off';
    
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
            
            sgtitle(fn_prefix, 'Interpreter', 'none');
            % Save variables and figure
            save(strcat(fn_prefix, '_phase_response_',string(bin_counts(iN)),'_bins'), 'out');
            print(this_fig, '-dpng', strcat(fn_prefix, '_phase_response_hist_', string(bin_counts(iN)),'_bins'));
%             savefig(this_fig, strcat(fn_prefix, '_phase_response_hist_', string(nbins),'_bins'));
            close;
        end
    end
end
