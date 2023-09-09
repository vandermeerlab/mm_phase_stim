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
    fbands = {[2 5], [6 10], [30 55]};
    c_list = {'red', 'blue', 'green'};
    nshufs = 100;
    bin_counts = 5; %[5, 7, 10, 15, 20]; % in msec
    

    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        % Load the stim_responses
        load(strcat(fn_prefix, '_stim_response.mat'));

        % Load the nonstim_phases
        load(strcat(fn_prefix, '_non_stim_phases.mat'));


        all_lat = od.trial_nonstim.latency;
        % Now subtracting baseline firing rate
        all_bfr = od.trial_nonstim.bfr;
        all_fr = od.trial_nonstim.fr;

        
        % Repeat for each bin count
        for iN = 1:length(bin_counts)   
            nbins = bin_counts(iN);
            phase_bins = -pi:2*pi/nbins:pi;
            x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));
            % Results field to be saved at the end
            [ns_out, ns_out.lat, ns_out.fr] = deal([]);
            ns_out.nshufs = nshufs;
            ns_out.bin_count = nbins;
            ns_out.overall_response = sum(~isnan(all_lat))/length(all_lat);
            ns_out.bfr = all_bfr;

            [ns_out.lat.bin, ns_out.fr.bin] = deal(nan(length(fbands),nbins));
            [ns_out.lat.ratio, ns_out.lat.zscore, ns_out.fr.ratio, ns_out.fr.zscore] = deal(nan(1,length(fbands)));
            [ns_out.lat.shufs, ns_out.fr.shufs] = deal(nan(length(fbands),nshufs,nbins)); 
    
            % Generate indices for shuffling
            shuf_idx  = zeros(nshufs, length(all_lat));
            for iShuf = 1:nshufs
                shuf_idx(iShuf,:) = randperm(length(all_lat));
            end
    
            this_fig = figure('WindowState', 'maximized'); 
            for iF = 1:length(fbands)
                this_phase = nonstim_causal_phase(iF,:);
                [this_count, ~, this_bin] = histcounts(this_phase, phase_bins);
                [this_lat, this_fr] = deal(zeros(size(this_count)));
                for iB = 1:nbins
                   this_lat(iB) = sum(~isnan(all_lat(this_bin==iB)))/this_count(iB);
                   this_fr(iB) = sum(all_fr(this_bin==iB))/this_count(iB);
                end
                
                % Using (max-min)/(max+min) as the ratio
                ns_out.lat.bin(iF,:)= this_lat;
                ns_out.lat.ratio(iF) = (max(this_lat) - min(this_lat))/(max(this_lat) + min(this_lat));
                ns_out.fr.bin(iF,:) = this_fr;
                ns_out.fr.ratio(iF) = (max(this_fr) - min(this_fr))/(max(this_fr) + min(this_fr));
                
                % Generate shuffles
                [shuf_lat_ratio, shuf_fr_ratio] = deal(zeros(nshufs, 1));
                for iShuf = 1:nshufs
                    % shuffle the bin indices
                    shuf_bin = this_bin(shuf_idx(iShuf,:));
                    [shuf_lat, shuf_fr] = deal(zeros(1,nbins));
                    for iB = 1:nbins
                       shuf_lat(iB) = sum(~isnan(all_lat(shuf_bin==iB)))/this_count(iB);
                       shuf_fr(iB) = sum(all_fr(shuf_bin==iB))/this_count(iB);
                    end
                    shuf_lat_ratio(iShuf) = (max(shuf_lat) - min(shuf_lat))/(max(shuf_lat) + min(shuf_lat));
                    shuf_fr_ratio(iShuf) = (max(shuf_fr) - min(shuf_fr))/(max(shuf_fr) + min(shuf_fr));
                    % Also save shufs for significance calculation and description in figures later
                    ns_out.lat.shufs(iF,iShuf,:) = shuf_lat'; 
                    ns_out.fr.shufs(iF,iShuf,:) = shuf_fr';
                end        
                % Saving Z-scored values
                ns_out.lat.zscore(iF) = (ns_out.lat.ratio(iF) - mean(shuf_lat_ratio))/std(shuf_lat_ratio);
                ns_out.fr.zscore(iF) = (ns_out.fr.ratio(iF) - mean(shuf_fr_ratio))/std(shuf_fr_ratio);
    
                % Load Stim-Phase Response
                load(strcat(fn_prefix, '_phase_response_', string(nbins), '_bins.mat'));

                % Load stim Response for side by side plotting
                ax = subplot(4,3,(iF-1)*3+1);
                bar(ax,x_ticks,this_count/sum(this_count),1,c_list{iF});
                ax.Title.String = 'NonStim Phase distribution';
                ax.YLabel.String = sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2));
                ax.Title.FontSize = 20;
                ax.XLim = [-3.25 3.25];
    
                ax = subplot(4,3,(iF-1)*3+2);
                bar(ax,x_ticks,this_fr,1,c_list{iF});
                ax.Title.String = sprintf('FR nonStim');
                ax.Title.FontSize = 20;
                ax.XLim = [-3.25 3.25];
                ax.YLim = [0 10];

                ax = subplot(4,3,(iF-1)*3+3);
                bar(ax,x_ticks,out.fr.bin(iF,:),1,c_list{iF});
                ax.Title.String = sprintf('FR OptoStim');
                ax.Title.FontSize = 20;
                ax.XLim = [-3.25 3.25];
                ax.YLim = [0 100];
            end
                     
            ax = subplot(4,3,11);
            ax.FontSize = 12;
            hold on
            for iF = 1:length(fbands)
                scatter(mean(fbands{iF}), ns_out.fr.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                    'MarkerEdgeColor', c_list{iF})
            end
            ax.XLim = [0 60];
            ax.YLim = [-3 12];
            ax.XAxis.Label.String = 'Freq (Hz)';
            ax.YAxis.Label.String = 'Z-score';
    
            ax = subplot(4,3,12);
            ax.FontSize = 12;
            hold on
            for iF = 1:length(fbands)
                scatter(mean(fbands{iF}), out.fr.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                    'MarkerEdgeColor', c_list{iF})
            end
            ax.XLim = [0 60];
            ax.YLim = [-3 12];
            ax.XAxis.Label.String = 'Freq (Hz)';
            ax.YAxis.Label.String = 'Z-score';
            
            % Save variables and figure
            save(strcat(fn_prefix, '_nonstim_phase_response_',string(bin_counts(iN)),'_bins'), 'ns_out');
            print(this_fig, '-dpng', strcat(fn_prefix, '_nonstim_phase_response_hist_', string(bin_counts(iN)),'_bins'));
%             savefig(this_fig, strcat(fn_prefix, '_phase_response_hist_', string(nbins),'_bins'));
            close;
        end
    end
end