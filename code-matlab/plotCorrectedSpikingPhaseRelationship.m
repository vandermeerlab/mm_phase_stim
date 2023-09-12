%% Script to generate the corrected spike phase relationship
% Assumes that *stim_phase_5_bins.mat and *nonstim_phase_5_bins.mat exits
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
    phase_bins = [-pi:2*pi/5:pi];
    x_ticks = 0.5*(phase_bins(1:end-1)+phase_bins(2:end));
    
    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        % Load the stim phase responses
        load(strcat(fn_prefix, '_phase_response_5_bins.mat'));

        % Load the nonstim phase responses
        load(strcat(fn_prefix, '_nonstim_phase_response_5_bins.mat'));
        corrected_fr_bin = out.fr.bin - ns_out.fr.bin;
        [corrected_fr_r, corrected_fr_z] = deal(zeros(length(fbands),1));

        fig = figure('WindowState', 'maximized');
        for iF = 1:length(fbands)
            this_fr = corrected_fr_bin(iF,:);
            corrected_fr_r(iF) = (max(this_fr) - min(this_fr))/(max(this_fr) + min(this_fr));
            q0 = squeeze(out.fr.shufs(iF,:,:)) - ns_out.fr.bin(iF,:);
            q1 = max(q0, [] , 2);
            q2 = min(q0, [], 2);
            q3 = (q1 - q2)./(q1 + q2); 
            corrected_fr_z(iF) = (corrected_fr_r(iF) - mean(q3))/std(q3);

            % Load stim Response for side by side plotting
            % Load PLV and bin those_phases
            ax = subplot(4,4,(iF-1)*4+1);
            fn2 = fn_prefix;
            fn2(end-1) = '-';
            load(strcat(fn2, '_spike_phaselock_causal_plv.mat'));
            if ~isempty(trial_spk_phase)  
                this_phase = trial_spk_phase(iF,:);
                [this_count, ~, ~] = histcounts(this_phase,phase_bins);
                bar(ax,x_ticks,this_count/sum(this_count),1,c_list{iF});
                ax.Title.String = sprintf('PLV phase distribution: %d', length(this_phase));
                ax.Title.FontSize = 20;
                ax.XLim = [-3.25 3.25];
            end

            ax = subplot(4,4,(iF-1)*4+2);
            bar(ax,x_ticks,ns_out.fr.bin(iF,:),1,c_list{iF});
            ax.Title.String = '{\Delta} FR nonStim';
            ax.Title.FontSize = 20;
            ax.XLim = [-3.25 3.25];
%                 ax.YLim = [-10 10];

            ax = subplot(4,4,(iF-1)*4+3);
            bar(ax,x_ticks,out.fr.bin(iF,:),1,c_list{iF});
            ax.Title.String = '{\Delta} FR OptoStim';
            ax.Title.FontSize = 20;
            ax.XLim = [-3.25 3.25];
%           ax.YLim = [0 100];

            ax = subplot(4,4,(iF-1)*4+4);
            bar(ax,x_ticks,corrected_fr_bin(iF,:),1,c_list{iF});
            ax.Title.String = '{\Delta} FR corrected';
            ax.Title.FontSize = 20;
            ax.XLim = [-3.25 3.25];
%           ax.YLim = [0 100];
        end
                     
        ax = subplot(4,4,14);
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

        ax = subplot(4,4,15);
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

        ax = subplot(4,4,16);
        ax.FontSize = 12;
        hold on
        for iF = 1:length(fbands)
            scatter(mean(fbands{iF}), corrected_fr_z(iF), 'MarkerFaceColor', c_list{iF}, ...
                'MarkerEdgeColor', c_list{iF})
        end
        ax.XLim = [0 60];
        ax.YLim = [-3 12];
        ax.XAxis.Label.String = 'Freq (Hz)';
        ax.YAxis.Label.String = 'Z-score';
        print(fig, '-dpng', strcat(fn_prefix, '_corrected_phase_response_hist_5_bins'));
        close;
    end
end