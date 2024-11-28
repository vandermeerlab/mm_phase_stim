%% Script to plot firing rate modulation by phase plots for opto cells and select co-recorded MSNs
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        this_label = {'M295-2022-01-06-TT06_5.t', 'M295-2022-01-06-TT06_4.t', ...
            'M295-2022-01-06-TT08_4.t',}; % n.s dStr FSI with clearly inhibited co-recorded cells
%         this_label = {'M295-2022-01-06-TT06_5.t', 'M295-2022-01-06-TT06_4.t', ...
%             'M295-2022-01-06-TT06_6.t','M295-2022-01-06-TT07_1.t', ...
%             'M295-2022-01-06-TT07_3.t', 'M295-2022-01-06-TT08_4.t',}; % n.s dStr FSI with clearly inhibited co-recorded cells
%         this_label = {'M320-2022-05-28-TT01_1.t', 'M320-2022-05-28-TT02_2.t', ...
%             'M320-2022-05-28-TT03_1.t','M320-2022-05-28-TT04_4.t', ...
%             'M320-2022-05-28-TT05_3.t'}; % sig vStr FSI with some inhibited co-recorded MSNS
        doStuff(this_label)
    end

end
%%
function doStuff(label)
    LoadExpKeys;
    % if not the required cell, skip
    if isempty(ExpKeys.goodCell) | ~strcmp(ExpKeys.goodCell, label{1})
        return
    end
    
    evs = LoadEvents([]);
    cfg = [];
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.min_cluster_quality = 3;
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
  
    S = LoadSpikes(cfg);
    % Select the opto cell and other 'good' co-recorded cell

    S_opto = SelectTS([],S,ismember(S.label, label(1)));
    S_others = SelectTS([],S,ismember(S.label, label(2:end)));

    % Set variables
    nbins = 5;
    fbands = {[2 5], [6 10], [12 28], [30 55]};
    c_list = {'red', 'blue','magenta', 'cyan'}; % One color for each frequency band
    % Alternate color scheme: 1 color for each bin
    temp = linspecer(nbins);
    c_rgb = {};
    for iB = 1:nbins
        c_rgb{iB} = temp(iB,:);
    end
%     c_rgb = {[1 0 0], [0 1 1], [1 0 1], [0 1 0], [0 0 1]};
    clear temp

    fig = figure('WindowState', 'maximized');
    nC = length(label);


    % Plot opto FSI first
    this_prefix = extractBefore(label{1}, '.t');
    load(strcat(this_prefix, '_phase_response_5_bins.mat'));
    
    % First plot the individual FR Mods
    for iF = 1:length(fbands)
        ax = subplot(length(fbands), nC,nC*(iF-1)+1);
        hold on
        b = bar(out.fr.bin(iF,:),1, 'FaceColor', c_list{iF});
        q0 = squeeze(out.fr.shufs(iF,:,:));
        q1 = max(q0,[],2);
        q2 =  min(q0,[],2);
        q3 = (q1 - q2)./(q2 + q1);
        qz2 = 2*std(q3)+ mean(q3);
        [qmin, min_idx] = min(out.fr.bin(iF,:));
        [qmax, max_idx] = max(out.fr.bin(iF,:));
        qoff = qmin*(qz2 + 1)/(1 - qz2);
        plot([min_idx, 5.65], [qmin, qmin], '--black')   
        plot([max_idx, 5.65], [qmax, qmax], '--black')
        errorbar(5.75, 0.5*(qmax+qmin), 0.5*(qmax-qmin), 'black', 'LineStyle', 'none');
        errorbar(6, 0.5*(qoff+qmin), 0.5*(qoff-qmin), 'red', 'LineStyle', 'none');
        xlim([0.5 6.5])
        ylim([0 50]) % Might need to be changed on a case to case basis
        xlabel('Phase bin')
        ylabel('{\Delta} FR (Hz)')
        ax.Title.String = sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2));
        ax.YAxis.FontSize = 12;
        ax.XAxis.FontSize = 12;
        ax.TickDir = 'out';
    end
    
    % Plot the scatter of z-scores
    ax = subplot(length(fbands)+1, nC, nC*length(fbands)+1);
    hold on;
    for iF = 1:length(fbands)
        scatter(mean(fbands{iF}), out.fr_ws.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                    'MarkerEdgeColor', c_list{iF})
    end
    ylim([-3 8])
    yticks([-3 0 2 8])
    yline(2, '--black')
    xlim([0 60])
    xticks(cellfun(@(x) mean(x), fbands))
    xlabel('Frequency (Hz)')
    ylabel('Z-score')
    ax.YAxis.FontSize = 12;
    ax.XAxis.FontSize = 12;
    ax.TickDir = 'out';
%     title ("")
    dummy = 1;

    clear out
    % Next plot the co-recorded cells
    for iC = 2:length(label)
        this_prefix = extractBefore(label{iC}, '.t');
        load(strcat(this_prefix, '_phase_response_5_bins.mat'));

        % First plot the individual FR Mods
        for iF = 1:length(fbands)
            ax = subplot(length(fbands)+1, nC,nC*(iF-1)+iC);
            hold on
            b = bar(out.fr.bin(iF,:),1, 'FaceColor', c_list{iF});
            q0 = squeeze(out.fr.shufs(iF,:,:));
            q1 = max(q0,[],2);
            q2 =  min(q0,[],2);
            q3 = abs((q1 - q2)./(q2 + q1));
            qz2 = 2*std(q3)+ mean(q3);
            [qmin, min_idx] = min(out.fr.bin(iF,:));
            [qmax, max_idx] = max(out.fr.bin(iF,:));
            % Think and fix later
            qoff = qmax*(qz2 + 1)/(1-qz2);
            plot([min_idx, 5.65], [qmin, qmin], '--black')   
            plot([max_idx, 5.65], [qmax, qmax], '--black')
            errorbar(5.75, 0.5*(qmax+qmin), 0.5*(qmax-qmin), 'black', 'LineStyle', 'none');
            errorbar(6, 0.5*(qoff+qmax), 0.5*(qoff-qmax), 'red', 'LineStyle', 'none');
            xlim([0.5 7.5])
            ylim([-20 2]) % Might need to be changed on a case to case basis
            xlabel('Phase bin')
            ylabel('{\Delta} FR (Hz)')
            ax.Title.String = sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2));
            ax.YAxis.FontSize = 12;
            ax.XAxis.FontSize = 12;
            ax.TickDir = 'out';
        end
        
        % Plot the scatter of z-scores
        ax = subplot(length(fbands)+1, nC, nC*length(fbands)+iC);
        hold on;
        for iF = 1:length(fbands)
            scatter(mean(fbands{iF}), out.fr_ws.zscore(iF), 'MarkerFaceColor', c_list{iF}, ...
                        'MarkerEdgeColor', c_list{iF})
        end
        ylim([-3 8])
        yticks([-3 0 2 8])
        yline(2, '--black')
        xlim([0 60])
        xticks(cellfun(@(x) mean(x), fbands))
        xlabel('Frequency (Hz)')
        ylabel('Z-score')
        ax.YAxis.FontSize = 12;
        ax.XAxis.FontSize = 12;
        ax.TickDir = 'out';
    %     title ("")
        clear out
    end

end

