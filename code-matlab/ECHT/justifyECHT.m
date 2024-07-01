%% Assumes that good LFPs have been picked out 'echt_baseline_evaluation.mat' exists

top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.mean, summary.std] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end

fbands = {[2 5], [6 10], [30 55]};
c_list = {'red', 'blue', 'green'};
c_list2 = {'black','red', 'blue', 'green', 'cyan'};
%% 
iF = 2;
for i = 1:53
    for j = 1:5
        scatter((i-1)*10+j,summary.mean(i,iF,j), c_list2{j})
        errorbar((i-1)*10+j,summary.mean(i,iF,j), summary.std(i,iF,j), 'Color', c_list2{j}, 'LineStyle', 'none')
        hold on
    end
end

%% For each frequency, what's the mean across all neurons for each time window (mean of means) and standard deviation

fig = figure('WindowState', 'maximized');
[all_mean, all_std] = deal(nan(length(fbands), 5));
for iF = 1:length(fbands)
    temp = squeeze(summary.mean(:,iF,:));
    all_mean(iF,:) = circ_mean(temp);
    all_std(iF,:) = circ_std(temp);

    ax = subplot(1,3,iF);
    hold on
%     scatter([0.5,0.75,1,1.25,1.5], all_mean(iF,:), c_list{iF});
    errorbar([0.5,0.75,1,1.25,1.5], all_mean(iF,:), all_std(iF,:), 'Color', ...
    c_list{iF});
% scatter([0.5,0.75,1,1.25,1.5], all_std(iF,:), c_list{iF});
    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
    xlabel('Window length (s)')
    ylabel('True Phase - causal phase')
    ylim([-pi-0.1 pi+0.1])
    xlim([0.4 1.6])
    yline(0, '--black');
    for iB = 1:6
        yline(-pi+(iB-1)*2*pi/5, 'black', 'LineWidth', 2)
    end
    xticks([0.5 0.75 1 1.25 1.5])
    ax.TickDir = 'out';
    ax.Box = 'off';
end

%%
function s_out = doStuff(s_in)
    s_out = s_in;
    idx = size(s_out.mean,1)+1; 
    LoadExpKeys;
    evs = LoadEvents([]);
    cfg_spk = [];
    cfg_spk.fc = ExpKeys.goodCell;
    if isempty(cfg_spk.fc)
        return
    end
    
    load('echt_baseline_evaluation.mat');
    for iF  = 1:3
        [temp_mean, temp_std] = deal([]);
        for iW = 1:5
            q0 = circ_stats(eval_acausal_phase(iF,:) - squeeze(eval_causal_phase(iW,iF,:))');
            temp_mean = [temp_mean, q0.mean];
            temp_std = [temp_std, q0.std];
        end
        s_out.mean(idx,iF,:) = temp_mean;
        s_out.std(idx,iF,:) = temp_std;
    end
end

