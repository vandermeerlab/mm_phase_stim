%% Assumes that good LFPs have been picked out

top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.depth, summary.freq, summary.og, ...
    summary.fooof, summary.irasa] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end
dStr_mask = summary.depth < 3.5;
% fband = {[2 5], [6 10], [12 30], [30, 55]};
fband = {[2 5], [6 10], [30, 55]};
%% Normalize all the psds
for i = 1:length(summary.depth)
    norm_og(i,:) = (summary.og(i,:) - min(summary.og(i,:)))/(max(summary.og(i,:)) - min(summary.og(i,:)));
    norm_fooof(i,:) = (summary.fooof(i,:) - min(summary.fooof(i,:)))/(max(summary.fooof(i,:)) - min(summary.fooof(i,:)));
    norm_irasa(i,:) = (summary.irasa(i,:) - min(summary.irasa(i,:)))/(max(summary.irasa(i,:)) - min(summary.irasa(i,:)));
end


%%


% Plot dStr
sel =  find(dStr_mask);
subplot(2,3,1)
imagesc(norm_og(dStr_mask,:))
hold on
for iF = 1:length(fband)
    F1 = fband{iF}(1);
    F2 = fband{iF}(2);
    rectangle('Position',[F1,0,F2-F1,length(sel)+1],'EdgeColor', 'red', 'LineWidth', 2)
end
title('dStr OG');

subplot(2,3,2)
imagesc(norm_fooof(dStr_mask,:))
hold on
for iF = 1:length(fband)
    F1 = fband{iF}(1);
    F2 = fband{iF}(2);
    rectangle('Position',[F1,0,F2-F1,length(sel)+1],'EdgeColor', 'red', 'LineWidth', 2)
end
title('dStr FOOOF');

subplot(2,3,3)
imagesc(norm_irasa(dStr_mask,:))
hold on
for iF = 1:length(fband)
    F1 = fband{iF}(1);
    F2 = fband{iF}(2);
    rectangle('Position',[F1,0,F2-F1,length(sel)+1],'EdgeColor', 'red', 'LineWidth', 2)
end
title('dStr IRASA')

% Plot vStr stuff
sel =  find(~dStr_mask);
subplot(2,3,4)
imagesc(norm_og(~dStr_mask,:))
hold on
for iF = 1:length(fband)
    F1 = fband{iF}(1);
    F2 = fband{iF}(2);
    rectangle('Position',[F1,0,F2-F1,length(sel)+1],'EdgeColor', 'red', 'LineWidth', 2)
end
title('vStr OG');

subplot(2,3,5)
imagesc(norm_fooof(~dStr_mask,:))
hold on
for iF = 1:length(fband)
    F1 = fband{iF}(1);
    F2 = fband{iF}(2);
    rectangle('Position',[F1,0,F2-F1,length(sel)+1],'EdgeColor', 'red', 'LineWidth', 2)
end
title('vStr FOOOF');

subplot(2,3,6)
imagesc(norm_irasa(~dStr_mask,:))
hold on
for iF = 1:length(fband)
    F1 = fband{iF}(1);
    F2 = fband{iF}(2);
    rectangle('Position',[F1,0,F2-F1,length(sel)+1],'EdgeColor', 'red', 'LineWidth', 2)
end
title('vStr IRASA')

%% Plot IRASA  
% grey out frequencies around 60 Hz 
norm_irasa(:,59:61) = nan;
fig = figure('WindowState', 'maximized');

sel =  find(dStr_mask);
% Plot dStr stuff
subplot(2,1,1)
imagesc(norm_irasa(dStr_mask,:))
hold on
for iF = 1:length(fband)
    F1 = fband{iF}(1);
    F2 = fband{iF}(2);
    rectangle('Position',[F1,0,F2-F1,length(sel)+1],'EdgeColor', 'red', 'LineWidth', 2)
end
title('dStr IRASA');
xticks([2 5 6 10 12 30 55, 60, 100])
xlim([1 100])
xlabel('Frequency (Hz)')
ylabel('Sessions')
yticks([])
ax = gca;
box off;
ax.TickDir = 'out';


% Plot vStr stuff
sel =  find(~dStr_mask);
subplot(2,1,2)
imagesc(norm_irasa(~dStr_mask,:))
hold on
for iF = 1:length(fband)
    F1 = fband{iF}(1);
    F2 = fband{iF}(2);
    rectangle('Position',[F1,0,F2-F1,length(sel)+1],'EdgeColor', 'red', 'LineWidth', 2)
end
title('vStr IRASA');
xticks([2 5 6 10 12 30 55, 60, 100])
xlim([1 100])
ylabel('Sessions')
xlabel('Frequency (Hz)')
yticks([])
ax = gca;
box off;
ax.TickDir = 'out';
fontname(fig, 'Helvetica');
fig.Renderer = 'painters';
fontsize(fig, 30, 'points');
%%
function s_out = doStuff(s_in)
    s_out = s_in;
    
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end

    % Go ahead only if this exists in a list of opto_cells
    
    % Load the list of final opto cells
    load('E:\Dropbox (Dartmouth College)\AnalysisResults\phase_stim_results\FinalOptoCells.mat');
    if (sum(contains(ExpKeys.goodCell, dStr_opto)) == 0) & (sum(contains(ExpKeys.goodCell, vStr_opto)) == 0)
        return
    end

    load('trial_psd.mat');
    s_out.freq  = [s_out.freq; psd.freq];
    s_out.og    = [s_out.og; 10*log10(psd.original)];
    s_out.fooof = [s_out.fooof; 10*log10(psd.original) - 10*log10(psd.fooof)];
    s_out.irasa = [s_out.irasa; 10*log10(psd.original) - 10*log10(psd.irasa)];
    s_out.depth = [s_out.depth; ExpKeys.probeDepth]; 
end