%% Script to collect and summarize PETHS for all isolated units
% Assumes that all_peths.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\'
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', ...
    'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', ...
    'M319', 'M321', 'M325'};
summary = [];
[summary.label, summary.peth, summary.zpeth, summary.shuf_mean, ...
    summary.shuf_zmean, summary.shuf_sd, summary.shuf_zsd, ...
    summary.tvec, summary.is_opto] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end
%% adding the 'is_really_opto' field
opto_labels = load('output\FinalOptoCells.mat');
temp0 = cellfun(@(x) any(strcmp(x, opto_labels.vStr_opto)), summary.label);
temp1 = cellfun(@(x) any(strcmp(x, opto_labels.dStr_opto)), summary.label);
summary.is_really_opto = temp0|temp1;

%% Plot mean peth for all z-scored putative opto cells
opto_keep = summary.is_opto == 1;
opto_peth = summary.zpeth(opto_keep,:);
opto_shuf_mean = summary.shuf_zmean(opto_keep,:);
opto_shuf_sd = summary.shuf_zsd(opto_keep,:);
non_opto_peth = summary.zpeth(~opto_keep,:);
non_opto_shuf_mean = summary.shuf_zmean(~opto_keep,:);
non_opto_shuf_sd = summary.shuf_zsd(~opto_keep,:);
%
figure;
subplot(1,2,1)
hold on
plot(summary.tvec, mean(opto_peth, 1), 'Color', 'blue', 'LineWidth', 2);
plot(summary.tvec, mean(opto_shuf_mean + 2*opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
plot(summary.tvec, mean(opto_shuf_mean - 2*opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
fill([summary.tvec, fliplr(summary.tvec)], [mean(opto_shuf_mean + 2*opto_shuf_sd, 1), ...
    fliplr(mean(opto_shuf_mean - 2*opto_shuf_sd, 1))], 'black', 'FaceAlpha', 0.1);

xline(0, '--black')
xlabel('Time from stim onset (s)')
ylabel('Z-scored firing-rate')
title('FSIs')
subplot(1,2,2)
hold on
plot(summary.tvec, mean(non_opto_peth, 1), 'Color', 'red', 'LineWidth', 2);
plot(summary.tvec, mean(non_opto_shuf_mean + 2*non_opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
plot(summary.tvec, mean(non_opto_shuf_mean - 2*non_opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);

% Alsp plot the clearly inhibited opto cells
plot(summary.tvec, summary.zpeth(2,:), 'Color', 'cyan', 'LineWidth',1);
plot(summary.tvec, summary.shuf_zmean(2,:) + 2*summary.shuf_zsd(2,:) , '--cyan')
plot(summary.tvec, summary.shuf_zmean(2,:) - 2*summary.shuf_zsd(2,:) , '--cyan')
% plot(summary.tvec, summary.zpeth(80,:), 'Color', 'green', 'LineWidth',1);
% plot(summary.tvec, summary.shuf_zmean(80,:) + 2*summary.shuf_zsd(80,:) , '--cyan')
% plot(summary.tvec, summary.shuf_zmean(80,:) - 2*summary.shuf_zsd(80,:) , '--cyan')
plot(summary.tvec, summary.zpeth(4,:), 'Color', 'magenta', 'LineWidth',1);
plot(summary.tvec, summary.shuf_zmean(4,:) + 2*summary.shuf_zsd(4,:) , '--magenta')
plot(summary.tvec, summary.shuf_zmean(4,:) - 2*summary.shuf_zsd(4,:) , '--magenta')

fill([summary.tvec, fliplr(summary.tvec)], [mean(non_opto_shuf_mean + 2*non_opto_shuf_sd, 1), ...
    fliplr(mean(non_opto_shuf_mean - 2*non_opto_shuf_sd, 1))], 'black', 'FaceAlpha', 0.1);
xline(0, '--black')
xlabel('Time from stim onset (s)')
ylabel('Z-scored firing-rate')
title('MSNs')
%% Plot mean peth after normalizing for all z-scored putative opto cells
opto_keep = summary.is_opto == 1;
opto_peth = summary.zpeth(opto_keep,:);
opto_shuf_mean = summary.shuf_zmean(opto_keep,:);
opto_shuf_sd = summary.shuf_zsd(opto_keep,:);

% Normalize all the peths
for i = 1:length(summary.depth)
    norm_og(i,:) = (summary.og(i,:) - min(summary.og(i,:)))/(max(summary.og(i,:)) - min(summary.og(i,:)));
    norm_fooof(i,:) = (summary.fooof(i,:) - min(summary.fooof(i,:)))/(max(summary.fooof(i,:)) - min(summary.fooof(i,:)));
    norm_irasa(i,:) = (summary.irasa(i,:) - min(summary.irasa(i,:)))/(max(summary.irasa(i,:)) - min(summary.irasa(i,:)));
end


non_opto_peth = summary.zpeth(~opto_keep,:);
non_opto_shuf_mean = summary.shuf_zmean(~opto_keep,:);
non_opto_shuf_sd = summary.shuf_zsd(~opto_keep,:);
%
figure;
subplot(1,2,1)
hold on
plot(summary.tvec, mean(opto_peth, 1), 'Color', 'blue', 'LineWidth', 2);
plot(summary.tvec, mean(opto_shuf_mean + 2*opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
plot(summary.tvec, mean(opto_shuf_mean - 2*opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
xline(0, '--black')
xlabel('Time from stim onset (s)')
ylabel('Z-scored firing-rate')
title('Opto cells')
subplot(1,2,2)
hold on
plot(summary.tvec, mean(non_opto_peth, 1), 'Color', 'red', 'LineWidth', 2);
plot(summary.tvec, mean(non_opto_shuf_mean + 2*non_opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
plot(summary.tvec, mean(non_opto_shuf_mean - 2*non_opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
xline(0, '--black')
xlabel('Time from stim onset (s)')
ylabel('Z-scored firing-rate')
title('Non Opto cells')

%%
function s_out = doStuff(s_in)
    s_out = s_in;
    % If all_peth.mat doesn't exist, skip
    if ~isfile('all_peths.mat')
        return
    end
    load("all_peths.mat");
    s_out.label = [s_out.label; string(out.labels)];
    s_out.peth = [s_out.peth; out.stim_peth];
    s_out.zpeth = [s_out.zpeth; out.stim_zpeth];
    s_out.shuf_mean = [s_out.shuf_mean; out.shuf_peth_mean];
    s_out.shuf_sd = [s_out.shuf_sd; out.shuf_peth_sd];
    s_out.shuf_zmean = [s_out.shuf_zmean; out.shuf_zpeth_mean];
    s_out.shuf_zsd = [s_out.shuf_zsd; out.shuf_zpeth_sd];
    s_out.is_opto = [s_out.is_opto; out.is_opto];
    s_out.tvec = out.tvec;
end