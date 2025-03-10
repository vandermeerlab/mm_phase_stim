%% Script to characterize opto_cells on the basis of firing rate changes and mean wave-form
% Assumes that stim_phases.mat and *stim_response.mat already exist in each folder
top_dir = 'data\'
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.nostim_delta, summary.opto_delta, ...
    summary.sham_delta, summary.poisson_delta, summary.nostim_delta_sd, summary.opto_delta_sd, ...
    summary.sham_delta_sd, summary.poisson_delta_sd, ...
    summary.p_val, summary.h, summary.p_val_t, summary.h_t, ...
    summary.isopto, summary.depth, summary.waveforms, summary.mfr,  ...
    summary.depth, summary.peth, summary.zpeth, summary.shuf_peth, ...
    summary.shuf_zpeth, summary.norm_shuf_zpeth, summary.is_opto, ...
    summary.tvec] = deal([]);

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
vStr_mask = summary.depth >= 3.5;
ks_mask = summary.h == 1; % That the sham and opto difference in firing-rate is different
og_opto_mask = summary.isopto == 1;
pct_mask = summary.h_t == 1; 

peak_to_trough = zeros(size(summary.depth));
for i = 1:length(summary.waveforms)
    norm_wf(i,:) = (summary.waveforms(i,:) - min(summary.waveforms(i,:)))/...
        (max(summary.waveforms(i,:)) - min(summary.waveforms(i,:)));
    % Reversing the waveform since input was inverted
    norm_wf(i,:) = -1 * norm_wf(i,:);
    [~, pidx] = max(norm_wf(i,:));
    [~, tidx] = min(norm_wf(i,:));
    peak_to_trough(i) = tidx - pidx;
end
%% Figure for showing Sham deltaFR vs Opto deltaFR
% Plot vStr stuff
fig = figure('WindowState','maximized');

% Plot dStr stuff
ax = subplot(2,4,[1 5]);
hold on
sel = find(ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'blue');
end

sel = find(~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'red');
end
xticks([1 2])
yticks([0 160 320])
xticklabels({'Sham stim', 'Opto stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-25 320])
title('All striatal neurons');
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';
ax.XAxis.FontSize = 30;
ax.XLabel.FontSize = 25;
ax.YAxis.FontSize = 30;
ax.YLabel.FontSize = 25;
xtickangle(0)

ax = subplot(2,4,2);
hold on
sel1 = find(ks_mask & summary.sham_delta < summary.opto_delta);
sel2 = find(ks_mask & summary.sham_delta >= summary.opto_delta);
for i = 1:length(sel1)
    h = plot((0:1:31)/32, norm_wf(sel1(i),:), 'blue');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
for i = 1:length(sel2)
    h = plot((0:1:31)/32, norm_wf(sel2(i),:), '--blue');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot((0:1:31)/32,mean(norm_wf(sel1,:)), 'black', 'LineWidth', 3);
xlabel('Time (ms)')
ylabel('Normalized Amplitude')
% legend({sprintf('Peak to trough %.2f +/- %.2f ms', ...
%     mean(peak_to_trough(sel1))/32, std(peak_to_trough(sel1)/32))}, 'FontSize', 12)
title('opto-responsive')
xticks([0 0.5 1]);
yticks([0 0.5 1]);
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';
ax.XAxis.FontSize = 30;
ax.XLabel.FontSize = 25;
ax.YAxis.FontSize = 30;
ax.YLabel.FontSize = 25;

ax = subplot(2,4,6);
hold on
sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    h = plot((0:1:31)/32, norm_wf(sel(i),:), 'red');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot((0:1:31)/32,mean(norm_wf(sel,:)), 'black', 'LineWidth', 3);
xlabel('Time (ms)')
ylabel('Normalized Amplitude')
% legend({sprintf('Peak to trough %.2f +/- %.2f ms', ...
%     mean(peak_to_trough(sel))/32, std(peak_to_trough(sel)/32))}, 'FontSize', 12)
title('non-responsive')
xticks([0 0.5 1]);
yticks([0 0.5 1]);
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';
ax.XAxis.FontSize = 30;
ax.XLabel.FontSize = 25;
ax.YAxis.FontSize = 30;
ax.YLabel.FontSize = 25;

% Plot the population peths
opto_keep = summary.is_opto == 1;
opto_peth = summary.zpeth(opto_keep,:);

% Norm opto_peth row_wise
row_min = min(opto_peth, [], 2);
row_max = max(opto_peth,[],2);
% Create matrices for broadcasting
min_matrix = repmat(row_min, 1, size(opto_peth, 2));
max_matrix = repmat(row_max, 1, size(opto_peth, 2));
% Normalize
opto_peth = (opto_peth - min_matrix) ./ (max_matrix - min_matrix);

% First smooth and then 

opto_shuf_dist = mean(summary.norm_shuf_zpeth(:,:,opto_keep),3);
opto_shuf_mean = mean(opto_shuf_dist,1);
opto_shuf_sd = std(opto_shuf_dist,1);

non_opto_peth = summary.zpeth(~opto_keep,:);
non_opto_shuf_dist = nanmean(summary.norm_shuf_zpeth(:,:,~opto_keep),3);
non_opto_shuf_mean = mean(non_opto_shuf_dist,1);
non_opto_shuf_sd = std(non_opto_shuf_dist,1);

% Parameters for gaussian kernel
sigma = 0.05;  % Standard deviation in seconds
timebase = 0.01;  % Your timebase in seconds
window = 0.2;  % Width of kernel in seconds (e.g., 4*sigma)

% Create gaussian kernel
kernel_t = -window:timebase:window;  % Kernel time vector
kernel = exp(-(kernel_t.^2)/(2*sigma^2));
kernel = kernel/sum(kernel);  % Normalize to sum to 1

ax = subplot(2,4,[3,4]);
hold on
plot(summary.tvec, mean(opto_peth, 1), 'Color', 'cyan', 'LineWidth', 2);
plot(summary.tvec, mean(opto_shuf_mean + 2*opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
plot(summary.tvec, mean(opto_shuf_mean - 2*opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
fill([summary.tvec, fliplr(summary.tvec)], [mean(opto_shuf_mean + 2*opto_shuf_sd, 1), ...
    fliplr(mean(opto_shuf_mean - 2*opto_shuf_sd, 1))], 'black', 'FaceAlpha', 0.1);

xline(0, '--black')
xlabel('Time from stim onset (s)')
ylabel('Z-scored firing rate')
title('FSIs')
ylim([-2 8])
yticks([-2 0 2 4 6 8])
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';
ax.XAxis.FontSize = 30;
ax.XLabel.FontSize = 25;
ax.YAxis.FontSize = 30;
ax.YLabel.FontSize = 25;

ax = subplot(2,4,[7,8]);
hold on
% plot(summary.tvec, conv(mean(non_opto_peth, 1),kernel,'same'), 'Color', 'magenta', 'LineWidth', 2);
plot(summary.tvec, mean(non_opto_shuf_mean + 2*non_opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);
plot(summary.tvec, mean(non_opto_shuf_mean - 2*non_opto_shuf_sd, 1), 'Color', 'black', 'LineWidth', 1);

fill([summary.tvec, fliplr(summary.tvec)], [mean(non_opto_shuf_mean + 2*non_opto_shuf_sd, 1), ...
    fliplr(mean(non_opto_shuf_mean - 2*non_opto_shuf_sd, 1))], 'black', 'FaceAlpha', 0.1);
xline(0, '--black')
xlabel('Time from stim onset (s)')
ylabel('Z-scored firing rate')
title('MSNs')
ylim([-0.02 0.02])
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';
ax.XAxis.FontSize = 30;
ax.XLabel.FontSize = 25;
ax.YAxis.FontSize = 30;
ax.YLabel.FontSize = 25;

% fontsize(fig, 30, 'points');
fontname(fig, 'Helvetica');
fig.Renderer = 'painters';

exportgraphics(fig, 'output\fig.eps', 'BackgroundColor', 'none', 'ContentType', 'vector')

%% Use this to write list of Final Opto Cells
dStr_opto_fsi = summary.labels(ks_mask & dStr_mask & summary.opto_delta > summary.sham_delta);
vStr_opto_fsi = summary.labels(ks_mask & vStr_mask & summary.opto_delta > summary.sham_delta);

dStr_opto_msn = summary.labels(ks_mask & dStr_mask & summary.opto_delta < summary.sham_delta);
vStr_opto_msn = summary.labels(ks_mask & vStr_mask & summary.opto_delta < summary.sham_delta);
%%
dStr_others = summary.labels(dStr_mask & ~(ks_mask & summary.opto_delta > summary.sham_delta));
vStr_others = summary.labels(vStr_mask & ~(ks_mask & summary.opto_delta > summary.sham_delta));

%%
opto_p = [summary.labels(ks_mask),summary.h(ks_mask),summary.p_val(ks_mask), summary.sham_delta(ks_mask), summary.opto_delta(ks_mask)];
other_p = [summary.labels(~ks_mask),summary.h(~ks_mask), summary.p_val(~ks_mask), summary.sham_delta(~ks_mask), summary.opto_delta(~ks_mask)];
%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
%     fprintf('%s\n',pwd);
    cfg = [];
%     cfg.fc = ExpKeys.goodCell;
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.min_cluster_quality = 3;
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
    S = LoadSpikes(cfg);

    for iC = 1:length(S.t)
        s_out.isopto = [s_out.isopto; contains(S.label{iC},ExpKeys.goodCell)];
        fn_prefix = extractBefore(S.label{iC}, '.t');
        
        % Load 20 ms if not goodcell, else load normal
        if contains(S.label{iC},ExpKeys.goodCell)
            load(strcat(fn_prefix,'_stim_response.mat'));
        else
            load(strcat(fn_prefix,'_stim_response20.mat'));
        end

        s_out.mfr = [s_out.mfr; od.trial_stim.mfr];
%         nostim_dfr = od.trial_nonstim.fr' - od.trial_nonstim.bfr';
        % Comapare to a similar distribution as before
        ridx = randi(1000000,1,10000);
        sham_dfr = od.sham_stim.fr(ridx)' - od.sham_stim.bfr(ridx)';
%         poisson_dfr = od.sham_stim.p_fr' - od.sham_stim.p_bfr';
        opto_idx = find(strcmp(S.label{iC},ExpKeys.goodCell));

        if isempty(opto_idx) % Non opto_cell
            opto_dfr = od.trial_stim.fr(min(min(ExpKeys.goodTrials)):max(max(ExpKeys.goodTrials))) - ...
                od.trial_stim.bfr(min(min(ExpKeys.goodTrials)):max(max(ExpKeys.goodTrials)));
        else
            opto_dfr = od.trial_stim.fr(ExpKeys.goodTrials(opto_idx,1):ExpKeys.goodTrials(opto_idx,2)) - ...
                od.trial_stim.bfr(ExpKeys.goodTrials(opto_idx,1):ExpKeys.goodTrials(opto_idx,2));
        end
        [h,p,~] = kstest2(sham_dfr, opto_dfr, 'Alpha',0.01);
        s_out.h = [s_out.h; h];
        s_out.p_val = [s_out.p_val; p];
        
        % Do a percentile-test based shuffle
        nshuf = 1000;
        sham_dfr_mean = zeros(1,nshuf);
        sham_dfr_sd = zeros(1,nshuf);
        for ishuf = 1:nshuf
            ridx = randi(1000000,1,(max(max(ExpKeys.goodTrials)) - ...
                min(min(ExpKeys.goodTrials)) + 1));
            this_shuf_dfr = od.sham_stim.fr(ridx) - od.sham_stim.bfr(ridx);
            sham_dfr_mean(ishuf) = mean(this_shuf_dfr);
            sham_dfr_sd(ishuf) = std(this_shuf_dfr);
        end        
        h2= (mean(opto_dfr) > prctile(sham_dfr_mean, 99)) | ...
            (mean(opto_dfr) < prctile(sham_dfr_mean, 1));
        p2 = max(sum(mean(opto_dfr) > sham_dfr_mean), ...
            sum(mean(opto_dfr) < sham_dfr_mean))/1000;
        s_out.h_t = [s_out.h_t; h2];
        s_out.p_val_t = [s_out.p_val_t; 1- p2];

        
%         % Plot x
        s_out.opto_delta = [s_out.opto_delta; mean(opto_dfr)];
        s_out.opto_delta_sd = [s_out.opto_delta_sd; std(opto_dfr)];
%         s_out.nostim_delta = [s_out.nostim_delta; mean(nostim_dfr)];
%         s_out.nostim_delta_sd = [s_out.nostim_delta_sd; std(nostim_dfr)];
        
        s_out.sham_delta = [s_out.sham_delta; mean(sham_dfr)];
        s_out.sham_delta_sd = [s_out.sham_delta_sd; std(sham_dfr)];
%         s_out.poisson_delta = [s_out.poisson_delta, mean(poisson_dfr)];
%         s_out.poisson_delta_sd = [s_out.poisson_delta, std(poisson_dfr)];

        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];

        fn_prefix = strrep(fn_prefix, '_', '-');
        % Load the mean waveform
        load(strcat(fn_prefix,'-wv.mat'));
        [~,midx]  = max(max(mWV)); % Getting the max amplitude channel waveform
        s_out.waveforms = [s_out.waveforms; mWV(:,midx)'];
    end
    % If all_peth.mat doesn't exist, skip
    if ~isfile('all_pethsv2.mat')
        return
    end


    load("all_pethsv2.mat");
    s_out.peth = [s_out.peth; out.stim_peth];
    s_out.zpeth = [s_out.zpeth; out.stim_zpeth];
    s_out.shuf_peth = cat(3, s_out.shuf_peth,out.shuf_peth);
    s_out.shuf_zpeth = cat(3, s_out.shuf_zpeth,out.shuf_zpeth);
    % Norming trial-wise or row_wise
    row_min = min(out.shuf_zpeth, [], 2);
    row_max = max(out.shuf_zpeth,[],2);
    % Create matrices for broadcasting
    min_matrix = repmat(row_min, 1, size(out.shuf_zpeth, 2));
    max_matrix = repmat(row_max, 1, size(out.shuf_zpeth, 2));
    % Normalize
    normalized_data = (out.shuf_zpeth - min_matrix) ./ (max_matrix - min_matrix);
    s_out.norm_shuf_zpeth = cat(3,s_out.norm_shuf_zpeth, normalized_data);
    s_out.is_opto = [s_out.is_opto; out.is_opto];
    s_out.tvec = out.tvec;
end