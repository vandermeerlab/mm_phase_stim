%% Script to characterize opto_cells on the basis of firing rate changes and mean wave-form
% Assumes that stim_phases.mat and *stim_response.mat already exist in each folder
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.nostim_delta, summary.opto_delta, ...
    summary.sham_delta, summary.poisson_delta, summary.nostim_delta_sd, summary.opto_delta_sd, ...
    summary.sham_delta_sd, summary.poisson_delta_sd, ...
    summary.p_val, summary.h, summary.p_val_t, summary.h_t, ...
    summary.isopto, summary.depth, summary.waveforms, ...
    summary.depth] = deal([]);
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
sel = find(dStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'blue');
end

sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'magenta');
end
xticks([1 2])
yticks([0 160 320])
xticklabels({'Sham stim', 'Opto stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-25 320])
title('dStr');
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';

ax = subplot(2,4,2);
hold on
sel1 = find(dStr_mask & ks_mask & summary.sham_delta < summary.opto_delta);
sel2 = find(dStr_mask & ks_mask & summary.sham_delta >= summary.opto_delta);
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
legend({sprintf('Peak to trough %.2f +/- %.2f ms', ...
    mean(peak_to_trough(sel1))/32, std(peak_to_trough(sel1)/32))}, 'FontSize', 12)
title('Opto Responsive')
xticks([0 0.5 1]);
yticks([0 0.5 1]);
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';

ax = subplot(2,4,6);
hold on
sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    h = plot((0:1:31)/32, norm_wf(sel(i),:), 'magenta');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot((0:1:31)/32,mean(norm_wf(sel,:)), 'black', 'LineWidth', 3);
xlabel('Time (ms)')
ylabel('Normalized Amplitude')
legend({sprintf('Peak to trough %.2f +/- %.2f ms', ...
    mean(peak_to_trough(sel))/32, std(peak_to_trough(sel)/32))}, 'FontSize', 12)
title('Opto non-responsive')
xticks([0 0.5 1]);
yticks([0 0.5 1]);
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';

% Plot vStr stuff
ax = subplot(2,4,[3 7]);
hold on
sel = find(vStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'blue');
end

sel = find(vStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'magenta');
end
xticks([1 2])
xticklabels({'Sham stim', 'Opto stim'})
ylabel('\Delta Firing rate (Hz)')
yticks([0 160 320])
xlim([0.85 2.15])
ylim([-25 320])
title('vStr')
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';

ax = subplot(2,4,4);
hold on
sel = find(vStr_mask & ks_mask);
for i = 1:length(sel)
    h = plot((0:1:31)/32, norm_wf(sel(i),:), 'blue');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot((0:1:31)/32,mean(norm_wf(sel,:)), 'black', 'LineWidth', 3);
xlabel('Time (ms)')
ylabel('Normalized Amplitude')
legend({sprintf('Peak to trough mean: %.2f +/- %.2f ms', ...
    mean(peak_to_trough(sel))/32, std(peak_to_trough(sel)/32))}, 'FontSize', 12)
title('Opto Responsive')
xticks([0 0.5 1]);
yticks([0 0.5 1]);
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';

ax = subplot(2,4,8);
hold on
sel = find(vStr_mask & ~ks_mask);
for i = 1:length(sel)
    h = plot((0:1:31)/32, norm_wf(sel(i),:), 'magenta');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot((0:1:31)/32,mean(norm_wf(sel,:)), 'black', 'LineWidth', 3);
xlabel('Time (ms)')
ylabel('Normalized Amplitude')
legend({sprintf('Peak to trough %.2f +/- %.2f ms', ...
    mean(peak_to_trough(sel))/32, std(peak_to_trough(sel)/32))}, 'FontSize', 12)
title('Opto non-responsive')
xticks([0 0.5 1]);
yticks([0 0.5 1]);
box off;
ax.TickLength(1) = 0.03;
ax.TickDir = 'out';

fontsize(fig, 30, 'points');
fontname(fig, 'Helvetica');
fig.Renderer = 'painters';

%% Use this to write list of Final Opto Cells
dStr_opto = summary.labels(ks_mask & dStr_mask & summary.opto_delta > summary.sham_delta);
vStr_opto = summary.labels(ks_mask & vStr_mask & summary.opto_delta > summary.sham_delta);

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
        
        % Load the stim_response
        load(strcat(fn_prefix,'_stim_response.mat'));

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
end