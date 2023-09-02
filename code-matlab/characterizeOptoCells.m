%% Script to characterize opto_cells on the basis of firing rate changes and mean wave-form
% Assumes that stim_phases.mat and *stim_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.sham_delta, summary.opto_delta, ...
    summary.nostim_delta, summary.rev_delta, ...
    summary.sham_delta_sd, summary.opto_delta_sd, ...
    summary.nostim_delta_sd, summary.rev_delta_sd, ...
    summary.stim_poisson, summary.nostim_poisson, ...
    summary.stim_poisson_sd, summary.nostim_poisson_sd, ...
    summary.sham_delta_sd, summary.opto_delta_sd, ...
    summary.p_val, summary.h, summary.isopto, summary.depth, ...
    summary.waveforms, summary.depth] = deal([]);
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


%% Look at p_values
dStr_mask = summary.depth < 3.5;
ks_mask = summary.h == 1; % That the sham and opto difference in firing-rate is different
og_opto_mask = summary.isopto == 1;
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
subplot(2,4,[1 5])
hold on
sel = find(dStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'cyan');
end

sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'magenta');
end
xticks([1 2])
xticklabels({'Sham stim', 'Opto stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-50 350])
title(sprintf('dStr : %d', sum(dStr_mask)))

subplot(2,4,2)
hold on
sel1 = find(dStr_mask & ks_mask & summary.sham_delta < summary.opto_delta);
sel2 = find(dStr_mask & ks_mask & summary.sham_delta >= summary.opto_delta);
for i = 1:length(sel1)
    h = plot((0:1:31)/32, norm_wf(sel1(i),:), 'cyan');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
for i = 1:length(sel2)
    h = plot((0:1:31)/32, norm_wf(sel2(i),:), '--cyan');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot((0:1:31)/32,mean(norm_wf(sel1,:)), 'black', 'LineWidth', 3);
xlabel('Time (ms)')
ylabel('Normalized Amplitude')
legend({sprintf('Peak to trough mean: %.2f ms, sd: %.2f ms', ...
    mean(peak_to_trough(sel1))/32, std(peak_to_trough(sel1)/32))}, 'FontSize', 12)
title(sprintf('dStr Opto Responsive: %d', length(sel1)))

subplot(2,4,6)
hold on
sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    h = plot((0:1:31)/32, norm_wf(sel(i),:), 'magenta');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot((0:1:31)/32,mean(norm_wf(sel,:)), 'black', 'LineWidth', 3);
xlabel('Time (ms)')
ylabel('Normalized Amplitude')
legend({sprintf('Peak to trough mean: %.2f ms, sd: %.2f ms', ...
    mean(peak_to_trough(sel))/32, std(peak_to_trough(sel)/32))}, 'FontSize', 12)
title(sprintf('dStr Opto non-responsive: %d', length(sel)))

% Plot vStr stuff
subplot(2,4,[3 7])
hold on
sel = find(~dStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'cyan');
end

sel = find(~dStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'magenta');
end
xticks([1 2])
xticklabels({'Sham stim', 'Opto stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-50 350])
title(sprintf('vStr : %d', sum(~dStr_mask)))

subplot(2,4,4)
hold on
sel = find(~dStr_mask & ks_mask);
for i = 1:length(sel)
    h = plot((0:1:31)/32, norm_wf(sel(i),:), 'cyan');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot((0:1:31)/32,mean(norm_wf(sel,:)), 'black', 'LineWidth', 3);
xlabel('Time (ms)')
ylabel('Normalized Amplitude')
legend({sprintf('Peak to trough mean: %.2f ms, sd: %.2f ms', ...
    mean(peak_to_trough(sel))/32, std(peak_to_trough(sel)/32))}, 'FontSize', 12)
title(sprintf('vStr Opto Responsive: %d', length(sel)))

subplot(2,4,8)
hold on
sel = find(~dStr_mask & ~ks_mask);
for i = 1:length(sel)
    h = plot((0:1:31)/32, norm_wf(sel(i),:), 'magenta');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot((0:1:31)/32,mean(norm_wf(sel,:)), 'black', 'LineWidth', 3);
xlabel('Time (ms)')
ylabel('Normalized Amplitude')
legend({sprintf('Peak to trough mean: %.2f ms, sd: %.2f ms', ...
    mean(peak_to_trough(sel))/32, std(peak_to_trough(sel)/32))}, 'FontSize', 12)
title(sprintf('vStr Opto non-responsive: %d', length(sel)))

%%
dStr_opto = summary.labels(ks_mask & dStr_mask & summary.opto_delta > summary.sham_delta);
vStr_opto = summary.labels(ks_mask & ~dStr_mask & summary.opto_delta > summary.sham_delta);

%%
opto_p = [summary.labels(ks_mask),summary.h(ks_mask),summary.p_val(ks_mask), summary.sham_delta(ks_mask), summary.opto_delta(ks_mask)];
other_p = [summary.labels(~ks_mask),summary.h(~ks_mask), summary.p_val(~ks_mask), summary.sham_delta(~ks_mask), summary.opto_delta(~ks_mask)];

%% Diagnostic_plot: Compare sham stim to opto stim, reversed spike train stim and no_Stim sham duration
fig = figure('WindowState','maximized');

% Plot dStr stuff
subplot(2,2,1)
hold on
sel = find(dStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'cyan');
end
sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'magenta');
end
xticks([1 2])
xticklabels({'Sham stim', 'Opto stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-50 350])
title(sprintf('dStr : %d', sum(dStr_mask)))


subplot(2,2,2)
hold on
sel = find(dStr_mask & ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-2.5, summary.depth(sel(i))+2.5], ...
        [summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], ...
        [summary.sham_delta_sd(sel(i)), summary.opto_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'cyan')
end
sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-2.5, summary.depth(sel(i))+2.5], ...
        [summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], ...
        [summary.sham_delta_sd(sel(i)), summary.opto_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'magenta')
end
xticks([0.5 5.5])
xticklabels({'Sham stim', 'Opto stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([-0.25 6.25])
ylim([-50 350])
title(sprintf('dStr : %d', sum(dStr_mask)))

subplot(2,2,3)
hold on
sel = find(dStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.nostim_delta(sel(i)), summary.sham_delta(sel(i))], 'cyan');
end
sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.nostim_delta(sel(i)), summary.sham_delta(sel(i))], 'magenta');
end
xticks([1 2])
xticklabels({'No stim period', 'Sham stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-20 50])
title(sprintf('dStr : %d', sum(dStr_mask)))

subplot(2,2,4)
hold on
sel = find(dStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.rev_delta(sel(i)), summary.sham_delta(sel(i))], 'cyan');
end
sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.rev_delta(sel(i)), summary.sham_delta(sel(i))], 'magenta');
end
xticks([1 2])
xticklabels({'Rev. Spikes', 'Sham stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-20 50])
title(sprintf('dStr : %d', sum(dStr_mask)))
%%
fig = figure('WindowState','maximized');

% Plot vStr stuff
subplot(2,2,1)
hold on
sel = find(vStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'cyan');
end
sel = find(vStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], 'magenta');
end
xticks([1 2])
xticklabels({'Sham stim', 'Opto stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-50 350])
title(sprintf('vStr : %d', sum(vStr_mask)))


subplot(2,2,2)
hold on
sel = find(vStr_mask & ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-3.5, summary.depth(sel(i))+3.5], ...
        [summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], ...
        [summary.sham_delta_sd(sel(i)), summary.opto_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'cyan')
end
sel = find(vStr_mask & ~ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-3.5, summary.depth(sel(i))+3.5], ...
        [summary.sham_delta(sel(i)), summary.opto_delta(sel(i))], ...
        [summary.sham_delta_sd(sel(i)), summary.opto_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'magenta')
end
xticks([0.5 7.5])
xticklabels({'Sham stim', 'Opto stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([-0.25 8.75])
ylim([-50 350])
title(sprintf('vStr : %d', sum(vStr_mask)))

subplot(2,2,3)
hold on
sel = find(vStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.nostim_delta(sel(i)), summary.sham_delta(sel(i))], 'cyan');
end
sel = find(vStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.nostim_delta(sel(i)), summary.sham_delta(sel(i))], 'magenta');
end
xticks([1 2])
xticklabels({'No stim period', 'Sham stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-20 50])
title(sprintf('vStr : %d', sum(vStr_mask)))

subplot(2,2,4)
hold on
sel = find(vStr_mask & ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.rev_delta(sel(i)), summary.sham_delta(sel(i))], 'cyan');
end
sel = find(vStr_mask & ~ks_mask);
for i = 1:length(sel)
    plot([1,2],[summary.rev_delta(sel(i)), summary.sham_delta(sel(i))], 'magenta');
end
xticks([1 2])
xticklabels({'Rev. Spikes', 'Sham stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([0.85 2.15])
ylim([-20 50])
title(sprintf('vStr : %d', sum(vStr_mask)))

%% Diagnostic plot to compare Poisson train to sham
fig = figure('WindowState', 'maximized');
% Plot dStr stuff
subplot(2,2,1)
hold on
sel = find(dStr_mask & ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-2.5, summary.depth(sel(i))+2.5], ...
        [summary.stim_poisson(sel(i)), summary.sham_delta(sel(i))], ...
        [summary.stim_poisson_sd(sel(i)), summary.sham_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'cyan')
end
sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-2.5, summary.depth(sel(i))+2.5], ...
        [summary.stim_poisson(sel(i)), summary.sham_delta(sel(i))], ...
        [summary.stim_poisson_sd(sel(i)), summary.sham_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'magenta')
end
xticks([0.5 5.5])
xticklabels({'Stim Poisson', 'Sham stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([-0.25 6.25])
ylim([-50 350])
title(sprintf('dStr : %d', sum(dStr_mask)))

subplot(2,2,2)
hold on
sel = find(dStr_mask & ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-2.5, summary.depth(sel(i))+2.5], ...
        [summary.nostim_poisson(sel(i)), summary.sham_delta(sel(i))], ...
        [summary.nostim_poisson_sd(sel(i)), summary.sham_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'cyan')
end
sel = find(dStr_mask & ~ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-2.5, summary.depth(sel(i))+2.5], ...
        [summary.nostim_poisson(sel(i)), summary.sham_delta(sel(i))], ...
        [summary.nostim_poisson_sd(sel(i)), summary.sham_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'magenta')
end
xticks([0.5 5.5])
xticklabels({'No Stim Poisson', 'Sham stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([-0.25 6.25])
ylim([-50 350])
title(sprintf('dStr : %d', sum(dStr_mask)))

% Plot vStr stuff
subplot(2,2,3)
hold on
sel = find(vStr_mask & ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-3.5, summary.depth(sel(i))+3.5], ...
        [summary.stim_poisson(sel(i)), summary.sham_delta(sel(i))], ...
        [summary.stim_poisson_sd(sel(i)), summary.sham_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'cyan')
end
sel = find(vStr_mask & ~ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-3.5, summary.depth(sel(i))+3.5], ...
        [summary.stim_poisson(sel(i)), summary.sham_delta(sel(i))], ...
        [summary.stim_poisson_sd(sel(i)), summary.sham_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'magenta')
end
xticks([0.5 7.5])
xticklabels({'Stim Poisson', 'Sham stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([-0.25 8.75])
ylim([-50 350])
title(sprintf('vStr : %d', sum(vStr_mask)))

subplot(2,2,4)
hold on
sel = find(vStr_mask & ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-3.5, summary.depth(sel(i))+3.5], ...
        [summary.nostim_poisson(sel(i)), summary.sham_delta(sel(i))], ...
        [summary.nostim_poisson_sd(sel(i)), summary.sham_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'cyan')
end
sel = find(vStr_mask & ~ks_mask);
for i = 1:length(sel)
    errorbar([summary.depth(sel(i))-3.5, summary.depth(sel(i))+3.5], ...
        [summary.nostim_poisson(sel(i)), summary.sham_delta(sel(i))], ...
        [summary.nostim_poisson_sd(sel(i)), summary.sham_delta_sd(sel(i))], ...
        'MarkerFaceColor','black', 'Marker', 'o', 'Color', 'magenta')
end
xticks([0.5 7.5])
xticklabels({'Stim Poisson', 'Sham stim'})
ylabel('\Delta Firing rate (Hz)')
xlim([-0.25 8.75])
ylim([-50 350])
title(sprintf('vStr : %d', sum(vStr_mask)))


%% Making a poisson Spike train and looking at sham stim delta FR 
% (The raw firing rate change can non-zero, because it is getting multiplied by a factor of 10)
dt = 0.001;
t = [0 3000]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);
 
pspike = 0.5; % probability of generating a spike in bin
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
num_sham = 10000;
max_delay = 0.01; % In seconds
this_on_events = sort(randsample(tvec, num_sham));
fr = zeros(size(this_on_events));
bfr = zeros(size(this_on_events));
for iStim = 1:length(this_on_events)
    this_post = sum((spk_poiss_t > this_on_events(iStim)) & ...
        (spk_poiss_t <= (this_on_events(iStim) + max_delay)));
    this_pre = sum((spk_poiss_t <= this_on_events(iStim)) & ...
        (spk_poiss_t > (this_on_events(iStim) - max_delay)));
    fr(iStim) = this_post/max_delay;
    bfr(iStim) = this_pre/max_delay;
end

% Uncomment to plot
% diag_fig = figure;
% hist(fr-bfr);
% hold on
% xline(mean(fr-bfr))
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
        
        sham_dfr = od.sham_stim.fr' - od.sham_stim.bfr';
        nostim_dfr = od.no_stim.fr' - od.no_stim.bfr';
        rev_dfr = od.rev_stim.fr' - od.rev_stim.bfr';
        opto_idx = find(strcmp(S.label{iC},ExpKeys.goodCell));

%         if sum(nostim_dfr < 0) > 0
%             disp('About time baseline makes some sense\n');
%         end
% 
%         if sum(rev_dfr < 0) > 0
%             disp('About time Reverse makes some sense\n');
%         end
% 
%         if sum(sham_dfr < 0) > 0
%             disp('About time sham makes some sense\n');
%         end
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
        s_out.sham_delta = [s_out.sham_delta; mean(sham_dfr)];
        s_out.opto_delta = [s_out.opto_delta; mean(opto_dfr)];
        s_out.sham_delta_sd = [s_out.sham_delta_sd; std(sham_dfr)];
        s_out.opto_delta_sd = [s_out.opto_delta_sd; std(opto_dfr)];

        s_out.rev_delta = [s_out.rev_delta; mean(rev_dfr)];
        s_out.nostim_delta = [s_out.nostim_delta; mean(nostim_dfr)];
        s_out.rev_delta_sd = [s_out.rev_delta_sd; std(rev_dfr)];
        s_out.nostim_delta_sd = [s_out.nostim_delta_sd; std(nostim_dfr)];

        % generate dfr based on poisson spike train with trial_stim mfr
        dt = 0.001;
        t = [0 3000]; % time interval (length) of spike train to generate
        tvec = t(1):dt:t(2);
         
        pspike = dt*od.trial_stim.mfr; % based on actual mean_firing_rate
        spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
        spk_poiss_t = tvec(spk_poiss < pspike)'; % generate spike train
        num_sham = 10000;
        max_delay = 0.01; % In seconds
        this_on_events = sort(randsample(tvec, num_sham));
        fr = zeros(size(this_on_events));
        bfr = zeros(size(this_on_events));
        for iStim = 1:length(this_on_events)
            this_post = sum((spk_poiss_t > this_on_events(iStim)) & ...
                (spk_poiss_t <= (this_on_events(iStim) + max_delay)));
            this_pre = sum((spk_poiss_t <= this_on_events(iStim)) & ...
                (spk_poiss_t > (this_on_events(iStim) - max_delay)));
            fr(iStim) = this_post/max_delay;
            bfr(iStim) = this_pre/max_delay;
        end

        s_out.stim_poisson = [s_out.stim_poisson, mean(fr - bfr)];
        s_out.stim_poisson_sd = [s_out.stim_poisson, std(fr - bfr)];

        % generate dfr based on poisson spike train with pre_stim_baseline mfr
        dt = 0.001;
        t = [0 3000]; % time interval (length) of spike train to generate
        tvec = t(1):dt:t(2);
         
        pspike = dt*od.no_stim.mfr; % based on actual mean_firing_rate
        spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
        spk_poiss_t = tvec(spk_poiss < pspike)'; % generate spike train
        num_sham = 10000;
        max_delay = 0.01; % In seconds
        this_on_events = sort(randsample(tvec, num_sham));
        fr = zeros(size(this_on_events));
        bfr = zeros(size(this_on_events));
        for iStim = 1:length(this_on_events)
            this_post = sum((spk_poiss_t > this_on_events(iStim)) & ...
                (spk_poiss_t <= (this_on_events(iStim) + max_delay)));
            this_pre = sum((spk_poiss_t <= this_on_events(iStim)) & ...
                (spk_poiss_t > (this_on_events(iStim) - max_delay)));
            fr(iStim) = this_post/max_delay;
            bfr(iStim) = this_pre/max_delay;
        end

        s_out.nostim_poisson = [s_out.nostim_poisson, mean(fr - bfr)];
        s_out.nostim_poisson_sd = [s_out.nostim_poisson, std(fr - bfr)];

        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];

        fn_prefix = strrep(fn_prefix, '_', '-');
        % Load the mean waveform
        load(strcat(fn_prefix,'-wv.mat'));
        [~,midx]  = max(max(mWV)); % Getting the max amplitude channel waveform
        s_out.waveforms = [s_out.waveforms; mWV(:,midx)']; 
    end
end