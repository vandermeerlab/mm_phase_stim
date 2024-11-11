%% Script to make plots that differentate between the opto-responsive MSNs and other opto-responsive FSIs
% Assumes that stim_phases.mat and *stim_response.mat already exist in each folder
top_dir = 'E:\Dartmouth College Dropbox\Manish Mohapatra\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.opto_delta, summary.sham_delta, ...
    summary.opto_delta_sd, summary.sham_delta_sd, summary.p_val, ...
    summary.h, summary.p_val_t, summary.h_t, ...
    summary.isopto, summary.depth, summary.waveforms, ...
    summary.depth, summary.stim_peth, summary.stim_tvec] = deal([]);
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
%%
msn_mask = ((summary.labels == 'M295-2022-01-06-TT06_4')| ...
    (summary.labels == 'M295-2022-01-06-TT08_4'));
opto_fsi = ks_mask & ~msn_mask;
%%
fig = figure('WindowState','maximized');
ax = subplot(2,3,1);
hold on

mean_fsi_wv = mean(norm_wf(opto_fsi,:),1);
plot((0:1:31)/32, mean_fsi_wv, 'cyan');
sel2 = find(msn_mask);
for i = 1:length(sel2)
    h = plot((0:1:31)/32, norm_wf(sel2(i),:), 'magenta');
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
title('Spike waveforms')
legend({'Mean opto-FSI', 'Opto-MSN 1', 'Opto-MSN 2' })
ax = subplot(2,3,[2,3]);
hold on;
mean_fsi_peth = mean(summary.stim_peth(opto_fsi,:));
plot(summary.stim_tvec(1,:), mean_fsi_peth, 'cyan');
for i = 1:length(sel2)
    h = plot(summary.stim_tvec(1,:), summary.stim_peth(sel2(i),:), 'magenta');
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
title('Normalized PETH')
legend({'Mean opto-FSI', 'Opto-MSN 1', 'Opto-MSN 2' })
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

        % Comapare to a similar distribution as before
        ridx = randi(1000000,1,10000);
        sham_dfr = od.sham_stim.fr(ridx)' - od.sham_stim.bfr(ridx)';
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

        s_out.opto_delta = [s_out.opto_delta; mean(opto_dfr)];
        s_out.opto_delta_sd = [s_out.opto_delta_sd; std(opto_dfr)];
 
        s_out.sham_delta = [s_out.sham_delta; mean(sham_dfr)];
        s_out.sham_delta_sd = [s_out.sham_delta_sd; std(sham_dfr)];

        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];

        fn_prefix = strrep(fn_prefix, '_', '-');
        % Load the mean waveform
        load(strcat(fn_prefix,'-wv.mat'));
        [~,midx]  = max(max(mWV)); % Getting the max amplitude channel waveform
        s_out.waveforms = [s_out.waveforms; mWV(:,midx)'];

        % Code to get PETH
        if contains(ExpKeys.light_source, 'LASER')
            start_delay = 0.0011;
            stop_delay = 0.0012; %0.0022 is the spike width in the case of 1 msec laser pulse
        else
            start_delay = 0;
            stop_delay = 0;
        end
        
        evs = LoadEvents([]);
        stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
        if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
            stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                              stim_on <= ExpKeys.stim_times(2));
        else
            stim_on = [];
        end

        % Load CSC
        cfg = []; cfg.fc = ExpKeys.goodLFP;
        if contains(cfg.fc, '-')
            temp = split(cfg.fc,'-');
            cfg.fc = {cat(2,temp{1},'.ncs')};
            csc = LoadCSC(cfg);
            cfg_temp.fc = {cat(2,temp{2},'.ncs')};
            ref = LoadCSC(cfg_temp);
            csc.data = csc.data - ref.data;    
            clear temp ref;
        else
            csc = LoadCSC(cfg);
        end
        
        % Restrict CSC to the stim_period
        csc = restrict(csc, iv(ExpKeys.stim_times(1), ExpKeys.stim_times(2)));

        % Sometimes csc.tvec can be have weird elements because of gaps in recording
        to_remove = find(diff(csc.tvec)<=0);
        while (~isempty(to_remove))
            csc.tvec(to_remove+1) = [];
            csc.data(to_remove+1) = [];
            to_remove = find(diff(csc.tvec)<=0);
        end
        
        % Calculate MUA
        cfg_MUA = []; 
        cfg_MUA.tvec = csc.tvec; % timebase to compute MUA on
        cfg_MUA.sigma = 0.001;

        % Parameters for PETH
        cfg_peth = []; % parameters for PETH
        cfg_peth.window = [-0.02 0.02];
        cfg_peth.dt = 0.001;
        cfg_peth.mode = 'interp';

        this_cell = SelectTS([],S,iC);
        this_MUA  = getMUA(cfg_MUA, this_cell);
        this_MUAz = zscore_tsd(this_MUA);
        this_out = TSDpeth_fast(cfg_peth, this_MUAz, stim_on);
        s_out.stim_peth = [s_out.stim_peth; (this_out.data-min(this_out.data))/(max(this_out.data) - min(this_out.data))];
        s_out.stim_tvec = [s_out.stim_tvec; this_out.tvec];
    end
end