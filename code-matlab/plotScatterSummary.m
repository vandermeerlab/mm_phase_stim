%% Script to generate the relationships between spiking and LFP Phase
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.response_p, summary.bfr, summary.depth,...
    summary.fr_r, summary.fr_z,  summary.phaselock_binned, summary.phaselock_phase, ...
    summary.excitable_phase, summary.phaselock_sig, summary.nbins] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end
fbands = {[2 5], [6 10], [12 30], [30 55]};
c_list = {'red', 'blue','magenta', 'green'};
dStr_mask = (summary.depth < 3.5)';

%% Sort by recording depth and fr_modulation_depth
[~, depth_sorted] = sort(summary.depth);

%% Scatter plot of exciatable phase vs intrinsic phase
%Plot for Dorsal Striatum
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_pl_sig = summary.phaselock_sig(:,iF) == 1;
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.phaselock_phase(dStr_mask,iF), summary.excitable_phase(dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on; 
    scatter(summary.phaselock_phase(dStr_mask & this_pl_sig,iF), summary.excitable_phase(dStr_mask & this_pl_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 25);
    scatter(summary.phaselock_phase(dStr_mask & this_ex_sig,iF), summary.excitable_phase(dStr_mask & this_ex_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0.5, 'Marker', '+','SizeData', 200);
    plot([-5 5],[-5 5], 'Color', c_list{iF})
    legend({'all ', 'Sig. Phase Locked', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    xlabel('Intrinsic Phase', 'FontSize', 12) 
    ylabel('Most Excitable Phase', 'FontSize', 12)
    xlim([-3.2 3.2])
    ylim([-3.2 3.2])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Dorsal Striatum')

%% Plot for Ventral Striatum
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_pl_sig = summary.phaselock_sig(:,iF) == 1;
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.phaselock_phase(~dStr_mask,iF), summary.excitable_phase(~dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
    hold  on; 
    scatter(summary.phaselock_phase(~dStr_mask & this_pl_sig,iF), summary.excitable_phase(~dStr_mask & this_pl_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SizeData', 25);
    scatter(summary.phaselock_phase(~dStr_mask & this_ex_sig,iF), summary.excitable_phase(~dStr_mask & this_ex_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0.5, 'Marker', '+','SizeData', 200);
    plot([-5 5],[-5 5], 'Color', c_list{iF})
    legend({'all ', 'Sig. Phase Locked', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    xlabel('Intrinsic Phase', 'FontSize', 12) 
    ylabel('Most Excitable Phase', 'FontSize', 12)
    xlim([-3.2 3.2])
    ylim([-3.2 3.2])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Ventral Striatum')
%% Scatter of most excitable phases for significant results
% dStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.excitable_phase(dStr_mask,iF), summary.fr_r(dStr_mask,iF),'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on; 
    scatter(summary.excitable_phase(dStr_mask & this_ex_sig,iF), summary.fr_r(dStr_mask & this_ex_sig,iF),  'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 25);
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    ylabel('Depth of modulation', 'FontSize', 12) 
    xlabel('Most Excitable Phase', 'FontSize', 12)
    ylim([-0.2 1.2])
    xlim([-3.2 3.2])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Dorsal Striatum')
%% vStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.excitable_phase(~dStr_mask,iF), summary.fr_r(~dStr_mask,iF),'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on; 
    scatter(summary.excitable_phase(~dStr_mask & this_ex_sig,iF), summary.fr_r(~dStr_mask & this_ex_sig,iF),  'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 25);
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    ylabel('Depth of modulation', 'FontSize', 12) 
    xlabel('Most Excitable Phase', 'FontSize', 12)
    ylim([-0.2 1.2])
    xlim([-3.2 3.2])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
end
sgtitle('Ventral Striatum')

%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    p_thresh = 0.99;
    fbands = {[2 5], [6 10], [12 30], [30 55]};
    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        % Load the phase_response
        load(strcat(fn_prefix, '_phase_response_10ms_bins.mat')); % Change this to what is decided to be the best binning option
        s_out.depth = [s_out.depth, ExpKeys.probeDepth];
        s_out.response_p = [s_out.response_p; out.overall_response];
        s_out.bfr = [s_out.bfr, out.mean_bfr];
        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.fr_z = [s_out.fr_z; out.fr.zscore];
        s_out.fr_r = [s_out.fr_r; out.fr.ratio];
        
        % Load phase_lock and shuf_spec
        fn_prefix = strrep(fn_prefix, '_', '-');
        load(strcat(fn_prefix, '_spike_phaselock.mat'));
        load(strcat(fn_prefix, '_shuf_spec.mat'));
           
        pct_ppc = sum(shuf_ppc < trial_ppc.vals', 1)/size(shuf_ppc,1); % Calculating p_val of PPC for each frequency
        [this_sig, this_pl_phase, this_pl_binned, this_nbins, this_ex_phase] = deal(zeros(1,length(fbands)));
        for iF = 1:length(fbands)
            f_idx = find(round(trial_sts.freqs) >= fbands{iF}(1) & ...
                round(trial_sts.freqs) <= fbands{iF}(2));
            mean_pct = mean(pct_ppc(f_idx));
            this_sig(iF) = mean_pct >= p_thresh;

            % Change accordingly
            this_nbins(iF) = length(out.fr.bin{iF}); %5;
            phase_bins = -pi:2*pi/this_nbins(iF):pi;
           
            % The mean angle at the maximum PPC
            [~, midx] = max(trial_ppc.vals(f_idx));
            this_ang = trial_ppc.ang(f_idx);
            this_pl_phase(iF) = this_ang(midx);
            [~, ~, this_bin] = histcounts(this_pl_phase(iF), phase_bins);
            this_pl_binned(iF) = mean(phase_bins(this_bin:this_bin+1)); % Bin the phase accordingly

            % The maximally excitable phase
            [~, midx] =  max(out.fr.bin{iF});% max(out.fr.bin);
            this_ex_phase(iF) = mean(phase_bins(midx:midx+1));
        end
        s_out.phaselock_sig = [s_out.phaselock_sig; this_sig];
        s_out.phaselock_phase = [s_out.phaselock_phase; this_pl_phase];
        s_out.phaselock_binned = [s_out.phaselock_binned, this_pl_binned];
        s_out.excitable_phase = [s_out.excitable_phase; this_ex_phase];
        s_out.nbins = [s_out.nbins; this_nbins];
    end
end
