%% Script to generate the relationships between spiking and LFP Phase
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.response_p, summary.bfr, summary.depth, summary.lat, summary.fr, ...
    summary.lat_r, summary.fr_r, summary.lat_z, summary.fr_z] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end

%% Sort by recording depth and fr_modulation_depth
[~, depth_sorted] = sort(summary.depth);
%%
delta_fr = squeeze(summary.fr(:,1,:));
norm_delta_fr = (delta_fr - min(delta_fr,[],2))./...
    repmat(max(delta_fr,[],2)-min(delta_fr,[],2),1,5);
imagesc(norm_delta_fr);

%% 
fbands = {[2 5], [6 10], [12 30], [30 55]};
c_list = {'red', 'blue','magenta', 'green'};
dStr_mask = (summary.depth < 3.5)';
%% Scatter plot of baseline firing rate vs depth of modulation
% dStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.bfr(dStr_mask), summary.fr_r(dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 200);
    hold  on; 
    scatter(summary.bfr(dStr_mask & this_ex_sig), summary.fr_r(dStr_mask & this_ex_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 50);
%     plot([-5 5],[-5 5], 'Color', c_list{iF})
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    xlabel('Baseline Firing Rate', 'FontSize', 12) 
    ylabel('Depth of Modulation', 'FontSize', 12)
    xlim([-5 35])
    ylim([-0.1 1.1])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
    sgtitle('Dorsal Striatum')
end
%% vStr
fig = figure('WindowState', 'maximized');
for iF = 1:length(fbands)
    this_ex_sig = summary.fr_z(:,iF) > 2;
    
    ax = subplot(2,2,iF);
    scatter(summary.bfr(~dStr_mask), summary.fr_r(~dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.5, 'SizeData', 200);
    hold  on; 
    scatter(summary.bfr(~dStr_mask & this_ex_sig), summary.fr_r(~dStr_mask & this_ex_sig,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SizeData', 50);
%     plot([-5 5],[-5 5], 'Color', c_list{iF})
    legend({'all ', 'Sig. Phase Modulated'}, 'Location', 'best', 'FontSize', 12);
    xlabel('Baseline Firing Rate', 'FontSize', 12) 
    ylabel('Depth of Modulation', 'FontSize', 12)
    xlim([-5 35])
    ylim([-0.1 1.1])
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
    sgtitle('Ventral Striatum')
end

%% Scatter plot of Z-scores
ax = subplot(2,2,1);
hold on
for iF = 1:length(fbands)
    scatter(summary.response_p(dStr_mask), summary.fr_z(dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 100);
end
ylim([-5 15])
yline(2)
title('Dorsal Striatum');
xlabel('Proportion of trials with any response')
ylabel('zscore')
legend({'Delta (2Hz - 5Hz)', 'Theta (6Hz - 10Hz)', 'Beta (12Hz - 30Hz)', ...
    'Low Gamma (30Hz - 55Hz)'}, 'Location','northwest');

ax = subplot(2,2,2);
hold on
for iF = 1:length(fbands)
        scatter(summary.response_p(~dStr_mask), summary.fr_z(~dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SizeData', 100);
end
ylim([-5 15])
yline(2)
title('Ventral Striatum');
xlabel('Proportion of trials with any response')
ylabel('zscore')
legend({'Delta (2Hz - 5Hz)', 'Theta (6Hz - 10Hz)', 'Beta (12Hz - 30Hz)', ...
    'Low Gamma (30Hz - 55Hz)'}, 'Location','northwest');

ax = subplot(2,2,3);
hold on
for iF = 1:length(fbands)
    scatter(summary.bfr(dStr_mask), summary.fr_z(dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Marker', 'd', 'SizeData', 100);
end
ylim([-5 15])
yline(2)
title('Dorsal Striatum');
xlabel('Mean Baseline Firing Rate')
ylabel('zscore')
legend({'Delta (2Hz - 5Hz)', 'Theta (6Hz - 10Hz)', 'Beta (12Hz - 30Hz)', ...
    'Low Gamma (30Hz - 55Hz)'}, 'Location','northwest');

ax = subplot(2,2,4);
hold on
for iF = 1:length(fbands)
        scatter(summary.bfr(~dStr_mask), summary.fr_z(~dStr_mask,iF), 'MarkerFaceColor', c_list{iF}, ...
        'MarkerEdgeColor', c_list{iF}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'SizeData', 100);
end
ylim([-5 15])
yline(2)
title('Ventral Striatum');
xlabel('Mean Baseline Firing Rate')
ylabel('zscore')
legend({'Delta (2Hz - 5Hz)', 'Theta (6Hz - 10Hz)', 'Beta (12Hz - 30Hz)', ...
    'Low Gamma (30Hz - 55Hz)'}, 'Location','northwest');




% The scatter plot looks too busy, think of alternate figures?
% Bars with each 
%% Scatter of pi
dStr_sig_only_1 = sum(summary.fr_z(dStr_mask,:) > 2)/sum(dStr_mask);
vStr_sig_only_1 = sum(summary.fr_z(~dStr_mask,:) > 2)/sum(~dStr_mask);

q0 = summary.fr_z(dStr_mask,:) > 2;
q0 = sum(q0,2);
dStr_no_sig = sum(q0<1)/length(q0);

q0 = summary.fr_z(dStr_mask,:) > 2;
q0 = sum(q0,2);
dStr_more_sig = sum(q0>1)/length(q0);

q0 = summary.fr_z(~dStr_mask,:) > 2;
q0 = sum(q0,2);
vStr_no_sig = sum(q0<1)/length(q0);

q0 = summary.fr_z(~dStr_mask,:) > 2;
q0 = sum(q0,2);
vStr_more_sig = sum(q0>1)/length(q0);
%%
ax = subplot(1,2,1);
pie(dStr_sig_only_1, {sprintf('2Hz - 5Hz, %.2f%%',100*dStr_sig_only_1(1)), ...
   sprintf('6Hz - 10Hz, %.2f%%',100*dStr_sig_only_1(2)), sprintf('12Hz - 30Hz, %.2f%%',100*dStr_sig_only_1(3)), ...
  sprintf('30Hz - 55Hz, %.2f%%',100*dStr_sig_only_1(4))});
ax.Colormap  = [1 0 0; 0 0 1; 1 0 1; 0 1 0];
title(sprintf('Dorsal Striatum (%d)', sum(dStr_mask)));

ax = subplot(1,2,2);
pie(vStr_sig_only_1, {sprintf('2Hz - 5Hz, %.2f%%',100*vStr_sig_only_1(1)), ...
   sprintf('6Hz - 10Hz, %.2f%%',100*vStr_sig_only_1(2)), sprintf('12Hz - 30Hz, %.2f%%',100*vStr_sig_only_1(3)), ...
  sprintf('30Hz - 55Hz, %.2f%%',100*vStr_sig_only_1(4))});
ax.Colormap  = [1 0 0; 0 0 1; 1 0 1; 0 1 0];
title(sprintf('Ventral Striatum (%d)', sum(~dStr_mask)));


%% Scatter plot of Z-scores 
fbands = {[2 5], [6 10], [12 30], [30 55]};
c_list = {'red', 'blue','magenta', 'green'};
dStr_mask = summary.depth < 3.5;
ax = subplot(2,1,[1 2]);
hold on
for iF = 1:length(fbands)
    text(summary.response_p(dStr_mask), summary.fr_z(dStr_mask,iF), summary.labels(dStr_mask), 'Color', 'red');
    text(summary.response_p(~dStr_mask), summary.fr_z(~dStr_mask,iF), summary.labels(~dStr_mask), 'Color', 'blue');
end
yline(2)
% The scatter plot looks too busy, think of alternate figures?
% Bars with each 
%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        % Load the phase_response
        load(strcat(fn_prefix, '_phase_response_5_bins.mat')); % Change this to what is decided to be the best binning option
        s_out.depth = [s_out.depth, ExpKeys.probeDepth];
        s_out.response_p = [s_out.response_p; out.overall_response];
        s_out.bfr = [s_out.bfr, out.mean_bfr];
        s_out.labels = [s_out.labels; string(fn_prefix)];
%         s_out.lat(size(s_out.lat,1)+1,:,:) = out.lat.bin;
%         s_out.lat_z = [s_out.lat_z; out.lat.zscore];
%         s_out.lat_r = [s_out.lat_r; out.lat.ratio];
%         s_out.fr(size(s_out.fr,1)+1,:,:) = out.fr.bin;
        s_out.fr_z = [s_out.fr_z; out.fr.zscore];
        s_out.fr_r = [s_out.fr_r; out.fr.ratio];
    end
end
