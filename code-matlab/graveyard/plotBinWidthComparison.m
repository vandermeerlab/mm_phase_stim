%% Script to generate the relationships between spiking and LFP Phase
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', ...
    'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.response_p, summary.depth, summary.fr_r, summary.fr_z] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end

%% Plot binning comparison results
fbands = {[2 5], [6 10], [12 30], [30 55]};
c_list = {'cyan', 'red','magenta', 'green'};
nshufs = 100;
bin_widths = [20, 15, 10, 5]; % in msec

for iF = 1:length(fbands)
    ax = subplot(3,4,iF); 
    hold on
    for iC = 1:size(summary.fr_r,1)
        plot(bin_widths, summary.fr_r(iC,:,iF), 'Color', c_list{iF});
    end
    mean_trend = squeeze(mean(summary.fr_r, 1));
    plot(bin_widths, mean_trend(:,iF), 'LineWidth', 2, 'Color', 'black');
    xlabel('Bin Width (msec)')
    ylabel('Depth of modulation')
    ax.Title.String = sprintf('%d Hz - %d Hz', fbands{iF}(1), fbands{iF}(2));
    ax.XAxis.FontSize = 13;
    ax.YAxis.FontSize = 13;
    ax.XDir = 'reverse';

    ax = subplot(3,4,iF+4); 
    hold on
    for iC = 1:size(summary.fr_z,1)
        plot(bin_widths, summary.fr_z(iC,:,iF), 'Color', c_list{iF});
    end
    mean_trend = squeeze(nanmean(summary.fr_z, 1));
    plot(bin_widths, mean_trend(:,iF), 'LineWidth', 2, 'Color', 'black');
    xlabel('Bin Width (msec)')
    ylabel('Zscore')
    ax.XAxis.FontSize = 13;
    ax.YAxis.FontSize = 13;
    ax.XDir = 'reverse';

    ax = subplot(3,4,iF+8);
    this_fr_z = squeeze(summary.fr_z(:,:,iF));
    this_sig = sum(this_fr_z >= 2);
    this_nan = sum(isnan(this_fr_z));
    for iN = 1:length(this_sig)
        text(0.2*(iN-1),0.8, string(this_sig(iN)), 'FontSize', 16);
        text(0.2*(iN-1),0.5, string(this_nan(iN)), 'FontSize', 16);
        text(0.2*(iN-1),0.2, string(ceil((1000/mean(fbands{iF}))/bin_widths(iN))), 'FontSize', 16);
    end
    yticks([0.2, 0.5, 0.8])
    yticklabels({'Bin Count', 'Ineligible', 'Count'});
    xticks([0:0.2:0.6])
    xticklabels({'20', '15', '10', '5'})
    xlim([-0.1 0.9])
    xlabel('Bin Width (msec)')
    ax.Box = 'off';
end   

%%
function s_out = doStuff(s_in)
    s_out = s_in;
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    bin_widths = [20, 15, 10, 5]; % in msec
    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        this_fr_r = zeros(length(bin_widths),4);
        this_fr_z = zeros(length(bin_widths),4);
        this_response = 0;
        for iN = 1:length(bin_widths)
            load(strcat(fn_prefix, strcat('_phase_response_', string(bin_widths(iN)),'ms_bins.mat')));
            this_response = out.overall_response;
            this_fr_r(iN,:) = out.fr.ratio;
            this_fr_z(iN,:) = out.fr.zscore;
        end
        s_out.depth = [s_out.depth; ExpKeys.probeDepth];
        s_out.response_p = [s_out.response_p; this_response];
        s_out.labels = [s_out.labels; string(fn_prefix)];
        s_out.fr_r(size(s_out.fr_r,1)+1,:,:) = this_fr_r;
        s_out.fr_z(size(s_out.fr_z,1)+1,:,:) = this_fr_z;
    end
end
