%%Plot PETH as a function of phase
load('E:\temp_phase_stim\ED\M19-2019-04-15_dStr_4p0_light_cells_TT6_min\latencies\M19_2019_04_15_lat.mat');
c_ord = linspecer(5);
set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesFontWeight','bold')


%% Values for 2-5 Hz
close all;
this_lat = out.M19_2019_04_15_TT06_01_Good.f_2_5.latency;
clean_lat = this_lat(:,~isnan(this_lat(2,:)));
bin_wise_latency = cell(1,5);
fig = figure;
hold on;

for i = 1:length(bin_wise_latency)
    bin_wise_latency{i} = clean_lat(2,find(clean_lat(1,:) == i));
    [values, edges] = histcounts(bin_wise_latency{i}, [-1:0.05:6],'Normalization', 'count');
    centers = (edges(1:end-1)+edges(2:end))/2;
%     plot(centers, values*(length(bin_wise_latency{i})/size(clean_lat,2)), 'Color', c_ord(i,:), 'LineWidth', 3);
    plot(centers, values/size(clean_lat,2), 'Color', c_ord(i,:), 'LineWidth', 3)
end

ax = gca;
ax.Box = 'off';
ax.TickDir = 'out';
ax.Box = 'off';
ax.TickDir = 'out';
ax.Box = 'off';
ax.TickDir = 'out';
ax.XLabel.String = 'Time (msec)';
ax.YLabel.String = 'Proportion of first responses to Opto stim';
ax.XAxis.FontSize = 32;
ax.YAxis.FontSize = 32;
ax.YLabel.FontSize = 24;
ax.XLabel.FontSize = 24;

%% code to plot p_spike
b = bar(cellfun(@(x) length(x)/size(clean_lat,2), bin_wise_latency));
b.BarWidth = 0.98;
ax = gca;
ax.Box = 'off';
ax.TickDir = 'out';
ax.XLabel.String = 'Phase Bins';
ax.YLabel.String = 'Proportion of responses';
ax.YAxis.FontSize = 32;
ax.XAxis.FontSize = 32;
ax.YLabel.FontSize = 24;
ax.XLabel.FontSize = 24;
ax.XTick = {};
ax.YTick = [0 0.1 0.2 0.3];


%%
