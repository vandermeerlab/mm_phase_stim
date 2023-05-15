clear;
close all;
cd('E:\Dropbox (Dartmouth College)\manish_data\PhotoSensor Tests\PhotoSensor-2021-03-12')
evs = LoadEvents([]);
cfg.fc = {'PhotoSensor.ncs'};
photo_signal = LoadCSC(cfg);
laser_idx = cellfun(@(x) contains(x, 'Laser'), evs.label, 'UniformOutput', 1);
led_idx = cellfun(@(x) contains(x, 'LED'), evs.label, 'UniformOutput', 1);
% Correct LED events were after the event 'Final'(evs.label{2})
after_final = cellfun(@(x) x(end) > evs.t{2}, evs.t, 'UniformOutput', 1);
led_idx = after_final & led_idx;
laser_labels = evs.label(laser_idx);
led_labels = evs.label(led_idx);

pulse_OFF_string = evs.label{52};
pulse_ON_string = evs.label{53};

pulse_OFF_times = evs.t{52};
pulse_ON_times = evs.t{53};

%% LED OFF train Analysis

led_temp = cellfun(@(x) split(x), led_labels, 'UniformOutput', false);
led_pw = cellfun(@(x) str2double(x{2}), led_temp);
led_signal_width = [];
led_start_delay = [];
led_end_delay = [];
[led_pw, led_sort_idx] = sort(led_pw);

led_labels = led_labels(led_sort_idx);
for i = 1:length(led_labels)
    start_idx = find(cellfun(@(x) contains(x,led_labels{i}), evs.label) == 1);
    start_time = evs.t{start_idx}(end); % 'end' because only the final LED label is correct
    if i ~= length(led_labels)
        end_idx = find(cellfun(@(x) contains(x,led_labels{i+1}), evs.label) == 1);
        end_time = evs.t{end_idx}(end);
    else
        end_time = evs.t{51}(end); % Because the correct set of LED recordings was the last
    end
    this_led_ON_times = pulse_ON_times(pulse_ON_times > start_time & pulse_ON_times < end_time);
    this_led_OFF_times= pulse_OFF_times(pulse_OFF_times > start_time & pulse_OFF_times < end_time);
    % We need to analyze the parts between off(1:end-1) to ON(2:end)
    this_led_ON_times = this_led_ON_times(2:end);
    this_led_OFF_times = this_led_OFF_times(1:end-1);
    snippet_starts = nearest_idx3(this_led_OFF_times, photo_signal.tvec);
    plot(photo_signal.tvec(snippet_starts(1)-1000:snippet_starts(end)+1000), ...
        photo_signal.data(snippet_starts(1)-1000:snippet_starts(end)+1000));
    title(led_labels{i})
    close;
end

%%
% start from 4 because < 800 msec looks weird
for i = 4:length(led_labels)
    start_idx = find(cellfun(@(x) contains(x,led_labels{i}), evs.label) == 1);
    start_time = evs.t{start_idx}(end); % 'end' because only the final LED label is correct
    if i ~= length(led_labels)
        end_idx = find(cellfun(@(x) contains(x,led_labels{i+1}), evs.label) == 1);
        end_time = evs.t{end_idx}(end);
    else
        end_time = evs.t{51}(end); % Because the correct set of LED recordings was the last
    end
    this_led_ON_times = pulse_ON_times(pulse_ON_times > start_time & pulse_ON_times < end_time);
    this_led_OFF_times= pulse_OFF_times(pulse_OFF_times > start_time & pulse_OFF_times < end_time);
    this_led_ON_times = this_led_ON_times(2:end);
    this_led_OFF_times = this_led_OFF_times(1:end-1);
    snippet_starts = nearest_idx3(this_led_OFF_times - 1e-2, photo_signal.tvec);
    snippet_ends = nearest_idx3(this_led_ON_times + 1e-2, photo_signal.tvec);
    snip_width = max(snippet_ends - snippet_starts);
    this_pw = this_led_ON_times - this_led_OFF_times;
    pulse_starts = nearest_idx3(this_led_OFF_times, photo_signal.tvec);
    pulse_ends = nearest_idx3(this_led_ON_times, photo_signal.tvec);
    pulse_starts = pulse_starts - snippet_starts;
    pulse_ends = pulse_ends - snippet_starts;
    
    % Calculating the average signal because there isn't much variation
    average_snippet = zeros(1,snip_width);
    for j = 1:length(snippet_starts)
        average_snippet = average_snippet + ...
            photo_signal.data(snippet_starts(j):snippet_starts(j)+snip_width-1);
    end
    average_snippet = (average_snippet/length(snippet_ends))+0.25; %Adding 0.25 to shift the whole thing above 0
    
    % Defining the actual signal start as where the signal reaches 95% of 
    % its max value and the signal end as where the signal drops below the
    % 95% of its max value

    threshold = 0.95*max(average_snippet);
    actual_start = find(average_snippet>threshold, 1);
    actual_width = find(average_snippet(actual_start+1:end) < threshold,1)+1;
    actual_end = actual_start+actual_width;
    actual_start_delay = actual_start - mean(pulse_starts);
    actual_end_delay = actual_end - mean(pulse_ends);
    led_signal_width = [led_signal_width actual_width];
    led_start_delay = [led_start_delay actual_start_delay];
    led_end_delay = [led_end_delay actual_end_delay];
    
    fig = figure('WindowState', 'maximized');
    hold on
    plot(average_snippet);
    xlabel("Time (microseconds)")
    ylabel("PhotoSensor Signal")
    vline(pulse_starts(1), 'r');
    vline(pulse_ends(1), '--r');
    vline(actual_start, 'g');
    vline(actual_end, '--g');
    
    txt = sprintf("Actual Width = %.2f usec", actual_width*1000000/32000);
    text(600, 0.245, txt, 'FontSize', 12);
    txt = sprintf("Start Delay = %.2f usec", actual_start_delay*1000000/32000);
    text(600, 0.235, txt, 'FontSize', 12);
    txt = sprintf("End Delay = %.2f usec", actual_end_delay*1000000/32000);
    text(600, 0.225, txt, 'FontSize', 12);
    txt = sprintf("LED Pulse Signal = %.2f usec", mean(this_pw)*1000000);
    title(txt);
    o_name = sprintf("LED_%.0f_Pulse",mean(this_pw)*1000000);
    o_name = char(o_name);
    WriteFig(fig, o_name, 1);
    close all;
end

% Restricting LED to its only 800 msec and above epochs
led_pw = led_pw(4:end);
%% Laser_Analysis
laser_temp = cellfun(@(x) split(x), laser_labels, 'UniformOutput', false);
laser_pw = cellfun(@(x) str2double(x{2}), laser_temp);
laser_signal_width = [];
laser_start_delay = [];
laser_end_delay = [];
[laser_pw, laser_sort_idx] = sort(laser_pw);

laser_labels = laser_labels(laser_sort_idx);
for i = 4:length(laser_labels)
    start_idx = find(cellfun(@(x) contains(x,laser_labels{i}), evs.label) == 1);
    start_time = evs.t{start_idx};
    if i ~= length(laser_labels)
        end_idx = find(cellfun(@(x) contains(x,laser_labels{i+1}), evs.label) == 1);
        end_time = evs.t{end_idx};
    else
        end_time = evs.t{51}(1); % Because Laser recording was first
    end
    this_laser_ON_times = pulse_ON_times(pulse_ON_times > start_time & pulse_ON_times < end_time);
    this_laser_OFF_times= pulse_OFF_times(pulse_OFF_times > start_time & pulse_OFF_times < end_time);
    snippet_starts = nearest_idx3(this_laser_ON_times - 1e-3, photo_signal.tvec);
    plot(photo_signal.tvec(snippet_starts(1)-1000:snippet_starts(end)+1000), ...
        photo_signal.data(snippet_starts(1)-1000:snippet_starts(end)+1000));
    title(laser_labels{i})
    close;
end
%%
% start from 4 because < 800 msec looks weird
for i = 4:length(laser_labels)
    start_idx = find(cellfun(@(x) contains(x,laser_labels{i}), evs.label) == 1);
    start_time = evs.t{start_idx};
    if i ~= length(laser_labels)
        end_idx = find(cellfun(@(x) contains(x,laser_labels{i+1}), evs.label) == 1);
        end_time = evs.t{end_idx};
    else
        end_time = evs.t{51}(1); % Because Laser recording was first
    end
    this_laser_ON_times = pulse_ON_times(pulse_ON_times > start_time & pulse_ON_times < end_time);
    this_laser_OFF_times= pulse_OFF_times(pulse_OFF_times > start_time & pulse_OFF_times < end_time);
    snippet_starts = nearest_idx3(this_laser_ON_times - 1e-2, photo_signal.tvec);
    snippet_ends = nearest_idx3(this_laser_OFF_times + 1e-2, photo_signal.tvec);
    snip_width = max(snippet_ends - snippet_starts);
    this_pw = this_laser_OFF_times - this_laser_ON_times;
    pulse_starts = nearest_idx3(this_laser_ON_times, photo_signal.tvec);
    pulse_ends = nearest_idx3(this_laser_OFF_times, photo_signal.tvec);
    pulse_starts = pulse_starts - snippet_starts;
    pulse_ends = pulse_ends - snippet_starts;

    % Snippets are pretty variable, consider only snippets that reach some
    % thresold of the max deivation
    snip_max = zeros(size(snippet_starts));
    for iS = 1:length(snippet_starts)
        snip_max(iS) = max(photo_signal.data(snippet_starts(iS):snippet_starts(iS)+snip_width-1)+0.25);
    end
    snip_thresh = 0.999*max(snip_max);

    % Calculating the average signal because there isn't much variation
    average_snippet = zeros(1,snip_width);
    for j = 1:length(snippet_starts)
        if snip_max(j) >= snip_thresh
            average_snippet = average_snippet + ...
                photo_signal.data(snippet_starts(j):snippet_starts(j)+snip_width-1);
        end
    end
    average_snippet = (average_snippet/sum(snip_max>=snip_thresh))+0.25; %Adding 0.25 to shift the whole thing above 0
    
    % Defining the actual signal start as where the signal reaches 93% of 
    % its max value and the signal end as where the signal drops below the
    % 93% of its max value

    threshold = 0.93*max(average_snippet);
    actual_start = find(average_snippet>threshold, 1);
    actual_width = find(average_snippet(actual_start+1:end) < threshold,1)+1;
    actual_end = actual_start+actual_width;
    actual_start_delay = actual_start - mean(pulse_starts);
    actual_end_delay = actual_end - mean(pulse_ends);
    laser_signal_width = [laser_signal_width actual_width];
    laser_start_delay = [laser_start_delay actual_start_delay];
    laser_end_delay = [laser_end_delay actual_end_delay];
    
    fig = figure('WindowState', 'maximized');
    plot(average_snippet);
    temp = xticklabels;
    temp = cellfun(@(x) char(string(str2double(x)*1000000/32000)), temp, 'UniformOutput', 0);
    xticklabels(temp);
    xlabel("Time (microseconds)")
    ylabel("PhotoSensor Signal")
    vline(pulse_starts(1), 'r');
    vline(pulse_ends(1), '--r');
    vline(actual_start, 'g');
    vline(actual_end, '--g');
    
    txt = sprintf("Actual Width = %.2f usec", actual_width*1000000/32000);
    text(600, 0.215, txt, 'FontSize', 12);
    txt = sprintf("Start Delay = %.2f usec", actual_start_delay*1000000/32000);
    text(600, 0.21, txt, 'FontSize', 12);
    txt = sprintf("End Delay = %.2f usec", actual_end_delay*1000000/32000);
    text(600, 0.205, txt, 'FontSize', 12);
    txt = sprintf("Laser Pulse Signal = %.2f usec", mean(this_pw)*1000000);
    title(txt);
    o_name = sprintf("Laser_%.0f_Pulse",mean(this_pw)*1000000);
    o_name = char(o_name);
    WriteFig(fig, o_name, 1);
    close all;
end

% Restricting LASER to its only 800 msec and above epochs
laser_pw = laser_pw(4:end);
%% Plotting aggregate statistics
% Plot Laser aggregate
fig = figure('WindowState', 'maximized');
scatter(laser_pw, laser_signal_width*1000000/32000, 'filled');
hold on;
scatter(laser_pw, laser_start_delay*1000000/32000, 'filled');
hold on;
scatter(laser_pw, laser_end_delay*1000000/32000, 'filled');
hold on;
plot((0:3000),(0:3000));
legend("laser signal width", "laser start delay", "laser end delay", "Y=X");
xlabel("Laser Pulse Width (microsec)");
ylabel("Time (microsec)");
grid on;
grid minor;
WriteFig(fig, char("Laser_aggregate"), 1);
close all;
%%
% Plot LED aggregate
fig = figure('WindowState', 'maximized');
scatter(led_pw, led_signal_width*1000000/32000, 'filled');
hold on;
scatter(led_pw, led_start_delay*1000000/32000, 'filled');
hold on;
scatter(led_pw, led_end_delay*1000000/32000, 'filled');
hold on;
plot((0:3000),(0:3000));
legend("led signal width", "led start delay", "led end delay", "Y=X");
xlabel("LED Pulse Width (microsec)");
ylabel("Time (microsec)");
grid on;
grid minor;
WriteFig(fig, char("LED_aggregate"), 1);
close all;

%% Extended Laser Analysis
laser_temp = cellfun(@(x) split(x), laser_labels, 'UniformOutput', false);
laser_pw = cellfun(@(x) str2double(x{2}), laser_temp);
laser_signal_width = [];
laser_start_delay = [];
laser_end_delay = [];
[laser_pw, laser_sort_idx] = sort(laser_pw);
all_actual_starts = {};
all_actual_ends = {};
all_actual_widths = {};
all_actual_start_delays = {};
all_actual_end_delays = {};
all_mid_starts = {};
all_mid_ends = {};
all_mid_widths = {};
all_mid_start_delay_1 = {};
all_mid_start_delay_2 = {};
all_last_starts = {};
all_last_ends = {};
all_last_widths = {};
all_last_start_delay_1 = {};
all_last_start_delay_2 = {};

laser_labels = laser_labels(laser_sort_idx);
for i = 1:length(laser_labels)
    start_idx = find(cellfun(@(x) contains(x,laser_labels{i}), evs.label) == 1);
    start_time = evs.t{start_idx};
    if i ~= length(laser_labels)
        end_idx = find(cellfun(@(x) contains(x,laser_labels{i+1}), evs.label) == 1);
        end_time = evs.t{end_idx};
    else
        end_time = evs.t{45}(1); % Because Laser recording was first
    end
    
    this_laser_ON_times = pulse_ON_times(pulse_ON_times > start_time & pulse_ON_times < end_time);
    this_laser_OFF_times= pulse_OFF_times(pulse_OFF_times > start_time & pulse_OFF_times < end_time);
    snippet_starts = nearest_idx3(this_laser_ON_times - 1e-3, photo_signal.tvec);
    mid_snippet_starts = nearest_idx3(this_laser_ON_times + 6e-3, photo_signal.tvec);
    mid_snippet_ends = nearest_idx3(this_laser_ON_times + 12e-3, photo_signal.tvec);
    
    last_snippet_starts = nearest_idx3(this_laser_ON_times + 14e-3, photo_signal.tvec);
    last_snippet_ends = nearest_idx3(this_laser_ON_times + 19e-3, photo_signal.tvec);
    
    
    snippet_ends = nearest_idx3(this_laser_ON_times + 2e-2, photo_signal.tvec);
    this_pw = this_laser_OFF_times - this_laser_ON_times;
    pulse_starts = nearest_idx3(this_laser_ON_times, photo_signal.tvec);
    pulse_ends = nearest_idx3(this_laser_OFF_times, photo_signal.tvec);
    pulse_starts = pulse_starts - snippet_starts;
    pulse_ends = pulse_ends - snippet_starts;
    
    this_actual_starts = [];
    this_actual_ends = [];
    this_actual_widths = [];
    this_actual_start_delays = [];
    this_actual_end_delays = [];
    this_mid_starts = [];
    this_mid_ends = [];
    this_mid_widths = [];
    this_mid_start_delay_1 = [];
    this_mid_start_delay_2 = [];
    this_last_starts = [];
    this_last_ends = [];
    this_last_widths = [];
    this_last_start_delay_1 = [];
    this_last_start_delay_2 = [];
    
    for j = 1:length(snippet_starts)
        snippet = photo_signal.data(snippet_starts(j):snippet_ends(j));
        mid_snippet = photo_signal.data(mid_snippet_starts(j):mid_snippet_ends(j));
        last_snippet = photo_signal.data(last_snippet_starts(j):last_snippet_ends(j));
        threshold = 0.9*max(snippet);
        actual_start = find(snippet>threshold);
        actual_start = actual_start(1);
        actual_end = find(snippet(actual_start+1:end) < threshold);
        actual_end = actual_start + actual_end(1);
        actual_width = actual_end - actual_start;
        actual_start_delay = actual_start - mean(pulse_starts);
        actual_end_delay = actual_end - mean(pulse_ends);
        this_actual_widths = [this_actual_widths actual_width];
        this_actual_start_delays = [this_actual_start_delays actual_start_delay];
        this_actual_end_delays = [this_actual_end_delays actual_end_delay];
        this_actual_starts = [this_actual_starts actual_start];
        this_actual_ends = [this_actual_ends actual_end];
        
        [mid_peak, mid_idx] = max(mid_snippet);
        % hacky way to go around peaks below 0
        if max(mid_snippet) < 0
            mid_threshold = 1.1 * max(mid_snippet);
        else
            mid_threshold = 0.9 * max(mid_snippet);
        end
        mid_start = find(mid_snippet(1:mid_idx) < mid_threshold);
        if isempty(mid_start)
            mid_start = mid_peak;
        else
            mid_start = mid_start(end);
        end
        mid_end = find(mid_snippet(mid_idx:end) < mid_threshold);
        if isempty(mid_end)
            mid_end = mid_peak;
        else
            mid_end = mid_end(1) + mid_idx;
        end
        mid_start = mid_start + mid_snippet_starts(j) - snippet_starts(j);
        mid_end = mid_end + mid_snippet_starts(j) - snippet_starts(j);
        mid_width = mid_end - mid_start;
        mid_start_delay_1 = mid_start - pulse_ends(j);
        mid_start_delay_2 = mid_start - actual_end;
        this_mid_starts = [this_mid_starts mid_start];
        this_mid_ends = [this_mid_ends mid_end];
        this_mid_widths = [this_mid_widths mid_width];
        this_mid_start_delay_1 = [this_mid_start_delay_1 mid_start_delay_1];
        this_mid_start_delay_2 = [this_mid_start_delay_2 mid_start_delay_2];
        
        % hacky way to go around peaks below 0
        [last_peak, last_idx] = max(last_snippet);
        if max(last_snippet) < 0
            last_threshold = 1.1 * max(last_snippet);
        else
            last_threshold = 0.9 * max(last_snippet);
        end
        last_start = find(last_snippet(1:last_idx) < last_threshold);
        if isempty(last_start)
            last_start = last_idx;
        else
            last_start = last_start(end);
        end
        last_end = find(last_snippet(last_idx:end) < last_threshold);
        if isempty(last_end)
            last_end = last_idx;
        else
            last_end = last_end(1) + last_idx;
        end
        last_start = last_start + last_snippet_starts(j) - snippet_starts(j);
        last_end = last_end + last_snippet_starts(j) - snippet_starts(j);
        last_width = last_end - last_start;
        last_start_delay_1 = last_start - pulse_ends(j);
        last_start_delay_2 = last_start - actual_end;
        this_last_starts = [this_last_starts last_start];
        this_last_ends = [this_last_ends last_end];
        this_last_widths = [this_last_widths last_width];
        this_last_start_delay_1 = [this_last_start_delay_1 last_start_delay_1];
        this_last_start_delay_2 = [this_last_start_delay_2 last_start_delay_2];
        
        % Uncomment to generate figures for each iteration
%         figure('WindowState', 'maximized');
%         plot(snippet);
%         vline(mid_start);
%         vline(mid_end);
%         vline(last_start);
%         vline(last_end);
%         txt = sprintf("Mid Width = %.2f usec", mid_width*1000000/32000);
%         text(600, 0.05, txt, 'FontSize', 12);
%         txt = sprintf("Mid Start Delay 1 = %.2f usec", mid_start_delay_1*1000000/32000);
%         text(600, 0.045, txt, 'FontSize', 12);
%         txt = sprintf("Mid Start Delay 2 = %.2f usec", mid_start_delay_2*1000000/32000);
%         text(600, 0.04, txt, 'FontSize', 12);
%         txt = sprintf("Last Width = %.2f usec", last_width*1000000/32000);
%         text(600, 0.035, txt, 'FontSize', 12);
%         txt = sprintf("Last Start Delay 1 = %.2f usec", last_start_delay_1*1000000/32000);
%         text(600, 0.03, txt, 'FontSize', 12);
%         txt = sprintf("Last Start Delay 2 = %.2f usec", last_start_delay_2*1000000/32000);
%         text(600, 0.025, txt, 'FontSize', 12);
%         close;
        
    end
    
    all_actual_starts{i} = this_actual_starts;
    all_actual_ends{i} = this_actual_ends;
    all_actual_widths{i} = this_actual_widths;
    all_actual_start_delays{i} = this_actual_start_delays;
    all_actual_end_delays{i} = this_actual_end_delays;
    all_mid_starts{i} = this_mid_starts;
    all_mid_ends{i} = this_mid_ends;
    all_mid_widths{i} = this_mid_widths;
    all_mid_start_delay_1{i} = this_mid_start_delay_1;
    all_mid_start_delay_2{i} = this_mid_start_delay_2;
    all_last_starts{i} = this_last_starts;
    all_last_ends{i} = this_last_ends;
    all_last_widths{i} = this_last_widths;
    all_last_start_delay_1{i} = this_last_start_delay_1;
    all_last_start_delay_2{i} = this_last_start_delay_2;    
end

%% Plot aggregate with error bars
fig = figure('WindowState', 'maximized');
all_pw = [500:100:2500];
mean_mid_widths = cellfun(@(x) mean(x)*1000000/32000, all_mid_widths, 'UniformOutput', 1);
error_mid_widths = cellfun(@(x) std(x)*1000000/32000, all_mid_widths, 'UniformOutput', 1);
errorbar(all_pw, mean_mid_widths, error_mid_widths);
hold on;
mean_mid_start_delay_1 = cellfun(@(x) mean(x)*1000000/32000, all_mid_start_delay_1, 'UniformOutput', 1);
error_mid_start_delay_1 = cellfun(@(x) std(x)*1000000/32000, all_mid_start_delay_1, 'UniformOutput', 1);
errorbar(all_pw, mean_mid_start_delay_1, error_mid_start_delay_1);
hold on;
mean_mid_start_delay_2 = cellfun(@(x) mean(x)*1000000/32000, all_mid_start_delay_2, 'UniformOutput', 1);
error_mid_start_delay_2 = cellfun(@(x) std(x)*1000000/32000, all_mid_start_delay_2, 'UniformOutput', 1);
errorbar(all_pw, mean_mid_start_delay_2, error_mid_start_delay_2);
xlim([0,3000]);
legend("Mid Bump Width", "Mid Bump onset delay from Pulse OFF", "Mid Bump onset delay from Signal OFF");
grid on;
grid minor;
ylabel("Time (microsec)");
xlabel("Laser Pulse Width (microsec)");
WriteFig(fig, char( "Mid_bump_aggregate"), 1);
close all;

fig = figure('WindowState', 'maximized');
all_pw = [500:100:2500];
mean_last_widths = cellfun(@(x) mean(x)*1000000/32000, all_last_widths, 'UniformOutput', 1);
error_last_widths = cellfun(@(x) std(x)*1000000/32000, all_last_widths, 'UniformOutput', 1);
errorbar(all_pw, mean_last_widths, error_last_widths);
hold on;
mean_last_start_delay_1 = cellfun(@(x) mean(x)*1000000/32000, all_last_start_delay_1, 'UniformOutput', 1);
error_last_start_delay_1 = cellfun(@(x) std(x)*1000000/32000, all_last_start_delay_1, 'UniformOutput', 1);
errorbar(all_pw, mean_last_start_delay_1, error_last_start_delay_1);
hold on;
mean_last_start_delay_2 = cellfun(@(x) mean(x)*1000000/32000, all_last_start_delay_2, 'UniformOutput', 1);
error_last_start_delay_2 = cellfun(@(x) std(x)*1000000/32000, all_last_start_delay_2, 'UniformOutput', 1);
errorbar(all_pw, mean_last_start_delay_2, error_last_start_delay_2);
xlim([0,3000]);
legend("Last Bump Width", "Last Bump onset delay from Pulse OFF", "Last Bump onset delay from Signal OFF");
grid on;
grid minor;
ylabel("Time (microsec)");
xlabel("Laser Pulse Width (microsec)");
WriteFig(fig, char( "Last_bump_aggregate"), 1);
close all;

    