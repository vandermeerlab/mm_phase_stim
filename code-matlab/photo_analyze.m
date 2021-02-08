clear;
close all;
cd('/Users/manishm/Work/vanDerMeerLab/Data/PhotoSensor_Tests/PhotoSensor-2021-02-05/')
evs = LoadEvents([]);
cfg.fc = {'PhotoSensor.ncs'};
photo_signal = LoadCSC(cfg);
laser_idx = cellfun(@(x) contains(x, 'Laser'), evs.label, 'UniformOutput', 1);
led_idx = cellfun(@(x) contains(x, 'LED'), evs.label, 'UniformOutput', 1);
laser_labels = evs.label(laser_idx);
led_labels = evs.label(led_idx);

pulse_OFF_string = evs.label{46};
pulse_ON_string = evs.label{47};

pulse_OFF_times = evs.t{46};
pulse_ON_times = evs.t{47};



%% Laser_Analysis
laser_temp = cellfun(@(x) split(x), laser_labels, 'UniformOutput', false);
laser_pw = cellfun(@(x) str2double(x{2}), laser_temp);
laser_signal_width = [];
laser_start_delay = [];
laser_end_delay = [];
[laser_pw, laser_sort_idx] = sort(laser_pw);

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
    snippet_ends = nearest_idx3(this_laser_ON_times + 2e-2, photo_signal.tvec);
    this_pw = this_laser_OFF_times - this_laser_ON_times;
    pulse_starts = nearest_idx3(this_laser_ON_times, photo_signal.tvec);
    pulse_ends = nearest_idx3(this_laser_OFF_times, photo_signal.tvec);
    pulse_starts = pulse_starts - snippet_starts;
    pulse_ends = pulse_ends - snippet_starts;
    
    % Calculating the average signal because there isn't much variation
    average_snippet = zeros(1,length(photo_signal.data(snippet_starts(1):snippet_ends(1))));
    for i = 2:length(snippet_starts)
        average_snippet = average_snippet + ...
            photo_signal.data(snippet_starts(i):snippet_ends(i));
    end
    average_snippet = average_snippet/(length(snippet_ends));
    
    % Defining the actual signal start as where the signal reaches 90% of 
    % its max value and the signal end as where the signal drops below the
    % 90% of its max value
    threshold = 0.9*max(average_snippet);
    actual_start = find(average_snippet>threshold);
    actual_start = actual_start(1);
    actual_end = find(average_snippet(actual_start+1:end) < threshold);
    actual_end = actual_start + actual_end(1);
%     actual_end = actual_end(1);
    actual_width = actual_end - actual_start;
    actual_start_delay = actual_start - mean(pulse_starts);
    actual_end_delay = actual_end - mean(pulse_ends);
    laser_signal_width = [laser_signal_width actual_width];
    laser_start_delay = [laser_start_delay actual_start_delay];
    laser_end_delay = [laser_end_delay actual_end_delay];
    
    fig = figure('WindowState', 'maximized');
    plot(average_snippet);
    vline(pulse_starts(1), 'r');
    vline(pulse_ends(1), 'r');
    vline(actual_start, 'g');
    vline(actual_end, 'g');
    
    txt = sprintf("Actual Width = %.2f usec", actual_width*1000000/32000);
    text(600, 0.05, txt, 'FontSize', 12);
    txt = sprintf("Start Delay = %.2f usec", actual_start_delay*1000000/32000);
    text(600, 0.045, txt, 'FontSize', 12);
    txt = sprintf("End Delay = %.2f usec", actual_end_delay*1000000/32000);
    text(600, 0.04, txt, 'FontSize', 12);
    txt = sprintf("Laser Pulse Signal = %.2f usec", mean(this_pw)*1000000);
    title(txt);
    o_name = sprintf("Laser_%.0f_Pulse",mean(this_pw)*1000000);
    o_name = char(o_name);
    WriteFig(fig, o_name, 1);
    close all;
end

%% LED Analysis
led_temp = cellfun(@(x) split(x), led_labels, 'UniformOutput', false);
led_pw = cellfun(@(x) str2double(x{2}), led_temp);
led_signal_width = [];
led_start_delay = [];
led_end_delay = [];
[led_pw, led_sort_idx] = sort(led_pw);

led_labels = led_labels(led_sort_idx);
for i = 1:length(led_labels)
    start_idx = find(cellfun(@(x) contains(x,led_labels{i}), evs.label) == 1);
    start_time = evs.t{start_idx};
    if i ~= length(led_labels)
        end_idx = find(cellfun(@(x) contains(x,led_labels{i+1}), evs.label) == 1);
        end_time = evs.t{end_idx};
    else
        end_time = evs.t{45}(2); % Because LED recording was second
    end
    this_led_ON_times = pulse_ON_times(pulse_ON_times > start_time & pulse_ON_times < end_time);
    this_led_OFF_times= pulse_OFF_times(pulse_OFF_times > start_time & pulse_OFF_times < end_time);
    snippet_starts = nearest_idx3(this_led_ON_times - 1e-3, photo_signal.tvec);
    snippet_ends = nearest_idx3(this_led_ON_times + 2e-2, photo_signal.tvec);
    this_pw = this_led_OFF_times - this_led_ON_times;
    pulse_starts = nearest_idx3(this_led_ON_times, photo_signal.tvec);
    pulse_ends = nearest_idx3(this_led_OFF_times, photo_signal.tvec);
    pulse_starts = pulse_starts - snippet_starts;
    pulse_ends = pulse_ends - snippet_starts;
    
    % Calculating the average signal because there isn't much variation
    average_snippet = zeros(1,length(photo_signal.data(snippet_starts(1):snippet_ends(1))));
    for i = 2:length(snippet_starts)
        average_snippet = average_snippet + ...
            photo_signal.data(snippet_starts(i):snippet_ends(i));
    end
    average_snippet = average_snippet/(length(snippet_ends));
    
    % Defining the actual signal start as where the signal reaches 90% of 
    % its max value and the signal end as where the signal drops below the
    % 90% of its max value
    threshold = 0.9*max(average_snippet);
    actual_start = find(average_snippet>threshold);
    actual_start = actual_start(1);
    actual_end = find(average_snippet(actual_start+1:end) < threshold);
    actual_end = actual_start + actual_end(1);
%     actual_end = actual_end(1);
    actual_width = actual_end - actual_start;
    actual_start_delay = actual_start - mean(pulse_starts);
    actual_end_delay = actual_end - mean(pulse_ends);
    led_signal_width = [led_signal_width actual_width];
    led_start_delay = [led_start_delay actual_start_delay];
    led_end_delay = [led_end_delay actual_end_delay];
    
    fig = figure('WindowState', 'maximized');
    plot(average_snippet);
    vline(pulse_starts(1), 'r');
    vline(pulse_ends(1), 'r');
    vline(actual_start, 'g');
    vline(actual_end, 'g');
    
    txt = sprintf("Actual Width = %.2f usec", actual_width*1000000/32000);
    text(600, 0.05, txt, 'FontSize', 12);
    txt = sprintf("Start Delay = %.2f usec", actual_start_delay*1000000/32000);
    text(600, 0.045, txt, 'FontSize', 12);
    txt = sprintf("End Delay = %.2f usec", actual_end_delay*1000000/32000);
    text(600, 0.04, txt, 'FontSize', 12);
    txt = sprintf("LED Pulse Signal = %.2f usec", mean(this_pw)*1000000);
    title(txt);
    o_name = sprintf("LED_%.0f_Pulse",mean(this_pw)*1000000);
    o_name = char(o_name);
    WriteFig(fig, o_name, 1);
    close all;
end

