%% Load the data from .evs and .ncs files
clear;
close all;
cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/WheelEncoderTests/WE-2021-05-10')
% cd('D:\Dropbox (Dartmouth College)\manish_data\WheelEncoderTests\WE-2021-05-10')
cfg.fc = {'WE1.ncs', 'WE2.ncs'};
we_csc = LoadCSC(cfg);

% Changing the sign of the data because input was inverted in the config
we_csc.data = - we_csc.data;


%% Plot Raw CSC Files to look at appropriate epochs
figure;
ax1 = subplot(2,1,1);
plot(we_csc.tvec - we_csc.tvec(1), we_csc.data(1,:))
title('WE1.ncs')
ax2 = subplot(2,1,2);
plot(we_csc.tvec - we_csc.tvec(1), we_csc.data(2,:), 'red')
title('WE2.ncs')
linkaxes([ax1 ax2],'x');

%% Find up and down transition times from raw_csc
% 'Good Forward': 0 to 11.5 sec 
% 'Good Backward' 28 to 34.2 sec
% 'Mixed' : 43.5 to 48.5 sec
this_csc = restrict(we_csc, we_csc.tvec(1)+28, we_csc.tvec(1) + 34.2);
we1_csc_up = find(diff(this_csc.data(1,:) > 0) > 0);
we1_csc_down = find(diff(this_csc.data(1,:) > 0) < 0);

we2_csc_up = find(diff(this_csc.data(2,:) > 0) > 0);
we2_csc_down = find(diff(this_csc.data(2,:) > 0) < 0);

%% Find up and down transition times from TTL
fn = '/Users/manishm/Dropbox (Dartmouth College)/manish_data/WheelEncoderTests/WE-2021-05-10/Events.nev';
% fn = 'D:\Dropbox (Dartmouth College)\manish_data\WheelEncoderTests\WE-2021-05-10\Events.nev';

[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);
keep1 = (EventIDs == 11);
keep2 = (EVTimeStamps <= this_csc.tvec(end)*1e6  & EVTimeStamps >= this_csc.tvec(1)*1e6);
keep = keep1 & keep2;
ts = EVTimeStamps(keep);
state = TTLs(keep);

% The first state will never be '0' because the TTL pins register a value 
% only when at least one pin is set, and after that every change.
we1_TTL_up = zeros(1, length(state));
we1_TTL_down = zeros(1, length(state));
we2_TTL_up = zeros(1, length(state));
we2_TTL_down = zeros(1, length(state));
for i = 1:length(state)
    if i == 1
        if state(i) == 16
            we2_TTL_up(i) = 1;
        elseif state(i) == 32
            we1_TTL_up(i) = 1;
        elseif state(i) == 48 % Unlikely case
            we1_TTL_up(i) = 1;
            we2_TTL_up(i) = 1;
        end
    else
        prev_state = state(i-1);
        cur_state = state(i);
        if prev_state == 0 && cur_state == 32
            we1_TTL_up(i) = 1;
        elseif prev_state == 16 && cur_state == 48
            we1_TTL_up(i) = 1;
        elseif prev_state == 16 && cur_state == 32 % Unlikely case
            we1_TTL_up(i) = 1;
            we2_TTL_down(i) = 1;
        elseif prev_state == 32 && cur_state == 0
            we1_TTL_down(i) = 1;
        elseif prev_state == 48 && cur_state == 16
            we1_TTL_down(i) = 1;
        elseif prev_state == 32 && cur_state == 16 % Unlikely case
            we1_TTL_down(i) = 1;
            we2_TTL_up(i) = 1;
        elseif prev_state == 0 && cur_state == 16
            we2_TTL_up(i) = 1;
        elseif prev_state == 32 && cur_state == 48
            we2_TTL_up(i) = 1;
        elseif prev_state == 16 && cur_state == 0
            we2_TTL_down(i) = 1;
        elseif prev_state == 48 && cur_state == 32
            we2_TTL_down(i) = 1;
        elseif prev_state == 48 && cur_state == 0 % Unlikely case
            we1_TTL_down(i) = 1;
            we2_TTL_down(i) = 1;
        elseif prev_state == 0 && cur_state == 48 % Unlikely case
            we1_TTL_up(i) = 1;
            we2_TTL_up(i) = 1;
        end
    end
    
end
total_we1_up = sum(we1_TTL_up);
total_we1_down = sum(we1_TTL_down);
total_we2_up = sum(we2_TTL_up);
total_we2_down = sum(we2_TTL_down);

s1 = sprintf("WE1 CSC UP: %d, WE1 TTL UP: %d", length(we1_csc_up), total_we1_up);
disp(s1)
s2 = sprintf("WE1 CSC DOWN: %d, WE1 TTL DOWN: %d", length(we1_csc_down), total_we1_down);
disp(s2)
s3 = sprintf("WE2 CSC UP: %d, WE2 TTL UP: %d", length(we2_csc_up), total_we2_up);
disp(s3)
s4 = sprintf("WE2 CSC DOWN: %d, WE2 TTL DOWN: %d", length(we2_csc_up), total_we2_up);
disp(s4)

%% Find indices of the TTL-derived transitions on the CSC timebase
we1_TTL_up_ts = ts(find(we1_TTL_up));
we1_TTL_up_idx = nearest_idx3(we1_TTL_up_ts, this_csc.tvec*1e6);

we1_TTL_down_ts = ts(find(we1_TTL_down));
we1_TTL_down_idx = nearest_idx3(we1_TTL_down_ts, this_csc.tvec*1e6);

we2_TTL_up_ts = ts(find(we2_TTL_up));
we2_TTL_up_idx = nearest_idx3(we2_TTL_up_ts, this_csc.tvec*1e6);

we2_TTL_down_ts = ts(find(we2_TTL_down));
we2_TTL_down_idx = nearest_idx3(we2_TTL_down_ts, this_csc.tvec*1e6);


%% Plot
figure
ax1 = subplot(2,1,1);
plot(this_csc.data(1,:));
hold on
% Hacky way to plot multiple xlines
for i = 1:length(we1_csc_up)
     xline(we1_csc_up(i), 'Color', 'g');
end
for i = 1:length(we1_TTL_up_idx)
     xline(we1_TTL_up_idx(i), 'Color', 'r');
end
for i = 1:length(we1_csc_down)
     xline(we1_csc_down(i), '--', 'Color', 'g');
end
for i = 1:length(we1_TTL_down_idx)
     xline(we1_TTL_down_idx(i), '--', 'Color', 'r');
end
title('WE1.ncs')

% 
% arrayfun(@(x) xline(x, 'Color', 'g'), we1_csc_up, 'UniformOutput', false);
% arrayfun(@(x) xline(x, 'Color', 'r'), we1_TTL_up_idx, 'UniformOutput', false);
% arrayfun(@(x) xline(x, '--','Color', 'g'), we1_csc_down, 'UniformOutput', false);
% arrayfun(@(x) xline(x, '--','Color', 'r'), we1_TTL_down_idx, 'UniformOutput', false);


ax2 = subplot(2,1,2);
plot(this_csc.data(2,:));
hold on
% Hacky way to plot multiple xlines
for i = 1:length(we2_csc_up)
     xline(we2_csc_up(i), 'Color', 'g');
end
for i = 1:length(we2_TTL_up_idx)
     xline(we2_TTL_up_idx(i), 'Color', 'r');
end
for i = 1:length(we2_csc_down)
     xline(we2_csc_down(i), '--', 'Color', 'g');
end
for i = 1:length(we2_TTL_down_idx)
     xline(we2_TTL_down_idx(i), '--', 'Color', 'r');
end

% arrayfun(@(x) xline(x, 'Color', 'g'), we2_csc_up, 'UniformOutput', false);
% arrayfun(@(x) xline(x, 'Color', 'r'), we2_TTL_up_idx, 'UniformOutput', false);
% arrayfun(@(x) xline(x, '--','Color', 'g'), we2_csc_down, 'UniformOutput', false);
% arrayfun(@(x) xline(x, '--','Color', 'r'), we2_TTL_down_idx, 'UniformOutput', false);
title('WE2.ncs')

linkaxes([ax1 ax2], 'xy')


