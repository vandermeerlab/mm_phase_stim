%% Load the data from .evs and .ncs files
clear;
close all;
% cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/WheelEncoderTests/WE-2021-05-10')
cd('D:\Dropbox (Dartmouth College)\manish_data\WheelEncoderTests\WE-2021-05-10')
cfg.fc = {'WE1.ncs', 'WE2.ncs'};
we_csc = LoadCSC(cfg);


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
% Restricting the data from 0 to 11.5 sec
this_csc = restrict(we_csc, we_csc.tvec(1), we_csc.tvec(1) + 11.5);
we1_csc_up = find(diff(this_csc.data(1,:) > 0) > 0);
we1_csc_down = find(diff(this_csc.data(1,:) > 0) < 0);

we2_csc_up = find(diff(this_csc.data(2,:) > 0) > 0);
we2_csc_down = find(diff(this_csc.data(2,:) > 0) < 0);

%% Find up and down transition times from TTL
% fn = '/Users/manishm/Dropbox (Dartmouth College)/manish_data/WheelEncoderTests/WE-2021-05-10/Events.nev';
fn = 'D:\Dropbox (Dartmouth College)\manish_data\WheelEncoderTests\WE-2021-05-10\Events.nev';

[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);
keep1 = (EventIDs == 11);
keep2 = (EVTimeStamps <= this_csc.tvec(end)*1e6);
keep = keep1 & keep2;
ts = EVTimeStamps(keep);
state = TTLs(keep);

we1_TTL_up = zeros(1, length(state)-1);
we1_TTL_down = zeros(1, length(state)-1);
we2_TTL_up = zeros(1, length(state)-1);
we2_TTL_down = zeros(1, length(state)-1);
for i = 1:length(state)-1
    prev_state = state(i);
    cur_state = state(i+1);
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

%%
we1_TTL_up_ts = ts(find(we1_TTL_up)+ 1);
we1_TTL_up_idx = nearest_idx3(we1_TTL_up_ts, this_csc.tvec*1e6);

we1_TTL_down_ts = ts(find(we1_TTL_down)+ 1);
we1_TTL_down_idx = nearest_idx3(we1_TTL_down_ts, this_csc.tvec*1e6);

we2_TTL_up_ts = ts(find(we2_TTL_up)+ 1);
we2_TTL_up_idx = nearest_idx3(we2_TTL_up_ts, this_csc.tvec*1e6);

we2_TTL_down_ts = ts(find(we2_TTL_down)+ 1);
we2_TTL_down_idx = nearest_idx3(we2_TTL_down_ts, this_csc.tvec*1e6);


%%


%State transitions are probably exchanged for WE1 and WE2


figure
plot(this_csc.data(1,:))
hold on
vline([0,10])
figure
plot(this_csc.data(1,:))
hold on
vline(we1_csc_up, 'green')
hold on
vline(we1_TTL_up_idx, 'red')