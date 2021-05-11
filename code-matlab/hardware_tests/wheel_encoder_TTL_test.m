%% Load the data from .evs and .ncs files
clear;
close all;
cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/WheelEncoderTests/WE-2021-05-10')
cfg.fc = {'WE1.ncs', 'WE2.ncs'};
we_csc = LoadCSC(cfg);
fn = '/Users/manishm/Dropbox (Dartmouth College)/manish_data/WheelEncoderTests/WE-2021-05-10/Events.nev';


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
up_thresh = 0.03;
down_thresh = -0.03;
some_vec = (this_csc.data(1,:) > 0);
some_dif = diff(some_vec);
some_up = find(some_dif >0);
some_down = find(some_dif < 0);

%% Find up and down transition times from TTL

[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);
keep1 = (EventIDs == 11);
keep2 = (EVTimeStamps <= this_csc.tvec(end)*1e6);
keep = keep1 & keep2;
ts = EVTimeStamps(keep);
state = TTLs(keep);

ON_state = zeros(1, length(state)-1);
OFF_state = zeros(1, length(state)-1);
for i = 2:length(state)
    prev_state = state(i-1);
    cur_state = state(i);
    if prev_state == 0 && cur_state == 32
        ON_state(i) = 1;
    elseif prev_state == 16 && cur_state == 48
        ON_state(i) = 1;
    elseif prev_state == 16 && cur_state == 32 % Unlikely case
        ON_state(i) = 1;
    elseif prev_state == 32 && cur_state == 0
        OFF_state(i) = 1;
    elseif prev_state == 48 && cur_state == 16
        OFF_state(i) = 1;
    elseif prev_state == 32 && cur_state == 16 % Unlikely case
        OFF_state(i) = 1;
    end
end
total_on = sum(ON_state);
total_off = sum(OFF_state);

%%
laser_idx = cellfun(@(x) contains(x, 'Laser'), evs.label, 'UniformOutput', 1);
led_idx = cellfun(@(x) contains(x, 'LED'), evs.label, 'UniformOutput', 1);
laser_labels = evs.label(laser_idx);
led_labels = evs.label(led_idx);

pulse_OFF_string = evs.label{46};
pulse_ON_string = evs.label{47};

pulse_OFF_times = evs.t{46};
pulse_ON_times = evs.t{47};