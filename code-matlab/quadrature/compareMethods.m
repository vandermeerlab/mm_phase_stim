%% Script to compare quadrature encoding via TTL and  Wheel Data
cd('E:\Dropbox (Dartmouth College)\manish_data\WheelEncoderTests\WE-2021-05-10\')
evs = LoadEvents([]);
%%
cfg = [];
cfg.chA='WE1.ncs';
cfg.chB='WE2.ncs';
updownTSD = getQEupdown(cfg);
state_tsd = ConvertQEUpDownToState(updownTSD);
[angle_tsd, wheel_tsd, bad_jumps] = ConvertQEStatesToAngle([], state_tsd);
% cfg = []; cfg.downsample = 1;
[d, speed, cfg] = ConvertWheeltoSpeed(cfg, wheel_tsd);
%% Script to test if quadrature encoding works on TTL data
  
% RR2 config
% str_AB = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0003).'; % both up
% str_A = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0002).'; % only chA up
% str_B = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0001).'; % only chB up
% str_d = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000).'; % both down 

% Surgery config
str_AB = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0030).'; % both up
str_A = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0020).'; % only chA up
str_B = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0010).'; % only chB up
str_d = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0000).'; % both down 

both_up = evs.t{strcmp(evs.label, str_AB)};
chA_up = evs.t{strcmp(evs.label, str_A)};
chB_up = evs.t{strcmp(evs.label, str_B)};
both_down = evs.t{strcmp(evs.label, str_d)};

all_events = [both_up; chA_up; chB_up; both_down]';

sparse_states = [ones(size(both_up)); 2*ones(size(chA_up)); ...
    3*ones(size(chB_up)); 4*ones(size(both_down))]';

[all_events, sort_idx] = sort(all_events);

sparse_states = sparse_states(sort_idx);

% Use existing timebase
state_tsd2 = state_tsd;
state_tsd2.data = zeros(size(state_tsd2.data));
cidx = nearest_idx3(all_events, state_tsd2.tvec);
state_tsd2.data(1:cidx(1)) = sparse_states(1);
state_tsd2.data(cidx(end):end) = sparse_states(end);
for i = 1:length(cidx)-1
    state_tsd2.data(cidx(i):cidx(i+1)-1) = sparse_states(i);
end 
assert(isempty(find(state_tsd2.data == 0)));
%%
% Create new timebase
% dt = min(diff(all_events));
% time = all_events(1):dt:all_events(end);
% 
% if(all_events(end) < time(end))
%     time = [time time(end)+dt];
% end
% 
% 
% data = zeros(size(time));
% state_idx = nearest_idx3(all_events, time);
% assert((state_idx(1) == 1) && state_idx(end)==length(data));
% for i = 1:length(state_idx)-1
%     data(state_idx(i):state_idx(i+1)-1) = sparse_states(i); 
% end    
% data(end) = sparse_states(end);
% assert(isempty(find(data == 0)));

%% Time taking step
[angle_tsd2, wheel_tsd2, bad_jumps2] = ConvertQEStatesToAngle([], state_tsd2);
cfg2 = []; cfg2.downsample = 13; % The above quantity is 16
[d2, speed2, cfg2] = ConvertWheeltoSpeed(cfg2, wheel_tsd2);
%% RR2 test plot
figure;
plot(speed2.tvec, speed2.data);
hold on
xline(evs.t{strcmp(evs.label, 'one rotation forward')}, '--g');
xline(evs.t{strcmp(evs.label, 'one rotation backward')}, '--r');
xline(evs.t{strcmp(evs.label, 'Backward motion fast')}, '--r');
xline(evs.t{strcmp(evs.label, 'Forward motion fast')}, '--g');
ax = gca;
ax.XAxis.TickDirection = 'out';
legend({'Running Speed', 'Forward', 'Backward'})


%% Surgery Continuous Sampling Plot
figure;
subplot(2,1,1)
hold on
plot(speed.tvec, speed.data)
xline(evs.t{strcmp(evs.label, 'Forward')}, '--g');
xline(evs.t{strcmp(evs.label, 'Backward')}, '--r');
xline(evs.t{strcmp(evs.label, 'Mixed')}, '--black');
legend({'Running Speed', 'Forward', 'Backward', 'Mixed'})
title('DC-Coupled Channel')
%% Surgery Digital IO Plot
subplot(2,1,2)
hold on
plot(speed2.tvec, speed2.data)
xline(evs.t{strcmp(evs.label, 'Forward')}, '--g');
xline(evs.t{strcmp(evs.label, 'Backward')}, '--r');
xline(evs.t{strcmp(evs.label, 'Mixed')}, '--black');
legend({'Running Speed', 'Forward', 'Backward', 'Mixed'})
title('Digital IO')
