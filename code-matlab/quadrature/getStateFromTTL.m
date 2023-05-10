%% Script to test if quadrature encoding works on TTL data

cd('E:\Dropbox (Dartmouth College)\manish_data\M322\M322-2022-07-22');

evs = LoadEvents([]);

str_AB = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0003).'; % both up
str_A = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0002).'; % only chA up
str_B = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0001).'; % only chB up
str_d = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000).'; % both down 

both_up = evs.t{strcmp(evs.label, str_AB)};
chA_up = evs.t{strcmp(evs.label, str_A)};
chB_up = evs.t{strcmp(evs.label, str_B)};
both_down = evs.t{strcmp(evs.label, str_d)};

all_events = [both_up; chA_up; chB_up; both_down]';

sparse_states = [ones(size(both_up)); 2*ones(size(chA_up)); ...
    3*ones(size(chB_up)); 4*ones(size(both_down))]';

[all_events, sort_idx] = sort(all_events);

sparse_states = sparse_states(sort_idx);

dt = min(diff(all_events));
time = all_events(1):dt:all_events(end);

if(all_events(end) < time(end))
    time = [time time(end)+dt];
end


data = zeros(size(time));
state_idx = nearest_idx3(all_events, time);
assert((state_idx(1) == 1) && state_idx(end)==length(data));
for i = 1:length(state_idx)-1
    data(state_idx(i):state_idx(i+1)-1) = sparse_states(i); 
end    
data(end) = sparse_states(end);
assert(isempty(find(data == 0)));
state_tsd = tsd(time,data);
%%



