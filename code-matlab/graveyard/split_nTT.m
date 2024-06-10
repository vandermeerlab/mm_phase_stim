cd('E:\Dropbox (Dartmouth College)\manish_data\M348\M348-2022-05-14\');
evs = LoadEvents([]);
clean_start = evs.t{10}(3); % Timestamp associated with the last clean Timestamp

%%
for i = 1:8
    old_file = strcat('oldTT',num2str(i),'.ntt');
    new_file = strcat('TT',num2str(i),'.ntt');
    [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = ...
    Nlx2MatSpike( old_file, [1 1 1 1 1], 1, 1, []);
    keep = (Timestamps * 1e-6) > clean_start;
    Timestamps_r = Timestamps(keep);
    ScNumbers_r = ScNumbers(keep);
    CellNumbers_r = CellNumbers(keep);
    Features_r = Features(:,keep);
    Samples_r = Samples(:,:,keep);
    Mat2NlxSpike(new_file, 0, 1, [], ...
    [1 1 1 1 1 1] , Timestamps_r, ScNumbers_r, CellNumbers_r,...
                  Features_r, Samples_r, Header);
end

%Clean start is 
%%
% The delay between in the time stamps of when the ttl reaches neuralynx
% (evs.t{9}) from Master-8 and after routing through cyclops + op-amp
% (evs.t{14}) is in the order of 10s of microseconds. Thus we can ignore
%  evs.t{14} and evs.t{15}

% All Laser off_events
laser_off = evs.t{11}; % evs.t{9} for M078-11-26 and M075-2021-26, evs.t{11} fot M074

%Extract all records from the Tetorde File:
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = ...
    Nlx2MatSpike( 'TT8.ntt', [1 1 1 1 1], 1, 1, []);

% Works, Uncomment to run
% % Writing back to a dummy ntt
% Mat2NlxSpike( 'new_TT05.ntt', 0, 1, [], ...
%     [1 1 1 1 1 1] , Timestamps, ScNumbers, CellNumbers,...
%                   Features, Samples, Header);

%%
% Get all Pulse_ON_events
all_on = zeros(1,length(laser_off));
k = 1;
for i = 12:15 %10:13 for M078-11-26 and M075-2021-01-26, 12:15 for M074
    for j = 1:length(evs.t{i})
        all_on(k) = evs.t{i}(j);
        k = k+1;
    end
end

% Find Timestamps within the first 2 msec (1 for the window, 1 for the
% delay)

% Converting to seconds from microseconds
temp_ts = Timestamps * 1e-6; %

to_be_removed = zeros(1,length(temp_ts));

% For M074 (Laser) and M075 (Laser), remove spikes between ON + 1 ms to ON + 2 ms
% For M078 (LED), remove spikes between ON to ON + 1 ms
for i = 1:length(all_on)
    this_tbr = find(temp_ts > all_on(i) + 1e-3 & temp_ts < all_on(i)+ 2e-3);
    to_be_removed(this_tbr) = 1;
end
%% Remove the to_be_removed stuff
Timestamps_r = Timestamps(~to_be_removed);
ScNumbers_r = ScNumbers(~to_be_removed);
CellNumbers_r = CellNumbers(~to_be_removed);
Features_r = Features(:,~to_be_removed);
Samples_r = Samples(:,:,~to_be_removed);

Timestamps_w = Timestamps(~(~to_be_removed));
ScNumbers_w = ScNumbers(~(~to_be_removed));
CellNumbers_w = CellNumbers(~(~to_be_removed));
Features_w = Features(:,~(~to_be_removed));
Samples_w = Samples(:,:,~(~to_be_removed));
%%
 Mat2NlxSpike( 'no_artifacts_TT08.ntt', 0, 1, [], ...
    [1 1 1 1 1 1] , Timestamps_r, ScNumbers_r, CellNumbers_r,...
                  Features_r, Samples_r, Header);
              
  Mat2NlxSpike( 'artifacts_TT08.ntt', 0, 1, [], ...
    [1 1 1 1 1 1] , Timestamps_w, ScNumbers_w, CellNumbers_w,...
                  Features_w, Samples_w, Header);
