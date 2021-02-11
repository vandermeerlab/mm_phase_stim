cd('D:\M078\M078-2020-11-26');
evs = LoadEvents([]);

% The delay between in the time stamps of when the ttl reaches neuralynx
% (evs.t{9}) from Master-8 and after routing through cyclops + op-amp
% (evs.t{14}) is in the order of 10s of microseconds. Thus we can ignore
%  evs.t{14} and evs.t{15}

% All Laser off_events
laser_off = evs.t{9};

%Extract all records from the Tetorde File:
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = ...
    Nlx2MatSpike( 'TT6.ntt', [1 1 1 1 1], 1, 1, []);

% Works, Uncomment to run
% % Writing back to a dummy ntt
% Mat2NlxSpike( 'new_TT05.ntt', 0, 1, [], ...
%     [1 1 1 1 1 1] , Timestamps, ScNumbers, CellNumbers,...
%                   Features, Samples, Header);

%%
% Get all Pulse_ON_events
all_on = zeros(1,length(evs.t{11}));
k = 1;
for i = 10:13
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
for i = 1:length(all_on)
    this_tbr = find(temp_ts > all_on(i) & temp_ts < all_on(i)+1e-3);
    to_be_removed(this_tbr) = 1;
end
disp(sum(to_be_removed));
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
 Mat2NlxSpike( 'no_artifacts_TT06.ntt', 0, 1, [], ...
    [1 1 1 1 1 1] , Timestamps_r, ScNumbers_r, CellNumbers_r,...
                  Features_r, Samples_r, Header);
              
  Mat2NlxSpike( 'artifacts_TT06.ntt', 0, 1, [], ...
    [1 1 1 1 1 1] , Timestamps_w, ScNumbers_w, CellNumbers_w,...
                  Features_w, Samples_w, Header);
