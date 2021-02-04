cd('/Users/manishm/Work/vanDerMeerLab/Data/M078/M078-2020-11-26/')
evs = LoadEvents([]);

% The delay between in the time stamps of when the ttl reaches neuralynx
% (evs.t{9}) from Master-8 and after routing through cyclops + op-amp
% (evs.t{14}) is in the order of 10s of microseconds. Thus we can ignore
%  evs.t{14} and evs.t{15}

% All Laser off_events
laser_off = evs.t{9};

%Extract all records from the Tetorde File:
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = ...
Nlx2MatSpike( './TT1.ntt', [1 1 1 1 1], 1, 1, []);


