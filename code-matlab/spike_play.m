cd('D:\M074\M074-2020-12-04');
cfg_in.uint= '64';
S = LoadSpikes([cfg_in]);
evs = LoadEvents([]);
all_OFF_times = evs.t{11};

%% Check pre-trial stim
p_start = evs.t{2};
p_end = evs.t{8};
p_on = evs.t{12};
p_off = all_OFF_times(all_OFF_times > p_start & all_OFF_times < p_end);
p_iv = iv(p_start, p_end);
p_S = restrict(S, p_iv);
p_spikes = p_S.t{:};
%%
% Need to find bad spikes, i.e spikes that occur witin