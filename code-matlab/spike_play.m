% cd('D:\M074\new_artifacts');
cd('D:\M075\artifacts');
cfg_in.uint= '64';
S = LoadSpikes([cfg_in]);
evs = LoadEvents([]);
all_OFF_times = evs.t{9};

%% Check pre-trial stim
p_start = evs.t{15};
p_end = evs.t{5};
p_on = evs.t{10};
p_off = all_OFF_times(all_OFF_times > p_start & all_OFF_times < p_end);
p_iv = iv(p_start, p_end);
p_S = restrict(S, p_iv);
%%

idx1 = nearest_idx3(p_S.t{1},p_on);
idx2 = nearest_idx3(p_S.t{2},p_on);

offsets_1 = p_S.t{1} - p_on(idx1);
offsets_2 = p_S.t{2} - p_on(idx2);
