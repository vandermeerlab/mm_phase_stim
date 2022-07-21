cd('E:\Dropbox (Dartmouth College)\manish_data\M321\M321-2022-07-13');
cfg_test.min_cluster_quality = 3;
cfg_test.uint = '64';

S = LoadSpikes(cfg_test);
evs = LoadEvents([]);
%%
this_on_events = evs.t{15};
for i = 1: length(S.t)
    myCell = SelectTS([], S, i);
    SpikePETHvdm([], myCell, this_on_events, S.label{i}, 0.05);
    close;
end