cd('data\IceCreamBet_Part2');
% cd('data/M074/M074-2020-12-04');
% cd('data/M078/M078-2020-11-26');
% cd ('data/M078/M078-2020-11-28');
% cd('data/M077/M077-2021-02-25');
cfg_in.fc = {'Sig1.ncs', 'CopyOfSig1.ncs'};
% cfg_in.fc = {'LFP3.ncs', 'LFP15.ncs', 'LFP18.ncs', 'LFP29.ncs'};
% cfg_in.fc = {'LFP3.ncs', 'LFP15.ncs'};

all_lfp  = LoadCSC(cfg_in);
evs = LoadEvents([]);
all_lfp1 = restrict(all_lfp, iv([evs.t{2}, evs.t{1}]));
all_lfp2 = restrict(all_lfp, iv([evs.t{1}, evs.t{3}]));

%%
Fs = all_lfp1.cfg.hdr{1}.SamplingFrequency;
wsize = 4*Fs;
figure;
for i = 1:length(all_lfp1.label)
    [Pxx, F] = pwelch(all_lfp1.data(i,:), rectwin(wsize), wsize/2, [], Fs);
    plot(F, 10*log10(Pxx));
    hold on
end
legend(all_lfp1.label)
xlim([0,120]);

%%
Fs = all_lfp2.cfg.hdr{1}.SamplingFrequency;
wsize = 4*Fs;
figure;
for i = 1:length(all_lfp2.label)
    [Pxx, F] = pwelch(all_lfp2.data(i,:), rectwin(wsize), wsize/2, [], Fs);
    plot(F, 10*log10(Pxx));
    hold on
end
legend(all_lfp2.label)
xlim([0,120]);