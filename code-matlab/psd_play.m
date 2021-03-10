cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/M075/M075-2021-01-26');
% cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/M074/M074-2020-12-04');
% cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/M078/M078-2020-11-26');
% cd ('/Users/manishm/Dropbox (Dartmouth College)/manish_data/M078/M078-2020-11-28');
% cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/M077/M077-2021-02-25');
cfg_in.fc = {'LFP4.ncs', 'LFP6.ncs', 'LFP28.ncs', 'LFP30.ncs'};
% cfg_in.fc = {'LFP3.ncs', 'LFP15.ncs', 'LFP18.ncs', 'LFP29.ncs'};
% cfg_in.fc = {'LFP3.ncs', 'LFP15.ncs'};

all_lfp  = LoadCSC(cfg_in);
%%
Fs = all_lfp.cfg.hdr{1}.SamplingFrequency;
wsize = 1024;
figure;
for i = 1:length(all_lfp.label)
    [Pxx, F] = pwelch(all_lfp.data(i,:), rectwin(wsize), wsize/2, [], Fs);
    plot(F, 10*log10(Pxx));
    hold on
end
legend(all_lfp.label)
xlim([0,120]);