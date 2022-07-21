cd('E:\Dropbox (Dartmouth College)\manish_data\M321\M321-2022-07-13\');
% cfg_in.fc ={'LFP29.ncs'};% {'LFP15_0001.ncs', 'LFP18_0001.ncs', 'LFP29_0001.ncs'};
cfg_in.fc = {'LFP29.ncs', 'LFP15.ncs', 'LFP18.ncs', 'LFP3.ncs'};%, 'LFP15.ncs', 'LFP18.ncs', 'LFP29.ncs'};
all_lfp  = LoadCSC(cfg_in);
%%
Fs = all_lfp.cfg.hdr{1}.SamplingFrequency;
wsize = 4*Fs;
figure;
for i = 1:length(all_lfp.label)    
    [Pxx, F] = pwelch(all_lfp.data(i,:), rectwin(wsize), wsize/2, [], Fs);
    plot(F, 10*log10(Pxx));
    hold on
    title(cfg_in.fc{i});
end
legend(all_lfp.label)
xlim([0,120]);