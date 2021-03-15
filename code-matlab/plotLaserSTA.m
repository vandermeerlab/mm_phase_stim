cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/M074/M074-2020-12-04');
evs = LoadEvents([]);
csc_fn = cell(1,32);
for iF = 1:32
   csc_fn{iF} = strcat('CSC',num2str(iF),'.ncs');
end
csc_fn(33:36) = {'LFP4.ncs', 'LFP6.ncs', 'LFP28.ncs', 'LFP30.ncs'};
t_start = evs.t{8}; % Fixed ISI protocol start
t_end = evs.t{7}; % post trial baseline recording statt
t_iv = iv(t_start, t_end);
for iF = 1:length{csc_fn}
    fig = figure('WindowState', 'maximized');
    this_cfg.fc = csc_fn(iF);
    this_lfp = LoadCSC(this_cfg);
    this_lfp = restrict(this_lfp, t_iv);
    Fs = trial_csc.cfg.hdr{1}.SamplingFrequency;
    wsize = Fs/2.5;  % arbitrary decision on Window Size
    [Pxx, F] = pwelch(this_lfp.data(:), rectwin(wsize), wsize/2, [], Fs);
    subplot(4,1,1)
    plot(F, 10*log10(Pxx));
    xlim([0,120]);
    title('PSD')
    
    
end