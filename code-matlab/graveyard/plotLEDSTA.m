cd('data\M078\M078-2020-11-28');
evs = LoadEvents([]);
csc_fn = cell(1,16); % cell(1,32);
for iF = 1:16
   csc_fn{iF} = strcat('CSC',num2str(iF),'.ncs');
end
csc_fn(17:18) = {'LFP3.ncs', 'LFP15.ncs'}; % {'LFP3.ncs', 'LFP15.ncs', 'LFP18.ncs', 'LFP29.ncs'};
t_start = evs.t{5}; % Fixed ISI protocol start
t_end = evs.t{3}; % post trial baseline recording statt
t_iv = iv(t_start, t_end);
led_on = evs.t{10};
control_on = led_on + 0.5;
for iF = 17:length(csc_fn)
    close all;
    if iF ==4
        continue
    end
    fig = figure('WindowState', 'maximized');
    this_cfg.fc = csc_fn(iF);
    this_csc = LoadCSC(this_cfg);
    this_csc = restrict(this_csc, t_iv);
    Fs = this_csc.cfg.hdr{1}.SamplingFrequency;
    wsize = round(Fs/2.5);  % arbitrary decision on Window Size
    [Pxx, F] = pwelch(this_csc.data(:), rectwin(wsize), round(wsize/2), [], Fs);
    subplot(4,1,1)
    plot(F, 10*log10(Pxx));
    xlim([0,120]);
    xlabel('Frequency (Hz)');
    title('PSD');
    subplot(4,1,2)
    [S,F,T,P] = spectrogram(this_csc.data,hanning(round(Fs/8)),round(Fs/16),1:120,Fs);
    imagesc(T,F,10*log10(P)); % converting to dB as usual
    yline(60, 'LineWidth', 1);
    axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)'); 
    title('Spectrogram')
    
    
    w = [-.1 .1]; % time window to compute STA over
    this_tvec = w(1):1/Fs:w(2); % time axis for STA
    for iEvt = 1:length(led_on) % for each stim ...
        on_sta_t = led_on(iEvt)+w(1);
        this_on_sta_idx = (nearest_idx3(on_sta_t,this_csc.tvec));
        % grab LFP snippet for this window
        this_on_toAdd = this_csc.data(this_on_sta_idx:this_on_sta_idx+ ...
            length(this_tvec)-1);
        this_on_sta(iEvt,:) = this_on_toAdd';
        
        control_on_sta_t = control_on(iEvt)+w(1);
        control_on_sta_idx = (nearest_idx3(control_on_sta_t,this_csc.tvec));
        % grab LFP snippet for this window
        control_on_toAdd = this_csc.data(control_on_sta_idx:control_on_sta_idx+ ...
            length(this_tvec)-1);
        control_on_sta(iEvt,:) = control_on_toAdd';
    end
    on_sta = mean(this_on_sta,1);
    control_sta = mean(control_on_sta,1);
    subplot(4,1,3);
    plot(this_tvec, on_sta);
    xline(0);
    title('STA aligned to LED ON')
    subplot(4,1,4);
    plot(this_tvec, control_sta);
    xline(0);
    title('STA aligned to Control ON')
    WriteFig(fig,csc_fn{iF});
    clear this_on_sta; clear control_on_sta;
end