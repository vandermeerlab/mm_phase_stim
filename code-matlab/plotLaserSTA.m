cd('D:\Dropbox (Dartmouth College)\manish_data\M074\M074-2020-12-04');
% cd('D:\Dropbox (Dartmouth College)\manish_data\M078\M078-2020-11-28');
% cd('D:\Dropbox (Dartmouth College)\manish_data\M077\M077-2021-02-26')
evs = LoadEvents([]);
csc_fn = cell(1,16);
csc_fn = cell(1,32);
for iF = 1:32
   csc_fn{iF} = strcat('CSC',num2str(iF),'.ncs');
end
% csc_fn(17:18) = {'LFP3.ncs', 'LFP15.ncs'}; % {'LFP3.ncs', 'LFP15.ncs', 'LFP18.ncs', 'LFP29.ncs'};
% t_start = evs.t{5}; % Fixed ISI protocol start
% t_end = evs.t{3}; % post trial baseline recording statt

csc_fn(33:36) = {'LFP4.ncs', 'LFP6.ncs', 'LFP28.ncs', 'LFP30.ncs'};
t_start = evs.t{8}; % Fixed ISI protocol start
t_end = evs.t{7}; % post trial baseline recording statt
t_iv = iv(t_start, t_end);
laser_on = evs.t{13};
actual_on = laser_on + 0.0011;
for iF = 1:length(csc_fn)
    close all;
    if iF ==4 % || iF == 9 || iF == 11 || iF == 14 || iF ==16
        continue
    end
    fig = figure('WindowState', 'maximized');
    this_cfg.fc = csc_fn(iF);
    this_csc = LoadCSC(this_cfg);
    this_csc = restrict(this_csc, t_iv);
    this_ref = str2num(this_csc.cfg.hdr{1, 1}.ReferenceChannel(end-2:end-1));
    this_ad = this_csc.cfg.hdr{1, 1}.ADChannel;
    Fs = this_csc.cfg.hdr{1}.SamplingFrequency;
    wsize = round(Fs*4);%round(Fs/2.5);  % Let Window size be 4 sec
    [Pxx, F] = pwelch(this_csc.data(:), rectwin(wsize), round(wsize/2), [], Fs);
    subplot(4,1,1)
    plot(F, 10*log10(Pxx));
    xlim([0,120]);
    xlabel('Frequency (Hz)');
    title('PSD');
    subplot(4,1,2)
    % Window size is 4 sec
    [S,F,T,P] = spectrogram(this_csc.data,hanning(round(Fs*4)),round(Fs*2),1:120,Fs);
    imagesc(T,F,10*log10(P)); % converting to dB as usual
    yline(60, 'LineWidth', 1);
    axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)'); 
    title('Spectrogram')
    
    
    w = [-.1 .1]; % time window to compute STA over
    this_tvec = w(1):1/Fs:w(2); % time axis for STA
    for iEvt = 1:length(laser_on) % for each stim ...
 
        on_sta_t = laser_on(iEvt)+w(1);
        actual_sta_t = actual_on(iEvt) + w(1);
        this_on_sta_idx = (nearest_idx3(on_sta_t,this_csc.tvec));
        this_actual_sta_idx = (nearest_idx3(actual_sta_t,this_csc.tvec));
        % grab LFP snippet for this window
        this_on_toAdd = this_csc.data(this_on_sta_idx:this_on_sta_idx+ ...
            length(this_tvec)-1);
        this_actual_toAdd = this_csc.data(this_actual_sta_idx:...
            this_actual_sta_idx+length(this_tvec)-1); 
        this_on_sta(iEvt,:) = this_on_toAdd';
        this_actual_sta(iEvt,:) = this_actual_toAdd'; 
    end
    on_sta = mean(this_on_sta,1);
    actual_sta = mean(this_actual_sta,1);
    subplot(4,1,3);
    plot(this_tvec,on_sta);
    xline(0);
    title('STA aligned to Laser ON')
    subplot(4,1,4);
    plot(this_tvec,actual_sta);
    xline(0);
    title('STA aligned to Actual ON')
    main_title = {this_csc.label{1}(1:end-4), ' AD Channel: ', num2str(this_ad), ' Reference AD: ', num2str(this_ref)};
    main_title = strjoin(main_title, ' ');
    suptitle(main_title)
    WriteFig(fig,csc_fn{iF});
    clear this_on_sta; clear this_actual_sta;
end