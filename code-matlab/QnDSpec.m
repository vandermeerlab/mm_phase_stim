cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/M175/M175-2021-04-10');
evs = LoadEvents([]);
csc_fn = {'LFP3.ncs', 'LFP15.ncs', 'LFP18.ncs', 'LFP29.ncs'};

%%
t_start = evs.t{2}; % Start recording
t_end = evs.t{3}; % Stop recording
t_iv = iv(t_start, t_end);
laser_on = evs.t{5};
laser_off = evs.t{4};
this_cfg.fc = {'R1.ncs'};
ref_sig = LoadCSC(this_cfg);
ref_sig = restrict(ref_sig, t_iv);

for iF = 1:length(csc_fn)
    close all;
    fig = figure('WindowState', 'maximized');
    this_cfg.fc = csc_fn(iF);
    this_csc = LoadCSC(this_cfg);
    this_csc = restrict(this_csc, t_iv);
    this_csc.data = this_csc.data - ref_sig.data;
    this_ref = str2num(this_csc.cfg.hdr{1, 1}.ReferenceChannel(end-2:end-1));
    this_ad = this_csc.cfg.hdr{1, 1}.ADChannel;
    Fs = this_csc.cfg.hdr{1}.SamplingFrequency;
    wsize = round(Fs/2.5);%round(Fs*4);%round(Fs/2.5);  % Let Window size be 4 sec
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
%%
    w = [-.1 .1]; % time window to compute STA over
    this_tvec = w(1):1/Fs:w(2); % time axis for STA
    for iEvt = 1:length(laser_on) % for each stim ...
 
        on_sta_t = laser_on(iEvt)+w(1);
        off_sta_t = laser_off(iEvt) + w(1);
        this_on_sta_idx = (nearest_idx3(on_sta_t,this_csc.tvec));
        this_off_sta_idx = (nearest_idx3(off_sta_t,this_csc.tvec));
        % grab LFP snippet for this window
        this_on_toAdd = this_csc.data(this_on_sta_idx:this_on_sta_idx+ ...
            length(this_tvec)-1);
        this_off_toAdd = this_csc.data(this_off_sta_idx:...
            this_off_sta_idx+length(this_tvec)-1); 
        this_on_sta(iEvt,:) = this_on_toAdd';
        this_off_sta(iEvt,:) = this_off_toAdd'; 
    end
    on_sta = mean(this_on_sta,1);
    off_sta = mean(this_off_sta,1);
    subplot(4,1,3);
    plot(this_tvec,on_sta);
    xline(0);
    title('STA aligned to Laser ON')
    subplot(4,1,4);
    plot(this_tvec,off_sta);
    xline(0);
    title('STA aligned to Laser OFF')
    main_title = {this_csc.label{1}(1:end-4), ' AD Channel: ', num2str(this_ad), ' Reference AD: ', num2str(this_ref)};
    main_title = strjoin(main_title, ' ');
    suptitle(main_title)
    WriteFig(fig,csc_fn{iF});
    clear this_on_sta; clear this_actual_sta;
end

%%
all_fn = {'LFP3.ncs', 'LFP15.ncs', 'LFP18.ncs', 'LFP29.ncs', 'R1.ncs'};
all_cscs = cell(size(all_fn));
for i = 1:5
    x.fc = all_fn(i);
    all_cscs{i} = LoadCSC(x);
end
%%
close all
t1_start = 6000000;
t1 = all_cscs{1}.tvec(t1_start);
t1_end = nearest_idx3(t1+5, all_cscs{1}.tvec);
tvec = all_cscs{1}.tvec(t1_start:t1_end) - all_cscs{1}.tvec(t1_start);
figure
plot(tvec, all_cscs{1}.data(t1_start:t1_end));
hold on;
plot(tvec, all_cscs{2}.data(t1_start:t1_end) - 0.003);
plot(tvec, all_cscs{3}.data(t1_start:t1_end) - 0.002);
plot(tvec, all_cscs{4}.data(t1_start:t1_end) - 0.004);
plot(tvec, all_cscs{5}.data(t1_start:t1_end) - 0.008);
legend(all_fn)


%%
tiledlayout(5,1)

% First plot
ax1 = nexttile;
plot(all_cscs{1}.data)
title(all_fn{1})

ax2 = nexttile;
plot(all_cscs{2}.data)
title(all_fn{2})

ax3 = nexttile;
plot(all_cscs{3}.data)
title(all_fn{3})

ax4 = nexttile;
plot(all_cscs{4}.data)
title(all_fn{4})

ax5 = nexttile;
plot(all_cscs{5}.data)
title(all_fn{5})

linkaxes([ax1 ax2 ax3 ax4 ax5],'xy')
linkaxes([ax1 ax2 ax3 ax4 ax5],'off')