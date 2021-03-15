cd('D:\Dropbox (Dartmouth College)\manish_data\M074\M074-2020-12-04');
cfg_in.fc = {'LFP30.ncs'};
all_lfp  = LoadCSC(cfg_in);
cfg_in.fc = {'CSC30.ncs', 'CSC11.ncs'};
all_csc = LoadCSC(cfg_in);
evs = LoadEvents([]);
psd_plot = 0; % set to 1 if you want to see psds
%%
t_start = evs.t{8}; % Fixed ISI protocol start
t_end = evs.t{7}; % post trial baseline recording statt
t_iv = iv(t_start, t_end);
trial_lfp  = restrict(all_lfp, t_iv);
trial_csc = restrict(all_csc, t_iv);
clear all_csc
clear all_lfp
%% Check PSD
if psd_plot
    Fs = trial_lfp.cfg.hdr{1}.SamplingFrequency;
    wsize = 1024;
    figure;
    for i = 1:length(trial_lfp.label)
        [Pxx, F] = pwelch(trial_lfp.data(i,:), rectwin(wsize), wsize/2, [], Fs);
        plot(F, 10*log10(Pxx));
        hold on
    end
    legend(trial_lfp.label)
    xlim([0,120]);
    Fs = trial_csc.cfg.hdr{1}.SamplingFrequency;
    wsize = 8192; %Need to fix this
    figure;
    for i = 1%:length(trial_lfp.label)
        [Pxx, F] = pwelch(trial_csc.data(i,:), rectwin(wsize), wsize/2, [], Fs);
        plot(F, 10*log10(Pxx));
        hold on
    end
    legend(trial_csc.label)
    xlim([0,120]);
end
%%
laser_on = evs.t{13};
lfp_on_idx = nearest_idx3(laser_on, trial_lfp.tvec);
csc_on_idx = nearest_idx3(laser_on, trial_csc.tvec);
actual_on = laser_on + 0.0011;
lfp_actual_on_idx = nearest_idx3(actual_on, trial_lfp.tvec);
csc_actual_on_idx = nearest_idx3(actual_on, trial_csc.tvec);

%% Find phases using FilterLFP + Hilbert transform

f_list = {[3 5], [7 10],[15 25], [30 40],[40 60], [60 80]};
f_list_label = {'3 - 5', '7 - 10', '15 - 25', '30 - 40', '40 - 60', '60 - 80'};
fstop_list = {[2.5 5.5], [5.5 10.5],[14 26], [28 42],[38 62], [55 85]};

% Save filtered CSC and all_phase for later diagnosis
csc_filtered = cell(length(f_list));
all_csc_phase = cell(length(f_list));
csc_phase_estimates = cell(length(f_list), 2);

% Save filtered LFP and all_phase for later diagnosis
lfp_filtered = cell(length(f_list));
all_lfp_phase = cell(length(f_list));
lfp_phase_estimates = cell(length(f_list), 2);


for iF = 1:length(f_list)
    % Process CSCs first
     cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
     this_filt = FilterLFP(cfg_filt, trial_csc);
     csc_filtered{iF} = this_filt;
     this_res = zeros(size(this_filt.data));
     for iC = 1:length(trial_csc.label)
        this_res(iC,:) = angle(hilbert(this_filt.data(iC,:)));
     end
     all_csc_phase{iF} = this_res;
     csc_phase_estimates{iF,1} = this_res(:,csc_on_idx);
     csc_phase_estimates{iF,2} = this_res(:,csc_actual_on_idx);
     
     % Process LFPs next
     cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
     this_filt = FilterLFP(cfg_filt, trial_lfp);
     lfp_filtered{iF} = this_filt;
     this_res = zeros(size(this_filt.data));
     for iC = 1:length(trial_lfp.label)
        this_res(iC,:) = angle(hilbert(this_filt.data(iC,:)));
     end
     all_lfp_phase{iF} = this_res;
     lfp_phase_estimates{iF,1} = this_res(:,lfp_on_idx);
     lfp_phase_estimates{iF,2} = this_res(:,lfp_actual_on_idx);
end
%% STA for artifact detection

w = [-.1 .1]; % time window to compute STA over
lfp_Fs = trial_lfp.cfg.hdr{1}.SamplingFrequency;
lfp_tvec = w(1):1/lfp_Fs:w(2); % time axis for STA

csc_Fs = trial_csc.cfg.hdr{1}.SamplingFrequency;
csc_tvec = w(1):1/csc_Fs:w(2); % time axis for STA

for iEvt = 1:length(laser_on) % for each stim ...
 
   on_sta_t = laser_on(iEvt)+w(1);
   actual_sta_t = actual_on(iEvt) + w(1);
   
   lfp_on_sta_idx = nearest_idx3(on_sta_t,trial_lfp.tvec); % find index of leading window edge
   lfp_actual_sta_idx = nearest_idx3(actual_sta_t,trial_lfp.tvec); % find index of leading window edge
   
   csc_on_sta_idx = nearest_idx3(on_sta_t,trial_csc.tvec); % find index of leading window edge
   csc_actual_sta_idx = nearest_idx3(actual_sta_t,trial_csc.tvec); % find index of leading window edge
   
    
   lfp_on_toAdd = trial_lfp.data(lfp_on_sta_idx:lfp_on_sta_idx+length(lfp_tvec)-1); % grab LFP snippet for this window
   lfp_actual_toAdd = trial_lfp.data(lfp_actual_sta_idx:lfp_actual_sta_idx+length(lfp_tvec)-1); % grab LFP snippet for this window

   csc_on_toAdd = trial_csc.data(2,csc_on_sta_idx:csc_on_sta_idx+length(csc_tvec)-1); % grab LFP snippet for this window
   csc_actual_toAdd = trial_csc.data(2,csc_actual_sta_idx:csc_actual_sta_idx+length(csc_tvec)-1); % grab LFP snippet for this window

 
   lfp_on_sta(iEvt,:) = lfp_on_toAdd';
   lfp_actual_sta(iEvt,:) = lfp_actual_toAdd'; 
   csc_on_sta(iEvt,:) = csc_on_toAdd';
   csc_actual_sta(iEvt,:) = csc_actual_toAdd'; 
end



%% PLOT STA
figure;
q1 = mean(lfp_on_sta,1);
q2 = mean(lfp_actual_sta,1);

subplot(2,1,1)
plot(q1);
hold on
vline(1334);
title('LFP Aligned with Laser ON')


subplot(2,1,2)
plot(q2);
hold on;
vline(1334);
title('Aligned with Actual ON')
suptitle('STA for laser stim')
%%
figure;
q1 = mean(csc_on_sta,1);
q2 = mean(csc_actual_sta,1);
subplot(2,1,1)
plot(q1);
hold on
vline(16001);
title(' CSC Aligned with Laser ON')


subplot(2,1,2)
plot(q2);
hold on;
vline(16001);
title('Aligned with Actual ON')

suptitle('STA for laser stim')

%% Eric's method (no spline, phase at laser_on)
figure;
k = 1;
for i = 1:3
    for j = 1:7
        switch k
            case 1
               disp('Skipped')
            case 2
                subplot(3,7,k)
                text(0.5,0.5,f_list_label(1), 'FontSize', 14)
                axis off
            case 3
                subplot(3,7,k)
                text(0.5,0.5,f_list_label(2), 'FontSize', 14)
                axis off
            case 4
                subplot(3,7,k)
                text(0.5,0.5,f_list_label(3), 'FontSize', 14)
                axis off
            case 5
                subplot(3,7,k)
                text(0.5,0.5,f_list_label(4), 'FontSize', 14)
                axis off
            case 6
                subplot(3,7,k)
                text(0.5,0.5,f_list_label(5), 'FontSize', 14)
                axis off
            case 7
                subplot(3,7,k)
                text(0.5,0.5,f_list_label(6), 'FontSize', 14)
                axis off
            case 8
                subplot(3,7,k)
                text(0.5,0.5,trial_csc.label{1}, 'FontSize', 14)
                axis off
            case 15
                subplot(3,7,k)
                text(0.5,0.5,trial_csc.label{2}, 'FontSize', 14)
                axis off
            otherwise
                subplot(3,7,k)
                hist(csc_phase_estimates{j-1,1}(i-1,:),5);
        end
        k = k+1;
    end
    suptitle('No spline, phase at laser ON')
end