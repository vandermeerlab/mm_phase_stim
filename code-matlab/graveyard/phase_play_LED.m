cd('data/M078/M078-2020-11-26');
cfg_in.fc = {'LFP3.ncs', 'LFP15.ncs', 'LFP18.ncs', 'LFP29.ncs'};
all_lfp  = LoadCSC(cfg_in);
evs = LoadEvents([]);
%%

t_start = evs.t{6}; % Fixed ISI protocol start
t_end = evs.t{5}; % post trial baseline recording statt
t_iv = iv(t_start, t_end);
trial_lfp  = restrict(all_lfp, t_iv);

%% Check PSD
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

%% 

led_on = evs.t{11};
on_idx = nearest_idx3(led_on, trial_lfp.tvec);


% spline_type1 (actual_on_idx to actual_on_idx+ 0.001 sec)
% spline_type2 (actual_on_idx to actual_on_idx+ 0.0015 sec)

spline_t1 = interpolateArtifacts('spline', led_on, 0.001, trial_lfp);
spline_t2 = interpolateArtifacts('spline', led_on, 0.0015, trial_lfp);

%% Find phases using FilterLFP + Hilbert transform

% make a cell array of the various processed LFPS
lfp_prefilt{1} = trial_lfp;
lfp_prefilt{2} = spline_t1;
lfp_prefilt{3} = spline_t2;

f_list = {[3 5], [7 10],[15 25], [30 40],[40 60], [60 80]};
f_list_label = {'3 - 5', '7 - 10', '15 - 25', '30 - 40', '40 - 60', '60 - 80'};
fstop_list = {[2.5 5.5], [5.5 10.5],[14 26], [28 42],[38 62], [55 85]};

% Save filtered LFP and all_phase for later diagnosis
lfp_filtered = cell(length(f_list), length(lfp_prefilt));
all_phase = cell(length(f_list), length(lfp_prefilt));

phase_estimates = cell(length(f_list), length(lfp_prefilt));
for iF = 1:length(f_list)
     cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
     for iL = 1:length(lfp_prefilt)
         this_filt = FilterLFP(cfg_filt, lfp_prefilt{iL});
         lfp_filtered{iF,iL} = this_filt;
         this_res = zeros(size(this_filt.data));
         for iC = 1:4
            this_res(iC,:) = angle(hilbert(this_filt.data(iC,:)));
         end
         all_phase{iF,iL} = this_res;
         phase_estimates{iF,iL} = this_res(:,on_idx);

     end
end

%% STA for artifact detection

w = [-.5 .5]; % time window to compute STA over
tvec = w(1):1/Fs:w(2); % time axis for STA

 
for iEvt = 1:length(led_on) % for each stim ...
 
   on_sta_t = led_on(iEvt)+w(1);
   on_sta_idx = nearest_idx3(on_sta_t,trial_lfp.tvec); % find index of leading window edge
    
   on_toAdd = trial_lfp.data(4,on_sta_idx:on_sta_idx+length(tvec)-1); % grab LFP snippet for this window
   % note this way can be dangerous if there are gaps in the data
 
   on_sta(iEvt,:) = on_toAdd';
end
figure;
q1 = mean(on_sta,1);
plot(q1);
hold on
vline(1334);
title('Aligned with LED ON')

suptitle('STA for LED stim')

%% histogram for all_phase
figure;
k = 1;
for i = 1:5
    for j = 1:7
        switch k
            case 1
               disp('Skipped')
            case 2
                subplot(5,7,k)
                text(0.5,0.5,f_list_label(1), 'FontSize', 14)
                axis off
            case 3
                subplot(5,7,k)
                text(0.5,0.5,f_list_label(2), 'FontSize', 14)
                axis off
            case 4
                subplot(5,7,k)
                text(0.5,0.5,f_list_label(3), 'FontSize', 14)
                axis off
            case 5
                subplot(5,7,k)
                text(0.5,0.5,f_list_label(4), 'FontSize', 14)
                axis off
            case 6
                subplot(5,7,k)
                text(0.5,0.5,f_list_label(5), 'FontSize', 14)
                axis off
            case 7
                subplot(5,7,k)
                text(0.5,0.5,f_list_label(6), 'FontSize', 14)
                axis off
            case 8
                subplot(5,7,k)
                text(0.5,0.5,trial_lfp.label{1}, 'FontSize', 14)
                axis off
            case 15
                subplot(5,7,k)
                text(0.5,0.5,trial_lfp.label{2}, 'FontSize', 14)
                axis off
            case 22
                subplot(5,7,k)
                text(0.5,0.5,trial_lfp.label{3}, 'FontSize', 14)
                axis off
            case 29
                subplot(5,7,k)
                text(0.5,0.5,trial_lfp.label{4}, 'FontSize', 14)
                axis off
            otherwise
                subplot(5,7,k)
                hist(all_phase{j-1,1}(i-1,:),5);
        end
        k = k+1;
    end
    suptitle('Phase across all LFP samples')
end
