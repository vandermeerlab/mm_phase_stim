cd('data/M074/M074-2020-12-04');
cfg_in.fc = {'LFP4.ncs', 'LFP6.ncs', 'LFP28.ncs', 'LFP30.ncs'};
all_lfp  = LoadCSC(cfg_in);
evs = LoadEvents([]);
%%

t_start = evs.t{8}; % Fixed ISI protocol start
t_end = evs.t{7}; % post trial baseline recording statt
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

laser_on = evs.t{13};
on_idx = nearest_idx3(laser_on, trial_lfp.tvec);
actual_on = laser_on + 0.0011;
actual_on_idx = nearest_idx3(actual_on, trial_lfp.tvec);


% spline_type1 (actual_on_idx to actual_on_idx+ 0.001 sec)
% spline_type2 (actual_on_idx to actual_on_idx+ 0.0022 sec)
% spline_type3 (actual_on_idx to actual_on_idx+ 0.0027 sec)

spline_t1 = interpolateArtifacts('spline', actual_on, 0.001, trial_lfp);
spline_t2 = interpolateArtifacts('spline', actual_on, 0.0022, trial_lfp);
spline_t3 = interpolateArtifacts('spline', actual_on, 0.0027, trial_lfp);

%% Find phases using FilterLFP + Hilbert transform

% make a cell array of the various processed LFPS
lfp_prefilt{1} = trial_lfp;
lfp_prefilt{2} = spline_t1;
lfp_prefilt{3} = spline_t2;
lfp_prefilt{4} = spline_t3;


f_list = {[3 5], [7 10],[15 25], [30 40],[40 60], [60 80]};
f_list_label = {'3 - 5', '7 - 10', '15 - 25', '30 - 40', '40 - 60', '60 - 80'};
fstop_list = {[2.5 5.5], [5.5 10.5],[14 26], [28 42],[38 62], [55 85]};

% Save filtered LFP and all_phase for later diagnosis
lfp_filtered = cell(length(f_list), length(lfp_prefilt));
all_phase = cell(length(f_list), length(lfp_prefilt));

phase_estimates = cell(length(f_list), length(lfp_prefilt), 2);

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
         phase_estimates{iF,iL,1} = this_res(:,on_idx);
         phase_estimates{iF,iL,2} = this_res(:,actual_on_idx);
     end
end

%% STA for artifact detection

w = [-.5 .5]; % time window to compute STA over
tvec = w(1):1/Fs:w(2); % time axis for STA

 
for iEvt = 1:length(laser_on) % for each stim ...
 
   on_sta_t = laser_on(iEvt)+w(1);
   actual_sta_t = actual_on(iEvt) + w(1);
   on_sta_idx = nearest_idx3(on_sta_t,trial_lfp.tvec); % find index of leading window edge
   actual_sta_idx = nearest_idx3(actual_sta_t,trial_lfp.tvec); % find index of leading window edge
    
   on_toAdd = trial_lfp.data(2,on_sta_idx:on_sta_idx+length(tvec)-1); % grab LFP snippet for this window
   actual_toAdd = trial_lfp.data(2,actual_sta_idx:actual_sta_idx+length(tvec)-1); % grab LFP snippet for this window
   % note this way can be dangerous if there are gaps in the data
 
   on_sta(iEvt,:) = on_toAdd';
   actual_sta(iEvt,:) = actual_toAdd'; 
end
figure;
q1 = mean(on_sta,1);
q2 = mean(actual_sta,1);

subplot(2,1,1)
plot(q1);
hold on
vline(1334);
title('Aligned with Laser ON')


subplot(2,1,2)
plot(q2);
hold on;
vline(1334);
title('Aligned with Actual ON')

suptitle('STA for laser stim')

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


%% Eric's method (no spline, phase at laser_on)
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
                hist(phase_estimates{j-1,1,1}(i-1,:),5);
        end
        k = k+1;
    end
    suptitle('No spline, phase at laser ON')
end
%% No spline + Stim_delay
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
                hist(phase_estimates{j-1,1,2}(i-1,:),5);
        end
        k = k+1;
    end
    suptitle('No spline, phase at laser ON + 1.1 msec')
end
%% 1 msec spline + stim_delay
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
                hist(phase_estimates{j-1,2,2}(i-1,:),5);
        end
        k = k+1;
    end
    suptitle('1 msec spline from laser ON + 1.1 msec')
end

%% 2.2 msec spline + stim_delay
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
                hist(phase_estimates{j-1,3,2}(i-1,:),5);
        end
        k = k+1;
    end
    suptitle('2.2 msec spline from laser ON + 1.1 msec')
end


%% 2.7 msec spline + stim_delay
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
                hist(phase_estimates{j-1,3,2}(i-1,:),5);
        end
        k = k+1;
    end
    suptitle('2.7 msec spline from laser ON + 1.1 msec')
end

