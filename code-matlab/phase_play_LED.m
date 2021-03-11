cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/M078/M078-2020-11-26');
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

phase_estimates = cell(length(f_list), length(lfp_prefilt));
all_phase = cell(length(f_list), length(lfp_prefilt));
for iF = 1:length(f_list)
     cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
     for iL = 1:length(lfp_prefilt)
         this_filt = FilterLFP(cfg_filt, lfp_prefilt{iL});
         this_res = angle(hilbert(this_filt.data));
         all_phase{iF,iL} = this_res;
         phase_estimates{iF,iL} = this_res(:,on_idx);

     end
end

%% Diagnosis to see what's wrong with Hilbert

% Take a few ON_indices to see how the raw, filtered and hilbert
% transformed signal looks like

% for i = 1:4
%     for j = 1:length(f_list)
%         figure;
%         cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{j};
%         this_filt = FilterLFP(cfg_filt, lfp_prefilt{1});
%         this_angle = angle(hilbert(this_filt.data));
%         t_angle = angle(hilbert(this_filt.data(1,:)));
%         subplot(3,1,1);
%         plot(trial_lfp.data(1,on_idx(1)-1000:on_idx(1)+1000))
%         hold on;
%         plot(this_filt.data(1,on_idx(1)-1000:on_idx(1)+1000))
%         subplot(3,1,2);
%         plot(this_angle(1,on_idx(1)-1000:on_idx(1)+1000))
%         subplot(3,1,3);
%         plot(t_angle(on_idx(1)-1000:on_idx(1)+1000))
%         close;
%     end
% end
%% Eric's method (no spline, phase at led_on)
figure;
k = 1;
for i = 1:4
    for j = 1:6
        subplot(4,6,k);
        hist(phase_estimates{j,1}(i,:),5);
        k = k+1;
    end
end
%%
figure;
k = 1;
for i = 1:4
    for j = 1:6
        subplot(4,6,k);
        hist(all_phase{j,1}(i,:),5);
        k = k+1;
    end
end


%% No spline + Stim_delay
figure;
k = 1;
for i = 1:4
    for j = 1:6
        subplot(4,6,k);
        hist(phase_estimates{j,1,2}(i,:),5);
        k = k+1;
    end
end
%% 1 msec spline + stim_delay
figure;
k = 1;
for i = 1:4
    for j = 1:6
        subplot(4,6,k);
        hist(phase_estimates{j,3,2}(i,:),5);
        k = k+1;
    end
end

%% 2.2 msec spline + stim_delay
figure;
k = 1;
for i = 1:4
    for j = 1:6
        subplot(4,6,k);
        hist(phase_estimates{j,4,2}(i,:),5);
        k = k+1;
    end
end

%% 2.7 msec spline + stim_delay
figure;
k = 1;
for i = 1:4
    for j = 1:6
        subplot(4,6,k);
        hist(phase_estimates{j,5,2}(i,:),5);
        k = k+1;
    end
end
