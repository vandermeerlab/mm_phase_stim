cd('/Users/manishm/Dropbox (Dartmouth College)/manish_data/M074/M074-2020-12-04');
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


% spline_type1 (on_idx to on_idx +1)
% spline_type2 (actual_on_idx to actual_on_idx+ 0.001 sec)
% spline_type3 (actual_on_idx to actual_on_idx+ 0.0022 msec)
% spline_type4 (actual_on_idx to actual_on_idx+ 0.0027 msec)

spline_t1 = interpolateArtifacts('spline', laser_on, 0.001, trial_lfp);
spline_t2 = interpolateArtifacts('spline', actual_on, 0.001, trial_lfp);
spline_t3 = interpolateArtifacts('spline', actual_on, 0.0022, trial_lfp);
spline_t4 = interpolateArtifacts('spline', actual_on, 0.0027, trial_lfp);

%% Find phases

f_list = {[3 5], [7 10],[15 25], [30 40],[40 60], [60 80]};
f_list_label = {'3 - 5', '7 - 10', '15 - 25', '30 - 40', '40 - 60', '60 - 80'};
fstop_list = {[2.5 5.5], [5.5 10.5],[14 26], [28 42],[38 62], [55 85]};

