%% Script to show the problem

% Example Eric Data
% 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M18\M18-2019-04-12_dStr_3p4_light_cells_TT4_min'% Laser
% 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M18\M18-2019-04-13_dStr_3p8_light_cells_TT6_min' % Laser
% 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M19\M19-2019-04-14_vStr_4p2_light_cells_TT5_min' % Laser

% Example Manish Data
% 'E:\Dropbox (Dartmouth College)\manish_data\M235\M235-2021-07-16' % Laser
% 'E:\Dropbox (Dartmouth College)\manish_data\M325\M325-2022-08-03' % LED

folder = 'E:\Dropbox (Dartmouth College)\manish_data\M325\M325-2022-08-03';
cd(folder);
LoadExpKeys;
evs = LoadEvents([]);

% var_stim_on = evs.t{strcmp(evs.label,ExpKeys.laser_on)} +  0.0011; % Eric data
var_stim_on = evs.t{strcmp(evs.label,ExpKeys.trial_stim_on)};%+ 0.0011; % Manish data

control_on = var_stim_on + 0.33;

%% Load LFP and look at STA
% cfg_lfp.fc = {ExpKeys.goodCSC}; % Eric Data
cfg_lfp.fc = ExpKeys.goodLFP; % Manish Data
if contains(cfg_lfp.fc, '-')
    temp = split(cfg_lfp.fc,'-');
    cfg_lfp.fc = {cat(2,temp{1},'.ncs')};
    this_lfp = LoadCSC(cfg_lfp);
    cfg_temp.fc = {cat(2,temp{2},'.ncs')};
    ref = LoadCSC(cfg_temp);
    this_lfp.data = this_lfp.data - ref.data;
    clear temp ref;
else
    this_lfp = LoadCSC(cfg_lfp);
end
%%
w = [-.1 .1]; % time window to compute STA over
Fs = this_lfp.cfg.hdr{1}.SamplingFrequency;
this_tvec = w(1):1/Fs:w(2); % time axis for STA
for iEvt = 1:length(var_stim_on) % for each stim ...
    on_sta_t = var_stim_on(iEvt)+w(1);
    this_on_sta_idx = (nearest_idx3(on_sta_t,this_lfp.tvec));
    % grab LFP snippet for this window
    this_on_toAdd = this_lfp.data(this_on_sta_idx:this_on_sta_idx+ ...
        length(this_tvec)-1);
    this_on_sta(iEvt,:) = this_on_toAdd';

    control_on_sta_t = control_on(iEvt)+w(1);
    control_on_sta_idx = (nearest_idx3(control_on_sta_t,this_lfp.tvec));
    % grab LFP snippet for this window
    control_on_toAdd = this_lfp.data(control_on_sta_idx:control_on_sta_idx+ ...
        length(this_tvec)-1);
    control_on_sta(iEvt,:) = control_on_toAdd';
end
on_sta = mean(this_on_sta,1);
control_sta = mean(control_on_sta,1);
%%
figure;
hold on;
plot(this_tvec, on_sta);
plot(this_tvec, control_sta);
xline(0);
legend({'STA', 'Control', 't=0'},'Location','southeast');
title('STA');
%%
figure;
sample_idx = randi(length(var_stim_on), 1, 5);
for i = 1:5
    subplot(5,1,i)
    hold on;
    plot(this_tvec, this_on_sta(sample_idx(i),:));
    xline(0);
end
suptitle('Sample Stimulus triggered Snippets')
%% Do filtering stuff
Hf = fdesign.lowpass('N,Fc', 32, 400/this_lfp.cfg.hdr{1}.SamplingFrequency);
Hd = design(Hf, 'window', 'window', @bartlett, 'systemobject', true);
b = Hd.Numerator;
filtered_lfp = this_lfp;
filtered_lfp.data = filter(b,1, this_lfp.data);
%% STA stuff on filtered_data
w = [-.1 .1]; % time window to compute STA over
Fs = filtered_lfp.cfg.hdr{1}.SamplingFrequency;
this_tvec = w(1):1/Fs:w(2); % time axis for STA
for iEvt = 1:length(var_stim_on) % for each stim ...
    filt_on_sta_t = var_stim_on(iEvt)+w(1);
    this_on_sta_idx = (nearest_idx3(filt_on_sta_t,filtered_lfp.tvec));
    % grab LFP snippet for this window
    this_on_toAdd = filtered_lfp.data(this_on_sta_idx:this_on_sta_idx+ ...
        length(this_tvec)-1);
    filt_this_on_sta(iEvt,:) = this_on_toAdd';

    filt_control_on_sta_t = control_on(iEvt)+w(1);
    control_on_sta_idx = (nearest_idx3(filt_control_on_sta_t,filtered_lfp.tvec));
    % grab LFP snippet for this window
    control_on_toAdd = filtered_lfp.data(control_on_sta_idx:control_on_sta_idx+ ...
        length(this_tvec)-1);
    filt_control_on_sta(iEvt,:) = control_on_toAdd';
end
filt_on_sta = mean(filt_this_on_sta,1);
filt_control_sta = mean(filt_control_on_sta,1);
%% Plot filtered STA
figure;
subplot(2,1,1)
hold on;
plot(this_tvec, on_sta);
plot(this_tvec, control_sta);
xline(0);
legend({'STA', 'Control', 't=0'},'Location','southeast');
title('STA');
subplot(2,1,2)
hold on;
plot(this_tvec, filt_on_sta);
plot(this_tvec, filt_control_sta);
xline(0);
legend({'STA', 'Control', 't=0'},'Location','southeast');
title('Filtered')
%% Plot filtered snippets
figure;
sample_idx = randi(length(var_stim_on), 1, 5);
for i = 1:5
    subplot(5,1,i)
    hold on;
    plot(this_tvec, this_on_sta(sample_idx(i),:));
    plot(this_tvec, filt_this_on_sta(sample_idx(i),:));
    xline(0);
end
suptitle('Sample Stimulus triggered Snippets')
%%  First let's look at the two PSDs
wsize = Fs;
[P_Og, F] = pwelch(this_lfp.data, hanning(wsize), wsize/2, [], Fs);
% [P_filt,~] = pwelch(filtered_lfp.data, hanning(wsize), wsize/2, [], Fs);
figure;
hold on;
plot(F, 10*log10(P_Og));
% plot(F, 10*log10(P_filt));
xlim([0 120]);

%% FOOOF (Doesn't work)

F = (F(F > 0 & F < 120))'; % very important to get rid of the '0' frequency for FOOOF to work, and shape it into 1 x N array
P_Og = (P_Og(1:length(F)))';
reshaped_P = reshape(P_Og,1 ,1, length(P_Og));
% All these parameters are borrowed from "process_fooof.m"
opt.freq_range = [F(1) F(end)];
opt.power_line = '60';
opt.peak_width_limits = [0.5,12];
opt.max_peaks = 3;
opt.min_peak_height = 0.3;
opt.aperiodic_mode = 'knee'; %Check with 'fixed' first
opt.peak_threshold = 2;
opt.return_spectrum = 1;
opt.border_threshold = 1;
opt.peak_type = 'best'; %There is an error in documenation where it says 'both'
opt.proximity_threshold = 2;
opt.guess_weight = 'none';
opt.thresh_after = 1;
opt.sort_type = 'param';
opt.sort_param = 'frequency';
%     opt.sort_bands = {{'delta'}, {'2', '4'}; {'theta'}, {'5', '7'}; ...
%         {'alpha'}, {'8','12'}; {'beta'}, {'15', '29'}; {'gamma1'}, {'30',' 59'};
%         {'gamma2'}, {'60','90'}};
opt.sort_bands = {{'delta'}, {'2', '5'}; {'theta'}, {'6', '10'}; ...
    {'beta'},{'15', '29'}; {'gamma1'}, {'30',' 55'}; ...
    {'gamma2'}, {'65','90'}};
[fs, fg] = process_fooof('FOOOF_matlab', reshaped_P, F, opt, 1);
powspctrm_f = cat(1, fg.ap_fit);
for k = 1:size(powspctrm_f,1)
    aperiodic_P(k,:) = interp1(fs, powspctrm_f(k,:), F, 'linear', nan);
end
figure;
plot(F, 10*log10(P_Og))
hold on;
plot(F, 10*log10(aperiodic_P))
for iR = 1:size(fg.peak_params,1)
    xline(fg.peak_params(iR,1), 'black');
end

%% Filter in delta and slow gamma-ish range
f_list = {[2 5], [25 50]};
bp_lfp = cell(length(f_list),1);
lfp_phase = cell(length(f_list),1);
% bp_filt = cell(length(f_list),1);
% filt_phase = cell(length(f_list),1);
for iB = 1:length(f_list)
    cfg_filt.type = 'fdesign'; 
    cfg_filt.f  = f_list{iB};
    bp_lfp{iB} = FilterLFP(cfg_filt, this_lfp);
%     bp_filt{iB} = FilterLFP(cfg_filt, filtered_lfp);
    lfp_phase{iB} = angle(hilbert(bp_lfp{iB}.data));
%     filt_phase{iB} = angle(hilbert(bp_filt{iB}.data));
end
%% Plot some phase scatters
figure
on_idx = nearest_idx3(var_stim_on, this_lfp.tvec);
control_idx = nearest_idx3(control_on, this_lfp.tvec);

subplot(2,2,1)
scatter(lfp_phase{1}(on_idx), filt_phase{1}(on_idx));
title('Delta Phase at Stim times')
xlabel('Unfiltered Phase')
ylabel('Filtered Phase')

subplot(2,2,2)
scatter(lfp_phase{1}(control_idx), filt_phase{1}(control_idx));
title('Delta Phase at Control times')
xlabel('Unfiltered Phase')
ylabel('Filtered Phase')

subplot(2,2,3)
scatter(lfp_phase{2}(on_idx), filt_phase{2}(on_idx));
title('Slowish Gamma Phase at Stim times')
xlabel('Unfiltered Phase')
ylabel('Filtered LFP')

subplot(2,2,4)
scatter(lfp_phase{2}(control_idx), filt_phase{2}(control_idx));
title('Slowish Gamma Phase at Control times')
xlabel('Unfiltered Phase')
ylabel('Filtered LFP')

%% Plot some phase hist
figure
on_idx = nearest_idx3(var_stim_on, this_lfp.tvec);
control_idx = nearest_idx3(control_on, this_lfp.tvec);

subplot(4,2,1)
hist(lfp_phase{1}(on_idx), 5);
title('Unfiltered Delta Phase distribution at Stim times')

subplot(4,2,2)
hist(filt_phase{1}(on_idx), 5);
title('Filtered Delta Phase distribution at Stim times')

subplot(4,2,3)
hist(lfp_phase{1}(control_idx), 5);
title('Unfiltered Delta Phase distribution at Control times')

subplot(4,2,4)
hist(filt_phase{1}(control_idx), 5);
title('Filtered Delta Phase distribution at Control times')

subplot(4,2,5)
hist(lfp_phase{2}(on_idx), 5);
title('Unfiltered slow gamma-ish Phase distribution at Stim times')

subplot(4,2,6)
hist(filt_phase{2}(on_idx), 5);
title('Filtered slow gamma-ish Phase distribution at Stim times')

subplot(4,2,7)
hist(lfp_phase{2}(control_idx), 5);
title('Unfiltered slow gamma-ish Phase distribution at Control times')

subplot(4,2,8)
hist(filt_phase{2}(control_idx), 5);
title('Filtered slow gamma-ish Phase distribution at Control times')

%% FOR Eric's data
figure
on_idx = nearest_idx3(var_stim_on, this_lfp.tvec);
control_idx = nearest_idx3(control_on, this_lfp.tvec);

subplot(2,2,1)
hist(lfp_phase{1}(on_idx), 5);
title('Delta Phase distribution at Stim times')

subplot(2,2,2)
hist(lfp_phase{1}(control_idx), 5);
title('Delta Phase distribution at Control times')

subplot(2,2,3)
hist(lfp_phase{2}(on_idx), 5);
title('Slowish Gamma Phase distribution at Stim times')

subplot(2,2,4)
hist(lfp_phase{2}(control_idx), 5);
title('Slowish Gamma Phase distribution at Control times')


%% Filter LFP where phase doesn't exist

f_list2 = {[6 10],[80 100]};
bp_lfp2 = cell(length(f_list2),1);
lfp_phase2 = cell(length(f_list2),1);
% bp_filt2 = cell(length(f_list2),1);
% filt_phase2 = cell(length(f_list2),1);
for iB = 1:length(f_list2)
    cfg_filt.type = 'fdesign'; 
    cfg_filt.f  = f_list2{iB};
    bp_lfp2{iB} = FilterLFP(cfg_filt, this_lfp);
%     bp_filt2{iB} = FilterLFP(cfg_filt, filtered_lfp);
    lfp_phase2{iB} = angle(hilbert(bp_lfp2{iB}.data));
%     filt_phase2{iB} = angle(hilbert(bp_filt2{iB}.data));
end

%% Plot Sanity check stuff
figure
on_idx = nearest_idx3(var_stim_on, this_lfp.tvec);
control_idx = nearest_idx3(control_on, this_lfp.tvec);

subplot(2,2,1)
scatter(lfp_phase2{1}(on_idx), filt_phase2{1}(on_idx));
title('Theta Phase at Stim times')
xlabel('Unfiltered Phase')
ylabel('Filtered Phase')

subplot(2,2,2)
scatter(lfp_phase2{1}(control_idx), filt_phase2{1}(control_idx));
title('Theta Phase at Control times')
xlabel('Unfiltered Phase')
ylabel('Filtered Phase')

subplot(2,2,3)
scatter(lfp_phase2{2}(on_idx), filt_phase2{2}(on_idx));
title('Fast Gamma Phase at Stim times')
xlabel('Unfiltered Phase')
ylabel('Filtered LFP')

subplot(2,2,4)
scatter(lfp_phase2{2}(control_idx), filt_phase2{2}(control_idx));
title('Fast Gamma Phase at Control times')
xlabel('Unfiltered Phase')
ylabel('Filtered LFP')

suptitle('Sanity Check')

figure
subplot(4,2,1)
hist(lfp_phase2{1}(on_idx), 5);
title('Unfiltered Theta Phase distribution at Stim times')
subplot(4,2,2)
hist(filt_phase2{1}(on_idx), 5);
title('Filtered Theta Phase distribution at Stim times')

subplot(4,2,3)
hist(lfp_phase2{1}(control_idx), 5);
title('Unfiltered Theta Phase distribution at Control times')

subplot(4,2,4)
hist(filt_phase2{1}(control_idx), 5);
title('Filtered Theta Phase distribution at Control times')

subplot(4,2,5)
hist(lfp_phase2{2}(on_idx), 5);
title('Unfiltered fast gamma-ish Phase distribution at Stim times')

subplot(4,2,6)
hist(filt_phase2{2}(on_idx), 5);
title('Filtered fast gamma-ish Phase distribution at Stim times')

subplot(4,2,7)
hist(lfp_phase2{2}(control_idx), 5);
title('Unfiltered fast gamma-ish Phase distribution at Control times')

subplot(4,2,8)
hist(filt_phase2{2}(control_idx), 5);
title('Filtered fast gamma-ish Phase distribution at Control times')
suptitle('Sanity Check')

%% Sanity FOR Eric's data
figure
on_idx = nearest_idx3(var_stim_on, this_lfp.tvec);
control_idx = nearest_idx3(control_on, this_lfp.tvec);

subplot(2,2,1)
hist(lfp_phase2{1}(on_idx), 5);
title('Theta Phase distribution at Stim times')

subplot(2,2,2)
hist(lfp_phase2{1}(control_idx), 5);
title('Theta Phase distribution at Control times')

subplot(2,2,3)
hist(lfp_phase2{2}(on_idx), 5);
title('Theta Gamma Phase distribution at Stim times')

subplot(2,2,4)
hist(lfp_phase2{2}(control_idx), 5);
title('Theta Gamma Phase distribution at Control times')

suptitle('Sanity Check')
