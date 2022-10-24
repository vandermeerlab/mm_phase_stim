%% Script to show the problem

% Example Eric Data
% 'E:\Dropbox (Dartmouth
% College)\EC_State_inProcess\M18\M18-2019-04-12_dStr_3p4_light_cells_TT4_min'% Laser
% 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M18\M18-2019-04-13_dStr_3p8_light_cells_TT6_min' % Laser
% 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M19\M19-2019-04-14_vStr_4p2_light_cells_TT5 _min' % Laser

% Example Manish Data
% 'E:\Dropbox (Dartmouth College)\manish_data\M235\M235-2021-07-16' % Laser
% 'E:\Dropbox (Dartmouth College)\manish_data\M325\M325-2022-08-03' % LED

folder = 'E:\Dropbox (Dartmouth College)\manish_data\M235\M235-2021-07-16';
cd(folder);
LoadExpKeys;
evs = LoadEvents([]);

% var_stim_on = evs.t{strcmp(evs.label,ExpKeys.laser_on)} +  0.0011; % Eric data
var_stim_on = evs.t{strcmp(evs.label,ExpKeys.trial_stim_on)} + 0.0011; % Manish data

control_on = var_stim_on + 0.33;

%% Load LFP and look at STA
% cfg_lfp.fc = ExpKeys.goodCSC; % Eric Data
cfg_lfp.fc = ExpKeys.goodLFP; % Manish Data
if contains(cfg_lfp.fc, '-')
    temp = split(cfg_lfp.fc1,'-');
    cfg_lfp.fc = {cat(2,temp{1},'.ncs')};
    this_lfp = LoadCSC(cfg_lfp);
    cfg_temp.fc = {cat(2,temp{2},'.ncs')};
    ref = LoadCSC(cfg_temp);
    this_lfp.data = this_lfp.data - ref.data;
    clear temp ref;
else
    this_lfp = LoadCSC(cfg_lfp);
end
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
%%
