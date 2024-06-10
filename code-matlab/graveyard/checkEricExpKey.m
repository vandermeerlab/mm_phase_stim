%% Script to diagnose Eric Data 
clear;clc;close all;
folders = {'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-15_vStr_4p2_light_cells_TT6_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-16_dStr_3p3_light_cells_TT5_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-17_vStr_4p0_light_cells_TT7_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-18_vStr_4p0_light_cells_TT4_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-19_vStr_3p9_light_cells_TT4_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-20_vStr_4p3_light_cells_TT5_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-22_vStr_4p4_light_cells_TT8_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-23_dStr_2p5_light_cells_TT4_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-25_dStr_3p1_light_cells_TT1_TT3_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-26_dStr_3p3_light_cells_TT5_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M16\\M16-2019-02-27_vStr_3p9_light_cells_TT2_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-15_dStr_3p0_light_cells_TT3_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-16_dStr_3p2_light_cells_TT5_lost_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-16_dStr_3p6_light_cells_TT4_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-17_dStr_3p0_light_cells_TT8_TT6', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-18_dStr_3p7_light_cells_TT6_TT5_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-19_dStr_2p5_light_cells_TT5_short_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-20_vStr_3p7_light_cells_TT5_TT6_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-21_vStr_4p2_light_cells_TT7_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-24_dStr_3p6_light_cells_TT1_TT3', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-25_dStr_3p1_light_cells_TT5', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M17\\M17-2019-02-25_dStr_3p9_light_cells_TT1_TT3', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M18\\M18-2019-04-10-dStr_3p8_light_cells_TT4_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M18\\M18-2019-04-11_vStr_4p2_light_cells_TT7_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M18\\M18-2019-04-12_dStr_3p4_light_cells_TT4_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M18\\M18-2019-04-12_dStr_3p8_light_cells_TT7_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M18\\M18-2019-04-13_dStr_3p8_light_cells_TT6_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M18\\M18-2019-04-14_dStr_4p0_light_cells_TT3_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M18\\M18-2019-04-15_dStr_3p3_light_cells_TT8', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M18\\M18-2019-04-15_vStr_4p0_light_cells_TT7_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M19\\M19-2019-04-12_vStr_4p2_light_cells_TT7_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M19\\M19-2019-04-13-vStr_4p7_light_cells_TT2_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M19\\M19-2019-04-13_vStr_4p2_light_cells_TT3_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M19\\M19-2019-04-14_dStr_3p3_light_cells_TT8_min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M19\\M19-2019-04-14_vStr_4p2_light_cells_TT5 _min', ...
 'E:\\Dropbox (Dartmouth College)\\EC_State_inProcess\\M19\\M19-2019-04-15_dStr_4p0_light_cells_TT6_min', ...
 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M20\M20-2019-06-07_dStr_3p8_light_cells_TT6_TT8_min', ...
 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M20\M20-2019-06-07_dStr_4p6_light_cells_TT6_TT8_min', ...
 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M20\M20-2019-06-08_vStr_4p2_light_cells_TT7_min', ...
 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M20\M20-2019-06-09_vStr_4p7_light_cells_TT7_min', ...
 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M20\M20-2019-06-10_vStr_4p2_light_cells_TT5_TT6_TT7_min'};

% folders = {'E:\Dropbox (Dartmouth College)\EC_State_inProcess\M17\M17-2019-02-16_dStr_3p6_light_cells_TT4_min'};
% get rid of fixed values
fixed = false(size(folders));

for i = 1:length(folders)
    cd(folders{i})
    
    % Load ExpKeys and Events
    LoadExpKeys;
    evs = LoadEvents([]);
    stim_offset = 0.0011;
    
    % Load LFP
    cfg_lfp.fc = {ExpKeys.goodCSC}; % Eric Data
    this_lfp = LoadCSC(cfg_lfp);
    w = [-.1 .1]; % time window to compute STA over
    Fs = this_lfp.cfg.hdr{1}.SamplingFrequency;
    this_tvec = w(1):1/Fs:w(2); % time axis for STA

    fig = figure('WindowState','Maximized');
    % Plot for laser_on
    subplot(2,2,1)
    hold on;
    this_on = evs.t{strcmp(evs.label,ExpKeys.laser_on)};
    this_on = this_on((this_on >= ExpKeys.timeOnWheel) & (this_on <= ExpKeys.timeOffWheel));
    this_on  = this_on + stim_offset;
    for iEvt = 1:length(this_on) % for each stim ...
        on_sta_t = this_on(iEvt)+w(1);
        this_on_sta_idx = (nearest_idx3(on_sta_t,this_lfp.tvec));
        % grab LFP snippet for this window
        this_on_toAdd = this_lfp.data(this_on_sta_idx:this_on_sta_idx+ ...
            length(this_tvec)-1);
        this_on_sta(iEvt,:) = this_on_toAdd';
    end
    on_sta = mean(this_on_sta,1);
    plot(this_tvec, on_sta);
    xline(0);
    xlim([-0.1 0.1])
    title(sprintf("Laser ON: %d events\n",length(this_on)));
    clear this_on this_on_sta this_on_toAdd on_sta_t this_on_sta_idx

    subplot(2,2,2)
    hold on;
    this_on = evs.t{strcmp(evs.label,ExpKeys.var_trig)};
    this_on = this_on((this_on >= ExpKeys.timeOnWheel) & (this_on <= ExpKeys.timeOffWheel));
    this_on  = this_on + stim_offset;
    for iEvt = 1:length(this_on) % for each stim ...
        on_sta_t = this_on(iEvt)+w(1);
        this_on_sta_idx = (nearest_idx3(on_sta_t,this_lfp.tvec));
        % grab LFP snippet for this window
        this_on_toAdd = this_lfp.data(this_on_sta_idx:this_on_sta_idx+ ...
            length(this_tvec)-1);
        this_on_sta(iEvt,:) = this_on_toAdd';
    end
    on_sta = mean(this_on_sta,1);
    plot(this_tvec, on_sta);
    xline(0);
    xlim([-0.1 0.1])
    title(sprintf("Var Trig: %d events\n",length(this_on)));
    clear this_on this_on_sta this_on_toAdd on_sta_t this_on_sta_idx

    subplot(2,2,3)
    hold on;
    this_on = evs.t{strcmp(evs.label,ExpKeys.short_laser_on)};
    this_on = this_on((this_on <= ExpKeys.timeOnWheel) | (this_on >= ExpKeys.timeOffWheel));
    this_on  = this_on + stim_offset;
    for iEvt = 1:length(this_on) % for each stim ...
        on_sta_t = this_on(iEvt)+w(1);
        this_on_sta_idx = (nearest_idx3(on_sta_t,this_lfp.tvec));
        % grab LFP snippet for this window
        this_on_toAdd = this_lfp.data(this_on_sta_idx:this_on_sta_idx+ ...
            length(this_tvec)-1);
        this_on_sta(iEvt,:) = this_on_toAdd';
    end
    on_sta = mean(this_on_sta,1);
    plot(this_tvec, on_sta);
    xline(0);
    xlim([-0.1 0.1])
    title(sprintf("Hundred Stim: %d events\n",length(this_on)));
    clear this_on this_on_sta this_on_toAdd on_sta_t this_on_sta_idx

    subplot(2,2,4)
    hold on;
    if i ~= 13
        this_on = evs.t{strcmp(evs.label,ExpKeys.wide_laser_on)};
        this_on = this_on((this_on >= ExpKeys.timeOffWheel) & (this_on <= this_lfp.tvec(end)));
        this_on  = this_on + stim_offset;
        for iEvt = 1:length(this_on) % for each stim ...
            on_sta_t = this_on(iEvt)+w(1);
            this_on_sta_idx = (nearest_idx3(on_sta_t,this_lfp.tvec));
            % grab LFP snippet for this window
            this_on_toAdd = this_lfp.data(this_on_sta_idx:this_on_sta_idx+ ...
                length(this_tvec)-1);
            this_on_sta(iEvt,:) = this_on_toAdd';
        end
        on_sta = mean(this_on_sta,1);
        plot(this_tvec, on_sta);
        xline(0);
        xlim([-0.1 0.1])
        title(sprintf("Long Stim: %d events\n",length(this_on)));
    end
    clear this_on this_on_sta this_on_toAdd on_sta_t this_on_sta_idx

    WriteFig(fig, 'STA', 1);
    close;

%     fixed(i) =  contains(ExpKeys.stim_mode, 'fixed');
%     subplot(1,2,1)
%     hist(diff(evs.t{strcmp(evs.label,ExpKeys.laser_on)}))
%     title('Laser On')
%     subplot(1,2,2)
%     hist(diff(evs.t{strcmp(evs.label,ExpKeys.var_trig)}))
%     title('Var Trig')
%     sgtitle(strcat(ExpKeys.subject,'-',ExpKeys.date))
%     disp('Edit the ExpKey');
%     close all;
end