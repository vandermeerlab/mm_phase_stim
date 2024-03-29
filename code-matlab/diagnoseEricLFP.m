%% Script to diagnose Eric Data 

folders = {'E:\temp_phase_stim\ED\M16-2019-02-15_vStr_4p2_light_cells_TT6_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-16_dStr_3p3_light_cells_TT5_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-17_vStr_4p0_light_cells_TT7_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-18_vStr_4p0_light_cells_TT4_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-19_vStr_3p9_light_cells_TT4_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-20_vStr_4p3_light_cells_TT5_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-22_vStr_4p4_light_cells_TT8_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-23_dStr_2p5_light_cells_TT4_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-25_dStr_3p1_light_cells_TT1_TT3_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-26_dStr_3p3_light_cells_TT5_min', ...
'E:\temp_phase_stim\ED\M16-2019-02-27_vStr_3p9_light_cells_TT2_min', ...
'E:\temp_phase_stim\ED\M17-2019-02-15_dStr_3p0_light_cells_TT3_min', ...
'E:\temp_phase_stim\ED\M17-2019-02-16_dStr_3p2_light_cells_TT5_lost_min', ...
'E:\temp_phase_stim\ED\M17-2019-02-17_dStr_3p0_light_cells_TT8_TT6', ...
'E:\temp_phase_stim\ED\M17-2019-02-18_dStr_3p7_light_cells_TT6_TT5_min', ...
'E:\temp_phase_stim\ED\M17-2019-02-19_dStr_2p5_light_cells_TT5_short_min', ...
'E:\temp_phase_stim\ED\M17-2019-02-21_vStr_4p2_light_cells_TT7_min', ...
'E:\temp_phase_stim\ED\M17-2019-02-24_dStr_3p6_light_cells_TT1_TT3', ...
'E:\temp_phase_stim\ED\M17-2019-02-25_dStr_3p1_light_cells_TT5', ...
'E:\temp_phase_stim\ED\M17-2019-02-25_dStr_3p9_light_cells_TT1_TT3', ...
'E:\temp_phase_stim\ED\M18-2019-04-10-dStr_3p8_light_cells_TT4_min', ...
'E:\temp_phase_stim\ED\M18-2019-04-11_vStr_4p2_light_cells_TT7_min', ...
'E:\temp_phase_stim\ED\M18-2019-04-12_dStr_3p4_light_cells_TT4_min', ...
'E:\temp_phase_stim\ED\M18-2019-04-13_dStr_3p8_light_cells_TT6_min', ...
'E:\temp_phase_stim\ED\M18-2019-04-14_dStr_4p0_light_cells_TT3_min', ...
'E:\temp_phase_stim\ED\M18-2019-04-15_dStr_3p3_light_cells_TT8', ...
'E:\temp_phase_stim\ED\M18-2019-04-15_vStr_4p0_light_cells_TT7_min', ...
'E:\temp_phase_stim\ED\M19-2019-04-12_vStr_4p2_light_cells_TT7_min', ...
'E:\temp_phase_stim\ED\M19-2019-04-13-vStr_4p7_light_cells_TT2_min', ...
'E:\temp_phase_stim\ED\M19-2019-04-14_dStr_3p3_light_cells_TT8_min', ...
'E:\temp_phase_stim\ED\M19-2019-04-14_vStr_4p2_light_cells_TT5 _min', ...
'E:\temp_phase_stim\ED\M19-2019-04-15_dStr_4p0_light_cells_TT6_min'};

for iD = 1%23:length(folders)
    doStuff(folders{iD});
end
% folder = folders{iD}; % Choose what session you want
% cd(folder);
% LoadExpKeys;
% evs = LoadEvents([]);


%%
function doStuff(cur_dir)
    cd(cur_dir);
    LoadExpKeys;
    evs = LoadEvents([]);
    stim_offset = 0.0011;
    % ExpKeys.laser_on looks flat, meaning the var_stim trigger is the correct
    % thing to look at in variable ISI sessions
    var_stim_on = evs.t{strcmp(evs.label,ExpKeys.laser_on)} + stim_offset;
    ISIs = [diff(var_stim_on);ExpKeys.ISI];
    control_offset = arrayfun(@(x) x*rand(), ISIs); %generate random delays
    control_on = var_stim_on + control_offset;

    %%  Load LFP
    cfg_lfp.fc = {ExpKeys.goodCSC}; % Eric Data

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


    % %% Load at the STA, a few random Snips and decide the artifact window for interpolation
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

    %% Genrerate some random split indices
    num_snips = 5;
    sample_idx = randi(length(var_stim_on), 1, num_snips);

    % %% Plot some random splits
    % figure;
    % subplot(num_snips+1,1,1);
    % hold on
    % plot(this_tvec, on_sta);
    % plot(this_tvec, control_sta);
    % xline(0);
    % xlim([-0.1 0.1])
    % legend({'STA', 'Control', 'Stim Time'},'Location','southeast');
    % title('STA');
    % 
    % for iS = 1:num_snips
    %     hold on
    %     subplot(num_snips+1,1,iS+1)
    %     plot(this_tvec, this_on_sta(sample_idx(iS),:));
    %     xline(0);
    %     xlim([-0.1 0.1])
    % end


%     %% Interpolate the signal
%     artifact_window = 0.01; %in seconds
%     int_lfp = interpolateArtifacts('spline', var_stim_on, artifact_window, this_lfp);

    %%  Calculate the the PSD of the original and the interpolated LFP
    wsize = Fs;
%     [P_OG, F] = pwelch(this_lfp.data, hanning(wsize), wsize/2, [], Fs);
%     [P_int,~] = pwelch(int_lfp.data, hanning(wsize), wsize/2, [], Fs);

    %% Plot the two PSDs

%     figure;
%     hold on
%     plot(F, 10*log10(P_OG));
%     plot(F, 10*log10(P_int));
%     xlim([0 120]);
%     legend({'Original Signal', 'Interpolated Signal'})
%     % There seems to be a slight difference in higher frequencies
%     WriteFig('Diag_PSD');
%     close;

%     %% Calculate STAs for the interpolated signal
%     for iEvt = 1:length(var_stim_on) % for each stim ...
%         int_on_sta_t = var_stim_on(iEvt)+w(1);
%         this_on_sta_idx = (nearest_idx3(int_on_sta_t,int_lfp.tvec));
%         % grab LFP snippet for this window
%         this_on_toAdd = int_lfp.data(this_on_sta_idx:this_on_sta_idx+ ...
%             length(this_tvec)-1);
%         int_this_on_sta(iEvt,:) = this_on_toAdd';
% 
%         int_control_on_sta_t = control_on(iEvt)+w(1);
%         control_on_sta_idx = (nearest_idx3(int_control_on_sta_t,int_lfp.tvec));
%         % grab LFP snippet for this window
%         control_on_toAdd = int_lfp.data(control_on_sta_idx:control_on_sta_idx+ ...
%             length(this_tvec)-1);
%         int_control_on_sta(iEvt,:) = control_on_toAdd';
%     end
%     int_on_sta = mean(int_this_on_sta,1);
%     int_control_sta = mean(int_control_on_sta,1);


%     %% Plot the STA and snippets of the original signal
%     fig = figure;
%     subplot(num_snips+1,2,1);
%     hold on
%     plot(this_tvec, on_sta);
%     plot(this_tvec, control_sta);
%     xline(0);
%     xlim([-0.1 0.1])
%     legend({'STA', 'Control', 'Stim Time'});
%     title('STA for original signal');
% 
%     subplot(num_snips+1,2,2);
%     hold on
%     plot(this_tvec, int_on_sta);
%     plot(this_tvec, int_control_sta);
%     xline(0);
%     xlim([-0.1 0.1])
%     legend({'STA', 'Control', 'Stim Time'});
%     title('STA for interpolated signal');
%     for iS = 1:5
%         hold on
%         subplot(num_snips+1,2,(iS*2)+1)
%         plot(this_tvec, this_on_sta(sample_idx(iS),:));
%         xline(0);
%         xlim([-0.1 0.1])
% 
%         hold on
%         subplot(num_snips+1,2,(iS*2)+2)
%         plot(this_tvec, int_this_on_sta(sample_idx(iS),:));
%         xline(0);
%         xlim([-0.1 0.1])
%     end
% 
%     WriteFig('Diag_STA');
    close;

    %% Filter the original and the interpolated signals in 2 actual bands and 2 sanity check bands and calculate the hilbert phases
    f_list = {[2 5], [6 10], [25 55], [65 90]};
    filt_lfp = cell(length(f_list),2);
    filt_phase = cell(length(f_list),2);
    for iB = 1:length(f_list)
        cfg_filt.type = 'fdesign'; 
        cfg_filt.f  = f_list{iB};
        filt_lfp{iB,1} = FilterLFP(cfg_filt, this_lfp);
%         filt_lfp{iB,2} = FilterLFP(cfg_filt, int_lfp);
        filt_phase{iB,1} = angle(hilbert(filt_lfp{iB,1}.data));
%         filt_phase{iB,2} = angle(hilbert(filt_lfp{iB,2}.data));
    end

%     %% Plot distribution of hilbert phases of orignal and interpolated signals at stim_times and control_times
%     figure;
%     on_idx = nearest_idx3(var_stim_on, this_lfp.tvec);
%     control_idx = nearest_idx3(control_on, this_lfp.tvec);
% 
%     subplot(4,4,1)
%     hist(filt_phase{1,1}(on_idx), 5);
%     title('OG Delta Phases at Stim times')
% 
%     subplot(4,4,2)
%     hist(filt_phase{1,2}(on_idx), 5);
%     title('Interpolated Delta Phases at Stim times')
% 
%     subplot(4,4,3)
%     hist(filt_phase{1,1}(control_idx), 5);
%     title('OG Delta Phases at Control times')
% 
%     subplot(4,4,4)
%     hist(filt_phase{1,2}(control_idx), 5);
%     title('Interpolated Delta Phases at Control times')
% 
%     subplot(4,4,5)
%     hist(filt_phase{2,1}(on_idx), 5);
%     title('OG Theta Phases at Stim times')
% 
%     subplot(4,4,6)
%     hist(filt_phase{2,2}(on_idx), 5);
%     title('Interpolated Theta Phases at Stim times')
% 
%     subplot(4,4,7)
%     hist(filt_phase{2,1}(control_idx), 5);
%     title('OG Theta Phases at Control times')
% 
%     subplot(4,4,8)
%     hist(filt_phase{2,2}(control_idx), 5);
%     title('Interpolated Theta Phases at Control times')
% 
%     subplot(4,4,9)
%     hist(filt_phase{3,1}(on_idx), 5);
%     title('OG slow Gamma-ish Phases at Stim times')
% 
%     subplot(4,4,10)
%     hist(filt_phase{3,2}(on_idx), 5);
%     title('Interpolated Slow Gamma-ish Phases at Stim times')
% 
%     subplot(4,4,11)
%     hist(filt_phase{3,1}(control_idx), 5);
%     title('OG Slow Gamma-ish at Control times')
% 
%     subplot(4,4,12)
%     hist(filt_phase{3,2}(control_idx), 5);
%     title('Interpolated Slow Gamma-ish Phases at Control times')
% 
%     subplot(4,4,13)
%     hist(filt_phase{4,1}(on_idx), 5);
%     title('OG Fast Gamma-ish Phases at Stim times')
% 
%     subplot(4,4,14)
%     hist(filt_phase{4,2}(on_idx), 5);
%     title('Interpolated Fast Gamma-ish Phases at Stim times')
% 
%     subplot(4,4,15)
%     hist(filt_phase{4,1}(control_idx), 5);
%     title('OG Fast Gamma-ish at Control times')
% 
%     subplot(4,4,16)
%     hist(filt_phase{4,2}(control_idx), 5);
%     title('Interpolated Fast Gamma-ish Phases at Control times')
% 
%     WriteFig('Diag_Phase_Distribution');
%     close;
% 
%     %% Uncomment and run if you want to change number and/or identity of snips
%     num_snips = 5;
%     sample_idx = randi(length(var_stim_on), 1, num_snips);
% 
%     %% Extract snippets of filtered signal and their hilbert phase near stim_times
% 
%     filt_snips = cell(num_snips,2);
%     filt_snip_phases = cell(num_snips,2);
%     filt_stim_start = nearest_idx3(var_stim_on(sample_idx) + w(1), this_lfp.tvec);
%     filt_control_start = nearest_idx3(control_on(sample_idx) + w(1), this_lfp.tvec);
% 
%     %% Add desctiption to how the data is organized in these cells
%     for iS = 1:num_snips
%         filt_snips{iS,1} = cellfun(@(x) [x.data(filt_stim_start(iS):filt_stim_start(iS) + length(this_tvec) - 1)], filt_lfp, 'UniformOutput', false);
%         filt_snip_phases{iS,1} = cellfun(@(x) [x(filt_stim_start(iS):filt_stim_start(iS) +  length(this_tvec) - 1)], filt_phase, 'UniformOutput', false);
%         filt_snips{iS,2} = cellfun(@(x) [x.data(filt_control_start(iS):filt_control_start(iS) + length(this_tvec) - 1)], filt_lfp, 'UniformOutput', false);
%         filt_snip_phases{iS,2} = cellfun(@(x) [x(filt_control_start(iS):filt_control_start(iS) +  length(this_tvec) - 1)], filt_phase, 'UniformOutput', false);
%     end

%     %% Plot a figure for each Snippet 
%     for iS = 1:num_snips
%         figure
% 
%         % Plot Delta stuff
% 
%         % Plot stim stuff
%         subplot(4,4,1)
%         hold on
%         plot(this_tvec, filt_snips{iS,1}{1,1}, 'blue');
%         plot(this_tvec, filt_snips{iS,1}{1,2}, 'red');
%         plot(this_tvec, this_on_sta(sample_idx(iS),:), 'green');
%         xline(0, 'black')
%         title('Delta Signal at stim')
% 
%         subplot(4,4,5)
%         hold on
%         plot(this_tvec, filt_snip_phases{iS,1}{1,1}, 'blue');
%         plot(this_tvec, filt_snip_phases{iS,1}{1,2}, 'red');
%         xline(0, 'black')
%         title('Delta Phase at stim')
% 
%         % Plot control stuff
%         subplot(4,4,9)
%         hold on
%         plot(this_tvec, filt_snips{iS,2}{1,1}, 'blue');
%         plot(this_tvec, filt_snips{iS,2}{1,2}, 'red');
%         plot(this_tvec, control_on_sta(sample_idx(iS),:), 'green');
%         xline(0, 'black')
%         title('Delta Signal at control')
% 
%         subplot(4,4,13)
%         hold on
%         plot(this_tvec, filt_snip_phases{iS,2}{1,1}, 'blue');
%         plot(this_tvec, filt_snip_phases{iS,2}{1,2}, 'red');
%         xline(0, 'black')
%         title('Delta Phase at control')
% 
%         % Plot Theta stuff
% 
%         % Plot stim stuff
%         subplot(4,4,2)
%         hold on
%         plot(this_tvec, filt_snips{iS,1}{2,1}, 'blue');
%         plot(this_tvec, filt_snips{iS,1}{2,2}, 'red');
%         plot(this_tvec, this_on_sta(sample_idx(iS),:), 'green');
%         xline(0, 'black')
%         title('Theta Signal at stim')
% 
%         subplot(4,4,6)
%         hold on
%         plot(this_tvec, filt_snip_phases{iS,1}{2,1}, 'blue');
%         plot(this_tvec, filt_snip_phases{iS,1}{2,2}, 'red');
%         xline(0, 'black')
%         title('Theta Phase at stim')
% 
%         % Plot control stuff
%         subplot(4,4,10)
%         hold on
%         plot(this_tvec, filt_snips{iS,2}{2,1}, 'blue');
%         plot(this_tvec, filt_snips{iS,2}{2,2}, 'red');
%         plot(this_tvec, control_on_sta(sample_idx(iS),:), 'green');
%         xline(0, 'black')
%         title('Theta Signal at control')
% 
%         subplot(4,4,14)
%         hold on
%         plot(this_tvec, filt_snip_phases{iS,2}{2,1}, 'blue');
%         plot(this_tvec, filt_snip_phases{iS,2}{2,2}, 'red');
%         xline(0, 'black')
%         title('Theta Phase at control')
% 
%         % Plot Slow Gamma-ish stuff
% 
%         % Plot stim stuff
%         subplot(4,4,3)
%         hold on
%         plot(this_tvec, filt_snips{iS,1}{3,1}, 'blue');
%         plot(this_tvec, filt_snips{iS,1}{3,2}, 'red');
%         plot(this_tvec, this_on_sta(sample_idx(iS),:), 'green');
%         xline(0, 'black')
%         title('Slow Gamma-ish Signal at stim')
% 
%         subplot(4,4,7)
%         hold on
%         plot(this_tvec, filt_snip_phases{iS,1}{3,1}, 'blue');
%         plot(this_tvec, filt_snip_phases{iS,1}{3,2}, 'red');
%         xline(0, 'black')
%         title('Slow Gamma-ish Phase at stim')
% 
%         % Plot control stuff
%         subplot(4,4,11)
%         hold on
%         plot(this_tvec, filt_snips{iS,2}{3,1}, 'blue');
%         plot(this_tvec, filt_snips{iS,2}{3,2}, 'red');
%         plot(this_tvec, control_on_sta(sample_idx(iS),:), 'green');
%         xline(0, 'black')
%         title('Slow Gamma-ish Signal at control')
% 
%         subplot(4,4,15)
%         hold on
%         plot(this_tvec, filt_snip_phases{iS,2}{3,1}, 'blue');
%         plot(this_tvec, filt_snip_phases{iS,2}{3,2}, 'red');
%         xline(0, 'black')
%         title('Slow Gamma-ish Phase at control')
% 
%         % Plot Fast Gamma-ish stuff
% 
%         % Plot stim stuff
%         subplot(4,4,4)
%         hold on
%         plot(this_tvec, filt_snips{iS,1}{4,1}, 'blue');
%         plot(this_tvec, filt_snips{iS,1}{4,2}, 'red');
%         plot(this_tvec, this_on_sta(sample_idx(iS),:), 'green');
%         xline(0, 'black')
%         title('Fast Gamma-ish Signal at stim')
% 
%         subplot(4,4,8)
%         hold on
%         plot(this_tvec, filt_snip_phases{iS,1}{4,1}, 'blue');
%         plot(this_tvec, filt_snip_phases{iS,1}{4,2}, 'red');
%         xline(0, 'black')
%         title('Fast Gamma-ish Phase at stim')
% 
%         % Plot control stuff
%         subplot(4,4,12)
%         hold on
%         plot(this_tvec, filt_snips{iS,2}{4,1}, 'blue');
%         plot(this_tvec, filt_snips{iS,2}{4,2}, 'red');
%         plot(this_tvec, control_on_sta(sample_idx(iS),:), 'green');
%         xline(0, 'black')
%         title('Fast Gamma-ish Signal at control')
% 
%         subplot(4,4,16)
%         hold on
%         plot(this_tvec, filt_snip_phases{iS,2}{4,1}, 'blue');
%         plot(this_tvec, filt_snip_phases{iS,2}{4,2}, 'red');
%         xline(0, 'black')
%         title('Fast Gamma-ish Phase at control')
% 
%         WriteFig(cat(2,'Diag_Snippet_', num2str(iS)));
%         close;
% 
%     end
    dummy = 1;
end
%%
