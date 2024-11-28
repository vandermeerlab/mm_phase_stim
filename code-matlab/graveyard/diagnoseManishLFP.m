%% Script to diagnose Manish Data 

% Example Manish Data
% 'data\M235\M235-2021-07-16' % Laser
% 'data\M325\M325-2022-08-03' % LED

folder = 'data\M320\M320-2022-05-30';
cd(folder);
LoadExpKeys;
evs = LoadEvents([]);
%%
stim_offset = 0;
if strcmp(ExpKeys.light_source, 'Shanghai Dream LASER')
    stim_offset = 0.0011;
end
var_stim_on = evs.t{strcmp(evs.label,ExpKeys.trial_stim_on)} + stim_offset;
ISIs = [diff(var_stim_on);ExpKeys.ISI];
control_offset = arrayfun(@(x) x*rand(), ISIs); %generate random delays
control_on = var_stim_on + control_offset;

%%  Load LFP
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


%% Load at the STA, a few random Snips and decide the artifact window for interpolation
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

%% Plot some random splits
figure;
subplot(num_snips+1,1,1);
hold on
plot(this_tvec, on_sta);
plot(this_tvec, control_sta);
xline(0);
xlim([-0.1 0.1])
legend({'STA', 'Control', 'Stim Time'},'Location','southeast');
title('STA');

for iS = 1:num_snips
    hold on
    subplot(num_snips+1,1,iS+1)
    plot(this_tvec, this_on_sta(sample_idx(iS),:));
    xline(0);
    xlim([-0.1 0.1])
end


%% Interpolate the signal
artifact_window = 0.0015; %in seconds
int_lfp = interpolateArtifacts('linear', var_stim_on, artifact_window, this_lfp);

%%  Calculate the the PSD of the original and the interpolated LFP
wsize = Fs;
[P_OG, F] = pwelch(this_lfp.data, hanning(wsize), wsize/2, [], Fs);
[P_int,~] = pwelch(int_lfp.data, hanning(wsize), wsize/2, [], Fs);

%% Plot the two PSDs
figure;
hold on
plot(F, 10*log10(P_OG));
plot(F, 10*log10(P_int));
xlim([0 120]);
legend({'Original Signal', 'Interpolated Signal'})
% There seems to be a slight difference in higher frequencies

%% Calculate STAs for the interpolated signal
for iEvt = 1:length(var_stim_on) % for each stim ...
    int_on_sta_t = var_stim_on(iEvt)+w(1);
    this_on_sta_idx = (nearest_idx3(int_on_sta_t,int_lfp.tvec));
    % grab LFP snippet for this window
    this_on_toAdd = int_lfp.data(this_on_sta_idx:this_on_sta_idx+ ...
        length(this_tvec)-1);
    int_this_on_sta(iEvt,:) = this_on_toAdd';

    int_control_on_sta_t = control_on(iEvt)+w(1);
    control_on_sta_idx = (nearest_idx3(int_control_on_sta_t,int_lfp.tvec));
    % grab LFP snippet for this window
    control_on_toAdd = int_lfp.data(control_on_sta_idx:control_on_sta_idx+ ...
        length(this_tvec)-1);
    int_control_on_sta(iEvt,:) = control_on_toAdd';
end
int_on_sta = mean(int_this_on_sta,1);
int_control_sta = mean(int_control_on_sta,1);


%% Plot the STA and snippets of the original signal
figure;
subplot(num_snips+1,2,1);
hold on
plot(this_tvec, on_sta);
plot(this_tvec, control_sta);
xline(0);
xlim([-0.1 0.1])
legend({'STA', 'Control', 'Stim Time'});
title('STA for original signal');

subplot(num_snips+1,2,2);
hold on
plot(this_tvec, int_on_sta);
plot(this_tvec, int_control_sta);
xline(0);
xlim([-0.1 0.1])
legend({'STA', 'Control', 'Stim Time'});
title('STA for interpolated signal');
for iS = 1:5
    hold on
    subplot(num_snips+1,2,(iS*2)+1)
    plot(this_tvec, this_on_sta(sample_idx(iS),:));
    xline(0);
    xlim([-0.1 0.1])
    
    hold on
    subplot(num_snips+1,2,(iS*2)+2)
    plot(this_tvec, int_this_on_sta(sample_idx(iS),:));
    xline(0);
    xlim([-0.1 0.1])
end

%% Filter the original and the interpolated signals in 2 actual bands and 2 sanity check bands and calculate the hilbert phases
f_list = {[2 5], [6 10], [25 55], [65 90]};
filt_lfp = cell(length(f_list),2);
filt_phase = cell(length(f_list),2);
for iB = 1:length(f_list)
    cfg_filt.type = 'fdesign'; 
    cfg_filt.f  = f_list{iB};
    filt_lfp{iB,1} = FilterLFP(cfg_filt, this_lfp);
    filt_lfp{iB,2} = FilterLFP(cfg_filt, int_lfp);
    filt_phase{iB,1} = angle(hilbert(filt_lfp{iB,1}.data));
    filt_phase{iB,2} = angle(hilbert(filt_lfp{iB,2}.data));
end

%% Plot distribution of hilbert phases of orignal and interpolated signals at stim_times and control_times
figure;
on_idx = nearest_idx3(var_stim_on, this_lfp.tvec);
control_idx = nearest_idx3(control_on, this_lfp.tvec);

subplot(4,4,1)
hist(filt_phase{1,1}(on_idx), 5);
title('OG Delta Phases at Stim times')

subplot(4,4,2)
hist(filt_phase{1,2}(on_idx), 5);
title('Interpolated Delta Phases at Stim times')

subplot(4,4,3)
hist(filt_phase{1,1}(control_idx), 5);
title('OG Delta Phases at Control times')

subplot(4,4,4)
hist(filt_phase{1,2}(control_idx), 5);
title('Interpolated Delta Phases at Control times')

subplot(4,4,5)
hist(filt_phase{2,1}(on_idx), 5);
title('OG Theta Phases at Stim times')

subplot(4,4,6)
hist(filt_phase{2,2}(on_idx), 5);
title('Interpolated Theta Phases at Stim times')

subplot(4,4,7)
hist(filt_phase{2,1}(control_idx), 5);
title('OG Theta Phases at Control times')

subplot(4,4,8)
hist(filt_phase{2,2}(control_idx), 5);
title('Interpolated Theta Phases at Control times')

subplot(4,4,9)
hist(filt_phase{3,1}(on_idx), 5);
title('OG slow Gamma-ish Phases at Stim times')

subplot(4,4,10)
hist(filt_phase{3,2}(on_idx), 5);
title('Interpolated Slow Gamma-ish Phases at Stim times')

subplot(4,4,11)
hist(filt_phase{3,1}(control_idx), 5);
title('OG Slow Gamma-ish at Control times')

subplot(4,4,12)
hist(filt_phase{3,2}(control_idx), 5);
title('Interpolated Slow Gamma-ish Phases at Control times')

subplot(4,4,13)
hist(filt_phase{4,1}(on_idx), 5);
title('OG Fast Gamma-ish Phases at Stim times')

subplot(4,4,14)
hist(filt_phase{4,2}(on_idx), 5);
title('Interpolated Fast Gamma-ish Phases at Stim times')

subplot(4,4,15)
hist(filt_phase{4,1}(control_idx), 5);
title('OG Fast Gamma-ish at Control times')

subplot(4,4,16)
hist(filt_phase{4,2}(control_idx), 5);
title('Interpolated Fast Gamma-ish Phases at Control times')


%% Uncomment and run if you want to change number and/or identity of snips
% num_snips = 5;
% sample_idx = randi(length(var_stim_on), 1, num_snips);

%% Extract snippets of filtered signal and their hilbert phase near stim_times

filt_snips = cell(num_snips,2);
filt_snip_phases = cell(num_snips,2);
filt_stim_start = nearest_idx3(var_stim_on(sample_idx) + w(1), this_lfp.tvec);
filt_control_start = nearest_idx3(control_on(sample_idx) + w(1), this_lfp.tvec);

%% Add desctiption to how the data is organized in these cells
for iS = 1:num_snips
    filt_snips{iS,1} = cellfun(@(x) [x.data(filt_stim_start(iS):filt_stim_start(iS) + length(this_tvec) - 1)], filt_lfp, 'UniformOutput', false);
    filt_snip_phases{iS,1} = cellfun(@(x) [x(filt_stim_start(iS):filt_stim_start(iS) +  length(this_tvec) - 1)], filt_phase, 'UniformOutput', false);
    filt_snips{iS,2} = cellfun(@(x) [x.data(filt_control_start(iS):filt_control_start(iS) + length(this_tvec) - 1)], filt_lfp, 'UniformOutput', false);
    filt_snip_phases{iS,2} = cellfun(@(x) [x(filt_control_start(iS):filt_control_start(iS) +  length(this_tvec) - 1)], filt_phase, 'UniformOutput', false);
end

%% Plot a figure for each Snippet 
for iS = 1:num_snips
    figure
    
    % Plot Delta stuff
    
    % Plot stim stuff
    subplot(4,4,1)
    hold on
    plot(this_tvec, filt_snips{iS,1}{1,1}, 'blue');
    plot(this_tvec, filt_snips{iS,1}{1,2}, 'red');
    plot(this_tvec, this_on_sta(sample_idx(iS),:), 'green');
    xline(0, 'black')
    title('Delta Signal at stim')

    subplot(4,4,5)
    hold on
    plot(this_tvec, filt_snip_phases{iS,1}{1,1}, 'blue');
    plot(this_tvec, filt_snip_phases{iS,1}{1,2}, 'red');
    xline(0, 'black')
    title('Delta Phase at stim')
    
    % Plot control stuff
    subplot(4,4,9)
    hold on
    plot(this_tvec, filt_snips{iS,2}{1,1}, 'blue');
    plot(this_tvec, filt_snips{iS,2}{1,2}, 'red');
    plot(this_tvec, control_on_sta(sample_idx(iS),:), 'green');
    xline(0, 'black')
    title('Delta Signal at control')
    
    subplot(4,4,13)
    hold on
    plot(this_tvec, filt_snip_phases{iS,2}{1,1}, 'blue');
    plot(this_tvec, filt_snip_phases{iS,2}{1,2}, 'red');
    xline(0, 'black')
    title('Delta Phase at control')
    
    % Plot Theta stuff
    
    % Plot stim stuff
    subplot(4,4,2)
    hold on
    plot(this_tvec, filt_snips{iS,1}{2,1}, 'blue');
    plot(this_tvec, filt_snips{iS,1}{2,2}, 'red');
    plot(this_tvec, this_on_sta(sample_idx(iS),:), 'green');
    xline(0, 'black')
    title('Theta Signal at stim')

    subplot(4,4,6)
    hold on
    plot(this_tvec, filt_snip_phases{iS,1}{2,1}, 'blue');
    plot(this_tvec, filt_snip_phases{iS,1}{2,2}, 'red');
    xline(0, 'black')
    title('Theta Phase at stim')
    
    % Plot control stuff
    subplot(4,4,10)
    hold on
    plot(this_tvec, filt_snips{iS,2}{2,1}, 'blue');
    plot(this_tvec, filt_snips{iS,2}{2,2}, 'red');
    plot(this_tvec, control_on_sta(sample_idx(iS),:), 'green');
    xline(0, 'black')
    title('Theta Signal at control')
    
    subplot(4,4,14)
    hold on
    plot(this_tvec, filt_snip_phases{iS,2}{2,1}, 'blue');
    plot(this_tvec, filt_snip_phases{iS,2}{2,2}, 'red');
    xline(0, 'black')
    title('Theta Phase at control')
    
    % Plot Slow Gamma-ish stuff
    
    % Plot stim stuff
    subplot(4,4,3)
    hold on
    plot(this_tvec, filt_snips{iS,1}{3,1}, 'blue');
    plot(this_tvec, filt_snips{iS,1}{3,2}, 'red');
    plot(this_tvec, this_on_sta(sample_idx(iS),:), 'green');
    xline(0, 'black')
    title('Slow Gamma-ish Signal at stim')

    subplot(4,4,7)
    hold on
    plot(this_tvec, filt_snip_phases{iS,1}{3,1}, 'blue');
    plot(this_tvec, filt_snip_phases{iS,1}{3,2}, 'red');
    xline(0, 'black')
    title('Slow Gamma-ish Phase at stim')
    
    % Plot control stuff
    subplot(4,4,11)
    hold on
    plot(this_tvec, filt_snips{iS,2}{3,1}, 'blue');
    plot(this_tvec, filt_snips{iS,2}{3,2}, 'red');
    plot(this_tvec, control_on_sta(sample_idx(iS),:), 'green');
    xline(0, 'black')
    title('Slow Gamma-ish Signal at control')
    
    subplot(4,4,15)
    hold on
    plot(this_tvec, filt_snip_phases{iS,2}{3,1}, 'blue');
    plot(this_tvec, filt_snip_phases{iS,2}{3,2}, 'red');
    xline(0, 'black')
    title('Slow Gamma-ish Phase at control')
    
    % Plot Fast Gamma-ish stuff
    
    % Plot stim stuff
    subplot(4,4,4)
    hold on
    plot(this_tvec, filt_snips{iS,1}{4,1}, 'blue');
    plot(this_tvec, filt_snips{iS,1}{4,2}, 'red');
    plot(this_tvec, this_on_sta(sample_idx(iS),:), 'green');
    xline(0, 'black')
    title('Fast Gamma-ish Signal at stim')

    subplot(4,4,8)
    hold on
    plot(this_tvec, filt_snip_phases{iS,1}{4,1}, 'blue');
    plot(this_tvec, filt_snip_phases{iS,1}{4,2}, 'red');
    xline(0, 'black')
    title('Fast Gamma-ish Phase at stim')
    
    % Plot control stuff
    subplot(4,4,12)
    hold on
    plot(this_tvec, filt_snips{iS,2}{4,1}, 'blue');
    plot(this_tvec, filt_snips{iS,2}{4,2}, 'red');
    plot(this_tvec, control_on_sta(sample_idx(iS),:), 'green');
    xline(0, 'black')
    title('Fast Gamma-ish Signal at control')
    
    subplot(4,4,16)
    hold on
    plot(this_tvec, filt_snip_phases{iS,2}{4,1}, 'blue');
    plot(this_tvec, filt_snip_phases{iS,2}{4,2}, 'red');
    xline(0, 'black')
    title('Fast Gamma-ish Phase at control')
end

