%%
cd('data/M295/M295-2022-01-06');
LoadExpKeys;
evs = LoadEvents([]);
cfg_spk = [];
cfg_spk.fc = {'M295-2022-01-06-TT06_4.t'};
if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg_spk.getRatings = 1;
        cfg_spk.uint = '64';
end
S = LoadSpikes(cfg_spk);
%%
if contains(ExpKeys.light_source, 'LASER')
    start_delay = 0.0011;
    stop_delay = 0.0012; %0.0022 is the spike width in the case of 1 msec laser pulse
else
    start_delay = 0;
    stop_delay = 0;
end
stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
    stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                      stim_on <= ExpKeys.stim_times(2));
else
    stim_on = [];
end

% Load CSC
cfg = []; cfg.fc = ExpKeys.goodLFP;
if contains(cfg.fc, '-')
    temp = split(cfg.fc,'-');
    cfg.fc = {cat(2,temp{1},'.ncs')};
    csc = LoadCSC(cfg);
    cfg_temp.fc = {cat(2,temp{2},'.ncs')};
    ref = LoadCSC(cfg_temp);
    csc.data = csc.data - ref.data;    
    clear temp ref;
else
    csc = LoadCSC(cfg);
end

% Sometimes csc.tvec can be have weird elements because of gaps
% in recording
to_remove = find(diff(csc.tvec)<=0);
while (~isempty(to_remove))
    csc.tvec(to_remove+1) = [];
    csc.data(to_remove+1) = [];
    to_remove = find(diff(csc.tvec)<=0);
end

% restrict csc
csc = restrict(csc, iv(ExpKeys.stim_times));
%%
% decimate csc
cfg_d = [];
cfg_d.decimateFactor = 16;
csc2 = decimate_tsd(cfg_d, csc);

%%
% Parameters for MUA
cfg_MUA = []; 
cfg_MUA.tvec = csc2.tvec; % timebase to compute MUA on
cfg_MUA.sigma = 0.001;

% Parameters for PETH
cfg_peth = []; % parameters for PETH
cfg_peth.window = [-0.2 0.2];
cfg_peth.dt = 0.001;
cfg_peth.mode = 'interp';

this_MUA = getMUA(cfg_MUA, S); % "MUA" for one cell is just that cell's firing rate
this_MUAz = zscore_tsd(this_MUA);
this_stim_peth = TSDpeth_fast(cfg_peth, this_MUA, ...
    stim_on(ExpKeys.goodTrials(1,1):ExpKeys.goodTrials(1,2)));
this_stim_zpeth = TSDpeth_fast(cfg_peth, this_MUAz, ...
    stim_on(ExpKeys.goodTrials(1,1):ExpKeys.goodTrials(1,2)));
%%
figure;
subplot(1,2,1)
plot (this_stim_peth)
title('Raw Peth')

subplot(1,2,2)
plot (this_stim_zpeth)
title('Z-scored Peth')

%%
fprintf('Sampling rate = %.2f \nMean of this_MUA.data = %.2f, Std of this_MUA.data = %.2f\nMax of this_MUAz.data = %.2f, Min of this_MUAz.data = %.2f\n', ...
    1/median(diff(csc.tvec)), mean(this_MUA.data), ...
    std(this_MUA.data), max(this_MUAz.data), min(this_MUAz.data));
%%
figure
ax1 = subplot(2,1,1);
plot(ax1,this_MUA.tvec, this_MUA.data)
ax1.Title.String = 'Raw MUA';
ax2 = subplot(2,1,2);
plot(ax2,this_MUAz.tvec, this_MUAz.data)
ax2.Title.String = 'Z-scored MUA';
sgtitle(sprintf('Bin-size = %.3f milliseconds',1000*median(diff(csc2.tvec))));
% linkaxes([ax1, ax2], 'x')


%%
% Parameters for gaussian kernel
sigma = 0.025;  % Standard deviation in seconds
timebase = median(diff(csc2.tvec));  % Your timebase in seconds
window = 0.025;  % Width of kernel in seconds (e.g., 4*sigma)

% Create gaussian kernel
kernel_t = -window:timebase:window;  % Kernel time vector
kernel = exp(-(kernel_t.^2)/(2*sigma^2));
kernel = kernel/sum(kernel);  % Normalize to sum to 1

norm_MUA = this_MUA;
norm_MUA.data = conv(norm_MUA.data,kernel,'same');

norm_MUAz = zscore_tsd(norm_MUA);

% Parameters for PETH
cfg_peth = []; % parameters for PETH
cfg_peth.window = [-0.2 0.2];
cfg_peth.dt = 0.001;
cfg_peth.mode = 'interp';

norm_zpeth = TSDpeth_fast(cfg_peth, norm_MUAz, ...
    stim_on(ExpKeys.goodTrials(1,1):ExpKeys.goodTrials(1,2)));

%
figure;
ax1 = subplot(3,2,1);
plot(ax1,this_MUA.tvec, this_MUA.data)
ax1.Title.String = 'Raw MUA';
ax2 = subplot(3,2,3);
plot(ax2,norm_MUA.tvec, norm_MUA.data)
ax2.Title.String = sprintf('smoothedMUA, gaussian sigma = %.4f sec, kernel width = %.4f sec', ...
    sigma, window);
ax3 = subplot(3,2,5);
plot(ax3,norm_MUAz.tvec, norm_MUAz.data)
ax3.Title.String = 'z-scored smoothedMUA';
ax4 = subplot(3,2,[2,4,6]);
plot(ax4, norm_zpeth.tvec, norm_zpeth.data);
ax4.Title.String = 'Peth for z-scored smoothedMUA';