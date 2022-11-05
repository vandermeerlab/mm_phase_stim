function [sdf] = MS_get_gau_sdf(cfg_in, tvec_in, data_in)
%% MS_get_gau_sdf: computes the spike density function (SDF) from time stamp data (data_in).
%
%
%  Inputs
%
%   - cfg_in [struct]  contains congifuration parameters
%
%
%
%
%
%  EC 2021-06-04  WIP


%% set up config
cfg_def = [];
cfg_def.dt = 0.001;
cfg_def.gauss_window = 1;
cfg_def.gauss_sd = 0.02;

cfg = ProcessConfig2(cfg_def, cfg_in);


%% compute the firing rate across tvec_in

tbin_edges = tvec_in(1):cfg.dt:tvec_in(end);
spk_count = histc(data_in,tbin_edges);
spk_count = spk_count(1:end-1);


% set up gau kernal
gauss_window = cfg.gauss_window./cfg.dt; % 1 second window
gauss_SD = cfg.gauss_sd./cfg.dt; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./cfg.dt; % normalize by binsize

sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
