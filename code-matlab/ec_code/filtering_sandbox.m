%% filtering_sandbox.m
%
% code for comparing true and estimated phases of last sample in various
% signals

%% generate signal and set up filter
nP = 1e4;
x = randn(nP, 1); % random signal
Fs = 1e3; dt = 1 ./ Fs;

tvec = 0:dt:(nP*dt)-dt; % or, sine wave plus white gaussian noise
x = sin(2*pi*8.1*tvec'); x = x + 0.3*rand(size(x)); % freq should not be a divisor of cut points

f = [2.5 5.5];

d = fdesign.bandpass('N,F3dB1,F3dB2', 4, f(1), f(2), Fs);
flt = design(d,'butter');

%% or, use an actual neural signal
cfg = [];
cfg.decimateByFactor = 30;
this_csc = LoadCSC(cfg);
Fs = 1 ./ median(diff(this_csc.tvec));
%%
x = this_csc.data(40001:80000)';
%% filter whole thing and get hilbertized version
xf = filter(flt, x);
xa = angle(hilbert(xf)); % uniform phase distribution
[~,~,~,~,xp] = InstFreq_v2(xf', 1./Fs); % non-uniform phase distribution...
xp = wrapToPi(xp-pi/2)'; % make phases consistent
%% loop over different cut points
% set up set of cut points (for plotting against length later)
cp = 1e3:1e3/10:4e4;

for iP = 1:length(cp)
    % filter up to cut point, keep track of filtered and hilbertized points
    this_x = x(1:cp(iP));
    
    this_xf = filter(flt, this_x);
    this_xa = angle(hilbert(this_xf));
    [~,~,~,~,this_xp] = InstFreq_v2(this_xf', 1./Fs);

    ALL_raw(iP) = this_x(end);
    ALL_f(iP) = this_xf(end);
    ALL_a(iP) = this_xa(end);
    ALL_p(iP) = this_xp(end-1); % note last sample phase is not defined with this method
    
end

%% plot
figure;
subplot(221);
plot(x(cp), ALL_raw, '.'); refline(1,0); axis square;
title('raw data');

subplot(222);
plot(xf(cp), ALL_f, '.'); refline(1,0);  axis square;
title(sprintf('filtered %.1f-%.1f Hz', f(1), f(2)));

subplot(223);
plot(xa(cp), ALL_a, '.'); refline(1,0);  axis square;
title(sprintf('hilbert angle', f(1), f(2)));

subplot(224);
plot(xp(cp-1), ALL_p, '.'); refline(1,0);  axis square;
title(sprintf('realtime angle', f(1), f(2)));


%% expectations:
% *filtering cut vs. whole thing shouldn't matter, because filter is causal
% (only uses points in the past, not the future)
% *not sure what hilbert transform would do
% *could use some other phase estimation method like linear interpolation
% between peaks (predict next peak based on last peak, running average,
% etc...)