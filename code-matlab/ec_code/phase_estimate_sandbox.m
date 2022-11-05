%% phase_estimate_sandbox.m
% some code for testing ways of estimating the phase of the last sample in
% a signal

%% hilbert phase is accurate up to last sample of perfect sine wave
dt = 0.001;
tvec = 0:dt:2;

y = sin(2*pi*4*tvec);

plot(tvec, y); hold on;

phi = angle(hilbert(y));
plot(tvec, phi ./ pi, 'r');

%% what about some actual data? (note requires statedep_sandbox to have been run)
offs = 300; % plot some example stim events starting at this event number
for iE = 1:9

temp = restrict(this_csc, laser_on.t{1}(iE+offs)-2.5, laser_on.t{1}(iE+offs)+0.5);
this_idx = nearest_idx3(laser_on.t{1}(iE+offs), temp.tvec);

subplot(3, 3, iE);
plot(temp); hold on;
plot(temp.tvec(this_idx), temp.data(this_idx), '.k', 'MarkerSize', 10);
axis tight; box off;

end

%% test with some synthetic data
cfg_filt = []; cfg_filt.order = 4; cfg_filt.f = [6 10]; cfg_filt.pad = 0;
d = fdesign.bandpass('N,F3dB1,F3dB2', cfg_filt.order, cfg_filt.f(1), cfg_filt.f(2), Fs);
Hd = design(d,'butter');
b = Hd.sosMatrix; a = Hd.scaleValues;

% generate some data
random_csc = this_csc; random_csc.data = randn(size(random_csc.data)); % random
%random_csc = this_csc; random_csc.data = sin(2*pi*8*random_csc.tvec)'; % pure oscillation

random_csc_f = random_csc; 
random_csc_f.data = filter(Hd, random_csc_f.data);
random_csc_f.data = angle(hilbert(random_csc_f.data));

% filter trialized data
trial_hist = zeros(size(hist(random_csc_f.data, 36)));
clear first_phases last_phases first_values last_values;
for iT = length(laser_on.t{1}):-1:1
    
    this_trial = restrict(random_csc, laser_on.t{1}(iT)-2.5, laser_on.t{1}(iT));
    
    raw_first_values(iT) = this_trial.data(1);
    raw_last_values(iT) = this_trial.data(end);
    
    if cfg_filt.pad % add some padding at start and end
        
        x = this_trial.data;

        Np = 500;

        %x1 = -flipud(x(2:Np + 1)) + 2*x(1);
        %x2 = -flipud(x(end - Np:end-1)) + 2*x(end);
        x1 = fliplr(diff(x(1:Np)));
        x1 = cumsum(x1);
        x1 = x1 + (x(1) - x1(end));
        
        x2 = fliplr(diff(x(end - Np + 1:end)));
        x2 = cumsum(cat(2, x(end), x2));
        
        this_trial.data = cat(2, x1(1:end-1), x, x2(2:end));
        
    end
    
    this_trial.data = filter(Hd, this_trial.data);
    f_first_values(iT) = this_trial.data(1);
    f_last_values(iT) = this_trial.data(end);
    
    this_trial.data = angle(hilbert(this_trial.data));
    
    if cfg_filt.pad
        this_trial.data = this_trial.data(:, Np + 1:end - Np);
    end
    
    trial_hist = trial_hist + hist(this_trial.data, 36);
    first_phases(iT) = this_trial.data(1);
    last_phases(iT) = this_trial.data(end);

end

subplot(221);
hist(random_csc_f.data, 36);
title('full data phase hist');

subplot(222);
bar(trial_hist);
title('trialized data phase hist');

subplot(223);
hist(first_phases, 36);
title('start phase hist');

subplot(224);
hist(last_phases, 36);
title('end phase hist');

figure;
subplot(221);
hist(raw_first_values, 36);
title('raw first values');

subplot(222);
hist(raw_last_values, 36);
title('raw last values');

subplot(223);
hist(f_first_values, 36);
title('filtered first values');

subplot(224);
hist(f_last_values, 36);
title('filtered last values');

%% idea: run forward filter up to time of stim to obtain phase estimate (works on actual data)
fs = 18;
fpass_list = {[3 5], [7 9], [30 40], [65 80]};
fstop_list = {[2.5 5.5], [6 10], [28 42], [60 85]};

for iF = 1:length(fpass_list) % loop across freqs
   
    % set up filter
    cfg_filt = [];
    cfg_filt.fpass = fpass_list{iF};
    cfg_filt.fstop = fstop_list{iF};
    cfg_filt.debug = 0;
    cfg_filt.filtfilt = 0;
    cfg_filt.pad = [];
    
    stim_phase = FindPreStimPhase(cfg_filt, laser_on, this_csc);
    
    % STIM PHASE HISTO THIS IS IMPORTANT
    figure(2); subplot(2, 2, iF);
    hist(stim_phase, 36); title(sprintf('stim phase histo (%.1f-%.1f Hz)', fpass_list{iF}(1), fpass_list{iF}(2)));
    
    figure(1)
    subplot(322); 
    plot(F, 10*log10(Pxx), 'k', 'LineWidth', 2);
    set(gca, 'XLim', [0 150], 'FontSize', fs); grid on;
    xlabel('Frequency (Hz)');
    
    for iP = 1 % loop across some phase splits
    
        phase_low_idx = find(stim_phase < 0);
        phase_high_idx = find(stim_phase >= 0);
        
        [this_ccf_low, tvec] = ccf(cfg, laser_on.t{1}(phase_low_idx), this_S.t{1});
        [this_ccf_high, tvec] = ccf(cfg, laser_on.t{1}(phase_high_idx), this_S.t{1});
        
        subplot(3, 2, 2 + iF);
        h(1) = plot(tvec, this_ccf_low, 'b', 'LineWidth', 2); hold on;
        h(2) = plot(tvec, this_ccf_high, 'r', 'LineWidth', 2);
        
        legend(h, {'phase < 0', 'phase >= 0'}, 'Location', 'Northwest'); legend boxoff;
        set(gca, 'FontSize', fs); xlabel('time (s)'); ylabel('spike count');
        title(sprintf('phase split %.1f-%.1f Hz', fpass_list{iF}(1), fpass_list{iF}(2)));
    
    end
    
    drawnow;
    
end % of freq loop