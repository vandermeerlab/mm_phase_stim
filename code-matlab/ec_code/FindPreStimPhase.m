function stim_phase = FindPreStimPhase(cfg_in, stim_times, csc_in)

cfg_def = [];
cfg_def.skip_time = 0.5; % time (in s) to start next trial
cfg_def.isi = 3; % inter-stim interval
cfg_def.fpass = [7 9];
cfg_def.fstop = [6 10];
cfg_def.att = 30; % stopband minimum attenuation
cfg_def.ripple = 2; % passband ripple
cfg_def.filtfilt = 0;
cfg_def.debug = 0;
cfg_def.pad = []; % pad trials before filtering; if non-empty, use this number of points on either end of signal
cfg_def.order = 4;
cfg_def.method = 'cicone'; % phase estimation method: 'cicone' or 'hilbert'
% addtion
cfg_def.filter_in = []; % is a prebuilt filter is specifed.  
cfg_def.filter_dir = cd; % place to save built filters


cfg = ProcessConfig(cfg_def, cfg_in);

% make ivs
tstart = cat(2, stim_times.t{1}(1) - cfg.isi + cfg.skip_time, stim_times.t{1}(1:end-1) + cfg.skip_time);
tend = stim_times.t{1};
trials = iv(tstart, tend);

% set up filter
fprintf('Creating filter...\n');

% Fs = 1 ./ median(diff(csc_in.tvec));
Fs = csc_in.cfg.hdr{1}.SamplingFrequency;

if ~isfield(cfg, 'filter_in') || isempty(cfg.filter_in) 

d = fdesign.bandpass(cfg.fstop(1), cfg.fpass(1), cfg.fpass(2), cfg.fstop(2), cfg.att, cfg.ripple, cfg.att, Fs);
flt = design(d, 'equiripple'); % equiripple method (FIR) gives linear phase delay

save([cfg.filter_dir 'Filt_' num2str(round(cfg.fpass(1))) '_' num2str(round(cfg.fpass(2))) '_Fs_' num2str(Fs) '.mat'], 'flt', '-v7.3')
else
    flt = cfg.filter_in;

end
if cfg.debug
   fvtool(flt);
   pause;
end

% loop over ivs, filter & get final phase
for iT = length(tstart):-1:1
    
    this_trial = restrict(csc_in, tstart(iT), tend(iT));
    this_trial.tvec = this_trial.tvec(1:end-1);
    this_trial.data = this_trial.data(1:end-1);
    
    if ~isempty(cfg.pad)
        
        x = this_trial.data;

        Np = cfg.pad;

        x1 = fliplr(diff(x(1:Np)));
        x1 = cumsum(x1);
        x1 = x1 + (x(1) - x1(end));
        
        x2 = fliplr(diff(x(end - Np + 1:end)));
        x2 = cumsum(cat(2, x(end), x2));
        
        this_trial.data = cat(2, x1(1:end-1), x, x2(2:end));
        
    end
    
    switch cfg.filtfilt
        case 0
            trial_f = filter(flt, this_trial.data); % note, this will generate a phase shift
        case 1
            trial_f = filtfilt(flt.Numerator, 1, this_trial.data);
    end
    
    switch cfg.method
        case 'hilbert'
            trial_f = angle(hilbert(trial_f));
        case 'cicone'
            [~,~,~,~,trial_f] = InstFreq_v2(trial_f, 1./Fs);
            trial_f = wrapToPi(trial_f - pi/2); % make consistent with hilbert
            trial_f = trial_f(1:end-1); % last sample is NaN
    end
            
    if cfg.pad
        trial_f = trial_f(:, Np + 1:end - Np);
    end
    
    stim_phase(iT) = trial_f(end);
    
end