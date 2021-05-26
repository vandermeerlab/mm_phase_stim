% cd('D:\Dropbox (Dartmouth College)\vandermeerlab_tutorial\data\R016-2012-10-03');
cd('/Users/manishm/Dropbox (Dartmouth College)/vandermeerlab_tutorial/data/R016-2012-10-03');
% goodGamma: {'R016-2012-10-03-CSC04d.ncs'
% goodSWR: {'R016-2012-10-03-CSC02b.ncs'}
% goodTheta: {'R016-2012-10-03-CSC02b.ncs'}
cfg_in.fc = {'R016-2012-10-03-CSC02b.ncs'};
trial_lfp  = LoadCSC(cfg_in);
evs = LoadEvents([]);

%% Check PSD
Fs = trial_lfp.cfg.hdr{1}.SamplingFrequency;
wsize = Fs;
for i = 1:length(trial_lfp.label)
    figure;
    subplot(2,1,1)
    [Pxx, F] = pwelch(trial_lfp.data(i,:), rectwin(wsize), wsize/2, [], Fs);
    plot(F, 10*log10(Pxx));
    xlim([0,120]);
    title('PSD')
    subplot(2,1,2)
    [S,F,T,P] = spectrogram(trial_lfp.data(i,:),hanning(round(Fs*4)),round(Fs*2),1:120,Fs);
    imagesc(T,F,10*log10(P)); % converting to dB as usual
%     yline(60, 'LineWidth', 1);
    axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)'); 
    title('Spectrogram')
    suptitle(trial_lfp.label{i})
end


%% Estimate ground truth phases using FilterLFP + Hilbert transform

% make a cell array of the various processed LFPS
lfp_prefilt{1} = trial_lfp;

f_list = {[2 5], [6.5 9.5],[14 25] ,[40 65], [70 100]};
f_list_label = {'2 - 5', '6.5 - 9.5', '14 - 25', '40 - 65', '70 - 100'};
fstop_list = {[1.5 5.5], [6 10],[13 26],[38 67], [68 102]};

% Save filtered LFP and all_phase for later diagnosis
lfp_filtered = cell(length(f_list), length(lfp_prefilt));
all_phase = cell(length(f_list), length(lfp_prefilt));


for iF = 1:length(f_list)
     cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
     for iL = 1:length(lfp_prefilt)
         this_filt = FilterLFP(cfg_filt, lfp_prefilt{iL});
         lfp_filtered{iF,iL} = this_filt;
         this_res = zeros(size(this_filt.data));
         for iCSC = 1:length(cfg_in.fc)
            this_res(iCSC,:) = angle(hilbert(this_filt.data(iCSC,:)));
         end
         all_phase{iF,iL} = this_res;
     end
end


%% Create random trials

trial_data = cell(1,25);
trial_phase = cell(1,25);
trial_filt = cell(1,25);
seeds = randi(length(lfp_prefilt{1}.data),1,25);
for i = 1:25
    t2 = lfp_prefilt{1}.tvec(seeds(i));
    t1_idx = nearest_idx3(t2 - 2.5, lfp_prefilt{1}.tvec);
    % If the trial length is not even, make it so!
    if rem(seeds(i)-t1_idx, 2) == 0
        t1_idx = t1_idx-1;
    end
    t1 = lfp_prefilt{1}.tvec(t1_idx);
    trial_data{i} = restrict(lfp_prefilt{1}, iv(t1,t2));
    phase_length = length(trial_data{i}.tvec) - round(length(trial_data{i}.tvec)*0.98);
    trial_filt{i} = zeros(length(trial_lfp.label),length(all_phase),phase_length);
    trial_phase{i} = zeros(length(trial_lfp.label),length(all_phase),phase_length);
    for j = 1:length(all_phase)
        trial_phase{i}(:,j,:) = all_phase{j}(:,seeds(i)+1-phase_length:seeds(i));
        trial_filt{i}(:,j,:) = lfp_filtered{j}.data(:,seeds(i)+1-phase_length:seeds(i));
    end
end

%% Specify various models

% Specify different models
all_freqs = [2, 7, 19.5, 56, 80];
all_bands = fstop_list;
all_sigmaFreqs = [1, 0.1, 0.01, 0.001, 0.0001];

% model specifcations
model_all.spec = true(1,5);
model_lo_gamma.spec = [false(1,3), true, false];
model_hi_gamma.spec = [false(1,4), true];
model_gamma.spec = [false(1,3), true(1,2)];
model_no_delta.spec = [false, true(1,4)];
model_no_theta.spec = [true, false, true(1,3)];
model_no_beta.spec = [true(1,2), false, true(1,2)];
model_no_gamma.spec = ~model_gamma.spec;
model_beta.spec = ~model_no_beta.spec;
model_theta.spec = ~model_no_theta.spec;
model_delta.spec = ~model_no_delta.spec;
model_delta_theta.spec = [true(1,2), false(1,3)];
model_gamma_delta.spec = model_gamma.spec | model_delta.spec;
model_gamma_theta.spec = model_gamma.spec | model_theta.spec;
model_gamma_beta.spec = model_gamma.spec | model_beta.spec;
model_beta_theta.spec = model_beta.spec | model_theta.spec;
model_hgamma_theta.spec = model_hi_gamma.spec | model_theta.spec;
model_lgamma_theta.spec = model_lo_gamma.spec | model_theta.spec;


% add more models if necessary


%% Use SSPE

% models for "good theta"
models = [model_all, model_no_gamma, model_gamma_theta, model_delta_theta, ...
    model_beta_theta, model_hgamma_theta, model_lgamma_theta, model_theta];

% models for "good-gamma"
% models = [model_all, model_gamma, model_hi_gamma, model_lo_gamma, model_no_delta, ...
%     model_no_theta, model_no_beta, model_gamma_delta, model_gamma_theta, model_gamma_beta];

sspe_pack = cell(1,25);
for iS = 1:25
    this_data = trial_data{iS}.data(1,:);
    this_res = cell(1, length(models));
    for iM = 1:length(models)
        this_model = models(iM);
        this_cfg.freqs = all_freqs(this_model.spec);
        this_cfg.freqBands = fstop_list(this_model.spec);
        this_cfg.Fs = Fs;
        this_cfg.ampVec = repmat(0.99, 1, length(this_cfg.freqs));
        this_cfg.sigmaObs = 1;
        this_cfg.window = round(length(this_data)*0.98);    
        this_cfg.sigmaFreqs = all_sigmaFreqs(this_model.spec);
        this_res{iM}.cfg_sspe = this_cfg;
        % The magic happens here
        [this_res{iM}.omega, this_res{iM}.phase, this_res{iM}.phase_bounds, ...
         this_res{iM}.fullX] =  causalPhaseEM_MKmdl_temp(this_data, ...
                                                    this_res{iM}.cfg_sspe);
    end
    % Calculate PSD for each segment
    wsize = 512;
    [Pxx, F] = pwelch(this_data, rectwin(wsize), wsize/2, [], Fs);
    % Save results to the struct
    sspe_pack{iS}.data = this_data;
    sspe_pack{iS}.res = this_res;
    sspe_pack{iS}.filt_data = trial_filt{iS}(1,:,:);
    sspe_pack{iS}.filt_phase = trial_phase{iS}(1,:,:);
    sspe_pack{iS}.psd.Pxx = Pxx;
    sspe_pack{iS}.psd.F = F;
    
end
%% Plot results




