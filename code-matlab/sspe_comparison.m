cd('/Users/manishm/Work/vanDerMeerLab/Tutorials/data/R016-2012-10-08');
% goodGamma: {'R016-2012-10-03-CSC04d.ncs'
% goodSWR: {'R016-2012-10-03-CSC02b.ncs'}
% goodTheta: {'R016-2012-10-03-CSC02b.ncs'}
% cfg_in.fc = {'LFP4.ncs', 'LFP6.ncs', 'LFP28.ncs', 'LFP30.ncs'};
cfg_in.fc = {'R016-2012-10-08-CSC02d.ncs'};
all_lfp  = LoadCSC(cfg_in);
evs = LoadEvents([]);
%%
% 
% t_start = evs.t{5}; % pre trial baseline recording start
% t_end = evs.t{4}; % pre trial baseline recording ended
% t_iv = iv(t_start, t_end);
% trial_lfp  = restrict(all_lfp, t_iv);
trial_lfp = all_lfp;
%% Check PSD
Fs = trial_lfp.cfg.hdr{1}.SamplingFrequency;
wsize = Fs;
figure;
for i = 1:length(trial_lfp.label)
    [Pxx, F] = pwelch(trial_lfp.data(i,:), rectwin(wsize), wsize/2, [], Fs);
    plot(F, 10*log10(Pxx));
    hold on
end
legend(trial_lfp.label)
xlim([0,120]);

%% Estimate ground truth phases using FilterLFP + Hilbert transform

% make a cell array of the various processed LFPS
lfp_prefilt{1} = trial_lfp;

f_list = {[2 5], [6 10],[15 25], [30 40],[40 60], [60 80]};
f_list_label = {'2 - 5', '6 - 10', '15 - 25', '30 - 40', '40 - 60', '60 - 80'};
fstop_list = {[1.5 5.5], [5.5 10.5],[14 26], [28 42],[38 62], [55 85]};

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
seeds = randi(length(lfp_prefilt{1}.data),1,25);
for i = 1:25
    t2 = lfp_prefilt{1}.tvec(seeds(i));
    t1 = nearest_idx3(t2 - 2.5, lfp_prefilt{1}.tvec);
    % If the trial length is not even, make it so!
    if rem(seeds(i)-t1, 2) == 0
        t1 = t1-1;
    end
    t1 = lfp_prefilt{1}.tvec(t1);
    trial_data{i} = restrict(lfp_prefilt{1}, iv(t1,t2));     
end

%% Use SSPE
for iD = 1:25
    close all;
    fig = figure;
    this_data = trial_data{iD}.data;
    this_tvec = trial_data{iD}.tvec;
    subplot(8,1,1);
    plot(this_tvec, this_data);
    subplot(8,1,2);
    [Pxx, F] = pwelch(this_data, rectwin(wsize), wsize/2, [], Fs);
    plot(F, 10*log10(Pxx));
    xlim([0,120]);
    
    
    
    for iF = 1%:length(f_list)
        close all;
        fig = figure;
        cfg_sspe.freqs = 60%mean(f_list{iF}); % initialization using above not great at identifying starting freq
        cfg_sspe.Fs = Fs;
        cfg_sspe.ampVec = 0.99;
        cfg_sspe.sigmaFreqs = 1;
        cfg_sspe.sigmaObs = 1;
        cfg_sspe.window = length(trial_data{1}.data)/2;
        cfg_sspe.lowFreqBand = [50,70]%f_list{iF};
        [phase,phaseBounds, fullX] = causalPhaseEM_MKmdl(trial_data{1}.data, cfg_sspe);
        phase = reshape(phase', size(phase,1) * size(phase,2),1);
        phaseBounds = reshape(permute(phaseBounds,[2,1,3]), size(phaseBounds,1) * size(phaseBounds,2),size(phaseBounds,3));
        fullX = reshape(permute(fullX,[2,1,3]), size(fullX,1) * size(fullX,2),size(fullX,3));
        wsize = 1024;
        [Pxx, F] = pwelch(trial_data{1}.data, rectwin(wsize), wsize/2, [], Fs);
        subplot(2,1,1)
        plot(F, 10*log10(Pxx));
        xlim([0,120]);
        subplot(2,1,2)
        plot(phase);
    end
end





