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
trial_phase = cell(1,25);
trial_filt = cell(1,25);
seeds = randi(length(lfp_prefilt{1}.data),1,25);
for i = 1:25
    t2 = lfp_prefilt{1}.tvec(seeds(i));
    t1_idx = nearest_idx3(t2 - 5, lfp_prefilt{1}.tvec);
    % If the trial length is not even, make it so!
    if rem(seeds(i)-t1_idx, 2) == 0
        t1_idx = t1_idx-1;
    end
    t1 = lfp_prefilt{1}.tvec(t1_idx);
    trial_data{i} = restrict(lfp_prefilt{1}, iv(t1,t2));
    trial_filt{i} = zeros(length(trial_lfp.label),length(all_phase),length(trial_data{i}.tvec));
    trial_phase{i} = zeros(length(trial_lfp.label),length(all_phase),length(trial_data{i}.tvec));
    for j = 1:length(all_phase)
        trial_phase{i}(:,j,:) = all_phase{j}(:,t1_idx:seeds(i));
        trial_filt{i}(:,j,:) = lfp_filtered{j}.data(:,t1_idx:seeds(i));
    end
end

%% Use SSPE
for iD = 1:25
    for iF = 2%1:length(f_list)
        this_data = trial_data{1}.data(1,:);
        this_filt_data = trial_filt{iD}(1,iF,:);
        this_filt_phase = trial_phase{iD}(1,iF,:);
        fig = figure;
        cfg_sspe.freqs = mean(f_list{iF}); % initialization using above not great at identifying starting freq
        cfg_sspe.Fs = Fs;
        cfg_sspe.ampVec = 0.99;
        cfg_sspe.sigmaFreqs = 10^(1-iF);;
        cfg_sspe.sigmaObs = 1;
        cfg_sspe.window = length(this_data)/2;
        cfg_sspe.lowFreqBand = fstop_list{iF};
        
%         [phase,phaseBounds, fullX] = causalPhaseEM_MKmdl(this_data, cfg_sspe);
        [~,~,~,~,stateVec,~ ] = fit_MKModel_multSines(this_data, ...
                                    cfg_sspe.freqs, cfg_sspe.Fs, ...
                                    cfg_sspe.ampVec, cfg_sspe.sigmaFreqs, ...
                                    cfg_sspe.sigmaObs);
        phase = angle(stateVec( 1, :) + 1i*stateVec( 2, :));
        wsize = 1024;
        [Pxx, F] = pwelch(this_data, rectwin(wsize), wsize/2, [], Fs);
        tiledlayout(3,1)
        ax1 = nexttile;
        plot(F, 10*log10(Pxx));
        xlim([0,120]);
        title('PSD')
        ax2 = nexttile;
        plot((1:length(phase))/Fs,phase);
        hold on;
        plot((1:length(phase))/Fs,this_filt_phase(:))
        ax3 = nexttile;
        plot((1:length(phase))/Fs,this_data);
        hold on
        plot((1:length(phase))/Fs,this_filt_data(:));
        linkaxes([ax2 ax3], 'x');
        dummy = 0;
    end
end
%%
