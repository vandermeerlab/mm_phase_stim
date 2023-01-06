%Testbed for comparing various phase estimation methods



%% Change folder and restrict data to epoch with no opto-stim
cd('E:\Dropbox (Dartmouth College)\EC_State_inProcess\M20\M20-2019-06-07_dStr_3p8_light_cells_TT6_TT8_min')
LoadExpKeys;
evs = LoadEvents([]);
cfg.fc = {ExpKeys.goodCSC};
csc = LoadCSC(cfg);

eval_csc = restrict(csc, iv(ExpKeys.PreRecord)); % Alternatively use ExpKeys.PostRecord
% try after downsampling the signal
cfg = []; cfg.decimateFactor = 16;
eval_csc = decimate_tsd(cfg, eval_csc); 
Fs = 1./median(diff(eval_csc.tvec));

%% Plot PSD and extract acausal phase
fig = figure('WindowState', 'maximized') ;
subplot(2,10,[1 2]);
wsize = floor(Fs);
[P, F] = pwelch(eval_csc.data, hanning(wsize), wsize/2, [], Fs);
plot(F, 10*log10(P));
xlim([0 120]); xlabel('Frequency (Hz)'); ylabel('Power'); title('PSD');

fbands = {[2 5], [6 10], [20 55], [55 95]};
true_phase = cell(length(fbands),1);
for iB = 1:length(fbands)
    cfg_filt.type = 'fdesign'; 
    cfg_filt.f  = fbands{iB};
    filt_lfp = FilterLFP(cfg_filt, eval_csc);
    true_phase{iB} = angle(hilbert(filt_lfp.data));
end

%% Add Path for the current method
addpath('D:\vstr_phase_stim\mm_phase_stim\code-matlab\phase_estimation\ECHT');

%% Run the current method on eval_data

nSamples = 1000;
win_length = 1.5; %in seconds
nEnds = randi(length(eval_csc.data), nSamples, 1);
nStarts = nearest_idx3(eval_csc.tvec(nEnds) - win_length, eval_csc.tvec);
for iB = 1:length(fbands)
    estimated_phase = zeros(1,nSamples);
    for iS = 1:nSamples
       this_echt = echt(eval_csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
       this_phase = angle(this_echt);
       estimated_phase(iS) = this_phase(end); % The last sample's phase
    end
    subplot(2,10, [(2*iB)+1, (2*iB)+2])
    scatter(estimated_phase, true_phase{iB}(nEnds), 4);
    hold on;
    plot([-pi pi], [-pi pi], 'k'); 
    l = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
    axis tight; grid on; set(gca, 'XTick', -pi:pi/2:pi, 'YTick', -pi:pi/2:pi, 'XTickLabel', l, 'YTickLabel', l);

    title(sprintf("%d Hz - %d Hz", fbands{iB}(1), fbands{iB}(2))); xlabel('estimated phase'); ylabel('true phase');
end

%% TODO: Plot the distribution of phases with just hilbert transform and the causal method

%% Obtain Hilbert-transfrom phases
test_csc = restrict(csc, iv([ExpKeys.timeOnWheel, ExpKeys.timeOffWheel]));
% try after downsampling the signal
cfg = []; cfg.decimateFactor = 16;
Fs = 1./median(diff(test_csc.tvec));

stim_times = evs.t{strcmp(evs.label,ExpKeys.laser_on)};
% This sanity check is necessary because of M020
stim_times = stim_times(stim_times > ExpKeys.timeOnWheel);
ISIs = [100 diff(stim_times)']; %100 is used as an arbitrarily large number so that the first stim is always included
keep = ISIs > win_length;

test_csc = decimate_tsd(cfg, test_csc);
ht_phase = zeros(length(fbands), sum(keep));

for iB = 1:length(fbands)
    cfg_filt = [];
    cfg_filt.type = 'fdesign'; 
    cfg_filt.f  = fbands{iB};
    filt_lfp = FilterLFP(cfg_filt, test_csc);
    filt_phase = angle(hilbert(filt_lfp.data));
    ht_phase(iB, :) = filt_phase(nearest_idx3(stim_times(keep), test_csc.tvec));
end

%% Obtain phases thorugh causal method
causal_phase = zeros(length(fbands), sum(keep));
nEnds = nearest_idx3(stim_times, test_csc.tvec);
nStarts = nearest_idx3(stim_times - win_length, test_csc.tvec);
for iB = 1:length(fbands)
    for iS = 1:sum(keep)
       this_echt = echt(test_csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
       this_phase = angle(this_echt);
%        plot(this_phase);
       causal_phase(iB,iS) = this_phase(end); % The last sample's phase
    end
end

%% Plot the distributions
for iB = 1:length(fbands)
   subplot(2,10,(2*iB)+ 11)
   histogram(ht_phase(iB,:), 5, 'FaceColor', 'Cyan');
   title('HT phases')
   subplot(2,10,(2*iB)+ 12)
   histogram(causal_phase(iB,:), 5, 'FaceColor', 'Magenta');
   title('Causal Phases')
end

%% Put some text
subplot(2,10, [11 12])
text(0.1, 0.6, strcat(ExpKeys.subject, '_', ExpKeys.date), 'Interpreter', 'none', 'FontSize', 16)
text(0.1, 0.5 , strcat('Window length Used: ', num2str(win_length), ' sec'), 'Interpreter', 'none', 'FontSize', 16)
text(0.1, 0.4 , strcat('Trials left:  ', num2str(sum(keep))), 'Interpreter', 'none', 'FontSize', 16)
text(0.1, 0.3, 'Method Used: ECHT','Interpreter', 'none', 'FontSize', 16)
box off
grid off
axis off
WriteFig(fig, 'echt', 1)

