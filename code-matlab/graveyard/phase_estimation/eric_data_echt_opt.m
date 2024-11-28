%Testbed for comparing various phase estimation methods

%TODO make a comparison between control offsets to see where the problem
%comes from

%% Change folder and restrict data to epoch with no opto-stim
cd('data\M20\M20-2019-06-07_dStr_4p6_light_cells_TT6_TT8_min')
LoadExpKeys;
evs = LoadEvents([]);
cfg.fc = {ExpKeys.goodCSC};
csc = LoadCSC(cfg);

eval_csc = restrict(csc, iv(ExpKeys.PreRecord)); % Alternatively use ExpKeys.PostRecord
% try after downsampling the signal
cfg = []; cfg.decimateFactor = 16;
eval_csc = decimate_tsd(cfg, eval_csc); 
Fs = 1./median(diff(eval_csc.tvec));
seed = 4994;
rng(seed);


%% Plot PSD and extract acausal phase
fig = figure('WindowState', 'maximized') ;
subplot(5,10,[1 2]);
wsize = floor(Fs);
[P, F] = pwelch(eval_csc.data, hanning(wsize), wsize/2, [], Fs);
plot(F, 10*log10(P));
xlim([0 120]); xlabel('Frequency (Hz)'); ylabel('Power'); title('PSD');
fbands = {[2 5], [6 10], [20 55], [55 95]};

% The filtering is a time taking step and is best done outside the
% optimization step
filt_phase = cell(length(fbands),1);
for iB = 1:length(fbands)
    cfg_filt.type = 'fdesign'; 
    cfg_filt.f  = fbands{iB};
    filt_lfp = FilterLFP(cfg_filt, eval_csc);
    filt_phase{iB} = angle(hilbert(filt_lfp.data));
end

%% Add Path for the current method
addpath('mm_phase_stim\code-matlab\phase_estimation\ECHT');

%% Try optimization method
delete(gcp('nocreate'));
parpool('local');
popsz = 20;
bounds_window = [10,20];
nSamples = 1000;
tic;
[opt_params, opt_out] = optimize_echt(eval_csc, filt_phase, fbands, Fs, nSamples, popsz, bounds_window);
toc;
delete(gcp('nocreate'));
%% Choose winning parameters and plot scatter and hist for data with no stim 

%Choose optimal parameters from opimization above
% win_length = opt_params.window_length/10; 
win_length = 1.5;
% Choose nEnds in a way such that the smallest sample is win_length long
min_start = ceil(win_length*Fs);
nEnds = randi(length(eval_csc.data) - min_start, nSamples, 1) + min_start;
nStarts = nearest_idx3(eval_csc.tvec(nEnds) - win_length, eval_csc.tvec);
for iB = 1:length(fbands)
    estimated_phase = zeros(1,nSamples);
    true_phase = filt_phase{iB}(nEnds);
    for iS = 1:nSamples
       this_echt = echt(eval_csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
       this_phase = angle(this_echt);
       estimated_phase(iS) = this_phase(end); % The last sample's phase
    end
    subplot(5,10, [(2*iB)+1, (2*iB)+2])
    scatter(estimated_phase, true_phase, 4);
    hold on;
    plot([-pi pi], [-pi pi], 'k'); 
    l = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
    axis tight; grid on; set(gca, 'XTick', -pi:pi/2:pi, 'YTick', -pi:pi/2:pi, 'XTickLabel', l, 'YTickLabel', l);

    title(sprintf("%d Hz - %d Hz", fbands{iB}(1), fbands{iB}(2))); xlabel('estimated phase'); ylabel('true phase');
    
    subplot(5,10,(2*iB)+ 11)
    histogram(true_phase, -pi:2*pi/5:pi, 'FaceColor', 'Cyan');
    title('True phases')
    subplot(5,10,(2*iB)+ 12)
    histogram(estimated_phase, -pi:2*pi/5:pi, 'FaceColor', 'Magenta');
    title('Causal Phases')
end

%% Sanity check: Run this method on white noise and sinusoid and look at phase distributions
nSamples = 1000;
sim_Fs = Fs;
sim_time = 0:1/sim_Fs:10000; %simulate data upto 1000 s
min_start = ceil(win_length*sim_Fs);
pink_data = pinknoise(1,length(sim_time));
nEnds = randi(length(pink_data) - min_start, nSamples, 1) + min_start;
nStarts = nearest_idx3(sim_time(nEnds) - win_length, sim_time);
pink_phase = zeros(length(fbands), nSamples);
sine_phase = zeros(length(fbands), nSamples);
for iB = 1:length(fbands)
    sine_data =  sin(2*pi*mean(fbands{iB})*sim_time);
    for iS = 1:nSamples
        pink_echt = echt(pink_data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), sim_Fs);
        this_phase = angle(pink_echt);
        pink_phase(iB,iS) = this_phase(end); % The last sample's phase
        sine_echt = echt(sine_data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), sim_Fs);
        this_phase = angle(sine_echt);
        sine_phase(iB,iS) = this_phase(end); % The last sample's phase
    end
    
end
clear sine_data
clear pink_data
clear sim_time
% Plot the distributions
for iB = 1:length(fbands)
   subplot(5,10,(2*iB)+ 21)
   histogram(sine_phase(iB,:), -pi:2*pi/5:pi, 'FaceColor', 'Cyan');
   title('Sine phases')
   subplot(5,10,(2*iB)+ 22)
   histogram(pink_phase(iB,:), -pi:2*pi/5:pi, 'FaceColor', 'Magenta');
   title('Pink Phases')
end


%% Obtain Hilbert-transfrom phases
test_csc = restrict(csc, iv([ExpKeys.timeOnWheel, ExpKeys.timeOffWheel]));
% try after downsampling the signal
cfg = []; cfg.decimateFactor = 16;
test_csc = decimate_tsd(cfg, test_csc); 
Fs = 1./median(diff(test_csc.tvec));

stim_times = evs.t{strcmp(evs.label,ExpKeys.laser_on)};
% This sanity check is necessary because of M020
stim_times = stim_times(stim_times > ExpKeys.timeOnWheel);


ISIs = [100 diff(stim_times)']; %100 is used as an arbitrarily large number so that the first stim is always included
keep = ISIs > win_length;

ht_phase = zeros(length(fbands), sum(keep));

for iB = 1:length(fbands)
    cfg_filt = [];
    cfg_filt.type = 'fdesign'; 
    cfg_filt.f  = fbands{iB};
    filt_lfp = FilterLFP(cfg_filt, test_csc); %Filtering the lfp data using bandpass
    filt_phase = angle(hilbert(filt_lfp.data)); %obtaining phases from the filtered data
    ht_phase(iB, :) = filt_phase(nearest_idx3(stim_times(keep), test_csc.tvec));
end

%% Obtain phases thorugh causal method
causal_phase = zeros(length(fbands), sum(keep));
nEnds = nearest_idx3(stim_times(keep), test_csc.tvec);
nStarts = nearest_idx3(stim_times(keep) - win_length, test_csc.tvec);

% To check for the source of discrepancy
c_off = 0.2; % in sec
c_times = stim_times(keep) - c_off;
%Control A (Original snippet, but look at an earlier section)
cEnds1 = nearest_idx3(c_times, test_csc.tvec) - nStarts + 1;
%Control B (Snippet that doesn't include the last bit)
cEnds2 = nearest_idx3(c_times, test_csc.tvec);
cStarts = nearest_idx3(c_times - win_length, test_csc.tvec);
controlA_phase = zeros(length(fbands), sum(keep));
controlB_phase = zeros(length(fbands), sum(keep));
for iB = 1:length(fbands)
    for iS = 1:sum(keep)
       this_echt = echt(test_csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
       this_phase = angle(this_echt);
%        plot(this_phase);
       causal_phase(iB,iS) = this_phase(end); % The last sample's phase
% Control A
       controlA_phase(iB,iS) = this_phase(cEnds1(iS));
% Control B
       %Replace cStarts with nStarts to rule out inclusion of stimuli
       c_echt = echt(test_csc.data(cStarts(iS):cEnds2(iS)), fbands{iB}(1), fbands{iB}(2), Fs); 
       c_phase = angle(c_echt);
       controlB_phase(iB,iS) = c_phase(end);
    end
end

% Plot the distributions
for iB = 1:length(fbands)
   subplot(5,10,(2*iB)+ 31)
   histogram(ht_phase(iB,:), -pi:2*pi/5:pi, 'FaceColor', 'Cyan');
   title('HT phases')
   subplot(5,10,(2*iB)+ 32)
   histogram(causal_phase(iB,:), -pi:2*pi/5:pi, 'FaceColor', 'Magenta');
   title('Causal Phases')
   subplot(5,10,(2*iB)+ 41)
   histogram(controlA_phase(iB,:),-pi:2*pi/5:pi, 'FaceColor', 'Green');
   title('Same Samples, fixed offset')
   subplot(5,10,(2*iB)+ 42)
   histogram(controlB_phase(iB,:), -pi:2*pi/5:pi, 'FaceColor', 'Blue');
   title('Shorter Offset')
end

%% Put some text
subplot(5,10, [11 12])
text(0.1, 0.6, strcat(ExpKeys.subject, '_', ExpKeys.date), 'Interpreter', 'none', 'FontSize', 16)
text(0.1, 0.5 , strcat('Window length Used: ', num2str(win_length), ' sec'), 'Interpreter', 'none', 'FontSize', 16)
text(0.1, 0.4 , strcat('Trials left:  ', num2str(sum(keep))), 'Interpreter', 'none', 'FontSize', 16)
text(0.1, 0.3, 'Method Used: ECHT','Interpreter', 'none', 'FontSize', 16)
box off
grid off
axis off
WriteFig(fig, 'echt', 1)

%% Helper functions
function y = pinknoise(m, n)
    % function: y = pinknoise(m, n)
    % m - number of matrix rows
    % n - number of matrix columns
    % y - matrix with pink (flicker) noise samples 
    %     with mu = 0 and sigma = 1 (columnwise)
    % The function generates a matrix of pink (flicker) noise samples
    % (columnwise). In terms of power at a constant bandwidth, pink  
    % noise falls off at 3 dB/oct, i.e. 10 dB/dec.
    % difine the length of the noise vector and ensure  
    % that M is even, this will simplify the processing
    m = round(m); n = round(n); N = m*n;
    if rem(N, 2)
        M = N+1;
    else
        M = N;
    end
    % generate white noise sequence
    x = randn(1, M);
    % FFT
    X = fft(x);
    % prepare a vector with frequency indexes 
    NumUniquePts = M/2 + 1;     % number of the unique fft points
    k = 1:NumUniquePts;         % vector with frequency indexes 
    % manipulate the left half of the spectrum so the PSD
    % is proportional to the frequency by a factor of 1/f, 
    % i.e. the amplitudes are proportional to 1/sqrt(f)
    X = X(1:NumUniquePts);      
    X = X./sqrt(k);
    % prepare the right half of the spectrum - a conjugate copy of the left
    % one except the DC component and the Nyquist component - they are unique,
    % and reconstruct the whole spectrum
    X = [X conj(X(end-1:-1:2))];
    % IFFT
    y = real(ifft(X));
    % ensure that the length of y is N
    y = y(1, 1:N);
    % form the noise matrix and ensure unity standard 
    % deviation and zero mean value (columnwise)
    y = reshape(y, [m, n]);
    y = bsxfun(@minus, y, mean(y));
    y = bsxfun(@rdivide, y, std(y));
end

