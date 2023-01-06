%Testbed for comparing various phase estimation methods



%% Change folder and restrict data to epoch with no opto-stim
cd('E:\temp_phase_stim\ED\M16-2019-02-15_vStr_4p2_light_cells_TT6_min')
LoadExpKeys;
cfg.fc = {ExpKeys.goodCSC};
csc = LoadCSC(cfg);

eval_csc = restrict(csc, iv(ExpKeys.PreRecord)); % Alternatively use ExpKeys.PostRecord
% try after downsampling the signal
cfg = []; cfg.decimateFactor = 16;
eval_csc = decimate_tsd(cfg, eval_csc); 
Fs = 1./median(diff(eval_csc.tvec));

%% Plot PSD and extract acausal phase
fig = figure;
subplot(2,5,1);
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
addpath('D:\vstr_phase_stim\mm_phase_stim\code-matlab\phase_estimation\Czrenner');

%% Run the current method on eval_data

%%Train

forders = [2 3 4 5]; % low gamma
D = {};
for iOrd = 1:length(forders)
    D{iOrd} = designfilt('bandpassfir', 'FilterOrder', round(forders(iOrd) * (Fs/33)), 'CutoffFrequency1', 19 , 'CutoffFrequency2', 56, 'SampleRate', Fs, 'DesignMethod', 'window'); % gamma
end
%%
csc_filtered = eval_csc;
csc_filtered.data = filtfilt(D{end}, csc_filtered.data);
%csc_filtered.data = sin(2*pi*csc_filtered.tvec*4);

data_len = length(csc_filtered.data); 
data_phase = angle(hilbert(csc_filtered.data));

nIter = 100; clear this_data; clear ALL;
for iI = nIter:-1:1
   
    this_sample = 10000 + ceil(rand(1).*(data_len-10000));
    
    ALL.data(:, iI) = csc_filtered.data(this_sample-10000:this_sample)';
    ALL.true_phase(iI) = data_phase(this_sample);
    
end

optimal_parameters = phastimate_optimize(ALL.data, ALL.true_phase, D, [1 length(forders)], [600 2000], [20 80], [20 100], 64)


%% Test  

csc_filtered = eval_csc;
csc_filtered.data = filtfilt(D{1}, csc_filtered.data);

data_len = length(csc_filtered.data); 
data_phase = angle(hilbert(csc_filtered.data));

nIter = 1000;
window_len = floor(Fs*3);
for iI = nIter:-1:1
   
    this_sample = window_len + ceil(rand(1).*(data_len-window_len));
    
    this_data = csc_filtered.data(this_sample-window_len:this_sample)';
    [ALL.phase(iI), ALL.amplitude(iI)] = phastimate(this_data, D{1}, 20, 96, 1988); % low-gamma
    %[ALL.phase(iI), ALL.amplitude(iI)] = phastimate(this_data, D{1}, 23, 20, 64); % delta
    ALL.true_phase(iI) = data_phase(this_sample);
    
end
%%
subplot(1,2,2)
scatter(ALL.phase, ALL.true_phase, 10, ALL.amplitude);
hold on;
plot([-pi pi], [-pi pi], 'k'); 
l = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
axis tight; grid on; set(gca, 'FontSize', 18, 'XTick', -pi:pi/2:pi, 'YTick', -pi:pi/2:pi, 'XTickLabel', l, 'YTickLabel', l);

title('Data - lowgamma'); xlabel('estimated phase'); ylabel('true phase');
