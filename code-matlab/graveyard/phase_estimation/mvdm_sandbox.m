%% set path
addpath('D:\vstr_phase_stim\mm_phase_stim\code-matlab\phase_estimation\Czrenner');
%addpath('C:\Users\mvdm\Documents\GitHub\phastimate');

cd('E:\temp_phase_stim\ED\M16-2019-02-15_vStr_4p2_light_cells_TT6_min');
%cd('C:\Data\R016-2012-10-03');

%% params
PARAMS.twin = [-2 0];

%% example from toolbox code
Fs = 200;
data = [sin([1:Fs]/Fs*5*2*pi)' sin([1:Fs]/Fs*6*2*pi)'] + randn(Fs,2);
D = designfilt('bandpassfir', 'FilterOrder', round(Fs/5), 'CutoffFrequency1', 4, 'CutoffFrequency2', 7, 'SampleRate', Fs);
phase = phastimate(data, D, 25, 20, 64)

%% extended to check output against acausal hilbert phase on sine waves
data_len = Fs * 10; 
data = sin([1:data_len]/Fs*5*2*pi)';
data_phase = angle(hilbert(data));
clear ALL

nIter = 100;
for iI = 1:nIter
   
    this_sample = 120 + ceil(rand(1).*(data_len-120));
    
    this_data = data(1:this_sample);
    ALL.phase(iI) = phastimate(this_data, D, 25, 20, 64); % MvdM: only gives one phase value per channel -- final phase?
    ALL.true_phase(iI) = data_phase(this_sample);
    
end

subplot(221); cla;
plot(ALL.phase, ALL.true_phase, '.');
hold on;
plot([-pi pi], [-pi pi], 'k'); 
l = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
axis tight; grid on; set(gca, 'FontSize', 18, 'XTick', -pi:pi/2:pi, 'YTick', -pi:pi/2:pi, 'XTickLabel', l, 'YTickLabel', l);

title('Sinusoid'); xlabel('estimated phase'); ylabel('true phase');

%%
LoadExpKeys;
please = [];
please.fc = {ExpKeys.goodCSC};
csc = LoadCSC(please);

Fs = 1./median(diff(csc.tvec));

%% restrict data to pre-record
clean_csc = restrict(csc, iv(ExpKeys.PreRecord));

%% train model

%% Create Filter
clear D ALL
peak_frequency = 3.5;
ord = [3 4 5]; % low gamma
%filter_orders = 10; % delta
filter_count = 1;
for iOrd = 1:length(ord)
D{filter_count} = designfilt('bandpassfir', 'FilterOrder', round(ord(iOrd) * (Fs/peak_frequency)), 'StopbandFrequency1', 1, 'PassbandFrequency1', 2, 'PassbandFrequency2', 5, 'StopbandFrequency2', 6, 'SampleRate', Fs, 'DesignMethod', 'ls');
filter_count = filter_count + 1;
end

% D{1} = designfilt('bandpassiir', ...       % Response type
%        'StopbandFrequency1',2, ...    % Frequency constraints
%        'PassbandFrequency1',3, ...
%        'PassbandFrequency2',5, ...
%        'StopbandFrequency2',6, ...
%        'StopbandAttenuation1',60, ...   % Magnitude constraints
%        'PassbandRipple',2, ...
%        'StopbandAttenuation2',60, ...
%        'DesignMethod','ellip', ...      % Design method
%        'SampleRate',2000);
%%
       
csc_filtered = csc;
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

optimal_parameters = phastimate_optimize(ALL.data, ALL.true_phase, D, [1 length(filter_orders)], [1000 2000], [20 80], [20 100], 128)

%% test
%clear D;
%D{1} = designfilt('bandpassfir', 'FilterOrder', filter_orders(3), 'CutoffFrequency1', 45, 'CutoffFrequency2', 65, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', 2000); % low-gamma
% D{1} = designfilt('bandpassiir', ...       % Response type
%        'StopbandFrequency1',2, ...    % Frequency constraints
%        'PassbandFrequency1',3, ...
%        'PassbandFrequency2',5, ...
%        'StopbandFrequency2',6, ...
%        'StopbandAttenuation1',30, ...   % Magnitude constraints
%        'PassbandRipple',2, ...
%        'StopbandAttenuation2',60, ...
%        'DesignMethod','ellip', ...      % Design method
%        'MatchExactly','passband', ...   % Design method options
%        'SampleRate',2000);
csc_filtered = csc;
csc_filtered.data = filtfilt(D{1}, csc_filtered.data);

data_len = length(csc_filtered.data); 
data_phase = angle(hilbert(csc_filtered.data));

nIter = 1000;
window_len = 1826;
for iI = nIter:-1:1
   
    this_sample = window_len + ceil(rand(1).*(data_len-window_len));
    
    this_data = csc_filtered.data(this_sample-window_len:this_sample)';
    [ALL.phase(iI), ALL.amplitude(iI)] = phastimate(this_data, D{1}, 20, 95, 1988); % low-gamma
    %[ALL.phase(iI), ALL.amplitude(iI)] = phastimate(this_data, D{1}, 23, 20, 64); % delta
    ALL.true_phase(iI) = data_phase(this_sample);
    
end
%%
figure;
scatter(ALL.phase, ALL.true_phase, 10, ALL.amplitude);
hold on;
plot([-pi pi], [-pi pi], 'k'); 
l = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
axis tight; grid on; set(gca, 'FontSize', 18, 'XTick', -pi:pi/2:pi, 'YTick', -pi:pi/2:pi, 'XTickLabel', l, 'YTickLabel', l);

title('Data - lowgamma'); xlabel('estimated phase'); ylabel('true phase');
