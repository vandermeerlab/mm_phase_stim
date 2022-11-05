function STATE_PETH_Check

%remove temp ExpKeys saves
these_files = dir;
for iF = 1:length(these_files)
    if length(these_files(iF).name) <=2
        continue
    else
        if strcmp(these_files(iF).name(1:2), '._')
            delete(these_files(iF).name)
        end
    end
end

font_size = 18;

Matlab_v = version('-release');

if str2double(Matlab_v(1:4)) < 2015
    error('Versions of Matlab before 2015 give add rasters and response rates. Should not be trusted ATM')
end

LoadExpKeys

%% load events
cfg = [];
cfg.eventList = {ExpKeys.laser_on};
cfg.eventLabel = {'laser on'};
laser_on = LoadEvents(cfg);

cfg = [];
cfg.eventList = {'Starting Recording', 'Stopping Recording'};
cfg.eventLabel = {'start', 'stop'};
start_stop = LoadEvents(cfg);

% find the longest recording
for ii = 1:length(start_stop.t{1})
    rec_times(ii) = start_stop.t{2}(ii)-start_stop.t{1}(ii);
end
[duration, main_rec_idx] = max(rec_times);
disp(['Longest Recording interval is ' num2str(duration/60) ' minutes in recording slot number ' num2str(main_rec_idx)])


laser_on = restrict(laser_on, start_stop.t{1}(main_rec_idx), start_stop.t{2}(main_rec_idx));

% check number of pulses.
if length(laser_on.t{1}) ~= 1000 && length(laser_on.t{1}) ~= 600 && length(laser_on.t{1}) ~= 1500
    warning('Wrong number of laser pulses. %0.2f when it should be 1000/1500 or in early sessions 600',length(laser_on.t{1}))
    
end

%% load spikes
cfg = []; cfg.getTTnumbers = 0;
if isfield(ExpKeys, 'goodCell')
    cfg.fc = ExpKeys.goodCell;
end
S = LoadSpikes(cfg);
figure
text(0.35, 0.98,[ExpKeys.subject '_' ExpKeys.date], 'fontsize', font_size)


for iC = 1:length(S.label)
    
    
    %pick this cell
    this_S = SelectTS([], S, iC);
    cell_id = this_S.label{1}(1:end-2);
    cell_id = strrep(cell_id, '_SS', '');
    
    
    %% make a PETH with color coding
    subplot(3,ceil(length(S.label)/3)+1, iC)        % Get the phase labels
    evt_count = 0;
    
    hold on
    % plot(linspace(0,0.001),-1*ones(1,length(linspace(0,0.001))))
    for iEvt = 1:length(laser_on.t{1})
        this_event_spikes = restrict(this_S, laser_on.t{1}(iEvt), laser_on.t{1}(iEvt)+0.025);
        if str2double(Matlab_v(1:4)) > 2015
            if ~isempty(this_event_spikes.t{1}) % versions beyond 2015b change the way rasters look (too small) so this is a hack around. 
                plot((this_event_spikes.t{1}-laser_on.t{1}(iEvt))*1000 ,iEvt, '.k', 'markersize', 5);
                evt_count = evt_count+1;
            end
        else % versions before 2016 which make 'classic' rasters using lines. 
            if ~isempty(this_event_spikes.t{1}) && length(this_event_spikes.t{1}) >2
                plot([(this_event_spikes.t{1}-laser_on.t{1}(iEvt))*1000 (this_event_spikes.t{1}-laser_on.t{1}(iEvt))*1000],[iEvt-0.5 iEvt+0.5], 'k', 'linewidth', 3);
                evt_count = evt_count+1;
            elseif ~isempty(this_event_spikes.t{1}) && length(this_event_spikes.t{1}) ==2 % for odd behaviour that plots two values on a diagonal
                for ii = 1:2
                    plot([(this_event_spikes.t{1}(ii)-laser_on.t{1}(iEvt))*1000 (this_event_spikes.t{1}(ii)-laser_on.t{1}(iEvt))*1000],[iEvt-0.5 iEvt+0.5],'k', 'linewidth',3);
                end
                evt_count = evt_count+1;
            elseif ~isempty(this_event_spikes.t{1}) && length(this_event_spikes.t{1}) ==1
                
                plot([(this_event_spikes.t{1}(1)-laser_on.t{1}(iEvt))*1000 (this_event_spikes.t{1}(1)-laser_on.t{1}(iEvt))*1000],[iEvt-0.5 iEvt+0.5],'k', 'linewidth', 3);
                evt_count = evt_count+1;
            else
                %disp([num2str(iEvt) 'n spikes ' num2str(length(this_event_spikes.t{1}))])
            end
            
        end
    end
    xlim([-1 10])
    ylim([1 length(laser_on.t{1})])
    set(gca, 'xtick',  -1:1:10)
    xlabel('time (ms)')
    ylabel('stim #')
    title(strrep(S.label{iC}, '_', ' '))
    text(7, length(laser_on.t{1})*.9, {[num2str(evt_count) '/' num2str(length(laser_on.t{1})) '  ' num2str(floor((evt_count/length(laser_on.t{1}))*100)) '%']});
    rectangle('position', [0, 0, 1, length(laser_on.t{1})], 'facecolor',[([4,172,218]./255) 0.5], 'edgecolor',[([4,172,218]./255) 0.5] )
    
    %% add the waveform
    if exist([S.label{iC}(1:4) 'AvgWaveforms.csv'], 'file')
        if verLessThan('matlab', '9.0.1')
            wave = csvread([S.label{iC}(1:4) 'AvgWaveforms.csv']);
        else
            wave = readmatrix([S.label{iC}(1:4) 'AvgWaveforms.csv']);
        end
        if size(wave,1) >iC*4
            this_wave = wave((iC*4)-2:(iC*4)+1,:);
            
            plot(9+wave(1,:), (this_wave(1,:)')*2500+length(laser_on.t{1})*.5)
            plot(9+wave(1,:), (this_wave(2,:)')*2500+length(laser_on.t{1})*.3)
            plot(7.5+wave(1,:), (this_wave(3,:)')*2500+length(laser_on.t{1})*.5)
            plot(7.5+wave(1,:), (this_wave(4,:)')*2500+length(laser_on.t{1})*.3)
        end
    end
end
%%
%% load running data
run_file = FindFiles('*run.csv');
if ~isempty(run_file)
    [Sys_time, Real_time, Encoder] = Load_Wheel(run_file{1});
    subplot(3,ceil(length(S.label)/3)+1, length(S.label)+1) ;
    tvec = (Real_time-Real_time(1))/60;
    plot(tvec, [0 ;diff(Encoder)]);
    xlabel('Time (min)');
    ylabel('Encoder position')
    [~,breaks] = findpeaks(diff(tvec), 'MinPeakHeight', 0.05);
    hold on
    plot(tvec([breaks]), zeros(1,length(breaks)), '*r', 'markersize', 4)
    xlim([tvec(1) tvec(end)])
end
%%
%     SetFigure([], gcf)
saveas(gcf, [ExpKeys.subject '_' ExpKeys.date 'check_peth.png']);
saveas_eps([ExpKeys.subject '_' ExpKeys.date 'check_peth'], cd)

