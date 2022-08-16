function STATE_Write_Exp(data_dir)
%% this will create an expkeys file by taking the information from the data files and putting it into an exicutable script
%
% put the path to the data folder as the input.  Or run it without an input
% while in the data folder.


%% initialize

if nargin == 0
    data_dir = cd;
end


%% get the basic information form the name

if isunix
    fname = strsplit(cd, '/');
else
    fname = strsplit(cd, '\');
end
fname = fname{end};
fname = fname(1:strfind(fname,'p')+1);
fname = strrep(fname, '-', '_');

subject_id = fname(1:3);
date_id = fname(5:14);
depth = str2double([fname(strfind(fname, 'p')-1) '.' fname(strfind(fname, 'p')+1)]); % convert the depth from the title p = . in file names used here

if depth <= 2; % if less than 2mm define as cortex
    target = 'crtx';
elseif depth > 2 && depth < 4 % define dorsal striatum
    target = 'dStr';
elseif depth >= 4 % define as vStr
    target = 'vStr';
end


% open the file to write
fid = fopen([fname(1:14) '_ExpKeys.m'], 'w');

%% fill in the consistent information
fprintf(fid, '%% ''This ExpKeys.m was generated using STATE_Write_Exp.m'';\n');

fprintf(fid, ['ExpKeys.version = ' num2str(1) ';\n']);
fprintf(fid, 'ExpKeys.species = ''mouse'';\n');
fprintf(fid, 'ExpKeys.experimenter = ''EC'';\n');
fprintf(fid, 'ExpKeys.behavior = ''wheel'';\n');
fprintf(fid, 'ExpKeys.probe = ''Buz32'';\n');

%% Get the flexible subject information

fprintf(fid, ['ExpKeys.target = ''' target ''';\n']);
fprintf(fid, ['ExpKeys.subject = ''' subject_id ''';\n']);
fprintf(fid, 'ExpKeys.genetics = ''PV-Cre:Ai32'';\n');
fprintf(fid, ['ExpKeys.date = ''' date_id ''';\n']);
fprintf(fid, ['ExpKeys.session_number = ''NaN'';\n']);
fprintf(fid, ['ExpKeys.tetrodeDepths = ' num2str(depth) ';\n']);
if strcmp(subject_id, 'M18') || strcmp(subject_id, 'M19')
    fprintf(fid, 'ExpKeys.goodCSC = ''CSC31.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
    fprintf(fid, 'ExpKeys.goodCSC2 = ''CSC1.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
    
else
    fprintf(fid, 'ExpKeys.goodCSC = ''CSC22.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
end
fprintf(fid, 'ExpKeys.quality = ''NaN''; %%0 is poor, 1 means cell faded, 2 means ok, 3 means great!, NaN means not yet filled in\n');

% to determine hemisphere based on mouse and date of switch.
if strcmp(subject_id, 'M13') && datetime(strrep(date_id, '_', '-'), 'format', 'yyyy-MM-dd') > datetime('2018-12-14', 'format', 'yyyy-MM-dd');
    fprintf(fid, 'ExpKeys.hemisphere = ''L'';\n');
elseif strcmp(subject_id, 'M14') && datetime(strrep(date_id, '_', '-'), 'format', 'yyyy-MM-dd') > datetime('2018-12-03', 'format', 'yyyy-MM-dd');
    fprintf(fid, 'ExpKeys.hemisphere = ''L'';\n');
elseif strcmp(subject_id, 'M15') % was always on R
    fprintf(fid, 'ExpKeys.hemisphere = ''R'';\n');
else
    fprintf(fid, 'ExpKeys.hemisphere = ''R'';\n');
end

fprintf(fid, '\n%%Notes\n');
fprintf(fid, 'ExpKeys.notes = '''';\n');


% laser information
fprintf(fid, '\n%%Laser information\n');
fprintf(fid, 'ExpKeys.laser_duration = %.1f; %% in ms.\n', 1);
fprintf(fid, 'ExpKeys.laser_mW = ''NaN''; %% in mW.\n');
fprintf(fid, ['ExpKeys.laser_wave = ''' num2str(473) ''';\n']);
fprintf(fid, 'ExpKeys.laser_type = ''Shanghai Laser'';\n');
fprintf(fid, 'ExpKeys.laser_shutter = ''yes'';\n');
fprintf(fid, 'ExpKeys.fibre_cable = %.0d ;%% microns\n', 125);
fprintf(fid, 'ExpKeys.fibre_probe = %.0d ;%% microns\n', 50);
fprintf(fid, 'ExpKeys.fibre_NA = %.2f ;%% \n', 0.39);


%% experimental variables: NLX digital I/O
fprintf(fid, '\n%%NLX digital I/O codes\n');

if strcmp(subject_id, 'M16') || strcmp(subject_id, 'M17') || strcmp(subject_id, 'M18') || strcmp(subject_id, 'M19')% different ids for variabel stim ISI version
    fprintf(fid, 'ExpKeys.laser_on = ''TTL Input on AcqSystem1_0 board 0 port 2 value (0x0002).''; %% when the 1ms laser came on \n');
    fprintf(fid, 'ExpKeys.var_trig = ''TTL Input on AcqSystem1_0 board 0 port 2 value (0x0004).''; %% when the 1ms laser came on \n');
    fprintf(fid, 'ExpKeys.wide_laser_on = ''TTL Input on AcqSystem1_0 board 0 port 2 value (0x0010).''; %% when the 50ms laser came on.  Only used in verfication recordings at the end of the session.\n');
    fprintf(fid, 'ExpKeys.short_laser_on = ''TTL Input on AcqSystem1_0 board 0 port 2 value (0x0008).''; %% when the 50ms laser came on.  Only used in verfication recordings at the end of the session.\n');
    fprintf(fid, 'ExpKeys.stim_mode = ''variable ISI''; %% if fixed then just used ISI, if variable then ISI is between 0-1500ms \n');
    fprintf(fid, 'ExpKeys.ISI = %.2f ;%% \n', 2);
    fprintf(fid, 'ExpKeys.stim_ran_range = [%.2f %.2f];%% \n', 0, 1500);
    
    
else
    fprintf(fid, 'ExpKeys.laser_on = ''TTL Input on AcqSystem1_0 board 0 port 2 value (0x000A).''; %% when the 1ms laser came on \n');
    fprintf(fid, 'ExpKeys.wide_laser_on = ''TTL Input on AcqSystem1_0 board 0 port 2 value (0x0012).''; %% when the 50ms laser came on.  Only used in verfication recordings at the end of the session.\n');
    fprintf(fid, 'ExpKeys.stim_mode = ''fixed ISI''; %% if fixed then just used ISI, if variable then ISI is between 0-1500ms \n');
    fprintf(fid, 'ExpKeys.ISI = %.2f ;%% \n', 3);
    fprintf(fid, 'ExpKeys.stim_ran_range = %.2f ;%% \n', 0);
end

fprintf(fid, 'ExpKeys.event_off = ''TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000).''; %% off signal for all events\n');


%% experimental timing variables;  Used to get the time of the main recording only at this time.
fprintf(fid, '\n%%Timing\n');

% find the main wheel section which will be the longest.
cfg_Evt = [];
cfg_Evt.eventList = {'Starting Recording', 'Stopping Recording'};
cfg_Evt.eventLabel = {'start', 'stop'};
Evt = LoadEvents(cfg_Evt);

for ii = 1:length(Evt.t{1})
    rec_times(ii) = Evt.t{2}(ii)-Evt.t{1}(ii);
end
[duration, main_rec_idx] = max(rec_times);
disp(['Longest Recording interval is ' num2str(duration/60) ' minutes in recording slot number ' num2str(main_rec_idx)])

fprintf(fid, 'ExpKeys.timeOnWheel = %.4d;\n',Evt.t{1}(main_rec_idx));
fprintf(fid, 'ExpKeys.timeOffWheel = %.4d;\n',Evt.t{2}(main_rec_idx) );

fprintf(fid, 'ExpKeys.PreRecord = [%.4d %.4d];\n',Evt.t{1}(1),Evt.t{2}(1));
if length(Evt.t{2}) > main_rec_idx
    fprintf(fid, 'ExpKeys.PostRecord = [%.4d %.4d];\n',Evt.t{1}(main_rec_idx+1),Evt.t{2}(main_rec_idx+1));
else
    fprintf(fid, 'ExpKeys.PostRecord = NaN;\n');
end
fprintf(fid, 'ExpKeys.Hundred_Stims = [%.4d %.4d];\n',Evt.t{1}(2),Evt.t{2}(2));


fclose(fid);
disp('ExpKeys written')


