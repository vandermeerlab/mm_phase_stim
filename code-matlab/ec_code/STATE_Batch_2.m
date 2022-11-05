%% initialize

STATE_init
global PARAMS

cd(PARAMS.raw_data_dir)

%% loop through sessions and preprocess the data. based on STATE_batch.m but with more utilitatian sturcture
sess_list = {};
d = dir;
d=d(~ismember({d.name},{'.','..', '._*'}));
for iSess = 1:length(d)
    if ~strcmp(d(iSess).name(1:2), '._')
        sess_list{end+1} = d(iSess).name;

    end
end


to_process = {}; % custom sessions to process


for iSess = 22:length(sess_list) % skip early sessions
    this_date = sess_list{iSess}(5:14);
    this_sub = sess_list{iSess}(1:3);
    this_str = sess_list{iSess}(16:19);
    this_depth = sess_list{iSess}(21:23);
    cd(sess_list{iSess})
    LoadExpKeys()
    for iGoodcell = 1:length(ExpKeys.goodCell)
        if sum(ismember(PARAMS.Good_cells, strrep([this_sub '-' this_date '_' this_depth '_' ExpKeys.goodCell{iGoodcell}(1:end-2)], '-', '_'))) ==1
            
            fprintf('STATE_Batch_2: Processing session: %s...\n', sess_list{iSess})
            
            % RUN DATA preprocessing and spike-phase relationships.
            statedep_latency_better; close all
            
        end
    end
    clear ExpKeys
    cd(PARAMS.raw_data_dir)
end

