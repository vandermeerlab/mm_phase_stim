%% Script to create PETHS for all isolated units
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\'
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', ...
    'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', ...
    'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff;
    end
end

%%
function doStuff
    out = [];
    [out.labels, out.stim_peth, out.stim_zpeth, ...
    out.shuf_peth, out.shuf_zpeth, out.tvec, out.is_opto] = deal([]);
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    evs = LoadEvents([]);
    cfg = [];
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.min_cluster_quality = 3;
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
%     cfg.fc = ExpKeys.goodCell;
  
    S = LoadSpikes(cfg);
    % Reject goodCell and continue if other cells in the session
    S_opto = SelectTS([],S,ismember(S.label, ExpKeys.goodCell));
    S_others = SelectTS([],S,~ismember(S.label, ExpKeys.goodCell));
    % Do this for all  
    if ~isempty(S_opto.label) | ~isempty(S_others.label)
        % Set Variables
        if contains(ExpKeys.light_source, 'LASER')
            start_delay = 0.0011;
            stop_delay = 0.0012; %0.0022 is the spike width in the case of 1 msec laser pulse
        else
            start_delay = 0;
            stop_delay = 0;
        end
        stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
        if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
            stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                              stim_on <= ExpKeys.stim_times(2));
        else
            stim_on = [];
        end
    
        % Load CSC
        cfg = []; cfg.fc = ExpKeys.goodLFP;
        if contains(cfg.fc, '-')
            temp = split(cfg.fc,'-');
            cfg.fc = {cat(2,temp{1},'.ncs')};
            csc = LoadCSC(cfg);
            cfg_temp.fc = {cat(2,temp{2},'.ncs')};
            ref = LoadCSC(cfg_temp);
            csc.data = csc.data - ref.data;    
            clear temp ref;
        else
            csc = LoadCSC(cfg);
        end
        
        % Sometimes csc.tvec can be have weird elements because of gaps
        % in recording
        to_remove = find(diff(csc.tvec)<=0);
        while (~isempty(to_remove))
            csc.tvec(to_remove+1) = [];
            csc.data(to_remove+1) = [];
            to_remove = find(diff(csc.tvec)<=0);
        end
        
        % restrict csc
        csc = restrict(csc, iv(ExpKeys.stim_times));
        % Select random stim times for sham
        num_sham = max(max(ExpKeys.goodTrials)) - min(min(ExpKeys.goodTrials)) +1;
        num_shuf = 1000;
        shuf_idx = zeros(num_shuf, num_sham);
        for iShuf = 1:num_shuf
            shuf_idx(iShuf,:) = sort(randsample(csc.tvec, num_sham));
        end

        % Parameters for MUA
        cfg_MUA = []; 
        cfg_MUA.tvec = csc.tvec; % timebase to compute MUA on
        cfg_MUA.sigma = 0.001;

        % Parameters for PETH
        cfg_peth = []; % parameters for PETH
        cfg_peth.window = [-0.2 0.2];
        cfg_peth.dt = 0.001;
        cfg_peth.mode = 'interp';

        for iC = 1:length(S_opto.t)
            fn_prefix = extractBefore(S_opto.label{iC}, '.t');
            this_cell = SelectTS([], S_opto, iC); 
            this_MUA = getMUA(cfg_MUA, this_cell); % "MUA" for one cell is just that cell's firing rate
            this_MUAz = zscore_tsd(this_MUA);
            this_stim_peth = TSDpeth_fast(cfg_peth, this_MUA, ...
                stim_on(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2)));
            this_stim_zpeth = TSDpeth_fast(cfg_peth, this_MUAz, ...
                stim_on(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2)));
            % Now we need to calculate PETHs for shuffles
            [big_peth, big_zpeth] = deal(zeros(num_shuf, length(this_stim_peth.data)));
            for iShuf = 1:num_shuf
                this_shuf = TSDpeth_fast(cfg_peth, this_MUA, shuf_idx(iShuf,:));
                big_peth(iShuf,:) = this_shuf.data;
                this_zshuf = TSDpeth_fast(cfg_peth, this_MUAz, shuf_idx(iShuf,:));
                big_zpeth(iShuf,:) = this_zshuf.data;
            end
            out.labels = [out.labels; fn_prefix];
            out.stim_peth = [out.stim_peth; this_stim_peth.data];
            out.stim_zpeth = [out.stim_zpeth; this_stim_zpeth.data];
            % Instead of saving mean save all of it
            out.shuf_peth = cat(3, out.shuf_peth, big_peth);
            out.shuf_zpeth = cat(3, out.shuf_zpeth, big_zpeth);
            out.tvec = this_stim_peth.tvec;
            out.is_opto = [out.is_opto; 1]; % This is the opto loop
        end
        for iC = 1:length(S_others.t)
            fn_prefix = extractBefore(S_others.label{iC}, '.t');
            this_cell = SelectTS([], S_others, iC); 
            this_MUA = getMUA(cfg_MUA, this_cell); % "MUA" for one cell is just that cell's firing rate
            this_MUAz = zscore_tsd(this_MUA);
            this_stim_peth = TSDpeth_fast(cfg_peth, this_MUA, ...
                stim_on(min(min(ExpKeys.goodTrials)):max(max(ExpKeys.goodTrials))));
            this_stim_zpeth = TSDpeth_fast(cfg_peth, this_MUAz,  ...
                stim_on(min(min(ExpKeys.goodTrials)):max(max(ExpKeys.goodTrials))));
            % Now we need to calculate PETHs for shuffles
            [big_peth, big_zpeth] = deal(zeros(num_shuf, length(this_stim_peth.data)));
            for iShuf = 1:num_shuf
                this_shuf = TSDpeth_fast(cfg_peth, this_MUA, shuf_idx(iShuf,:));
                big_peth(iShuf,:) = this_shuf.data;
                this_zshuf = TSDpeth_fast(cfg_peth, this_MUAz, shuf_idx(iShuf,:));
                big_zpeth(iShuf,:) = this_zshuf.data;
            end
            out.labels = [out.labels; fn_prefix];
            out.stim_peth = [out.stim_peth; this_stim_peth.data];
            out.stim_zpeth = [out.stim_zpeth; this_stim_zpeth.data];
            out.shuf_peth = cat(3, out.shuf_peth, big_peth);
            out.shuf_zpeth = cat(3, out.shuf_zpeth, big_zpeth);
            out.tvec = this_stim_peth.tvec;
            out.is_opto = [out.is_opto; 0]; % This is the non-opto loop
        end

        % Debug here
        dummy = 1;

        % Save the output variable in each session folder
        save('all_pethsv2', 'out');
        close;
    end
end


