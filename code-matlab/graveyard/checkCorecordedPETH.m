%% Script to generate various scatter summary plots
% Assumes that *phase_response.mat already exist in each folder
rng(2023); % Setting the seed for reproducibility
top_dir = 'data\'
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M295'}; %, ...
    %'M074', 'M075', 'M077', 'M078', 'M235', 'M265', ...
    %'M295', 'M320', 'M319', 'M321', 'M325'};
summary = [];
[summary.labels, summary.stim_mode, summary.short_sw, ...
    summary.long_sw, summary.depth, ...
    summary.pre_stim, summary.trial_stim, summary.post_stim, ...
    summary.long_stim, summary.ntrials] = deal([]);
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        summary = doStuff(summary);
    end
end

%%
function s_out = doStuff(s_in)
    s_out = s_in;
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
    if length(S_others.label) ~= 0
        this_fig = figure('WindowState', 'maximized');
        % Set Variables
        max_delay = 0.01; % sec (window for the first response since stimulus)
        max_long_delay = 0.4; % sec
        short_win = 0.02; % sec
        long_win = 0.2; % sec
        short_dt = 0.001; % sec
        long_dt = 0.01; % sec
        fbands = {[2 5], [6 10], [30 55]};
    
        if contains(ExpKeys.light_source, 'LASER')
            start_delay = 0.0011;
            stop_delay = 0.0012; %0.0022 is the spike width in the case of 1 msec laser pulse
        else
            start_delay = 0;
            stop_delay = 0;
        end
    
        pre_stim_on = evs.t{strcmp(evs.label, ExpKeys.pre_trial_stim_on)} + start_delay;
        if ~isempty(pre_stim_on) && ~isempty(ExpKeys.pre_stim_times)
            pre_stim_on = pre_stim_on(pre_stim_on >= ExpKeys.pre_stim_times(1) & ...
                                      pre_stim_on <= ExpKeys.pre_stim_times(2));
        else
            pre_stim_on = [];
        end
    
        stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
        if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
            stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                              stim_on <= ExpKeys.stim_times(2));
        else
            stim_on = [];
        end
    
        post_stim_on = evs.t{strcmp(evs.label, ExpKeys.post_trial_stim_on)} + start_delay;
        if ~isempty(post_stim_on) && ~isempty(ExpKeys.post_stim_times)
            post_stim_on = post_stim_on(post_stim_on >= ExpKeys.post_stim_times(1) & ...
                                        post_stim_on <= ExpKeys.post_stim_times(2));
        else
            post_stim_on = [];
        end
    
        if sum(strcmp(evs.label, ExpKeys.long_stim_on)) ~= 0
            long_stim_on = evs.t{strcmp(evs.label, ExpKeys.long_stim_on)} + start_delay;
        else
            long_stim_on = [];
        end
        if ~isempty(long_stim_on) && ~isempty(ExpKeys.long_stim_times)
            long_stim_on = long_stim_on(long_stim_on >= ExpKeys.long_stim_times(1) & ...
                    long_stim_on <= ExpKeys.long_stim_times(2));
            % Take only the first 25 long stim if special stim present
            if strcmp(ExpKeys.hasSpecialStim, 'Yes')
                long_stim_on = long_stim_on(1:25); 
            end
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

        % Calculate MUA
        cfg_MUA = []; 
        cfg_MUA.tvec = csc.tvec; % timebase to compute MUA on
        cfg_MUA.sigma = 0.001;
        cfg_MUA2 = cfg_MUA;
        cfg_MUA2.sigma = 0.01;

        ax1 = subplot(4,1,1); % Pre-stim
        hold on;
        ax2 = subplot(4,1,2); % Trial-stim
        hold on;
        ax3 = subplot(4,1,3); % Post-stim
        hold on;
        ax4 = subplot(4,1,4); % Long-stim
        hold on;

        % Pre-stim
        
        legend_names = {};
        hold on;
        cfg_peth = []; % parameters for PETH
        cfg_peth.window = [-0.02 0.02];
        cfg_peth.dt = 0.001;
        cfg_peth.mode = 'interp';
        cfg_peth2 = cfg_peth;
        cfg_peth2.window = [-0.2 0.2];
        for iC = 1:length(S_opto.t)
            fn_prefix = extractBefore(S_opto.label{iC}, '.t');
            this_cell = SelectTS([], S_opto, iC); 
            this_MUA = getMUA(cfg_MUA, this_cell); % "MUA" for one cell is just that cell's firing rate
            this_MUAz = zscore_tsd(this_MUA);
            this_out = TSDpeth_fast(cfg_peth, this_MUAz, pre_stim_on);
            plot(ax1, this_out.tvec, this_out.data, 'LineWidth', 1.5);
            this_out = TSDpeth_fast(cfg_peth, this_MUAz, ...
                stim_on(ExpKeys.goodTrials(iC,1):ExpKeys.goodTrials(iC,2)));
            plot(ax2, this_out.tvec, this_out.data, 'LineWidth', 1.5);
            this_out = TSDpeth_fast(cfg_peth, this_MUAz, post_stim_on);
            plot(ax3, this_out.tvec, this_out.data, 'LineWidth', 1.5);
            this_out = TSDpeth_fast(cfg_peth2, this_MUAz, long_stim_on);
            plot(ax4, this_out.tvec, this_out.data, 'LineWidth', 1.5);
            legend_names{length(legend_names)+1} = strcat('Opto cell#', num2str(iC));
        end
        for iC = 1:length(S_others.t)
            fn_prefix = extractBefore(S_others.label{iC}, '.t');
            this_cell = SelectTS([], S_others, iC); 
            this_MUA = getMUA(cfg_MUA, this_cell); % "MUA" for one cell is just that cell's firing rate
            this_MUAz = zscore_tsd(this_MUA);
            this_out = TSDpeth_fast(cfg_peth, this_MUAz, pre_stim_on);
            plot(ax1, this_out.tvec, this_out.data, 'LineWidth', 1.5);
            this_out = TSDpeth_fast(cfg_peth, this_MUAz, ...
                stim_on(min(min(ExpKeys.goodTrials)): max(max(ExpKeys.goodTrials))));
            plot(ax2, this_out.tvec, this_out.data, 'LineWidth', 1.5);
            this_out = TSDpeth_fast(cfg_peth, this_MUAz, post_stim_on);
            plot(ax3, this_out.tvec, this_out.data, 'LineWidth', 1.5);
            this_out = TSDpeth_fast(cfg_peth2, this_MUAz, long_stim_on);
            this_MUA = getMUA(cfg_MUA2, this_cell); % "MUA" for one cell is just that cell's firing rate
            this_MUAz = zscore_tsd(this_MUA);
            plot(ax4, this_out.tvec, this_out.data, 'LineWidth', 1.5);
            legend_names{length(legend_names)+1} = strcat('Non-opto cell#', num2str(iC));
        end

        ax1.YLabel.String = 'Z-score';
        ax1.XLabel.String = 'Time(s)';
        ax1.XTick = [-0.02:0.005:0.02];
        ax1.Title.String = 'Pre-stim';

        ax2.YLabel.String = 'Z-score';
        ax2.XLabel.String = 'Time(s)';
        ax2.XTick = [-0.02:0.005:0.02];
        ax2.Title.String = 'Trial-stim';

        ax3.YLabel.String = 'Z-score';
        ax3.XLabel.String = 'Time(s)';
        ax3.XTick = [-0.02:0.005:0.02];
        ax3.Title.String = 'Post-stim';

        ax4.YLabel.String = 'Z-score';
        ax4.XLabel.String = 'Time(s)';
        ax4.XTick = [-0.2:0.05:0.2];
        ax4.Title.String = 'Long-stim';

        xline(ax1, 0, '--black')
        xline(ax2, 0, '--black')
        xline(ax3, 0, '--black')
        xline(ax4, 0, '--black')

        legend(ax1, legend_names, 'Interpreter', 'none', ...
            'FontSize', 10, 'Location', 'westoutside');

        sgtitle(strcat(ExpKeys.subject_id,'-',ExpKeys.date), 'Interpreter', 'none');
        ax1.Title.FontSize = 13;
        ax2.Title.FontSize = 13;
        ax3.Title.FontSize = 13;
        ax4.Title.FontSize = 13;
        ax1.XAxis.FontSize = 12;
        ax2.XAxis.FontSize = 12;
        ax3.XAxis.FontSize = 12;
        ax4.XAxis.FontSize = 12;


        ax1.Legend.Box = 'off';
        pause(5);
        ax1.Legend.Position(2) = ax1.Position(2);
        ax1.Legend.Position(1) = ax2.Position(1) - ax1.Legend.Position(3) - 0.05;
        
        % Save png file
        dest_path = 'output\';
        print(this_fig, '-dpng',  ...
            strcat(dest_path,ExpKeys.subject_id,'-',ExpKeys.date, '-inhibPETH'));
        close all;
    end
end


