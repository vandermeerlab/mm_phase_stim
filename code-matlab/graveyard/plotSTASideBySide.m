%% Assumes that good LFPs have been picked out

top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};

for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff
    end
end

%%
function doStuff
    % Declaring variables
    % Setting up parameters
    fbands = {[2 5], [6 10],[30 55]};
    c_list = {'red', 'blue','green'};
    phase_bins = -pi:2*pi/5:pi;

    LoadExpKeys;
    evs = LoadEvents([]);
    cfg_spk = [];
    cfg_spk.fc = ExpKeys.goodCell;
    if isempty(cfg_spk.fc)
        return
    end

        
    cell_ID = 'M017-2019-02-16-TT04_2.t'; % Great oscillations
    if isempty(find(strcmp(cell_ID, ExpKeys.goodCell)))
       return
    end

    fn_prefix = split(pwd, '\');
    fn_prefix = fn_prefix{end};

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
    Fs = csc.cfg.hdr{1}.SamplingFrequency;

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
    control_on = stim_on + rand(size(stim_on)) - 0.5;

    w =[-0.25,0.25]; % time window to compute STA over
    this_tvec = w(1):1/Fs:w(2); % time axis for STA
    [this_on_snip, this_control_snip] = deal(zeros(length(stim_on), length(this_tvec)));
    for iEvt = 1:length(stim_on) % for each stim ...
        this_start = nearest_idx3(stim_on(iEvt)+w(1), csc.tvec);
        control_start = nearest_idx3(control_on(iEvt)+w(1), csc.tvec);
        % grab LFP snippet for this window
        this_on_toAdd = csc.data(this_start:this_start+length(this_tvec)-1);
        this_on_snip(iEvt,:) = this_on_toAdd';
        this_on_toAdd = csc.data(control_start:control_start+length(this_tvec)-1);
        this_control_snip(iEvt,:) = this_on_toAdd';   
    end

    fig = figure('WindowState', 'Maximized');
    plot(this_tvec, mean(this_on_snip), 'Color', 'magenta', 'LineWidth', 2);
    hold on
    plot(this_tvec, mean(this_control_snip), 'Color', [0.6, 0.6, 0.6], 'LineWidth', 2);
    xline(0, '--black')
    legend({'STA centered at stim-onset','Control'})
    sgtitle('DSP filter (0.1 - 400 Hz) applied during Acquisition', 'Interpreter', 'none');
%     sgtitle('DSP filter OFF during Acquistion')
    xlim([-0.1 0.1])
    xticks([- 0.1 0 0.1]);
    xlabel('Time (sec)')
    yticks([])
    box off
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength = [0.03 0.02];
    fontname(fig, 'Helvetica');
    fontsize(fig, 45, 'points');
    fig.Renderer = 'painters';
    dummy = 1;
end
