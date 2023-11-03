%% Assumes that good LFPs have been picked out

top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
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

%
function doStuff
    % Declaring variables
    % Setting up parameters
    fbands = {[2 5], [6 10], [30 55]};

    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    evs = LoadEvents([]);

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
    Fs = 1/median(diff(csc.tvec));

    csc = restrict(csc, iv(ExpKeys.stim_times));

    % Plot phase distribution of stim before and after rejecting plots
    if contains(ExpKeys.light_source, 'LASER')
        start_delay = 0.0011;
    else
        start_delay = 0;
    end

    stim_on = evs.t{strcmp(evs.label, ExpKeys.trial_stim_on)} + start_delay;
    if ~isempty(stim_on) && ~isempty(ExpKeys.stim_times)
        stim_on = stim_on(stim_on >= ExpKeys.stim_times(1) & ...
                          stim_on <= ExpKeys.stim_times(2));
    else
        stim_on = [];
    end

    causal_phase = nan(length(fbands), length(stim_on));
    nEnds = nearest_idx3(stim_on, csc.tvec);

    for iB = 2%1:length(fbands)
        win_length  = 0.5; % we are using this window length because we don't see much of a difference in phase estimation % TODO : Find a way to communicate this to the readers
        nStarts = nearest_idx3(stim_on - win_length, csc.tvec);
        for iS = 1:length(stim_on)
            this_echt = echt(csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
            this_phase = angle(this_echt);
            causal_phase(iB,iS) = this_phase(end); % The last sample's phase
            % Diagnostic to check how are angles assined (Uncomment to run)
            diag_fig = figure;
            ax1 = subplot(2,1,1);
            plot(abs(this_echt))
            ax2 = subplot(2,1,2);
            plot(angle(this_echt))
            hold on
            yline(0, 'red')
            yline(pi/2, 'green')
            yline(pi, 'black')
            linkaxes([ax1,ax2],'x')
            close(diag_fig)
        end
    end
    % Assume you are in the correct folder
    save('stim_phases','causal_phase'); % should add option to save in specified output dir
    close;
end