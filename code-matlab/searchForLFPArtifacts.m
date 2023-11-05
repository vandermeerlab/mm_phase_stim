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

%     cell_ID = 'M019-2019-04-15-TT06_1.t'; % Great oscillations
%     cell_ID = 'M020-2019-06-07-TT08_1.t'; % Great artifact
%     if isempty(find(strcmp(cell_ID, ExpKeys.goodCell)))
%        return
%     end

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

    % Compute filtered LFP and Hilbert Phases
    [filt_csc, filt_phase] = deal(cell(size(fbands)));

    for iF = 1:length(fbands)
         cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = fbands{iF};
         filt_csc{iF} = FilterLFP(cfg_filt, csc);
         filt_phase{iF} = angle(hilbert(filt_csc{iF}.data));
    end

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
    [filt_snip, phase_snip] = deal(zeros(length(fbands), length(stim_on), length(this_tvec)));
    [stim_phase, control_phase] = deal(zeros(length(fbands), length(stim_on)));
    for iEvt = 1:length(stim_on) % for each stim ...
        this_start = nearest_idx3(stim_on(iEvt)+w(1), csc.tvec);
        control_start = nearest_idx3(control_on(iEvt)+w(1), csc.tvec);
        % grab LFP snippet for this window
        this_on_toAdd = csc.data(this_start:this_start+length(this_tvec)-1);
        this_on_snip(iEvt,:) = this_on_toAdd';
        for iF = 1:3
            filt_snip(iF,iEvt,:) = filt_csc{iF}.data(this_start:this_start+length(this_tvec)-1)';
            phase_snip(iF,iEvt,:) = filt_phase{iF}(this_start:this_start+length(this_tvec)-1)';
            stim_phase(iF,iEvt) = filt_phase{iF}(nearest_idx3(stim_on(iEvt), csc.tvec));
            control_phase(iF,iEvt) = filt_phase{iF}(nearest_idx3(control_on(iEvt), csc.tvec));
        end
        this_on_toAdd = csc.data(control_start:control_start+length(this_tvec)-1);
        this_control_snip(iEvt,:) = this_on_toAdd';   
    end

    % Load the phases at stim_on in various frequency bands
    load('stim_phases.mat');
    % Get rid of all the 3rd band stuff IF there are 4 bands
    if size(causal_phase, 1) == 4 causal_phase(3,:) = []; end

    fig = figure('WindowState', 'Maximized');
    for iF = 1:length(fbands)
        subplot(3,4, (iF-1)*4+1)
        histogram(stim_phase(iF,:),phase_bins,'FaceColor',c_list{iF});
        title(sprintf('Hilbert at Stim %d - % d Hz',fbands{iF}(1), fbands{iF}(2)));

        subplot(3,4, (iF-1)*4+2)
        histogram(control_phase(iF,:),phase_bins,'FaceColor',c_list{iF});
        title(sprintf('Hilbert at Control %d - % d Hz',fbands{iF}(1), fbands{iF}(2)));

        subplot(3,4, (iF-1)*4+3)
        histogram(causal_phase(iF,:),phase_bins,'FaceColor',c_list{iF});
        title(sprintf('Causal at Stim %d - % d Hz',fbands{iF}(1), fbands{iF}(2)));
    end

    subplot(3,4,[4,8,12])
    plot(this_tvec, mean(this_on_snip), 'cyan');
    hold on
    plot(this_tvec, mean(this_control_snip), 'magenta');
    xline(0, '--black')
    legend({'Stim','Control'})
    sgtitle(pwd, 'Interpreter', 'none');
    mkdir('STA_dump');
    cd('STA_dump');
    print(fig, '-dpng', 'stim_phase_dist' );
    close;
    % Code to find examples of artifacts
    first = 0;
    normfun = @(x) (x - min(x))/(max(x) - min(x));
    while first < length(stim_on)
        fig = figure('WindowState', 'maximized');
        for iP = 1:5
            if first + iP > length(stim_on)
                break
            end
            for iF = 1:3
                subplot(5,3,(iP-1)*3+iF)
                hold on
                plot(this_tvec, normfun(this_on_snip(first+iP,:)), 'Color', 'black');
                plot(this_tvec, squeeze(normfun(filt_snip(iF,first+iP,:))), 'Color', c_list{iF}, 'LineStyle', '--');
                plot(this_tvec, squeeze(normfun(phase_snip(iF,first+iP,:))), 'Color', c_list{iF});
                xline(0, '--black')
                xticks([-0.25 0 0.25])
                yticks([])
                if iP == 1
                    title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
                end
                if iF == 1
                    ylabel(sprintf('Stim #%d', first+iP))
                end
                if iF == 3
                    xticks([-0.1 0 0.1]);
                    xlim([-0.1 0.1]);
                end
            end
        end
        first = first+5;
        print(fig, '-dpng',  num2str(first));
        close;
    end       
end

%% Function to return optimal parameters for echt
