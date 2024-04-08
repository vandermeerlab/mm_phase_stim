%% Assumes that good LFPs have been picked out

top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};

for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        rng(491994) % Random seed set for each session
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff
    end
end

%
function doStuff
    % Declaring variables
    % Setting up parameters
    fbands = {[2 5], [6 10], [12 28], [30 55]};
    c_list = {'red', 'blue','magenta', 'cyan'};
    % we are using this window length because we don't 
    % see much of a difference in phase estimation 
    bin_width = 0.01; % seconds
    stim_win = 0.25; % seconds around stim to ignore spikes
    nbins = 25;
    ns_bins =  -pi:2*pi/nbins:pi;
    ns_ticks = 0.5*(ns_bins(1:end-1) + ns_bins(2:end));
    
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
    cfg.fc = ExpKeys.goodCell;
    S = LoadSpikes(cfg);

    cfg = [];
    cfg.fc = ExpKeys.goodLFP;
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

    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        this_cell = SelectTS([], S, iC);
        this_cell = restrict(this_cell, iv(ExpKeys.stim_times));
        
        % Get rid of all spikes that happen within 0.25 around stim
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

        [all_sc_phase, all_fr_phase] = deal(nan(length(stim_on), length(fbands), nbins));
        [all_phase,all_spk_phase] = deal(cell(size(fbands)));
        sanity_spk = 0; sanity_t = 0;
        for iStim = 2:length(stim_on)
            end_t = stim_on(iStim) - 0.01;
            start_t = stim_on(iStim-1) + stim_win;
            if end_t - start_t > 0.5
                this_spks = this_cell.t{1}((this_cell.t{1} >= start_t) & ...
                    (this_cell.t{1} < end_t));
                sanity_spk = sanity_spk + length(this_spks);
                sanity_t = sanity_t + end_t - start_t;
                end_idx = nearest_idx3(end_t, csc.tvec);
                start_idx = nearest_idx3(start_t, csc.tvec);
                actual_t = csc.tvec(start_idx:end_idx);
                [~,edges,~] = histcounts(actual_t,start_t:bin_width:end_t);
                t_center = 0.5*(edges(1:end-1)+edges(2:end));
                t_idx = nearest_idx3(t_center, actual_t);
                spk_idx = nearest_idx3(this_spks, t_center);
                real_idx = nearest_idx3(this_spks, actual_t);
                for iF = 1:length(fbands)
                    this_echt = echt(csc.data(start_idx:end_idx), fbands{iF}(1), fbands{iF}(2), Fs);
                    this_phase = angle(this_echt);
                    % Inverting the phases here because data was recorded with Input inverted
                    this_phase = -1*this_phase;
                    all_phase{iF} = [all_phase{iF}, this_phase];
                    all_spk_phase{iF} = [all_spk_phase{iF}, this_phase(real_idx)];
                    this_phase = this_phase(t_idx); % Only the relevant phases
                    [pcounts, ~, pbins] = histcounts(this_phase, ns_bins);
                    [sc_phase, fr_phase] = deal(zeros(length(fbands), nbins));
                    for iBin = 1:nbins
                        sc_phase(iF, iBin) =  sum(sum(spk_idx == find(pbins==iBin)));
                        fr_phase(iF, iBin) = sc_phase(iF, iBin)/(pcounts(iBin)*bin_width);
                    end
                    all_sc_phase(iStim,iF,:) = sc_phase(iF,:);
                    all_fr_phase(iStim,iF,:) = fr_phase(iF,:);
                end
            end
        end

        phase_out = [];
        phase_out.binned_spk_count = squeeze(sum(all_sc_phase,1,'omitnan'));
        phase_out.binned_fr = squeeze(mean(all_fr_phase,1,'omitnan'));
        % Saving the lookup
        phase_out.ns_lookup = phase_out.binned_fr;
        phase_out.ns_lookup(:,1:end-1) = diff(phase_out.binned_fr,[],2);
        phase_out.ns_lookup(:,end) = phase_out.binned_fr(:,1) - phase_out.binned_fr(:,end);
    
        phase_out.real_mfr = sanity_spk/sanity_t;
        phase_out.spk_count = sanity_spk;
        phase_out.all_spk_phase = all_spk_phase;
        phase_out.all_phase = all_phase;
    
    
    %     % Sanity plot: Uncomment to debug
    %     fig = figure('WindowState', 'maximized');
    %     for iF = 1:length(fbands)
    %         subplot(3,3,3*(iF-1)+1)
    %         bar(ns_ticks,phase_out.binned_spk_count(iF,:),1, c_list{iF});
    %         xlabel('Phase Bins')
    %         ylabel('Spike Count')
    %         title(sprintf('%d - %d Hz', fbands{iF}(1), fbands{iF}(2)))
    % 
    %         subplot(3,3,3*(iF-1)+2)
    %         bar(ns_ticks,phase_out.binned_fr(iF,:),1, c_list{iF});
    %         xlabel('Phase Bins')
    %         ylabel('Firing rate (Hz)')
    %         title(sprintf('Mean FR: %.2f Hz', mean(phase_out.binned_fr(iF,:))))
    % 
    %         subplot(3,3,3*(iF-1)+3)
    %         bar(ns_ticks,phase_out.ns_lookup(iF,:),1, c_list{iF});
    %         xlabel('Phase Bins')
    %         ylabel('{\Delta} FR')
    %         title('NS Lookup')
    %     end
    %     sgtitle(sprintf('Real Mean Firing rate: %.2f', phase_out.real_mfr));
    %     close(fig)
    %     
    
    % Assume you are in the correct folder
    save(strcat(fn_prefix, '_nonstim_spk_phases'),'phase_out'); % should add option to save in specified output dir
    end
end
