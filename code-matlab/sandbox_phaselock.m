%% Script to generate spike-LFP phase locking

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

function doStuff
    LoadExpKeys;
    evs = LoadEvents([]);
    if isempty(ExpKeys.goodCell)
        return
    end
    
    cfg = []; cfg.fc = ExpKeys.goodLFP;
    if contains(cfg.fc, '-')
        temp = split(cfg.fc,'-');
        cfg.fc = {cat(2,temp{1},'.ncs')};
        csc = LoadCSC(cfg);
        ft_csc = ft_read_neuralynx_interp(cfg.fc);
        cfg_temp.fc = {cat(2,temp{2},'.ncs')};
        ref = LoadCSC(cfg_temp);
        csc.data = csc.data - ref.data;
        ft_ref = ft_read_neuralynx_interp(cfg_temp.fc);
        ft_csc.trial{1} = ft_csc.trial{1} - ft_ref.trial{1};
        clear temp ref ft_ref;
    else
        csc = LoadCSC(cfg);
        ft_csc = ft_read_neuralynx_interp(cfg.fc);
    end
    % Because of the way ft_read_data has been rewritten, the sampleinfo
    % field needs to be changed accordingly
    if ft_csc.sampleinfo(2) ~= length(ft_csc.time{1})
        fprintf('Adjusting sample info manually to avoid running into problems later\n');
        ft_csc.sampleinfo(2) = length(ft_csc.time{1});
    end
    overall_max = csc.cfg.hdr{1}.InputRange * 1e-6;
    overall_min = -overall_max;
    all_saturated = (csc.data(1,:) == overall_max | csc.data(1,:) == overall_min);
    all_saturated = csc.tvec(all_saturated);
    % The units of ft_csc are not the same as LoadCSC
    
    % Modify ft_csc to enable easy restricting
    temp_tvec = [0:length(ft_csc.time{1})-1];
    temp_offset = (double(ft_csc.hdr.LastTimeStamp)/1e6 - double(ft_csc.hdr.FirstTimeStamp)/1e6)/(length(temp_tvec) - 1);
    temp_tvec = temp_tvec * temp_offset;
     
    ft_csc.time{1} = temp_tvec;
    ft_csc.fsample = 1/temp_offset;
    % clean_spikes
    cfg = []; cfg.fc = ExpKeys.goodCell;
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
    S = LoadSpikes(cfg);
    S.ft_spikes = cell(size(S.t));
    
    % Remove spikes during short stim times
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
    
    long_stim_on = evs.t{strcmp(evs.label, ExpKeys.long_stim_on)} + start_delay;
    if ~isempty(long_stim_on) && ~isempty(ExpKeys.long_stim_times)
        long_stim_on = long_stim_on(long_stim_on >= ExpKeys.long_stim_times(1) & ...
                                    long_stim_on <= ExpKeys.long_stim_times(2));
    else
        long_stim_on = [];
    end
    
    all_stim = [pre_stim_on-0.5, pre_stim_on+0.5; stim_on-0.5,stim_on+0.5; ...
        post_stim_on-0.5,post_stim_on+0.5; long_stim_on-0.5,long_stim_on+0.5];
    
    clean_iv = MergeIV([], iv(all_stim));
    clean_iv = InvertIV(clean_iv, ExpKeys.recording_times(1), ExpKeys.recording_times(2));
    restricted_S = restrict(S, clean_iv);
    
    for iC = 1:length(S.t)
        this_fig = figure('WindowState', 'maximized');
        if ~strcmp(ExpKeys.experimenter, 'EC')
            S.ft_spikes{iC} =  ft_read_spike(S.label{iC}, 'encoding', 64);
        else
            S.ft_spikes{iC} =  ft_read_spike(S.label{iC}, 'encoding', 32);
        end
        % For clean_spikes, use restricted_S
        S.ft_clean_spikes{iC} = S.ft_spikes{iC};
        S.ft_clean_spikes{iC}.timestamp{1} = restricted_S.t{iC}*1e6; % TODO: Check for EC data too

        all_spike_data = ft_appendspike([], ft_csc, S.ft_spikes{iC});
        clean_spike_data = ft_appendspike([], ft_csc, S.ft_clean_spikes{iC});
    
        % Find the spike indices of the clean spikes among all_spikes
        all_idx = find(all_spike_data.trial{1}(2,:));
        clean_idx = find(clean_spike_data.trial{1}(2,:)); 
        clean_spk_idx = nearest_idx3(all_spike_data.time{1}(1,clean_idx), all_spike_data.time{1}(1,all_idx));

        % Find the largest Unsaturated block and calculate everything for
        % that and then subset
        % code to separate out all stim duration
        temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;
        big_start = ExpKeys.recording_times(1);
        big_stop = ExpKeys.recording_times(2);
        % Find the maximum unsaturated period in this epoch
        if ~isempty(all_saturated)
            this_saturated = (all_saturated > big_start) & (all_saturated < big_stop);
            this_saturated = all_saturated(this_saturated);
            if ~isempty(this_saturated)
                this_segments = [big_start; this_saturated; big_stop];
                [max_len, midx] = max(diff(this_segments));
                fprintf('Taking the longest unsaturated segment of %f sec out of total %f sec\n', ...
                    max_len, big_stop-big_start);
                big_start = this_segments(midx);
                big_stop = this_segments(midx+1);
            end
        end
        temp_start = nearest_idx3(big_start, temp_tvec);
        temp_end = nearest_idx3(big_stop, temp_tvec);
        % This temp_end has to be chosen carefully to avoid no Nan values at
        % the end of all_data.trial{1}(1,:)
        cfg=[];
        cfg.begsample = temp_start;
        cfg.endsample = temp_end;
        all_data = ft_redefinetrial(cfg, all_spike_data);
        all_clean_data = all_data;
        all_clean_data.trial{1}(2,:) = zeros(size(all_clean_data.trial{1}(2,:)));
        all_clean_data.trial{1}(2,clean_idx) = 1;
        
        all_spk_count = sum(all_data.trial{1}(2,:));
        all_clean_spk_count = sum(all_clean_data.trial{1}(2,:));
        
        if all_spk_count ~= 0
            % Calculate STA
            cfg = [];
            cfg.timwin = [-0.5 0.5];
            cfg.spikechannel = S.ft_spikes{iC}.label{1};
            cfg.channel = ft_csc.label(1);
            this_sta = ft_spiketriggeredaverage(cfg, all_data);
            all_sta = [];
            all_sta.time = this_sta.time;
            all_sta.vals = this_sta.avg(:,:)';
            if all_clean_spk_count ~= 0
                this_sta = ft_spiketriggeredaverage(cfg, all_clean_data);
                all_clean_sta.time = this_sta.time;
                all_clean_sta.vals = this_sta.avg(:,:)';
            else
                all_clean_sta = [];
            end
        
            % Calculate STS
            cfg = [];
            cfg.method = 'mtmconvol';
            cfg.foi = 1:1:100;
            cfg.t_ftimwin = 5./cfg.foi;
            cfg.taper = 'hanning';
            cfg.spikechannel = S.ft_spikes{iC}.label{1};
            cfg.channel = ft_csc.label{1};
            cfg.rejectsaturation = 'no'; % This is important because we are manually selecting the largest unsaturated gap
            big_sts = ft_spiketriggeredspectrum(cfg, all_data);
            all_sts.hasnan = ~isempty(find(isnan(big_sts.fourierspctrm{1}),1));
            all_sts.freqs = big_sts.freq;
            all_sts.vals = nanmean(sq(abs(big_sts.fourierspctrm{1})));
            if all_clean_spk_count ~= 0 
                % Instead of recalulating, reuse the values
                clean_sts = big_sts;
                clean_sts.fourierspctrm{1} = big_sts.fourierspctrm{1}(clean_spk_idx,:,:);
                clean_sts.time{1} = big_sts.time{1}(clean_spk_idx,:);
                clean_sts.trial{1} = big_sts.trial{1}(clean_spk_idx,:);
                all_clean_sts.hasnan = ~isempty(find(isnan(clean_sts.fourierspctrm{1}),1));
                all_clean_sts.freqs = clean_sts.freq;
                all_clean_sts.vals = nanmean(sq(abs(clean_sts.fourierspctrm{1})));
            else
                all_clean_sts = [];
            end
        
            % Calculate PPC
            cfg               = [];
            cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
            cfg.spikechannel  = big_sts.label;
            cfg.channel       = big_sts.lfplabel; % selected LFP channels
            cfg.avgoverchan   = 'weighted';
            cfg.timwin        = 'all'; % compute over all available spikes in the window
            this_ppc          = ft_spiketriggeredspectrum_stat(cfg,big_sts);
            all_ppc.hasnan = isempty(find(isnan(this_ppc.ppc0),1));
            all_ppc.vals = this_ppc.ppc0';
            if all_clean_spk_count ~= 0
                this_ppc = ft_spiketriggeredspectrum_stat(cfg, clean_sts);
                all_clean_ppc.hasnan = isempty(find(isnan(this_ppc.ppc0),1));
                all_clean_ppc.vals = this_ppc.ppc0';
            else
                all_clean_sts = [];
            end
        else
            all_sta = [];
            all_clean_sta = [];
            all_sts = [];
            all_clean_sts = [];
            all_ppc =  [];
            all_clean_ppc = [];
        end
        
        subplot(3,4,4)
        hold off
        plot(all_sta.time, all_sta.vals)
        if all_clean_spk_count ~= 0
            hold on
            plot(all_sta.time, all_clean_sta.vals)
        end
        xlabel('Time')
        ylabel('STA')
        yticks([])
        title('Whole Recording')
        
        subplot(3,4,8)
        hold off
        plot(all_sts.freqs, all_sts.vals)
        if all_clean_spk_count ~= 0
            hold on
            plot(all_sts.freqs, all_clean_sts.vals)
        end
        xlim([0 100])
        xlabel('Freqs')
        yticks([])
        ylabel('STS')
        
        ax = subplot(3,4,12);
        hold off
        plot(all_sts.freqs, all_ppc.vals)
        if all_clean_spk_count ~= 0
            hold on
            plot(all_sts.freqs, all_clean_ppc.vals)
            legend({sprintf('All: %d', all_spk_count), sprintf('Clean: %d', all_clean_spk_count)}, 'FontSize', 12, 'Location','best')
        else
            legend({sprintf('All: %d', all_spk_count)}, 'FontSize', 12, 'Location','best')
        end
        xlim([0 100])
        xlabel('Freqs')
        ylabel('PPC')
        ax.YAxis.Exponent = 0;
    
%         clear all_data all_clean_data temp_tvec % to avoid running out of space

        % code to separate out pre-stim duration
        this_start = max(ExpKeys.recording_times(1), big_start);
        this_stop = min(ExpKeys.pre_stim_times(2), big_stop);
        if this_start >= this_stop 
            % The longest overall unsaturated segment does not include this epoch, so find
            % the longest such segment in this epoch and recalculate stuff
            this_start = ExpKeys.recording_times(1);
            this_stop = ExpKeys.pre_stim_times(2);
            if ~isempty(all_saturated)
                this_saturated = (all_saturated > this_start) & (all_saturated < this_stop);
                this_saturated = all_saturated(this_saturated);
                if ~isempty(this_saturated)
                    this_segments = [this_start; this_saturated; this_stop];
                    [max_len, midx] = max(diff(this_segments));
                    fprintf('Taking the longest unsaturated segment of %f sec out of total %f sec\n', ...
                        max_len, this_stop-this_start);
                    this_start = this_segments(midx);
                    this_stop = this_segments(midx+1);
                end
            end
            temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;
            temp_start = nearest_idx3(this_start, temp_tvec);
            temp_end = nearest_idx3(this_stop, temp_tvec);
            cfg=[];
            cfg.begsample = temp_start;
            cfg.endsample = temp_end;
            pre_data = ft_redefinetrial(cfg, all_spike_data);
            pre_clean_data = ft_redefinetrial(cfg, clean_spike_data);

            pre_spk_count = sum(pre_data.trial{1}(2,:));
            pre_clean_spk_count = sum(pre_clean_data.trial{1}(2,:));
        
            if pre_spk_count ~= 0
                % Calculate STA
                cfg = [];
                cfg.timwin = [-0.5 0.5];
                cfg.spikechannel = S.ft_spikes{iC}.label{1};
                cfg.channel = ft_csc.label(1);
                this_sta = ft_spiketriggeredaverage(cfg, pre_data);
                pre_sta = [];
                pre_sta.time = this_sta.time;
                pre_sta.vals = this_sta.avg(:,:)';
                if pre_clean_spk_count ~= 0
                    this_sta = ft_spiketriggeredaverage(cfg, pre_clean_data);
                    pre_clean_sta.time = this_sta.time;
                    pre_clean_sta.vals = this_sta.avg(:,:)';
                else
                    pre_clean_sta = [];
                end
            
                % Calculate STS
                cfg = [];
                cfg.method = 'mtmconvol';
                cfg.foi = 1:1:100;
                cfg.t_ftimwin = 5./cfg.foi;
                cfg.taper = 'hanning';
                cfg.spikechannel = S.ft_spikes{iC}.label{1};
                cfg.channel = ft_csc.label{1};
                cfg.rejectsaturation = 'no'; % This is important because we are manually selecting the largest unsaturated gap
                this_sts = ft_spiketriggeredspectrum(cfg, pre_data);
                pre_sts.hasnan = ~isempty(find(isnan(this_sts.fourierspctrm{1}),1));
                pre_sts.freqs = this_sts.freq;
                pre_sts.vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                if pre_clean_spk_count ~= 0
                    clean_sts = ft_spiketriggeredspectrum(cfg, pre_clean_data);
                    pre_clean_sts.hasnan = ~isempty(find(isnan(clean_sts.fourierspctrm{1}),1));
                    pre_clean_sts.freqs = clean_sts.freq;
                    pre_clean_sts.vals = nanmean(sq(abs(clean_sts.fourierspctrm{1})));
                else
                    pre_clean_sts = [];
                end
            
                % Calculate PPC
                cfg               = [];
                cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg.spikechannel  = this_sts.label;
                cfg.channel       = this_sts.lfplabel; % selected LFP channels
                cfg.avgoverchan   = 'weighted';
                cfg.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc          = ft_spiketriggeredspectrum_stat(cfg,this_sts);
                pre_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                pre_ppc.vals = this_ppc.ppc0';
                if pre_clean_spk_count ~= 0
                    this_ppc = ft_spiketriggeredspectrum_stat(cfg, clean_sts);
                    pre_clean_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                    pre_clean_ppc.vals = this_ppc.ppc0';
                else
                    pre_clean_sts = [];
                end
            else
                pre_sta = [];
                pre_clean_sta = [];
                pre_sts = [];
                pre_clean_sts = [];
                pre_ppc =  [];
                pre_clean_ppc = [];
            end

        else % Use clever indexing
            temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;
            temp_start = nearest_idx3(this_start, temp_tvec);
            temp_end = nearest_idx3(this_stop, temp_tvec);
            cfg=[];
            cfg.begsample = temp_start;
            cfg.endsample = temp_end;
            pre_data = ft_redefinetrial(cfg, all_spike_data);
            pre_clean_data = ft_redefinetrial(cfg, clean_spike_data);
            pre_spk_count = sum(pre_data.trial{1}(2,:));
            pre_clean_spk_count = sum(pre_clean_data.trial{1}(2,:));
            
            if pre_spk_count ~= 0
                % Calculate STA
                cfg = [];
                cfg.timwin = [-0.5 0.5];
                cfg.spikechannel = S.ft_spikes{iC}.label{1};
                cfg.channel = ft_csc.label(1);
                this_sta = ft_spiketriggeredaverage(cfg, pre_data);
                pre_sta = [];
                pre_sta.time = this_sta.time;
                pre_sta.vals = this_sta.avg(:,:)';
                if pre_clean_spk_count ~= 0
                    this_sta = ft_spiketriggeredaverage(cfg, pre_clean_data);
                    pre_clean_sta.time = this_sta.time;
                    pre_clean_sta.vals = this_sta.avg(:,:)';
                else
                    pre_clean_sta = [];
                end
                % Find where the spikes and 'clean' spikes in this epoch lie among all spikes
                this_spk_idx  = nearest_idx3(find(all_data.trial{1}(2,pre_data.sampleinfo(1):pre_data.sampleinfo(2))) + ...
                    pre_data.sampleinfo(1) - 1,  find(all_data.trial{1}(2,:)));
                assert(length(this_spk_idx) == length(unique(this_spk_idx)), 'Clever indexing failed for pre stim epoch');
                this_clean_spk_idx = nearest_idx3(find(all_clean_data.trial{1}(2,pre_data.sampleinfo(1):pre_data.sampleinfo(2))) + ...
                    pre_data.sampleinfo(1) - 1, find(all_data.trial{1}(2,:)));
                assert(length(this_clean_spk_idx) == length(unique(this_clean_spk_idx)), 'Clever indexing failed for pre stim epoch');
                this_sts = big_sts;
                this_sts.fourierspctrm{1} = big_sts.fourierspctrm{1}(this_spk_idx,:,:);
                this_sts.time{1} = big_sts.time{1}(this_spk_idx,:);
                this_sts.trial{1} = big_sts.trial{1}(this_spk_idx,:);
                this_sts.trialtime(1) = (pre_data.sampleinfo(1) - all_data.sampleinfo(1))'/diff(all_data.sampleinfo)*big_sts.trialtime(2);
                this_sts.trialtime(2) = (pre_data.sampleinfo(2) - all_data.sampleinfo(1))'/diff(all_data.sampleinfo)*big_sts.trialtime(2);
                pre_sts.hasnan = ~isempty(find(isnan(this_sts.fourierspctrm{1}),1));
                pre_sts.freqs = this_sts.freq;
                pre_sts.vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                if pre_clean_spk_count ~= 0
                    clean_sts = big_sts;
                    clean_sts.fourierspctrm{1} = big_sts.fourierspctrm{1}(this_clean_spk_idx,:,:);
                    clean_sts.time{1} = big_sts.time{1}(this_clean_spk_idx,:);
                    clean_sts.trial{1} = big_sts.trial{1}(this_clean_spk_idx,:);
                    clean_sts.trialtime = this_sts.trialtime;
                    pre_clean_sts.hasnan = ~isempty(find(isnan(clean_sts.fourierspctrm{1}),1));
                    pre_clean_sts.freqs = clean_sts.freq;
                    pre_clean_sts.vals = nanmean(sq(abs(clean_sts.fourierspctrm{1})));
                else
                    pre_clean_sts = [];
                end           
           
                % Calculate PPC
                cfg               = [];
                cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg.spikechannel  = this_sts.label;
                cfg.channel       = this_sts.lfplabel; % selected LFP channels
                cfg.avgoverchan   = 'weighted';
                cfg.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc          = ft_spiketriggeredspectrum_stat(cfg,this_sts);
                pre_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                pre_ppc.vals = this_ppc.ppc0';
                if pre_clean_spk_count ~= 0
                    this_ppc = ft_spiketriggeredspectrum_stat(cfg, clean_sts);
                    pre_clean_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                    pre_clean_ppc.vals = this_ppc.ppc0';
                else
                    pre_clean_sts = [];
                end
            else
                pre_sta = [];
                pre_clean_sta = [];
                pre_sts = [];
                pre_clean_sts = [];
                pre_ppc =  [];
                pre_clean_ppc = [];
            end
        end

        subplot(3,4,1)
        hold off
        plot(pre_sta.time, pre_sta.vals)
        if pre_clean_spk_count ~= 0
            hold on
            plot(pre_sta.time, pre_clean_sta.vals)
        end
        xlabel('Time')
        ylabel('STA')
        yticks([])
        title('Pre-stim')
        
        subplot(3,4,5)
        hold off
        plot(pre_sts.freqs, pre_sts.vals)
        if pre_clean_spk_count ~= 0
            hold on
            plot(pre_sts.freqs, pre_clean_sts.vals)
        end
        xlim([0 100])
        xlabel('Freqs')
        yticks([])
        ylabel('STS')
        
        ax = subplot(3,4,9);
        hold off
        plot(pre_sts.freqs, pre_ppc.vals)
        if pre_clean_spk_count ~= 0
            hold on
            plot(pre_sts.freqs, pre_clean_ppc.vals)
            legend({sprintf('All: %d', pre_spk_count), sprintf('Clean: %d', pre_clean_spk_count)}, 'FontSize', 12, 'Location','best')
        else
            legend({sprintf('All: %d', pre_spk_count)}, 'FontSize', 12, 'Location','best')
        end
        xlim([0 100])
        xlabel('Freqs')
        ylabel('PPC')
        ax.YAxis.Exponent = 0;

        clear pre_data pre_clean_data temp_tvec % to avoid running out of space
    
        % code to separate out trial-stim duration
        this_start = max(ExpKeys.stim_times(1), big_start);
        this_stop = min(ExpKeys.stim_times(2), big_stop);
        if this_start >= this_stop 
            % The longest overall unsaturated segment does not include this epoch, so find
            % the longest such segment in this epoch and recalculate stuff
            this_start = ExpKeys.stim_times(1);
            this_stop = ExpKeys.stim_times(2);
            if ~isempty(all_saturated)
                this_saturated = (all_saturated > this_start) & (all_saturated < this_stop);
                this_saturated = all_saturated(this_saturated);
                if ~isempty(this_saturated)
                    this_segments = [this_start; this_saturated; this_stop];
                    [max_len, midx] = max(diff(this_segments));
                    fprintf('Taking the longest unsaturated segment of %f sec out of total %f sec\n', ...
                        max_len, this_stop-this_start);
                    this_start = this_segments(midx);
                    this_stop = this_segments(midx+1);
                end
            end
            temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;
            temp_start = nearest_idx3(this_start, temp_tvec);
            temp_end = nearest_idx3(this_stop, temp_tvec);
            cfg=[];
            cfg.begsample = temp_start;
            cfg.endsample = temp_end;
            trial_data = ft_redefinetrial(cfg, all_spike_data);
            trial_clean_data = ft_redefinetrial(cfg, clean_spike_data);

            trial_spk_count = sum(trial_data.trial{1}(2,:));
            trial_clean_spk_count = sum(trial_clean_data.trial{1}(2,:));
        
            if trial_spk_count ~= 0
                % Calculate STA
                cfg = [];
                cfg.timwin = [-0.5 0.5];
                cfg.spikechannel = S.ft_spikes{iC}.label{1};
                cfg.channel = ft_csc.label(1);
                this_sta = ft_spiketriggeredaverage(cfg, trial_data);
                trial_sta = [];
                trial_sta.time = this_sta.time;
                trial_sta.vals = this_sta.avg(:,:)';
                if trial_clean_spk_count ~= 0
                    this_sta = ft_spiketriggeredaverage(cfg, trial_clean_data);
                    trial_clean_sta.time = this_sta.time;
                    trial_clean_sta.vals = this_sta.avg(:,:)';
                else
                    trial_clean_sta = [];
                end
            
                % Calculate STS
                cfg = [];
                cfg.method = 'mtmconvol';
                cfg.foi = 1:1:100;
                cfg.t_ftimwin = 5./cfg.foi;
                cfg.taper = 'hanning';
                cfg.spikechannel = S.ft_spikes{iC}.label{1};
                cfg.channel = ft_csc.label{1};
                cfg.rejectsaturation = 'no'; % This is important because we are manually selecting the largest unsaturated gap
                this_sts = ft_spiketriggeredspectrum(cfg, trial_data);
                trial_sts.hasnan = ~isempty(find(isnan(this_sts.fourierspctrm{1}),1));
                trial_sts.freqs = this_sts.freq;
                trial_sts.vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                if trial_clean_spk_count ~= 0
                    clean_sts = ft_spiketriggeredspectrum(cfg, trial_clean_data);
                    trial_clean_sts.hasnan = ~isempty(find(isnan(clean_sts.fourierspctrm{1}),1));
                    trial_clean_sts.freqs = clean_sts.freq;
                    trial_clean_sts.vals = nanmean(sq(abs(clean_sts.fourierspctrm{1})));
                else
                    trial_clean_sts = [];
                end
            
                % Calculate PPC
                cfg               = [];
                cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg.spikechannel  = this_sts.label;
                cfg.channel       = this_sts.lfplabel; % selected LFP channels
                cfg.avgoverchan   = 'weighted';
                cfg.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc          = ft_spiketriggeredspectrum_stat(cfg,this_sts);
                trial_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                trial_ppc.vals = this_ppc.ppc0';
                if trial_clean_spk_count ~= 0
                    this_ppc = ft_spiketriggeredspectrum_stat(cfg, clean_sts);
                    trial_clean_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                    trial_clean_ppc.vals = this_ppc.ppc0';
                else
                    trial_clean_sts = [];
                end
            else
                trial_sta = [];
                trial_clean_sta = [];
                trial_sts = [];
                trial_clean_sts = [];
                trial_ppc =  [];
                trial_clean_ppc = [];
            end

        else % Use clever indexing
            temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;
            temp_start = nearest_idx3(this_start, temp_tvec);
            temp_end = nearest_idx3(this_stop, temp_tvec);
            cfg=[];
            cfg.begsample = temp_start;
            cfg.endsample = temp_end;
            trial_data = ft_redefinetrial(cfg, all_spike_data);
            trial_clean_data = ft_redefinetrial(cfg, clean_spike_data);
            trial_spk_count = sum(trial_data.trial{1}(2,:));
            trial_clean_spk_count = sum(trial_clean_data.trial{1}(2,:));
            
            if trial_spk_count ~= 0
                % Calculate STA
                cfg = [];
                cfg.timwin = [-0.5 0.5];
                cfg.spikechannel = S.ft_spikes{iC}.label{1};
                cfg.channel = ft_csc.label(1);
                this_sta = ft_spiketriggeredaverage(cfg, trial_data);
                trial_sta = [];
                trial_sta.time = this_sta.time;
                trial_sta.vals = this_sta.avg(:,:)';
                if trial_clean_spk_count ~= 0
                    this_sta = ft_spiketriggeredaverage(cfg, trial_clean_data);
                    trial_clean_sta.time = this_sta.time;
                    trial_clean_sta.vals = this_sta.avg(:,:)';
                else
                    trial_clean_sta = [];
                end
                % Find where the spikes and 'clean' spikes in this epoch lie among all spikes
                this_spk_idx  = nearest_idx3(find(all_data.trial{1}(2,trial_data.sampleinfo(1):trial_data.sampleinfo(2))) + ...
                    trial_data.sampleinfo(1) - 1,  find(all_data.trial{1}(2,:)));
                assert(length(this_spk_idx) == length(unique(this_spk_idx)), 'Clever indexing failed for trial stim epoch');
                this_clean_spk_idx = nearest_idx3(find(all_clean_data.trial{1}(2,trial_data.sampleinfo(1):trial_data.sampleinfo(2))) + ...
                    trial_data.sampleinfo(1) - 1, find(all_data.trial{1}(2,:)));
                assert(length(this_clean_spk_idx) == length(unique(this_clean_spk_idx)), 'Clever indexing failed for trial stim epoch');
                this_sts = big_sts;
                this_sts.fourierspctrm{1} = big_sts.fourierspctrm{1}(this_spk_idx,:,:);
                this_sts.time{1} = big_sts.time{1}(this_spk_idx,:);
                this_sts.trial{1} = big_sts.trial{1}(this_spk_idx,:);
                this_sts.trialtime(1) = (trial_data.sampleinfo(1) - all_data.sampleinfo(1))'/diff(all_data.sampleinfo)*big_sts.trialtime(2);
                this_sts.trialtime(2) = (trial_data.sampleinfo(2) - all_data.sampleinfo(1))'/diff(all_data.sampleinfo)*big_sts.trialtime(2);
                trial_sts.hasnan = ~isempty(find(isnan(this_sts.fourierspctrm{1}),1));
                trial_sts.freqs = this_sts.freq;
                trial_sts.vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                if trial_clean_spk_count ~= 0
                    clean_sts = big_sts;
                    clean_sts.fourierspctrm{1} = big_sts.fourierspctrm{1}(this_clean_spk_idx,:,:);
                    clean_sts.time{1} = big_sts.time{1}(this_clean_spk_idx,:);
                    clean_sts.trial{1} = big_sts.trial{1}(this_clean_spk_idx,:);
                    clean_sts.trialtime = this_sts.trialtime;
                    trial_clean_sts.hasnan = ~isempty(find(isnan(clean_sts.fourierspctrm{1}),1));
                    trial_clean_sts.freqs = clean_sts.freq;
                    trial_clean_sts.vals = nanmean(sq(abs(clean_sts.fourierspctrm{1})));
                else
                    trial_clean_sts = [];
                end           
           
                % Calculate PPC
                cfg               = [];
                cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg.spikechannel  = this_sts.label;
                cfg.channel       = this_sts.lfplabel; % selected LFP channels
                cfg.avgoverchan   = 'weighted';
                cfg.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc          = ft_spiketriggeredspectrum_stat(cfg,this_sts);
                trial_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                trial_ppc.vals = this_ppc.ppc0';
                if trial_clean_spk_count ~= 0
                    this_ppc = ft_spiketriggeredspectrum_stat(cfg, clean_sts);
                    trial_clean_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                    trial_clean_ppc.vals = this_ppc.ppc0';
                else
                    trial_clean_sts = [];
                end
            else
                trial_sta = [];
                trial_clean_sta = [];
                trial_sts = [];
                trial_clean_sts = [];
                trial_ppc =  [];
                trial_clean_ppc = [];
            end
        end

        subplot(3,4,2)
        hold off
        plot(trial_sta.time, trial_sta.vals)
        if trial_clean_spk_count ~= 0
            hold on
            plot(trial_sta.time, trial_clean_sta.vals)
        end
        xlabel('Time')
        ylabel('STA')
        yticks([])
        title('Trial-stim')
        
        subplot(3,4,6)
        hold off
        plot(trial_sts.freqs, trial_sts.vals)
        if trial_clean_spk_count ~= 0
            hold on
            plot(trial_sts.freqs, trial_clean_sts.vals)
        end
        xlim([0 100])
        xlabel('Freqs')
        yticks([])
        ylabel('STS')
        
        ax=subplot(3,4,10);
        hold off
        plot(trial_sts.freqs, trial_ppc.vals)
        if trial_clean_spk_count ~= 0
            hold on
            plot(trial_sts.freqs, trial_clean_ppc.vals)
            legend({sprintf('All: %d', trial_spk_count), sprintf('Clean: %d', trial_clean_spk_count)}, 'FontSize', 12, 'Location','best')
        else
            legend({sprintf('All: %d', trial_spk_count)}, 'FontSize', 12, 'Location','best')
        end
        xlim([0 100])
        xlabel('Freqs')
        ylabel('PPC')
        ax.YAxis.Exponent = 0;

        clear trial_data trial_clean_data temp_tvec % to avoid running out of space
        
        if ~isempty(ExpKeys.post_baseline_times)
            % code to separate out post-stim duration
            this_start = max(ExpKeys.post_baseline_times(1), big_start);
            this_stop = min(ExpKeys.recording_times(2), big_stop);
            if this_start >= this_stop 
                % The longest overall unsaturated segment does not include this epoch, so find
                % the longest such segment in this epoch and recalculate stuff
                this_start = ExpKeys.post_baseline_times(1);
                this_stop = ExpKeys.recording_times(2);
                if ~isempty(all_saturated)
                    this_saturated = (all_saturated > this_start) & (all_saturated < this_stop);
                    this_saturated = all_saturated(this_saturated);
                    if ~isempty(this_saturated)
                        this_segments = [this_start; this_saturated; this_stop];
                        [max_len, midx] = max(diff(this_segments));
                        fprintf('Taking the longest unsaturated segment of %f sec out of total %f sec\n', ...
                            max_len, this_stop-this_start);
                        this_start = this_segments(midx);
                        this_stop = this_segments(midx+1);
                    end
                end
                temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;
                temp_start = nearest_idx3(this_start, temp_tvec);
                temp_end = nearest_idx3(this_stop, temp_tvec);
                cfg=[];
                cfg.begsample = temp_start;
                cfg.endsample = temp_end;
                post_data = ft_redefinetrial(cfg, all_spike_data);
                post_clean_data = ft_redefinetrial(cfg, clean_spike_data);
    
                post_spk_count = sum(post_data.trial{1}(2,:));
                post_clean_spk_count = sum(post_clean_data.trial{1}(2,:));
            
                if post_spk_count ~= 0
                    % Calculate STA
                    cfg = [];
                    cfg.timwin = [-0.5 0.5];
                    cfg.spikechannel = S.ft_spikes{iC}.label{1};
                    cfg.channel = ft_csc.label(1);
                    this_sta = ft_spiketriggeredaverage(cfg, post_data);
                    post_sta = [];
                    post_sta.time = this_sta.time;
                    post_sta.vals = this_sta.avg(:,:)';
                    if post_clean_spk_count ~= 0
                        this_sta = ft_spiketriggeredaverage(cfg, post_clean_data);
                        post_clean_sta.time = this_sta.time;
                        post_clean_sta.vals = this_sta.avg(:,:)';
                    else
                        post_clean_sta = [];
                    end
                
                    % Calculate STS
                    cfg = [];
                    cfg.method = 'mtmconvol';
                    cfg.foi = 1:1:100;
                    cfg.t_ftimwin = 5./cfg.foi;
                    cfg.taper = 'hanning';
                    cfg.spikechannel = S.ft_spikes{iC}.label{1};
                    cfg.channel = ft_csc.label{1};
                    cfg.rejectsaturation = 'no'; % This is important because we are manually selecting the largest unsaturated gap
                    this_sts = ft_spiketriggeredspectrum(cfg, post_data);
                    post_sts.hasnan = ~isempty(find(isnan(this_sts.fourierspctrm{1}),1));
                    post_sts.freqs = this_sts.freq;
                    post_sts.vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                    if post_clean_spk_count ~= 0
                        clean_sts = ft_spiketriggeredspectrum(cfg, post_clean_data);
                        post_clean_sts.hasnan = ~isempty(find(isnan(clean_sts.fourierspctrm{1}),1));
                        post_clean_sts.freqs = clean_sts.freq;
                        post_clean_sts.vals = nanmean(sq(abs(clean_sts.fourierspctrm{1})));
                    else
                        post_clean_sts = [];
                    end
                
                    % Calculate PPC
                    cfg               = [];
                    cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                    cfg.spikechannel  = this_sts.label;
                    cfg.channel       = this_sts.lfplabel; % selected LFP channels
                    cfg.avgoverchan   = 'weighted';
                    cfg.timwin        = 'all'; % compute over all available spikes in the window
                    this_ppc          = ft_spiketriggeredspectrum_stat(cfg,this_sts);
                    post_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                    post_ppc.vals = this_ppc.ppc0';
                    if post_clean_spk_count ~= 0
                        this_ppc = ft_spiketriggeredspectrum_stat(cfg, clean_sts);
                        post_clean_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                        post_clean_ppc.vals = this_ppc.ppc0';
                    else
                        post_clean_sts = [];
                    end
                else
                    post_sta = [];
                    post_clean_sta = [];
                    post_sts = [];
                    post_clean_sts = [];
                    post_ppc =  [];
                    post_clean_ppc = [];
                end
    
            else % Use clever indexing
                temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;
                temp_start = nearest_idx3(this_start, temp_tvec);
                temp_end = nearest_idx3(this_stop, temp_tvec);
                cfg=[];
                cfg.begsample = temp_start;
                cfg.endsample = temp_end;
                post_data = ft_redefinetrial(cfg, all_spike_data);
                post_clean_data = ft_redefinetrial(cfg, clean_spike_data);
                post_spk_count = sum(post_data.trial{1}(2,:));
                post_clean_spk_count = sum(post_clean_data.trial{1}(2,:));
                
                if post_spk_count ~= 0
                    % Calculate STA
                    cfg = [];
                    cfg.timwin = [-0.5 0.5];
                    cfg.spikechannel = S.ft_spikes{iC}.label{1};
                    cfg.channel = ft_csc.label(1);
                    this_sta = ft_spiketriggeredaverage(cfg, post_data);
                    post_sta = [];
                    post_sta.time = this_sta.time;
                    post_sta.vals = this_sta.avg(:,:)';
                    if post_clean_spk_count ~= 0
                        this_sta = ft_spiketriggeredaverage(cfg, post_clean_data);
                        post_clean_sta.time = this_sta.time;
                        post_clean_sta.vals = this_sta.avg(:,:)';
                    else
                        post_clean_sta = [];
                    end
                    % Find where the spikes and 'clean' spikes in this epoch lie among all spikes
                    this_spk_idx  = nearest_idx3(find(all_data.trial{1}(2,post_data.sampleinfo(1):post_data.sampleinfo(2))) + ...
                        post_data.sampleinfo(1) - 1,  find(all_data.trial{1}(2,:)));
                    assert(length(this_spk_idx) == length(unique(this_spk_idx)), 'Clever indexing failed for post stim epoch');
                    this_clean_spk_idx = nearest_idx3(find(all_clean_data.trial{1}(2,post_data.sampleinfo(1):post_data.sampleinfo(2))) + ...
                        post_data.sampleinfo(1) - 1, find(all_data.trial{1}(2,:)));
                    assert(length(this_clean_spk_idx) == length(unique(this_clean_spk_idx)), 'Clever indexing failed for post stim epoch');
                    this_sts = big_sts;
                    this_sts.fourierspctrm{1} = big_sts.fourierspctrm{1}(this_spk_idx,:,:);
                    this_sts.time{1} = big_sts.time{1}(this_spk_idx,:);
                    this_sts.trial{1} = big_sts.trial{1}(this_spk_idx,:);
                    this_sts.trialtime(1) = (post_data.sampleinfo(1) - all_data.sampleinfo(1))'/diff(all_data.sampleinfo)*big_sts.trialtime(2);
                    this_sts.trialtime(2) = (post_data.sampleinfo(2) - all_data.sampleinfo(1))'/diff(all_data.sampleinfo)*big_sts.trialtime(2);
                    post_sts.hasnan = ~isempty(find(isnan(this_sts.fourierspctrm{1}),1));
                    post_sts.freqs = this_sts.freq;
                    post_sts.vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                    if post_clean_spk_count ~= 0
                        clean_sts = big_sts;
                        clean_sts.fourierspctrm{1} = big_sts.fourierspctrm{1}(this_clean_spk_idx,:,:);
                        clean_sts.time{1} = big_sts.time{1}(this_clean_spk_idx,:);
                        clean_sts.trial{1} = big_sts.trial{1}(this_clean_spk_idx,:);
                        clean_sts.trialtime = this_sts.trialtime;
                        post_clean_sts.hasnan = ~isempty(find(isnan(clean_sts.fourierspctrm{1}),1));
                        post_clean_sts.freqs = clean_sts.freq;
                        post_clean_sts.vals = nanmean(sq(abs(clean_sts.fourierspctrm{1})));
                    else
                        post_clean_sts = [];
                    end           
               
                    % Calculate PPC
                    cfg               = [];
                    cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                    cfg.spikechannel  = this_sts.label;
                    cfg.channel       = this_sts.lfplabel; % selected LFP channels
                    cfg.avgoverchan   = 'weighted';
                    cfg.timwin        = 'all'; % compute over all available spikes in the window
                    this_ppc          = ft_spiketriggeredspectrum_stat(cfg,this_sts);
                    post_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                    post_ppc.vals = this_ppc.ppc0';
                    if post_clean_spk_count ~= 0
                        this_ppc = ft_spiketriggeredspectrum_stat(cfg, clean_sts);
                        post_clean_ppc.hasnan = ~isempty(find(isnan(this_ppc.ppc0),1));
                        post_clean_ppc.vals = this_ppc.ppc0';
                    else
                        post_clean_sts = [];
                    end
                else
                    post_sta = [];
                    post_clean_sta = [];
                    post_sts = [];
                    post_clean_sts = [];
                    post_ppc =  [];
                    post_clean_ppc = [];
                end
            end
    
            subplot(3,4,3)
            hold off
            plot(post_sta.time, post_sta.vals)
            if post_clean_spk_count ~= 0
                hold on
                plot(post_sta.time, post_clean_sta.vals)
            end
            xlabel('Time')
            ylabel('STA')
            yticks([])
            title('Post-stim')
            
            subplot(3,4,7)
            hold off
            plot(post_sts.freqs, post_sts.vals)
            if post_clean_spk_count ~= 0
                hold on
                plot(post_sts.freqs, post_clean_sts.vals)
            end
            xlim([0 100])
            xlabel('Freqs')
            yticks([])
            ylabel('STS')
            
            ax = subplot(3,4,11);
            hold off
            plot(post_sts.freqs, post_ppc.vals)
            if post_clean_spk_count ~= 0
                hold on
                plot(post_sts.freqs, post_clean_ppc.vals)
                legend({sprintf('All: %d', post_spk_count), sprintf('Clean: %d', post_clean_spk_count)}, 'FontSize', 12, 'Location','best')
            else
                legend({sprintf('All: %d', post_spk_count)}, 'FontSize', 12, 'Location','best')
            end
            xlim([0 100])
            xlabel('Freqs')
            ylabel('PPC')
            ax.YAxis.Exponent = 0;
    
            clear post_data post_clean_data temp_tvec % to avoid running out of space
        end
        
        clear all_data all_clean_data % to avoid running out of space

        sgtitle(S.label{iC}, 'Interpreter', 'None')
        
        % Save figure
        fn_prefix = extractBefore(S.label{iC}, '.t');
        fn_prefix = strrep(fn_prefix, '_', '-');
        print(this_fig, '-dpng', strcat(fn_prefix,'-SpikePhaseLock'));
        print(this_fig, '-dpng',  strcat('E:\Dropbox (Dartmouth College)\AnalysisResults\phase_stim_results\', fn_prefix, '-SpikePhaseLock'));
        % Save variables
        save(strcat(fn_prefix, '_spike_phaselock'), 'pre_spk_count', 'pre_clean_spk_count', 'pre_sta', 'pre_clean_sta', 'pre_sts', 'pre_clean_sts', 'pre_ppc', 'pre_clean_ppc', ...
            'trial_spk_count', 'trial_clean_spk_count', 'trial_sta', 'trial_clean_sta', 'trial_sts', 'trial_clean_sts', 'trial_ppc', 'trial_clean_ppc', ...
            'post_spk_count', 'post_clean_spk_count', 'post_sta', 'post_clean_sta', 'post_sts', 'post_clean_sts', 'post_ppc', 'post_clean_ppc', ...
            'all_spk_count', 'all_clean_spk_count', 'all_sta', 'all_clean_sta', 'all_sts', 'all_clean_sts', 'all_ppc', 'all_clean_ppc');
        close
    end
end
%% Raster to check if cleanup is working

% Pre-stim raster
% if ~isempty(ExpKeys.pre_stim_times)
%     fig2 = figure;
%     ax = subplot(1,3,1);
%     this_cell = SelectTS([], restricted_S, iC);
%     this_cell = SelectTS([], S, iC);
%     this_on_events = long_stim_on;
%     [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm([], ...
%         this_cell, this_on_events, '', 0.5);
%     hold on
%     plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
%     plot([0 0], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
%     plot([ExpKeys.short_stim_pulse_width+stop_delay ExpKeys.short_stim_pulse_width+stop_delay], [0 length(this_on_events)], 'color', 'red', 'linewidth', 1);
%     ylabel('Pre-Stim #')
%     ylim([1 length(this_on_events)])
%     xlim([-0.01 0.01]);
%     xlabel("Time (sec)")
% end