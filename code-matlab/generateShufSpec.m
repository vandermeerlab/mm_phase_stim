%% Script to generate fake spike-phase locking data for significance testing (PLV version)
% Assumes that surrogate_sts.mat exists in each session's folder
rng(491994); % Setting the seed for reproducibility
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
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
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    nshufs = 1000;
    
    % Load the surrogate spectral stuff
    load('surrogate_sts.mat'); 

    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        fn_prefix = strrep(fn_prefix, '_', '-');
        % Load the spike_phaselock
        load(strcat(fn_prefix, '_spike_phaselock.mat'));
        [shuf_sts, shuf_ppc] = deal(zeros(nshufs, length(pool_sts.freq)));
        for iShuf = 1:nshufs
            pool_count = length(pool_sts.time{1});
            keep = randperm(pool_count);
            keep = keep(1:trial_clean_spk_count);
            this_shuf = pool_sts;
            this_shuf.label{1} = fn_prefix; % Because the pooled STS might have a differnt spike channel label
            this_shuf.fourierspctrm{1} = pool_sts.fourierspctrm{1}(keep,:,:);
            this_shuf.time{1} = pool_sts.time{1}(keep,:);
            this_shuf.trial{1} = pool_sts.trial{1}(keep,:);
            shuf_sts(iShuf,:) = nanmean(sq(abs(this_shuf.fourierspctrm{1})));
    
            cfg               = [];
            cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
            cfg.spikechannel  = this_shuf.label;
            cfg.channel       = this_shuf.lfplabel; % selected LFP channels
            cfg.avgoverchan   = 'weighted';
            this_shuf_ppc     = ft_spiketriggeredspectrum_stat(cfg, this_shuf);
            shuf_ppc(iShuf,:) = this_shuf_ppc.ppc0;
        end
        save(strcat(fn_prefix, '_shuf_spec'), 'shuf_sts', 'shuf_ppc', 'nshufs');
    end
end