%% Assumes that good LFPs have been picked out

top_dir = 'data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};

for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        rng(491994);
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff
    end
end

%%
function doStuff
    
    % declaring variables
    nsamples = 1000;

    LoadExpKeys;
    evs = LoadEvents([]);
    cfg_spk = [];
    cfg_spk.fc = ExpKeys.goodCell;
    if isempty(cfg_spk.fc)
        return
    end
    
%     % For debugging
%     if ~(ExpKeys.hasWheelData) | contains(ExpKeys.light_source, 'LASER')
%         return
%     end

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
    fbands = {[2 5], [6 10], [30 55]}; % Hz
    all_win = [0.5, 0.75, 1, 1.25, 1.5]; % sec, size of windows to be used for ECHT
    
   % Include both pre and post baseline if available
   eval_csc{1} = restrict(csc, iv(ExpKeys.pre_baseline_times));
   true_phase = cell(length(fbands), 2);
   for iF = 1:length(fbands)
        cfg_filt = [];
        cfg_filt.type = 'fdesign';
        cfg_filt.f = fbands{iF};
        this_filt = FilterLFP(cfg_filt, eval_csc{1});
        true_phase{iF,1} = angle(hilbert(this_filt.data));
   end

    if ~isempty(ExpKeys.post_baseline_times)
        eval_csc{2} = restrict(csc, iv(ExpKeys.post_baseline_times));
        for iF = 1:length(fbands)
            cfg_filt = [];
            cfg_filt.type = 'fdesign';
            cfg_filt.f = fbands{iF};
            this_filt = FilterLFP(cfg_filt, eval_csc{2});
            true_phase{iF,2} = angle(hilbert(this_filt.data));
       end
    end
    eval_acausal_phase = nan(length(fbands), nsamples);
    eval_causal_phase = nan(length(all_win), length(fbands), nsamples);
    for iS = 1:nsamples
        for iF = 1:length(fbands)
            min_start = ceil(max(all_win)*Fs); % making sure the end remains the same for all cases
            % If post and pre both exist, use a coin toss to decide which one to pick from
            if ~isempty(ExpKeys.post_baseline_times)
                choice = rand();
                if choice < 0.5
                   choice = 1;
                else
                    choice = 2;
                end
            else
                choice = 1;
            end
            this_end = randi(length(eval_csc{choice}.tvec) - min_start) + min_start;
            eval_acausal_phase(iF,iS) = true_phase{iF,choice}(this_end);
            for iW = 1:length(all_win)
                wsz = all_win(iW);
                assert(eval_csc{choice}.tvec(this_end) - eval_csc{choice}.tvec(1) >= wsz); 
                this_start = nearest_idx3(eval_csc{choice}.tvec(this_end) - wsz, eval_csc{choice}.tvec);
                this_echt = echt(eval_csc{choice}.data(this_start:this_end), ...
                    fbands{iF}(1), fbands{iF}(2), Fs);
                this_phase = angle(this_echt);
                eval_causal_phase(iW,iF,iS) = this_phase(end); % The last sample's phase
            end
        end
    end
    save('echt_baseline_evaluation', 'eval_causal_phase', 'eval_acausal_phase');
end

