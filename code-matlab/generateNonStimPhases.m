%% Assumes that good LFPs have been picked out

top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
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
    fbands = {[2 5], [6 10], [30 55]};
    % we are using this window length because we don't 
    % see much of a difference in phase estimation 
    % TODO : Find a way to communicate this to the readers
    win_length  = 0.5; 

    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end

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
        % Load the stim_responses
        load(strcat(fn_prefix, '_stim_response.mat'));
        
        stim_on = od.trial_nonstim.on_events;
        nonstim_causal_phase = nan(length(fbands), length(stim_on));
        nEnds = nearest_idx3(stim_on, csc.tvec);

         for iB = 1:length(fbands)
            nStarts = nearest_idx3(stim_on- win_length, csc.tvec);
            for iS = 1:length(nStarts)
                this_echt = echt(csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
                this_phase = angle(this_echt);
                nonstim_causal_phase(iB,iS) = this_phase(end); % The last sample's phase
%                 % Diagnostic to check how are angles assigned (Uncomment to run)
%                 diag_fig = figure;
%                 ax1 = subplot(2,1,1);
%                 plot(abs(this_echt))
%                 ax2 = subplot(2,1,2);
%                 plot(angle(this_echt))
%                 hold on
%                 yline(0, 'red')
%                 yline(pi/2, 'green')
%                 yline(pi, 'black')
%                 linkaxes([ax1,ax2],'x')
%                 close(diag_fig)
            end
        end
        % Assume you are in the correct folder
        save(strcat(fn_prefix, '_non_stim_phases'),'nonstim_causal_phase'); % should add option to save in specified output dir
        close;  
    end
end