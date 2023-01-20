%% Script to generate spike-LFP phase locking
rng(491994); % Setting up seed for reproducibility
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
    % Declare parameters and variables
    fbands = {[2 5], [6 10], [30 55]};
    wsz = 0.5; % seconds, window to use for causal phase estimation
    min_spikes = 200;
    num_subamples = 1000;
 
    LoadExpKeys;
    evs = LoadEvents([]);
    if isempty(ExpKeys.goodCell)
        return
    end

    cfg = []; cfg.fc = ExpKeys.goodCell;
    if ~strcmp(ExpKeys.experimenter, 'EC')
        cfg.getRatings = 1;
        cfg.uint = '64';
    end
    S = LoadSpikes(cfg);

    for iC = 1:length(S.t)
        this_cell = SelectTS([], S, iC);
        fn_prefix = extractBefore(this_cell.label{1}, '.t');
        
        % Load non_stim lookup at various frequency bands
        load(strcat(fn_prefix, '_nonstim_spk_phases.mat'));     

        [trial_unsampled_mean_phase, trial_unsampled_plv, ...
            trial_subsampled_mean_phase, trial_subsampled_plv] = deal(nan(size(fbands)));
        trial_spk_phase = [];
        for iF = 1:length(fbands)
            trial_spk_phase(iF,:) = phase_out.all_spk_phase{iF};
            trial_unsampled_plv(iF) = resultantlength(trial_spk_phase(iF,:));
            trial_unsampled_mean_phase(iF) = circ_mean(trial_spk_phase(iF,:), [], 2);       
            if phase_out.spk_count > min_spikes % Subsample
                [temp_mean_phase, temp_plv] = deal(nan(1,num_subamples));
                for iSub = 1:num_subamples
                    this_idx = randperm(phase_out.spk_count);
                    temp_mean_phase(iSub) = circ_mean(trial_spk_phase(iF,this_idx(1:200)), [], 2);
                    temp_plv(iSub) = resultantlength(trial_spk_phase(iF,this_idx(1:200)));   
                end
                trial_subsampled_mean_phase(iF) = circ_mean(temp_mean_phase, [], 2); 
                trial_subsampled_plv(iF) = mean(temp_plv);
                clear temp_plv temp_mean_phase
            end
        end
        % Save variables
        save(strcat(fn_prefix, '_spike_phaselock_causal_plv'), ...
            'trial_unsampled_plv', 'trial_unsampled_mean_phase', ...
            'trial_subsampled_plv', 'trial_subsampled_mean_phase', ...
            'trial_spk_phase');
    end
end

%% Helper functions 
% Borrowed from FieldTrip
function [resLen] = resultantlength(angles)

n = sum(~isnan(angles));
resLen = abs(nansum(angles))./n;

end
