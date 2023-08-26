%% Script to generate circularly shifted shuffles for signifcance testing of spike-phase locking (PLV version)
% Circularly shifting but not subsampling
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
    rng(491994); % Setting the seed for reproducibility (needs to be set for session at the start)
    LoadExpKeys;
    if isempty(ExpKeys.goodCell)
        return
    end
    nshufs = 1000; % This has to be less than min_spikes in this case
    min_spikes = 200;

    for iC = 1:length(ExpKeys.goodCell)
        fn_prefix = extractBefore(ExpKeys.goodCell{iC}, '.t');
        fn_prefix = strrep(fn_prefix, '_', '-');
        % Load the spike_phaselock data
        load(strcat(fn_prefix, '_spike_phaselock_plv.mat'));
        shuf_circ2_plv = nan(nshufs, (size(trial_spk_phase,1)));
        to_shift = randperm(trial_spk_count);
        if trial_spk_count >  min_spikes
            for iShuf = 1:nshufs
                % circularly shift spikes
                this_spk_phase = cellfun(@(x) circshift(x, to_shift(iShuf)), trial_spk_phase, 'UniformOutput', false);
                shuf_circ2_plv(iShuf,:) = cellfun(@(x) resultantlength(x), this_spk_phase);
            end
        end
        save(strcat(fn_prefix, '_shuf_spec_circ2_plv'), 'shuf_circ2_plv');
    end
end


%% Helper functions 
% Borrowed from FieldTrip
function [resLen] = resultantlength(angles)
    n = sum(~isnan(angles));
    resLen = abs(nansum(angles))./n;
end
