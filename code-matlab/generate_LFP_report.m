% Assumes that references and all LFP channels have the same sampling
% frequency

cd('');
LoadExpKeys;
evs = LoadEvents([]);

if strcmp(ExpKeys.isReferenceRecordedSeparately,'Yes')
    cfg_lfp.fc = ExpKeys.LFP_channels;
    all_lfp = LoadCSC(cfg_lfp);
    cfg_refs.fc = ExpKeys.ref_channels;
    all_refs = LoadCSC(cfg_refs);
    % Create reference subtracted versions 
    all_ref_labels = cell(length(cfg_lfp.fc) * length(cfg_refs.fc));
    all_ref_ = ;
end