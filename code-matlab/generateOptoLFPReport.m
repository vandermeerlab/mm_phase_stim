%% Assumes that good LFPs have been picked out
% top_dir = 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\';
% mice = {'M016', 'M017', 'M018', 'M019', 'M020'};
top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff
    end
end

