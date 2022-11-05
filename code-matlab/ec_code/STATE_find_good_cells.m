function STATE_find_good_cells(data_dir)
%% Simple function for collecting all of the "ExpKeys.goodCell" from each recording 



if nargin ~= 1
    data_dir = cd;
end

cd(data_dir)

sess_list = {};
d = dir;
d=d(~ismember({d.name},{'.','..', '._*'}));
for iSess = 1:length(d)
    if ~strcmp(d(iSess).name(1:2), '._')
        sess_list{end+1} = d(iSess).name;

    end
end

count = 0;
for iSess = 1:length(sess_list)
    if iSess  == 30
        continue
    else
        
        if strcmp(sess_list{iSess}(1:3), 'M13') || strcmp(sess_list{iSess}(1:3), 'M14') %|| strcmp(sess_list{iSess}(1:3), 'M16')
            continue
        else
            cd(sess_list{iSess})
            
            LoadExpKeys
            
            for iCell = 1:length(ExpKeys.goodCell)
                
                fprintf('''%s_%s_%s'',...\n', strrep(sess_list{iSess}(1:14),'-','_'),strrep(num2str(ExpKeys.tetrodeDepths),'.','p'),  ExpKeys.goodCell{iCell}(1:end-2))
                count = count+1;
            end
            
            cd(data_dir)
            clear ExpKeys
            
        end
    end
end
    disp(count)

