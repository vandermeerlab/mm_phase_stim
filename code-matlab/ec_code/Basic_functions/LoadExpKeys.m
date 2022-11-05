%% LoadExpkeys: loading function for the 'Keys.m' files.  Keys files are
% used to track all experimental variables like subject, tragets, behaviour...
%
%
%
%   Inputs:
%       - none: will look for anything ending in '*Keys.m'
%      
%   Outputs:
%       - none.  Instead it will run the *Keys.m file resulting in the Keys
%       variable.
%
%
%   EC 2019-12-30: initial version based on ExpKeys/LoadExpKeys/FindFile by
%   MvdM and AD Redish
%   v1.0 can now deal with temporary or '._' files.
%
%% extract inputs and run or

    file_name = dir('*Keys.m');
    if isempty(file_name)
        error('No *Keys.m files found here')
    elseif length(file_name) >1
%         fprintf('%d *Keys.m files found. Finding the correct file\n', length(file_name))
        for iF = 1:length(file_name)
            if ~strcmp(file_name(iF).name(1:2), '._') && ~strcmp(file_name(iF).name(end-1:end), '~')
                good_ExpKeys(iF) = 1;
%                 fprintf('<strong>%s</strong> accepted\n', file_name(iF).name);
            else
                good_ExpKeys(iF) = 0;
%                 fprintf('<strong>%s</strong> rejected\n', file_name(iF).name);
            end
        end
        if sum(good_ExpKeys) > 1
            error('Too many *Keys.m" files in this dir')
        else
            run(file_name(find(good_ExpKeys==1)).name)
%             fprintf('\nKeys file: <strong>%s</strong> loaded.\n',file_name(find(good_ExpKeys==1)).name);
        end
        
    else
        run(file_name.name);
%             fprintf('\nKeys file: <strong>%s</strong> loaded.\n',file_name.name);
    end
    
clear good_ExpKeys file_name iF