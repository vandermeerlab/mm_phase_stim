function S = LoadSpikes(tfilelist, encoding)

% adapted from M-clust function LoadSpikes


%-------------------
% Check input type
%-------------------
if ~isa(tfilelist, 'cell')
   error('LoadSpikes: tfilelist should be a cell array.');
end

nFiles = length(tfilelist);


S = cell(nFiles, 1);
for iF = 1:nFiles
	tfn = tfilelist{iF};
	if ~isempty(tfn)
		tfp = fopen(tfn, 'rb','b');
		if (tfp == -1)
			warning([ 'Could not open tfile ' tfn]);
		end
		
		ReadHeader(tfp);
        
        if encoding == 64
            S{iF} = fread(tfp,inf,'uint64');	% read as 64 bit ints
        else
            S{iF} = fread(tfp,inf,'uint32'); % read as 32 bit ints     
        end	
	    S{iF} = double(S{iF}*100);

		fclose(tfp);		
	end 		% if tfn valid
end		% for all files
fprintf(2,'\n');

