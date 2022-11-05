function [h, hStrings] = add_num_imagesc(h, vals, n_dec, font_size)
%% add_num_imagesc:  adds the numerical values of a matrix used to generate an imagesc.
%
%          Inputs:
%           - h: axis handle
%           - vals: matrix
%           - font_size: default is what ever the figure already had
%          Outputs:
%           - h: axis handle
%
% EC - 2016-10-31
if nargin <3
    n_dec = 3;
end
if nargin <4
    font_size = get(gca, 'fontsize');
end
if n_dec == 2
    textStrings = num2str(vals(:),'%0.2f');  %# Create strings from the matrix values
elseif n_dec == 3
    textStrings = num2str(vals(:),'%0.3f');  %# Create strings from the matrix values
elseif n_dec == 4
    textStrings = num2str(vals(:),'%0.4f');  %# Create strings from the matrix values
elseif n_dec == 5
    textStrings = num2str(vals(:),'%0.5f');  %# Create strings from the matrix values
end
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x, y]=meshgrid(1:size(vals,2), 1:size(vals,1));   %# Create x and y coordinates for the strings
flag =0; % flag if only contains NaNs
for ii = 1:length(textStrings)
    if ~strcmp(textStrings(ii), 'NaN')
        hStrings = text(x(ii),y(ii),textStrings(ii),'HorizontalAlignment','center');
        flag = 1; 
    end
end
if flag ~=0
set(hStrings,'fontsize',font_size)
else
    disp('Input vals only contain NaNs')
end