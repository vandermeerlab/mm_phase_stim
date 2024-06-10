% Script to generate correlations between wideband signal
cd('E:\Dropbox (Dartmouth College)\manish_data\M321\M321-2022-07-13');
LoadExpKeys;
evs = LoadEvents([]);
%Ordered from shank A top to bottom, then Shank B top to bottom)
order1 = [3,16,5,12,9,8,11,6,7,2,1,4,13,14,10,15]; 
order2 =  [29,17,28,21,30,22,32,25,31,27,26,24,20,19,23,18];
cfg1.fc = cell(size(order1));
for i = 1:length(cfg1.fc)
    cfg1.fc{i} = strcat('CSC', num2str(order1(i)), '.ncs');
end
cfg2.fc = cell(size(order2));
for i = 1:length(cfg2.fc)
    cfg2.fc{i} = strcat('CSC', num2str(order2(i)), '.ncs');
end

%%
shank1_csc = LoadCSC(cfg1);

%%
cfg_tt4.fc = {'CSC13.ncs', 'CSC14.ncs', 'CSC10.ncs', 'CSC15.ncs'};
tt4_csc = LoadCSC(cfg_tt4);
corr_tt4 = corrcoef(tt4_csc.data');
%%
cfg_tt3.fc = {'CSC7.ncs', 'CSC2.ncs', 'CSC1.ncs', 'CSC4.ncs'};
tt3_csc = LoadCSC(cfg_tt3);
corr_tt3 = corrcoef(tt3_csc.data');
%%
