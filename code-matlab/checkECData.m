top_dir = 'D:\Dropbox (Dartmouth College)\EC_State_inProcess';

good_cells =  {'M16_2019_02_15_4p2_TT6_SS_02_Good',...
'M16_2019_02_16_3p3_TT5_SS_01_Great',...
'M16_2019_02_17_4_TT7_SS_01_Good',...
'M16_2019_02_18_4_TT4_SS_01_Good',...
'M16_2019_02_19_3p9_TT4_SS_01_Good',...
'M16_2019_02_20_4p3_TT5_SS_01_Good',...
'M16_2019_02_22_4p4_TT8_SS_01_Good',...
'M16_2019_02_23_2p5_TT4_SS_01_Good',...
'M16_2019_02_25_3p1_TT1_SS_01_Good',...
'M16_2019_02_26_3p3_TT5_SS_01_Good',...
'M16_2019_02_27_3p9_TT2_SS_01_Good',...
'M17_2019_02_15_3_TT3_SS_02_OK',...
'M17_2019_02_16_3p2_TT5_SS_01_Good',...
'M17_2019_02_16_3p6_TT4_SS_02_Good',...
'M17_2019_02_17_3_TT8_SS_04_OK',...
'M17_2019_02_18_3p7_TT5_SS_02_OK',...
'M17_2019_02_19_2p5_TT5_SS_01_OK',...
'M17_2019_02_20_3p7_TT5_SS_01_Good',...
'M17_2019_02_21_4p2_TT7_SS_01_Good',...
'M17_2019_02_24_3p6_TT1_SS_01_Good',...
'M17_2019_02_25_3p1_TT5_SS_02_OK',...
'M17_2019_02_25_3p9_TT1_SS_01_Good',...
'M18_2019_04_10_3p8_TT4_SS_01_Good',...
'M18_2019_04_11_4p2_TT7_SS_01_Good',...
'M18_2019_04_12_3p4_TT4_SS_02_Good',...
'M18_2019_04_12_3p8_TT7_SS_01_Good',...
'M18_2019_04_13_3p8_TT6_SS_02_OK',...
'M18_2019_04_14_4_TT3_SS_01_Good',...
'M18_2019_04_15_3p3_TT8_SS_08_Good',...
'M18_2019_04_15_4_TT7_SS_01_Good',...
'M19_2019_04_12_4p2_TT7_SS_01_OK',...
'M19_2019_04_13_4p7_TT2_SS_01_OK',...
'M19_2019_04_13_4p2_TT3_SS_01_OK',...
'M19_2019_04_14_3p3_TT8_SS_01_OK',...
'M19_2019_04_14_4p2_TT5_SS_01_Good',...
'M20_2019_06_07_3p8_TT8_SS_01_Good',...
'M20_2019_06_07_4p6_TT8_SS_01_Good',...
'M20_2019_06_08_4p2_TT7_SS_02_Good',...
'M20_2019_06_09_4p7_TT7_SS_01_Good',...
'M20_2019_06_10_4p2_TT6_SS_01_OK'};

% Add after find the ExpKeys
% 'M19_2019_04_15_4_TT6_SS_01_OK',...

ntt_exists = false(size(good_cells));
folder_paths = {};

% For each of the cells, generate Folder Path and check if NTT files exist 
for iC = 1:length(good_cells)
    tokens = split(good_cells{iC},'_');
    this_mouse = tokens{1};
    this_session = join(tokens(1:4),'-');
    this_session = this_session{1};
    search_string = join({top_dir, this_mouse, this_session}, '\');
    search_string = strcat(search_string{1},'*');
    this_tetrode = tokens{6};
    this_tetrode(length(this_tetrode)+1) = this_tetrode(end);
    this_tetrode(end-1) = '0';
    this_ntt_file = strcat(this_tetrode, '.ntt');
    this_dir = dir(search_string);
    for iD = 1:length(this_dir)
        this_filename = join({this_dir(iD).folder, this_dir(iD).name, [this_session, '_', this_ntt_file]}, '\');
        if isfile(this_filename)
            ntt_exists(iC) = true;
            this_folder =  join({this_dir(iD).folder, this_dir(iD).name}, '\');
            folder_paths{length(folder_paths)+1} = this_folder{1};
            clear this_dir this_mouse this_session this_ntt_file this_tetrode this_filename search_string
            break
        else
            dummy = 1;
        end
    end
end
intact_good_cells = good_cells(ntt_exists);

%% Separate session into variable ISI and fixed ISI
is_variable_ISI = false(size(folder_paths));
pushdir(top_dir);
for iF = 1:length(folder_paths)
   pushdir(folder_paths{iF});
   LoadExpKeys;
   if ExpKeys.stim_mode == 'variable ISI'
       is_variable_ISI(iF) = true;
   end
   popdir;
end
% All sessions have variable ISI
%% 
pushdir(top_dir);
artifact_main = cell(size(folder_paths));
artifact_pre = cell(size(folder_paths));
artifact_post = cell(size(folder_paths));
count_main = zeros(size(folder_paths));
count_pre = zeros(size(folder_paths));
count_post = zeros(size(folder_paths));
for iF = 1:length(folder_paths)
   pushdir(folder_paths{iF});
   LoadExpKeys;
   evs = LoadEvents([]);
   main_on = evs.t{find(cellfun(@(x) strcmp(x, ExpKeys.var_trig), evs.label))};
   short_on = evs.t{find(cellfun(@(x) strcmp(x, ExpKeys.short_laser_on), evs.label))};
   pre_on = short_on(short_on <= ExpKeys.Hundred_Stims(2));
   post_on = short_on(short_on > ExpKeys.Hundred_Stims(2));
   count_main(iF) = length(main_on);
   count_pre(iF) = length(pre_on);
   count_post(iF) = length(post_on);
   S = LoadSpikes([]);
   iv_pre = iv(pre_on+0.001, pre_on+0.002);
   iv_post = iv(post_on+0.0011, post_on+0.0021);
   iv_main = iv(main_on+0.0011, main_on+0.0021);
   artifact_pre{iF} = restrict(S, iv_pre);
   artifact_post{iF} = restrict(S, iv_post);
   artifact_main{iF} = restrict(S, iv_main);
   dummy = 1;
   popdir;
end
bad_pre = count_pre ~= 100;
bad_post = count_post ~= 100;
bad_main = count_main ~= 1500;
%%
%Look at 'bad' sessions later
all_good = ~(bad_pre | bad_post | bad_main);
good_main = artifact_main(all_good);
good_pre = artifact_pre(all_good);
good_post = artifact_post(all_good);
spk_main = zeros(size(good_main));
spk_pre = zeros(size(good_pre));
spk_post = zeros(size(good_post));
for iS = 1:length(good_main)
    spk_main(iS) = length(good_main{iS}.t{1});
    spk_pre(iS) = length(good_pre{iS}.t{1});
    spk_post(iS) = length(good_post{iS}.t{1});
end


