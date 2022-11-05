% /Volumes/Fenrir/State_dep/EC_State/M17-2019-02-21_vStr_4p2_light_cells_TT7_min
%  CSC22

data_temp = restrict(this_csc, start_stop.t{1}(main_rec_idx)+10, start_stop.t{1}(main_rec_idx)+20);

plot(data_temp.tvec, data_temp.data, 'k', 'linewidth', 2);

xlim([start_stop.t{1}(main_rec_idx)+11, start_stop.t{1}(main_rec_idx)+12])

hold on

    f_list = {[1 2], [3 5], [7 10], [30 40],[40 60]};
    f_list_label = {'3 - 5', '7 - 10', '15 - 25', '40 - 60'};
    
for iF = 2:length(f_list) % loop across freqs
cfg_filt = []; cfg_filt.type = 'fdesign'; cfg_filt.f  = f_list{iF};
csc_f = FilterLFP(cfg_filt, data_temp);

plot(csc_f.tvec, csc_f.data-(iF*0.0001), 'linewidth', 2)
xlim([start_stop.t{1}(main_rec_idx)+11, start_stop.t{1}(main_rec_idx)+12])
end