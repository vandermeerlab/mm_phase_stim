cd('E:\Dropbox (Dartmouth College)\manish_data\M325\')
sess_list = dir;
sess_list(1:2) = [];
main_dir = cd;

good_sess_list = {'M320-2022-05-28', ...
    'M319-2022-06-28', ...
    'M321-2022-07-08', ...
    'M321-2022-07-09', ...
    'M321-2022-07-13', ...
    'M322-2022-07-19', ...
    'M322-2022-07-21', ...
    'M322-2022-07-22', ...
    'M325-2022-07-28', ...
    };
% good_sess_list = {'M16_2019_02_18',...
%     'M16_2019_02_19',...
%     'M16_2019_02_20',...
%     'M16_2019_02_22',...
%     'M16_2019_02_23',...
%     'M16_2019_02_25',...
%     'M16_2019_02_27',...
%     'M17_2019_02_15',...
%     'M17_2019_02_18',...
%     'M17_2019_02_19',...
%     'M17_2019_02_20',...
%     'M17_2019_02_21',...
%     'M18_2019_04_10',....
%     'M18_2019_04_11',....
%     'M18_2019_04_12',....
%     'M18_2019_04_13',....
%     'M18_2019_04_14',....
%     'M19_2019_04_12',....
%     'M19_2019_04_13',....
%     'M19_2019_04_14',....
%     };
% 
% This_list = {'M17_2019_02_15', 'M16_2019_02_17', 'M16_2019_02_23', 'M16_2019_02_25', 'M16_2019_02_27'};
for iSess = 1:length(sess_list) % made it to 30
%         if ismember(strrep(sess_list(iSess).name(1:14), '-', '_'), This_list)%~strcmp(sess_list(iSess).name(1:3), 'M16')  && ~strcmp(sess_list(iSess).name(1:3), 'M17') && 
    if ~strcmp(sess_list(iSess).name(1), 'M') ||   ~ismember(sess_list(iSess).name(1:15), good_sess_list)
        continue
    elseif ~strcmp(sess_list(iSess).name(end-3:end), '.zip')

        disp(sess_list(iSess).name)
        cd(sess_list(iSess).name)

        % remove odd hidden '._MXX..._ExpKeys.m' that causes problems
        these_files = dir;
        for iF = 1:length(these_files)
            if length(these_files(iF).name) <=2
                continue
            else
                if strcmp(these_files(iF).name(1:2), '._')
                    delete(these_files(iF).name)
                end
            end
        end

        % actually run the
%             statedep_all_phase_sandbox; close all
        statedep_latency_manish; close all
%             statedep_PETH_Phase
        close all
        cd(main_dir)

    end
end


%%
% main_dir = '/Volumes/Fenrir/State_dep/EC_State'
% good_list = {'M13-2018-12-16_dStr_3p6_light_cells_TT3_min', 'M13-2018-12-11_vStr_4p1_light_cells_TT7_min', ...
%     'M13-2018-12-17_dStr_3p7_light_cells_TT1_2_min', 'M14-2018-12-08_dStr_2p8_light_cells_tt1_min', 'M14-2018-12-09_dStr_2p6_light_cells_TT6_8_post_min', ...
%     'M14-2018-12-10_vStr_2p9_light_cells_TT1_TT2_min', 'M14-2018-12-15_dStr_3p7_light_cells_TT1_3_min', 'M14-2018-12-17_dStr_3p6_light_cells_TT2_3_min'};
% 
% for iSess =good_list
%         if ~strcmp(iSess{1}(1), 'M')
%             continue
%         elseif ~strcmp(iSess{1}(end-3:end), '.zip')
%             
%             disp(iSess{1})
%             cd(iSess{1})
%             
%             % remove odd hidden '._MXX..._ExpKeys.m' that causes problems
%             these_files = dir;
%             for iF = 1:length(these_files)
%                 if length(these_files(iF).name) <=2
%                     continue
%                 else
%                     if strcmp(these_files(iF).n ame(1:2), '._')
%                         delete(these_files(iF).name)
%                     end
%                 end
%             end
%             
%             % actually run the
%             statedep_all_phase_sandbox
%             statedep_latency_phase_sandbox
%             close all
%             cd(main_dir)
%             
%         end
% end
% 
% 
% %%
% % get the firing rate
% figure(1)
% for iPhase = 1:n_phases
%     [~, this_phase_idx] = find(all_lat_2(1,:)==iPhase);
% subplot(2,3,iPhase)
% evt_count = 0;
%         hold on
%         counts(1,this_phase_idx) = iPhase;
% % plot(linspace(0,0.001),-1*ones(1,length(linspace(0,0.001))))
% for iEvt = this_phase_idx
%     evt_count = evt_count+1;
%     binsize = 0.01;
%     this_event_spikes = restrict(this_S, laser_on.t{1}(iEvt), laser_on.t{1}(iEvt)+0.025);
%     % this_event_spikes = restrict(this_S, this_peak_time{iPhase}(1)+laser_on.t{1}(iEvt), this_peak_time{iPhase}(2)+laser_on.t{1}(iEvt));
%     %     tbin_edges = laser_on.t{1}(iEvt):binsize:laser_on.t{1}(iEvt)+1; % vector of time bin edges (for histogram)
%     
% %     tbin_edges = laser_on.t{1}(iEvt):binsize:laser_on.t{1}(iEvt)+0.06; % vector of time bin edges (for histogram)
% %     tbin_centers = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)
% %     
% %     spk_count = histc(this_event_spikes.t{1},tbin_edges); % get spike counts for each bin
% %     spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
%     %     this_tvec = tbin
%     %     spk_count = spk_count/binsize;  % correct for time.
%     if ~isempty(this_event_spikes.t{1}) && length(this_event_spikes.t{1}) >2
%         plot([this_event_spikes.t{1}-laser_on.t{1}(iEvt) this_event_spikes.t{1}-laser_on.t{1}(iEvt)],[evt_count-0.5 evt_count+0.5],'Color',[0 0 0], 'linewidth', 4);
%     elseif ~isempty(this_event_spikes.t{1}) && length(this_event_spikes.t{1}) ==2 % for odd behaviour that plots two values on a diagonal
%         for ii = 1:2
%             plot([this_event_spikes.t{1}(ii)-laser_on.t{1}(iEvt) this_event_spikes.t{1}(ii)-laser_on.t{1}(iEvt)],[evt_count-0.5 evt_count+0.5],'Color',[0 0 0], 'linewidth', 4);
%         end
%     elseif ~isempty(this_event_spikes.t{1}) && length(this_event_spikes.t{1}) ==1
%         plot([this_event_spikes.t{1}(1)-laser_on.t{1}(iEvt) this_event_spikes.t{1}(1)-laser_on.t{1}(iEvt)],[evt_count-0.5 evt_count+0.5],'Color',[0 0 0], 'linewidth', 4);
%     else
%         disp([num2str(iEvt) 'n spikes ' num2str(length(this_event_spikes.t{1}))])
%     end
%     
%     if ~isempty(this_event_spikes.t{1})
%                 disp([num2str(iEvt) 'n spikes ' num2str(length(this_event_spikes.t{1}))])
%         counts(2, iEvt) = length(this_event_spikes.t{1});
%     else
%         counts(2, iEvt) = NaN;
%     end
% %     xlim([0 0.005])
% %     set(gca, 'xtick', [0 0.005 0.01])
% %     pause(0.1)
% %     if ~isnan(all_lat(iPhase,iEvt))
% %         all_spk_counts(iEvt,:) = spk_count;
% %     else
% %         all_spk_counts(iEvt,:) = nan(1, length(tbin_centers)) ;
% %     end
% end
% x_tick_val = get(gca, 'xtick');
% set(gca, 'xtick', x_tick_val*1000);
% end
