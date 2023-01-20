%% Assumes that good LFPs have been picked out
top_dir = 'E:\Dropbox (Dartmouth College)\EC_State_inProcess\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020'};
% top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
% mice = {'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        doStuff
    end
end


function doStuff
    LoadExpKeys;
    evs = LoadEvents([]);
    cfg_spk = [];
    cfg_spk.fc = ExpKeys.goodCell;
    if isempty(cfg_spk.fc)
        return
    end
    
%     % For debugging
%     if ~(ExpKeys.hasWheelData) | contains(ExpKeys.light_source, 'LASER')
%         return
%     end

    cfg = []; cfg.fc = ExpKeys.goodLFP;
    if contains(cfg.fc, '-')
        temp = split(cfg.fc,'-');
        cfg.fc = {cat(2,temp{1},'.ncs')};
        csc = LoadCSC(cfg);
        cfg_temp.fc = {cat(2,temp{2},'.ncs')};
        ref = LoadCSC(cfg_temp);
        csc.data = csc.data - ref.data;
        clear temp ref;
    else
        csc = LoadCSC(cfg);
    end
    
    fig = figure('WindowState', 'maximized');
    ax1 = subplot(3, 1, [1 2]);
    
    if csc.cfg.hdr{1}.SamplingFrequency > 30000
        cfg = []; cfg.decimateFactor = 12;
        csc = decimate_tsd(cfg, csc); 
    end
    
    % Save the original time_base for speed_data
    og_tvec = csc.tvec;
    csc = restrict(csc, iv(ExpKeys.stim_times));

    Fs = 1/median(diff(csc.tvec));
    % 1 sec window works for most sessions except some
    weird_spec_8 = {'2021-12-27'};
    weird_spec_4 = {'2022-05-28'};
    weird_spec_2 = {'2022-08-03'};
    if contains(pwd, weird_spec_8)
        wsize = pow2(floor(log2(8*Fs)));
    elseif contains(pwd, weird_spec_4)
        wsize = pow2(floor(log2(4*Fs)));
    elseif contains(pwd, weird_spec_2)
        wsize = pow2(floor(log2(2*Fs)));
    else
        wsize = pow2(floor(log2(Fs)));
    end
    % Plot spectrogram
    [S F T P] = spectrogram(csc.data, hanning(wsize), wsize/2, 1:120, Fs);
%     imagesc(T,F,10*log10(P));
    % Saving spectrogram stuff that was calculated for the trial stimperiod
    spec_data = [];
    spec_data.wsize  = wsize;
    spec_data.Fs = Fs;
    spec_data.S = S;
    spec_data.T = T;
    spec_data.P = P;
    spec_data.F = F;
    spec_data.firstTimeStamp = csc.tvec(1);
    save('spectrogram_data', 'spec_data');
    title('Spectrogram');
    ylabel('Freq (Hz)');
    xlabel('Time (sec)');
    colorbar
    ax1.FontSize = 14;
    if isfield(ExpKeys, 'hasWheelData') & ExpKeys.hasWheelData
        % If running wheel data present, process below
        if ~strcmp(ExpKeys.experimenter, 'EC')
            if contains(ExpKeys.light_source, 'LASER')  % Surgery config
                str_AB = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0030).'; % both up
                str_A = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0020).'; % only chA up
                str_B = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0010).'; % only chB up
                str_d = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0000).'; % both down
                cfg_w = [];
                cfg_w.chA='WE1.ncs';
                cfg_w.chB='WE2.ncs';
                updownTSD = getQEupdown(cfg_w);
                state_tsd = ConvertQEUpDownToState(updownTSD);
                [angle_tsd, wheel_tsd, bad_jumps] = ConvertQEStatesToAngle([], state_tsd);
                [d, speed, cfg_w] = ConvertWheeltoSpeed(cfg_w, wheel_tsd);
            else    % RR2 config                     
                str_AB = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0003).'; % both up
                str_A = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0002).'; % only chA up
                str_B = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0001).'; % only chB up
                str_d = 'TTL Input on AcqSystem1_0 board 0 port 2 value (0x0000).'; % both down
                both_up = evs.t{strcmp(evs.label, str_AB)};
                chA_up = evs.t{strcmp(evs.label, str_A)};
                chB_up = evs.t{strcmp(evs.label, str_B)};
                both_down = evs.t{strcmp(evs.label, str_d)};
               
                all_events = [both_up; chA_up; chB_up; both_down]';
                sparse_states = [ones(size(both_up)); 2*ones(size(chA_up)); ...
                    3*ones(size(chB_up)); 4*ones(size(both_down))]';
                [all_events, sort_idx] = sort(all_events);
                sparse_states = sparse_states(sort_idx);
                % Use LFP timebase
                state_tsd.tvec = og_tvec;
                state_tsd.data = zeros(size(state_tsd.tvec));
                cidx = nearest_idx3(all_events, state_tsd.tvec);
                state_tsd.data(1:cidx(1)) = sparse_states(1);
                state_tsd.data(cidx(end):end) = sparse_states(end);
                for i = 1:length(cidx)-1
                    state_tsd.data(cidx(i):cidx(i+1)-1) = sparse_states(i);
                end 
                assert(isempty(find(state_tsd.data == 0)));
                cfg_w = [];
                [angle_tsd, wheel_tsd, bad_jumps] = ConvertQEStatesToAngle(cfg_w, state_tsd);
                [d, speed, cfg_w] = ConvertWheeltoSpeed(cfg_w, wheel_tsd);
            end
                % Assume you are in the correct folder
                save('speed_data','speed');
                keep = (speed.tvec >= ExpKeys.stim_times(1) & speed.tvec <= ExpKeys.stim_times(2));
                sub_speed = [];
                sub_speed.tvec = speed.tvec(keep);
                sub_speed.data = speed.data(keep);
                sub_speed.tvec = sub_speed.tvec - sub_speed.tvec(1); % converting to 0 time base
                % For now no anti-aliasing filter, just subsampling the timepoints from the spectrogram timebase
                keep = nearest_idx3(T, sub_speed.tvec);
                sub_speed.tvec = sub_speed.tvec(keep);
                sub_speed.data = sub_speed.data(keep);
                ax2 = subplot(3, 1, 3);
                plot(sub_speed.tvec, sub_speed.data);
                ylabel('Running Speed (cm/sec)');
                xlabel('Time (sec)');
                ax2.FontSize = 14;
                linkaxes([ax1, ax2], 'x');
        else
            % Depends on what Eric says about the running data
%             run_file = FindFiles('*run.csv');
%             if ~isempty(run_file)
%                 [Sys_time, Real_time, Encoder] = Load_Wheel_EC(run_file{1});
%                 figure;
% %                 subplot(3,ceil(length(S.label)/3)+1, length(S.label)+1) ;
%                 tvec = (Real_time-Real_time(1))/60;
%                 plot(tvec, [0 ;diff(Encoder)]);
%                 xlabel('Time (min)');
%                 ylabel('Encoder position')
%                 [~,breaks] = findpeaks(diff(tvec), 'MinPeakHeight', 0.05);
%                 hold on
%                 plot(tvec([breaks]), zeros(1,length(breaks)), '*r', 'markersize', 4)
%                 xlim([tvec(1) tvec(end)])
%             end
        end
    end
    fn_prefix = split(pwd, '\');
    fn_prefix = fn_prefix{end};
    WriteFig(fig, strcat('E:\Dropbox (Dartmouth College)\EC_State_inProcess\', fn_prefix, '-specNspeed'),1);
    close;
end