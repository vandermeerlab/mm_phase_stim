
function statedep_summary
%% collect the outputs from statedep_all_phase_sandbox and statedep_latency_phase_sandbox
%
%
%
STATE_init
%% defaults


global PARAMS
if isunix
    addpath(genpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'));
    addpath('/Users/jericcarmichael/Documents/GitHub/EC_state/Basic_functions');
    
    all_fig_dir = '/Volumes/Fenrir/State_dep/all_checks/';
    all_lat_dir = '/Volumes/Fenrir/State_dep/all_lat/';
    all_ccf_dir = '/Volumes/Fenrir/State_dep/all_ccf/';
    summary_dir = '/Volumes/Fenrir/State_dep/Summary/';
else
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\EC_State'))
    
    all_fig_dir = 'G:\State_data\all_checks\';
    all_lat_dir = 'G:\State_data\all_lat\';
    
end

mkdir(all_lat_dir); mkdir(all_fig_dir);



font_size = 18;
%% generate a latency and count summary
cd(all_lat_dir);
sess_list = dir(all_lat_dir);
sess_list = sess_list(3:end);

for iSess = 1:length(sess_list)
    load(sess_list(iSess).name);
    
    these_cells = fieldnames(out);
    for iC = 1:length(these_cells)
        
        this_cell = these_cells{iC};
        if ismember([sess_list(iSess).name(1:14) '_' this_cell], PARAMS.Good_cells) % check if this is a 'good' cell from above approvedd list
            
            freq_list = fieldnames(out.(this_cell));
            
            
            all_cells.([sess_list(iSess).name(1:14) '_' this_cell]) = out.(this_cell);
            
            % subject hdr
            %             all_cells.(this_cell).hdr.name = sess_list(iSess).name(1:14);
            %             all_cells.(this_cell).hdr.date = sess_list(iSess).name(5:14);
            %             all_cells.(this_cell).hdr.subject = sess_list(iSess).name(1:3);
        else
            disp(['Skippking ' sess_list(iSess).name(1:14) '_' this_cell '.  Not a ''good'' cell'])
        end
        
    end
    
end

%%
nShuf = 10000;
cell_list = fieldnames(all_cells);
all_lat =[];
all_count = [];
all_resp_ratio = [];
for iC = 1:length(cell_list)
    this_cell = cell_list{iC};
    f_list = fieldnames(all_cells.(this_cell));
    cell_depths(iC) = all_cells.(this_cell).ExpKeys.tetrodeDepths;
    for iF = 1:length(f_list)
        if strcmp(f_list{iF}, 'ExpKeys') || strcmp(f_list{iF}, 'hdr')
            continue
        else
            x_phase_resp(iC, iF) = nanmean(all_cells.(this_cell).(f_list{iF}).resp(2,:));
            
            
            %% get the shuffle for the response max-min
            for iShuf = nShuf:-1:1
                this_shuf = [];
                mix = randperm(length(all_cells.(this_cell).(f_list{iF}).latency(1,:)));
                this_shuf(1,:) = all_cells.(this_cell).(f_list{iF}).latency(1,mix);
                
                for iPhase = unique(all_cells.(this_cell).(f_list{iF}).latency(1,:))
                    this_phase_idx = find(this_shuf(1,:) == iPhase);
                    this_phase_mean(iPhase) = nanmean(all_cells.(this_cell).(f_list{iF}).resp(2,this_phase_idx));
                end
                all_shuf_ratio(iShuf) = max(this_phase_mean)/min(this_phase_mean);
                % check for inf from shuffle
                all_shuf_ratio(isinf(all_shuf_ratio)) = NaN;
                
                
            end
            % get the phase response averages
            for iPhase = unique(all_cells.(this_cell).(f_list{iF}).latency(1,:))
                this_phase_idx = find(all_cells.(this_cell).(f_list{iF}).latency(1,:) == iPhase);
                this_resp(iPhase) = nanmean(all_cells.(this_cell).(f_list{iF}).resp(2,this_phase_idx));
            end
            %
            all_resp_phase(iC,iF) = max(this_resp)/min(this_resp);
            all_resp_shuffle(iC, iF) = nanmean(all_shuf_ratio);
            all_resp_ratio(iC, iF) = ((max(this_resp)/min(this_resp)) - nanmean(all_shuf_ratio))/nanstd(all_shuf_ratio);
            [~, all_resp_p_phase_max(iC, iF)] = max(this_resp);
            [~, all_resp_p_phase_min(iC, iF)] = min(this_resp);
            
            if isnan(all_resp_ratio(iC, iF))
                disp([num2str(iC), 'f' num2str(iF)])
            end
            
            
            
            
            %% get the
            
            for iPhase = unique(all_cells.(this_cell).(f_list{iF}).latency(1,:)) % get the phase segment numbers (should be 1:5 if phase split by 5
                % labels for x axes
                phase_labels{iPhase} = ['p' num2str(iPhase)];
                
                % get the latencies
                this_phase_idx = find(all_cells.(this_cell).(f_list{iF}).latency(1,:) == iPhase);
                all_lat(iC, iPhase,iF)  = nanmean(all_cells.(this_cell).(f_list{iF}).latency(2,this_phase_idx)); % in ms get the mean value for this cell in this freq at this phase bin
                % same for spike count
                all_count(iC, iPhase, iF) = nanmean(all_cells.(this_cell).(f_list{iF}).count(2,this_phase_idx)); % in ms get the mean value for this cell in this freq at this phase bin
                
                all_resp(iC,iPhase,iF) = nanmean(all_cells.(this_cell).(f_list{iF}).resp(2,this_phase_idx));
                
                % get all the reponses
                
                % maybe shuffles?
                all_lat_shuf_mean(iC, iPhase, iF) = nanmean(all_cells.(this_cell).(f_list{iF}).latency_shuf(1,this_phase_idx)); % in ms get the mean value for this cell in this freq at this phase bin
                all_lat_shuf_std(iC, iPhase, iF) = nanstd(all_cells.(this_cell).(f_list{iF}).latency_shuf(1,this_phase_idx)); % in ms get the mean value for this cell in this freq at this phase bin
                
                all_resp_shuf_mean(iC,iPhase, iF) = nanmean(all_cells.(this_cell).(f_list{iF}).resp_shuf(1,this_phase_idx));
                all_resp_shuf_std(iC,iPhase, iF) = nanstd(all_cells.(this_cell).(f_list{iF}).resp_shuf(1,this_phase_idx));
                
                
                if all_lat(iC, iPhase,iF) > all_lat_shuf_mean(iC,iPhase,iF)+ 1.96*(all_lat_shuf_std(iC,iPhase,iF)) || all_lat(iC, iPhase,iF) < all_lat_shuf_mean(iC,iPhase,iF)- 1.96*(all_lat_shuf_std(iC,iPhase,iF))
                    all_lat_v_shuf(iC,iPhase,iF) = all_lat(iC, iPhase, iF);
                else
                    all_lat_v_shuf(iC, iPhase, iF) = NaN;
                end
                
                %exceeds 2sd of shuffle for response
                %                  if all_resp(iC, iPhase,iF) > all_resp_shuf_mean(iC,iPhase,iF)+ 1.96*(all_resp_shuf_std(iC,iPhase,iF)) || all_resp(iC, iPhase,iF) < all_resp_shuf_mean(iC,iPhase,iF)- 1.96*(all_resp_shuf_std(iC,iPhase,iF))
                all_resp_v_shuf(iC,iPhase,iF) = (all_resp(iC,iPhase,iF) - all_resp_shuf_mean(iC,iPhase,iF))/all_resp_shuf_std(iC,iPhase,iF);
                %                 else
                %                     all_resp_v_shuf(iC, iPhase, iF) = NaN;
                %                 end
                
            end
            z_lat(iC, :,iF) = zscore(all_lat(iC, :,iF));
            z_resp(iC, :,iF) = zscore(all_resp(iC, :,iF));
            all_resp_v_all(iC,:,iF) = all_resp(iC,:,iF)./x_phase_resp(iC,iF);
            
            sess_id{iC} = this_cell;
            freq_vals(iC, iF) = iF;
        end
    end
    
end
%% sort base on depth
[depth_sort, depth_ord] = sort(cell_depths, 'descend');
x_phase_resp_sort = x_phase_resp(depth_ord,:);
all_resp_v_shuf_sort = all_resp_v_shuf(depth_ord,:,:);
all_lat_sort = all_lat(depth_ord,:,:);
all_resp_sort = all_resp(depth_ord,:,:);
all_resp_v_all_sort = all_resp_v_all(depth_ord,:,:);
all_lat_v_shuf_sort = all_lat_v_shuf(depth_ord,:,:);
all_count_sort = all_count(depth_ord,:,:);
all_resp_ratio_sort = all_resp_ratio(depth_ord, :,:);
all_resp_p_phase_max_sort = all_resp_p_phase_max(depth_ord, :,:);
all_resp_p_phase_min_sort = all_resp_p_phase_min(depth_ord, :,:);

for ii =1:length(depth_sort)
    
    labels{ii} = strcat( num2str(floor(x_phase_resp_sort(ii)*100)), '%'); %num2str(depth_sort(ii)), '_m_m ',
    all_sess_id_sort{ii} = sess_id{depth_ord(ii)};
    
end

% get the freq list for legend
for iF = 1:length(f_list)
    if strcmp(f_list{iF}, 'ExpKeys') || strcmp(f_list{iF}, 'hdr')
        continue
    else
        leg_freq{iF} = [ strrep(f_list{iF}(3:end), '_', '-') 'Hz'];
    end
    
    
end

%% remove cells with a response less than 20%
cut_off = 0.2;
for iC = length(x_phase_resp_sort):-1:1
    if x_phase_resp_sort(iC) < cut_off
        all_resp_v_shuf_sort(iC,:,:) = [];
        all_lat_sort(iC,:,:) = [];
        all_resp_sort(iC,:,:) = [];
        all_resp_v_all_sort(iC,:,:) = [];
        all_lat_v_shuf_sort(iC,:,:) = [];
        all_count_sort(iC,:,:) = [];
        all_resp_ratio_sort(iC,:,:) = [];
        all_resp_p_phase_max_sort(iC,:,:) = [];
        all_resp_p_phase_min_sort(iC,:,:) = [];
        x_phase_resp_sort(iC,:,:) = [];
        
        
        freq_vals(iC,:) = [];
        labels{iC} = [];
        %         all_sess_id_sort{iC} = [];
        depth_sort(iC) = [];
    end
end
labels = labels(~cellfun('isempty',labels));
% all_sess_id_sort = all_sess_id_sort(~cellfun('isempty',all_sess_id_sort));
%%
% relative to shuffle
% zscore
figure(8)
for iF = 1:length(f_list)
    if strcmp(f_list{iF}, 'hdr') || strcmp(f_list{iF}, 'ExpKeys')
        continue
    else
        subplot(2,ceil(length(f_list)-1)/2, iF);
        imagesc(all_resp_v_shuf_sort(:,:,iF));
        axis xy
        xlabel('phase');
        ylabel('cell');
        set(gca, 'xticklabel', []);
        set(gca,'ytick', [1:length(all_resp_v_shuf_sort(:,:,iF))], 'yticklabel',labels, 'fontsize',14); % [num2str(cell_depths) ' - ' num2str(floor(x_phase_resp(:,iF)*100))]
        text(floor(length(phase_labels)/2),length(all_resp_v_shuf_sort)+1, f_list{iF}, 'fontsize', font_size)
        %         if strcmp(f_list{iF}, 'f_60_80')
        %         colorbar
        caxis([-2.5 2.5])
        %         end
    end
    sig_cells = double(all_resp_v_shuf_sort(:,:,iF)>1.96);
    temp_all_resp_v_shuf = all_resp_v_shuf_sort(:,:,iF);
    temp_all_resp_v_shuf(sig_cells==0)=NaN;
    %     add_num_imagesc(gca, temp_all_resp_v_shuf, 2, 12)
end
SetFigure([], gcf)
% set(gca, 'fontsize', 8);
% for iSub = 1:6
%     subplot(2,3,iSub)
% set(gca, 'yticklabel',labels, 'fontsize',12)
% end
set(gcf, 'position', [282   50  1100  800])
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.35, 0.98,['Zscore v Shuf per stim across cells'], 'fontsize', font_size)
saveas(gcf, [summary_dir 'Summary_resp_v_shuf.png']);
saveas_eps('Summary_resp_v_shuf', summary_dir)

% make a plot just for the color bar
figure(888)
subplot(2,3,6)
imagesc([-2.5 2.5])
xlabel('phase');
ylabel('cell');
set(gca, 'xticklabel', []);
set(gca,'ytick', [1:length(all_resp_v_shuf_sort(:,:,3))], 'yticklabel',labels, 'fontsize',14 ); % [num2str(cell_depths) ' - ' num2str(floor(x_phase_resp(:,iF)*100))]
text(floor(length(phase_labels)/2),length(all_resp_v_shuf_sort)+1, f_list{3}, 'fontsize', font_size)
colorbar
caxis([-2.5 2.5])
SetFigure([], gcf)
set(gcf, 'position', [282   50  1100  800])
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
saveas_eps('Summary_resp_v_shuf_Cbar', summary_dir)
%% make a response plot



% get all responses p(Spike|stim)

% x_phase_resp should be the same across f_list since it is the number of
% responses and doesn't care about phase or freq.
f_cor = linspecer(size(all_resp_ratio,2));
figure(9)
freq_colors = [];
for iC =1:size(all_resp_ratio_sort,1)
    for iF = 1:size(all_resp_ratio_sort,2)
        % % get max phase
        % % norm to lowest
        %
        % [max_val, max_idx] = max(all_resp_v_shuf(iC,:,iF))%./min(all_resp_v_shuf(iC,:,iF)));
        %
        % % get min phase
        % [min_val, min_idx] = min(all_resp_v_shuf(iC,:,iF));
        %
        %
        % % get ratio
        % resp_ratio(iC,iF) = abs(max_val-min_val);
        freq_colors{iC,iF} = f_cor(iF,:);
        freq(iC, iF) = iF;
        cell_id(iC,iF) = iC;
    end
end
resp_ratio_1d = reshape(all_resp_ratio_sort, 1, numel(all_resp_ratio_sort));
resp_all_1d = reshape(x_phase_resp_sort, 1, numel(x_phase_resp_sort));
freq_colors_1d = reshape(freq_colors,1, numel(freq_colors));
freq_id_1d = reshape(freq_vals,1, numel(freq_vals));
depth_all_1d = reshape(repmat(depth_sort',1,size(x_phase_resp_sort,2)), 1, numel(x_phase_resp_sort));
cell_id_1d = reshape(cell_id,1, numel(cell_id));

% generate figure
for iC = 1:length(resp_all_1d)
    hold on
    if depth_all_1d(iC) >3.5
        S(iC) = scatter(resp_all_1d(iC), resp_ratio_1d(iC),100,freq_colors_1d{iC}, 'filled');
    else
        scatter(resp_all_1d(iC), resp_ratio_1d(iC),100,freq_colors_1d{iC}, 'filled', 'd')
    end
end
xlabel('p(spike|stim)')
ylabel('zscore response (max phase/min phase)')
%ylabel('ratio of max phase / min phase response')
h = zeros(size(all_resp_ratio,2), 1);
for iF = 1:size(all_resp_ratio,2)
    h(iF) = plot(NaN,NaN,'o', 'color', f_cor(iF, :), 'MarkerFaceColor', f_cor(iF,:));
end
lgd = legend(h, leg_freq, 'location', 'NorthEastOutside');
legend('boxoff')
lgd.FontSize = 18;

xlim([cut_off .7])
ylim([-2 10])
line([0 0.7], [1.96, 1.96], 'color', [0.3 0.3 0.3])
line([0 0.7], [3.1, 3.1], 'color', [0.3, 0.3, 0.3])

SetFigure([], gcf)
saveas(gcf, [summary_dir 'Summary_resp_ratio.png']);
saveas_eps('Summary_resp_ratio', summary_dir)

% extra legend figure
figure(111)
hold on
h(1) = plot(NaN,NaN,'o', 'color', 'k', 'MarkerFaceColor', 'k');
h(2) = plot(NaN,NaN,'d', 'color', 'k', 'MarkerFaceColor', 'k');
axis off
lgd = legend({'vStr', 'dStr'});
legend('boxoff')
lgd.FontSize = 18;
saveas_eps('Summary_Str_legend', summary_dir)
% %% make a polar of the prefered phase
%
%  ver = version;
%     ver = str2double(ver(end-5:end-2));
% n_phases = 5;
%         fprintf('\nPhase split into %1d bins\n', n_phases)
%         [~, edges, bin] = histcounts(-pi:pi, n_phases, 'BinLimits', [-pi, pi]);
%         phase_center = (edges(1:end-1) + edges(2:end))/2;
%         for iPhase = 1:max(all_resp_p_phase_max_sort)
%         all_resp_p_phase_max_sort(all_resp_p_phase_max_sort==iPhase) = phase_center(iPhase);
%         end
%
%         figure(123)
%         for iFreq = 1:size(all_resp_p_phase_max_sort,2)
%
%             phase_mean(iFreq) = circ_mean(all_resp_p_phase_max_sort(:,iFreq));
%             vector_r(iFreq) = circ_r(all_resp_p_phase_max_sort(:,iFreq));
%
%             if ver >= 2017
%                 h = polarplot( [0 phase_mean(iFreq)],[0 vector_r(iFreq)]);
%                 rlim([0 1])
%                 set(gca, 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top','RTick' ,[0.5 1], 'FontSize', 22, 'FontName', 'Helvetica')
%                 cfg_polar.axes_labels = 0; % turns off any scaling of the x and y label which don't work in polar plots
%             else
%                 h_temp = polar([0 0], [0 1]);
%                 hold on
%                 h = polar([0 phase_mean(iFreq)], [0 vector_r(iFreq)]);
%                 set(h_temp, 'visible', 'off') ;
%                 cfg_polar.axes_labels = 1; % turns off any scaling of the x and y label which don't work in polar plots
%             end
%
%
% %         [t,r ] = rose(all_resp_p_phase_max_sort(:,iFreq));
% %         subplot(2,3,iFreq)
% %         polar(t,r)
% %         polarhistogram(all_resp_p_phase_max_sort(:,iFreq),100)
% %         title(f_list{iFreq})
%         end
%% get p-values for each cell and use to build FDR

all_p = [];

for iC = 1:length(resp_ratio_1d)
    p_val = [];
    p_val = normcdf(resp_ratio_1d(iC),0, 1, 'upper');
all_p = [all_p, p_val];


end 

[fdr_out, q] = mafdr(all_p,'showplot', true)

Keep_idx = resp_all_1d >0.05;

keep_resp_ratio_1d = resp_ratio_1d(Keep_idx);
keep_resp_all_1d = resp_all_1d(Keep_idx);
keep_depth_1d = depth_all_1d(Keep_idx);
keep_freq_1d = freq_id_1d(Keep_idx);
keep_cell_1d = cell_id_1d(Keep_idx);

dStr_count = 0; vStr_count = 0;
d_freq_dist =NaN(1,length(keep_resp_all_1d));
v_freq_dist =NaN(1,length(keep_resp_all_1d));


fprintf('\nNumber of responsive cells: %.2d/%.2d = %.2d %%\n',length(unique(keep_cell_1d(fdr_out <0.05))),length(unique(keep_cell_1d)),round((length(unique(keep_cell_1d(fdr_out <0.05)))/length(unique(keep_cell_1d)))*100) )
for iC = 1:length(keep_resp_all_1d)
            this_f_str = strsplit(freq_list{keep_freq_1d(iC)},'_');
    if fdr_out(iC) <0.05 && depth_all_1d(iC) >3.5
        fprintf(['vStr Responsive cell ' all_sess_id_sort{keep_cell_1d(iC)} ' at ' num2str(keep_resp_ratio_1d(iC),2) ' SD, ' freq_list{keep_freq_1d(iC)} 'Hz. Overall Resp:' num2str((keep_resp_all_1d(iC)*100),2) '\n'])
        vStr_count = vStr_count+1;
        v_freq_dist(iC) = str2double(this_f_str{2});
            
    elseif fdr_out(iC) <0.05 && depth_all_1d(iC) <=3.5
        fprintf(['dStr Responsive cell ' all_sess_id_sort{keep_cell_1d(iC)} ' at ' num2str(keep_resp_ratio_1d(iC),2) ' SD, ' freq_list{keep_freq_1d(iC)} 'Hz. Overall Resp:' num2str((keep_resp_all_1d(iC)*100),2) '\n'])
        dStr_count = dStr_count+1;
        d_freq_dist(iC) = str2double(this_f_str{2});
    elseif depth_all_1d(iC) <=3.5
         d_freq_dist(iC) = 0;
         
    else
        v_freq_dist(iC) = 0;
% 
    end
    
end


fprintf('dStr: %.0d  vStr: %.0d\n',dStr_count, vStr_count )
fprintf('\n')
fprintf('Depths')
depth_sort;

%% make a table of responsive cells and which frequencies
Keep_idx = resp_all_1d >0.05;

keep_resp_ratio_1d = resp_ratio_1d(Keep_idx);
keep_resp_all_1d = resp_all_1d(Keep_idx);
keep_depth_1d = depth_all_1d(Keep_idx);
keep_freq_1d = freq_id_1d(Keep_idx);
keep_cell_1d = cell_id_1d(Keep_idx);

dStr_count = 0; vStr_count = 0;
d_freq_dist =NaN(1,length(keep_resp_all_1d));
v_freq_dist =NaN(1,length(keep_resp_all_1d));


fprintf('\nNumber of responsive cells: %.2d/%.2d = %.2d %%\n',length(unique(keep_cell_1d(keep_resp_ratio_1d >1.96))),length(unique(keep_cell_1d)),round((length(unique(keep_cell_1d(keep_resp_ratio_1d >1.96)))/length(unique(keep_cell_1d)))*100) )
for iC = 1:length(keep_resp_all_1d)
            this_f_str = strsplit(freq_list{keep_freq_1d(iC)},'_');
    if keep_resp_ratio_1d(iC) >1.96 && depth_all_1d(iC) >3.5
        fprintf(['vStr Responsive cell ' all_sess_id_sort{keep_cell_1d(iC)} ' at ' num2str(keep_resp_ratio_1d(iC),2) ' SD, ' freq_list{keep_freq_1d(iC)} 'Hz. Overall Resp:' num2str((keep_resp_all_1d(iC)*100),2) '\n'])
        vStr_count = vStr_count+1;
        v_freq_dist(iC) = str2double(this_f_str{2});
            
    elseif keep_resp_ratio_1d(iC) >1.96 && depth_all_1d(iC) <=3.5
        fprintf(['dStr Responsive cell ' all_sess_id_sort{keep_cell_1d(iC)} ' at ' num2str(keep_resp_ratio_1d(iC),2) ' SD, ' freq_list{keep_freq_1d(iC)} 'Hz. Overall Resp:' num2str((keep_resp_all_1d(iC)*100),2) '\n'])
        dStr_count = dStr_count+1;
        d_freq_dist(iC) = str2double(this_f_str{2});
    elseif depth_all_1d(iC) <=3.5
         d_freq_dist(iC) = 0;
         
    else
        v_freq_dist(iC) = 0;
% 
    end
    
end


fprintf('dStr: %.0d  vStr: %.0d\n',dStr_count, vStr_count )
fprintf('\n')
fprintf('Depths')
depth_sort;

%% make apie chart of the phase mod for vStr and dStr
% cell_mod = NaN(numel(all_resp_ratio_sort));
% for ii = 1:numel(all_resp_ratio_sort,1)
%     if all_resp_ratio_sort(ii) >1.96
%         
% %     for jj = size(all_resp_ratio_sort,2)
% %         if x_phase_resp_sort (ii, :) >1.96
%         all_resp_ratio_sort(ii, :)
%     end
% end
v_freq_dist_2d = reshape(v_freq_dist,size(all_resp_ratio_sort));
d_freq_dist_2d = reshape(d_freq_dist,size(all_resp_ratio_sort));

v_freq_dist_2d(any(isnan(v_freq_dist_2d), 2), :) = [];
d_freq_dist_2d(any(isnan(d_freq_dist_2d), 2), :) = [];

for ii =1:size(v_freq_dist_2d,1)
    if sum(v_freq_dist_2d(ii,:) ~=0) ==1
        v_cells(ii) = v_freq_dist_2d(ii,v_freq_dist_2d(ii,:) ~=0);
    elseif sum(v_freq_dist_2d(ii,:) ~=0) >1
        temp_vals = v_freq_dist_2d(ii,v_freq_dist_2d(ii,:) ~=0);
        v_cells(ii) = str2num(strrep(num2str(temp_vals), ' ', '.'));
    else
        v_cells(ii) = 0;
    end
end

for ii =1:size(d_freq_dist_2d,1)
    if sum(d_freq_dist_2d(ii,:) ~=0) ==1
        d_cells(ii) = d_freq_dist_2d(ii,d_freq_dist_2d(ii,:) ~=0);
    elseif sum(d_freq_dist_2d(ii,:) ~=0) >1
        temp_vals = d_freq_dist_2d(ii,d_freq_dist_2d(ii,:) ~=0);
        d_cells(ii) = str2num(strrep(num2str(temp_vals), '  ', '.'));
    else
        d_cells(ii) = 0;
    end
end




figure(9999)
subplot(2,1,1)
[a,b]=hist(d_cells,unique(d_cells));
p_a = a/sum(a);
for ii = 1:length(b)
pie_labels{ii} = strcat(num2str(b(ii)), ': ');
end
p = pie(p_a,[1, zeros(1,length(p_a)-1)]);

pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
combinedtxt = strcat(pie_labels,percentValues'); 
for ii = 1:length(pText)
    pText(ii).String = strrep(combinedtxt(ii), ':', ':  ');
    pText(ii).FontSize = 18;
    pText(ii).Position = [pText(ii).Position(1)*1.25,pText(ii).Position(2)*1, pText(ii).Position(3)] ;
end

title('dStr')

subplot(2,1,2)
[a,b]=hist(v_cells,unique(v_cells));
p_a = a/sum(a);
for ii = 1:length(b)
    
pie_labels{ii} = strcat(num2str(b(ii)), ': ');
end
p = pie(p_a,[1, zeros(1,length(p_a)-1)]);

pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
combinedtxt = strcat(pie_labels,percentValues'); 
for ii = 1:length(pText)
    pText(ii).String = strrep(combinedtxt(ii), ':', ':  ');
    pText(ii).FontSize = 18;
    pText(ii).Position = [pText(ii).Position(1)*1.25,pText(ii).Position(2)*1, pText(ii).Position(3)] ;
end
title('vStr')

SetFigure([], gcf)
saveas(gcf, [summary_dir 'Summary_Pie.png']);
saveas_eps('Summary_Pie', summary_dir)
%% format for GLM
cell_list = fieldnames(all_cells);
All_GLM_cells = zeros(length(cell_list), 6); % make an array for whether or not the phase or phase x amp is a predictor
% cycle cells



% initialize models
Phase_GLME.model.baseline.modelspec = 'Response ~ 1 + Pulse_num';
Phase_GLME.model.dphase.modelspec = 'Response ~ 1 + Pulse_num + Phase_delta';
Phase_GLME.model.tphase.modelspec = 'Response ~ 1 + Pulse_num + Phase_theta';
Phase_GLME.model.bphase.modelspec = 'Response ~ 1 + Pulse_num + Phase_beta';
Phase_GLME.model.lGphase.modelspec = 'Response ~ 1 + Pulse_num + Phase_lG';
Phase_GLME.model.mGphase.modelspec = 'Response ~ 1 + Pulse_num + Phase_mG';
Phase_GLME.model.hGphase.modelspec = 'Response ~ 1 + Pulse_num + Phase_hG';

% list the amp model
Phase_GLME.model.dphase_amp.modelspec = 'Response ~ 1 + Pulse_num + Phase_delta * Amp_delta';
Phase_GLME.model.tphase_amp.modelspec = 'Response ~ 1 + Pulse_num + Phase_theta *Amp_theta';
Phase_GLME.model.bphase_amp.modelspec = 'Response ~ 1 + Pulse_num + Phase_beta *Amp_beta';
Phase_GLME.model.lGphase_amp.modelspec = 'Response ~ 1 + Pulse_num + Phase_lG * Amp_lG';
Phase_GLME.model.mGphase_amp.modelspec = 'Response ~ 1 + Pulse_num + Phase_mG *Amp_mG';
Phase_GLME.model.hGphase_amp.modelspec = 'Response ~ 1 + Pulse_num + Phase_hG *Amp_hG';

%     Phase_GLME.model.all_model.modelspec ='Response ~ 1 + Pulse_num + Phase_delta + Phase_theta + Phase_beta + Phase_lG + Phase_mG +Phase_hG + Amp_delta+ Amp_theta+ Amp_beta + Amp_lG+ Amp_mG+ Amp_hG';

models = fieldnames(Phase_GLME.model);

% initialize outputs
for iMod = 1:length(models);
    Phase_GLME.model.(models{iMod}).err = zeros(length(cell_list), 1500); % needs to be zeros because error output will be added to this
    %     Phase_GLME.model.(models{iMod}).tstat = nan(length(cell_list), 11);
end

for iC = 1:length(cell_list)
    clear this_depth this_subject this_response this_phase this_amp
    
    this_cell = cell_list{iC};
    f_list = fieldnames(all_cells.(this_cell));
    % clear and setup GLM variables
    %clearvars('GLM*', 'Phase_GLME')
    
    GLM_depth = [];
    GLM_subject = [];
    GLM_response = [];
    % phases
    GLM_phase_delta = [];
    GLM_phase_theta = [];
    GLM_phase_beta = [];
    GLM_phase_lG = [];
    GLM_phase_mG =[];
    GLM_phase_hG = [];
    % amp
    GLM_amp_delta = [];
    GLM_amp_theta = [];
    GLM_amp_beta = [];
    GLM_amp_lG = [];
    GLM_amp_mG =[];
    GLM_amp_hG = [];
    
    
    
    % for cells with less than complete responses (faded or
    % unstable after some point
    cut_idx = all_cells.(this_cell).ExpKeys.goodTrials(1,1:2);
    
    % get the probe depth
    this_depth(1,1:(cut_idx(2)-cut_idx(1))+1) = all_cells.(this_cell).ExpKeys.tetrodeDepths;
    % subect ID
    this_subject(1,1:(cut_idx(2)-cut_idx(1))+1) = str2double(all_cells.(this_cell).ExpKeys.subject(2:3));
    
    %get the overall response.  Note frequency doesnt matter
    %since the overal responses were just copied for each freq
    %range.  This is leftover from an earlier version
    this_response = all_cells.(this_cell).(f_list{3}).resp(2,cut_idx(1):cut_idx(2));
    
    for iF = 1:length(f_list)
        if strcmp(f_list{iF}, 'ExpKeys') || strcmp(f_list{iF}, 'hdr')
            continue
        else
            % get the phase and amp for this frequency
            %             this_phase.(f_list{iF}) = all_cells.(this_cell).(f_list{iF}).amp(1,cut_idx(1):cut_idx(2));
            
            this_phase.(f_list{iF}) = categorical(all_cells.(this_cell).(f_list{iF}).amp(1,cut_idx(1):cut_idx(2)),'Ordinal',false);
            this_amp.(f_list{iF}) = all_cells.(this_cell).(f_list{iF}).amp(2,cut_idx(1):cut_idx(2));
            
        end
    end % end freq loop
    
    % collect and concatinate all the cells for the GLM
    if isequal(length(this_amp.(f_list{3})),length(this_phase.(f_list{3})),length(this_response),length(this_subject),length(this_depth))
        
        GLM_pulse_num = cut_idx(1):cut_idx(2);
        GLM_depth = this_depth;
        GLM_subject = this_subject;
        GLM_response = this_response;
        
        GLM_phase_delta =  this_phase.f_3_5;
        GLM_phase_theta =  this_phase.f_7_10;
        GLM_phase_beta = this_phase.f_15_25;
        GLM_phase_lG = this_phase.f_30_40;
        GLM_phase_mG = this_phase.f_40_60;
        GLM_phase_hG = this_phase.f_60_80;
        
        GLM_amp_delta =  this_amp.f_3_5;
        GLM_amp_theta =  this_amp.f_7_10;
        GLM_amp_beta = this_amp.f_15_25;
        GLM_amp_lG =  this_amp.f_30_40;
        GLM_amp_mG =  this_amp.f_40_60;
        GLM_amp_hG =  this_amp.f_60_80;
    else
        error('these are not equal')
    end
    
    
    Phase_GLME.tbl = table(GLM_subject', GLM_depth', GLM_response',GLM_pulse_num',...
        GLM_phase_delta', GLM_phase_theta', GLM_phase_beta', GLM_phase_lG', GLM_phase_mG', GLM_phase_hG', ...
        GLM_amp_delta', GLM_amp_theta', GLM_amp_beta', GLM_amp_lG', GLM_amp_mG', GLM_amp_hG', ...
        'VariableNames',{'SubjectID','Depth', 'Response','Pulse_num' 'Phase_delta', 'Phase_theta', 'Phase_beta', 'Phase_lG', 'Phase_mG', 'Phase_hG',...
        'Amp_delta', 'Amp_theta', 'Amp_beta', 'Amp_lG', 'Amp_mG', 'Amp_hG'});
    
    Phase_GLME.tbl.SubjectID = nominal(Phase_GLME.tbl.SubjectID);
    Phase_GLME.tbl.Response = logical(Phase_GLME.tbl.Response);
    
    
    %     amp_models = fieldnames(Amp_GLM.model);
    
    % make CV partitions
    nPleats = 1;
    nFolds = 10;
    for iPleat = nPleats:-1:1
        C{iPleat} = cvpartition(length(Phase_GLME.tbl.Response), 'KFold', nFolds);
    end
    
    %% run the GLMS w/ CV
    
    for iPleat = 1:nPleats
        
        for iFold = 1:C{iPleat}.NumTestSets
            
            fprintf('Pleat %d fold %d...\n', iPleat, iFold);
            
            % get idxs for training and testing set
            tr_idx = C{iPleat}.training(iFold); te_idx = C{iPleat}.test(iFold);
            
            %loop models
            for iMod = 1:length(models);
                if strcmp(models{iMod}, 'baseline')
                    base_line = fitglme(Phase_GLME.tbl(tr_idx,:), Phase_GLME.model.(models{iMod}).modelspec,'link','logit','Distribution','binomial');
                end
                model_out =  fitglme(Phase_GLME.tbl(tr_idx,:), Phase_GLME.model.(models{iMod}).modelspec,'link','logit','Distribution','binomial');
                %                 model_out.(models{iMod}); % [temp] just to print the outputs
                
                % hold the tstat
                Phase_GLME.model.(models{iMod}).tstat(iC,:,iFold) = model_out.Coefficients.tStat; % should initialize this, but that's a pain
                Phase_GLME.model.(models{iMod}).varnames{iC,:} = model_out.CoefficientNames;
                
                % hold the Fstat
                model_anova = anova(model_out);
                
                Phase_GLME.model.(models{iMod}).fstat = model_anova.FStat(end);
                Phase_GLME.model.(models{iMod}).fstat_var = model_anova.Term{end};
                % compare models
                if ~strcmp(models{iMod}, 'baseline')
                    model_comp = compare(model_out, base_line)
                    
                    
                    %                 this_err = model_out.predict(Phase_GLME.tbl(te_idx,:)); % gives same output as below
                    model_pred = predict(model_out, Phase_GLME.tbl(te_idx,:));
                    % other version
                    mse.(models{iMod})(:,iFold,iPleat) = nanmean((Phase_GLME.tbl.Response(te_idx,:) - model_pred).^2);
                    
                    model_pred = (model_pred - Phase_GLME.tbl.Response(te_idx)).^2;
                    Phase_GLME.model.(models{iMod}).err(iC,te_idx) = Phase_GLME.model.(models{iMod}).err(iC,te_idx) + (model_pred ./ nPleats)';
                    
                    
                end % models
                %             mse_mean
            end
            
        end
        for iMod = 1:length(models);
            % get all the tstats
            all_tstat.(models{iMod}){iC} = nanmean(Phase_GLME.model.(models{iMod}).tstat(iC,end,:));
            all_tstat_name{iMod} = Phase_GLME.model.(models{iMod}).varnames{iC}{end};
            
            all_mse.(models{iMod}){iC} = nanmean(sqrt(reshape(mse.(models{iMod}),1,numel(mse.(models{iMod})))));
            all_mse_relative.(models{iMod}){iC} = all_mse.baseline{iC} - all_mse.(models{iMod}){iC};
        end
    end
end

% plot out the MSE
for iMod = 1:length(models);
    all_mse_relative.mat(:,iMod) = cell2mat(all_mse_relative.(models{iMod}));
    all_tstat.mat(:,iMod) = cell2mat(all_tstat.(models{iMod}));
    
end
figure(210)
imagesc(all_mse_relative.mat)
set(gca,'xtick', 1:length(models), 'xticklabel', strrep(models,'_', ' '),'XTickLabelRotation', 45)
set(gca, 'yticklabel', [])
ylabel('cell number')
title('RMSE relative to baseline model')
saveas(gcf, [summary_dir 'GLM_MSE.png']);
saveas_eps('GLM_MSE', summary_dir)

figure(211)
imagesc(all_tstat.mat)
set(gca,'xtick', 1:length(all_tstat_name), 'xticklabel', strrep(all_tstat_name, '_', ' '),'XTickLabelRotation', 45)
title('mean t-stat')
saveas(gcf, [summary_dir 'GLM_tstat.png']);
saveas_eps('GLM_tstat', summary_dir)

figure(212)
plot(1:length(all_tstat_name), nanmean(all_tstat.mat))
set(gca,'xtick', 1:length(all_tstat_name), 'xticklabel', strrep(all_tstat_name, '_', ' '),'XTickLabelRotation', 45)


%% Compare across models
% check tStats
figure(234)
imagesc(Phase_GLME.model.(models{iMod}).tstat(:,:,1))
set(gca,'xtick', 1:length(model_out.CoefficientNames), 'xticklabel', strrep(model_out.CoefficientNames, '_', ' '),'XTickLabelRotation', 45)



for iMod = 1:length(models);
    Phase_GLME.model.(models{iMod}).err(Phase_GLME.model.(models{iMod}).err==0) =NaN;
    if strcmp(models{iMod},'baseline')
        continue;
    end
    
    model_v_baseline_err = Phase_GLME.model.baseline.err - Phase_GLME.model.(models{iMod}).err;
    celldiffmean = nanmean(model_v_baseline_err ,2);
    
    figure(iMod)
    plot(celldiffmean,'LineWidth',1); hold on;
    ylabel('Prediction improvement'); xlabel('cell #');
    title(models{iMod});
    
    
    
end

%         All_GLMS.model.(models{iMod}).pValue( 1:length(model_out.(models{iMod}).Coefficients.pValue)) = model_out.(models{iMod}).Coefficients.pValue; % should initialize this, but that's a pain
%         Phase_GLME.model.(models{iMod}).varnames = model_out.(models{iMod}).PredictorNames;
%
%Identify cells that pass and flag them phase mod
%                 if sum(model_out.(models{iMod}).Coefficients.pValue(2:end)<0.05) >0 && iMod > 1
%                     All_GLM_cells(iC, iMod-1) = 1;
%
%                     % add in LFP amplitude
%                     model_w_amp = fitglm(Phase_GLME.tbl, Amp_GLM.model.(amp_models{iMod}).modelspec,'link','logit','Distribution','binomial');
%
%                     % compare models, if better flag as phase x amp
%                     if model_w_amp.ModelCriterion.AIC <  model_out.(models{iMod}).ModelCriterion.AIC
%                         All_GLM_cells(iC, iMod-1) = 2;
%
%                     end % end compare
%                 end % end phase x amp




%% summary of all modulated cells
figure(111)
imagesc(All_GLM_cells)
set(gca, 'xticklabels', {'D', 'Th', 'B', 'lG', 'mG', 'hG'})
ylabel('cell number')



%%


%                    % run the GLM
%                    for iF = 1:length(f_list)
%                        if strcmp(f_list{iF}, 'ExpKeys') || strcmp(f_list{iF}, 'hdr')
%                            continue
%                        else
%
%                            % check for model vs intercept
%                                    Phase_GLME.glm = fitglm(Phase_GLME.tbl,'Response~1+Condition+(1|SubjectID)+(1|SessID)');
%                                 b = glmfit(x,y,'binomial','link','logit')
%
%
%                        end % check for expkeys field
%
%
%
%                    end % end freq loop
%


%% Statistical model

% Loop cells




% loop frequencys


% check baseline model with just an intercept


% glm(response ~ phase at freq)

% if sig, then compare to interaction with amplitude at that
% frequency




%%

%                     GLM_depth = [GLM_depth, this_depth];
%                     GLM_subject = [GLM_subject, this_subject];
%                     GLM_response = [GLM_response, this_response];
%
%                     GLM_phase_delta = [GLM_phase_delta, this_phase.f_3_5];
%                     GLM_phase_theta = [GLM_phase_theta, this_phase.f_7_10];
%                     GLM_phase_beta = [GLM_phase_beta, this_phase.f_15_25];
%                     GLM_phase_lG = [GLM_phase_lG, this_phase.f_30_40];
%                     GLM_phase_mG = [GLM_phase_mG, this_phase.f_40_60];
%                     GLM_phase_hG = [GLM_phase_hG, this_phase.f_60_80];
%
%                     GLM_amp_delta = [GLM_amp_delta, this_amp.f_3_5];
%                     GLM_amp_theta = [GLM_amp_theta, this_amp.f_7_10];
%                     GLM_amp_beta = [GLM_amp_beta, this_amp.f_15_25];
%                     GLM_amp_lG = [GLM_amp_lG, this_amp.f_30_40];
%                     GLM_amp_mG = [GLM_amp_mG, this_amp.f_40_60];
%                     GLM_amp_hG = [GLM_amp_hG, this_amp.f_60_80];

%%


% %% make a table of responsive cells, where they are, subject ID, ....
% Sub_list = {'M13', 'M14'};
% summary_table = cell(6,length(Sub_list)+1);
% summary_table{2,1} = 'n Resp Cells';
% for iSub = 1:length(Sub_list)
%     summary_table{1, iSub+1} = Sub_list{iSub};
%
% end
% for iC = 1:length(cell_list)
%     this_cell = cell_list{iC};
%     f_list = fieldnames(all_cells.(this_cell));
%
%     summary_table(
%
%
%
%
% end


%% plots
% figure(1)
% for iF = 1:length(f_list)
%     if strcmp(f_list{iF}, 'hdr') || strcmp(f_list{iF}, 'ExpKeys')
%         continue
%     else
%         subplot(2,ceil(length(f_list)-1)/2, iF)
%         imagesc(all_lat_sort(:,:,iF))
%         axis xy
%         xlabel('phase')
%         ylabel('cell id')
%         set(gca, 'xticklabel', phase_labels)
%         text(floor(length(phase_labels)/2),length(all_lat_sort)+1, f_list{iF}, 'fontsize', font_size)
%         colorbar
%
%     end
% end
% SetFigure([], gcf)
% set(gcf, 'position', [282   50  1200  720])
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.35, 0.98,['Latency to first spike across cells'], 'fontsize', font_size)
%
% saveas(gcf, [summary_dir 'Summary_latency_raw.png']);
% saveas_eps('Summary_latency_raw', summary_dir)
% %%
% figure(2)
% for iF = 1:length(f_list)
%     if strcmp(f_list{iF}, 'hdr') || strcmp(f_list{iF}, 'ExpKeys')
%         continue
%     else
%         subplot(2,ceil(length(f_list)-1)/2, iF)
%
%         imagesc(all_lat_v_shuf_sort(:,:,iF))
%         axis xy
%         xlabel('phase')
%         ylabel('cell id')
%         set(gca, 'xticklabel', phase_labels)
%         text(floor(length(phase_labels)/2),length(all_lat_v_shuf_sort)+1, f_list{iF}, 'fontsize', font_size)
%         colorbar
%
%     end
% end
% SetFigure([], gcf)
% set(gcf, 'position', [282   50  1200  720])
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.35, 0.98,['Latency to first spike across cells (above 2d of shuf)'], 'fontsize', font_size)
%
% saveas(gcf, [summary_dir 'Summary_latency_zshuf.png']);
% saveas_eps('Summary_latency_z_shuf', summary_dir)
% %%
% % same but for zscore
% figure(3)
% for iF = 1:length(f_list)
%     if strcmp(f_list{iF}, 'hdr') || strcmp(f_list{iF}, 'ExpKeys')
%         continue
%     else
%         subplot(2,ceil(length(f_list)-1)/2, iF)
%         imagesc(z_lat(:,:,iF))
%         axis xy
%         xlabel('phase')
%         ylabel('cell id')
%         set(gca, 'xticklabel', phase_labels)
%         text(floor(length(phase_labels)/2),length(all_lat)+1, f_list{iF}, 'fontsize', font_size)
%
%         colorbar
%         caxis([-2.5 2.5])
%     end
% end
% SetFigure([], gcf)
% set(gcf, 'position', [282   50  1200  720])
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.35, 0.98,['Zscore first spike latency across cells'], 'fontsize', font_size)
% saveas(gcf, [summary_dir 'Summary_zscore.png']);
% saveas_eps('Summary_zscore', summary_dir)
%
% %% count version
%
% figure(3)
% for iF = 1:length(f_list)
%     if strcmp(f_list{iF}, 'hdr') || strcmp(f_list{iF}, 'ExpKeys')
%         continue
%     else
%         subplot(2,ceil(length(f_list)-1)/2, iF)
%         imagesc(all_count_sort(:,:,iF))
%         axis xy
%         xlabel('phase')
%         ylabel('cell id')
%         set(gca, 'xticklabel', phase_labels)
%         text(floor(length(phase_labels)/2),length(all_count_sort)+1, f_list{iF}, 'fontsize', font_size)
%
%         colorbar
%         %         caxis([-2.5 2.5])
%     end
% end
% SetFigure([], gcf)
% set(gcf, 'position', [282   50  1200  720])
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.35, 0.98,['nSpikes per stim across cells'], 'fontsize', font_size)
% saveas(gcf, [summary_dir 'Summary_count.png']);
% saveas_eps('Summary_count', summary_dir)
%
%
% %% response version
% figure(5)
% for iF = 1:length(f_list)
%     if strcmp(f_list{iF}, 'hdr') || strcmp(f_list{iF}, 'ExpKeys')
%         continue
%     else
%         subplot(2,ceil(length(f_list)-1)/2, iF)
%         imagesc(all_resp_sort(:,:,iF))
%         axis xy
%         xlabel('phase')
%         ylabel('cell id')
%         set(gca, 'xticklabel', phase_labels)
%         text(floor(length(phase_labels)/2),length(all_resp_sort)+1, f_list{iF}, 'fontsize', font_size)
%
%         colorbar
%         %         caxis([-2.5 2.5])
%     end
% end
% SetFigure([], gcf)
% set(gcf, 'position', [282   50  1200  720])
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.35, 0.98,['Response per stim across cells'], 'fontsize', font_size)
% saveas(gcf, [summary_dir 'Summary_resp.png']);
% saveas_eps('Summary_resp', summary_dir)
%
% % zscore
% figure(6)
% for iF = 1:length(f_list)
%     if strcmp(f_list{iF}, 'hdr') || strcmp(f_list{iF}, 'ExpKeys')
%         continue
%     else
%         subplot(2,ceil(length(f_list)-1)/2, iF)
%         imagesc(z_resp(:,:,iF))
%         axis xy
%         xlabel('phase')
%         ylabel('cell id')
%         set(gca, 'xticklabel', phase_labels)
%         text(floor(length(phase_labels)/2),length(all_lat)+1, f_list{iF}, 'fontsize', font_size)
%
%         colorbar
%         %         caxis([-2.5 2.5])
%     end
% end
% SetFigure([], gcf)
% set(gcf, 'position', [282   50  1200  720])
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.35, 0.98,['Zscore response per stim across cells'], 'fontsize', font_size)
% saveas(gcf, [summary_dir 'Summary_resp_z.png']);
% saveas_eps('Summary_resp_z', summary_dir)
%
% % relative to all
% % zscore
% figure(7)
% for iF = 1:length(f_list)
%     if strcmp(f_list{iF}, 'hdr') || strcmp(f_list{iF}, 'ExpKeys')
%         continue
%     else
%         subplot(2,ceil(length(f_list)-1)/2, iF)
%         imagesc(all_resp_v_all_sort(:,:,iF))
%         axis xy
%         xlabel('phase')
%         ylabel('cell id')
%         set(gca, 'xticklabel', phase_labels)
%         text(floor(length(phase_labels)/2),length(all_resp_v_all_sort)+1, f_list{iF}, 'fontsize', font_size)
%
%         colorbar
%         %         caxis([-2.5 2.5])
%     end
% end
% SetFigure([], gcf)
% set(gcf, 'position', [282   50  1200  720])
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.35, 0.98,['Relative response per stim across cells'], 'fontsize', font_size)
% saveas(gcf, [summary_dir 'Summary_resp_v_all.png']);
% saveas_eps('Summary_resp_v_all', summary_dir)
%%


