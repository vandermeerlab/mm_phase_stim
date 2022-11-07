
function statedep_summary

all_lat_dir = 'E:\temp_phase_stim\test';
% all_fig_dir = 'E:\temp_phase_stim\ED\M16-2019-02-18_vStr_4p0_light_cells_TT4_min\figures';
summary_dir = 'E:\temp_phase_stim\';

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesFontWeight','bold')

PARAMS.Good_cells = {'M16_2019_02_18_TT04_01_Good'};
font_size = 18;
%% generate a latency and count summary
cd(all_lat_dir);
sess_list = dir(all_lat_dir);
sess_list = sess_list(3:end);

for iSess = 1:length(sess_list)
    load(sess_list(iSess).name);
    
    these_cells = fieldnames(out);
    for iC = 1:length(these_cells)
        %For now assume all cells are good
        this_cell = these_cells{iC};
        freq_list = fieldnames(out.(this_cell));
        all_cells.(this_cell) = out.(this_cell);
%         all_cells.([sess_list(iSess).name(1:14) '_' this_cell]) = out.(this_cell);
%         if ismember([sess_list(iSess).name(1:14) '_' this_cell], PARAMS.Good_cells) % check if this is a 'good' cell from above approvedd list
%             
%             freq_list = fieldnames(out.(this_cell));
%             
%             
%             all_cells.([sess_list(iSess).name(1:14) '_' this_cell]) = out.(this_cell);
%             
%             % subject hdr
%             %             all_cells.(this_cell).hdr.name = sess_list(iSess).name(1:14);
%             %             all_cells.(this_cell).hdr.date = sess_list(iSess).name(5:14);
%             %             all_cells.(this_cell).hdr.subject = sess_list(iSess).name(1:3);
%         else
%             disp(['Skippking ' sess_list(iSess).name(1:14) '_' this_cell '.  Not a ''good'' cell'])
%         end     
    end 
end

%%
nShuf = 10000;
cell_list = fieldnames(all_cells);
isMD = false(size(cell_list));
all_lat =[];
all_count = [];
all_resp_ratio = [];
for iC = 1:length(cell_list)
    this_cell = cell_list{iC};
    f_list = fieldnames(all_cells.(this_cell));
    if isfield(all_cells.(this_cell).ExpKeys, 'tetrodeDepths')
       cell_depths(iC) = all_cells.(this_cell).ExpKeys.tetrodeDepths;
    else
       cell_depths(iC) = all_cells.(this_cell).ExpKeys.probeDepth;
       isMD(iC) = true;
    end
     
    for iF = 1:length(f_list)
        if strcmp(f_list{iF}, 'ExpKeys') || strcmp(f_list{iF}, 'hdr')
            continue
        else
            x_phase_resp(iC, iF) = nanmean(all_cells.(this_cell).(f_list{iF}).resp(2,:));
               
            %get the shuffle for the response max-min
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
            
            
            
            
            % get the
            
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
cell_names = fieldnames(all_cells);
cell_names = cellfun(@(x) extractAfter(x,15), cell_names, 'UniformOutput', false);
[depth_sort, depth_ord] = sort(cell_depths, 'ascend');
cell_names = cell_names(depth_ord);
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
isMD = isMD(depth_ord);

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
cut_off = 0;
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

%% Manish Results Figure
fig = figure;
%Plot 2-5 Hz
ax1 = subplot(3,2,[1,3]);
iF = find(strcmp('f_2_5', f_list));
p1 = imagesc(all_resp_v_shuf_sort(:,:,iF));
caxis([-2.5 2.5])
ax1.Box = 'off';
ax1.XTick = {};
ax1.TickDir = 'out';
ax1.YAxis.Label.String = 'Cells';
ax1.YAxis.FontSize = 16;
ax1.XAxis.Label.String = 'Phase';
ax1.XAxis.FontSize = 16;
ax1.YTick = [1:length(all_resp_v_shuf_sort(:,:,iF))];
ax1.YTickLabel = labels;
ax1.Title.String = '2 - 5 Hz';
ax1.Title.FontSize = 24;

%Plot 6-10 Hz
ax2 = subplot(3,2,[2,4]);
iF = find(strcmp('f_6_10', f_list));
p2 = imagesc(all_resp_v_shuf_sort(:,:,iF));
caxis([-2.5 2.5])
ax2.Box = 'off';
ax2.XTick = {};
ax2.TickDir = 'out';
ax2.YAxis.Label.String = 'Cells';
ax2.YAxis.FontSize = 16;
ax2.XAxis.Label.String = 'Phase';
ax2.XAxis.FontSize = 16;
ax2.YTick = [1:length(all_resp_v_shuf_sort(:,:,iF))];
ax2.YTickLabel = labels;
ax2.Title.String = ' 6 - 10 Hz';
ax2.Title.FontSize = 24;

%Plot 25-55 Hz
ax3 = subplot(3,2,5);
iF = find(strcmp('f_25_55', f_list));
p3 = imagesc(all_resp_v_shuf_sort(~isMD,:,iF));
caxis([-2.5 2.5])
ax3.Box = 'off';
ax3.XTick = {};
ax3.TickDir = 'out';
ax3.YAxis.Label.String = 'Cells';
ax3.YAxis.FontSize = 16;
ax3.XAxis.Label.String = 'Phase';
ax3.XAxis.FontSize = 16;
ax3.YTick = [1:length(all_resp_v_shuf_sort(~isMD,:,iF))];
ax3.YTickLabel = labels;
ax3.Title.String = '25 - 55 Hz';
ax3.Title.FontSize = 24;

%Plot 65-100 Hz
ax4 = subplot(3,2,6);
iF = find(strcmp('f_65_100', f_list));
p4 = imagesc(all_resp_v_shuf_sort(~isMD,:,iF));
caxis([-2.5 2.5])
ax4.Box = 'off';
ax4.XTick = {};
ax4.TickDir = 'out';
ax4.YAxis.Label.String = 'Cells';
ax4.YAxis.FontSize = 16;
ax4.XAxis.Label.String = 'Phase';
ax4.XAxis.FontSize = 16;
ax4.YTick = [1:length(all_resp_v_shuf_sort(~isMD,:,iF))];
ax4.YTickLabel = labels;
ax4.Title.String = '65 - 100 Hz';
ax4.Title.FontSize = 24;

%% Manish summary response plot


f_cor = linspecer(2);
figure(9)
freq_colors = [];
%First scatter Manish Data
for iC =1:size(all_resp_ratio_sort(isMD,:),1)
    for iF = 1:2
        freq_colors{iC,iF} = f_cor(iF,:);
        freq(iC, iF) = iF;
        cell_id(iC,iF) = iC;
    end
end
resp_ratio_1d = reshape(all_resp_ratio_sort(isMD,:), 1, numel(all_resp_ratio_sort(isMD,:)));
resp_all_1d = reshape(x_phase_resp_sort(isMD,:), 1, numel(x_phase_resp_sort(isMD,:)));
freq_colors_1d = reshape(freq_colors,1, numel(freq_colors));
% freq_id_1d = reshape(freq_vals(isMD),1, numel(freq_vals(isMD)));
depth_all_1d = reshape(repmat(depth_sort(isMD)',1,size(x_phase_resp_sort(isMD,:),2)), 1, numel(x_phase_resp_sort(isMD,:)));
% cell_id_1d = reshape(cell_id,1, numel(cell_id));

% generate figure 

for iC = 1:34 % this is really_hacky %length(resp_all_1d)
    hold on
    if depth_all_1d(iC) >3.5
        S(iC) = scatter(resp_all_1d(iC), resp_ratio_1d(iC),300,freq_colors_1d{iC}, 'filled');
    else
        scatter(resp_all_1d(iC), resp_ratio_1d(iC),300,freq_colors_1d{iC}, 'filled', 'd');
    end
end
% Then scatter Eric Data
f_cor = linspecer(4);
freq_colors = [];
for iC =1:size(all_resp_ratio_sort(~isMD,:),1)
    for iF = 1:4
        freq_colors{iC,iF} = f_cor(iF,:);
        freq(iC, iF) = iF;
        cell_id(iC,iF) = iC;
    end
end
resp_ratio_1d = reshape(all_resp_ratio_sort(~isMD,:), 1, numel(all_resp_ratio_sort(~isMD,:)));
resp_all_1d = reshape(x_phase_resp_sort(~isMD,:), 1, numel(x_phase_resp_sort(~isMD,:)));
freq_colors_1d = reshape(freq_colors,1, numel(freq_colors));
% freq_id_1d = reshape(freq_vals(isMD),1, numel(freq_vals(isMD)));
depth_all_1d = reshape(repmat(depth_sort(~isMD)',1,size(x_phase_resp_sort(~isMD,:),2)), 1, numel(x_phase_resp_sort(~isMD,:)));
% cell_id_1d = reshape(cell_id,1, numel(cell_id));

% generate figure 

for iC = 1:length(resp_all_1d)
    if depth_all_1d(iC) >3.5
        S(iC) = scatter(resp_all_1d(iC), resp_ratio_1d(iC),300,freq_colors_1d{iC}, 'filled');
    else
        scatter(resp_all_1d(iC), resp_ratio_1d(iC),300,freq_colors_1d{iC}, 'filled', 'd');
    end
end

%
xlabel('Overall Response Probability')
ylabel('Z-scored response (max phase/min phase)')
%ylabel('ratio of max phase / min phase response')
h = zeros(size(all_resp_ratio,2), 1);
for iF = 1:size(all_resp_ratio,2)
    h(iF) = plot(NaN,NaN,'o', 'color', f_cor(iF, :), 'MarkerFaceColor', f_cor(iF,:));
end
lgd = legend(h, leg_freq, 'location', 'NorthEastOutside');
legend('boxoff')
lgd.FontSize = 18;
q0 = yline(1.96, 'color', 'black');
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
h(1) = plot(NaN,NaN,'o', 'color', 'k', 'MarkerFaceColor', 'k');
h(2) = plot(NaN,NaN,'d', 'color', 'k', 'MarkerFaceColor', 'k');
ax = gca;
ax.Legend.String = {'2 - 5 Hz', '6 - 10 Hz', '25 - 55 Hz', '65 - 100 Hz', 'vStr', 'dStr'};
ax.TickDir = 'out';
ax.Box = 'off';
ax.XTick = [0 0.25 0.5 0.75 1];
ax.YTick = [-2 2 6 12];
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
ax.XLim =[0 1];
xlim([cut_off 1])
ylim([-2 12])
%%
q0 = xline(1.96, 'color', 'black');
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
line([0 1.5], [1.96, 1.96], 'color', [0.3 0.3 0.3])
line([0 1.5], [3.1, 3.1], 'color', [0.3, 0.3, 0.3])

SetFigure([], gcf)
saveas(gcf, [summary_dir 'Summary_resp_ratio.fig']);
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
%         set(gca,'ytick', [1:length(all_resp_v_shuf_sort(:,:,iF))], 'yticklabel',labels, 'fontsize',14); % [num2str(cell_depths) ' - ' num2str(floor(x_phase_resp(:,iF)*100))]
        set(gca,'ytick', [1:length(all_resp_v_shuf_sort(:,:,iF))], 'yticklabel',strcat(cell_names,'_',labels'), 'fontsize',10); % [num2str(cell_depths) ' - ' num2str(floor(x_phase_resp(:,iF)*100))]
%         text(floor(length(phase_labels)/2),length(all_resp_v_shuf_sort)+1, f_list{iF}, 'fontsize', font_size)
        title(f_list{iF});
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
% SetFigure([], gcf)
% set(gca, 'fontsize', 8);
% for iSub = 1:6
%     subplot(2,3,iSub)
% set(gca, 'yticklabel',labels, 'fontsize',12)
% end
set(gcf, 'position', [282   50  1100  800])
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.35, 0.98,['Zscore v Shuf per stim across cells'], 'fontsize', font_size)
saveas(gcf, [summary_dir 'Summary_resp_v_shuf.fig']);
% saveas_eps('Summary_resp_v_shuf', summary_dir)
%%
% make a plot just for the color bar
figure(888)
subplot(3,1,1:3)
imagesc([-2.5 2.5])
% xlabel('phase');
% ylabel('cell');
% set(gca, 'xticklabel', []);
% set(gca,'ytick', [1:length(all_resp_v_shuf_sort(:,:,3))], 'yticklabel',labels, 'fontsize',14 ); % [num2str(cell_depths) ' - ' num2str(floor(x_phase_resp(:,iF)*100))]
% text(floor(length(phase_labels)/2),length(all_resp_v_shuf_sort)+1, f_list{3}, 'fontsize', font_size)
hp = colorbar;
caxis([-2.5 2.5])
% SetFigure([], gcf)
% set(gcf, 'position', [282   50  1100  800])
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
hp.TickDirection = 'out';
hp.FontSize = 18;
hp.Ticks = [-2.5 0 2.5];
% saveas_eps('Summary_resp_v_shuf_Cbar', summary_dir)
%% make a response plot
f_cor = linspecer(size(all_resp_ratio,2));
figure(10)
freq_colors = [];
for iC =1:size(all_resp_ratio_sort,1)
    for iF = 1:size(all_resp_ratio_sort,2)
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

xlim([cut_off 1])
ylim([-10 10])
line([0 1.5], [1.96, 1.96], 'color', [0.3 0.3 0.3])
line([0 1.5], [3.1, 3.1], 'color', [0.3, 0.3, 0.3])

SetFigure([], gcf)
saveas(gcf, [summary_dir 'Summary_resp_ratio.fig']);
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
