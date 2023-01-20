%% Script to play around with correlograms of spectrograms

top_dir = 'E:\Dropbox (Dartmouth College)\manish_data\';
mice = {'M016', 'M017', 'M018', 'M019', 'M020', 'M074', 'M075', 'M077', 'M078', 'M235', 'M265', 'M295', 'M320', 'M319', 'M321', 'M325'};

all_cm = {};
all_labels = {};
for iM  = 1:length(mice)
    all_sess = dir(strcat(top_dir, mice{iM}));
    sid = find(arrayfun(@(x) contains(x.name, mice{iM}), all_sess));
    for iS = 1:length(sid)
        this_dir = strcat(top_dir, mice{iM}, '\', all_sess(sid(iS)).name);
        cd(this_dir);
        LoadExpKeys;
        if isempty(ExpKeys.goodCell)
            continue
        end
        clear ExpKeys
        all_cm{length(all_cm)+1} = doStuff;
        all_labels{length(all_labels)+1} = all_sess(sid(iS)).name;
    end
end

%%
for i = 1:length(all_cm)
    imagesc(all_cm{i});
    title(all_labels{i}, 'Interpreter','none');
    caxis([0 0.5]);
    dummy = 1;
end
%%
gcm = zeros(length(all_cm), 120,120);
for i = 1:length(all_cm)
    gcm(i,:,:) = all_cm{i};
end
%%
gmean = squeeze(mean(gcm,1));
gsd = squeeze(std(gcm,0,1));
gpct = squeeze(prctile(gcm,50,1));
imagesc(gmean);
f_list = {[2 5], [6 10], [12 30], [30, 55], [65 100]};
for iF = 1:length(f_list)
    F1 = f_list{iF}(1);
    F2 = f_list{iF}(2);
    rectangle('Position',[F1,F1,F2-F1,F2-F1],'EdgeColor', 'red', 'LineWidth', 3)
end
caxis([0 0.3]);

%%
odir = 'E:\Dropbox (Dartmouth College)\AnalysisResults\phase_stim_results\spectrocorrelogram\';
for iS = 1:length(all_cm)
    this_fig = figure;
    cm = all_cm{iS};
    imagesc(cm);
    caxis([0 0.3])
    hold on
    for iF = 1:length(f_list)
        F1 = f_list{iF}(1);
        F2 = f_list{iF}(2);
        rectangle('Position',[F1,F1,F2-F1,F2-F1],'EdgeColor', 'red', 'LineWidth', 3)
        title(all_labels{iS}, 'Interpreter','none');
        total = 0;
        hi = 0;
        for i = F1:F2-1
            for j = i+1:F2
                total = total+1;
                if (cm(i,j) >= gpct(i,j))
                    hi = hi+1;
                    scatter(i,j, 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0, 'SizeData', 10); 
                end
            end
        end
        crossed(iF) = hi/total;
    end
    print(this_fig, '-dpng', strcat(odir,all_labels{iS},'-SC'));
    close;
end

%%
function cm = doStuff
    load('spectrogram_data.mat');
    P = spec_data.P;
    cm = corr(P',P');
end
%%
% %%
% % cd('E:\Dropbox (Dartmouth College)\EC_State_inProcess\M018\M018-2019-04-10-dStr_3p8_light_cells_TT4_min\')
% % cd('E:\Dropbox (Dartmouth College)\manish_data\M325\M325-2022-07-28')
% load('spectrogram_data.mat')
% P = spec_data.P;
% P2 = 10*log10(P);
% F = spec_data.F;
% 
% % Create covariance matrices
% cm = corr(P',P');
% 
% % Generate 1000 shuffles and their correlation matrices
% nshuf = 1000;
% rcm = zeros(nshuf, size(cm,1), size(cm,2));
% for i = 1:nshuf
%     Pr = P;
%     for iR = 1:size(Pr,1)
%         Pr(iR,:) = Pr(iR, randperm(size(Pr,2))); % shuffle all rows
%     end
%     rcm(i,:,:) = corr(Pr',Pr');
% end
% %%
% rmean = squeeze(mean(rcm,1));
% rsd = squeeze(std(rcm,0,1));
% rthresh = rmean + 3*rsd;
% p_thresh = 99;
% rpct = squeeze(prctile(rcm,p_thresh,1));
% 
% %%
% figure;
% imagesc(cm);
% colorbar;
% hold on;
% caxis([0 0.3]);
% f_list = {[2 5], [6 10], [12 30], [30, 55], [65 100]};
% crossed = zeros(size(f_list));
% for iF = 1:length(f_list)
%     F1 = F(nearest_idx3(f_list{iF}(1),F));
%     F2 = F(nearest_idx3(f_list{iF}(2),F));
%     rectangle('Position',[F1,F1,F2-F1,F2-F1],'EdgeColor', 'red', 'LineWidth', 3)
%     total = 0;
%     hi = 0;
%     for i = F1:F2-1
%         for j = i+1:F2
%             total = total+1;
%             if (cm(i,j) >= rthresh(i,j))
%                 hi = hi+1;
%                 scatter(i,j, 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0, 'SizeData', 10); 
%             end
%         end
%     end
%     crossed(iF) = hi/total;
% end
%% Attempt at unsupervised_clustering
% threshold = 0.2;
% cm2 = cm1>0.2;  
% for i = 1:119
%     for j = i+1:120
% 
%         this_row = cm1(i, )
%     end
% end
