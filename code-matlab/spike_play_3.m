cd('D:\Dropbox (Dartmouth College)\manish_data\M075\artifacts');
cfg_in.uint= '64';
S = LoadSpikes([cfg_in]);
evs = LoadEvents([]);
all_OFF_times = evs.t{9};

%% pre-trial stim
p_start = evs.t{15};
p_end = evs.t{5};
p_on = evs.t{10};
p_off = all_OFF_times(all_OFF_times > p_start & all_OFF_times < p_end);
p_iv = iv(p_start, p_end);
p_S = restrict(S, p_iv);
pre_offsets = cell(1, length(p_S.t));

for i = 1:length(p_S.t)
   idx = nearest_idx3(p_S.t{i},p_on);
   this_offsets = 1e3*(p_S.t{i} - p_on(idx)); % Convert to milliseconds
   % Some offsets are less than 0 and let's ignore them for now
   this_offsets = this_offsets(this_offsets >= 0);
   pre_offsets{i} = this_offsets;
end

%% trial stim
p_start = evs.t{5};
p_end = evs.t{2};
p_on = evs.t{11};
p_off = all_OFF_times(all_OFF_times > p_start & all_OFF_times < p_end);
p_iv = iv(p_start, p_end);
p_S = restrict(S, p_iv);
trial_offsets = cell(1, length(p_S.t));

for i = 1:length(p_S.t)
   idx = nearest_idx3(p_S.t{i},p_on);
   this_offsets = 1e3*(p_S.t{i} - p_on(idx)); % Convert to milliseconds
   % Some offsets are less than 0 and let's ignore them for now
   this_offsets = this_offsets(this_offsets >= 0);
   trial_offsets{i} = this_offsets;
end

%% post-trial stim
p_start = evs.t{6};
p_end = evs.t{14};
p_on = evs.t{12};
p_off = all_OFF_times(all_OFF_times > p_start & all_OFF_times < p_end);
p_iv = iv(p_start, p_end);
p_S = restrict(S, p_iv);
post_offsets = cell(1, length(p_S.t));

for i = 1:length(p_S.t)
   idx = nearest_idx3(p_S.t{i},p_on);
   this_offsets = 1e3*(p_S.t{i} - p_on(idx)); % Convert to milliseconds
   % Some offsets are less than 0 and let's ignore them for now
   this_offsets = this_offsets(this_offsets >= 0);
   post_offsets{i} = this_offsets;
end

%% post_long stim

p_start = evs.t{14};
p_end = evs.t{8};
p_on = evs.t{13};
p_off = all_OFF_times(all_OFF_times > p_start & all_OFF_times < p_end);
p_iv = iv(p_start, p_end);
p_S = restrict(S, p_iv);
long_offsets = cell(1, length(p_S.t));

for i = 1:length(p_S.t)
   idx = nearest_idx3(p_S.t{i},p_on);
   this_offsets = 1e3*(p_S.t{i} - p_on(idx)); % Convert to milliseconds
   % Some offsets are less than 0 and let's ignore them for now
   this_offsets = this_offsets(this_offsets >= 0);
   long_offsets{i} = this_offsets;
end

%% Plotting
for i = 1:length(p_S.t)
    close all;
    fig = figure('WindowState', 'maximized');
    subplot(4,1,1);
    hist(pre_offsets{i});
    title('100 Pre trial');
    subplot(4,1,2);
    hist(trial_offsets{i});
    title('1000 trial');
    subplot(4,1,3);
    hist(post_offsets{i});
    title('100 Post trial');
    subplot(4,1,4);
    hist(long_offsets{i});
    title('25 Post trial long');
    sgtitle(p_S.label{i}(1:4));
    WriteFig(fig, p_S.label{i}(1:4), 1);
    close;
end