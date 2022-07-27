%Script to fix funky LFP recordings due to filter settings
cd('E:\Dropbox (Dartmouth College)\manish_data\M320\M320-2022-05-26');
file_list = {'LFP3_OG.ncs', 'LFP15_OG.ncs', 'LFP18_OG.ncs', ...
    'R1_OG.ncs', 'R4_OG.ncs'};

%%
% Extract All the records from the CSC File
for iF = 1:length(file_list)
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, ...
        Samples, Header] = Nlx2MatCSC(file_list{iF}, [1 1 1 1 1], 1, 1, [] );

    % LoadCSC only reads the 'valid samples' of a signal, so we should subtract
    % only the mean of these valid samples from the data;
    all_valid = (NumberOfValidSamples == 512);
    all_valid_sum = sum(sum(Samples(:,all_valid)));
    some_valid = find(~all_valid);
    some_valid_sum = 0;
    for iS = 1:length(some_valid)
        some_valid_sum = some_valid_sum + ...
            sum(Samples(NumberOfValidSamples(some_valid(iS)), some_valid(iS)));
    end
    signal_mean = (all_valid_sum + some_valid_sum)/sum(NumberOfValidSamples);
    
    % Subtract the mean from the signal
    new_Samples = Samples - repmat(signal_mean, size(Samples));

    % Write the file back after modifications
    new_filename = strcat(extractBefore(file_list{iF}, '_OG.ncs'), '.ncs');
    Mat2NlxCSC(new_filename, 0, 1, 1, [1 1 1 1 1 1], Timestamps, ChannelNumbers, ...
        SampleFrequencies, NumberOfValidSamples, new_Samples, Header);
end