function [data_all] = convert_tsd_to_ft(csc_in)

% FT_READ_NEURALYNX_INTERP reads a cell-array of NCS files and performs interpolation
% and gap-filling of the timestamp axis to correct for potentially different offsets
% across channels, and potential gaps within the recordings. All samples are being 
% read out.
%
% Inputs:
%  input FNAME is the list of CSC filenames
%  All of these channels should have the sample sampling frequency.
% Outputs:
%  data_all is a raw data structure containing interpolated data and NaNs at the gaps, 
%  based on all the available samples in a recording.
%

ts = csc_in.tvec;  
min_all = ts(1);
max_all = ts(end); 
ts1   = ts(1);
md = mode(diff(double(ts-ts1)));

% take the mode of the modes and construct the interpolation axis
% the interpolation axis should be casted in doubles
mode_dts  = mode(md);
rng       = double(max_all-min_all); % this is small num, can be double
tsinterp  = [0:mode_dts:rng]; % the timestamp interpolation axis

cfg         = [];
cfg.dataset = csc_in.label{1};
data        = ft_preprocessing(cfg); % Let's see if this has any effect on the pipeline 
% Replace the required fields
data.time{1} = ts';
data.trial{1} = csc_in.data;
data.hdr.Fs =  csc_in.cfg.hdr{1}.SamplingFrequency;

% original timestamaps in doubles, with the minimum ts subtracted
ts          = double(ts-ts(1)); 
  
% check if there are gaps to correct
gaps     = find(diff(ts)>2*mode_dts); % skips at least a sample
if isempty(gaps)
    fprintf('there are no gaps and all channels start and end at same time, no interpolation performed\n');
else  
    % interpolate the data
    [unique_ts,keep_idx] = unique(ts);
    if length(ts) ~= length(unique_ts)
    disp('Duplicate timestamps removed!'); 
        ts = unique(ts);
        data.trial{1} = data.trial{1}(keep_idx);
    end
    datinterp   = interp1(ts, data.trial{1}, tsinterp);

    % Interpolate the signal in the gaps
    gaps = find(diff(ts)>2*mode_dts); % skips at least a sample
    if ~isempty(gaps)
            fprintf('Largest gap is %.2f sec\n', max(diff(ts))/data.hdr.Fs);
        for igap = 1:length(gaps)
          sel = find(tsinterp < ts(gaps(igap)+1) & tsinterp > ts(gaps(igap)));
          datinterp(sel) = interp1(tsinterp(1:sel(end)-1),datinterp(1:sel(end)-1), tsinterp(sel), 'spline', 'extrap');
        end
    end 
    % update the FieldTrip data structure
    data.trial{1} = datinterp; clear datinterp
    data.time{1}  = [0:length(tsinterp)-1].*(1./data.hdr.Fs);
end

data_all = data;
clear data

% correct the header information and the sampling information
data_all.hdr.FirstTimeStamp     = min_all;
data_all.hdr.LastTimeStamp      = uint64(tsinterp(end)) + min_all;
data_all.hdr.TimeStampPerSample = mode_dts;
len = length(tsinterp);
data_all.hdr.nSamples           = len;
data_all.sampleinfo             = [1 len];
data_all.cfg.trl(2)             = len;
data_all.fsample                = data_all.hdr.Fs;