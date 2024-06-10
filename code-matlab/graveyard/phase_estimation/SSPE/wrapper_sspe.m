function phases = wrapper_sspe(wsz, fq_exp, csc, fbands, filt_phase, Fs, nSamples)
% This is a wrapper function that returns true and estimated phase values
% as obtained by the resonant oscillator method
% Let's try to optimize without caring about initial frequencies
%     Input:
%         wsz: window size in tenths of a second (for faster optimization)
%         fq_exp: nx1 array of sigma_freq exponentials
%         csc: mvdmlab TSD struct
%         fbands: nx2 array of frequency bands [[lfq1 hfq];[lfq1 hfq2]; ..., lfqn hfqn]]
%         filt_phase: cell array where each cell contians hilbert transformed phases of 'csc'
%         Fs: sampling frequency of csc
%         nSamples: Number of samples to be test this method on
%     Output:
%         phases: 1x2 cell array, phases{1} has the true phases and phases{2} has the estimated phases

    wsz = wsz/10;
    min_start = ceil(wsz*Fs);
    nEnds = randi(length(csc.data) - min_start, nSamples, 1) + min_start;
    nStarts = nearest_idx3(csc.tvec(nEnds) - wsz, csc.tvec);
    % hack to make sure the 99% of the smallest window size is at least Fs
    if ceil(0.99*min(nEnds - nStarts)) <= round(Fs)
        nEnds = nEnds + ceil(round(Fs) - 0.99*min(nEnds - nStarts));
    end
    true_phase = cell2mat(cellfun(@(x) x(nEnds), filt_phase, 'UniformOutput', false));
    output_phase = zeros(length(fq_exp), nSamples);
    for iS = 1:nSamples
       this_sample = csc.data(nStarts(iS):nEnds(iS));  
       this_cfg.freqs = cellfun(@mean, fbands);
       this_cfg.freqBands = fbands;
       this_cfg.Fs = Fs;
       this_cfg.ampVec = repmat(0.99, 1, length(this_cfg.freqs));
       this_cfg.sigmaObs = 1;
       this_cfg.window = round(length(this_sample)*0.99);
       this_cfg.sigmaFreqs = fq_exp;
       [omega, phase, phase_bounds, fullX] = causalPhaseEM_MKmdl_all(this_sample, this_cfg);
       output_phase(:,iS) = phase(:,end); % The last sample's phase
    end
    phases = cell(1,2);
    phases{1} = true_phase;
    phases{2} = output_phase;
end
