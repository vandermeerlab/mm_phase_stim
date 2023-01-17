function phases = wrapper_echt(wsz, csc, filt_phase, fbands, Fs, nSamples)
% This is a wrapper function that returns true and estimated phase values
% as obtained by the echt method
%     Input:
%         wsz: window size in tenths of a second (for faster optimization)
%         csc: mvdmlab TSD struct
%         filt_phase: cell array where each cell contians hilbert transformed phases of 'csc'
%         fbands: nx2 array of frequency bands [[lfq1 hfq];[lfq1 hfq2]; ..., lfqn hfqn]]
%         Fs: sampling frequency of csc
%         nSamples: Number of samples to be test this method on
%     Output:
%         phases: 1x2 cell array, phases{1} has the true phases and phases{2} has the estimated phases
% 

    wsz = wsz/10;
    min_start = ceil(wsz*Fs);
    nEnds = randi(length(csc.data) - min_start, nSamples, 1) + min_start;
    nStarts = nearest_idx3(csc.tvec(nEnds) - wsz, csc.tvec);
    true_phase = zeros(length(fbands), nSamples);
    output_phase = zeros(length(fbands), nSamples);
    for iB = 1:length(fbands)
        true_phase(iB,:) = filt_phase{iB}(nEnds);
        for iS = 1:nSamples
           this_echt = echt(csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
           this_phase = angle(this_echt);
           output_phase(iB,iS) = this_phase(end); % The last sample's phase
        end
    end
    phases = cell(1,2);
    phases{1} = true_phase;
    phases{2} = output_phase;
end
