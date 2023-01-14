function phases = wrapper_echt(csc, wsz, fbands, Fs, nSamples, seed)
% This is a wrapper function that returns true and estimated phase values
% as obtained by the echt method
    rng(seed);
    min_start = ceil(wsz*Fs);
    nEnds = randi(length(csc.data) - min_start, nSamples, 1) + min_start;
    nStarts = nearest_idx3(csc.tvec(nEnds) - wsz, csc.tvec);
    true_phase = zeros(length(fbands), nSamples);
    output_phase = zeros(length(fbands), nSamples);
    for iB = 1:length(fbands)
        cfg_filt = [];
        cfg_filt.type = 'fdesign'; 
        cfg_filt.f  = fbands{iB};
        filt_lfp = FilterLFP(cfg_filt, csc);
        this_phase = angle(hilbert(filt_lfp.data));
        true_phase(iB,:) = this_phase(nEnds);
        clear filt_lfp this_phase;
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
