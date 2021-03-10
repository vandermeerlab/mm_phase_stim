function csc_out = interpolateArtifacts(method, stim_times, stim_window,  csc_in)
    stim_start = stim_times;
    stim_end = stim_start + stim_window;
    art_idx = false(size(csc_in.tvec));
    % set indices for the  artifact samples
    for i = 1:length(stim_times)
        art_idx(csc_in.tvec > stim_start(i) & csc_in.tvec < stim_end(i)) = 1;
    end
    midx = find(art_idx);
    pidx = find(~art_idx);
    csc_out = csc_in;
    for i = 1:length(csc_in.label)
        mval = interp1(pidx, csc_out.data(i,pidx), midx, method);
        csc_in.data(i,midx) = mval;
    end
end