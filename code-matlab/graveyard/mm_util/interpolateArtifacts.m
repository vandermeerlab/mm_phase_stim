function csc_out = interpolateArtifacts(method, stim_times, artifact_window,  csc_in)
    stim_start = nearest_idx3(stim_times, csc_in.tvec);
    stim_end = nearest_idx3(stim_times + artifact_window, csc_in.tvec);
    to_be_int = false(size(csc_in.tvec));
    for iS = 1:length(stim_start)
        this_idx = stim_start(iS):stim_end(iS);
        to_be_int(this_idx) = true;
    end
    midx = find(to_be_int);
    pidx = find(~to_be_int);
    csc_out = csc_in;
    for iC = 1:length(csc_out.label)
        mval = interp1(pidx, csc_out.data(iC,pidx), midx, method);
        csc_out.data(iC,midx) = mval;
    end
end