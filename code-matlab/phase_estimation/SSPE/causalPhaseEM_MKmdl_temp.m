% causal phase estimates using the SP model and EM with fixed interval
% smoothing across windows of data
% primary benefit: assume a transitory burst of oscillatory activity in the
% range of your bandpass filter/ assume a peak shift in the data towards the edge
% of the band pass filter. These are problems unaddressed for instantaneous
% phase estimation right now
% Potential extension: a diffuse prior on x to estimate all frequencies? 

% Algorithm:
% after estimating reasonable initialization points we need to run EM on
% the data - at whatever rate makes it possible to run it again before
% getting the next window of data. So while it might be slow here, a C++
% implementation is likely going to be much faster meaning we have have
% small windows (up to the frequency limits of course)

% for the prototype all I need to do i run it on data already collected.
% And split it into whatever sized segments make sense. And then run three
% methods on it to show that the EM approach works best given a certain
% type of data. 

% 10/19/2020
% The final algorithm estimates parameters on an initial window and never
% updates parameters.

% all the pieces should be run together:
% 1. Initialization - Using the MK initialization approach
% 2. Use the EM to estimate parameters from the first window
% 3 Kalman filter with the latest parameter estiamtes
% Last edit: Ani Wodeyar 10/19/2020

% Modify code to check if omega estimates lie within the band-stop
% frequencies, and use that instead of the pre-defined omaega values.

function [omega, phase,phaseBounds,allX_full,circstd,cosineSum] = causalPhaseEM_MKmdl_temp(y,initParams)

    freqs = initParams.freqs;
    Fs = initParams.Fs;
    ampVec = initParams.ampVec;
    sigmaFreqs = initParams.sigmaFreqs;
    sigmaObs = initParams.sigmaObs;
    windowSize = initParams.window;
    freqBands = initParams.freqBands;

    if windowSize < Fs
        disp('The window size needs to be different. Setting it equal to sampling rate')
        windowSize = Fs;
    end


    ang_var2dev = @(v) sqrt(-2*log(v)); % note the difference in definition (ie not (1-v))

    data = y(1:windowSize);
    % first run to set up parameters
    [omega_est, ampEst, allQ, R, stateVec, stateCov] = fit_MKModel_multSines(data,freqs, Fs,ampVec, sigmaFreqs,sigmaObs);
%     lowFreqLoc = find((omega>lowFreqBand(1)) & (omega<lowFreqBand(2)),1);

%     if isempty(lowFreqLoc)
%         disp('Low freq band limits incorrect OR there is no low freq signal; retaining initial params')
%         omega = freqs;
%         ampEst = ampVec;
%         allQ = sigmaFreqs;
%         [~,lowFreqLoc] = min(abs(freqs-mean(lowFreqBand))); % pick frequency closest to middle of low frequency range
%     end
    
    omega = zeros(size(freqs));
    ampEst_temp = zeros(size(ampVec));
    allQ_temp = zeros(size(sigmaFreqs));
    for iF = 1:length(freqs)
        if omega_est(iF) >= freqBands{iF}(1) && omega_est(iF) <= freqBands{iF}(2)
            omega(iF) = omega_est(iF);
            ampEst_temp(iF) = amp_est(iF);
            allQ_temp(iF) = allQ_temp(iF);
        else
            disp_str = sprintf('For frequency band %.0f Hz - %.0f Hz, initial estimates are used\n', ...
                        freqBands{iF}(1), freqBands{iF}(2));
            fprintf(disp_str);
            omega(iF) = freqs(iF);
            ampEst_temp(iF) = ampVec(iF);
            allQ_temp(iF) = sigmaFreqs(iF);
        end
    end


    ampEst = ampEst_temp;
    allQ = allQ_temp;

    % for loop that runs through rest of the data reestimating parameters after
    % generating phase estimates for the whole period using past parameter ests
    % and the kalman filter
    [phi, Q, M] = genParametersSoulatMdl_sspp(omega, Fs, ampEst, allQ);
    y_rest = y(windowSize+1:end);
    wz_rest = length(y_rest);
    phase = zeros(length(freqs), wz_rest);
    phaseBounds = zeros(length(freqs), wz_rest,2);
    allX_full = zeros(length(freqs)*2, wz_rest);
    circstd = zeros(1,wz_rest);
    cosineSum = zeros(1,wz_rest);
    allP = zeros(length(freqs)*2,length(freqs)*2, wz_rest);

    x = stateVec(:,end);
    P = squeeze(stateCov(:,:,end));

    for i = 1:wz_rest
        [x_new,P_new] = oneStepKFupdate_sspp(x,y_rest(i),phi,M,Q,R,P);
        allX_full(:,i) = x_new;
        P_new = (P_new + P_new') /2; % forcing symmetry to kill off rounding errors
        allP(:,:,i) = P_new; 

        for iF = 1:length(freqs)
            % estimate phase
            phase(iF, i) = angle(x_new(iF*2-1) + 1i* x_new(iF*2));
            samples = mvnrnd(x_new(iF*2-1:iF*2), P_new(iF*2-1:iF*2,iF*2-1:iF*2),2000);
            sampleAngles = (angle(exp(1i*angle(samples(:,1) + 1i*samples(:,2)) - 1i*phase(iF,i)))); % removing mean
            tmpAngleStd = (wrapTo2Pi((prctile(sampleAngles,97.5) - prctile(sampleAngles,2.5)))/2); %[0,pi/2]
            phaseBounds(iF,i,:) = sort([(phase(iF,i)- tmpAngleStd), (phase(iF,i) + tmpAngleStd)]); % can have a range of [0,pi]
            circstd(iF,i) = ang_var2dev(abs(mean(exp(1i*sampleAngles))));
            cosineSum(iF,i) = mean(cos(2*(sampleAngles)));
        end

        % update state and state cov
        P = P_new;
        x = x_new;
    end
end




