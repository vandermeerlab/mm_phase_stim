% this is a demo of the causalPhaseEM_MK
% first generate some data
% creating a pure sine with some amplitude noise and added pink noise at
% sampling frequency of 1000 Hz for 10 seconds

Fs = 1000;
time = 10; 
Vlo = (10).*cos(2*pi*(6).*[1/Fs:1/Fs:time]); % random movement in amplitude, creates a little bump in PSD
V1 = (2.5).*cos(2*pi*(50).*[1/Fs:1/Fs:time] + pi/3);
[pn] = make_pink_noise(1,1e4,1/Fs);
pn = 10*pn;
data = Vlo + V1;
truePhase =  wrapToPi((2*pi*(6).*[1/Fs:1/Fs:time]))';

%%

% following helps to initialize more reasonable estimates for the scaling
% parameter 'a' and the variances for the state vector
% [init_f, init_a,init_sigma,R0] = initializeParams(data,1,1,2000); % s

% setting up initial parameters to start the causalPhaseEM code
initParams.freqs = [4,45, 80]; % initialization using above not great at identifying starting freq
initParams.Fs = 1000;
initParams.ampVec = [.99,.99, 0.99]; % in a pinch this can be initialized to 0.99 to start
initParams.sigmaFreqs = [10,1, 0.1]; % its important to use a value of this that is in the same ballpark scale
initParams.sigmaObs = 1;
initParams.window = 8000;
initParams.lowFreqBand = [4,8];

numFreqs = length(initParams.freqs);
[phase,phaseBounds, fullX] = causalPhaseEM_MKmdl_temp(data, initParams);

figure
for i = 1:numFreqs
    subplot(numFreqs,1,i)
    plot(phase(i,:))
end

% TODO: Fix the code below this
phase = reshape(phase', size(phase,1) * size(phase,2),1);
phaseBounds = reshape(permute(phaseBounds,[2,1,3]), size(phaseBounds,1) * size(phaseBounds,2),size(phaseBounds,3));
fullX = reshape(permute(fullX,[2,1,3]), size(fullX,1) * size(fullX,2),size(fullX,3));

figure
err_Spcausal = angle(mean(exp(1i*(phase(2001:end) - truePhase(2001:end)))))*(180/pi);
imagesc( 1:10000,-3:3, (abs(fullX(:,1) + 1i*fullX(:,2)))', 'AlphaData', .5)
colormap(summer)
hold on; plot(1:10000, phase, 'Linewidth', 2, 'Color', 'b')
plot(squeeze(phaseBounds(:,1)), 'Linewidth', 2, 'color','r')
hold on
plot(squeeze(phaseBounds(:,2)), 'Linewidth', 2,'color','r')
set(gca,'Fontsize', 16)
xlabel('Time (in samples)')
ylabel('Phase')
h = colorbar;
ylabel(h,'Amplitude')
title(['Causal Phase Estimate with SP-EM, Error = ',num2str(err_Spcausal)])
xlim([0,10000])
