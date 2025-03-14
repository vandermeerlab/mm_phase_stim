  %Testbed for comparing various phase estimation methods



%% Change folder and restrict data to epoch with no opto-stim
cd('data\M18\M18-2019-04-11_vStr_4p2_light_cells_TT7_min')
LoadExpKeys;
evs = LoadEvents([]);
cfg.fc = {ExpKeys.goodCSC};
csc = LoadCSC(cfg);

eval_csc = restrict(csc, iv(ExpKeys.PreRecord)); % Alternatively use ExpKeys.PostRecord
% try after downsampling the signal
cfg = []; cfg.decimateFactor = 16;
eval_csc = decimate_tsd(cfg, eval_csc); 
Fs = 1./median(diff(eval_csc.tvec));

%% Plot PSD and extract acausal phase
fig = figure('WindowState', 'maximized') ;
subplot(2,10,[1 2]);
wsize = floor(Fs);
[P, F] = pwelch(eval_csc.data, hanning(wsize), wsize/2, [], Fs);
plot(F, 10*log10(P));
xlim([0 120]); xlabel('Frequency (Hz)'); ylabel('Power'); title('PSD');

fbands = {[2 5], [6 10], [20 55], [55 95]};
true_phase = cell(length(fbands),1);
for iB = 1:length(fbands)
    cfg_filt.type = 'fdesign'; 
    cfg_filt.f  = fbands{iB};
    filt_lfp = FilterLFP(cfg_filt, eval_csc);
    true_phase{iB} = angle(hilbert(filt_lfp.data));
end

%% Add Path for the current method
% addpath('mm_phase_stim\code-matlab\phase_estimation\Rosenblum');
% This method has been implemented by slight modification of the source
% code (see helper functions at the bottom)

%% Run the current method on eval_data

nSamples = 1000;
win_length = 1.5; %in seconds
nEnds = randi(length(eval_csc.data), nSamples, 1);
nStarts = nearest_idx3(eval_csc.tvec(nEnds) - win_length, eval_csc.tvec);

for iB = 1:length(fbands)
    estimated_phase = zeros(1,nSamples);
    for iS = 1:nSamples
        this_sample = eval_csc.data(nStarts(iS):nEnds(iS));
        npt = length(this_sample);  % the number of points to be measured and processed
              % (put here your own value)
        fs = Fs;    % sampling frequency in Hz (put here your own value)
        dt=1/fs;      % sampling interval
        nu=mean(fbands{iB});        % main frequency of the beta-band peak (put your own value)
              % only rough estimation required

        s=zeros(1,npt);      % reserve space for the filtered input signal
        Dphase=s; Dampl=s;   % reserve space for the "device" amplitude and phase

        % parameters of the measurement "device"
        om0=5*nu;   om0_2=om0*om0;         % oscillator frequency
        alpha_a=80; gama=alpha_a/2;        % damping for the "amplitude device"     
        alpha_p=10; gamp=alpha_p/2;        % damping for the "phase device"
        factor=sqrt((om0_2-nu*nu)^2+(alpha_a*nu)^2);
        % precomputed coefficients for solving oscillators's equations
        % for the "amplitude device":
        [C1a,C2a,C3a,enuDela,ealDela,etaa]=ExSolCoefs(om0,dt,gama);  
        % for the "phase device":
        [C1p,C2p,C3p,enuDelp,ealDelp,etap]=ExSolCoefs(om0,dt,gamp);  
        % Initialization
        x=0; y=0; u=0; v=0;              
        %
        % Now comes the causal phase and amplitude estimation
        % We need three time-points to start estimation
        s(1)=this_sample(1); % here shall be your function that provides the
                        % new measurement; it is assumed that s(k) is 
                        % already causally bandpassed, if required
        s(2)=this_sample(2);
        for k=3:npt             % loop over the number of measurements
            s(k)=this_sample(k); 
            % the next two lines provide the amplitude for s(k);
            % comment them out for speed, if you do not need the amplitude
%             [x,y]=OneStep(x,y,gama,etaa,enuDela,ealDela,C1a,C2a,C3a,s(k-2),s(k-1),s(k)); 
%             z=y/nu;  Dampl(k)=factor*sqrt(z*z+x*x); % new amplitude value
            % the next two lines provide the phase for s(k);                         
            % comment them out for speed, if you need only the amplitude
            [u,v]=OneStep(u,v,gamp,etap,enuDelp,ealDelp,C1p,C2p,C3p,s(k-2),s(k-1),s(k)); 
            z=v/nu;  Dphase(k)=atan2(-z,u);            % new phase value
        end
       estimated_phase(iS) = Dphase(end); % The last sample's phase
    end
    
    subplot(2,10, [(2*iB)+1, (2*iB)+2])
    scatter(estimated_phase, true_phase{iB}(nEnds), 4);
    hold on;
    plot([-pi pi], [-pi pi], 'k'); 
    l = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
    axis tight; grid on; set(gca, 'XTick', -pi:pi/2:pi, 'YTick', -pi:pi/2:pi, 'XTickLabel', l, 'YTickLabel', l);

    title(sprintf("%d Hz - %d Hz", fbands{iB}(1), fbands{iB}(2))); xlabel('estimated phase'); ylabel('true phase');
end

%% TODO: Plot the distribution of phases with just hilbert transform and the causal method

%% Obtain Hilbert-transfrom phases
test_csc = restrict(csc, iv([ExpKeys.timeOnWheel, ExpKeys.timeOffWheel]));
% try after downsampling the signal
cfg = []; cfg.decimateFactor = 16;
Fs = 1./median(diff(test_csc.tvec));

stim_times = evs.t{strcmp(evs.label,ExpKeys.laser_on)};
% This sanity check is necessary because of M020
stim_times = stim_times(stim_times > ExpKeys.timeOnWheel);
ISIs = [100 diff(stim_times)']; %100 is used as an arbitrarily large number so that the first stim is always included
keep = ISIs > win_length;

test_csc = decimate_tsd(cfg, test_csc);
ht_phase = zeros(length(fbands), sum(keep));

for iB = 1:length(fbands)
    cfg_filt = [];
    cfg_filt.type = 'fdesign'; 
    cfg_filt.f  = fbands{iB};
    filt_lfp = FilterLFP(cfg_filt, test_csc);
    filt_phase = angle(hilbert(filt_lfp.data));
    ht_phase(iB, :) = filt_phase(nearest_idx3(stim_times(keep), test_csc.tvec));
end

%% Obtain phases thorugh causal method
causal_phase = zeros(length(fbands), sum(keep));
nEnds = nearest_idx3(stim_times, test_csc.tvec);
nStarts = nearest_idx3(stim_times - win_length, test_csc.tvec);
for iB = 1:length(fbands)
    for iS = 1:sum(keep)
       this_echt = echt(test_csc.data(nStarts(iS):nEnds(iS)), fbands{iB}(1), fbands{iB}(2), Fs);
       this_phase = angle(this_echt);
%        plot(this_phase);
       causal_phase(iB,iS) = this_phase(end); % The last sample's phase
    end
end

%% Plot the distributions
for iB = 1:length(fbands)
   subplot(2,10,(2*iB)+ 11)
   histogram(ht_phase(iB,:), 5, 'FaceColor', 'Cyan');
   title('HT phases')
   subplot(2,10,(2*iB)+ 12)
   histogram(causal_phase(iB,:), 5, 'FaceColor', 'Magenta');
   title('Causal Phases')
end

%% Put some text
subplot(2,10, [11 12])
text(0.1, 0.6, strcat(ExpKeys.subject, '_', ExpKeys.date), 'Interpreter', 'none', 'FontSize', 16)
text(0.1, 0.5 , strcat('Window length Used: ', num2str(win_length), ' sec'), 'Interpreter', 'none', 'FontSize', 16)
text(0.1, 0.4 , strcat('Trials left:  ', num2str(sum(keep))), 'Interpreter', 'none', 'FontSize', 16)
text(0.1, 0.3, 'Method Used: Non-Rez','Interpreter', 'none', 'FontSize', 16)
box off
grid off
axis off
WriteFig(fig, 'nonrez', 1)

%% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C1,C2,C3,eetadel,ealdel,eta]=ExSolCoefs(om0,dt,alpha) 
% alpha is the half of the damping: x''+2*alpha*x'+om0^2*x=input
eta2=om0^2-alpha^2; eta=sqrt(eta2);
eetadel=exp(1i*eta*dt); a=1/eetadel; ealdel=exp(alpha*dt);
I1=1i*(a-1)/eta; I2=(a*(1+1i*eta*dt)-1)/eta2;
I3=(a*(dt*eta*(2+1i*dt*eta)-2*1i)+2*1i)/eta2/eta;
C1=(I3-I2*dt)/2/dt^2/ealdel;   C2=I1-I3/dt^2;
C3=ealdel*(I2*dt+I3)/2/dt^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,xd]=OneStep(x,xd,alpha,eta,eetaDel,ealDel,C1,C2,C3,sprev,s,snew) 
A=x-1i*(xd+alpha*x)/eta;
A=A-1i*(C1*sprev+C2*s+C3*snew)/eta;
d=A*eetaDel;
y=real(d); yd=1i*eta*(d-conj(d))/2;
x=y/ealDel; xd=(yd-alpha*y)/ealDel;
end
