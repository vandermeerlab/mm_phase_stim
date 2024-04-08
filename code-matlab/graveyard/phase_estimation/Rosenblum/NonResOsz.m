function NonResOsz
% This code performs real-time measurement of phase and amplitude
% using a non-resonant oscillator as a "measurement" device
%
% For a description see:
%    Scientific Reports, M. Rosenblum, A. Pikovsky, A. Kuehn2, 
%    and J.L. Busch, Real-time estimation of phase and amplitude with
%    application to neural data (2021, submitted)
%    arXiv:2105.10404 [q-bio.NC]
%
%    For an example of application, see beta-band brain activity example
%    and Fig. 4 in the mentioned publication
%    

npt = 10000;  % the number of points to be measured and processed
              % (put here your own value)
fs = 1000;    % sampling frequency in Hz (put here your own value)
dt=1/fs;      % sampling interval
nu=18;        % main frequency of the beta-band peak (put your own value)
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
s(1)=GetNewMeasurement; % here shall be your function that provides the
                        % new measurement; it is assumed that s(k) is 
                        % already causally bandpassed, if required
s(2)=GetNewMeasurement;
for k=3:npt             % loop over the number of measurements
    s(k)=GetNewMeasurement; 
    % the next two lines provide the amplitude for s(k);
    % comment them out for speed, if you do not need the amplitude
    [x,y]=OneStep(x,y,gama,etaa,enuDela,ealDela,C1a,C2a,C3a,s(k-2),s(k-1),s(k)); 
    z=y/nu;  Dampl(k)=factor*sqrt(z*z+x*x); % new amplitude value
    % the next two lines provide the phase for s(k);                         
    % comment them out for speed, if you need only the amplitude
    [u,v]=OneStep(u,v,gamp,etap,enuDelp,ealDelp,C1p,C2p,C3p,s(k-2),s(k-1),s(k)); 
    z=v/nu;  Dphase(k)=atan2(-z,u);            % new phase value
end
end
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
