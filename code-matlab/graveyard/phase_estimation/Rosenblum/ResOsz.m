function ResOsz
% This code performs real-time measurement of phase and amplitude
% using a resonant oscillator as a "measurement" device
%
% For a description see:
%    Scientific Reports, M. Rosenblum, A. Pikovsky, A. Kuehn2, 
%    and J.L. Busch, Real-time estimation of phase and amplitude with
%    application to neural data (2021, submitted)
%    arXiv:2105.10404 [q-bio.NC]
%
%    For an example of application, see tremor data example
%    and Fig. 3 in the mentioned publication
%    
%    This version includes a frequency-adaptation algorithm
%    This version removes low-frequency trend (baseline fluctuations)

npt = 10000;  % the number of points to be measured and processed
              % (put here your own value)
fs = 1000;    % sampling frequency in Hz (put here your own value)
dt=1/fs;      % sampling interval
nu=4.5;       % rough estimate of the tremor frequency (put your own value)
om0=2*pi*nu;  % angular frequency
             
s=zeros(1,npt);      % reserve space for the filtered input signal
Dphase=s; Dampl=s;   % reserve space for the "device" amplitude and phase
sdemean=s;           % reserve space for the baseline-corrected signal

% parameters of the measurement "device"
alpha=0.3*om0; mu=500;  gam=alpha/2; % device parameters 
% parameters of adaptation algorithm
nperiods=1;   % number of previous periods for frequency correction
npt_period=round(2*pi/om0/dt);  % number of points per average period
npbuf=nperiods*npt_period;   % number of points for frequency correction buffer
update_factor=20;            % number of frequency updates per period
update_step=round(npt_period/update_factor);
% precomputed quantities for linear fit for frequency adaptation
tbuf=(1:npbuf)*dt; Sx=sum(tbuf); denom=npbuf*sum(tbuf.^2)-Sx*Sx; 
% precomputed coefficients for the device
[C1,C2,C3,enuDel,ealDel,eta]=ExSolCoefs(om0,dt,gam);  % for the oscillator
edelmu=exp(-dt/mu);                                   % for the integrator      
% Initialization
runav=-0.1;                        % initial guess for the dc-component
x=0; y=0; z=0; ypp=0; yp=0;            
updatepoint=3*npbuf;  
%
% Now comes the causal phase and amplitude estimation
% We need three time-points to start estimation
s(1)=GetNewMeasurement; % here shall be your function that provides the
                        % new measurement; it is assumed that s(k) is 
                        % already causally bandpassed, if required
s(2)=GetNewMeasurement;
for k=3:npt             % loop over the number of measurements
    s(k)=GetNewMeasurement; 
    sdemean(k)=s(k)-runav;  % baseline correction
    [x,y]=OneStep(x,y,gam,eta,enuDel,ealDel,C1,C2,C3,...
        sdemean(k-2),sdemean(k-1),sdemean(k)); 
    z=OneStepInt(z,edelmu,mu,dt,ypp,yp,y);
    ypp=yp;   yp=y;  v=mu*z*om0;
    Dphase(k)=atan2(v,y);            % new phase value
    Dampl(k)=alpha*sqrt(y*y+v*v);    % new amplitude value
    %
    if(k>updatepoint)
        buffer=unwrap(Dphase(k-npbuf+1:k))'; % buffer for frequency estimation
            % it contains phases of npbuf previous points 
            % frequency via the linear fit of the phases in the buffer: 
        om0=(npbuf*sum(tbuf.*buffer)-Sx*sum(buffer))/denom;  
            % recompute integration parameters for the new frequency:
        [C1,C2,C3,enuDel,ealDel,eta]=ExSolCoefs(om0,dt,gam); 
        runav=mean(s(k-npbuf+1:k));              % update running average
        updatepoint=updatepoint+update_step;     % point for the next update
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C1,C2,C3,eetadel,ealdel,eta]=ExSolCoefs(om0,dt,alpha) 
% alpha is the half of the dampling: x''+2*alpha*x''+om0^2*x=input
eta2=om0^2-alpha^2; eta=sqrt(eta2);
eetadel=exp(1i*eta*dt); a=1/eetadel; ealdel=exp(alpha*dt);
I1=1i*(a-1)/eta; I2=(a*(1+1i*eta*dt)-1)/eta2;
I3=(a*(dt*eta*(2+1i*dt*eta)-2*1i)+2*1i)/eta2/eta;
C1=(I3-I2*dt)/2/dt^2/ealdel;   C2=I1-I3/dt^2;
C3=ealdel*(I2*dt+I3)/2/dt^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z=OneStepInt(z,edelmu,mu,dt,ypp,yp,y)
dt2=dt*dt;
a=yp; b=(y-ypp)/dt/2;  c=(ypp-2*yp+y)/dt2/2;
d=-a+b*mu-2*c*mu^2;    C0=z+d; 
z=C0*edelmu-d+b*dt-2*c*mu*dt+c*dt2;
end

