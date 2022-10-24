function ModelDataAllAlg
% Model data, all three algorithms, the code reproduces Fig 1.
% Notice: in this example all frequencies are angular frequencies

[s,dt,npt]=GenerateMultiCompDataFrMod;  % generate test data
time=(1:npt)*dt; HTampl=abs(hilbert(s)); HTphase=angle(hilbert(s));
figure(1); 

%%%%% panel (a) the signal and its Hilbert amplitude  %%%%%%%%%%%
subplot(4,1,1); plot(time,s,time,HTampl); 
ylabel('$s,a_H$','Interpreter','Latex','FontSize',20);
set(gca,'Xticklabel',[]); ylim([-2,3]); 
title('phases');

%%%%% panel (b): phase locking approach  %%%%%%%%%%%%%%%%%    
LBphase=LockingBasedPhase(s,dt,npt);  
phidif=PhaseDifference(HTphase,LBphase);
subplot(4,1,2); plot(time,phidif); 
ylabel('$\varphi_H-\varphi_L$','Interpreter','Latex','FontSize',20);
set(gca,'Xticklabel',[]); ylim([-pi pi]);

%%%%% panel (c): nonresonant oscillator     %%%%%%%%%%%%%%%%%    
[Dphase,NRDampl]=NonResonantOscillator(s,dt,npt);
phidif=PhaseDifference(HTphase,Dphase);
subplot(4,1,3); plot(time,phidif); 
ylabel('$\varphi_H-\varphi_N$','Interpreter','Latex','FontSize',20);
set(gca,'Xticklabel',[]); ylim([-pi pi]);

%%%%% panel (d): resonant oscillator     %%%%%%%%%%%%%%%%%    
[Dphase,RDampl]=ResonantOscillator(s,dt,npt);
phidif=PhaseDifference(HTphase,Dphase);
subplot(4,1,4); plot(time,phidif);  ylim([-pi pi]);
ylabel('$\varphi_H-\varphi_R$','Interpreter','Latex','FontSize',20);
xlabel('time');

%%%%%%%%% amplitudes  %%%%%%%%%%%%
figure(2); 
subplot(2,1,1); plot(time,s,time,NRDampl); 
ylabel('$s,a_N$','Interpreter','Latex','FontSize',20);
set(gca,'Xticklabel',[]); title('amplitudes'); ylim([-2,4]);
subplot(2,1,2); plot(time,s,time,RDampl); 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dphase,Dampl]=ResonantOscillator(s,dt,npt)
Dphase=zeros(1,npt);  % reserve space for the "device" phase
Dampl=Dphase;         % reserve space for the "device" amplitude
om0=1.1;     % initial guess for frequency (10% large than the true value)
% parameters of the measurement "device"
alpha=0.3*om0; mu=500;  gam=alpha/2; 
% precomputed coefficients for the device
[C1,C2,C3,enuDel,ealDel,eta]=ExSolCoefs(om0,dt,gam);  % for the oscillator
edelmu=exp(-dt/mu);                                   % for the integrator      
% parameters of adaptation algorithm
update_factor=5;            % frequency correction 20 times per period
npbuf=round(2*pi/om0/dt);  % number of points for frequency estimation
updatepoint=2*npbuf;
update_step=round(npbuf/update_factor);
% precomputed quantities for linear fit for frequency adaptation
tbuf=(1:npbuf)*dt;  Sx=sum(tbuf); denom=npbuf*sum(tbuf.^2)-Sx*Sx; 
% Initialization
x=0; y=0; z=0; ypp=0; yp=0;            
for k=3:npt    % starting the "real-time" estimation
    [x,y]=OneStep(x,y,gam,eta,enuDel,ealDel,C1,C2,C3,s(k-2),s(k-1),s(k)); 
    z=OneStepInt(z,edelmu,mu,dt,ypp,yp,y);
    ypp=yp;   yp=y;  v=mu*z*om0;
    Dphase(k)=atan2(v,y);            % new phase value
    Dampl(k)=alpha*sqrt(y*y+v*v);    % new amplitude value
    if(k>updatepoint)    % frequency adaptation
        buffer=unwrap(Dphase(k+1-npbuf:k));     % buffer for frequency estimation
        om0=(npbuf*sum(tbuf.*buffer)-Sx*sum(buffer))/denom;  % new frequency estimate
        [C1,C2,C3,enuDel,ealDel,eta]=ExSolCoefs(om0,dt,gam); 
        updatepoint=updatepoint+update_step; % point for the next frequency correction
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dphase,Dampl]=NonResonantOscillator(s,dt,npt)
Dphase=zeros(1,npt);  % reserve space for the "device" phase
Dampl=Dphase;         % reserve space for the "device" amplitude
nu=1.1;     % initial guess for frequency (10% large than the true value)
% parameters of the non-resonant "device"
om0=5*nu;  om0_2=om0*om0;         % oscillator frequency
alpha_a=6; gama=alpha_a/2;       % damping for the "amplitude device"     
alpha_p=0.2;  gamp=alpha_p/2;       % damping for the "phase device"
factor=sqrt((om0_2-nu*nu)^2+(alpha_a*nu)^2);
% precomputed coefficients for oscillators
[C1a,C2a,C3a,enuDela,ealDela,etaa]=ExSolCoefs(om0,dt,gama);  % for the "amplitude device"
[C1p,C2p,C3p,enuDelp,ealDelp,etap]=ExSolCoefs(om0,dt,gamp);  % for the "phase device"
% parameters of adaptation algorithm
update_factor=5;            % frequency correction 20 times per period
npbuf=round(2*pi/nu/dt);  % number of points for frequency estimation
updatepoint=2*npbuf;
update_step=round(npbuf/update_factor);
% precomputed quantities for linear fit for frequency adaptation
tbuf=(1:npbuf)*dt;  Sx=sum(tbuf); denom=npbuf*sum(tbuf.^2)-Sx*Sx; 
% Initialization
x=0; y=0; u=0; v=0;          
for k=3:npt    % starting the "real-time" estimation
    % the next two lines provide the amplitude for s(k);
    % comment them out for speed, if you do not need the amplitude
     [x,y]=OneStep(x,y,gama,etaa,enuDela,ealDela,C1a,C2a,C3a,s(k-2),s(k-1),s(k)); 
    z=y/nu;  Dampl(k)=factor*sqrt(z*z+x*x); % new amplitude value
    % the next two lines provide the phase for s(k);                         
    % comment them out for speed, if you need only the amplitude
    [u,v]=OneStep(u,v,gamp,etap,enuDelp,ealDelp,C1p,C2p,C3p,s(k-2),s(k-1),s(k)); 
    z=v/nu;  Dphase(k)=atan2(-z,u);            % new phase value
    if(k>updatepoint)    % frequency adaptation
        buffer=unwrap(Dphase(k+1-npbuf:k));     % buffer for frequency estimation
        nu=(npbuf*sum(tbuf.*buffer)-Sx*sum(buffer))/denom;  % new frequency estimate
        updatepoint=updatepoint+update_step; % point for the next frequency correction
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LBphase=LockingBasedPhase(s,dt,npt)
LBphase=zeros(1,npt);             % reserve space for the device phase
% parameters of the measurement "device"
epsilon=0.8; K=1;
nRKsteps=5;                       % number of Runge-Kutta steps per dt
% parameters of adaptation algorithm
nuest=1.1;     % initial guess for frequency (10% large than the true value)
update_factor=20;            % frequency correction 20 times per period
npbuf=round(2*pi/nuest/dt);  % number of points for frequency estimation
updatepoint=2*npbuf;
update_step=round(npbuf/update_factor);
% precomputed quantities for linear fit for frequency adaptation
tbuf=(1:npbuf)*dt;  Sx=sum(tbuf); denom=npbuf*sum(tbuf.^2)-Sx*Sx; 
for k=3:npt    % starting the "real-time" estimation
    LBphase(k)=rk(LBphase(k-1),dt,nRKsteps,nuest,epsilon,s(k),s(k-1),s(k));  
    if(k>updatepoint)    % frequency adaptation
        buffer=LBphase(k+1-npbuf:k);     % buffer for frequency estimation
        nu_quasi=(npbuf*sum(tbuf.*buffer)-Sx*sum(buffer))/denom;  %slope
        nuest=nuest+K*(nu_quasi-nuest);
        updatepoint=updatepoint+update_step; % point for the next frequency correction
    end
end
end
% %%%%%%%%%%%%%%%%%%%%%%%%  RK-solver, adapted  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi=rk(phi,dt,nsteps,omega,epsilon,sp,s,snew)    % x is time,
h=dt/nsteps;
b=(snew-sp)/dt/2; c=(sp-2*s+snew)/dt^2;
hh=h/2.; h6=h/6.;                        % y is phi (dependent variable),
t=0;                                     % initial time
for i=1:nsteps                           % h is the step Delta,
   th=t+hh;                              % nsteps is the number of steps
   dy0=omega - epsilon*(s+b*t+c*t*t)*sin(phi);
   phit=phi+hh*dy0;
   dyt=omega - epsilon*(s+b*th+c*th*th)*sin(phit);
   phit=phi+hh*dyt;
   dym=omega - epsilon*(s+b*th+c*th*th)*sin(phit);
   phit=phi+h*dym;
   dym=dym+dyt;
   t=t+h;
   dyt=omega - epsilon*(s+b*t+c*t*t)*sin(phit);
   phi=phi+h6*(dy0+2*dym+dyt);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,xd]=OneStep(x,xd,alpha,eta,eetaDel,ealDel,C1,C2,C3,sprev,s,snew) 
A=x-1i*(xd+alpha*x)/eta;
A=A-1i*(C1*sprev+C2*s+C3*snew)/eta;
d=A*eetaDel;
y=real(d); yd=1i*eta*(d-conj(d))/2;
x=y/ealDel; xd=(yd-alpha*y)/ealDel;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C1,C2,C3,eetadel,ealdel,eta]=ExSolCoefs(om0,del,alpha) 
% del is the time step dt
% alpha is the half of the dampling: x''+2*alpha*x'...
eta2=om0^2-alpha^2; eta=sqrt(eta2);
eetadel=exp(1i*eta*del); a=1/eetadel; ealdel=exp(alpha*del);
I1=1i*(a-1)/eta; I2=(a*(1+1i*eta*del)-1)/eta2;
I3=(a*(del*eta*(2+1i*del*eta)-2*1i)+2*1i)/eta2/eta;
C1=(I3-I2*del)/2/del^2/ealdel;   C2=I1-I3/del^2;
C3=ealdel*(I2*del+I3)/2/del^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z=OneStepInt(z,edelmu,mu,dt,ypp,yp,y)
dt2=dt*dt;
a=yp; b=(y-ypp)/dt/2;  c=(ypp-2*yp+y)/dt2/2;
d=-a+b*mu-2*c*mu^2;    C0=z+d; 
z=C0*edelmu-d+b*dt-2*c*mu*dt+c*dt2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phidif=PhaseDifference(phi1,phi2)
cosphidif=cos(phi1).*cos(phi2)+sin(phi1).*sin(phi2);
sinphidif=sin(phi1).*cos(phi2)-cos(phi1).*sin(phi2);
phidif=atan2(sinphidif,cosphidif);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,dt,npt]=GenerateMultiCompDataFrMod
npt=100000;
dt=1/100; t=(1:npt)*dt;
omega1=sqrt(2)/30;
omega2=sqrt(5)/60;
amp=1+0.95*cos(omega1*t);
p=t+5*sin(omega2*t);
s=amp.*(cos(p)+0.2*cos(2*p+pi/6)+0.1*cos(3*p+pi/3));
end

