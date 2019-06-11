clearvars
addpath('..\functions');
addpath('..\');
%% generate 1D phase profile
Fe = 50e6; % sampling frequency
dt = 1/Fe;
N = 2^14;
t = (-N/2:N/2-1)*dt;
t0 = min(t) ; % zero for screening
f0 = 6e6;  % US frequency in Hz

Tmax = max(t) - min(t);

F = TF1D(N,Fe);

% structuring frequency 
Tjump = 1e-6 ;
N = floor(Tmax/Tjump);
Np = floor(Tjump*Fe) ;
nu = 5/(max(F.t)-min(F.t));
t_r = round(rand(1,N));

%% AO tagged photons

%% bascules de phase aléatoires
phi = pi + 0*F.t;
As  = 1 + 0*F.t;
for i=0:(N-1)
phi( (i*Np+1):((i+1)*Np)) = phi( (i*Np+1):((i+1)*Np))*t_r(i+1) ;
end
%phi = phi + 2*pi*f0*F.t ;


 %% frequency ramps
Trep = 10e-6;
Nrep = 10 ;
tau_c = 2*Nrep*Trep;   % camera integration time
PHI2 = 0*(0.9e6)/tau_c; % 0.185

% phi = 2*pi*f0*F.t + df*(F.t).^2;
 
% periodic ramp

% phi = 2*pi*f0*Trep*sawtooth(2*pi*t/Trep) ;
% phi = 2*pi*f0*Trep*sawtooth(2*pi*t/Trep) + 2*pi*df*(Trep*sawtooth(2*pi*t/Trep)).^2;
phi = 2*pi*0*t + 2*pi*PHI2*( Trep*( sawtooth(2*pi*t/Trep) + 1 )/2 ).^2 ;


figure(1); 
hold on
%plot(F.t/Trep, phi)
plot(F.t/Trep, phi )
xlabel('t( \mu s )')
title('\phi(t)')

%% modulation d'amplitude
As = (1+cos(2*pi*t/Trep))/2;
Bs = (1+cos(2*pi*t/Trep+pi))/2; 

%% command AO
Es = As.*exp(1i*( phi ) );
Es_fft = F.fourier(Es);

delay = (10e-3)/1540 ; % delay = z/v_us 
Er = (abs(F.t) < tau_c/2 ).*interp1(t,Es,F.t-delay,'linear',0);

% figure(2)
% subplot(211)
% plot(F.t*1e6,real(Es))
% hold on
% plot(F.t*1e6,real(Er))
% xlabel('t( \mu s )')
% hold off
% subplot(212)
% plot(F.f/f0,abs(Es_fft))
% xlim([0 2])

%% autocorrelation
[acor,lag] = xcorr(Er,Bs);
[~,I] = max(abs(acor));
lagDiff = lag(I);
timeDiff = lagDiff/Fe ;
tau = lag/Fe;

%% manual evaluation

T_fwhm = FWHM( abs(acor)/max(abs(acor)) , 1e6*tau );

figure(3);
hold on
plot( 1e6*tau,abs(acor)/max(abs(acor)) )
title(['xCorr maximum at ',num2str(1e6*timeDiff),'\mu s'])
xlabel('time(\mu s)')
ylabel('xcorr (a.u)')
xlim([-50 50])
legend(['fwhm = ',num2str(T_fwhm),'\mu s'])


%  hold on
%  Ith = exp(2*1i*pi*PHI2*Trep*(tau-Trep)).*tau.*sinc(2*PHI2*tau.*(tau-Trep)) ...
%      + exp(2*1i*pi*PHI2*Trep*tau).*(Trep-tau).*sinc(2*PHI2*tau.*(tau-Trep)) ;
% %  Ith(tau<0) = 0 ;
% %  Ith(tau>Trep) = 0;
%  plot(1e6*tau,abs(Ith)/max(abs(Ith)))
%  %xlim([0 1e6*Trep])
%  hold off



%%
