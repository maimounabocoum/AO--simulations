%% clasical modelization of noise:
clearvars ;
param = [1:10] ;

for loop = 1:length(param)
P           = param(loop) ;       % power in W
lambda0     = 800e-9;      % m
c           = 3e8;         % m/s
nu          = c/lambda0;   % Hz
h           = 6.626*1e-34; % J.s
Ephoton     = h*nu ;       % J
eta         = 1 ;

Fs = 1e6 ; % sampling frequency
T = 1/Fs ; % (condition for Df T = 1)

% generation of random trace 

% Fourier structure : 
N = 2^12 ;       % number of point for FFT
F = TF_t(N,Fs) ; % Fourier time structure

for i = 1:50
I = (P/Ephoton)*T + F.t*0 ; % unit = Number of photon
I_FFT = F.fourier(I) ; 

eta = 1e6 ;
Lindewidth = 2*eta./(F.f.^2 + (eta)^2).*exp(1i*2*pi*rand(1,F.N));
I = abs(F.ifourier(Lindewidth)) ; 
I = (P/Ephoton)*T*I/mean(I) ; % renormalize on average power 



% renormalize
trapz(F.f,abs(I_FFT).^2)
trapz(F.t,abs(I).^2)

I = PoissonNoise(I) ; 
I = I*Ephoton/T ;

%% get the mean value ; 

mesure(i) = trapz(F.t,I) ;

end

mu(loop) = mean(mesure) ;
sigma(loop) = var(mesure,1,'all') ;

end

figure(4) ; hold on
plot(param,mu)
%% plot the signal as a function of time

figure(2); clf;
subplot(211)
plot(F.t*1e3,I*1e3,'o-')
% ylim([0 2])
xlabel('time(ms)')
ylabel('power(mW)')

% psd in [unit]/sqrt(Hz)

[f_psd,PowerSpectDens] = F.psd(I-mean(I)) ;
subplot(212)
loglog( f_psd(2:end)*1E-3 , 10*log10(sqrt(PowerSpectDens(2:end))) ,'o-') ; hold on
loglog( f_psd(2:end)*1E-3, 10*log10(sqrt(2*Ephoton*P)) + 0*f_psd(2:end) , 'red') ; 
ylabel('dB [$W/\sqrt{Hz}]$',Interpreter='latex')
xlabel('frequency (kHz)')
















