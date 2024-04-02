addpath('..\..\..\AO--commons\shared functions folder');
addpath('..\..\..\AO--commons\common subfunctions');

%%
clearvars :
%% generate poissoninan statistics on laser using "real" parameters:

P           = 1e-3 ;       % power in W
lambda0     = 800e-9;      % m
c           = 3e8;         % m/s
nu          = c/lambda0;   % Hz
h           = 6.626*1e-34; % J.s
Ephoton     = h*nu ;       % J
eta         = 1 ;


% sampling of the pulse : integration time for each point = T


Fs = 1e6 ; % sampling frequency
T = 1/Fs ; % (condition for Df T = 1)

%% generation of random trace 

% Fourier structure : 
N = 2^10 ;       % number of point for FFT
F = TF_t(N,Fs) ; % Fourier time structure


I = (P/Ephoton)*T + F.t*0 ; % unit = Number of photon
I = PoissonNoise(I) ; 
I = I*Ephoton/T ;

% measurement of the variance od power
var(I,1)
sigmaP2 = (Ephoton/T)*P
mean((I-mean(I)).^2)
trapz( F.t/T , ( I-mean(I) ).^2 )/F.N

I_FFT = F.fourier(I-mean(I)) ; 

% trapz(F.f,abs(I_FFT).^2)
% trapz(F.t,abs(I).^2)
% mean(abs(I).^2)*F.T

%% plot the signal as a function of time

figure(1); clf ;
subplot(211)
plot(F.t*1e3,I*1e3,'o-')
ylim([0 2])
xlabel('time(ms)')
ylabel('power(mW)')

% psd in [unit]/sqrt(Hz)
PowerSpectDens = 2*abs( I_FFT( ((F.N)/2+1):end ) ).^2 ;
f_psd = F.f( ((F.N)/2+1):end );

subplot(212) ;
loglog( f_psd , PowerSpectDens ,'o-') ; hold on
loglog( f_psd, 2*Ephoton*P + 0*f_psd , 'red') ; 
ylabel('W/sqrt(Hz)')

















