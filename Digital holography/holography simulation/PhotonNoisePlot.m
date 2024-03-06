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


% sampling of the pulse : integration time for each point = T

T = 1e-9 ;
Fs = 1e3 ; % sampling frequency

%% generation of random trace 

% Fourier structure : 
N = 2^10 ;       % number of point for FFT
F = TF_t(N,Fs) ; % Fourier time structure


I = (P/Ephoton)*T + F.t*0 ; % unit = Number of photon
I = PoissonNoise(I) ; 
I_FFT = F.fourier(I) ; 

trapz(F.f,abs(I_FFT).^2)
trapz(F.t,abs(I).^2)
mean(abs(I).^2)*F.T

%% plot the signal as a function of time

figure(1); clf ;
subplot(211)
plot(F.t,I,'o-')

% psd in [unit]/sqrt(Hz)
PowerSpectDens = 2*abs( I_FFT( ((F.N)/2+1):end ) ).^2 ;
f_psd = F.f( ((F.N)/2+1):end );

subplot(212)
loglog( f_psd , PowerSpectDens ,'o-') ; hold on
loglog( f_psd, mean(PowerSpectDens) + 0*f_psd , 'red') ; 

















