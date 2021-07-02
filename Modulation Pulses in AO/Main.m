%% addpath
addpath('Q:\AO--commons\shared functions folder')

%% load TF structure:

F = TF_t(2^12,20e6);

% 100e-6 long pulse

% enveloppe
Pulse = ones(1,F.N);
Pulse(abs(F.t)>= 50e-6) = 0;

% carrier
f_AO = 0.3e6 ;
w0 = 2*pi*f_AO;
Pulse = Pulse.*exp(1i*w0*F.t).*cos(2*pi*(0.03e6)*F.t);
%Pulse = Pulse.*exp(1i*w0*(F.t - (F.t).^2/(200e-6)));
PulseSpectra = F.fourier(Pulse.*exp(-1i*w0*F.t));

figure(1);
% subplot(121)
plot(1e6*F.t,abs(Pulse),'linewidth',5);
xlim([-60 60])
xlabel('time (\mu s)')
ylabel('a.u')
% subplot(122)
% plot(1e-6*(F.f+f_AO),abs(PulseSpectra),'linewidth',2);
% xlim(f_AO*1e-6+[-1 1])
% xlabel('frequency (MHz)')
% ylabel('a.u')

%% properties of AP modulator
Cs = 4200 ;         % sound velocity in m/s in Te02
F = 80e6 ;          % center frequency for Bragg
lambda = 700e-9 ;   % illumination laser wavelength
n = 1 ;             % optical index of the material

Theta_Bragg = (lambda*F/(2*Cs*n));

Theta_spread = 2*Theta_Bragg

%%









