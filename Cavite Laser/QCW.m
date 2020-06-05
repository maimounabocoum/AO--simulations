%% QCW laser amplifier
clearvars;
%% input laser
addpath('Q:\AO--commons\shared functions folder');

F = TF_t(1024,5e6);
taus_fwhm = 100e-6;
E0_s      = 0.1e-3;
Pulse0    = exp(-log(2)*(2*F.t/taus_fwhm).^6);
Pulse0    = E0_s*Pulse0/trapz(F.t,Pulse0);
 
figure(1)
plot(F.t*1e6,Pulse0)
xlabel('time(\mu s)')
ylabel('peak power (W)')

%%