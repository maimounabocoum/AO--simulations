%% field fourier transform
clearvars

lambda0 = 800e-9 ;
c = 3e8 ;
w0 = 2*pi*(c/lambda0) ;
wus = 2*pi*6e6;
L = 20e-6;
tau = 0.5e-12;
x = linspace(-L/2,L/2,100);

% temporal profile 
N = 2^13 ;
dt = 1e-15 ;
Fmax = 1/dt;
F1 = TF1D(N,Fmax) ;
F2 = TF1D(N,500e6) ;

PulseLOA = exp(-(F1.t).^2/(5e-15)^2);
SpectraLOA = F1.fourier(PulseLOA);


SpectraIL = exp(-(F2.f).^2/(3e6)^2)+ exp(-(F2.f-6e6).^2/(1e6)^2);
PulseIL = F2.ifourier(SpectraIL);
% PHI_t = exp(-(F2.t).^2/(2e-6)^2).*cos(wus*F2.t);
% PulseIL = exp(1i*w0*F2.t+PHI_t);
% SpectraIL = F1.fourier(PulseLOA);


%  figure; plot(F1.t*1e15,PulseLOA)
%  figure; plot(F1.f,abs(SpectraLOA))
 figure; plot(F2.t*1e6,real(PulseIL))
 xlabel('\mu s')
 figure; plot(F2.f*1e-6,abs(SpectraIL))
 xlabel('MHz')
 
fwhm_t = FWHM(PulseLOA,F1.t*1e15)
fwhmt_f = FWHM(abs(SpectraLOA),F1.f)


(fwhm_t*fwhmt_f)/8.8326e+14

fwhm_t = FWHM(abs(PulseIL),F2.t*1e6)
coherenceLength = (3e8)*(fwhm_t*1e-6)
fwhmt_f = FWHM(SpectraIL,F2.f*1e-6)
