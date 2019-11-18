%% resolving power test
clearvars
% spatial profile 
N = 2^13 ;
dx = 10e-6 ;
Fmax = 1/dx;
F1 = TF1D(N,Fmax) ;


%% two exponential
P = exp(-(F1.t-2e-3).^2/(1e-3)^2) + exp(-(F1.t+2e-3).^2/(1e-3)^2);

%% fft
Spectre = F1.fourier(P);

figure;
plot(F1.t*1e3,P)
xlabel('x(mm)')

figure;
plot(F1.f*1e-3,abs(Spectre))
xlabel('f(mm^{-1})')
