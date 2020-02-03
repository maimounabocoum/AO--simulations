%% simulate reconstruction with periodic fringes autocorrelation
clearvars
addpath('..\functions');
addpath('..\..\AO--commons\shared functions folder');

%% definition of pression profile
Fz = 1/(10e-6);
Fs = 50e6;
MyF = TF2D(2^10,2^10,Fs,Fz);
omega_us =  2*pi*6e6;
omega_mod =  2*pi*9*50e3;
c    = 1540 ;
K_us    =  omega_us/c ;
K_mod    =  omega_mod/c ;
[T,Z] = meshgrid(MyF.x,MyF.z);

phi = 0.1 ;

P = (sin(omega_mod*T - K_mod*Z)>0).*...
    sin(omega_us*T - K_us*Z).*...
    sign(sin(omega_mod*T - K_mod*Z)) ;
                
%P = sin(omega_mod*T - K_mod*Z).*sin(omega_us*T - K_us*Z);
%P = sin(omega_mod*T - K_mod*Z);
%P = sin(omega_us*T - K_us*Z);
P = exp(1i*2*pi*phi*P);
P_fft = MyF.fourier(P);

figure(1)
subplot(121)
imagesc(real(P))
subplot(122)
imagesc(abs(P_fft)) ; colorbar
xlabel('f_t')
ylabel('f_z')
caxis([0 3e-8])
axis([610 670 450 490])

