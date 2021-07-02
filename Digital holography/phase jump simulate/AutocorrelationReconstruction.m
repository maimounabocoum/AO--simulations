%% simulate reconstruction with periodic fringes autocorrelation
clearvars
 addpath('D:\AO--commons\shared functions folder')
 addpath('sequences');
 addpath('subfunctions');
 addpath('C:\Program Files (x86)\Gage\CompuScope\CompuScope MATLAB SDK\CsMl')
 addpath('D:\AO--commons\read and write files')



%% generate 1D phase profile
N = 2^10;
S = TF2D(N,N,10000,10000);

[X,Y] = meshgrid(S.x,S.y);
Io = exp(-X.^2/(0.01)^2-Y.^2/(0.01)^2);
I_oFFT = S.fourier(Io);


figure(1)
imagesc(abs(Io))
xlabel('x (mm)')
ylabel('z (mm)')

figure(2)
imagesc(S.fx,S.fy,abs(I_oFFT))
colorbar
axis([-50 50 -50 50])
xlabel('\nu_x (mm^{-1})')
ylabel('\nu_z (mm^{-1})')

M   = 2^9;
nx  = [-10:-1,1:10];
ny  = 1:20;
nu0y = 10; % sampling frequency mm^-1
nu0x = 10; % sampling frequency mm^-1
phase = [0,0.25,0.5,0.75];
[NBX,NBY] = meshgrid(nx,ny) ;
R = TF2D(M,(M-1)*nu0x,(M-1)*nu0y);






I_FFT = zeros(M,M);

for i = 1:length(NBX(:))

SIG = 0;

% ideal projection
Mnm = exp(1i*(2*pi*nu0x*NBX(i)*(X-10) + 2*pi*nu0y*NBY(i)*(Y-15) ));
SIG = SIG + trapz(S.x,trapz(S.y,Io.*Mnm) );

% 4-phase projection
% for j = 1:length(phase)
% Mnm = ( 1 + cos( 2*pi*nu0x*NBX(i)*X + 2*pi*nu0y*NBY(i)*Y + 2*pi*phase(j) ) )/2 ;
% sum(sum(Io.*Mnm))
% 
% SIG = SIG + sum(sum(Io.*Mnm))*exp(1i*2*pi*phase(j)) ;
% %  figure(1);
% % imagesc(Io.*Mnm)
% %  drawnow
%     end

 I_FFT(M/2+1 + NBY(i),M/2+1 + NBX(i))  =  SIG ;
 I_FFT(M/2+1 - NBY(i),M/2+1 - NBX(i))  =  conj(SIG) ;
 

end

Ioriginal = R.ifourier(I_FFT);


figure(3);
imagesc(Ioriginal)
figure(4);
imagesc(R.fx,R.fy,abs(I_FFT))
axis([-50 50 -50 50])
colorbar



