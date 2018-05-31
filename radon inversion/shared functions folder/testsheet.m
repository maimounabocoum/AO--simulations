%% 
clearvars

x = (-19:0.2:19)*1e-3 ;
z = 0:(150e-6):(30e-3);

sigma = 5e-3;
 I = 0*exp(-x.^2/sigma^2) ;
%  I = I + 1;
 I = I + 0.1*rand(1,length(I));


Fmax = 1/(x(2)-x(1));
N = 2^7;
MyImage = TF_t(N,Fmax);

Iinterp = interp1(x,I,MyImage.t,'linear',0);
Ifft = MyImage.fourier(Iinterp);

figure(1);
hold on
plot(MyImage.t*1e3,Iinterp,'o-')

figure(2);
hold on
semilogy(MyImage.f/1000,abs(Ifft),'-o')