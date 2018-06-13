%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017

N = 2^11 ;
Fmax_x = 1/(0.1e-3);
Fmax_z = 1/(0.1e-3);



MyTest = TF2D(N,Fmax_x,Fmax_z) ;

Fs = 20/( (N-1)*MyTest.dx );
theta = pi/5 ;

[X,Z] = meshgrid(MyTest.x,MyTest.z);
[FX, FZ] = meshgrid(MyTest.fx,MyTest.fz);
[X_phantom,Z_phantom] = meshgrid(x_phantom,z_phantom) ;
I = interp2(X_phantom,Z_phantom,MyTansmission,X,Z,'linear',0);



% fourier transform :
I_ktkx = MyTest.fourier(I) ;
Hsup = ( FZ.^2 + FX.^2 >= (Fs)^2 );
I_ktkx = I_ktkx.*Hsup ;

I = MyTest.ifourier(I_ktkx);

figure;
imagesc(MyTest.x*1e3,MyTest.z*1e3,I)

figure;
subplot(211) ; imagesc(MyTest.fx,MyTest.z*1e3,abs(I_ktkx)) ;
subplot(212) ; imagesc(MyTest.fx,MyTest.z*1e3,angle(I_ktkx)) ;

s_ctks = I_kxkz(:,(N/2+1)+ floor(Fs/MyTest.dfx));

s_ktks = MyTest.ifourierz(s_ctks ) ;

figure
imagesc(MyTest.x,MyTest.z,abs(I))
% axis([-100 100 -100 100])

xi = [Fs*sin(theta)+ Fs*cos(theta),-Fs*sin(theta)- Fs*cos(theta)] ;
zi = [Fs*cos(theta)- Fs*sin(theta),-Fs*cos(theta)+ Fs*sin(theta)] ;
% 
% improfile(abs(I_kxkz),xi,zi)











