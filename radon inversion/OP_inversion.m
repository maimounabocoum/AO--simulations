
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% load experiemental data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
load('experiement images - JB - test\OP0deg-2016-02-01_13-11.mat');
c = 1540 ; % sound velocity in m/s
% update of the OP structure
% X : theta 
% Y : monotonic t variable

MyImage = OP(data(:,:,1),X,Y,Param,c); 
MyImage.Show_R();
MyImage.Fmax() % maximum frequency sampling = 1/dt
MyImage.Get_Lsample()

%% Nyist principle states the sampling of the object to reconstruct to be such that w > w_max/2 

%%% Fourier Tranform of Radon input image with respect to t
N = 2^11;
% creating a fourier parameter set using the measure sample size in m
FourierParam = TF_t(N,MyImage.Fmax());


% fourier transform of image
R = interp1(MyImage.t,MyImage.R,FourierParam.t,'linear',0);
F_R = FourierParam.fourier(R);

%% image show
figure;
subplot(221)
imagesc(MyImage.theta,FourierParam.t,R)
xlabel('\theta (°)')
ylabel('\omega (m^{-1})')
title('interpolation of Radon transform')
subplot(222)
imagesc(MyImage.theta,FourierParam.w(N/2:end),log(abs(F_R(N/2:end,:))))
xlabel('\theta (°)')
ylabel('\omega (m^{-1})')
title('Fourier Transform of Radon')

% representation in polar coordinates:

subplot(223)
[THETA, W] = meshgrid(MyImage.theta*pi/180,FourierParam.w(N/2:end));
[X,Y] = pol2cart(THETA, W);
surfc(X,Y,log(abs(F_R(N/2:end,:))))
view(0,90)
shading interp
xlabel('\omega\it_{x} (\itm^{-1})')
ylabel('\omega\it_{y} (\itm^{-1})')
title('Fourier Transform of Radon')


% filtered inverse fourier transform :
% defninition of the filter :

FILTER = abs(FourierParam.w)'*ones(1,length(MyImage.theta));


I = FourierParam.ifourier(F_R.*FILTER);


subplot(224)
imagesc(MyImage.theta,FourierParam.t,I)
xlabel('\theta (°)')
ylabel('\omega (m^{-1})')
title('identity op')





