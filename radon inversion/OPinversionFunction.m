function Ireconstruct = OPinversionFunction( Angle_Rad , Z_m , I , SampleRate_Hz, c)
% function copied from OP_inversion script 
% maimouna bocoum - 03/07/2017
% clearvars;
addpath('functions')
MyImage = OP(I , Angle_Rad , Z_m , SampleRate_Hz , c );

N       = 2^12;
Lobject = 3e-3;
Fc      = 1/Lobject;  % Lobject is the size of the object to detect. Using simple model (sinc function)
                      % we set it to kc = 100/Lobject 
MyImage = MyImage.InitializeFourier(N,10*Fc);
MyImage.F_R = MyImage.fourier(MyImage.R) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtered inverse fourier transform :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter options : 'ram-lak' (default) , 'cosine', 'hamming' , 'hann'
FilterType = 'ram-lak';%'ram-lak' 

filt = FilterRadon(MyImage.f, MyImage.N ,FilterType , Fc);
filt = filt(:);
FILTER = filt*ones(1,length(MyImage.theta));
I = MyImage.ifourier(MyImage.F_R.*FILTER);
 
% extract image back to initial size :
[I,z_out] = ReduceDataSize( I,'y',MyImage.t,[0 50e-3]);%MyImage.L max(Z_m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% reconstruction BOX initialization (retroprojection):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %xsonde = (1:size(DelayLAWS,1))*0.2e-3; %linspace(0,192*0.2e-3,128);
 xsonde = linspace(0,192*0.2e-3,192);
 xsonde = xsonde - mean(xsonde) ;
 
 % retreive t0 correction :
 Zref = 0*mean(MyImage.L) ; % Zref : position supposed to be invariant in rotation

 

