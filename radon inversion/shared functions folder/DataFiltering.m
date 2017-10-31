function [I,z_out] = DataFiltering(MyImage)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N       = 2^12;
Lobject = 1e-3;
Fc      = 1/Lobject;  % Lobject is the size of the object to detect. Using simple model (sinc function)
                      % we set it to kc = 100/Lobject 
MyImage = MyImage.InitializeFourier(N,10*Fc);
%MyImage.Show_R();    % show Radon transform (ie interpolated raw data)
MyImage.Fmax()        % maximum frequency sampling = 1/dt
                 
%% Nyist principle states the sampling of the object to reconstruct to be such that w > w_max/2 

%%% Fourier Tranform of Radon input image with respect to t
% creating a fourier parameter set using the measure sample size in m

MyImage.F_R = MyImage.fourier(MyImage.R) ;
% MyImage = MyImage.PhaseCorrection(Fc);
% MyImage.Show_F_R(Fc); % Fc : cut-off frequency used for screening


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtered inverse fourier transform :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter options : 'ram-lak' (default) , 'cosine', 'hamming' , 'hann'
FilterType = 'ram-lak';%'ram-lak' 

filt = FilterRadon(MyImage.f, MyImage.N ,FilterType , Fc);
filt = filt(:);
FILTER = filt*ones(1,length(MyImage.theta));
% hold on
% plot(MyImage.f,filt)

%p = bsxfun(@times, p, H); % faster than for-loop
%  I = MyImage.ifourier(MyImage.F_R);
 %I = MyImage.ifourier(MyImage.F_R.*FILTER);
 I = MyImage.ifourier(MyImage.F_R);
% extract image back to initial size :
 [I,z_out] = ReduceDataSize( I,'y',MyImage.t,MyImage.L);%MyImage.L



end

