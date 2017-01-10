
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% load experiemental data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
load('experiement images - JB - test\OP0deg-2016-02-01_13-11.mat');

% update of the OP structure
% X : theta 
% Y : monotonic t variable

MyImage = OP(data(:,:,1),X,Y); 
MyImage.Show_R();

%% Nyist principle states the sampling of the object to reconstruct to be such that w > w_max/2 


%%% Fourier Tranform of Radon input image with respect to t
N = 2^12;

MyImage = MyImage.t_Fourier_R();
MyImage.Show_TF_R();