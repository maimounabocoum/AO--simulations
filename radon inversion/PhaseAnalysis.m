%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% load experiemental data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
addpath('functions')
addpath('..\Scan routines')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update of the OP structure
                            % X : theta 
                            % Y : monotonic t variable in points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% experiemental input datas :
% load('experiement images - JB - test\OP0deg-2016-02-01_13-11.mat');
% c = 1540 ; % sound velocity in m/s
% MyImage = OP(data(:,:,1),X,Y,Param.SamplingRate*1e6,c); 

%% simulation traces 
load('saved images\Simulation.mat');
load('saved images\SimulationTransmission.mat');

N = 2^12;
MyImage = MyImage.InitializeFourier(N);
MyImage.Fmax()       % maximum frequency sampling = 1/dt
Lobject = 5e-3;
Fc = 100/Lobject;    % Lobject is the size of the object to detect. Using simple model (sinc function)
                     % we set it to kc = 100/Lobject
                   
%% Nyist principle states the sampling of the object to reconstruct to be such that w > w_max/2 

%%% Fourier Tranform of Radon input image with respect to t
% creating a fourier parameter set using the measure sample size in m

MyImage.F_R = MyImage.fourier(MyImage.R) ;
MyImage.Show_F_R(Fc); % Fc : cut-off frequency used for screening

% find index for angle = 0
N_theta0 = find(MyImage.theta == 0);

% extract spectrum for this value :
% plotyy(MyImage.f,abs(MyImage.F_R(:,N_theta0)),MyImage.f,unwrap( angle(MyImage.F_R(:,N_theta0)) ) )
% hold on 
%plot(MyImage.theta*180/pi,angle(MyImage.F_R(MyImage.N/2,:)))

































