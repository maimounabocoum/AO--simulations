%%%%%%%%%%%%%%%%% using conversion interation %%%%%%%%%%%%%%%%%%%

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
MyImage.Show_R();    % show Radon transform (ie interpolated raw data)
MyImage.Fmax()       % maximum frequency sampling = 1/dt
Lobject = 5e-3;
Fc = 100/Lobject;    % Lobject is the size of the object to detect. Using simple model (sinc function)

MyImage.F_R = MyImage.fourier(MyImage.R);
MyImage.Show_F_R(Fc); % Fc : cut-off frequency used for screening

% initial iteration :
 xsonde = linspace(0,128*0.2e-3,128);
 xsonde = xsonde - mean(xsonde) ;
 [X,Z]= meshgrid(xsonde,z_out);
 img = zeros(size(X,1),size(X,2),'like',X);
 
% constraint in fourier domain 
C = abs ( MyImage.fourier(MyImage.R) );   

for iteration = 1:10
    
  % fourier transform  
 FT_img = fft(img) ;
 MyImage.F_R = abs(MyImage.F_R).*angle() ;
    

    
    
end























%rmpath('..\Scan routines')
%rmpath('functions')

