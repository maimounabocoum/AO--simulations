function [E1,E2]=Error_Analysis_3(N,RectoPolar,Polar)
%====================================================================
% This function performs an error analysis for the Polar Fourier Transform, computed 
% based on Recto-Polar coordinates. The idea is the same as that described in 
% Error_Analysis_1.m and Error_Analysis_2.m. The difference from Error_Analysis_2.m 
% is that the RectoPolar FT are computed exactly (assuming we have a fast algorithm for 
% their computation), and then we interpolate to get the Polar points. The interpolation
% used is the simplest possible - Nearest-Neighbor.
%
% Synopsis: Error_Analysis_3(N,RectoPolar,Polar)
%
% Inputs - 
%       N    - The number of samples in the signal to be treated (N*N) (default=8)
%       RectoPolar - The RectoPolar coordinates to use for exact computation 
%                [Nr - number of points on each ray, Nt - number of angles, and 
%                  Range - the radius to use] (default=[8,32,pi])
%       Polar - The parameters defining the polar grid to use for computing the error
%                [Nr - number of points on each ray, Nt - number of angles, and 
%                  Range - the radius to use] (default=[8,32,pi])
%
% Output - 
%       E1 - the relative error, where the ratio is computed to the frequency domain energy
%       E2 - the relative error, where the ratio is computed to the spatial domain energy
%
% Example: 
%
%     for r=1:1:10, [E1,E2]=Error_Analysis_3(6,[8*r,32*r,pi],[8,32,pi]); end;
%
%     A 6*6 signal in the spatial domain is used and we are interested in computing its 
%     Polar FT using 8 points along each ray and 32 rays. Using a varying size recto-polar 
%     grid and with exact evaluations on it, we get the errors
% 
% Written by Miki Elad on March 20th, 2005. 
%====================================================================

if nargin<1,
    N=8; % size of the input signal 
    RectoPolar=[8,32,pi]; % size of oversampling of the 2D-FFT
    Polar=[8,32,pi]; % Parameters of the destination grid to use
elseif nargin<2,
    RectoPolar=[8,32,pi]; % size of oversampling of the 2D-FFT
    Polar=[8,32,pi]; % Parameters of the destination grid to use
elseif nargin<3,    
    Polar=[8,32,pi]; % Parameters of the destination grid to use
end;

figure(1); clf;
% Create_Grid('C',[N,N,-pi,pi,-pi,pi],'*r'); % The Cartesian grid for regular 2D-FFT
[RPX,RPY]=Create_Grid('R',RectoPolar,'c.'); % The destination grid being Polar
hold on;
[PX,PY]=Create_Grid('P',Polar,'ob'); % The destination grid being Polar
title('The involved grids: Red - Input, Cyan - Upsampled, Blue - Estimated');

FRP=Transform_Matrix(N,N,RPX,RPY); % The RectoPolar transform 
FP = Transform_Matrix(N,N,PX,PY); % The RectoPolar transform 

Bint=zeros(Polar(1)*Polar(2),RectoPolar(1)*RectoPolar(2));
count=1;
for kk=1:1:size(PX,1), % angle
    for jj=1:1:size(PX,2), % radius
        Dist=(RPX-PX(kk,jj)).^2+(RPY-PY(kk,jj)).^2;
        [kk0,jj0]=find(Dist==min(Dist(:)));
        kk0=kk0(1); jj0=jj0(1);
        Bint(count,(kk0-1)*RectoPolar(1)+jj0)=1;
        count=count+1;
    end; 
end;

% Error relative to the frequency domain energy
[U,D]=eig((FP-Bint*FRP)'*(FP-Bint*FRP),FP'*FP);
dd1=diag(real(D));
Pos1=find(dd1==max(dd1));
Pos1=Pos1(1);
Xopt1=U(:,Pos1);
Xopt1=Xopt1/sqrt(Xopt1'*Xopt1);
E1=dd1(Pos1);

% Error relative to the spatial domain energy
[U,D]=eig((FP-Bint*FRP)'*(FP-Bint*FRP));
dd2=diag(real(D));
Pos2=find(dd2==max(dd2));
Pos2=Pos2(1);
Xopt2=U(:,Pos2);
Xopt2=Xopt2/sqrt(Xopt2'*Xopt2);
E2=dd2(Pos2);

% results
fprintf('Frequency relative and spatial relative errors:  ');
fprintf(' %10.5f  %10.5f \n',E1,E2);

return;