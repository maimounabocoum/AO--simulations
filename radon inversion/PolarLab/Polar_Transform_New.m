function [Out]=Polar_Transform_New(In,method1,method2,OS1,OS2)

%====================================================================
% This function performs a Fourier Transfrom over the Polar grid, using either spline or 
% Hermite interpolation for the two stages. As opposed to the XPolar_Transform routine, 
% this code starts by regular Pseudo-Polar transform with oversampling by zero padding. 
% Then the interpolations are done by working on the rows and the columns of the 
% produced array. A shortcoming of this algorithm is that the oversampling on the two 
% axes (r,theta) is forced to be the same.
%
% Synopsis: [Out]=Polar_Transform_New(In,method1,method2,OS1,OS2)
%
% Input:     In            - The Input as a N*N array holding the signal.
%               method1 - kind of interpolation for the concentric squares (1-spline, 2-Hermite).
%               method2 - kind of interpolation for the rays (1-spline, 2-Hermite).
%               OS1        - Processing oversampling factor.
%               OS2        - Output oversampling factor.
% Output:   Out        - The Output as a 2N*OS2 by 2N*OS2 new values that correpsonds 
%                                 to a Polar grid (Polar grid with 2N*OS2 samples along each ray, 
%                                 and 2N*OS2 rays).
%
% Examples:
%
%   1. In this example we compare run-time and results between several ways to 
%       compute the Polar transform (mainly for debugging means). 
%       N=32; X=randn(N,N); 
%       a. Fast method
%           tic; Y1=Polar_Transform_New(X,1,1,4,1); toc;
%           tic; Y2=Polar_Transform_New(X,1,2,4,1); toc;
%           tic; Y3=Polar_Transform_New(X,2,1,4,1); toc;
%           tic; Y4=Polar_Transform_New(X,2,2,4,1); toc;
%       b. Brute force transform to the destination grid.
%           [Xc,Yc]=Create_Grid('X',[N,pi],'.r');
%           tic; Yref=Brute_Force_Transform(X,Xc,Yc); toc;
%       disp(['Error between exact and approx.: ',num2str(max(abs(Y1(:)-Yref(:))))]);
%       disp(['Error between exact and approx.: ',num2str(max(abs(Y2(:)-Yref(:))))]);
%       disp(['Error between exact and approx.: ',num2str(max(abs(Y3(:)-Yref(:))))]);
%       disp(['Error between exact and approx.: ',num2str(max(abs(Y4(:)-Yref(:))))]);
%
% Written by Miki Elad on September 12th, 2002. 
%===============================================================================

N=size(In,1); % assuming square input
if size(In,2)~=N,
    disp('The program expects square signals');
    return;
end;
ii=sqrt(-1);

%------------------------------------------------------------------------------------------------------------------------------
%                       Stage 1 - Starting by Pseudo-Polar FFT with zero-padded oversampling
%------------------------------------------------------------------------------------------------------------------------------

In=[In, zeros(N,N*(OS1-1))];
In=[In; zeros(N*(OS1-1),N*OS1)];
Out=RectoPolar_Trans_New(In);
% Out is the Pseudo-Polar FFT result on a grid with 2N*OS1 rays each containing 2N*OS1 samples

%------------------------------------------------------------------------------------------------------------------------------
%                            Stage 2 - Rotating the rays to become equaly spaced angelwise
%------------------------------------------------------------------------------------------------------------------------------

Out1=zeros(2*N*OS1,2*N);

% spline interpolation for the basically vertical part
Cor1=2*pi/(N*OS1)*(-OS1*N/2:1:OS1*N/2-1)/OS1/N;
Cor2=pi/(N*OS1)*tan(pi*(-N/2:1:N/2-1)/N/2);
for ll=-N*OS1:1:N*OS1-1,
    Vec=Out(ll+OS1*N+1,1:N*OS1).';
    if ll~=0,
        Out1(ll+OS1*N+1,1:N)=interp1(Cor1*ll,Vec,Cor2*ll,'spline');
    else, 
        Out1(ll+OS1*N+1,1:N)=Vec(1:OS1:end).';
    end;
end;

% spline interpolation for the basically horizntal part
Cor1=2*pi/(N*OS1)*(-OS1*N/2+1:1:OS1*N/2)/OS1/N;
Cor2=pi/(N*OS1)*tan(pi*(-N/2+1:1:N/2)/N/2);
for ll=-N*OS1:1:N*OS1-1,
    Vec=Out(ll+OS1*N+1,1+N*OS1:2*N*OS1).';
    if ll~=0,
        Out1(ll+OS1*N+1,1+N:2*N)=interp1(Cor1*ll,Vec,Cor2*ll,'spline');
    else, 
        Out1(ll+OS1*N+1,1+N:2*N)=Vec(1:OS1:end).';
    end;
end;

% spline interpolation along the basically vertical rays 
Cor=-pi:2*pi/(2*N*OS1):pi-pi/(2*N*OS1); 
for k=1:1:N,
    Vec=Out1(:,k);
    Factor=cos((k-N/2-1)*pi/(2*N));
    Out2(:,k)=interp1(Cor/Factor,Vec,Cor(1:OS1:end),'spline').';
end;

% spline interpolation along the basically horizontal rays 
Cor=-pi:2*pi/(2*N*OS1):pi-pi/(2*N*OS1); 
for k=N+1:1:2*N,
    Vec=Out1(:,k);
    Factor=cos((k-N/2)*pi/(2*N));
    Out2(:,k)=interp1(Cor/Factor,Vec,Cor(1:OS1:end),'spline').';
end;

Out=Out2;

return;
