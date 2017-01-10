function [Out]=PFFT(In,OS1,OS2)

%=====================================================================
% This function performs a Fourier Transfrom over the polar grid. It uses
% Hermite interpolation along the angles, and spline one for the rays
%
% Synopsis: [Out]=PFFT(In,OS1,OS2)
% 
% Input:     In           - The Input as a N*N array holding the signal.
%               OS1        - Oversampling factor for the S-Polar transform
%                                 (equiangle but still with squares)
%               OS2        - Oversampling factor for resampling of the rays.
% Output:   Out        - The Output as a 2N by 2N new values that correpsonds to a new 
%                                 X-Polar grid (Polar grid with 2N samples along each ray).
%
% Example:
%   1. In this example we compare run-time and results between several ways to 
%       compute the X-Polar transform (mainly for debugging means). 
%       a. Fastest method - Starting with fast recto-polar and then conversion to S-polar.
%           N=100; X=randn(N,N); 
%           tic; Y=PFFT(X,5,5); toc;
%       b. Brute force transform to the destination grid.
%           [Xc,Yc]=Create_Grid('X',[N,pi],'.r');
%           tic; Yref=Brute_Force_Transform(X,Xc,Yc); toc;
%           imagesc(abs([Yref,Yref-y])); axis image; axis off; 
%
% Written by Miki Elad on March 20th, 2005. 
%=====================================================================

if nargin==1, 
    OS1=3;
    OS2=8;
end;

N=size(In,1); % assuming square input
if size(In,2)~=N,
    disp('The program expects square signals');
    return;
end;
ii=sqrt(-1);

%------------------------------------------------------------------------------------------------------------------------------
%                                        Stage 1 - Basically Vertical Rays
%------------------------------------------------------------------------------------------------------------------------------

h=waitbar(0,'Applying the transform');

% A. FFT on the columns with zero padding
f_tilde=fft([In; zeros((2*OS2-1)*N,N)],[],1); 
f_tilde=fftshift(f_tilde,1);

% B. FFT on the rows while putting the grid points on the S-Polar
Xcor1=2*pi/(N*OS2)*(-OS1*N/2:1:OS1*N/2-1)/OS1/N;
Xcor2=pi/(N*OS2)*tan(pi*(-N/2:1:N/2-1)/N/2);
F=zeros(2*OS2*N,N);
for ll=-N*OS2:1:N*OS2-1,
    waitbar((ll+N*OS2)/(6*N*OS2));
    Vec=[f_tilde(ll+OS2*N+1,:),zeros(1,OS1*N-N)];
    alpha=ll/N^2/OS1/OS2;
    Temp1=My_FRFT_Centered(Vec,alpha).';
    Vec=Vec.*(-ii*[0:1:OS1*N-1]);
    Temp2=My_FRFT_Centered(Vec,alpha).'; % first derivative
    Vec=Vec.*(-ii*[0:1:OS1*N-1]);
    Temp3=My_FRFT_Centered(Vec,alpha).'; % second derivative
    if ll~=0,
        F(ll+OS2*N+1,1:N)=interp_SPLINE(Xcor1*ll,Temp1,Temp2,Temp3,Xcor2*ll);
    else,
        F(ll+OS2*N+1,1:N)=Temp1(1:OS1:end);
    end;
end;

% C. resampling the rays to the X-Polar locations
Fout=zeros(2*N,N);
Xcor=-pi:2*pi/(2*N*OS2):pi-pi/(2*N*OS2); % these are the projected locations
for k=1:1:N,
    waitbar((k*OS2+2*N*OS2)/(6*N*OS2));
    Ray=F(:,k);
    Factor=cos((k-N/2-1)*pi/(2*N));
    Fout(:,k)=interp1(Xcor/Factor,Ray,Xcor(1:OS2:end),'spline').';
end;

%------------------------------------------------------------------------------------------------------------------------------
%                                      Stage 2 - Basically Horizontal Rays
%------------------------------------------------------------------------------------------------------------------------------

% A. FFT on the rows with zero padding
f_tilde=fft([In, zeros(N,(2*OS2-1)*N)],[],2); 
f_tilde=fftshift(f_tilde,2);
f_tilde=f_tilde.';

% B. FFT on the columns while putting the grid points on the S-Polar
Ycor1=2*pi/(N*OS2)*(-OS1*N/2+1:1:OS1*N/2)/OS1/N;
Ycor2=pi/(N*OS2)*tan(pi*(-N/2+1:1:N/2)/N/2);
G=zeros(2*OS2*N,N);
for ll=-N*OS2:1:N*OS2-1,
    waitbar((ll+N*OS2+3*N*OS2)/(6*N*OS2));
    alpha=ll/N^2/OS1/OS2;
    Factor=exp(i*2*pi*(0:1:OS1*N-1)*(OS1*N/2-1)*ll/N^2/OS1/OS2);
    Vec=[f_tilde(ll+OS2*N+1,:),zeros(1,OS1*N-N)].*Factor;
    Temp1=My_FRFT(Vec,alpha).';
    Vec=Vec.*(-ii*[0:1:OS1*N-1]);
    Temp2=My_FRFT(Vec,alpha).'; % first derivative
    Vec=Vec.*(-ii*[0:1:OS1*N-1]);
    Temp3=My_FRFT(Vec,alpha).'; % second derivative
    if ll~=0,
        G(ll+OS2*N+1,1:N)=interp_SPLINE(Ycor1*ll,Temp1,Temp2,Temp3,Ycor2*ll);
    else,
        G(ll+OS2*N+1,1:N)=Temp1(1:OS1:end);
    end;
end;

% C. resampling the rays to the X-Polar locations
Gout=zeros(2*N,N);
Xcor=-pi:2*pi/(2*N*OS2):pi-pi/(2*N*OS2); % these are the projected locations
for k=1:1:N,
    waitbar((OS2*k+5*N*OS2)/(6*N*OS2));
    Ray=G(:,k);
    Factor=cos((k-N/2)*pi/(2*N));
    Gout(:,k)=interp1(Xcor/Factor,Ray,Xcor(1:OS2:end),'spline').';
end;

%------------------------------------------------------------------------------------------------------------------------------
%                             Stage 3 - Merging the two Transform sections
%------------------------------------------------------------------------------------------------------------------------------

Out=[fliplr(Fout),Gout]; % we have to flip in order to get a smoothly rotating rays
close(h);

return;
