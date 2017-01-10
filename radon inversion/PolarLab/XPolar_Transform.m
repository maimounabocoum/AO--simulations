function [Out]=XPolar_Transform(In,method1,OS1,method2,OS2)

%===============================================================================
% This function performs a Fourier Transfrom over the X-polar grid. Several methods are 
% tested for the S-Polar transform and then on the final resampling stage.
%
% Synopsis: [Out]=XPolar_Transform(In,method1,OS1,method2,OS2)
% 
% Input:     In           - The Input as a N*N array holding the signal.
%               method1 - for the S-Polar transform (see 'Spolar_Transform.m')
%               OS1        - Oversampling factor for the S-Polar transform
%               method2 - for the final resampling of the rays (not applicable yet).
%               OS2        - Oversampling factor for resampling of the rays.
% Output:   Out        - The Output as a 2N by 2N new values that correpsonds to a new 
%                                   X-Polar grid (Polar grid with 2N samples along each ray).
%
% Examples:
%
%   1. In this example we compare run-time and results between several ways to 
%       compute the X-Polar transform (mainly for debugging means). 
%       a. Fastest method - Starting with fast recto-polar and then conversion to S-polar.
%           N=100; X=randn(N,N); 
%           tic; Y1=XPolar_Transform(X,1,5,1,5); toc;
%           tic; Y2=XPolar_Transform(X,2,5,1,5); toc;
%           tic; Y3=XPolar_Transform(X,3,5,1,5); toc;
%       b. Brute force transform to the destination grid.
%           [Xc,Yc]=Create_Grid('X',[N,pi],'.r');
%           tic; Yref=Brute_Force_Transform(X,Xc,Yc); toc;
%       disp(['Error between exact and approx.: ',num2str(max(abs(Y1(:)-Y2(:))))]);
%       disp(['Error between exact and approx.: ',num2str(max(abs(Y1(:)-Y3(:))))]);
%       disp(['Error between fast and slow (accurate): ',num2str(max(abs(Y1(:)-Yref(:))))]);
%
% Written by Miki Elad on March 20th, 2005. 
%===============================================================================

N=size(In,1); % assuming square input
if size(In,2)~=N,
    disp('The program expects square signals');
    return;
end;
ii=sqrt(-1);

%--------------------------------------------------------------------------------------------------------
%                                        Stage 1 - Basically Vertical Rays
%--------------------------------------------------------------------------------------------------------

% A. FFT on the columns with zero padding
f_tilde=fft([In; zeros((2*OS2-1)*N,N)],[],1); 
f_tilde=fftshift(f_tilde,1);

% B. FFT on the rows while putting the grid points on the S-Polar
switch method1,
case 1, % brute-force exact transform for the rows
    k=(0:1:N-1);
    F=zeros(2*OS2*N,N);
    h=waitbar(0,'sweeping the rows');
    for ll=-N*OS2:1:N*OS2-1,
        waitbar((ll+N*OS2)/(2*N*OS2));
        for mm=-N/2:1:N/2-1,
            Map=exp(-ii*k*pi*ll/(N*OS2)*tan(pi*mm/N/2));
            F(ll+OS2*N+1,mm+N/2+1)=sum(f_tilde(ll+N*OS2+1,:).*Map);
        end;    
    end;
    close(h);
case 2, % oversampling by fast FRFTand then interpolation
    Xcor1=2*pi/(N*OS2)*(-OS1*N/2:1:OS1*N/2-1)/OS1/N;
    Xcor2=pi/(N*OS2)*tan(pi*(-N/2:1:N/2-1)/N/2);
    F=zeros(2*OS2*N,N);
    h=waitbar(0,'sweeping the rows');
    for ll=-N*OS2:1:N*OS2-1,
        waitbar((ll+N*OS2)/(2*N*OS2));
        Temp=My_FRFT_Centered([f_tilde(ll+OS2*N+1,:),zeros(1,OS1*N-N)],ll/N^2/OS1/OS2).';
        if ll~=0,
            F(ll+OS2*N+1,1:N)=interp1(Xcor1*ll,Temp,Xcor2*ll,'spline');
        else, 
            F(ll+OS2*N+1,1:N)=Temp(1:OS1:end);
        end;
    end;        
    close(h);
case 3, % oversampling by fast FRFTand then interpolation using derivatives
    Xcor1=2*pi/(N*OS2)*(-OS1*N/2:1:OS1*N/2-1)/OS1/N;
    Xcor2=pi/(N*OS2)*tan(pi*(-N/2:1:N/2-1)/N/2);
    F=zeros(2*OS2*N,N);
    h=waitbar(0,'sweeping the rows');
    for ll=-N*OS2:1:N*OS2-1,
        waitbar((ll+N*OS2)/(2*N*OS2));
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
    close(h);
end;

% C. resampling the rays to the X-Polar locations
Fout=zeros(2*N,N);
Xcor=-pi:2*pi/(2*N*OS2):pi-pi/(2*N*OS2); % these are the projected locations
for k=1:1:N,
    Ray=F(:,k);
    Factor=cos((k-N/2-1)*pi/(2*N));
    Fout(:,k)=interp1(Xcor/Factor,Ray,Xcor(1:OS2:end),'spline').';
end;

%--------------------------------------------------------------------------------------------------------
%                                      Stage 2 - Basically Horizontal Rays
%--------------------------------------------------------------------------------------------------------

% A. FFT on the rows with zero padding
f_tilde=fft([In, zeros(N,(2*OS2-1)*N)],[],2); 
f_tilde=fftshift(f_tilde,2);
f_tilde=f_tilde.';

% B. FFT on the columns while putting the grid points on the S-Polar
switch method1,
case 1, % brute-force exact trabsform for the rows
    n=(0:1:N-1);
    G=zeros(2*OS2*N,N);
    h=waitbar(0,'sweeping the columns');
    for ll=-N*OS2:1:N*OS2-1,
        waitbar((ll+N*OS2)/(2*N*OS2));
        for mm=-N/2+1:1:N/2,
            Map=exp(-ii*pi*n*ll/(N*OS2)*tan(pi*mm/N/2));
            G(ll+OS2*N+1,mm+N/2)=sum(f_tilde(ll+OS2*N+1,:).*Map);
        end;    
    end;
    close(h);
case 2, % oversampling by fast FRFTand then interpolation
    Ycor1=2*pi/(N*OS2)*(-OS1*N/2+1:1:OS1*N/2)/OS1/N;
    Ycor2=pi/(N*OS2)*tan(pi*(-N/2+1:1:N/2)/N/2);
    G=zeros(2*N*OS2,N);
    h=waitbar(0,'sweeping the columns');
    for ll=-N*OS2:1:N*OS2-1,
        waitbar((ll+N*OS2)/(2*N*OS2));
        Factor=exp(i*2*pi*(0:1:OS1*N-1)*(OS1*N/2-1)*ll/N^2/OS1/OS2);
        Temp=My_FRFT([f_tilde(ll+OS2*N+1,:),zeros(1,OS1*N-N)].*Factor,ll/N^2/OS1/OS2).';
        if ll~=0,
            G(ll+OS2*N+1,1:N)=interp1(Ycor1*ll,Temp,Ycor2*ll,'spline');
        else,
            G(ll+OS2*N+1,1:N)=Temp(1:OS1:end);
        end;
    end;        
    close(h);    
case 3, % oversampling by fast FRFTand then interpolation using derivatives
    Ycor1=2*pi/(N*OS2)*(-OS1*N/2+1:1:OS1*N/2)/OS1/N;
    Ycor2=pi/(N*OS2)*tan(pi*(-N/2+1:1:N/2)/N/2);
    G=zeros(2*OS2*N,N);
    h=waitbar(0,'sweeping the columns');
    for ll=-N*OS2:1:N*OS2-1,
        waitbar((ll+N*OS2)/(2*N*OS2));
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
    close(h);
end;

% C. resampling the rays to the X-Polar locations
Gout=zeros(2*N,N);
Xcor=-pi:2*pi/(2*N*OS2):pi-pi/(2*N*OS2); % these are the projected locations
for k=1:1:N,
    Ray=G(:,k);
    Factor=cos((k-N/2)*pi/(2*N));
    Gout(:,k)=interp1(Xcor/Factor,Ray,Xcor(1:OS2:end),'spline').';
end;

%--------------------------------------------------------------------------------------------------------
%                             Stage 3 - Merging the two Transform sections
%--------------------------------------------------------------------------------------------------------

Out=[fliplr(Fout),Gout]; % we have to flip in order to get a smoothly rotating rays

return;
