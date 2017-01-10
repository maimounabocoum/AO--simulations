function [Out]=SPolar_Transform(In,method,OS)

%====================================================================
% This function performs a Fourier Transfrom over the s-polar grid (which is a half-way 
% towrads the polar-FFT). Several methods are tested:
%             a. exact
%             b. oversampling (direct and through FRFT) and bilinear interpolation
%             c. oversampling and spline interpolation
%             d. oversampling, computing derivatives, and spline interpolation
%
% Synopsis: [Out]=SPolar_Transform(In,method,OS)
% 
% Input:     In - The Input as a N*N array holding the signal.
%               method - as above from 1 to 4.
%               OS - Oversampling factor
% Output:   Out - The Output as a 2N by 2N new values that correpsonds to a new 
%                         S-Polar grid (having uniform angle sampling).
%
% Examples:
%
%   1. In this example we compare run-time and results between several ways to 
%       compute the S-Polar transform (mainly for debugging means). We see that the 
%       derivative-spline method is best, and higher oversampling is needed in order to 
%       get a similar result using regular spline.
%       a. Fastest method - Starting with fast recto-polar and then conversion to S-polar.
%           N=32; X=randn(N,N); 
%           tic; Y1=SPolar_Transform(X,1,[]); toc;
%           tic; Y2=SPolar_Transform(X,2,12); toc;
%           tic; Y3=SPolar_Transform(X,3,12); toc;
%           tic; Y4=SPolar_Transform(X,4,12); toc;
%       b. Brute force transform to the destination grid.
%           [Xc,Yc]=Create_Grid('S',[N,pi],'.r');
%           tic; Yref=Brute_Force_Transform(X,Xc,Yc); toc;
%       disp(['Error between exact and approx.: ',num2str(max(abs(Y1(:)-Y2(:))))]);
%       disp(['Error between exact and approx.: ',num2str(max(abs(Y1(:)-Y3(:))))]);
%       disp(['Error between exact and approx.: ',num2str(max(abs(Y1(:)-Y4(:))))]);
%       disp(['Error between fast and slow (accurate): ',num2str(max(abs(Y1(:)-Yref(:))))]);
%
% Written by Miki Elad on March 20th, 2005. 
%===================================================================

N=size(In,1); % assuming square input
if size(In,2)~=N,
    disp('The program expects square signals');
    return;
end;
ii=sqrt(-1);

%--------------------------------------------------------------------------------------------------------
%                                        Stage 1 - Basically Vertical Rays
%--------------------------------------------------------------------------------------------------------

% FFT on the columns with zero padding
f_tilde=fft([In; zeros(N,N)],[],1); 
f_tilde=fftshift(f_tilde,1);

switch method,
case 1, % brute-force exact transform for the rows
    k=(0:1:N-1);
    F=zeros(2*N,N);
    for ll=-N:1:N-1,
        for mm=-N/2:1:N/2-1,
            Map=exp(-ii*k*pi*ll/N*tan(pi*mm/N/2));
            % Map=exp(-ii*2*pi*k*mm/N*(ll/N)); % in the Pseudo-Polar case
            F(ll+N+1,mm+N/2+1)=sum(f_tilde(ll+N+1,:).*Map);
        end;    
    end;
case 2, % oversampling by brute force and then interpolation
    k=(0:1:N-1);
    Xcor1=2*pi/N*(-OS*N/2:1:OS*N/2-1)/OS/N;
    Xcor2=pi/N*tan(pi*(-N/2:1:N/2-1)/N/2);
    F=zeros(2*N,N);
    for ll=-N:1:N-1,
        Temp=zeros(1,OS*N); % the oversampled vector
        for mm=-OS*N/2:1:OS*N/2-1,
            Map=exp(-ii*2*pi*k*mm/N*(ll/N)/OS); 
            Temp(mm+OS*N/2+1)=sum(f_tilde(ll+N+1,:).*Map);
        end;
        if ll~=0,
            F(ll+N+1,1:N)=interp1(Xcor1*ll,Temp,Xcor2*ll,'spline');
        else, 
            F(ll+N+1,1:N)=Temp(1:OS:end);
        end;
    end;
case 3, % oversampling by fast FRFTand then interpolation
    Xcor1=2*pi/N*(-OS*N/2:1:OS*N/2-1)/OS/N;
    Xcor2=pi/N*tan(pi*(-N/2:1:N/2-1)/N/2);
    F=zeros(2*N,N);
    for ll=-N:1:N-1,
        Temp=My_FRFT_Centered([f_tilde(ll+N+1,:),zeros(1,OS*N-N)],ll/N^2/OS).';
        if ll~=0,
            F(ll+N+1,1:N)=interp1(Xcor1*ll,Temp,Xcor2*ll,'spline');
        else, 
            F(ll+N+1,1:N)=Temp(1:OS:end);
        end;
    end;        
case 4, % oversampling by fast FRFTand then interpolation using derivatives
    Xcor1=2*pi/N*(-OS*N/2:1:OS*N/2-1)/OS/N;
    Xcor2=pi/N*tan(pi*(-N/2:1:N/2-1)/N/2);
    F=zeros(2*N,N);
    for ll=-N:1:N-1,
        Vec=[f_tilde(ll+N+1,:),zeros(1,OS*N-N)];
        Temp1=My_FRFT_Centered(Vec,ll/N^2/OS).';
        Vec=Vec.*(-ii*[0:1:OS*N-1]); 
        Temp2=My_FRFT_Centered(Vec,ll/N^2/OS).'; % first derivative
        Vec=Vec.*(-ii*[0:1:OS*N-1]);
        Temp3=My_FRFT_Centered(Vec,ll/N^2/OS).'; % second derivative
        if ll~=0,
            F(ll+N+1,1:N)=interp_SPLINE(Xcor1*ll,Temp1,Temp2,Temp3,Xcor2*ll);
        else, 
            F(ll+N+1,1:N)=Temp1(1:OS:end);
        end;
    end;            
end;

%--------------------------------------------------------------------------------------------------------
%                                      Stage 2 - Basically Horizontal Rays
%--------------------------------------------------------------------------------------------------------

% FFT on the rows with zero padding
f_tilde=fft([In, zeros(N,N)],[],2); 
f_tilde=fftshift(f_tilde,2);
f_tilde=f_tilde.';

switch method,
case 1, % brute-force exact trabsform for the rows
    n=(0:1:N-1);
    G=zeros(2*N,N);
    for ll=-N:1:N-1,
        for mm=-N/2+1:1:N/2,
            Map=exp(-ii*pi*n*ll/N*tan(pi*mm/N/2));
            % Map=exp(-ii*2*pi*n*mm/N*(ll/N)); % in the Pseudo-Polar case
            G(ll+N+1,mm+N/2)=sum(f_tilde(ll+N+1,:).*Map);
        end;    
    end;
case 2, % oversampling by brute force and then bilinear interpolation
    n=(0:1:N-1);
    Ycor1=2*pi/N*(-OS*N/2+1:1:OS*N/2)/OS/N;
    Ycor2=pi/N*tan(pi*(-N/2+1:1:N/2)/N/2);
    G=zeros(2*N,N);
    for ll=-N:1:N-1,
        Temp=zeros(1,OS*N); % the oversampled vector
        for mm=-OS*N/2+1:1:OS*N/2,
            Map=exp(-ii*2*pi*n*mm/N*(ll/N)/OS); 
            Temp(mm+OS*N/2)=sum(f_tilde(ll+N+1,:).*Map);
        end;
        if ll~=0,
            G(ll+N+1,1:N)=interp1(Ycor1*ll,Temp,Ycor2*ll,'spline');
        else,
            G(ll+N+1,1:N)=Temp(1:OS:end);
        end;
    end;    
case 3, % oversampling by fast FRFTand then interpolation
    n=(0:1:N-1);
    Ycor1=2*pi/N*(-OS*N/2+1:1:OS*N/2)/OS/N;
    Ycor2=pi/N*tan(pi*(-N/2+1:1:N/2)/N/2);
    G=zeros(2*N,N);
    for ll=-N:1:N-1,
        Factor=exp(i*2*pi*(0:1:OS*N-1)*(OS*N/2-1)*ll/N^2/OS);
        Temp=My_FRFT([f_tilde(ll+N+1,:),zeros(1,OS*N-N)].*Factor,ll/N^2/OS).';
        if ll~=0,
            G(ll+N+1,1:N)=interp1(Ycor1*ll,Temp,Ycor2*ll,'spline');
        else,
            G(ll+N+1,1:N)=Temp(1:OS:end);
        end;
    end;        
case 4, % oversampling by fast FRFTand then interpolation using derivatives
    Ycor1=2*pi/N*(-OS*N/2+1:1:OS*N/2)/OS/N;
    Ycor2=pi/N*tan(pi*(-N/2+1:1:N/2)/N/2);
    G=zeros(2*N,N);
    for ll=-N:1:N-1,
        Factor=exp(i*2*pi*(0:1:OS*N-1)*(OS*N/2-1)*ll/N^2/OS);
        Vec=[f_tilde(ll+N+1,:),zeros(1,OS*N-N)].*Factor;
        Temp1=My_FRFT(Vec,ll/N^2/OS).';
        Vec=Vec.*(-ii*[0:1:OS*N-1]); 
        Temp2=My_FRFT(Vec,ll/N^2/OS).'; % first derivative
        Vec=Vec.*(-ii*[0:1:OS*N-1]);
        Temp3=My_FRFT(Vec,ll/N^2/OS).'; % second derivative
        if ll~=0,
            G(ll+N+1,1:N)=interp_SPLINE(Ycor1*ll,Temp1,Temp2,Temp3,Ycor2*ll);
        else,
            G(ll+N+1,1:N)=Temp1(1:OS:end);
        end;
    end;        
end;

%--------------------------------------------------------------------------------------------------------
%                             Stage 3 - Merging the two Transform sections
%--------------------------------------------------------------------------------------------------------

Out=[fliplr(F),G]; % we have to flip in order to get a smoothly rotating rays

return;
