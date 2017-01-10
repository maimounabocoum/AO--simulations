function []=Test_FRFT(N,alpha)

%===================================================================
% Here we analyse the proper way to perform the Fractional Fourier
% Transform, by showing an evolution of methods, from brute-force to the
% fastest way. There are two parts to this program - one for the case where 
% the output is required for n=0,1,2, .... ,N-1, and the other for the
% centered case.
% 
% This program was written for debugging purpose only.
% 
% Synopsis:  Test_FRFT(N)
% 
% Input - N - length of the vector
%            alpha - shrink factor (~0.3/N).
%
% Written by Michael Elad on March 20th 2005.
%===================================================================

if nargin==0,
    N=24;
    alpha=0.2/N;
end;
x=randn(N,1);
i=sqrt(-1);

% Direct FRFT
y1=zeros(N,1);
for n=0:1:N-1,
    for k=0:1:N-1,
        y1(n+1)=y1(n+1)+x(k+1)*exp(-i*2*pi*k*n*alpha);
    end;
end;

% Still direct but with vector multiplication
y2=zeros(N,1);
k=(0:1:N-1);
for n=0:1:N-1,
    Map=exp(-i*2*pi*k'*n*alpha);
    y2(n+1)=sum(x.*Map);
end;
disp(max(abs(y1-y2)));

% The convolution term 
n=0:1:N-1;
Factor=exp(-i*pi*alpha*n'.^2);
x_tilde=x.*Factor;

y3=zeros(N,1);
k=(0:1:N-1);
for n=0:1:N-1,
    Map=exp(i*pi*(k-n)'.^2*alpha);
    y3(n+1)=sum(x_tilde.*Map);
end;
y3=y3.*Factor;
disp(max(abs(y1-y3)));

% The convolution via FFT 
n=[0:1:N-1, -N:1:-1]';
Factor=exp(-i*pi*alpha*n.^2);

x_tilde=[x; zeros(N,1)];
x_tilde=x_tilde.*Factor;    

XX=fft(x_tilde);
YY=fft(conj(Factor));
y4=ifft(XX.*YY);
y4=y4.*Factor;
y4=y4(1:N);

disp(max(abs(y1-y4)));

% Applying My_FRFT
y5=My_FRFT(x,alpha);
disp(max(abs(y1-y5)));

%=========================================================
% Here we analyse the proper way to perform the Fractional Fourier Transform
% ===> This part assumes that the output is required for n=-N/2, ... ,N/2-1
%=========================================================

if nargin==0,
    N=24;
    alpha=0.2/N;
end;
x=randn(N,1);
i=sqrt(-1);

% Direct FRFT
y1=zeros(N,1);
for n=-N/2:1:N/2-1,
    for k=0:1:N-1,
        y1(n+1+N/2)=y1(n+1+N/2)+x(k+1)*exp(-i*2*pi*k*n*alpha);
    end;
end;

% Still direct but with vector multiplication
y2=zeros(N,1);
k=(0:1:N-1);
for n=-N/2:1:N/2-1,
    Map=exp(-i*2*pi*k'*n*alpha);
    y2(n+1+N/2)=sum(x.*Map);
end;
disp(max(abs(y1-y2)));

% The convolution term 
n=-N/2:1:N/2-1;
Factor_N=exp(-i*pi*alpha*n'.^2);
k=0:1:N-1;
Factor_K=exp(-i*pi*alpha*k'.^2);
x_tilde=x.*Factor_K;

y3=zeros(N,1);
for n=-N/2:1:N/2-1,
    Map=exp(i*pi*(k-n)'.^2*alpha);
    y3(n+1+N/2)=sum(x_tilde.*Map);
end;
y3=y3.*Factor_N;
disp(max(abs(y1-y3)));

% The convolution via FFT 
n=[0:1:N-1, -N:1:-1]';
Factor1=exp(-i*pi*alpha*n.^2);
Factor2=exp(i*pi*(0:1:N-1)'*N*alpha);

x_tilde=x.*Factor2;
x_tilde=[x_tilde; zeros(N,1)];
x_tilde=x_tilde.*Factor1;    

XX=fft(x_tilde);
YY=fft(conj(Factor1));
y4=ifft(XX.*YY);
y4=y4.*Factor1;
y4=y4(1:N);

disp(max(abs(y1-y4)));

% Applying My_FRFT
Factor2=exp(i*pi*(0:1:N-1)'*N*alpha);
x_tilde=x.*Factor2;
y5=My_FRFT(x_tilde,alpha);
disp(max(abs(y1-y5)));

% Applying My_FRFT_Centered
y6=My_FRFT_Centered(x,alpha);
disp(max(abs(y1-y6)));
