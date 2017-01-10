function []=Ray_Behavior(N,OS,signal)

%====================================================================
% This function displays one ray of the FFT function in order to show that it is indeed 
% band limitted. This builds the justification for the polynomial interpolation. The program 
% prompts the user to choose a ray and then shows the sampling of its values with growing 
% sampling rate, while showing the accuracy obtained using a polynomial 
% approximation (spline).
%
% Synopsis: []=Ray_Behavior(N,OS)
% 
% Inputs: N   - Number of entries per row/column in the signal. Original X-Polar grid 
%                     will have  2N rays with 2N equi-spaced points along it.
%               OS - Oversampling factor along the rays. If we expect 2N points, there will be 
%                       2N*OS points after the oversampling. OS shopuld be a power of 2.
%               signal - the kind of signal to use: 1 - random, 2- square, 3 - checker board.
%
% Output - None
%
% Example: Ray_Behavior(16,256,3);
%
% Written by Miki Elad on March 20th, 2005.
%====================================================================

if nargin<1,
    N=16;
    OS=128;
end;

% stage 1 - creating a signal
if signal==1,
    X=randn(N)+sqrt(-1)*randn(N);
elseif signal==2,
    X=zeros(N,N);
    X(N/2-4:N/2+5,N/2-3:N/2+2)=1;
elseif signal==3,
    [kk,jj]=meshgrid(1:1:N);
    X=(-1).^(kk+jj)+0.01*sqrt(-1)*(-1).^(kk+jj);
end;

% stage 2 - building the overall S-Polar grid
VX=zeros(2*OS*N,N);
VY=zeros(2*OS*N,N);
ll=-N*OS:1:N*OS-1;
theta_y=pi*ll/N/OS;
m=-N/2:1:N/2-1;
for k=1:1:length(theta_y),
    theta_x=theta_y(k)*tan(m*pi/N/2);    
    VX(k,:)=theta_x;
    VY(k,:)=theta_y(k)*ones(1,N);
end;

HX=zeros(2*OS*N,N);
HY=zeros(2*OS*N,N);
ll=-N*OS:1:N*OS-1;
theta_x=pi*ll/N/OS;
m=-N/2+1:1:N/2;
for k=1:1:2*N*OS,
    theta_y=theta_x(k)*tan(m*pi/N/2);    
    HX(k,:)=theta_x(k)*ones(1,N);
    HY(k,:)=theta_y;
end;

GridX=[fliplr(VX),HX]; 
GridY=[fliplr(VY),HY];

figure(1); clf; 
plot(GridX(1:OS:end,:),GridY(1:OS:end,:),'.b');
axis equal;
axis([-4 4 -4 4]);
xlabel('The involved grid - please choose a ray');

% stage 3 - choose one ray
[x,y]=ginput(1);
Func=abs(GridX-x).^2+abs(GridY-y).^2;
[kk,jj]=find(min(Func(:))==Func);
RayX=GridX(:,jj(1));
RayY=GridY(:,jj(1));
hold on; 
plot(RayX(1:OS:end),RayY(1:OS:end,:),'.r');
xlabel('Red is the chosen ray');

% stage 4 - computing the exact FFT for this ray oversampled points
Y=Brute_Force_Transform(X,RayX,RayY);
figure(2); clf;
YY=fftshift(abs(fft(Y)));
YY=YY/max(YY);
semilogy(YY);

% stage 5 - compute the destination tranbsform values along an X-Polar grid ray with 2N points
Factor=sqrt(RayX(1)^2+RayY(1)^2)/pi;
Index_dest=-pi:2*pi/(2*N):pi;
Index_dest=Index_dest/Factor;
Index_dest=Index_dest(1:end-1);
Ydest=Brute_Force_Transform(X,RayX(1:OS:end)/Factor,RayY(1:OS:end)/Factor);
figure(1); plot(RayX(1:OS:end)/Factor,RayY(1:OS:end)/Factor,'g*');

% stage 6 - displaying the oversampling effects and interpolation error
figure(3); clf;
Index=-pi:2*pi/(2*N*OS):pi;
Index=Index(1:end-1);
subplot(2,1,1); plot(Index,real(Y),'r'); hold on; xlabel('real part');
subplot(2,1,2); plot(Index,imag(Y),'r'); hold on; xlabel('imaginary part');

for S=log(OS)/log(2):-1:0,
    Ysub=Y(1:2^S:end);
    Index=-pi:2*pi*2^S/(2*N*OS):pi;
    Index=Index(1:end-1);
    Yint=interp1(Index,Ysub,Index_dest,'spline');
    Error=sqrt(mean(abs(Yint-Ydest.').^2));
    subplot(2,1,1); plot(Index,real(Ysub),'.b');
    h1=plot(Index_dest,real(Yint),'go');
    subplot(2,1,2); plot(Index,imag(Ysub),'.b');
    h2=plot(Index_dest,imag(Yint),'go');
    subplot(2,1,1); 
    title(['Oversampling factor:',num2str(OS/2^S),'. Interpolation error: ',num2str(Error)]);
    pause;
    subplot(2,1,1); set(h1,'Visible','off');
    subplot(2,1,2); set(h2,'Visible','off');
end;
subplot(2,1,1); set(h1,'Visible','on');
subplot(2,1,2); set(h2,'Visible','on');

return;
