function [EE,TT,T1,T2]=Pseudo_Polar_Evolution(f);

%====================================================================
% Here we develop an evolution of stages from the brute-force grid based Transform and 
% to the fast Pseudo-Polar FT. The target is getting a fast PseudoPolar transform that will 
% match perfectly the results obtained by brute force.
% The tested methods are:
% 1. Brute force transform  with no vector multiplications - (like in C)
% 2. Brute force transform but with vector multiplications (exploiting Matlab capabilities)
% 3. First stage using regular 1D FFT and second stage using brute-force summations
% 4. First stage using regular 1D FFT and second stage using FR-FFT
%
% Note that the resulting fast transform is different from the one proposed by 
% Averbuch-Donoho in several minor asspects.
%
% Synopsis: [EE,TT,T1,T2]=Pseudo_Polar_Evolution(f);
%
% Input -    f -   The signal to work on. If its size is N*N, the transform will produce 2N*2N 
%                        values (actually, 2N rays and 2N elements on each of them).
% Output - EE - The errors between the various methods
%                 TT - The run-time measured per each method
%                 T1 - The full slow transform result
%                 T2 - The full fast transform result
% 
% Example:  
%           N=32; X=zeros(N,N); X(N/4:3*N/4,N/4:3*N/4)=1;
%           [EE,TT,T1,T2]=Pseudo_Polar_Evolution(X);
% 
% Written by Miki Elad on March 20th, 2005.
%====================================================================

%--------------------------------------------------------------------------------------------------------
%                                                  Stage 1 - Creating the signal
%--------------------------------------------------------------------------------------------------------

if nargin==0,   
    N=32; 
    f=randn(N,N);
else,
    N=size(f,1);
    if size(f,2)~=N,
        disp('The program expects square signals');
        return;
    end;
end;

%--------------------------------------------------------------------------------------------------------
%                                               Stage 2 - Basically Vertical Rays
%--------------------------------------------------------------------------------------------------------

% Creating the grid for the basically vertical rays
VX=zeros(2*N,N);
VY=zeros(2*N,N);
ll=-N:1:N-1;
theta_y=pi*ll/N;
m=-N/2:1:N/2-1;
for k=1:1:length(theta_y),
    theta_x=2*m*theta_y(k)/N;    
    VX(k,:)=theta_x;
    VY(k,:)=theta_y(k)*ones(1,N);
end;
[Nr,Nt]=size(VX);
figure(1); clf; plot(VX,VY,'+r'); hold on;
ii=sqrt(-1);

% Brute force transform for the basically vertical rays - detailed
tic; 
F1=zeros(Nr,Nt);
for ll=1:1:Nr,
    for mm=1:1:Nt,
        F1(ll,mm)=0;
        for k=0:1:N-1,
            for n=0:1:N-1,
                F1(ll,mm)=F1(ll,mm)+f(n+1,k+1)*exp(-ii*k*VX(ll,mm)-ii*n*VY(ll,mm));
            end;
        end;
    end;
end;
Time1=toc; 

% Brute force transform for the basically vertical rays as we do it
tic; 
[XX,YY]=meshgrid(0:1:N-1,0:1:N-1);
F2=zeros(Nr,Nt);
for k=1:1:Nr,
    for n=1:1:Nt,
        Map=exp(-ii*(VX(k,n).*XX+VY(k,n).*YY));
        F2(k,n)=sum(sum(Map.*f));
    end;
end;
Time2=toc; 

% Employing the FFT on the columns with zero padding - the basically vertical rays
tic; 
f_tilde=fft([f; zeros(N,N)],[],1); 
f_tilde=fftshift(f_tilde,1);
        % The above could also be done by:
        %       f_tilde=fft(f,2*N);
        %       f_tilde=[f_tilde(N+1:end,:); f_tilde(1:N,:)];

k=(0:1:N-1);
F3=zeros(Nr,Nt);
for ll=-N:1:N-1,
    for mm=-N/2:1:N/2-1,
        Map=exp(-ii*2*pi*k*mm/N*(ll/N));
        F3(ll+N+1,mm+N/2+1)=sum(f_tilde(ll+N+1,:).*Map);
    end;    
end;
Time3=toc;

% Employing the Fractional FT as well - the basically vertical rays
tic; 
f_tilde=fft([f; zeros(N,N)],[],1);
f_tilde=fftshift(f_tilde,1); 

F4=zeros(Nr,Nt);
for ll=-N:1:N-1,
     F4(ll+N+1,:)=My_FRFT_Centered(f_tilde(ll+N+1,:),ll/N^2).';
end;
Time4=toc;

%--------------------------------------------------------------------------------------------------------
%                                              Stage 3 - Basically Horizontal Rays
%--------------------------------------------------------------------------------------------------------


% Creating the grid for the basically horizontal rays
HX=zeros(2*N,N);
HY=zeros(2*N,N);
ll=-N:1:N-1;
theta_x=pi*ll/N;
m=-N/2+1:1:N/2;
for k=1:1:2*N,
    theta_y=2*m*theta_x(k)/N;    
    HX(k,:)=theta_x(k)*ones(1,N);
    HY(k,:)=theta_y;
end;
figure(1); clf; 
plot(VX,VY,'+r'); hold on;
plot(HX,HY,'+b'); 

% Brute force transform for the basically vertical rays - detailed
tic;
G1=zeros(Nr,Nt);
for ll=1:1:Nr,
    for mm=1:1:Nt,
        G1(ll,mm)=0;
        for k=0:1:N-1,
            for n=0:1:N-1,
                G1(ll,mm)=G1(ll,mm)+f(n+1,k+1)*exp(-ii*k*HX(ll,mm)-ii*n*HY(ll,mm));
            end;
        end;
    end;
end;
Time5=toc;

% Brute force transform for the basically horizontal rays as we do it
tic; 
[XX,YY]=meshgrid(0:1:N-1,0:1:N-1);
G2=zeros(Nr,Nt);
for k=1:1:Nr,
    for n=1:1:Nt,
        Map=exp(-ii*(HX(k,n).*XX+HY(k,n).*YY));
        G2(k,n)=sum(sum(Map.*f));
    end;
end;
Time6=toc;

% Employing the FFT on the columns with zero padding - the basically vertical rays
tic; 
f_tilde=fft([f, zeros(N,N)],[],2); 
f_tilde=fftshift(f_tilde,2);
f_tilde=f_tilde.';

n=(0:1:N-1);
G3=zeros(Nr,Nt);
for ll=-N:1:N-1,
    for mm=-N/2+1:1:N/2,
        Map=exp(-ii*2*pi*n*mm/N*(ll/N));
        G3(ll+N+1,mm+N/2)=sum(f_tilde(ll+N+1,:).*Map);
    end;    
end;
Time7=toc; 

% Employing the Fractional FT as well - the basically vertical rays
tic; 
f_tilde=fft([f, zeros(N,N)],[],2); 
f_tilde=fftshift(f_tilde,2);
f_tilde=f_tilde.';

G4=zeros(Nr,Nt);
for ll=-N:1:N-1,
     Factor=exp(i*2*pi*(0:1:N-1)*(N/2-1)*ll/N^2);
     G4(ll+N+1,:)=My_FRFT(f_tilde(ll+N+1,:).*Factor,ll/N^2).';
end;
Time8=toc;

%--------------------------------------------------------------------------------------------------------
%                                              Stage 4 - Summarizing Results
%--------------------------------------------------------------------------------------------------------

EE=[0; max(max(abs(F1-F2))); max(max(abs(F1-F3))); max(max(abs(F1-F4))); 
         0; max(max(abs(G1-G2))); max(max(abs(G1-G3))); max(max(abs(G1-G4)))]; 
TT=[Time1; Time2; Time3; Time4; Time5; Time6; Time7; Time8];

disp('    ');
disp('Treating the basically vertical rays');
disp('-----------------------------------------------------------------------------');
fprintf('Method:                                             Error:                         time:\n');
disp('-----------------------------------------------------------------------------');
fprintf('Brute force approach:                    REFERENCE            %12.8f \n ',Time1);
fprintf('Brute force explicit approach:  %12.8e     %12.8f \n ',EE(2),Time2');
fprintf('Method and using 1D FFT:        %12.8e     %12.8f \n ',EE(3),Time3');
fprintf('using 1D FFT and FFRFT:          %12.8e     %12.8f \n ',EE(4),Time4');
disp('-----------------------------------------------------------------------------');
disp('   ');
disp('Treating the basically horizontal rays');
disp('-----------------------------------------------------------------------------');
fprintf('Method:                                             Error:                         time:\n');
disp('-----------------------------------------------------------------------------');
fprintf('Brute force approach:                    REFERENCE           %12.8f \n ',Time5);
fprintf('Brute force explicit approach:  %12.8e     %12.8f \n ',EE(6),Time6');
fprintf('Method and using 1D FFT:        %12.8e     %12.8f \n ',EE(7),Time7');
fprintf('using 1D FFT and FFRFT:          %12.8e     %12.8f \n ',EE(8),Time8');

%--------------------------------------------------------------------------------------------------------
%                                    Stage 5 - Merging the two Transform sections
%--------------------------------------------------------------------------------------------------------

% Joining the coordinate systems into one continuous pair and visualizing the order
GX=[fliplr(VX),HX]; % we have to flip in order to get a smoothly rotating rays
GY=[fliplr(VY),HY];

figure(2); clf; 
subplot(1,2,1); 
h=plot(GX(:,1),GY(:,1),'.');
set(h,'Color',[0 0.1 0.1]);
xlabel('Progressing via the rows of the grid points');
axis image; 
axis([-pi pi -pi pi]*1.1);
hold on;

subplot(1,2,2);
h=plot(GX(1,:),GY(1,:),'.');
set(h,'Color',[0.1 0  0.1]);
xlabel('Progressing via the columns of the grid points');
axis image; 
axis([-pi pi -pi pi]*1.1);
hold on;

for k=1:1:2*N,
    subplot(1,2,1); 
    h=plot(GX(:,k),GY(:,k),'.');
    set(h,'Color',[k/(2*N+1),0.1,0.1]);
    
    subplot(1,2,2);
    h=plot(GX(k,:),GY(k,:),'.');
    set(h,'Color',[0.1,k/(2*N+1),0.1]);    
    pause;
end;

% Brute force full transform
tic; 
T1=zeros(2*N,2*N);
for ll=1:1:2*N,
    for mm=1:1:2*N,
        F1(ll,mm)=0;
        for k=0:1:N-1,
            for n=0:1:N-1,
                T1(ll,mm)=T1(ll,mm)+f(n+1,k+1)*exp(-ii*k*GX(ll,mm)-ii*n*GY(ll,mm));
            end;
        end;
    end;
end;
Time9=toc; 

% Fast Algorithm for both the basically vertical and horizontal rays
tic; 
T2=zeros(2*N,2*N);

f_tilde=fft([f; zeros(N,N)],[],1);
f_tilde=fftshift(f_tilde,1); 
for ll=-N:1:N-1,
     T2(ll+N+1,N:-1:1)=My_FRFT_Centered(f_tilde(ll+N+1,:),ll/N^2).';
end;

f_tilde=fft([f, zeros(N,N)],[],2); 
f_tilde=fftshift(f_tilde,2);
f_tilde=f_tilde.';
for ll=-N:1:N-1,
     Factor=exp(i*2*pi*(0:1:N-1)*(N/2-1)*ll/N^2);
     T2(ll+N+1,N+1:2*N)=My_FRFT(f_tilde(ll+N+1,:).*Factor,ll/N^2).';
end;
Time10=toc;

disp('             ');
disp('The merged Transform results:');
disp(['Error between brute force and fast algorithm: ',...
        num2str(max(max(abs(T1-T2))))]);
disp(['Run time results: Slow takes ',num2str(Time9),...
        ' sec. and Fast takes ',num2str(Time10), ' sec.']);

return;


