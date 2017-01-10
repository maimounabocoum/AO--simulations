function [Err1,Err2,Err3,Err4]=Error_Analysis_1(N,R,NN,HH,VV)
%====================================================================
% This function performs an error analysis for the Unequaly Spaced Fourier Transform.
% The basic idea is that given a 2D discrete singal X with N*N samples, its exact Fourier
% transform at an arbitrary place is linear w.r.t. the signal samples. When approximating
% this value using oversampled 2D FFT and interpolation we get again a linear operation. 
% Defining these two operations, we can find the worst signal to maximize the error. Four
% kinds of errors are discussed:
% 1. Error per location that finds the worst signal per location.
% 2. The same as above but forcing the worst-case signal to be real.
% 3. Error averaged over (NN+1)^2 locations, where one signal that is the worst on 
%     the average is found. 
% 4. As above with forcing the signal to be real.
% In all these cases we resort to an eigenvalue problem. In this routine we explore the 
% bilinear interpolation as the way to go from the oversampled Cartesian grid to the 
% desired USFT.
%
% Synopsis: [Err1,Err2,Err3,Err4]=Error_Analysis_1(N,R,NN,HH,VV)
%
% Inputs - 
%       N    - The number of samples in the signal to be treated (N*N) (default=8)
%       R    - The oversampling ratio in the USFFT (default=4)
%       NN - The number of samples to measure the error function. 
%       [HH,VV] - The block to work on in terms of numbers in the range [0,N]. For 
%                example, HH=[0,1] and VV=[0,1] work on the block [-pi:pi+2*pi/N]^2. 
%
% Output - 
%      Err1  - The error per location (varying worst case signal)
%      Err2  - The error per location (varying worst case signal while forcing it to be real)
%      Err3  - The error averaged (one worst case signal)
%      Err4  - The error averaged (one worst case real signal)
%
% Example:    
%     [E1,E2,E3,E4]=Error_Analysis_1(16,4,50,[0,1],[0,1]);
%                     This instruction takes a 16 by 16 signal, performs an upsampling by 4 and 
%                     the error is checked for the 1/16-by-1/16 block in the lower left corner. 
%                     This block is sampled 51^2 times and the error is computed per each of 
%                     this locaitions
%     [E1,E2,E3,E4]=Error_Analysis_1(8,4,100,[0,8],[0,8]);
%                     This instruction takes a 16 by 16 signal, performs an upsampling by 4 and 
%                     the error is checked for the 1/16-by-1/16 block in the lower left corner. 
%                     This block is sampled 51^2 times and the error is computed per each of 
%                     this locaitions
%
% Written by Miki Elad on March 20th, 2005. 
%===================================================================

if nargin<1,
    N=8; % size of the input signal 
    R=4; % size of oversampling of the 2D-FFT
    NN=20; % number of location inside the Delta*Delta block to sample for error measure
    HH=[0,8]; % horizontal span of the points to test
    VV=[0,8]; % vertical span of the points to test
elseif nargin<2,
    R=4; 
    NN=20; 
    HH=[0,8]; % horizontal span of the points to test
    VV=[0,8]; % vertical span of the points to test
elseif nargin<3,    
    NN=20; 
    HH=[0,8]; % horizontal span of the points to test
    VV=[0,8]; % vertical span of the points to test
end;

figure(1); clf;
[XC1,YC1]=Create_Grid('C',[N,N,-pi,pi,-pi,pi],'*r'); % The input Cartesian signal grid
XC1=[XC1,XC1; XC1,XC1]; % Extension to overcome boundary problems
YC1=[YC1,YC1; YC1,YC1];
hold on;
[XC2,YC2]=Create_Grid('C',[N*R,N*R,-pi,pi,-pi,pi],'.c'); % The oversampling Cart. grid
Delta1=2*pi/N;
Delta2=2*pi/N/R;
if XC1(1,HH(1)+1)<XC1(1,HH(2)+1), % no boundary problem
    HorPos=XC1(1,HH(1)+1):(HH(2)-HH(1))*Delta1/NN:XC1(1,HH(2)+1);
else,  % boundary problem
    HorPos=XC1(1,HH(1)+1):(HH(2)-HH(1))*Delta1/NN:XC1(1,HH(2)+1)+2*pi;
end;
if YC1(VV(1)+1,1)<YC1(VV(2)+1,1),
    VerPos=YC1(VV(1)+1,1):(VV(2)-VV(1))*Delta1/NN:YC1(VV(2)+1,1);
else,
    VerPos=YC1(VV(1)+1,1):(VV(2)-VV(1))*Delta1/NN:YC1(VV(2)+1,1)+2*pi;
end;
[x,y]=meshgrid(HorPos,VerPos); % The error assesment points grid
plot(x,y,'ob','Markersize',2);
title('The involved grids: Red - Input, Cyan - Upsampled, Blue - Estimated');
i=sqrt(-1);

Frn=Transform_Matrix(N,N,XC2,YC2); % The upsampled FFT transform 
                                                                                % [R^2N^2 by N^2] matrix
Err1=zeros(NN+1,NN+1);
Err2=zeros(NN+1,NN+1);
FFC=zeros((NN+1)^2,N^2);
FFM=zeros((NN+1)^2,N^2);
count=1;
[kk,jj]=meshgrid(0:1:N-1); % The input grid in integers
h=waitbar(0,'Sweeping through the estimation locations');
for Location1=1:1:NN+1,
    h=waitbar(Location1/(NN+1));
    for Location2=1:1:NN+1,

        Fc=exp(-i*kk.*x(Location2,Location1)-i*jj.*y(Location2,Location1));
        Fc=reshape(Fc,[1,N^2]); % The exact FT computation

        Bi=zeros(N*R,N*R); % The bilinear weighting vector
        Horizontal=max(find(x(Location2,Location1)-XC2(1,:)>=0));
        alpha=(x(Location2,Location1)-XC2(1,Horizontal))/Delta2;
        Vertical=max(find(y(Location2,Location1)-YC2(:,1)>=0));
        beta=(y(Location2,Location1)-YC2(Vertical,1))/Delta2;        
        Bi(Horizontal,Vertical)=(1-alpha)*(1-beta);
        if Vertical==R*N, % treating boundaries
            if Horizontal==R*N,
                Bi(Horizontal,1)=(1-alpha)*beta;
                Bi(1,Vertical)=alpha*(1-beta);
                Bi(1,1)=alpha*beta;
            else,
                Bi(Horizontal,1)=(1-alpha)*beta;
                Bi(Horizontal+1,Vertical)=alpha*(1-beta);
                Bi(Horizontal+1,1)=alpha*beta;
            end;
        else,
            if Horizontal==R*N,
                Bi(Horizontal,Vertical+1)=(1-alpha)*beta;
                Bi(1,Vertical)=alpha*(1-beta);
                Bi(1,Vertical+1)=alpha*beta;
            else,
                Bi(Horizontal,Vertical+1)=(1-alpha)*beta;
                Bi(Horizontal+1,Vertical)=alpha*(1-beta);
                Bi(Horizontal+1,Vertical+1)=alpha*beta;
            end;
        end;
        Bi=reshape(Bi,[1,R^2*N^2]);
        
        % Problem #1 - maximal error per each location separately - complex signal        
        Xopt=(Bi*Frn-Fc)';
        if norm(Xopt)>1e-6,
            Xopt=Xopt/sqrt(real(Xopt'*Xopt));
        end;
        Err1(Location2,Location1)=abs((Bi*Frn-Fc)*Xopt).^2;
  
        % Problem #2 - As above but while foceing Xopt to be real
        Matrix=[real(Bi*Frn-Fc); imag(Bi*Frn-Fc)];
        [U,D]=eig(Matrix'*Matrix);
        dd=diag(real(D));
        Pos=find(dd==max(dd));
        Pos=Pos(1);
        Xopt=U(:,Pos);
        Xopt=Xopt/sqrt(Xopt'*Xopt);
        Err2(Location2,Location1)=sum((Matrix*Xopt).^2);
        
        % Problem #3 & 4 - maximal error for all locations jointly
        FFC(count,:)=Fc;         % supposed to have (NN+1)^2 by N^2 entries
        FFM(count,:)=Bi*Frn; % supposed to have (NN+1)^2 by N^2 entries
        count=count+1;
        
    end;
end;
close(h);

% Problem 3 - finalyzing the error estimation
[U,D]=eig((FFM-FFC)'*(FFM-FFC));
dd=diag(real(D));
Pos=find(dd==max(dd));
Pos=Pos(1);
Xopt3=U(:,Pos);
Xopt3=Xopt3/sqrt(Xopt3'*Xopt3);
Err3=abs((FFM-FFC)*Xopt3).^2;
Err3=reshape(Err3,[NN+1,NN+1]);

% Problem 4 - finalyzing the error estimation
Matrix=real(FFM-FFC)'*real(FFM-FFC)+imag(FFM-FFC)'*imag(FFM-FFC);
[U,D]=eig(Matrix);
dd=diag(real(D));
Pos=find(dd==max(dd));
Pos=Pos(1);
Xopt4=U(:,Pos);
Xopt4=Xopt4/sqrt(Xopt4'*Xopt4);
Err4=abs((FFM-FFC)*Xopt4).^2;
Err4=reshape(Err4,[NN+1,NN+1]);

% Presenting the results
figure(2); clf; 
imagesc([Err1,zeros(NN+1,2),Err2,zeros(NN+1,2),Err3,zeros(NN+1,2),Err4]);
axis image; colormap(gray(256));
title('The errors as function of the location');

figure(3); clf; 
imagesc([reshape(real(Xopt3),[N,N]),reshape(imag(Xopt3),[N,N])]);
axis image; colormap(gray(256));
title('The worst image in terms of average error - real and imaginary parts');

figure(4); clf; 
imagesc(reshape(real(Xopt4),[N,N]));
axis image; colormap(gray(256));
title('The worst real image in terms of average error');
    
figure(5); clf; 
imagesc(fftshift(abs(fft2(reshape(Xopt4,[N,N]),256,256))));
colormap(jet);
title('The worst real image in terms of average error - Frequency domain');

disp('Error summary: ');
disp('=========================');
disp(['Per location (complex signal) ',num2str(sqrt(max(Err1(:))))]);
disp(['Per location (real signal)         ',num2str(sqrt(max(Err2(:))))]);
disp(['Averaged (complex signal)     ',num2str(sqrt(max(Err3(:))))]);
disp(['Averaged (real signal)             ',num2str(sqrt(max(Err4(:))))]);

return;