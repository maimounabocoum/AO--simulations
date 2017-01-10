function [Err1,Err2,Err3,Err4]=Error_Analysis_2(N,R,Polar)
%===================================================================
% This function performs an error analysis for the Polar Fourier Transform. The idea is 
% the same as that described in Error_Analysis_1.m. The difference is that the locations
% that the error is computed in are chosen as a polar grid. 
%
% Synopsis: [Err1,Err2,Err3,Err4]=Error_Analysis_2(N,R,Polar)
%
% Inputs - 
%       N    - The number of samples in the signal to be treated (N*N) (default=8)
%       R    - The oversampling ratio in the USFFT (default=4)
%       Polar - The parameters defining the polar grid to use for computing the error
%                Polar=[Kind - 1 for Polar and 0 for Recto-Polar; 
%                             Nr - number of points on each ray, Nt - number of angles, and 
%                             Range - the radius to use] (see Create_Grid for more details)
%
% Output - 
%      Err1  - The error per location (varying worst case signal)
%      Err2  - The error per location (varying worst case signal while forcing it to be real)
%      Err3  - The error averaged (one worst case signal)
%      Err4  - The error averaged (one worst case real signal)
%
% Example:    
%     [E1,E2,E3,E4]=Error_Analysis_2(16,4,[1,8,8,pi]);
%                     This instruction takes a 16 by 16 signal, performs an upsampling by 4 and 
%                     the error is checked for a polar grid that has
%     [E1,E2,E3,E4]=Error_Analysis_2(8,4,[1,8,8,pi]);
%                     This instruction takes a 16 by 16 signal, performs an upsampling by 4 and 
%                     the error is checked for the 1/16-by-1/16 block in the lower left corner. 
%                     This block is sampled 51^2 times and the error is computed per each of 
%                     this locaitions
%
% Written by Miki Elad on March 20th, 2004. 
%===================================================================

if nargin<1,
    N=8; % size of the input signal 
    R=4; % size of oversampling of the 2D-FFT
    Polar=[1,8,8,pi]; % Parameters of the destination grid to use
elseif nargin<2,
    R=4; 
    Polar=[1,8,8,pi]; % Parameters of the destination grid to use
elseif nargin<3,    
    Polar=[1,8,8,pi]; % Parameters of the destination grid to use
end;

figure(1); clf;
[XC1,YC1]=Create_Grid('C',[N,N,-pi,pi,-pi,pi],'*r'); % The input Cartesian signal grid
hold on;
[XC2,YC2]=Create_Grid('C',[N*R,N*R,-pi,pi,-pi,pi],'.c'); % The oversampling Cart. grid
Delta=2*pi/N/R;
if Polar(1)==1,
    [x,y]=Create_Grid('P',Polar(2:end),''); % The destination grid being Polar
else,
    [x,y]=Create_Grid('R',Polar(2:end),''); % The destination grid being Recto-Polar    
end;  
plot(x,y,'ob');
title('The involved grids: Red - Input, Cyan - Upsampled, Blue - Estimated');
i=sqrt(-1);

Frn=Transform_Matrix(N,N,XC2,YC2); % The upsampled FFT transform 
                                                                                % [R^2N^2 by N^2] matrix
Nr=Polar(2);
Nt=Polar(3);
Err1=zeros(Nt,Nr);
Err2=zeros(Nt,Nr);
FFC=zeros(Nr*Nt,N^2);
FFM=zeros(Nr*Nt,N^2);
count=1;
[kk,jj]=meshgrid(0:1:N-1); % The input grid in integers
h=waitbar(0,'Sweeping through the estimation locations');
for rrr=1:1:Nr,
    h=waitbar(rrr/Nr);
    for ttt=1:1:Nt,

        Fc=exp(-i*kk.*x(ttt,rrr)-i*jj.*y(ttt,rrr));
        Fc=reshape(Fc,[1,N^2]); % The exact FT computation

        Bi=zeros(N*R,N*R); % The bilinear weighting vector
        Horizontal=max(find(x(ttt,rrr)-XC2(1,:)>=0));
        alpha=(x(ttt,rrr)-XC2(1,Horizontal))/Delta;
        Vertical=max(find(y(ttt,rrr)-YC2(:,1)>=0));
        beta=(y(ttt,rrr)-YC2(Vertical,1))/Delta;        
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
        Err1(ttt,rrr)=abs((Bi*Frn-Fc)*Xopt).^2;
  
        % Problem #2 - As above but while foceing Xopt to be real
        Matrix=[real(Bi*Frn-Fc); imag(Bi*Frn-Fc)];
        [U,D]=eig(Matrix'*Matrix);
        dd=diag(real(D));
        Pos=find(dd==max(dd));
        Pos=Pos(1);
        Xopt=U(:,Pos);
        Xopt=Xopt/sqrt(Xopt'*Xopt);
        Err2(ttt,rrr)=sum((Matrix*Xopt).^2);
        
        % Problem #3 & 4 - maximal error for all locations jointly
        FFC(count,:)=Fc;         % supposed to have Nt*Nr by N^2 entries
        FFM(count,:)=Bi*Frn; % supposed to have Nt*Nr by N^2 entries
        count=count+1;
        
    end;
end;
close(h);

% Problem 3 - finalyzing the error estimation

% Note: A generalized eigenvalue problem would have been more suitable so as
%            to constrain the energy of the transformed input to be constant. Thus
%            solve : [U,D]=eig((FFM-FFC)'*(FFM-FFC),FFC'*FFC);

% [U,D]=eig((FFM-FFC)'*(FFM-FFC));
[U,D]=eig((FFM-FFC)'*(FFM-FFC),FFC'*FFC);
dd=diag(real(D));
Pos=find(dd==max(dd));
Pos=Pos(1);
Xopt3=U(:,Pos);
Xopt3=Xopt3/sqrt(Xopt3'*Xopt3);
Err3=abs((FFM-FFC)*Xopt3).^2;
Err3=reshape(Err3,[Nt,Nr]);

% Problem 4 - finalyzing the error estimation
Matrix=real(FFM-FFC)'*real(FFM-FFC)+imag(FFM-FFC)'*imag(FFM-FFC);
[U,D]=eig(Matrix);
dd=diag(real(D));
Pos=find(dd==max(dd));
Pos=Pos(1);
Xopt4=U(:,Pos);
Xopt4=Xopt4/sqrt(Xopt4'*Xopt4);
Err4=abs((FFM-FFC)*Xopt4).^2;
Err4=reshape(Err4,[Nt,Nr]);

% Presenting the results
figure(2); clf; 
imagesc([Err1,zeros(Nt,2),Err2,zeros(Nt,2),Err3,zeros(Nt,2),Err4]);
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