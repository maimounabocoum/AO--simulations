function [TransRef,TransOut,Err]=Cartesian_2_Polar(Func,Ratio,Polarsize,Show)

%====================================================================
% This function gets a 2D discrete signal on a Cartesian grid and applies Polar-FT
% on it using two methods:
%    a. Direct transform that may take long, and 
%    b. Computing an oversampled 2D Cartesian FFT and bilinear interpolating.
% These two methods are compared.
%
% Synopsis: [TransRef,TransOut,Err]=Cartesian_2_Polar(Func,Ratio,Param)
%
% Inputs:
%   Func - a 2D discrete signal over a Cartesian grid
%   Ratio - amount of oversampling in the frequency domain of the Cartesian grid
%   Polarsize - parameters defining the polar grid to use - (Nr and Nt).
%   Show - 1 for showing the grids
%
% Outputs:
%   TransRef - The polar transform computed in a brute force way accurately
%   TransOut - The approximated polar transform by bilinear interpolation
%   Err - L1, L2, and L_infinity errors between TransRef and TransOut
%
% Example: This example shows that my code (Cartesian_2_Polar) performs 
%                  like Donohos' code (AFTUSF_NGP_1).
%
%           In=randn(32,32);
%           [Ref,Out1,Err]=Cartesian_2_Polar(In,4,[16,64],1);
%           figure(1);
%           [XP,YP]=Create_Grid('P',[16,64,pi],'.c');
%           Out2=AFTUSF_NGP_1(In,XP,YP,4);
%           disp(mean(abs(Ref(:)-Out1(:))));
%           disp(mean(abs(Ref(:)-Out2(:))));
%           disp(mean(abs(Out1(:)-Out2(:))));
%
% Written by Miki Elad on March 20th, 2005.
%====================================================================

[Ny,Nx]=size(Func);
disp(['Number of distinct samples in the input signal grid is: ',num2str(Nx*Ny)]);

% Creation of the grids involved
if Show==1, 
    figure(1); clf; 
    NNy=Ratio*Ny;
    NNx=Ratio*Nx;
    [GxIn,GyIn]=Create_Grid('C',[NNx,NNy,-pi,pi,-pi,pi],'.r'); 
    hold on;
    [GxOut,GyOut]=Create_Grid('P',[Polarsize,pi],'+b');
else,
    NNy=Ratio*Ny;
    NNx=Ratio*Nx;
    [GxIn,GyIn]=Create_Grid('C',[NNx,NNy,-pi,pi,-pi,pi],''); 
    [GxOut,GyOut]=Create_Grid('P',[Polarsize,pi],'');
end;

% Transforming the given signal using the Polar grid for reference 
TransRef=Brute_Force_Transform(Func,GxOut,GyOut);
disp(['Number of distinct samples in the Polar grid is: ',...
                                                       num2str(Distinct([GxOut(:),GyOut(:)]))]);

% Transforming the given signal with an oversampled Cartesian grid and FFT2
TransIn=fftn(Func,[NNx,NNy]); 
TransIn=fftshift(TransIn);
disp(['Number of distinct samples in the Cartesian grid is: ',num2str(NNx*NNy)]);

% Using TransIn in order to compute TransOut by Bilinear interpolation
[Ny,Nx]=size(GxOut);
Dx=(GxIn(1,2)-GxIn(1,1));
Dy=(GyIn(2,1)-GyIn(1,1));

TransOut=zeros(size(GxOut));
h=waitbar(0,'Interpolating ...');
for k=1:1:Ny,
    waitbar(k/Ny)
    for j=1:1:Nx, 
        dist=GxOut(k,j)-GxIn(1,:);
        PosX=max(find(dist>=0));
        alpha=dist(PosX)/Dx;
        dist=GyOut(k,j)-GyIn(:,1);
        PosY=max(find(dist>=0));
        beta=dist(PosY)/Dy;
        if PosX==NNx,
            if PosY==NNy,
                TransOut(k,j)=TransIn(PosY,PosX);                
            else,
                TransOut(k,j)=(1-beta)*TransIn(PosY,PosX)+beta*TransIn(PosY+1,PosX);                
            end;
        else,
            if PosY==NNy,
                TransOut(k,j)=(1-alpha)*TransIn(PosY,PosX)+...
                                            alpha*TransIn(PosY,PosX+1);                                                       
            else,
                TransOut(k,j)=(1-alpha)*(1-beta)*TransIn(PosY,PosX)+...
                                            alpha*(1-beta)*TransIn(PosY,PosX+1)+...
                                               (1-alpha)*beta*TransIn(PosY+1,PosX)+...                                   
                                                  alpha*beta*TransIn(PosY+1,PosX+1);                                           
            end;
        end;
    end;
end;
close(h);

% Results summary
if Show==1,
    figure(2); clf; 
    imagesc(log(abs([TransOut,TransRef,TransOut-TransRef])+1));
    colormap(gray(256));
end;

L1=sum(sum(abs(TransOut-TransRef)))/sum(sum(abs(TransRef)));
Li=max(max(abs(TransOut-TransRef)))/max(abs(TransRef(:)));
L2=sqrt(sum(sum(abs(TransOut-TransRef).^2)))/...
                            sqrt(sum(sum(abs(TransRef).^2)));
Err=[L1,L2,Li];
disp(['L_1, L_inf, and L_2 normalized norm errors are: ',num2str(Err)]);


return;