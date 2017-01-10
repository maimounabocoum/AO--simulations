% This script compares Donoho's bilinear interpolation to mine. 

In=zeros(8,8); In(3:7,2:5)=1; % the signal
In=randn(8,8);
%----------------------------------------------------------------------------------------------------
%                                                           My Interpolation
%----------------------------------------------------------------------------------------------------

% original grid
N=8; % size of the input signal 
figure(1); clf;
[Xo,Yo]=Create_Grid('C',[N,N,-pi,pi,-pi,pi],'*r'); % The input Cartesian signal grid
hold on;

% oversampling grid
R=3; % size of oversampling of the 2D-FFT
[Xus,Yus]=Create_Grid('C',[N*R,N*R,-pi,pi,-pi,pi],'.c'); % The oversampling Cart. grid

% destination grid
NN=17;
[Xd,Yd]=Create_Grid('C',[NN,NN,-pi,pi,-pi,pi],'.b'); % The input Cartesian signal grid
% [Xd,Yd]=Create_Grid('P',[NN,4*NN,pi],'.b'); % The input Cartesian signal grid
title('The involved grids: Red - Input, Cyan - Upsampled, Blue - Estimated');

% Building the involved matrices
Ftrue=Transform_Matrix(N,N,Xd,Yd); 

Fupsamp=Transform_Matrix(N,N,Xus,Yus); % The upsampled transform
Bilinear=Bilinear_Interp_Matrix(N,R,Xd,Yd);

In_vec=In(:);
Out1_mine=Ftrue*In_vec;
Out1_mine=reshape(Out1_mine,[size(Xd,2),size(Xd,1)]).';

Out1_ref=Brute_Force_Transform(In,Xd,Yd);
Error1=Out1_mine-Out1_ref;
disp('The error between the exact transform computed by brute force ');
disp('    and the one computed by the matrix manipulations:');
disp(mean(abs(Error1(:))));

Out2_mine=Bilinear*Fupsamp*In_vec;
Out2_mine=reshape(Out2_mine,[size(Xd,2),size(Xd,1)]).';

%----------------------------------------------------------------------------------------------------
%                                                       Donoho's Interpolation
%----------------------------------------------------------------------------------------------------

Out2_donoho=AFTUSF_NGP_1(In,Xd,Yd,R); 
Error2=Out2_mine-Out2_donoho;
disp('The error between the approximated transform computed by Donohos code ');
disp('    and the one computed by my matrix manipulations:');
disp(mean(abs(Error2(:))));

