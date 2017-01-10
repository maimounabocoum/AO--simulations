function [Mat]=Transform_Matrix(NN1,NN2,GridX,GridY)

%=====================================================================
% This function builds the Discrete Fourier Transform matrix for a cartesian 
% signal of size NN1*NN2  on a required grid defined by [GridX,GridY].
%
% Synopsis: [Mat]=Transform_Matrix(NN1,NN2,GridX,GridY)
%
% Input: 
%    In - The input discrete signal given on the Cartesian grid (spatial domain)
%    GridX,GridY - two 2D arrays containing the grid points for the output 
%                   (frequency domain) 
%
% Output:
%    Mat - The matrix that multiplies the input in order to create the transform
%
% Example:
%           [XC,YC]=Create_Grid('C',[24,24,-pi,pi,-pi,pi],''); 
%           In=zeros(32,32); In(14:30,12:23)=1;
%           D=Transform_Matrix(32,32,XC,YC);
%           Out1=reshape(D*In(:),[24,24]).';
%           Out2=Brute_Force_Transform(In,XC,YC);
%           disp(mean(mean(abs(Out1-Out2))));
%
%           [XP,YP]=Create_Grid('P',[32,64,pi],''); 
%           In=zeros(12,12);
%           In(1:4,8:11)=1;
%           D=Transform_Matrix(12,12,XP,YP);
%           Out1=reshape(D*In(:),[32,64]).';
%           Out2=Brute_Force_Transform(In,XP,YP);
%           disp(mean(mean(abs(Out1-Out2))));
%
% Showing the condition number as a function of N for the pseudo-polar grid
%           for N=1:1:25, 
%               disp(N);
%               [XC,YC]=Create_Grid('D',[N,pi],''); 
%               D=Transform_Matrix(N,N,XC,YC); 
%               CC(N)=cond(D); 
%           end;
%           plot(CC);
%
% Written by Michael Elad on March 20th, 2005.
%=====================================================================

[N1,N2]=size(GridX);
[XX,YY]=meshgrid(0:1:NN2-1,0:1:NN1-1);
ii=sqrt(-1);

if N1*N2*NN1*NN2>20*1e6,
    disp('Too big  matrix.');
    Mat=[];
    return;
end;

Mat=zeros(N1*N2,NN1*NN2);
count=1;
for k=1:1:N1,
    for j=1:1:N2,
        Map=exp(-ii*(GridX(k,j).*XX+GridY(k,j).*YY));
        Mat(count,:)=reshape(Map,[1,NN1*NN2]);
        count=count+1;
    end;
end;

return;