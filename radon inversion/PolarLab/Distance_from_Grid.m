function []=Distance_from_Grid(N,M,Im)

%====================================================================
% This function computes the distance between the polar and its
% approximating grids (caresian and pseudo-polar). The distance is computed
% as a weighted Euclidean distance, matching for every polar coordinate its
% nearest neighbor in the approximating grid. The weights are chosen based
% on the relative energy that exists in each polar coordinate. 
%
% Synopsis: []=Distance_from_Grid(N,M,Im)
% 
% Inputs: N - size of blocks (N-by-N pixels) taken from the image
%              M - maximal oversampling 
%              Im - an input image to work on
%
% Output - None
%
% Example: Distance_from_Grid(32,20,Im);
%
% Written by Miki Elad on March 20th, 2005.
%====================================================================


if nargin==0,
    N=32; M=20;
    A=imread('Example_Image1.jpg');
    A=rgb2gray(A);
    A=double(A(1:10*N,1:10*N));
end;
[Nx,Ny]=size(A);
Nx=floor(Nx/N)*N;
Ny=floor(Ny/N)*N;
A=A(1:Nx,1:Ny);

[X0,Y0]=Create_Grid('X',[N,pi],'');
C0=[X0(:),Y0(:)];
L0=size(C0,1);

% Computing the relative weights based on energy
Weight=zeros(2*N,2*N);
h=waitbar(0,'Sweeping through blocks in the image');
count=0;
for k=1:N:Nx,
    for j=1:N:Ny,
        waitbar(count*N^2/(Nx*Ny));
        Block=A(k:N+k-1,j:N+j-1);
        Weight=Weight+abs(Brute_Force_Transform(Block,X0,Y0));
        count=count+1;
    end;
end;
close(h);
Weight=Weight/sum(Weight(:));
figure(1); clf; imagesc(Weight); axis image; colorbar;

% Computing the distance between polar and Cartesian grids (weighted)
Result1=zeros(M,1);
for S=1:1:M,
    [XC,YC]=Create_Grid('C',[N*S,N*S,-pi,pi,-pi,pi],'');
    CC=[XC(:),YC(:)];
    Lc=size(CC,1);
    DDD=[];
    for jj=1:1:L0,
        Dist=sum((CC-ones(Lc,1)*C0(jj,:)).^2,2);
        DDD=[DDD,min(Dist)];
    end;
    Result1(S)=sum(sqrt(DDD)*Weight(:));
    disp([S,Result1(S)]);
end;

% Computing the distance between psudopolar and Cartesian grids (weighted)
Result2=zeros(M,1);
for S=1:1:M,
    [XC,YC]=Create_Grid('D',[N*S/2,pi],''); 
    CC=[XC(:),YC(:)];
    Lc=size(CC,1);
    DDD=[];
    for jj=1:1:L0,
        Dist=sum((CC-ones(Lc,1)*C0(jj,:)).^2,2);
        DDD=[DDD,min(Dist)];
    end;
    Result2(S)=sum(sqrt(DDD)*Weight(:));
    disp([S,Result2(S)]);
end;
    
semilogy(1:1:M,Result1,'g',1:1:M,Result2,'b');



