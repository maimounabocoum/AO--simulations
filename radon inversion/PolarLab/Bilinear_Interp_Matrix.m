function [Matrix]=Bilinear_Interp_Matrix(N,R,Xd,Yd)

%====================================================================
% This function builds a matrix that takes an over-sampled uniform grid (N mulitplied 
% by R), and performs the bilinear interpolation of the values of the FFT for the 
% destination grid [Xd,Yd].
% 
% Remark: This function is shown to be exactly the same as Donoho's bilinear 
%                 interpolation (see "Compare_Bilinear_Interpolations").
%
% Synopsis: [Matrix]=Bilinear_Interp_Matrix(N,R,Xd,Yd)
%
% Inputs:     N      size of input signal grid is N*N
%                   R      oversampling factor to a grid NR*NR
%                   Xd    X positions of the new grid to evaluate FT values
%                   Yd    Y positions of the new grid to evaluate FT values
% Outputs:   Matrix - the bilinear interpolation matrix (very sparse)
%
% Written by Miki Elad on March 20th, 2005.
%====================================================================

[Xus,Yus]=Create_Grid('C',[N*R,N*R,-pi,pi,-pi,pi],''); % The oversampling Cart. grid

Matrix=zeros(length(Xd(:)),R^2*N^2);
[kk,jj]=meshgrid(0:1:N-1); % The input grid in integers
count=1; delta=2*pi/R/N; i=sqrt(-1);

for Location1=1:1:size(Xd,1),
    for Location2=1:1:size(Xd,2),        

        Bi=zeros(N*R,N*R); % The bilinear weighting vector
        Horizontal=max(find(Xd(Location1,Location2)-Xus(1,:)>=0));
        alpha=(Xd(Location1,Location2)-Xus(1,Horizontal))/delta;
        Vertical=max(find(Yd(Location1,Location2)-Yus(:,1)>=0));
        beta=(Yd(Location1,Location2)-Yus(Vertical,1))/delta;        
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
        
        Matrix(count,:)=reshape(Bi,[1,R^2*N^2]);
        count=count+1;
        
    end;
end;

return;