function [TT,CondRes]=Comparison_4_Polar_Corners

%====================================================================
% Comparison_4_Polar_Corners - ccomputing the condition number of the polar
%               FFT matrix to the same one with augmentation of zeros in
%               the corners. The idea is to show that the condition number
%               is very stable and low for such an additional constraint. 
%
% Synopsis: Comparison_4_Polar_Corners
% 
% Outputs:
%     TT - The polar transfrom matrices of varying sizes, along with the
%             nulling of the corners cartesian transform matrices 
%     CondRes - Condition number
%
% Example:  Comparison_4_Polar_Corners;
%
% Written by Miki Elad on March 20th, 2005.
%====================================================================

count=1;
TT=cell(11,2);
CondRes=zeros(11,2);

for N=4:2:26
    disp(N);
    
    % Building the transform matrices for the exact polar-fft methods.
    [Xc,Yc]=Create_Grid('X',[N,pi],''); 
    TT{count,1}=Transform_Matrix(N,N,Xc,Yc);
        
    % Compute Additional Part to refer to the NULLED corners 
    [Xc,Yc]=Create_Grid('C',[4*N,4*N,-pi,pi,-pi,pi],'');
    TransCart=Transform_Matrix(N,N,Xc,Yc);
    Weight=Xc.^2+Yc.^2>pi^2;
    Pos=find(Weight(:));
    TransCart=TransCart(Pos,:);
    TT{count,2}=[TT{count,1}; TransCart];
        
    % Compute Condition Number
    CondRes(count,1)=cond(TT{count,1});
    CondRes(count,2)=cond(TT{count,2});    
    count=count+1;
    
end;

figure(1); clf; 
semilogy(4:2:26,CondRes(:,1),'r');
hold on;
semilogy(4:2:26,CondRes(:,2),'b');
grid on;
xlabel('N [Input array size is N \times N]');
ylabel('Condition-Number');

return;