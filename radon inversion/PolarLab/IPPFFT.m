function [Xt]=IPPFFT(Y,accuracy)

%====================================================================
% This function performs a Recto (Pseudo) Polar inverse transform. It uses
% a simple steepest-descent algorithm with pre-conditioning. 
% The preconditioner is used in the following way: while T (the transform)
% is not so well-behaved, W*T is better. Thus, instead of minimizing
% ||Tx-y||^2 we minimize ||WTx-Wy||^2 which leads to the iterative
% procedure x_{k+1}=x_k - \mu T'W'[WTx_k-Wy]=x_k - \mu T'W^2*[Tx_k-y] 
%   
% Synopsis: X=IPPFFT(Y)
%
% Inputs -    Y     2N*2N matrix in pseudo-polar grid, (N is assumed to be even)           
%                 accuracy - the norm of the gradient in which the algorithm stops
% Outputs - X      N*N matrix (Cartesian grid)
%
% Example: 
%       N=128; X=randn(N)+sqrt(-1)*randn(N); 
%       Y=PPFFT(X,1,1);
%       X1=IPPFFT(Y);
%       disp(norm(X1-X));
%
% Written on March 20th, 2005 by Michael Elad.
%====================================================================

if nargin==1,
    accuracy=1e-3;
end;

[N,M]=size(Y);
N=N/2;

W=sqrt(abs(-N:1:N-1)/2)/N;
W(N+1)=sqrt(1/8)/N;
W=ones(2*N,1)*W;

Xt=zeros(N,N); 
Delta=1;
count=0;
method=2;
while Delta>accuracy, 
    Err=W.*(PPFFT(Xt).'-Y.'); % transpose is needed - see PPFFT comments
    D=APPFFT(W.*Err).';
    Delta=norm(D);
    disp([count,Delta]);
    mu=1/N;
    Xt_New=Xt-mu*D; 
    Xt=Xt_New;
    count=count+1;
end;
disp(['Number of required iterations is ',num2str(count)]);

return;

%====================================================================
% This is an alternative code for debugging, using the direct transform matrix
%====================================================================

N=24; X=randn(N)+sqrt(-1)*randn(N);
[xc,yc]=Create_Oversampled_Grid('D',[N,pi,1,1],'');
T=Transform_Matrix(N,N,xc,yc);

% First, running without pre-conditioner
Xc=X(:);
Yc=T*Xc;
Xt=Xc*0;
mu=real(1.5/max(eig(T'*T)));
Er1=zeros(50,1);
for k=1:1:50,
    disp(k);
    Xt=Xt -mu*T'*(T*Xt-Yc);
    Er1(k)=sqrt(mean(abs(Xc-Xt).^2));
end;
Er1=[sqrt(mean(abs(Xc).^2)); Er1];
figure(1); clf; semilogy(0:1:50,Er1);

% Second, introducing a pre-conditioner
W=sqrt(abs(-N:1:N-1)/2)/N;
W(N+1)=sqrt(1/8)/N;
W=ones(2*N,1)*W;
W=diag(W(:));

Xc=X(:);
Yc=T*Xc;
Xt=Xc*0;
mu=real(1.5/max(eig(T'*W*W*T)));
Er2=zeros(50,1);
for k=1:1:50,
    disp(k);
    Xt=Xt -mu*T'*W*(W*(T*Xt-Yc));
    Er2(k)=sqrt(mean(abs(Xc-Xt).^2));
end;
Er2=[sqrt(mean(abs(Xc).^2)); Er2];
figure(1); hold on; semilogy(0:1:50,Er2,'r');















