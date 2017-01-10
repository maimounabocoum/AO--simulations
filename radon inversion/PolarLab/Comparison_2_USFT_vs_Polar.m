function []=Comparison_2_USFT_vs_Polar(N,OS1,OS2)

%====================================================================
%                                Comparison_2_USFT_vs_Polar
%
% Here we explicitly derive the involved matrices that perform the exact and the
% approximate transformations. Then we use these to look at the worst-case
% errors, found by eigenvalue problems.
%
% Synopsis: Comparison_2_USFT_vs_Polar(N,OS1,OS2)
%
% Input: N - size of signal
%           OS1 - Oversampling along the angles
%           OS1 - Oversampling along the rays
%
% Example: Comparison_2_USFT_vs_Polar(32,4,20);
%====================================================================

if nargin==0,
    N=16;
    OS1=4;
    OS2=20;
end;
OS=ceil(sqrt(OS1*OS2));

Exists=1;
if Exists==1,
    load TransMatrices
else,
    [Xc,Yc]=Create_Grid('X',[N,pi],'');
    Trans=zeros(N^2*4,N^2);
    count=1;
    h=waitbar(0,'Building the exact Transform matrix');
    for k=1:1:N,
        waitbar(k/N);
        for j=1:1:N,
            x=zeros(N,N);
            x(k,j)=1;
            y=Brute_Force_Transform(x,Xc,Yc);
            Trans(:,count)=y(:);
            count=count+1;
        end;
    end;
    close(h);

    [Xc,Yc]=Create_Grid('X',[N,pi],'');
    Approx_Spline=zeros(N^2*4,N^2);
    count=1;
    h=waitbar(0,'Building the spline approximation matrix');
    for k=1:1:N,
        waitbar(k/N);
        for j=1:1:N,
            x=zeros(N,N);
            x(k,j)=1;
            y=AFTUSF_Spline_1(x,Xc,Yc,OS);
            Approx_Spline(:,count)=y(:);
            count=count+1;
        end;
    end;
    close(h);

    [Xc,Yc]=Create_Grid('X',[N,pi],'');
    Approx_Polar=zeros(N^2*4,N^2);
    count=1;
    h=waitbar(0,'Building the fast Polar approximation matrix');
    for k=1:1:N,
        waitbar(k/N);
        for j=1:1:N,
            x=zeros(N,N);
            x(k,j)=1;
            y=XPolar_Transform(x,3,OS1,1,OS2);
            Approx_Polar(:,count)=y(:);
            count=count+1;
        end;
    end;
    close(h);
end;

% Error analysis
E2=(Trans-Approx_Spline)'*(Trans-Approx_Spline);
E3=(Trans-Approx_Polar)'*(Trans-Approx_Polar);

% Result # 1 - Worst case error:
disp(['The Spline maximal error is:          ',num2str(max(abs(eig(E2))))]);
disp(['The Fast Polar maximal error is:   ',num2str(max(abs(eig(E3))))]);

% Result # 2 - Worst case signals:
[U2,D2]=eig(E2);
pos2=find(diag(D2)==max(diag(D2)));
disp(sqrt(abs(D2(pos2,pos2))));
[U3,D3]=eig(E3);
pos3=find(diag(D3)==max(diag(D3)));
disp(sqrt(abs(D3(pos3,pos3))));

figure(1); clf; 
subplot(2,3,1); imagesc(reshape(real(U2(:,pos2)),[N,N])); axis image;  axis off;
title('Worst image - real part');
subplot(2,3,2); imagesc(reshape(imag(U2(:,pos2)),[N,N])); axis image;  axis off;
title('Worst image - imaginary part');
subplot(2,3,3); imagesc(abs(freqz2(reshape(U2(:,pos2),[N,N]))));  axis image;  axis off;
title('Worst image - Frequency content');
subplot(2,3,4); imagesc(reshape(real(U3(:,pos3)),[N,N])); axis image;  axis off;
title('Worst image - real part');
subplot(2,3,5); imagesc(reshape(imag(U3(:,pos3)),[N,N])); axis image;  axis off;
title('Worst image - imaginary part');
subplot(2,3,6); imagesc(abs(freqz2(reshape(U3(:,pos3),[N,N]))));  axis image;  axis off;
title('Worst image - Frequency content');
colormap(gray(256));

% Result # 3 - Relative worst case error:
disp(['The Spline maximal error is:          ',num2str(sqrt(max(abs(eig(E2,Trans'*Trans)))))]);
disp(['The Fast Polar maximal error is:   ',num2str(sqrt(max(abs(eig(E3,Trans'*Trans)))))]);

% Result # 4 - Relative worst case signals:
[U5,D5]=eig(E2,Trans'*Trans);
pos5=find(diag(D5)==max(diag(D5)));
disp(sqrt(abs(D5(pos5,pos5))));
[U6,D6]=eig(E3,Trans'*Trans);
pos6=find(diag(D6)==max(diag(D6)));
disp(sqrt(abs(D6(pos6,pos6))));

figure(2); clf; 
subplot(2,3,1); imagesc(reshape(real(U5(:,pos5)),[N,N])); axis image;  axis off;
title('Worst image - real part');
subplot(2,3,2); imagesc(reshape(imag(U5(:,pos5)),[N,N])); axis image;  axis off;
title('Worst image - imaginary part');
subplot(2,3,3); imagesc(abs(freqz2(reshape(U5(:,pos5),[N,N]))));  axis image;  axis off;
title('Worst image - Frequency content');
subplot(2,3,4); imagesc(reshape(real(U6(:,pos6)),[N,N])); axis image;  axis off;
title('Worst image - real part');
subplot(2,3,5); imagesc(reshape(imag(U6(:,pos6)),[N,N])); axis image;  axis off;
title('Worst image - imaginary part');
subplot(2,3,6); imagesc(abs(freqz2(reshape(U6(:,pos6),[N,N]))));  axis image;  axis off;
title('Worst image - Frequency content');
colormap(gray(256));

% Result # 5 - Errors versus subspaces - Polar as reference
[U,D]=eig(E3);
dd=sqrt(abs(diag(D)));
[dd1,Pos]=sort(dd);
dd2=dd1;
for k=1:1:N^2,
    dd2(k)=sqrt(abs(diag(U(:,Pos(k))'*E2*U(:,Pos(k)))));
end;

figure(3); clf; 
semilogy(dd(Pos),'--');
hold on; axis on;
semilogy(dd2,'-');
grid on;
xlabel('The eigen-subspaces');
ylabel('MSE error');
axis([0 length(dd) 1e-10 1]);
legend({'The Polar-FFT method errors','The USFFT method errors'},4)
title('Direct Eigenvalue Analysis')

% Result # 6 - Errors versus subspaces - Polar as reference - Relative
[U,D]=eig(E3,Trans'*Trans);
dd=sqrt(abs(diag(D)));
[dd1,Pos]=sort(dd);
dd2=dd1;
for k=1:1:N^2,
    dd2(k)=sqrt(abs((U(:,Pos(k))'*E2*U(:,Pos(k)))/(U(:,Pos(k))'*Trans'*Trans*U(:,Pos(k)))));
end;

figure(4); clf; 
semilogy(dd(Pos),'--');
hold on; axis on;
semilogy(dd2,'-');
grid on;
xlabel('The eigen-subspaces');
ylabel('MSE error')
axis([0 length(dd) 1e-12 0.01]);
legend({'The Polar-FFT method errors','The USFFT method errors'},4)
title('Relative Eigenvalue Analysis')

%=========================================================================

% Result # 7 - Relative worst case signals within the circle:
lambda=1000;

[XC,YC]=Create_Grid('C',[4*N,4*N,-pi,pi,-pi,pi],''); 
F=Transform_Matrix(N,N,XC,YC);
Mask=(XC.^2+YC.^2)>=pi^2;
Mask=Mask(:);
Pos=find(Mask);
R=F(Pos,:);

[U7,D7]=eig(E2,Trans'*Trans+lambda*R'*R);
pos7=find(diag(D7)==max(diag(D7)));
pos7=pos7(1);
disp(sqrt(abs(D7(pos7,pos7))));
[U8,D8]=eig(E3,Trans'*Trans+lambda*R'*R);
pos8=find(diag(D8)==max(diag(D8)));
pos8=pos8(1);
disp(sqrt(abs(D8(pos8,pos8))));

figure(5); clf; 
subplot(2,3,1); imagesc(reshape(real(U7(:,pos7)),[N,N])); axis image;  axis off;
title('Worst image - real part');
subplot(2,3,2); imagesc(reshape(imag(U7(:,pos7)),[N,N])); axis image;  axis off;
title('Worst image - imaginary part');
subplot(2,3,3); imagesc(abs(freqz2(reshape(U7(:,pos7),[N,N]))));  axis image;  axis off;
title('Worst image - Frequency content');
subplot(2,3,4); imagesc(reshape(real(U8(:,pos8)),[N,N])); axis image;  axis off;
title('Worst image - real part');
subplot(2,3,5); imagesc(reshape(imag(U8(:,pos8)),[N,N])); axis image;  axis off;
title('Worst image - imaginary part');
subplot(2,3,6); imagesc(abs(freqz2(reshape(U8(:,pos8),[N,N]))));  axis image;  axis off;
title('Worst image - Frequency content');
colormap(gray(256));

% Result # 8 - Errors versus subspaces - Polar as reference - Relative
[U,D]=eig(E3,Trans'*Trans+lambda*R'*R);
dd=sqrt(abs(diag(D)));
[dd1,Pos]=sort(dd);
dd2=dd1;
for k=1:1:N^2,
    dd2(k)=sqrt(abs((U(:,Pos(k))'*E2*U(:,Pos(k)))/(U(:,Pos(k))'*Trans'*Trans*U(:,Pos(k))+...
        lambda*U(:,Pos(k))'*R'*R*U(:,Pos(k)))));
end;

figure(6); clf; 
semilogy(dd(Pos),'--');
hold on; axis on;
semilogy(dd2,'-');
grid on;
xlabel('The eigen-subspaces');
ylabel('MSE error')
legend({'The Polar-FFT method errors','The USFFT method errors'},4)
title('Relative Eigenvalue Analysis with Constraint')
axis([0 length(dd) 1e-12 0.01]);

return;



