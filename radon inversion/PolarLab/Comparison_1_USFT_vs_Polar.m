function []=Comparison_1_USFT_vs_Polar(N,signal)

%====================================================================
%                                Comparison_1_USFT_vs_Polar
%
% For several specific signals (simple rectangle, and random signals), this script
% computes the exact XPolar transform and its approximation using the USFT
% methods (bilinear and spline+derivatives) and our fast XPolar transform. 
% L_1, L_2, and L_inf errors are then computed to show clear advantage to 
% the fast algorithm proposed. A graph that show this error as a function of 
% the Over-sampling factor is given. 
% 
% Synopsis: Comparison_1_USFT_vs_Polar(N,signal)
%
% Input: N - size of signal
%           Signal - type of signal to use (1 - random smooth, 2 for  -
%           random, 3 - square, 4 - random complex)
%
% Example: Comparison_1_USFT_vs_Polar(32,3);
%====================================================================

% Creating the signal
switch signal
case 1,
    X=randn(500);
    for k=1:1:3,
        X=conv2(X,ones(15,15),'valid');
    end;
    X=X(251:251+N-1,251:251+N-1);
    X=X/sqrt(sum(X(:).^2));
case 2,
    X=randn(N);
case 3,
    X=zeros(N); X(N/2-2:N/2+3,N/2-3:N-1)=1;
case 4,
    X=randn(N)+sqrt(-1)*randn(N);
case 5,
    X=imread('lena.tif');
    X=double(X); 
    X=imresize(X,[N,N]);
end;

figure(2); imagesc(X); colormap(gray(256)); axis image; colorbar;

% The ground truth computed brute-force
[Xc,Yc]=Create_Grid('X',[N,pi],'');
Yref=Brute_Force_Transform(X,Xc,Yc); 

% The aproximation using USFT and bilinear interpolation
% Err1=zeros(20,3);
% disp('USFT method with bilinear interpolation');
% for OS=1:1:20,
%     disp(OS);
%     Y=AFTUSF_NGP_1(X,Xc,Yc,OS); 
%     Err1(OS,1)=mean(abs(Yref(:)-Y(:)));
%     Err1(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2));
%     Err1(OS,3)=max(abs(Yref(:)-Y(:)));
% end;

% The aproximation using USFT and derivative-spline interpolation
Err2=zeros(20,3);
disp('USFT method with spline interpolation');
for OS=1:20,
    disp(OS);
    Y=AFTUSF_Spline_1(X,Xc,Yc,OS); 
    Err2(OS,1)=mean(abs(Yref(:)-Y(:)))/mean(abs(Yref(:)));
    Err2(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2))/sqrt(mean(abs(Yref(:)).^2));
    Err2(OS,3)=max(abs(Yref(:)-Y(:)))/max(abs(Yref(:)));
end;

% The approximation using our Polar FFT
disp('Polar FFT method with OS1=1');
Err3=zeros(10,3);
for OS=1:20,
    disp(OS);
    Y=XPolar_Transform(X,3,1,1,OS); 
    Err3(OS,1)=mean(abs(Yref(:)-Y(:)))/mean(abs(Yref(:)));
    Err3(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2))/sqrt(mean(abs(Yref(:)).^2));
    Err3(OS,3)=max(abs(Yref(:)-Y(:)))/max(abs(Yref(:)));
end;

disp('Polar FFT method with OS1=2');
Err4=zeros(10,3);
for OS=1:20,
    disp(OS);
    Y=XPolar_Transform(X,3,2,1,OS); 
    Err4(OS,1)=mean(abs(Yref(:)-Y(:)))/mean(abs(Yref(:)));
    Err4(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2))/sqrt(mean(abs(Yref(:)).^2));
    Err4(OS,3)=max(abs(Yref(:)-Y(:)))/max(abs(Yref(:)));
end;

disp('Polar FFT method with OS1=3');
Err5=zeros(10,3);
for OS=1:20,
    disp(OS);
    Y=XPolar_Transform(X,3,3,1,OS); 
    Err5(OS,1)=mean(abs(Yref(:)-Y(:)))/mean(abs(Yref(:)));
    Err5(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2))/sqrt(mean(abs(Yref(:)).^2));
    Err5(OS,3)=max(abs(Yref(:)-Y(:)))/max(abs(Yref(:)));
end;

disp('Polar FFT method with OS1=4');
Err6=zeros(20,3);
for OS=1:20,
    disp(OS);
    Y=XPolar_Transform(X,3,4,1,OS); 
    Err6(OS,1)=mean(abs(Yref(:)-Y(:)))/mean(abs(Yref(:)));
    Err6(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2))/sqrt(mean(abs(Yref(:)).^2));
    Err6(OS,3)=max(abs(Yref(:)-Y(:)))/max(abs(Yref(:)));
end;

disp('Polar FFT method with OS1=5');
Err7=zeros(20,3);
for OS=1:20,
    disp(OS);
    Y=XPolar_Transform(X,3,5,1,OS); 
    Err7(OS,1)=mean(abs(Yref(:)-Y(:)))/mean(abs(Yref(:)));
    Err7(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2))/sqrt(mean(abs(Yref(:)).^2));
    Err7(OS,3)=max(abs(Yref(:)-Y(:)))/max(abs(Yref(:)));
end;

disp('Polar FFT method with OS1=6');
Err8=zeros(20,3);
for OS=1:20,
    disp(OS);
    Y=XPolar_Transform(X,3,6,1,OS); 
    Err8(OS,1)=mean(abs(Yref(:)-Y(:)))/mean(abs(Yref(:)));
    Err8(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2))/sqrt(mean(abs(Yref(:)).^2));
    Err8(OS,3)=max(abs(Yref(:)-Y(:)))/max(abs(Yref(:)));
end;

disp('Polar FFT method with OS1=7');
Err9=zeros(20,3);
for OS=1:20,
    disp(OS);
    Y=XPolar_Transform(X,3,7,1,OS); 
    Err9(OS,1)=mean(abs(Yref(:)-Y(:)))/mean(abs(Yref(:)));
    Err9(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2))/sqrt(mean(abs(Yref(:)).^2));
    Err9(OS,3)=max(abs(Yref(:)-Y(:)))/max(abs(Yref(:)));
end;

disp('Polar FFT method with OS1=8');
Err10=zeros(20,3);
for OS=1:20,
    disp(OS);
    Y=XPolar_Transform(X,3,8,1,OS); 
    Err10(OS,1)=mean(abs(Yref(:)-Y(:)))/mean(abs(Yref(:)));
    Err10(OS,2)=sqrt(mean(abs(Yref(:)-Y(:)).^2))/sqrt(mean(abs(Yref(:)).^2));
    Err10(OS,3)=max(abs(Yref(:)-Y(:)))/max(abs(Yref(:)));
end;

figure(1); clf;
subplot(2,1,1);
semilogy((1:1:20).^2,Err2(:,1),(1:1:20)*1,Err3(:,1),...
                 (1:1:20)*2,Err4(:,1),(1:1:20)*3,Err5(:,1),(1:1:20)*4,Err6(:,1),...
                 (1:1:20)*5,Err7(:,1),(1:1:20)*6,Err8(:,1),(1:1:20)*7,Err9(:,1),...
                 (1:1:20)*8,Err10(:,1));
legend({'USFFT','Fast Polar with OS1=1',...
               'Fast Polar with OS1=2','Fast Polar with OS1=3',...
               'Fast Polar with OS1=4','Fast Polar with OS1=5',...
               'Fast Polar with OS1=6','Fast Polar with OS1=7',...
               'Fast Polar with OS1=8'});
title('L_1 norm approx. error');
xlabel('Over-Sampling factor');
grid on;

subplot(2,1,2);
semilogy((1:1:20).^2,Err2(:,2),(1:1:20)*1,Err3(:,2),...
                 (1:1:20)*2,Err4(:,2),(1:1:20)*3,Err5(:,2),(1:1:20)*4,Err6(:,2),...
                 (1:1:20)*5,Err7(:,2),(1:1:20)*6,Err8(:,2),(1:1:20)*7,Err9(:,2),...
                 (1:1:20)*8,Err10(:,2));
legend({'USFFT','Fast Polar with OS1=1',...
               'Fast Polar with OS1=2','Fast Polar with OS1=3',...
               'Fast Polar with OS1=4','Fast Polar with OS1=5',...
               'Fast Polar with OS1=6','Fast Polar with OS1=7',...
               'Fast Polar with OS1=8'});
title('L_{2} norm (square-root of MSE) approx. error');
xlabel('Over-Sampling factor');
grid on;

save Comparison_1_Results.mat

% subplot(2,1,2);
% semilogy((1:1:20).^2,Err2(:,3),(1:1:20)*1,Err3(:,3),...
%                  (1:1:20)*2,Err4(:,3),(1:1:20)*3,Err5(:,3),(1:1:20)*4,Err6(:,3));
% legend({'Bilinear USFT','Spline USFT','Fast XPOlar with OS1=1',...
%                'Fast XPOlar with OS1=2','Fast XPOlar with OS1=3','Fast XPOlar with OS1=4'});
% title('L_{\infty} norm approx. error');
% xlabel('Over-Sampling factor');
% grid on;

