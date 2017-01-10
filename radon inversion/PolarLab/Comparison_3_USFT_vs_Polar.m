function []=Comparison_3_USFT_vs_Polar(N,OS1,OS2,signal,delta)

%====================================================================
% Comparing both the NGP and the spline USFT methods to the direct polar FFT using 
% the pseudo-polar approach. The comparison's objective is to show that the error in the 
% frequency domain is spatially dependent. We show this by first implementing the
% involved transforms on example images, and later by worst-case scenario analysis - 
% per each frequency point we find the signal (image) that brings the error to its maximum. 
%
% Synopsis: Comparison_3_USFT_vs_Polar(N,OS1,OS2,signal,delta)
% 
% Inputs:
%     N: size of the input array is N*N (Default - 16) 
%     OS1,OS2: Oversampling to use in both axes (default - OS1*OS2=80)
%                        In the USFFT case we use sqrt(OS1*OS2)
%     signal: Test signal to use (default - 3 for rectangular one part)
%     delta: Sampling distance in the frequency domain (default - 0.005) 
%
% Example:  Comparison_3_USFT_vs_Polar(16,4,20,3,0.05);
%
% Written by Miki Elad on March 20th, 2005.
%====================================================================

if nargin==0,
    N=16;
    OS1=4;
    OS2=20;
    signal==3;
    delta=0.01;
end;

% Building the transform matrices f or the exact the NGP, the spline, and
% the new polar-fft methods.

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

    % [Xc,Yc]=Create_Grid('X',[N,pi],'');
    % Approx_NGP=zeros(N^2*4,N^2);
    % count=1;
    % h=waitbar(0,'Building the NGP approximation matrix');
    % for k=1:1:N,
    %     waitbar(k/N);
    %     for j=1:1:N,
    %         x=zeros(N,N);
    %         x(k,j)=1;
    %         y=AFTUSF_NGP_1(x,Xc,Yc,ceil(sqrt(OS1*OS2)));
    %         Approx_NGP(:,count)=y(:);
    %         count=count+1;
    %     end;
    % end;
    % close(h);

    [Xc,Yc]=Create_Grid('X',[N,pi],'');
    Approx_Spline=zeros(N^2*4,N^2);
    count=1;
    h=waitbar(0,'Building the spline approximation matrix');
    for k=1:1:N,
        waitbar(k/N);
        for j=1:1:N,
            x=zeros(N,N);
            x(k,j)=1;
            y=AFTUSF_Spline_1(x,Xc,Yc,ceil(sqrt(OS1*OS2)));
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

%=================================================================
%                            W O R S T   C A S  E   A N A  L Y S I S 
%=================================================================

% Err_NGP=Trans-Approx_NGP; Y1=zeros(4*N^2,1);
Err_Spline=Trans-Approx_Spline; Y2=zeros(4*N^2,1);
Err_Polar=Trans-Approx_Polar; Y3=zeros(4*N^2,1);
for k=1:1:size(Err_Spline,1),
    %     vec=Err_NGP(k,:);
    %     if norm(vec)<1e-10, x=x*0; else, x=vec/norm(vec); end;
    %     Y1(k)=abs(vec*x');
    vec=Err_Spline(k,:);
    if norm(vec)<1e-10, x=x*0; else, x=vec/norm(vec); end;
    Y2(k)=sqrt(abs(vec*x'));
    vec=Err_Polar(k,:); 
    if norm(vec)<1e-10, x=x*0; else, x=vec/norm(vec); end;
    Y3(k)=sqrt(abs(vec*x'));
end;

data=[Xc(:),Yc(:),Y2,Y3];
L=length(-pi:delta:pi);
Result4=zeros(L,L);
Result5=zeros(L,L);
Result6=zeros(L,L);
ch=1; 
h=waitbar(0,'sweeping through 2D locations');
for hh=-pi:delta:pi,
    waitbar(ch/L);
    cv=1;
    for vv=-pi:delta:pi, 
        dist=(hh-data(:,1)).^2+(vv-data(:,2)).^2;
        pos=find(dist==min(dist));
        pos=pos(1);
        Result4(cv,ch)=data(pos,3);
        Result5(cv,ch)=data(pos,4);
        cv=cv+1;
    end;
    ch=ch+1;
end;
close(h);
[vv,hh]=meshgrid(-pi:delta:pi,-pi:delta:pi);
Mask=vv.^2+hh.^2<pi^2;

figure(1); clf; 
subplot(1,2,1); imagesc(Mask.*Result4); colorbar; axis off; axis image; colormap(gray(256));
title('Worst Error: USFFT vs. the exact');
xlabel('\xi_x'); ylabel('\xi_y'); 
subplot(1,2,2); imagesc(Mask.*Result5); colorbar; axis off; axis image; 
title('Worst Error: New method vs. exact');
xlabel('\xi_x'); ylabel('\xi_y'); 

figure(2); clf; 
imagesc([Mask.*Result4,ones(315,10)*max(Result5(:)),Mask.*Result5]); 
colorbar; axis off; axis image; colormap(gray(256)); 

return;