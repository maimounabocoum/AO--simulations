%%=========================== define camera screen
x_cam = 1:1024;
y_cam = 1:620;
dpixel = 5.5e-6;

% np = number of points / pixels
np = 1 ;
Fmax = np/(dpixel); % Fourier cut_off (x2)
N = 2^(nextpow2( np*length(x_cam) )+ 1); % number of point in Fourier

clearvars np
%% =========================== generate fourier structure

F       = TF2D(N,Fmax,Fmax);
[X,Y]   = meshgrid(F.x,F.y);
[FX,FY] = meshgrid(F.fx,F.fy);

clearvars N Fmax

%% ========= initial gaussian field in (x,y) plane + plane wave tilted

S0   = exp(-(X.^2+Y.^2)/(2e-6)^2); % initial Main laser beam 
Sref = exp(1i*9e7*X);              % initial Ref laser beam 

%figure; imagesc(F.x,F.y,S0)
%% ====================== Phase randomization induced by scattering phantom

IS0 = F.fourier(S0);
% generation of tagged photons
PHASE = pi*rand(size(X)) ;
IS0_tag = IS0.*exp(1i*PHASE);
% generation of untagged photons
PHASE = pi*rand(size(X)) ;
IS0_untag = IS0.*exp(1i*PHASE);

%% ============================= field out of Phantom
S0_tag = F.ifourier(IS0_tag);
S0_untag = F.ifourier(IS0_untag);

% % Iris beam
D_iris = 50e-3;
IRIS = (X.^2 + Y.^2) < (D_iris/2)^2;
S1_tag    = S0_tag.*IRIS;
S1_untag  = S0_untag.*IRIS;

% Fourier transform;
IS1_tag = F.fourier(S1_tag);
IS1_untag = F.fourier(S1_untag);

% Propagation in Free Space (methode 1)
%     d = 100000e-3;
%     lambda = 780e-9;
%     %H = exp(-2*pi*d*(FX.^2 + FY.^2 - 1/lambda^2)^(1/2));
%     H = exp(2*1i*pi*d/lambda).*exp(1i*pi*lambda*d*(FX.^2 + FY.^2));
%     IS1 = IS1.*H ;
%     S2 = F.ifourier(IS1);

% Fresnel approximation (methode 2)
        d = 5000e-3;
        lambda = 780e-9;
        S2_tag = interp2(FX,FY,IS1_tag,X/(d*lambda),Y/(d*lambda)) ;
        S2_untag = interp2(FX,FY,IS1_untag,X/(d*lambda),Y/(d*lambda)) ;





%% ============================  inverse fourier transform
alpha = max(abs(S2_tag(:)));

Icam = 0*( real(S2_untag).^2 + real((S2_tag + 0.3*alpha*Sref)).^2 ) + 0e-20*rand(size(X)) ;

Ifour0 =  F.fourier(Icam)  ;

%  =========================== selection of first order ===============
fx_c = 40e3;
fy_c = 0e3;
Filter0 = ((FX-fx_c).^2 + (FY-fy_c).^2 <= (20e3)^2);
Ifour =  Ifour0.*Filter0;

figure(2)
subplot(221)
imagesc(F.x*1e3,F.y*1e3,Icam)
colorbar
xlabel('x(mm)')
ylabel('y(mm)')
title('tagged photon profile')

subplot(223)
imagesc(F.x*1e3,F.y*1e3,real(S2_tag).^2)
colorbar
xlabel('x(mm)')
ylabel('y(mm)')
title('tagged photon profile')

subplot(222)
imagesc(F.fx*1e-3,F.fy*1e-3, abs(Ifour) )
axis([-100 100 -100 100])
caxis([0 max(abs(Ifour(:)))/100])
colorbar
xlabel('f_x(mm^{-1})')
ylabel('f_y(mm^{-1})')
title('Fourier plane Icam')

% inverse fourier transform
Icam1 = F.ifourier(Ifour);
Icam1 = Icam1.*exp(-1i*2*pi*fx_c.*X).*exp(-1i*2*pi*fy_c.*Y);

subplot(224)
imagesc(F.x*1e3,F.y*1e3, abs(Icam1) )
caxis([0 max(abs(Icam1(:)))/1])
colorbar
xlabel('f_x(mm^{-1})')
ylabel('f_y(mm^{-1})')
title('Fourier plane Icam')
