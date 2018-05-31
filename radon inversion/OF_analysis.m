%% OF analysis %%
% MyImage = OF(x*1e-3,z*1e-3,Datas,SampleRate*1e6,c) ;
% MyImage = OP(x*1e-3,z*1e-3,Datas,SampleRate*1e6,c) ;
% addpath('shared functions folder')

%%
    MyImage = OP(Datas,(-20:20)*pi/180,z,SampleRate*1e6,c) ;
    R_FT = MyImage.fourier(MyImage.R) ;
    FILTER = MyImage.GetFILTER(1e-3) ;
    I = MyImage.ifourier(R_FT.*FILTER) ;
    %[I,z_out] = DataFiltering(MyImage) ;
    
    [I,z_out] = ReduceDataSize( I,'y',MyImage.t,MyImage.L);%MyImage.L

    X_m = (1:system.probe.NbElemts)*(0.2*1e-3) ;
    [theta,M0,X0,Z0] = EvalDelayLaw_shared(X_m,DelayLAWS,ActiveLIST,c); 

    Hresconstruct = figure;
    set(Hresconstruct,'WindowStyle','docked');
    Ireconstruct = Retroprojection_shared(I , X_m , z_out ,theta ,M0,Hresconstruct);
    ylim([0 Prof])
    cb = colorbar;
    ylabel(cb,'a.u')
    colormap(parula)
    set(findall(Hresconstruct,'-property','FontSize'),'FontSize',15) 

%%
figure;
imagesc(x-mean(x),z*1e3,1e3*Datas)
xlabel('\theta(°)')
ylabel('z(mm)')
colormap(parula)
cb = colorbar ;
ylabel(cb, 'AC tension mV')
set(findall(gcf,'-property','FontSize'),'FontSize',15) 



%% view image FFT
FT = MyImage.fourier(MyImage.Data) ;

% convolutive filter
[MconvX,MconvY] = meshgrid(-500:500,-500:500);
Mconv = exp(-(MconvX.^2+MconvY.^2)/(0.5*5^2));
FilterData = conv2(Datas,Mconv,'same'); 
FilterData = sqrt(trapz(x,trapz(z,Datas.^2)))...
             *FilterData/sqrt(trapz(x,trapz(z,FilterData.^2)));



% filter high frequencies
a = 1 ;

% Mask1 = 1*(abs(MyImage.fz/1000) < a)'*1*(abs(MyImage.fx/1000) < a);
% Mask2 = 1*(abs(MyImage.fz/1000) > 0.002)'*1*(abs(MyImage.fx/1000) > 0.002);

figure;
imagesc(x,z,FilterData)
colormap(parula)
colorbar
figure;
imagesc(x,z,Datas-FilterData)
colormap(parula)
colorbar

% figure;
% imagesc(MyImage.fx/1000,MyImage.fz/1000,abs(FT))
% axis(a*[-1 1 -1 1])
% xlabel('f_x (mm^{-1})')
% ylabel('f_z (mm^{-1})')
% title('Fourier Transform')
% cb = colorbar;
% ylabel(cb,'a.u')
% colormap(parula)