%% OF analysis %%
addpath('..\shared functions folder');
% MyImage = OF(x*1e-3,z*1e-3,Datas,SampleRate*1e6,c) ;
% MyImage = OP(x*1e-3,z*1e-3,Datas,SampleRate*1e6,c) ;
% addpath('shared functions folder')

%%
    MyImage = OP(Datas,Alphas,z,SampleRate*1e6,c) ;
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
imagesc(Alphas*180/pi,z*1e3,1e3*Datas)
xlabel('x(mm)')
ylabel('z(mm)')
colormap(parula)
cb = colorbar ;
ylabel(cb, 'AC tension mV')
set(findall(gcf,'-property','FontSize'),'FontSize',15) 

[cx,cy,c] = improfile;
figure;
plot(sqrt((cx-cx(1)).^2 + (cy-cy(1)).^2),c)

