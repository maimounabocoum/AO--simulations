%% OF analysis %%
% addpath('..\shared functions folder');
% MyImage = OF(x*1e-3,z*1e-3,Datas,SampleRate*1e6,c) ;
% MyImage = OP(x*1e-3,z*1e-3,Datas,SampleRate*1e6,c) ;
% addpath('shared functions folder')

%%
figure;
imagesc(x,z,1e3*Datas)
xlabel('x(mm)')
ylabel('z(mm)')
colormap(parula)
cb = colorbar ;
ylabel(cb, 'AC tension mV')
set(findall(gcf,'-property','FontSize'),'FontSize',15) 

 [cx,cy,c] = improfile;
 figure;
 plot(cx(1) + sqrt((cx-cx(1)).^2 + (cy-cy(1)).^2),c/max(abs(c)))

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