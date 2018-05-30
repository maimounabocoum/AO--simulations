%% OF analysis %%
% addpath('shared functions folder') ;


%% view image FFT
FT = MyImage.fourier(MyImage.Data) ;

figure;
imagesc(MyImage.fx/1000,MyImage.fz/1000,abs(FT))
xlabel('f_x (mm^{-1})')
ylabel('f_z (mm^{-1})')