


figure('DefaultAxesFontSize',18);  

for nloop = 1:size(Field_Profile,3)

nloop = 2
figure('DefaultAxesFontSize',18);  
imagesc(x_phantom*1e3,z_phantom*1e3,Field_Profile(:,:,nloop))
colorbar
title(['N_x =',num2str(MyImage.decimation(nloop))])
xlabel('x (mm)')
ylabel('y (mm)')
caxis([(1-1e-3)*min(min(Field_Profile(:,:,nloop))) (1+1e-3)*max(max(Field_Profile(:,:,nloop)))])
drawnow   

saveas(gcf,'C:\Users\mbocoum\Dropbox\untitled.png')

end