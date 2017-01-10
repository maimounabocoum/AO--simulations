%% Gaussian tumors sript

clearvars;
close all

WindowSizeInPixels = [100,100]; % x, y in mm
PositionOfTumors = [0,0;]; % list of tumor positions (x,y) in mm
                                 % each line corresponds to a tumor
SizeOfTumors = [10];          % size at \sigma in mm pixels
%TypeOfTumors = 'gaussian'; % gaussian, plain, square
TypeOfTumors = 'square';

I = TumorImage(WindowSizeInPixels,PositionOfTumors,SizeOfTumors,TypeOfTumors);
I.ShowTumor();

%% sample de tumor with sampling frequency fs :

fs = 1/(5*SizeOfTumors); % in mm-1

I.sampleTumor(fs);
I.FourierTumor();


% %%% radon transform of the image %%%
% theta = 0:1:20;
% 
% [R,xp] = radon(I.I,theta);
% 
% %% inverse of radon transform
% %I_R = iradon(R,theta,'linear','ram-lak');%frequency_scaling,output_size)
% %I_R = iradon(R,theta,'linear','cosine');%frequency_scaling,output_size)
% %I_R = iradon(R,theta,'linear','hamming');%frequency_scaling,output_size)
% %I_R = iradon(R,theta,'linear','hann');%frequency_scaling,output_size)
% %I_R = iradon(R,theta,'linear','shepp-logan');%frequency_scaling,output_size)
% I_R = iradon(R,theta,'linear','none')%frequency_scaling,output_size)
% 
% %% radon transform of inverse Radon
% [R_IR,xp_IR] = radon(I_R,theta);
% 
% figure;
%     subplot(2,2,1)
% imagesc(I.x,I.y,I.I);
%     subplot(2,2,2)
% imagesc(theta,xp,R);
% title('R_{\theta} (X\prime)');
% xlabel('\theta (degrees)');
% ylabel('X\prime');
% set(gca,'XTick',0:20:180);
% colormap(hot);
% colorbar
%     subplot(2,2,3)
% imagesc(I_R)
%    subplot(2,2,4)
% imagesc(theta,xp_IR,R_IR);
% title('R_{\theta} (X\prime)');
% xlabel('\theta (degrees)');
% ylabel('X\prime');
% set(gca,'XTick',0:20:180);
% colormap(hot);
% colorbar
