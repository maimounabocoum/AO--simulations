%% Algo de rétroprojetction filtrée
% clear
 load('PressureScan_OP_vis15mm_50V_6MhZ-2017-05-03_17-18.mat')
 DataOP0=data;
%  load('OP90deg-2016-02-01_13-38.mat')
%  DataOP90=data;
clear data

theta=X;
Proj=DataOP0(:,:,1);

SamplingRate=Param.SamplingRate*1e6;

thetar=pi*theta/180;

xsonde=linspace(0,192*0.2,193);
yt=(0:size(Proj,1)-1)*1540/(SamplingRate)*1000;

FFTt=fftshift(fft(Proj));
f=linspace(-SamplingRate/(2*1540),SamplingRate/(2*1540),length(yt));

[x,y]=meshgrid(xsonde,yt);

BPt=1.5; %MHz
BP=BPt*1e6/1540; %m-1

imgfilt=zeros(length(yt),length(xsonde),length(theta));
FFTtfilt=FFTt;
% LP=ones(1,length(yt));
% LP=abs(sinc(f/(2*BP)));
LP=cos(pi*f/(2*BP));
LP(abs(f)>BP)=0;

for i=1:size(FFTt,2)
    FFTtfilt(:,i)=abs(f)'.*FFTt(:,i).*LP';
end
 
Projfilt=abs(ifft(fftshift(FFTtfilt)));

for i=1:length(theta)
        t = x.*sin(thetar(i)) + y.*cos(thetar(i))...
        -0.5*(1+sign(thetar(i)))*abs(xsonde(end).*tan(thetar(i)))...
        +abs(xsonde(end).*tan(max(abs(X))*pi/180));
    projContrib = interp1(yt',Projfilt(:,i),t(:),'linear');
    projContrib(isnan(projContrib))=0;%Projfilt(end,i);
    imgfilt(:,:,i) = reshape(projContrib,length(yt),length(xsonde));
        figure(1)
        imagesc(xsonde,yt,imgfilt(:,:,i))
    drawnow
end
% 
% Proj=DataOP90(:,:,1);
% FFTt2=fftshift(fft(Proj));
% 
% imgfilt2=zeros(length(yt),length(xsonde),length(theta));
% FFTtfilt=FFTt;
% 
% for i=1:size(FFTt,2)
%     FFTtfilt(:,i)=abs(f)'.*FFTt(:,i).*LP';
% end
%  
% Projfilt=abs(ifft(fftshift(FFTtfilt)));
% 
% for i=1:length(theta)
%         t = x.*sin(thetar(i)) + y.*cos(thetar(i))...
%         -0.5*(1+sign(thetar(i)))*abs(xsonde(end).*tan(thetar(i)))...
%         +abs(xsonde(end).*tan(max(abs(X))*pi/180));
%     projContrib = interp1(yt',Projfilt(:,i),t(:),'linear');
%     projContrib(isnan(projContrib))=0;%Projfilt(end,i);
%     imgfilt2(:,:,i) = reshape(projContrib,length(yt),length(xsonde));
%         figure(1)
%         imagesc(xsonde,yt,imgfilt2(:,:,i))
%     
% end
% 
 IMG = sum(imgfilt(:,:,:),3)/size(imgfilt,3);
 figure
 imagesc(IMG)
% 
% IMGt=IMG-mean(mean((IMG(1:20,:))));
% IMGt=IMGt-min(min((IMGt)));
% IMGt=IMGt/(max(max(IMGt)));
% % IMGt(x.*sin(thetar(end)) + y.*cos(thetar(end))>40)=0;
% % IMGt(x.*sin(thetar(1)) + y.*cos(thetar(1))+abs(xsonde(end).*tan(thetar(1)))>40)=0;
% IMGt=IMGt(yt<40,:,:);
% IMG50V=IMGt;
% 
% clear IMGt
% 
% Yt=yt(yt<40);
% yt=Yt;
% clear Yt
% 
% f2=figure(2);
% imagesc(xsonde,yt(yt<40),IMG)
% axis image
% title('Reconstructed image',...
%     'Fontsize',30)
% xlabel('Lateral position (mm)','fontsize',36)
% ylabel('Longitudinal position (mm)','fontsize',36)
% set(gca,'fontsize',36,'linewidth',3,'box','on','color','w')
% set(gcf,'color','none')
% set(f2,'units','normalized','position',[0 0 1 1],'paperpositionmode','auto')
% % caxis([0 1])
