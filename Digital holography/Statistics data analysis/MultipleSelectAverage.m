% clearvars

[Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
Nfiles_ = length(Filename);
if Nfiles_==1
   Filename = {Filename};
end
Nfiles = 24;
%% unfolder average

load('LogFile.mat');
%LogFile = table2array(LogFile);

%% 


    
sumImage = zeros(1,Nfiles);
Ny = 2^10;
Nx = 2^11;
dx = 3.5e-6;
dy = 3.5e-6;
F = TF2D(Nx,Ny,1/((Nx-1)*dx),1/((Ny-1)*dy));
Io = zeros(Ny,Nx,Nfiles);
% define filter in FFT
H = ones(Ny,Nx);
BOX = [350 750 100 260];
H( : ,(-Nx/2:1:Nx/2-1) > BOX(2) )  = 0;
H( :,(-Nx/2:1:Nx/2-1) < BOX(1)  )  = 0;
H( (-Ny/2:1:Ny/2-1) > BOX(4) , : ) = 0;
H( (-Ny/2:1:Ny/2-1) < BOX(3), :  ) = 0;

%% filter for tagged photons
Htagged = ones(Ny,Nx);
Htagged( : ,(1:Nx) > 990 | (1:Nx) < 100 )  = 0 ;
Htagged( (1:Ny) < 270 | (1:Ny) > 650  , : )  = 0;

%% loop
for n_file  = 1:Nfiles
    
REF(:,:,n_file) = 0*double( importdata([Foldername,Filename{1}]) ) ;

Iscan = n_file:24:Nfiles_;
%Io_temp = zeros(F.Ny,F.Nx,20);

for Iaverage =  1:length(Iscan)
    
REF_ = double( importdata([Foldername,Filename{Iscan(Iaverage)}]) ) ;

% Ref_temp(:,:,Iaverage) = REF_;
% Io_temp(1:960,1:1280,Iaverage)= REF_;
% Ifft_temp(:,:,Iaverage) = F.fourier(Io_temp(:,:,Iaverage));
% Io_taggedtemp(:,:,Iaverage) = F.ifourier(Ifft_temp(:,:,Iaverage).*H);

REF(:,:,n_file)  = REF(:,:,n_file) + REF_;

end

% ToProbe = abs(Io_taggedtemp) ;
% 
% MU = mean(ToProbe,3);
% VAR = var(ToProbe,0,3);
% VAR(VAR==0) = 1e23;
% figure(5)
% imagesc(Htagged.*( MU./sqrt(VAR) ) )
% colorbar
% caxis([0 10])

REF(:,:,n_file) = REF(:,:,n_file)/20;



% sumImage(n_file) = sum(sum(REF(:,:,n_file)))/(960*1280);

Io(1:size(REF(:,:,n_file),1),1:size(REF(:,:,n_file),2),n_file) = REF(:,:,n_file) ;

InputSum(n_file) = sum(sum(Io(:,:,n_file))) ;

% nomalization image

% alpha = trapz(F.x,trapz(F.y,( Io(:,:,n_file) ).^2));

Io(:,:,n_file) = Io(:,:,n_file)/InputSum(n_file)   ;

InputSum2(n_file) = sum(sum(Io(:,:,n_file))) ;


Ifft(:,:,n_file) = F.fourier(Io(:,:,n_file));

% filter
Io_taggedPhotons(:,:,n_file) = F.ifourier(Ifft(:,:,n_file).*H);
Io_taggedPhotons(:,:,n_file) = Io_taggedPhotons(:,:,n_file).*Htagged;
end       






figure(2)
imagesc(abs(Io_taggedPhotons(:,:,1)).*Htagged)
colorbar
caxis([0 0.1e-8])
%% fourier reconstruction

Ny = 2^10;
Nx = 2^11;
dFx = 13.0208;
dFy = 32.4675;
G = TF2D(Nx,Ny,(Nx-1)*dFx,(Ny-1)*dFy);
Gfft = zeros(Ny,Nx);





%% result analysis

ResultImage = zeros(size(Io_taggedPhotons,1),size(Io_taggedPhotons,2),size(Io_taggedPhotons,3)/4);
ResultProjection = zeros(1,size(Io_taggedPhotons,3)/4);
RawSum = zeros(1,Nfiles);

for n_file  = 1:Nfiles
    
 ResultImage(:,:,LogFile(n_file,1)) = ResultImage(:,:,LogFile(n_file,1)) + ...
     abs(Io_taggedPhotons(:,:,n_file)).*exp(1i*2*pi*LogFile(n_file,4));

 RawSum(n_file) = sum(sum(abs(Io_taggedPhotons(:,:,n_file))));
 
 ResultProjection(LogFile(n_file,1)) = ResultProjection(LogFile(n_file,1)) + ...
     RawSum(n_file).*exp(1i*2*pi*LogFile(n_file,4));
 
 
 Gfft(Ny/2+1+LogFile(n_file,3),Nx/2+1+LogFile(n_file,2)) = ResultProjection(LogFile(n_file,1));
 Gfft(Ny/2+1-LogFile(n_file,3),Nx/2+1-LogFile(n_file,2)) = conj( ResultProjection(LogFile(n_file,1)) );
 

end

% figure(3)
% hold on
% plot(  RawSum )



% for i=1:(Nfiles/4)
% 
% figure(1);
% subplot(1,6,i)
% imagesc( imag( ResultImage(:,:,i) ) )
% axis([0 1280 0 960])
% colorbar
% caxis([-50 50])
% drawnow   
% 
% end

figure(1);
imagesc(abs(Gfft))
colorbar
axis([Nx/2+1-15, Nx/2+1+15, Ny/2+1-15, Ny/2+1+15])
%caxis([0 1000000])
title('NbZ = 5:10, NbX = 0')


%% view results
Irecons = G.ifourier( Gfft ) ;

figure(5)
imagesc(Irecons)
title('object reconstruction')
% %%
% for i=1:6
% figure(5)
% imagesc( real(ResultImage(:,:,i)) )
% colorbar
% %caxis([-1e-9 1e-9])
% axis([0 1280 0 960])
% pause(1)
% end























