% clearvars

[Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end

%% unfolder average

load('LogFile.mat');
%LogFile = table2array(LogFile);

%% 
  
sumImage = zeros(1,Nfiles);
Ny = 2^10;
Nx = 2^11;
dx = 3.5e-6; % calibration using 3 periods at 10NBz
dy = 3.5e-6;
F = TF2D(Nx,Ny,1/dx,1/dy);
Io = zeros(Ny,Nx);

% define filter in FFT
H = ones(Ny,Nx);
BOX = [1100 1300 700 800];
H( : ,1:Nx > BOX(2) )  = 0;
H( :,1:Nx < BOX(1)  )  = 0;
H( 1:Ny > BOX(4) , : ) = 0;
H( 1:Ny < BOX(3), :  ) = 0;

% filter for tagged photons
Htagged = ones(Ny,Nx);
Htagged( : ,(1:Nx) > 1200 | (1:Nx) < 0 )  = 0 ;
Htagged( (1:Ny) < 0 | (1:Ny) > 950  , : )  = 0;

Io         = zeros(Ny,Nx,Nfiles);
Ifft       = zeros(Ny,Nx,Nfiles);
Io_tagged_ = zeros(Ny,Nx,Nfiles);

%% loop
for n_file  = 1:Nfiles
  
REF = double( importdata([Foldername, Filename{n_file} ]) ) ;

Io(1:size(REF,1),1:size(REF,2),n_file) = REF ;

Ifft(:,:,n_file) = F.fourier(Io(:,:,n_file));

Io_taggedPhotons(:,:,n_file) = F.ifourier(Ifft(:,:,n_file).*H);

end

%%

  MU( : , : )  =  mean(  Io/16 , 3 )   ;
  VAR( : , : ) =  var( Io/16 , 0 ,3 ) ;
  VAR(VAR==0) = 1e-23;
  


figure(1)
imagesc(MU)  
caxis([min(REF(:))/16 max(REF(:))/16])
axis([0 1280 0 960])
cb = colorbar ;
ylabel(cb,'counts (12 bit)');
title(['Average \mu over ',num2str(Nfiles)])
saveas(gcf,'figure1_mu.png')

figure(2)
imagesc(sqrt(VAR))
axis([0 1280 0 960])
cb = colorbar ;
ylabel(cb,'counts (12 bit)');
title(['STD \sigma over ',num2str(Nfiles)])
saveas(gcf,'figure2_std.png')

figure(3)
imagesc(16*VAR./MU)
axis([0 1280 0 960])
caxis([0 20])
colorbar
title(['Poisson compare = sigma^2 / \mu over ',num2str(Nfiles)])
saveas(gcf,'figure3_poisson.png')

figure(4)
imagesc(MU./sqrt(VAR))
axis([0 1280 0 960])
caxis([0 1500])
colorbar
title(['SNR = \mu / \sigma over ',num2str(Nfiles)])
saveas(gcf,'figure4_snr.png')
%%

