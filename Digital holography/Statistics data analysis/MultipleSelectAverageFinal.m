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
dx = 3*0.0031/242; % calibration using 3 periods at 10NBz
dy = 3*0.0031/242;
dFx = 13.0208; % m^-1
dFy = 32.4675; % m^-1
G = TF2D(Nx,Ny,(Nx-1)*dFx,(Ny-1)*dFy);
F = TF2D(Nx,Ny,1/dx,1/dy);
Io = zeros(Ny,Nx);

%% define filter in FFT
H = ones(Ny,Nx);
BOX = [1100 1300 700 800];
H( : ,1:Nx > BOX(2) )  = 0;
H( :,1:Nx < BOX(1)  )  = 0;
H( 1:Ny > BOX(4) , : ) = 0;
H( 1:Ny < BOX(3), :  ) = 0;

%% filter for tagged photons
Htagged = ones(Ny,Nx);
Htagged( : ,(1:Nx) > 1200 | (1:Nx) < 0 )  = 0 ;
Htagged( (1:Ny) < 0 | (1:Ny) > 950  , : )  = 0;

Io_tagged_ = zeros(Ny,Nx,10,20);

% loop
for n_file  = 1:40

%   Io_tagged_ = 0*Io_tagged_;
    Iscan = n_file:40:Nfiles;
    
    for Iaverage  = 1:length(Iscan)
    
REF = double( importdata([Foldername, Filename{Iscan(Iaverage)} ]) ) ;

Io(1:size(REF,1),1:size(REF,2)) = REF ;

InputSum = sum(sum(Io)) ;

Io = Io/InputSum   ;

Ifft = F.fourier(Io);
% filter
Io_taggedPhotons = G.ifourier(Ifft.*H);

Io_taggedPhotons = Io_taggedPhotons.*Htagged;

Io_tagged_(:,:,LogFile(n_file,1),Iaverage) = Io_tagged_(:,:,LogFile(n_file,1),Iaverage)...
    + abs(Io_taggedPhotons).*exp(1i*2*pi*LogFile(n_file,4)) ;

RawSum(n_file) = sum(sum( abs( Io_taggedPhotons ) ));


 
    end

end

%%
for i = 1: 10
  MU( : , : )  =  mean(  Io_tagged_(:,:,i,:)  , 4 ) ;
  
  for j = 1:size(Io_tagged_,4)
  cf(i,j) =    sum(sum( Io_tagged_(:,:,i,j)));
  end
  
  cf_mean(i) = mean(abs(cf(i,:))); % average over sum
  cf_std(i) = sqrt( var(abs(cf(i,:))) ); % sum over average
  
  VAR( : , : ) =  var( abs(Io_tagged_(:,:,i,:)), 0 , 4 ) ;
 
 
  figure(4)
   % plot(F.y*1e3,real(MU( : , 550 ))/max(real(MU( : , 550 ))))
   test_fft = fftshift(fft(real(MU( : , 550 ))/max(real(MU( : , 550 )))));
   max_fft(i) = test_fft( max(find(abs(test_fft) == max(abs(test_fft)))) );
   max_freq(i) =  F.fy(max(find(abs(test_fft) == max(abs(test_fft)))));
   hold on
   plot(F.y*1e3,angle(test_fft))
   title('NBz = 3 , \phi_2 = 0.15 [2 \pi Rad]')
%  title(['SNR (real(4-phases)) over 20 samples NBz =',num2str(i)])
%  %axis([1 1300 250 700])
%  colorbar
%  caxis([0 5])
%  drawnow
%  %saveas(gcf,['mu_NBz',num2str(i),'.png'])

end  

% figure(2)
% errorbar(cf_mean,cf_std)

%% result analysis
%% fourier reconstruction


Gfft = zeros(Ny,Nx);

for n_file  = 1:size(Io_tagged_,3)

 MU( : , : )  =  mean( Io_tagged_(:,:,n_file,:)  , 4 ) ;
 cf(n_file) =    sum(sum( MU( : , : ) )).*exp(1i*4*n_file) ;
 Gfft(Ny/2+1+n_file,Nx/2+1) = cf(n_file);
 Gfft(Ny/2+1-n_file,Nx/2+1) = conj( cf(n_file) );
 
 % view results
Irecons = G.ifourier( Gfft ) ;


end

%%
figure(2)
imagesc(G.x*1e3,G.y*1e3,Irecons)

figure(5)
hold on
plot(G.y*1e3,Irecons(:,Nx/2+1)/max(Irecons(:,Nx/2+1)))
title('object reconstruction')
xlabel('z(mm)')

figure(8)
hold on
phase = unwrap(angle( Gfft(:,Nx/2+1) - Gfft(Ny/2+1,Nx/2+1)));
phase = phase-phase(Ny/2+1);
plot(G.fy,  phase ,'-o')

plot(G.fy, 0.3*abs( Gfft(:,Nx/2+1) )/max(abs( Gfft(:,Nx/2+1) ) ) ,'-o')


figure(1);
imagesc(abs(Gfft))
colorbar
axis([Nx/2+1-15, Nx/2+1+15, Ny/2+1-15, Ny/2+1+15])
%caxis([0 1000000])
title('NbZ = 5:10, NbX = 0')




% for i=1:6
% figure(5)
% imagesc( real(ResultImage(:,:,i)) )
% colorbar
% %caxis([-1e-9 1e-9])
% axis([0 1280 0 960])
% pause(1)
% end























