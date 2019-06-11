clearvars
REF = double( importdata('NBz1-10 object/I2000mA_US0.tiff') );
        NbZ         = 1:30;        % 8; % Nb de composantes de Fourier en Z, 'JM'
        NbX         = 0;        % 20 Nb de composantes de Fourier en X, 'JM'
        [NBX,NBZ] = meshgrid(NbX,NbZ);

        
figure(1)
imagesc(REF)
colorbar

%% definition of Fourier transform in filter plane (pixel 3.5 um)
Ny = 2^10;
Nx = 2^11;
dx = 3*0.0031/242; % calibration using 3 periods at 10NBz
dy = 3*0.0031/242;
dFx = 13.0208; % m^-1
dFy = 32.4675; % m^-1
G = TF2D(Nx,Ny,(Nx-1)*dFx,(Ny-1)*dFy);
F = TF2D(Nx,Ny,1/dx,1/dy);


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


%% get object of Fourier dimensions
Io = zeros(Ny,Nx);
Io(1:size(REF,1),1:size(REF,2)) = REF - REF(1,1) ;


figure(2)
imagesc(F.x,F.y,Io)
colorbar
caxis([min(REF(:)) max(REF(:))])

Ifft = F.fourier(Io);

% figure
% test = 1 - exp(-(F.y-2e-3).^2/(1e-3)^2)- exp(-(F.y+2e-3).^2/(1e-3)^2) ;
% plot(F.y,test)
% figure
% plot(F.fy,abs( fftshift(fft(test)) ) )

figure(3)
plot(F.y*1e3,Io(:,600))
colorbar
caxis([0 15000])

Ifft(Ny/2+1,Nx/2+1 ) = 0;

figure(8)
hold on
phaseref = unwrap(angle(Ifft(:,Nx/2+1 ))-angle(Ifft(Ny/2+1,Nx/2+1 )));
phaseref = phaseref-phaseref(Ny/2+1);

plot(F.fy,phaseref-phaseref(Ny/2+1),'-o')
%plot(F.fy,abs(Ifft(:,Nx/2+1 ))./max(abs(Ifft(:,Nx/2+1 ))),'-o')
%xlim([-1000 1000])


Io_taggedPhotons = F.ifourier(Ifft.*H);

figure(4)
imagesc(abs(Io_taggedPhotons).*Htagged)
%axis([0 1028 0 960])
colorbar
%caxis([0 100])

%%
figure(2)
imagesc((-Nx/2:1:Nx/2-1),(-Ny/2:1:Ny/2-1),abs(Ifft).*H)
colorbar
%caxis([0 500])
%axis([-10 10 -10 10])






%% filtrage dans fourier
 
%  Ifft = G.fourier(Io);
%  I_fft_filtered = 0*Ifft;
%  
%  
%  for i =1:length(NBZ(:))
%  I_fft_filtered(Ny/2+1 + NBZ(i) ,Nx/2+1 + NBX(i)) = Ifft(Ny/2+1 + NBZ(i) , Nx/2+1 + NBX(i));
%  I_fft_filtered(Ny/2+1 - NBZ(i) ,Nx/2+1 - NBX(i)) = conj(Ifft(Ny/2+1 + NBZ(i) ,Nx/2+1 + NBX(i))) ;
%  end
% 
% figure(8)
% hold on
% plot(F.fy,abs(I_fft_filtered(:,Nx/2+1 ))./max(abs(I_fft_filtered(:,Nx/2+1 ))),'-o')
% xlim([-500 400])
% % imagesc((-Nx/2:1:Nx/2-1),(-Ny/2:1:Ny/2-1),abs(I_fft_filtered))
% % axis([-10 10 -10 10])
% 
% figure(9)
% Io_filtered = G.ifourier(I_fft_filtered);
% imagesc(Io_filtered.*Htagged)
% caxis([0 2000])
% axis([0 1280 0 960])
% colorbar
