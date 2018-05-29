disp(' ');
disp(['fdct_usfft_demo_basic.m -- This demo displays a curvelet'])
disp (['both in the spatial and frequency domains.']);
disp(' ');
disp(['This is achieved by setting all the coefficients in the curvelet'])
disp(['domain to zero except that at the required location (which'])
disp(['is set to one). The curvelet is obtained by taking the'])
disp(['adjoint curvelet transform. Notice how the curvelet is sharply '])
disp(['localized in both space and frequency.']); 
disp(' ');

% fdct_usfft_demo_basic.m -- This demo displays a curvelet
clearvars;

m = 512;
n = 512;

X = zeros(m,n);

%forward curvelet transform
disp('Take curvelet transform: fdct_usfft');
tic; C = fdct_usfft(X,0); toc;

%specify one curvelet
j = 5;
%2^(ceil(j/2)+1)
l = 1;
[A,B] = size(C{j}{l});
k1 = ceil((A+1)/2);
k2 = ceil((B+1)/2);
C{j}{l}(k1,k2) = 1;

%adjoint curvelet transform
disp('Take adjoint curvelet transform: afdct_usfft');
tic; Y = afdct_usfft(C,0); toc;

%display the curvelet
F = ifftshift(fft2(fftshift(Y)));

figure;
subplot(1,2,1); colormap gray; imagesc(real(Y)); axis('image'); ...
    title('a curvelet: spatial viewpoint');
subplot(1,2,2); colormap gray; imagesc(abs(F)); axis('image'); ...
    title('a curvelet: frequency viewpoint');

%get parameters
[SX,SY,FX,FY,NX,NY] = fdct_usfft_param(C);

