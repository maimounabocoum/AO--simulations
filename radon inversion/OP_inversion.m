
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% load experiemental data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
addpath('functions')
load('experiement images - JB - test\OP0deg-2016-02-01_13-11.mat');
c = 1540 ; % sound velocity in m/s
% update of the OP structure
% X : theta 
% Y : monotonic t variable
N = 2^11;
MyImage = OP(N,data(:,:,1),X,Y,Param,c); 
MyImage.Show_R();    % show Radon transform (ie raw data)
MyImage.Fmax()       % maximum frequency sampling = 1/dt
Lobject = 1e-3;
Fc = 10/Lobject; % Lobject is the size of the object to detect. Using simple model (sinc function)
                   % we set it to kc = 100/Lobject
                   
%% Nyist principle states the sampling of the object to reconstruct to be such that w > w_max/2 

%%% Fourier Tranform of Radon input image with respect to t

% creating a fourier parameter set using the measure sample size in m
% fourier transform of image

MyImage.F_R = MyImage.fourier(MyImage.R);
MyImage.Show_F_R(Fc); % Fc : cut-off frequency used for screening

%% image show

% representation in polar coordinates:

figure;
[THETA, W] = meshgrid(MyImage.theta*pi/180,MyImage.w(N/2:end));
[X,Y] = pol2cart(THETA, W);
surfc(X,Y,log(abs(MyImage.F_R(N/2:end,:))))
view(0,90)
shading interp
xlabel('\omega\it_{x} (\itm^{-1})')
ylabel('\omega\it_{y} (\itm^{-1})')
title('Fourier Transform of Radon')


% filtered inverse fourier transform :
% defninition of the filter :
% 
FILTER = abs(MyImage.w)'*ones(1,length(MyImage.theta));
 
 I = MyImage.ifourier(MyImage.F_R);
% 
% 
% I = MyImage.ifourier(F_R.*FILTER);
% 
xsonde=linspace(0,192*200e-6,193);
[x,y]=meshgrid(MyImage.t,xsonde);
img = zeros(size(x,1),size(x,2),'like',x);

  for i=1:length(MyImage.theta)
      t = x.*sin(MyImage.theta(i)*pi/180) + y.*cos(MyImage.theta(i)*pi/180);
      projContrib = interp1(MyImage.t',I(:,i),t(:),'linear',0);
      img = img + reshape(projContrib,length(xsonde),length(MyImage.t)); 
      
      clf;figure(15)
      imagesc(xsonde*1e3,MyImage.t*1e3,img)
      title(['angle',num2str(MyImage.theta(i))])
      %caxis([-400 400])
      colorbar
      drawnow 

  end


rmpath('functions')

