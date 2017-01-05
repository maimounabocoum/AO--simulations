%  Read data to make a CFM image from simulated data
%
%  Version 1.0, March 24, 1996, Joergen Arendt Jensen

%  Physical data

f0=3.5e6;                %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
Ncfm=10;                 %  Number of pulses in one direction
fprf=10e3;               %  Pulse emissions frequency  [Hz]
D=5;                     %  Sampling frequency decimation rate
est_dist=1/1000;         %  Distance between velocity estimates [m]
no_lines=20;             %  Number of imaging directions
dy=2/1000;               %  Distance between imaging lines [m]
v_est=0;                 %  Set estimates to zero [m/s]

%  Load the data for one image line and a number of
%  pulse emissions

for i=1:no_lines
  min_sample=0;
  data=0;

  for k=1:Ncfm
    cmd=['load sim_flow/rft',num2str(k),'l',num2str(i),'.mat']
    eval(cmd)
  
    %  Decimate the data and store it in data
  
    if (tstart>0)
      rf_sig = [zeros(tstart*fs-min_sample,1); rf_data];
    else
      rf_sig = rf_data( abs(tstart*fs):max(size(rf_data)) ) ;
      end
     
    rf_sig=hilbert(rf_sig(1:D:max(size(rf_sig))));
    data(1:max(size(rf_sig)),k)=rf_sig;
    end

  %  Make the velocity estimation

  Ndist=floor(2*est_dist/c*fs/D);  %  Rf samples between velocity estimates
  [Nsamples,M]=size(data);
  index=1;

  RMS=std(data(:,1));
  for k=1:Ndist:Nsamples

    %  Find the proper data
  
    vdata=diff(data(k,:));

    %  Calculate the autocorrelation and the velocity

    if (std(vdata) > RMS/5)
      auto  = vdata(2:(M-1)) * vdata(1:(M-2))' ;
      v_est(index,i) = c*fprf/(4*pi*f0) * atan2(imag(auto),real(auto));
    else
      v_est(index,i)=0;
      end
    index=index+1;
    end
  end

[Nx,Ny]=size(v_est);  
imagesc(((0:Ny-1)-Ny/2)*dy*1000,(1:Nx)*Ndist*D/fs*c/2*1000,v_est)
map=[1:64; zeros(2,64)]/64;
colormap(map')
colorbar
ylabel('Depth in tissue [mm]')
xlabel('Lateral distance [mm]')
drawnow

%  Make an interpolated image

ID=25;
[n,m]=size(v_est)
new_est1=zeros(n,m*ID);
for i=1:n
i
  new_est1(i,:)=abs(interp(v_est(i,:),ID));
  end
[n,m]=size(new_est1)
new_est=zeros(n*5,m);
Ndist=Ndist/5;
for i=1:m
i
  new_est(:,i)=abs(interp(new_est1(:,i),5));
  end

[Nx,Ny]=size(new_est);  
new_est=new_est/max(max(new_est))*64;
imagesc(((0:Ny-1)-Ny/2)*dy/ID*1000,(1:Nx)*Ndist*D/fs*c/2*1000,new_est)
map=[1:64; zeros(2,64)]/64;
colormap(map')
colorbar
ylabel('Depth in tissue [mm]')
xlabel('Lateral distance [mm]')
axis('image')

