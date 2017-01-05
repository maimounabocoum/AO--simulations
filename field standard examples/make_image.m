%  Compress the data to show 40 dB of
%  dynamic range for the CFM phantom image
%
%  version 1.2 by Joergen Arendt Jensen, April 6, 1997.

clearvars ;
addpath('..\Field_II')
field_init(0);

f0=3.5e6;                 %  Transducer center frequency [Hz]
fs=100e6;                 %  Sampling frequency [Hz]
c=1540;                   %  Speed of sound [m/s]
no_lines=50;              %  Number of lines in image
d_x=40/1000/no_lines;     %  Increment for image

%  Read the data and adjust it in time 

min_sample=0;
for i=1:no_lines

  %  Load the result

  cmd=['load sim_bmd/rf_ln',num2str(i),'.mat']
  eval(cmd)
  
  %  Find the envelope
  
  if (tstart>0)
    rf_env=abs(hilbert([zeros(tstart*fs-min_sample,1); rf_data]));
  else
    rf_env=abs(hilbert( rf_data( abs(tstart*fs):max(size(rf_data)) ) ));
    end
  env(1:max(size(rf_env)),i)=rf_env;
  end

%  Do logarithmic compression

D=20;   %  Sampling frequency decimation factor

log_env=env(1:D:max(size(env)),:)/max(max(env));
log_env=log(log_env+0.01);
log_env=log_env-min(min(log_env));
log_env=64*log_env/max(max(log_env));

%  Make an interpolated image

ID_bmode=10;
[n,m]=size(log_env)
new_env=zeros(n,m*ID_bmode);
for i=1:n
  if(rem(i,100) == 0)
    i
    end
  new_env(i,:)=abs(interp(log_env(i,:),ID_bmode));
  end
[n,m]=size(new_env)
  
fn_bmode=fs/D;
clg
image(((1:(ID_bmode*no_lines-1))*d_x/ID_bmode-no_lines*d_x/2)*1000,((1:n)/fn_bmode+min_sample/fs)*1540/2*1000,new_env)
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
colormap(gray(64))
%brighten(-0.35)
%axis([-20 20 30 90])
axis('image')

