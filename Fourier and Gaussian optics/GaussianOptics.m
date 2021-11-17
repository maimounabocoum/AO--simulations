%% some Gaussian function
clearvars;


%% calculate divergence from waist:
lambda = 1064e-9;   % optical wavelength
w0     = 3e-6;      % waist in focus

Zr = pi*w0^2/lambda; % Rayleigh length
theta = atan( lambda/(pi*w0) );

180*theta/pi
NA = 1*sin(theta)

4*30+67.69+4*26.56+23.21+192.80+195.74+4*22.13+4*28.03+3*64.42+118.03+ 63.44