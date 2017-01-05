%  Create a computer model of the tissue surrounding a
%  vessel
%
%  Calling: [positions, amp] = tissue_pht (N);
%
%  Parameters:  N - Number of scatterers in the phantom
%
%  Output:      positions  - Positions of the scatterers.
%               amp        - amplitude of the scatterers.
%
%  Version 1.0, March 26, 1997 by Joergen Arendt Jensen

function [positions, amp] = tissue_pht (N)

%  Dimensions of the phantom

x_size = 40/1000;   %  Width of phantom [mm]
y_size = 10/1000;   %  Transverse width of phantom [mm]
z_size = 90/1000;   %  Height of phantom [mm]
z_plus = z_size/2 + 10/1000;  %  Start of phantom surface [mm];

%  Initialize the ranges for the vessel

R=0.005;         %  Radius of blood vessel [m]
y_range=2*R;     %  y range for the scatterers  [m]
z_range=2*R;     %  z range for the scatterers  [m]
z_offset=0.06;   %  Offset of the mid-point of the scatterers [m]

%  Creat the general scatterers

x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = (rand (N,1)-0.5)*z_size + z_plus - z_offset;

%  Generate the amplitudes with a Gaussian distribution

amp=randn(N,1);

%  Generate the rotated and offset block of sample

theta=-45/180*pi;
xnew=x*cos(theta)+z*sin(theta);
znew=z*cos(theta)-x*sin(theta);

%  Make the vessel and set the amplitudes to -40 dB below inside

inside = (( y.^2 + (znew).^2) < R^2);
amp = amp .* (1-inside) + amp .* inside/100*0; 

%  Generate the rotated and offset block of sample

theta=45/180*pi;
x=xnew*cos(theta)+znew*sin(theta);
z=znew*cos(theta)-xnew*sin(theta) + z_offset;

%  Return the variables

positions=[x y z];
