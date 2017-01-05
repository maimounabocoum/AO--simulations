clear all
close all
clc

field_init(0);
CurrentDirectory = cd;

% Set initial parameters for Scatterers
f0=5.5e6;           % Transducer center frequency [Hz]
fs=200e6;           % Sampling frequency [Hz]
c=1540;             % Speed of sound [m/s]
lambda=c/f0;        % Wavelength [m]

N_elements = 1;
focus = [0 0 40]/1000;
R_elem = 3.5/1000;
el_size = 0.4/1000;

Th = xdc_concave(R_elem, focus(3), el_size);

impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
figure; plot(impulse)

excitation = sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (Th, impulse);
figure; plot(excitation)

% Scatterers Definition

x_phantom_size = 3/1000;
y_phantom_size = 3/1000;
z_phantom_size = 20/1000;
z_phantom_start = 30/1000;
N_scat = 2000;

x = (rand(N_scat,1)-0.5)*x_phantom_size;
y = (rand(N_scat,1)-0.5)*y_phantom_size;
z = rand(N_scat,1)*z_phantom_size + z_phantom_start;

phantom_amplitudes = randn(N_scat,1);
phantom_positions = [x y z];
clear x y z



% Diplacement field

%load('Dpl_tot.mat')
load('Dir1.mat');
load('Dir2.mat');
load('Dir3.mat');
load('InfoPressureField3D_6.mat');
InfoBeam.pasX = InfoBeam.pasX/1000;
InfoBeam.pasY = InfoBeam.pasY/1000;
InfoBeam.pasZ = InfoBeam.pasZ/1000;
InfoBeam.XDebut = InfoBeam.XDebut/1000;
InfoBeam.YDebut = InfoBeam.YDebut/1000;
InfoBeam.ZDebut = InfoBeam.ZDebut/1000;

%j=0;
%test1 = 0;
disp('Data Loaded')
for timing = 1:1:size(Dpl_dir1,4)
    timing
    phantom_wave(:,:,timing) = phantom_positions(:,:);
    for X0 = 1:1:size(Dpl_dir1,1)
        X = InfoBeam.XDebut + (X0 - 1)*InfoBeam.pasX;
        for Y0 = 1:1:size(Dpl_dir1,2)
            Y = InfoBeam.YDebut + (Y0 - 1)*InfoBeam.pasY;
            for Z0 = 1:1:size(Dpl_dir1,3)
                Z = InfoBeam.ZDebut + (Z0 - 1)*InfoBeam.pasZ;
                %j = j + 1;
                A = find((phantom_positions(:,1)>=(X-InfoBeam.pasX/2)) & (phantom_positions(:,1)<(X+InfoBeam.pasX/2)) & (phantom_positions(:,2)>=(Y-InfoBeam.pasY/2)) & (phantom_positions(:,2)<(Y+InfoBeam.pasY/2)) & (phantom_positions(:,3)>=(Z-InfoBeam.pasZ/2)) & (phantom_positions(:,3)<(Z+InfoBeam.pasZ/2)));
                if length(A)>0
                    %test1 = test1 + 1;
                    %tail(j)=length(A);
                    for i=1:1:length(A)
                        phantom_wave(i,1,timing) =  phantom_positions(i,1) + Dpl_dir1(X0,Y0,Z0,timing);
                        phantom_wave(i,2,timing) =  phantom_positions(i,2) + Dpl_dir2(X0,Y0,Z0,timing);
                        phantom_wave(i,3,timing) =  phantom_positions(i,3) + Dpl_dir3(X0,Y0,Z0,timing);
                    end
                end
            end
        end
    end
end
                


 












