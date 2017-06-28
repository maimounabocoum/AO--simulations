% example from user manual p.27
clearvars;
addpath(genpath('D:\GIT\AO---softwares-and-developpement\KWAVE'));

%% create the computational grid

Nx = 128;       % number of grid points in the x (row) direction
Ny = 256;       % number of grid points in the y (column) direction
dx = 50e-6;     % grid point spacing in the x direction [m]
dy = 50e-6;     % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);



%% define the medium properties
medium.sound_speed = 1500*ones(Nx, Ny); % [m/s]
medium.sound_speed(1:50, :) = 1800;     % [m/s]
medium.density = 1040;                  % [kg/m^3]

%% define an initial pressure using makeDisc
% predefined option for sources :
% makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);
% source mask for time varying pressure initial condition
disc_x_pos = 75;     % [grid points]
disc_y_pos = 120;    % [grid points]
disc_radius = 8;     % [grid points]
disc_mag = 0.01e6;   % [Pa]
source.p_mask = false(Nx, Ny);
source.p_mask(50,:) = true;
sampling_freq = 50e6;        % [Hz]
tone_burst_freq = 1e6;       % [Hz]
tone_burst_cycles = 3;
source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles);

% Hinitsource = figure;
%     set(Hinitsource,'WindowStyle','docked');
%     imagesc(source.p0 )
%     xlabel('x (px)')
%     ylabel('y (px)')
% %     axis equal
% %     axis tight
%     title('initalization of pressure field')
%     cb = colorbar;
%     ylabel(cb,'Pascal')
%     colormap(parula)
%     set(findall(Hinitsource,'-property','FontSize'),'FontSize',15) 
    
%% define a Cartesian sensor mask of a centered circle with 50 sensor elements
sensor_radius = 2.5e-3; % [m]
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

%% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);