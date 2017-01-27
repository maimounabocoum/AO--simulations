%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clear all
clearvars ;
parameters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Simulation BOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation box initialization : 

    Nx = 200;
    Ny = 1;
    Nz = 200;

    Xrange = [-15 15]/1000; % in m
    Yrange = 0/1000;%[-0.1 0.1]/1000;
    Zrange = [0 70]/1000; % in m
    
    % field 
    % gaussian : w0 , f0 (transducer freq) , fs (sampling frequency)
    % plane : angle of propagation : theta

SimulationBox = AO_FieldBox(Xrange,Yrange,Zrange,Nx,Ny,Nz);
%focus = [0 0 40]/1000; % Initial electronic focus
%Rfocus = 40/1000; % Elevation focus
SimulationBox = SimulationBox.GenerateField(N_elements*width,focus(3),f0,fs,c);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% calculation of the emitted Field :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SimulationBox.ShowMaxField('XZ');
