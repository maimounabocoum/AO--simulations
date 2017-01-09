%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main  program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clear all
clearvars ;

addpath('..\Field_II')
field_init(0);

parameters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %definition of the actuator array %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P = ActuatorProbe(N_elements,element_height,width,no_sub_x,no_sub_y,kerf,ActiveList);
    %P.ShowActuatorCenter();

% type of waves for delays line : 1. focused to Z , 2. plane waves with
% angle theta 3. user defined delay law 

    %P = P.Set_ActuatorDelayLaw('focus',[0 0 40]/1000,c);
    P = P.Set_ActuatorDelayLaw('plane',10*180/pi,c);
    %P.ShowDelay();

    %Probe = xdc_linear_array (N_elements, width, element_height, kerf,no_sub_x,no_sub_y, focus);
    Probe = xdc_rectangles(P.rect,[0 0 0], focus);
    % show_xdc (Probe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Setting the impulse response field (Green function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Emission signal : I think by befault, the green function is calculated
% for a different frequency

    t_impulseResponse = (0:1/fs:2/f0);
    impulse = sin(2*pi*f0*t_impulseResponse);
    impulse=impulse.*hanning(length(impulse))'; 
    xdc_impulse (Probe, impulse);
    % figure;
    % plot(1e6*t_impulseResponse,impulse);
    % xlabel('time(\mu s)')
    % ylabel('a. u')
    % title('Probe impulse Response setting')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of excitation function for actuator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %excitation= Actuators.Excitation(1,:);%sin(2*pi*f0*(0:1/fs:2/f0));
    Noc = 2;
    t_excitation = (0:1/fs:Noc*1.5/f0);
    excitation =  sin(2*pi*f0*t_excitation);
    excitation = excitation.*hanning(length(excitation))';
    % addind delay law
    xdc_focus_times(Probe,-1,P.DelayLaw);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Simulation BOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation box initialization : 

    Nx = 50;
    Ny = 1;
    Nz = 40;

    Xrange = [-5 5]/1000; % in m
    Yrange = 0;%[-0.1 0.1]/1000;
    Zrange = [1 45]/1000; % in m

SimulationBox = AO_FieldBox(Xrange,Yrange,Zrange,Nx,Ny,Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% calculation of the emitted Field :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,t] = calc_hp(Probe,SimulationBox.Points());
h = h/max(h(:));
SimulationBox = Get_SimulationResults(t,h,fs);
SimulationBox.ShowMaxField();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdc_free(Probe);
field_end;

