%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main  program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clear all
clearvars ;

addpath('..\Field_II')
field_init(0);

parameters;


 P = ActuatorProbe(N_elements,element_height,width,no_sub_x,no_sub_y,kerf,ActiveList);
% rect=[1 0/1000 0/1000 0 2/1000 0/1000 0 2/1000 5/1000 0 0/1000 5/1000 0 ...
% 1 2/1000 5/1000 1/1000 2.5/1000 0
% 1 -2/1000 0/1000 0 0/1000 0/1000 0 0/1000 5/1000 0 -2/1000 5/1000 0 ...
% 1 2/1000 5/1000 -1/1000 2.5/1000 0
% 1 -2/1000 -5/1000 0 0/1000 -5/1000 0 0/1000 0/1000 0 -2/1000 0/1000 0 ...
% 1 2/1000 5/1000 -1/1000 -2.5/1000 0
% 1 0/1000 -5/1000 0 2/1000 -5/1000 0 2/1000 0/1000 0 0/1000 0/1000 0 ...
% 1 2/1000 5/1000 1/1000 -2.5/1000 0];
% center=[0 0 0];
% focus=[0 0 70]/1000;
%Th= xdc_rectangles (rect, center, focus);


%Probe = xdc_linear_array (N_elements, width, element_height, kerf,no_sub_x,no_sub_y, focus);
Probe = xdc_rectangles(P.rect,[0 0 0], focus);

show_xdc (Probe);