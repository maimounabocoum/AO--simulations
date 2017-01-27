%%% this is a test routine %%%
%% maimouna bocoum 04-01-2017
clear all
clearvars ;

addpath('..\Field_II')
field_init(0);
parameters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate aperture for emission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = ActuatorProbe(N_elements,element_height,width,no_sub_x,no_sub_y,kerf,ActiveList,Rfocus);
%P.ShowProbe()
Probe = xdc_rectangles(P.rect,[0 0 0], focus);
%Probe = xdc_linear_array (N_elements, width, element_height, kerf,no_sub_x,no_sub_y, focus);
%Probe = xdc_focused_array(N_elements,width,element_height,kerf,Rfocus,no_sub_x,no_sub_y,focus);
%data = xdc_get(Probe,'rect');
show_xdc (Probe)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%xdc_free(Probe);
field_end;
