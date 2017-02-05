%%% this is a test routine %%%
%% maimouna bocoum 04-01-2017
clear all
%clearvars ;

addpath('..\Field_II')
field_init(0);
parameters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate aperture for emission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = ActuatorProbe(param.N_elements,param.element_height,param.width,...
                  param.no_sub_x,param.no_sub_y,param.kerf,param.ActiveList,param.Rfocus);
%P.ShowProbe()
%Probe = xdc_rectangles(P.rect,[0 0 0], param.focus);
Probe = xdc_focused_array(param.N_elements,param.width,param.element_height,...
         param.kerf,param.Rfocus,param.no_sub_x,param.no_sub_y,param.focus);
%data2 = xdc_get(Probe,'rect')
show_xdc (Probe)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%xdc_free(Probe);
field_end;
