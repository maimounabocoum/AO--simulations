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

%============ plot impulse response
            % calculate impulse response in FIELD II
                t_impulseResponse = (0:1/param.fs:2/param.f0);
                impulse = sin(2*pi*param.f0*t_impulseResponse);
                impulse=impulse.*hanning(length(impulse))'; 
                
                plot(t_impulseResponse*1e6,impulse)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%xdc_free(Probe);
field_end;
