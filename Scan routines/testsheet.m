%% test sheet
clearvars
parameters;



MyProbe = ActuatorProbe(param.N_elements,param.element_height,param.width,...
                          param.no_sub_x,param.no_sub_y,param.kerf,param.ActiveList,param.Rfocus);
                      
MyProbe.ShowProbe();

MyProbe.GetIndex([1:20]/1000)