classdef ExcitationField
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        Probe
    end
    
    properties
        t % time
        Excitation % Field for each actuator (1 actuator = 1 line)   
    end
    
    methods
        
        %f0 : defines the central frequency of the accoustique wave:
        function obj =  ExcitationField()
            % default probe
        obj.Probe = ActuatorProbe(192,6/1000,0.2/1000,1,10,0,1:192,35/1000);
            % default excitation field
        fs = 100e6 ;
        f0 = 2e6 ;
        obj.t = (0:1/fs:6/f0);
        obj.Excitation=  sin(2*pi*f0*obj.t).*hanning(length(obj.t)).^2';  

        end
        
        
    end
    
end

