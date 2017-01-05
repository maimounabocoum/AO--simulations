classdef ExcitationField
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N % number of actuator taken into account
        t % time
        Excitation % Field for each actuator (1 actuator = 1 line)
        
    end
    
    methods
        %f0 : defines the central frequency of the accoustique wave:
        function obj =  ExcitationField(f0,fs,N)
            % by default : sinusoidal wave
            
            % single excitation :
            obj.t = 0:1/fs:2/f0; % dt = 1/fs : sampling time
            
              % identical excitation wave for N actuators 
            obj.Excitation = repmat(sin(2*pi*f0*obj.t),N);
            
            
        end
        
                function obj =  Gaussian(f0,fs,N)
            % by default : sinusoidal wave
            
            % single excitation :
            obj.t = 0:1/fs:2/f0; % dt = 1/fs : sampling time
            
              % identical excitation wave for N actuators 
            obj.Excitation = repmat(sin(2*pi*f0*obj.t),N);
            
            
        end
        
    end
    
end

