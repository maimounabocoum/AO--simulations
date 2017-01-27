classdef LaserBeam
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        w0
        
    end
    
    methods
        function obj = LaserBeam()
            w0 = 10e-3;
        end
        
        function [] = ShowLaserBeam(obj)
%            figure;
%             x = linspace(-10*obj.w0,10*obj.w0,500);
%             y = linspace(-10*obj.w0,10*obj.w0,500);
%            imagesc(x*1e3,y*1e3,obj.I)
%            xlabel('x (mm)')
%            ylabel('y (mm) : US probe position')
%            title(['Laser intensity Profile'])
%            colorbar
        end
    
    end

end

