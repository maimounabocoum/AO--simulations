classdef AO_FieldBox
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x
        y
        z
        
        X
        Y
        Z
        
        time
        Field
       
    end
    
    methods
        function obj = AO_FieldBox(Xrange,Yrange,Zrange,Nx,Ny,Nz)
            
            obj.x = linspace(Xrange(1),Xrange(end),Nx);
            obj.y = linspace(Yrange(1),Yrange(end),Ny);
            obj.z = linspace(Zrange(1),Zrange(end),Nz);
            
          [X,Y,Z] = meshgrid(obj.x,obj.y,obj.z);  
          
          obj.X = X(:);
          obj.Y = Y(:);
          obj.Z = Z(:);

        
        end
        
        function List = Points(obj)
           
            List = [obj.X,obj.Y,obj.Z];
     
        end
        
        function obj = Get_SimulationResults(obj,t,h,fs)
            
            Field = h;
            time = t + [0:(size(h,1)-1)]/fs;
            
        end
        
        function [] = ShowMaxField(obj)
        
        % get maximum field value (with respect to time)    
        Field_max= reshape(max(obj.Field,[],1),[length(obj.x),length(obj.y),length(obj.z)]);
        
        Hf3 = figure(3);
        set(Hf3,'name','maximum field values')
        %imagesc(SimulationBox.x*1e3,SimulationBox.z*1e6,squeeze(Field_max(1,:,:)));
        imagesc(obj.x*1e3,obj.z*1e3,squeeze(Field_max(:,1,:))');
        xlabel('x (mm)')
        ylabel('z (mm)')
        colorbar
        
        end
        
    end
    
end

