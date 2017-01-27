classdef AO_FieldBox
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        x
        y
        z
        
        time
        Field
       
    end
    
    methods
        function obj = AO_FieldBox(Xrange,Yrange,Zrange,Nx,Ny,Nz)
            
            obj.x = linspace(Xrange(1),Xrange(end),Nx);
            obj.y = linspace(Yrange(1),Yrange(end),Ny);
            obj.z = linspace(Zrange(1),Zrange(end),Nz);
            
            
            
        end
        
        function obj = GenerateField(obj,w0i,z_x,f0,fs,c)
            
            % set simulation time length :
            dt = 1/fs; % sampling time 
            tmax = max(abs(obj.z))/c;
            obj.time = 0:dt:tmax;
            
             Noc = 10; % number of optical cycles
             excitation =  sin(2*pi*f0*obj.time);
             excitation = excitation.*hanning(length(excitation))';
            
            % defining w0 :
            lambda = c/f0 ;
            w0 = sqrt( (1/2)*( w0i^2 + sqrt(w0i^4 + 4*(z_x^2*lambda^2)/pi^2  ) ) ) ;
            Zr = pi*w0^2/lambda;
            
            [X,Z] = meshgrid(obj.x,obj.z);
            
            % define spatial profile :
            
            W = w0*sqrt(1 + Z.^2./Zr) ;
            obj.Field = (w0./W).*exp( - 2*X.^2./W.^2) ;
            
            for i = 1:10:length(obj.time)
            clf;figure(10);
            excitation =  sin(2*pi*f0*(obj.time(i) - Z/c)).*exp(-(obj.time(i) - Z/c).^2/(8/f0)^2);
            imagesc(1e3*obj.x,1e3*obj.z,abs(obj.Field.*excitation).^2)
            title(['time = ',num2str(1e6*obj.time(i)),'\mu s'])
            xlabel('x (mm)')
            ylabel('z (mm)')
            colorbar
            drawnow           
            end

            
        end
        
        function List = Points(obj)
            
           [X,Y,Z] = meshgrid(obj.x,obj.y,obj.z);  

            List = [X(:),Y(:),Z(:)];
     
        end
        
        function [Nx,Ny,Nz] = SizeBox(obj)
            Nx = length(obj.x);
            Ny = length(obj.y);
            Nz = length(obj.z);
        end
        
        function obj = Get_SimulationResults(obj,t,h,fs)
            
            obj.Field = h;
            obj.time = t + [0:(size(h,1)-1)]/fs;
            
        end
        
        function [] = ShowMaxField(obj,plane)
            [Nx,Ny,Nz] = SizeBox(obj);
            % get maximum field value (with respect to time)   
            Field_max= reshape(max(obj.Field,[],1),[Ny,Nx,Nz]);
            switch plane
                case 'XZ'
                    % selection of the interpolation plane:
                if (Ny == 1)
                    I_plane = 1;
                else
                    prompt = {'Enter Y coordinate (program will look for closest value):'};
                    dlg_title = 'Y plane select (mm)';
                    num_lines = 1;
                    answer = inputdlg(prompt,dlg_title,num_lines,{'0'});
                    V_plane = str2num(answer{1})*1e-3;                    
                    I_plane = Closest(V_plane,obj.y); 
                end
 


            Hf3 = figure(3);
            set(Hf3,'name','(XZ) maximum field values')
            %imagesc(SimulationBox.x*1e3,SimulationBox.z*1e6,squeeze(Field_max(1,:,:)));
            size(squeeze(Field_max(:,I_plane,:))')
            imagesc(obj.x*1e3,obj.z*1e3,squeeze(Field_max(I_plane,:,:))');
            shading interp
            xlabel('x (mm)')
            ylabel('z (mm)')
            title(['Maximum Field in plane Y = ',num2str(obj.y(I_plane)*1e3),'mm'])
            colorbar
            
               case 'YZ'
                   if (Nx == 1)
                   I_plane = 1;
                else
                    prompt = {'Enter X coordinate (program will look for closest value):'};
                    dlg_title = 'X plane select (mm)';
                    num_lines = 1;
                    answer = inputdlg(prompt,dlg_title,num_lines,{'0'});
                    V_plane = str2num(answer{1})*1e-3;                    
                    I_plane = Closest(V_plane,obj.x); 
                   end
                            Hf4 = figure(4);
            set(Hf4,'name','(YZ) maximum field values')
            %imagesc(SimulationBox.x*1e3,SimulationBox.z*1e6,squeeze(Field_max(1,:,:)));
            size(squeeze(Field_max(:,I_plane,:))')
            imagesc(obj.y*1e3,obj.z*1e3,squeeze(Field_max(:,I_plane,:))');
            shading interp
            xlabel('y (mm)')
            ylabel('z (mm)')
            title(['Maximum Field in plane X = ',num2str(obj.x(I_plane)*1e3),'mm'])
            colorbar

            end
        
        end
        
    end
    
end

function I_plane = Closest(V_plane,x)
    Distance = abs(x - V_plane);
    I_planeList = find(Distance == min(Distance));
    I_plane = I_planeList(1);

end

