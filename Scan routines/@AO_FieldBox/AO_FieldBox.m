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
        
        function I_pointIndex = GetIndexPoint(obj,point)
            % retreive value using normal distance
            % center is the closest value to gepmetrical center
   
            List = Points(obj);
            distance = sum(abs(  List - repmat( point , size(List,1),1 ) ).^2,2);
            I_pointIndex = find( distance == min(distance) ) ;
            I_pointIndex = I_pointIndex(1); % get single value if exist equidistant points
            
            
            
        end
        
        function center = GetCenter(obj)
           center = [mean(obj.x), mean(obj.y),mean(obj.z)];           
        end
        
        function [Nx,Ny,Nz] = SizeBox(obj)
            Nx = length(obj.x);
            Ny = length(obj.y);
            Nz = length(obj.z);
        end
        
        function obj = Get_SimulationResults(obj,t,h,fs)
            %t: min time
            obj.Field = h;
            obj.time = t + (0:(size(h,1)-1))/fs;
            
        end
                
        function Iout = ShowMaxField(obj,plane,FigHandle)
            
            if ~ishandle(FigHandle)
                FigHandle = figure ;
            end            
            % set(FigHandle,'WindowStyle','docked');      
            % set(FigHandle,'NextPlot', 'replace') ;           
            
            [Nx,Ny,Nz]       = SizeBox(obj);
            % Field : column = time (Npoint) same as obj.time
            %         line = number of point
            [Field_max,Tmax] = max(abs(obj.Field),[],1);
            % max(obj.Field,[],1) : returns for each colulm
            % the maximum field pressure.
            Field_max = reshape(Field_max,[Ny,Nx,Nz]);
            Tmax = reshape(obj.time(Tmax),[Ny,Nx,Nz]);
            
            switch plane
                
                case 'Xt'
            %set(FigHandle,'name','(XT) field values for fixed z');
            
            
            output = envelope(obj.Field,300).^2;
            %Enveloppe = envelope(obj.MySimulationBox.Field,300).^2;
            
            % size(output,1) = number of temporal points
            % size(output,2) = number of spatial points
            % creation of a 3D matrix
            Field_max = reshape(output',[Ny,Nx,Nz,length(obj.time)]);    
            
            Iout = squeeze(Field_max(1,:,1,:))';
            
                imagesc(obj.x*1e3,obj.time*1e6,Iout)
                xlabel('x (mm)')
                ylabel('t (\mu s)')
                title(['P(t), z(t)= ',num2str((obj.z(2))*1e3),'mm']);
                
                case 'Zt'
            %set(FigHandle,'name','(XT) field values for fixed z');
            
            F = TF_t( 2^13 , 1/(obj.time(2)-obj.time(1)) );
            output = obj.Field;
            %Enveloppe = envelope(obj.MySimulationBox.Field,300).^2;
            
            % size(output,1) = number of temporal points
            % size(output,2) = number of spatial points
            % creation of a 3D matrix
            myField = reshape(output',[Ny,Nx,Nz,length(obj.time)]);    
            
            % get index of field in the center of the probe
            Nextract = 92;
            Iout = squeeze( myField(1,Nextract,:,:))';
            
            subplot(121)
            
                imagesc(obj.time*1e6,obj.z*1e3,Iout')
                ylabel('z (mm)')
                xlabel('t (\mu s)')
                title(['P(t), z(t)= ',num2str((obj.x(Nextract))*1e3),'mm']);
                cb = colorbar;
                ylabel(cb,'Pressure a.u.')
                
                % interpolation minus propagation
                for i=1:size(Iout,2)
                   Pressure(:,i) = interp1(obj.time - mean(obj.time),hilbert(Iout(:,i)),F.t,'linear',0);  
                   Spectrum(:,i) = F.fourier( Pressure(:,i) );
                end
                
           subplot(122)
           
                imagesc(F.f/(3e6),obj.z*1e3,abs(Spectrum)')
                %xlim([5 20])
                ylabel('z (mm)')
                xlabel('t (\mu s)')
                title(['P(t), z(t)= ',num2str((obj.x(Nextract))*1e3),'mm']);
                cb = colorbar;
                ylabel(cb,'temporal |FFT| a.u.')
                
                case 'XZt'

            set(FigHandle,'name','(XZ) maximum field (t) values');
            D=[1 1 1;0 0 1;0 1 0;1 1 0;1 0 0;];
            F=[0 0.25 0.5 0.75 1];
            G=linspace(0,1,256);
            cmap = interp1(F,D,G);
            
            %output = envelope(obj.Field,300).^2 ;
            output = obj.Field;
            myField = reshape(output',[Ny,Nx,Nz,length(obj.time)]);     
            
            Nskip = max(1,floor(size(obj.Field,1)/100)) ;
            m = max(abs( squeeze(myField(:))'));
            for loop = 1:Nskip:size(obj.Field,1) % loop over timefloor(size(obj.Field,1)/2) %
  
                %imagesc(obj.x*1e3,obj.z*1e3,abs( squeeze(myField(1,:,:,loop))')/m );
                imagesc(obj.x*1e3,obj.z*1e3,real( squeeze(myField(1,:,:,loop))')/m );
                xlabel('x (mm)')
                ylabel('z (mm)')
                ylim([min(obj.z*1e3) max(obj.z*1e3)])
                % title(['P(t) on XZ, z(t)= ',num2str((obj.time(i))*1540*1e3),'mm']) 
                colormap(cmap)
                colorbar
                %caxis([0.2*min(Field_max(:))/m 0.7*max(Field_max(:))/m])
                drawnow
               
               % caxis('auto')                
               % saveas(gcf,['gif folder\image',num2str(i),'.png'],'png')
            end
                  
               
                case 'XZ'
                    % selection of the interpolation plane:
                if (Ny == 1)
                    I_plane = 1;
                else
                    prompt = {'Enter Y coordinate (program will look for closest value):'};
                    dlg_title = 'Y plane select (mm)';
                    num_lines = 1;
                    answer = inputdlg(prompt,dlg_title,num_lines,{'0'});
                    V_plane = str2double(answer{1})*1e-3;                    
                    I_plane = Closest(V_plane,obj.y); 
                end
 
                    set(FigHandle,'name','(XZ) maximum field values') ;
                    subplot(121)
                    imagesc(obj.x*1e3,obj.z*1e3,squeeze(Field_max(I_plane,:,:))');
                    %view([0,90])
                    shading interp
                    xlabel('x (mm)')
                    ylabel('z (mm)')
                    title(['Maximum Field in plane Y = ',num2str(obj.y(I_plane)*1e3),'mm'])
                    colorbar
                    drawnow
                    subplot(122)
                    imagesc(obj.x*1e3,obj.z*1e3,1e6*squeeze(Tmax(I_plane,:,:))');
                    %view([0,90])
                    shading interp
                    xlabel('x (mm)')
                    ylabel('z (mm)')
                    title(['Time at Maximum Y = ',num2str(obj.y(I_plane)*1e3),'mm'])
                    cb = colorbar ;
                    ylabel(cb,'arrival time (\mu s)')
                    
                    drawnow
           
               case 'YZ'
                   if (Nx == 1)
                   I_plane = 1;
                else
                    prompt = {'Enter X coordinate (program will look for closest value):'};
                    dlg_title = 'X plane select (mm)';
                    num_lines = 1;
                    answer = inputdlg(prompt,dlg_title,num_lines,{'0'});
                    V_plane = str2double(answer{1})*1e-3;                    
                    I_plane = Closest(V_plane,obj.x); 
                   end

                    set(FigHandle,'name','(YZ) maximum field values')
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

