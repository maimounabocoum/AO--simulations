function [] = Retroprojection( I , X_m, z_out , theta, M0)
% function created by maimouna bocoum 13/09/2017

z_out = z_out(:)';

% check consistancy of data dimensions : 
if size(I,1)~=length(z_out)
    error('inconsistent data size')
end




% retroprojection : 


 [X,Z]= meshgrid(X_m,z_out);
 Ireconstruct = zeros(size(X,1),size(X,2),'like',X);
 
 H = figure;
 A = axes ;
 set(H,'WindowStyle','docked');
 
  for i= 1:length(theta)
      
        T =   (X - M0(i,1)).*sin( theta(i) ) ...
            + (Z - M0(i,2)).*cos( theta(i) ) ;
      % common interpolation:  
        projContrib = interp1(z_out,I(:,i),T(:),'linear',0); %13.9765
     
     % retroprojection:  
        Ireconstruct = Ireconstruct + reshape(projContrib,length(z_out),length(X_m)); 
        
      %%% real time monitoring %%%   
       imagesc( X_m*1e3,z_out*1e3,Ireconstruct,'parent',A)
       colormap(parula)
       cb = colorbar ;
       title(['angle(°): ',num2str(theta(i)*180/pi)])
       xlabel('x (mm)')
       ylabel('z (mm)')
       caxis( [ min(Ireconstruct(:)) , max(Ireconstruct(:)) ] )
       drawnow 

        
  end

    
    %title('Reconstruction')
    ylabel(cb,'AC tension (mV)')
    colormap(parula)
    set(findall(H,'-property','FontSize'),'FontSize',15) 
  

end









