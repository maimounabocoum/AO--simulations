function Ireconstruct = Retroprojection_shared(I,X_m,z_out,theta, M0 ,H )
% function created by maimouna bocoum 13/09/2017

z_out = z_out(:)';

size(I,1)
length(z_out)
% check consistancy of data dimensions : 
if size(I,1)~=length(z_out)
    error('inconsistent data size')
end

% generation of mask function M + t ut
% followed by interpolation on grid 


% retroprojection : 


 [X,Z]= meshgrid(X_m,z_out);
 Ireconstruct = zeros(size(X,1),size(X,2),'like',X);

%  H = figure;
%  A = axes ;
 set(H,'WindowStyle','docked');
 
  for i= 1:length(theta)
       
         
        T =   (X - M0(i,1)).*sin( theta(i) ) ...
            + (Z - M0(i,2)).*cos( theta(i) ) ;
      % common interpolation:  
        %Mask = double( interp1(X_m,ActiveLIST(:,i),X,'linear',0) );
        
        projContrib = interp1(z_out,I(:,i),T(:),'linear',0);
         
     
       % retroprojection:  
        Ireconstruct = Ireconstruct + reshape(projContrib,length(z_out),length(X_m)); 
        %%% real time monitoring %%%   
       imagesc( X_m*1e3,z_out*1e3,Ireconstruct)
       colormap(parula)
       cb = colorbar ;
       title(['angle(°): ',num2str(theta(i)*180/pi)])
       xlabel('x (mm)')
       ylabel('z (mm)')
       caxis( [ min(Ireconstruct(:)) , max(Ireconstruct(:)) ] )
       
       %saveas(gcf,['Q:\AO---softwares-and-developpement\radon inversion\gif folder/image',num2str(i),'.png'])
       drawnow 
  end
  

    
    %title('Reconstruction')
    ylabel(cb,'AC tension (mV)')
    colormap(parula)
    set(findall(H,'-property','FontSize'),'FontSize',15) 



end









