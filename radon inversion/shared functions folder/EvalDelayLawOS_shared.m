function [angle, M0 , X0 , Z0 , C] = EvalDelayLawOS_shared( X_m , DelayLAWS , ActiveLIST , c )
%creation of function 11/09/2017 for delay law retreival using sequence
%parameters. 

% on fabrique la matrice les coordonnée initial du front d'onde 
% à partir de DelayLAWS de la taille de la sonde et des X_m, coordonnée en
% x local de la position des piezos élements coincident avec les coordonée
% du labo sur l'axe des x : (X_m,0)

% C is the fixed point. Value [0,0] is returned for theta = 0

% DelayLAWS : each column represents the delay law in s for all probe
% each line the angle of shoot
% Last dimension is the décimate value

% view result of not
ScreenResult = 0;

% convert seconds to distance
ct =    DelayLAWS*c ;
 
% ct = M0 M(t)
Nangle = size(ct ,  2 ) ;

% initialization:
angle = zeros(1,Nangle); % angle emission of wavefront inittialization
M0    = zeros(Nangle,2); % initial position of wavefront inittialization
C     = zeros(Nangle,2);
% 0 is defined by the (0,0) on probe linear plane

if ScreenResult ==1
  Hf = figure;
  set(Hf,'WindowStyle','docked');
  subplot(211)
    cc = jet(Nangle);
    colormap(jet);
end
 
    
    for i = 1:Nangle     
       % find index of minimum active element :
   
       Nmin = find(ActiveLIST(:,i) == 1, 1 );
       Nmax = find(ActiveLIST(:,i) == 1, 1, 'last' );
       
 
       angle(i) = asin( (ct(Nmax,i)-ct(Nmin,i))/(X_m(Nmax) - X_m(Nmin))  );
       % truncation of angle by 6 unit decimate
       angle(i) = round(angle(i)*100000)/100000 ;
       
       % definition ut vector orthogonal to initial (linear) wavefront:
       % equation of previous wavefront : y = 0
       % equation of new wavefront : y = ax + b
       % intersection equation : x = -b/a if a not equals to 0 , set to 0
       % otherwise
       
       u = [sin(angle(i)),cos(angle(i))];
       
       % affichage des coordonnée initiales :
       % (X_m,0) - ct*ut
       X0   = X_m - u(1)*DelayLAWS(:,i)'*c ;
       Z0   = 0   - u(2)*DelayLAWS(:,i)'*c;
       M0(i,1) = 0 - u(1)*DelayLAWS(1,i)'*c; 
       M0(i,2) = 0   - u(2)*DelayLAWS(1,i)'*c;
       if ScreenResult ==1
       plot( X0(:)*1e3 , Z0(:)*1e3 ,'linewidth',3,'color',cc(i,:)) 
       end
       if ~(Z0(end)-Z0(1)) == 0

%        a = (Z0(end)-Z0(1))/(X0(end)-X0(1));
%        b = -(Z0(end)*X0(1)-Z0(1)*X0(end))/(X0(end)-X0(1));
       C(i,1) = (Z0(end)*X0(1)-Z0(1)*X0(end))/(Z0(end)-Z0(1));
       C(i,2) = 0;
       end
       % evaluate fixed point over rotation
       if ScreenResult ==1
       hold on
       plot( C(i,1)*1e3 , C(i,2)*1e3 ,'x','color','black') 
       end
       
    end
    
    if ScreenResult ==1
        
    cb = colorbar ;
    ylabel(cb,'angular index (a.u.)')
    xlabel('X (mm)')
    ylabel('Z(mm)')
    set(gca,'YDir','reverse')
    
subplot(212)
plot(180/pi*angle,'-o','linewidth',3)
xlabel('shoot index')
ylabel('angle (°)')

    end
    






end

