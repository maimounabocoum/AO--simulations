function [angle, M0 , X0 , Z0] = EvalDelayLawOS_shared( X_m , DelayLAWS , ActiveLIST , c )
%creation of function 11/09/2017 for delay law retreival using sequence
%parameters. 

% on fabrique la matrice les coordonnée initial du front d'onde 
% à partir de DelayLAWS de la taille de la sonde et des X_m, coordonnée en
% x local de la position des piezos élements coincident avec les coordonée
% du labo sur l'axe des x : (X_m,0)


% DelayLAWS : each column represents the delay law in s for all probe
% each line the angle of shoot
% Last dimension is the décimate value

% convert seconds to distance
ct =    DelayLAWS*c ;
 
% ct = M0 M(t)
Nangle = size(ct ,  2 ) ;

% initialization:
angle = zeros(1,Nangle); % angle emission of wavefront inittialization
M0    = zeros(Nangle,2); % initial position of wavefront inittialization

% 0 is defined by the (0,0) on probe linear plane

  Hf = figure;
  set(Hf,'WindowStyle','docked');
  subplot(211)
    cc = jet(Nangle);
    colormap(jet);
 
    
    for i = 1:Nangle     
       % find index of minimum active element :
   
       Nmin = find(ActiveLIST(:,i) == 1, 1 );
       Nmax = find(ActiveLIST(:,i) == 1, 1, 'last' );
        
       angle(i) = asin( (ct(Nmax,i)-ct(Nmin,i))/(X_m(Nmax) - X_m(Nmin)) );
              % definition ut vector :
       u = [sin(angle(i)),cos(angle(i))];


       
       
       % affichage des coordonnée initiales :
       % (X_m,0) - ct*ut
       X0(:)   = X_m - u(1)*DelayLAWS(:,i)'*c ;
       Z0(:)   = 0   - u(2)*DelayLAWS(:,i)'*c;
       M0(i,1) = 0 - u(1)*DelayLAWS(1,i)'*c; 
       M0(i,2) = 0   - u(2)*DelayLAWS(1,i)'*c;
       plot( X0(:)*1e3 , Z0(:)*1e3 ,'linewidth',3,'color',cc(i,:)) 
       hold on
       
    end
    
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

