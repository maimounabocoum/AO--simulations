function [angle, M0] = EvalDelayLaw( X_m , Z_m , ActiveLIST  )
%creation of function 11/09/2017 for delay law retreival using sequence
%parameters. 

% ct = M0 M(t)
Nangle = size(Z_m , 1) ;
% 0 is defined by the (0,0) on probe linear plane

  Hf = figure;
  set(Hf,'WindowStyle','docked');
  subplot(211)
    cc = jet(Nangle);
    colormap(jet);
    
    for i = 1:Nangle     
       % find index of minimum active element :
       Nmin = min(find(ActiveLIST(:,i) == 1));
       Nmax = max(find(ActiveLIST(:,i) == 1));
        

       angle(i) = atan( (Z_m(i,Nmin)-Z_m(i,Nmax))/(X_m(Nmax) - X_m(Nmin)) );
       % definition ut vector :
       M0(i,:) = [ X_m(Nmin)  , Z_m(i,Nmin)  ];
       plot( X_m*1e3 , Z_m(i,:)*1e3 ,'linewidth',3,'color',cc(i,:)) 
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
xlabel('angle (°)')
    








end

