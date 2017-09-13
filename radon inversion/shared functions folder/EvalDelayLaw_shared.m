function [angle, M0] = EvalDelayLaw( X_m , Z_m  )
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
       hold on
       angle(i) = atan( (Z_m(i,1)-Z_m(i,end))/(X_m(end) - X_m(1)) );
       % definition ut vector :
      % ut(i,:)= [sin(angle(i)) , cos(angle(i))] ; % (x,z)
       M0(i,:) = [ X_m(1)  , Z_m(i,1)  ];
       plot( X_m*1e3 , Z_m(i,:)*1e3 ,'linewidth',3,'color',cc(i,:)) 
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

