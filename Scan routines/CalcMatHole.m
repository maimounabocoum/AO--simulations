function [nuX,nuZ,t,Mat] = CalcMatHole(f0,nbX,nbZ,nuX0,nuZ0,x,Fe,c,Bascule)
% c input in m/s
% f0  (en MHz) est la fréquence de porteuse
% nuZ (en mm-1)
% nuX (en mm-1)
% Fe : sampling frequency of Aixplorer in MHz ; 
% x  : coordinate vector in mm

% code modifié le 27-07-2019 par Maïmouna :
% objectif : n'envoyer qu'une seule période de modulation en z

% conversion of all input to SI units : 
Fe = Fe*1e6 ;           % MHz->Hz
f0 = f0*1e6;            % MHz->Hz
nuZ = 1e3*(nbZ*nuZ0) ;  % mm-1->m-1
nuX = 1e3*(nbX*nuX0);   % mm-1->m-1
x = x*1e-3;             % mm->m

% conversion to temporal frequency along z
fz = c*nuZ ; %Hz

dt = 1/Fe ; % in s
% Tmax = (20.1*1e-6);   % periode maximale d'un cycle élémentaire (en s)

Tz = 1/fz;          % periode de l'envelloppe A in s
T0 = 1/f0;          % periode de la porteuse A in s

% si la période de la modulation de phase est plus grande que Tmax --> Manip impossible
% (à remplacer par une exception..)
% if (Tz>Tmax) 
%    Mat = zeros(size(X));
%    return;
% end;

Nrep = nbZ;
%Nrep = floor(Tmax/Tz);    % on essaie de se rapprocher au mieux de Tmax 
                          % (essentially N = 1 for fundamental)
                          
% N = round(Nrep*(Tz/dt));  % Nombre de point effectifs
N = nbZ*(Fe/fz);

Tot = N*dt;        % durée totale de la séquence

Tz = Tot/Nrep;    % ré-ajustement de la période ??
fz = 1/Tz   ;     % ré-ajustement de la fréquence (en MHz) ??
nuZ = fz/c ;

k = round(Tot/T0); % nombre de cycles porteuse
f0 = k/Tot; % ré-ajustement de la fréquence porteuse

t = (0:N-1)*dt;  % time in us

[X,T] = meshgrid(x-mean(x),t);

alpha = nuX/fz;
carrier = sin(2*pi*f0*T);

if strcmp(Bascule,'on')
    
        Am          =  mod(ceil(2*fz*(T-alpha*X)),4);
        Am(Am==2)   =  0;
        Am(Am==3)   =  -1;
        Mat = sign(carrier).*Am;
        
else
    
        if nbZ==0
           Mat = sign(carrier) ;
        else
           Mat = sign(carrier).*( sin( 2*pi*fz*(T-alpha*X) )> 0 );   
        end    
        
end
% convert m-1->mm-1
nuZ = 1e-3*nuZ;
nuX = 1e-3*nuX;

%% print matrix
% figure(100)
%         imagesc(Mat)
%         xlabel('x (mm)')
%         ylabel('z(mm)')
%         drawnow
end
    