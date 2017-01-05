
% Definition des sources
% Input: Structure contenant les limites de la zone a etudier
% Output : Position = Matrice des points de la zone
%          distance = distance max couverte 

function [Positions,Amplitudes,InfoBeam] = zoneImageRemote(InfoBeam)

XDebut  = InfoBeam.XDebut;
XFin    = InfoBeam.XFin;
pasX    = InfoBeam.pasX;
YDebut  = InfoBeam.YDebut;
YFin    = InfoBeam.YFin;
pasY    = InfoBeam.pasY;
ZDebut  = InfoBeam.ZDebut;
ZFin    = InfoBeam.ZFin;
pasZ    = InfoBeam.pasZ;



Sourcexi = [XDebut:pasX:XFin]*1e-3; NbXs = length(Sourcexi);
Sourceyi = [YDebut:pasY:YFin]*1e-3; NbYs = length(Sourceyi);
Sourcezi = [ZDebut:pasZ:ZFin]*1e-3; NbZs = length(Sourcezi);
NbPointss = NbXs*NbYs*NbZs;

Sourcez = ones(NbXs*NbYs,1)*Sourcezi; 
Sourcez = reshape(Sourcez,NbPointss,1);
Sourcex = ones(NbZs,1)*Sourcexi; 
Sourcex = reshape(Sourcex',NbXs*NbZs,1);
Sourcex = ones(NbYs,1)*Sourcex'; 
Sourcex = reshape(Sourcex,NbPointss,1);
Sourcey = ones(NbZs*NbXs,1)*Sourceyi; 
Sourcey = reshape(Sourcey',NbPointss,1);

Positions = [Sourcex';Sourcey'; Sourcez']'
Amplitudes = ones(1,NbPointss);
InfoBeam.NbX = NbXs; InfoBeam.NbZ = NbZs; InfoBeam.NbY = NbYs;