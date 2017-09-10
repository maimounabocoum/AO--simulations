%% creation 03-07-2017 by maimouna Bocoum :
clearvars
filename = uigetfile('\\arsenal\mbocoum\datas\2017-07-01\*OP*.mat');
load(['\\arsenal\mbocoum\datas\2017-07-01\',filename])
c = 1540 ;
SampleRate = 10 ;

%% data initialization
%AlphaM= -AlphaM:dA:AlphaM ;
Nlines = length(AlphaM );
z = (1:size(raw,1))*(c/(1e6*SampleRate));
Datas = RetreiveDatas(raw,NTrig,Nlines,MedElmtList);

figure;
imagesc(AlphaM,z,Datas)

%% radon inversion
       Iradon = OPinversionFunction(AlphaM*pi/180,z,Datas,SampleRate*1e6,c);
