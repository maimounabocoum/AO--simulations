%% Gaussian tumors sript

clearvars;
close all

WindowSize = [40,50]*1e-3; % x, y in mm
PositionOfTumors = [20,40]*1e-3; % list of tumor positions (x,y) in mm
                                 % each line corresponds to a tumor
SizeOfTumors = [10]*1e-3;        % size at \sigma in mm pixels
TypeOfTumors = 'gaussian'; % gaussian, plain, square
%TypeOfTumors = 'square';

ph = Phantom();
laser = LaserBeam();
%I = Phantom(WindowSizeInPixels,PositionOfTumors,SizeOfTumors,TypeOfTumors);
ph.ShowTumor();
laser.ShowLaserBeam();
% simulate the radon scan in angles : 
%I = I.ScanTumor(linspace(-20,20,50));

%% sample de tumor with sampling frequency fs :




