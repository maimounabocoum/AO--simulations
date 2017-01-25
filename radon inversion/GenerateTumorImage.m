%% Gaussian tumors sript

clearvars;
close all

WindowSizeInPixels = [40,300]; % x, y in mm
PositionOfTumors = [20,40]; % list of tumor positions (x,y) in mm
                                 % each line corresponds to a tumor
SizeOfTumors = [10];          % size at \sigma in mm pixels
TypeOfTumors = 'gaussian'; % gaussian, plain, square
%TypeOfTumors = 'square';

I = TumorImage(WindowSizeInPixels,PositionOfTumors,SizeOfTumors,TypeOfTumors);
%I.ShowTumor();
% simulate the radon scan in angles : 
I = I.ScanTumor(linspace(-20,20,50));

%% sample de tumor with sampling frequency fs :


