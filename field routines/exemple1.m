%% Example of transducer divided with a  2 x 3 sampling
% Summary of example objective
%% maimouna bocoum 04-01-2017
clear all
clearvars ;

f0=6e6; % Transducer center frequency [Hz]
fs=200e6; % Sampling frequency [Hz]
N = 3;
% box initialization : 
Nx = 1;
Ny = 5;
Nz = 5;

Xrange = 0;
Yrange = [-10 10]/1000; % in m
Zrange = [0 50]/1000; % in m

SimulationBox = AO_FieldBox(Xrange,Yrange,Zrange,Nx,Ny,Nz);

excitation = ExcitationField(f0,fs,N);

