%%% Load simulated Datas

%% add path
addpath('..\Field_II')
addpath('..\radon inversion')
addpath('subscripts')
addpath('..\..\AO--commons\shared functions folder')

%% load datas

load('Q:\datas\simulated datas\2020-06-01\TaboltEffect_off_type_JM_19h33_23.mat')

%% view correlated datas
Hf = figure;
CurrentExperiement.ShowFieldCorrelation('XZ', Hf , 20e-6 , 1);